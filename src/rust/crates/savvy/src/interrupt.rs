//! User-interrupt polling that never unwinds through Rust.
//!
//! The solvers poll between iterations, on the thread that entered the FFI and
//! never inside a parallel region, and they only need an answer: the core
//! unwinds normally, reports an `interrupted` flag, and `balance()` re-signals
//! with `rlang::interrupt()` once the call returns. Detection is therefore the
//! whole job here. `R_CheckUserInterrupt` delivers an interrupt by longjmp,
//! which would skip Rust destructors, so it is never called bare.
//!
//! On Unix the detection is two reads, mirroring R's own test.
//! `R_interrupts_pending` is the flag R's `SIGINT` handler sets, and
//! `R_interrupts_suspended` gates it the way R gates it: while interrupts are
//! suspended, R records the interrupt and returns rather than jumping, so a
//! poll must report nothing pending. Reporting one anyway would stop the solve
//! and then hand `balance()` a stopped fit that `rlang::interrupt()` declines
//! to signal, which would fall through to the non-convergence path. R declares
//! both flags `LibExtern`, so both are exported. Both are read volatile and
//! neither is ever written.
//!
//! Reading the flags pushes nothing on R's context stack, which is the reason
//! for the platform split. The earlier implementation wrapped
//! `R_CheckUserInterrupt` in `R_ToplevelExec` on every platform to contain the
//! longjmp, and R's profiler does not survive a `CTXT_TOPLEVEL` context pushed
//! inside a `.Call`. Under `profvis` on macOS a `SIGPROF` sample landed in the
//! pop window, at the instruction after `Rf_endcontext` returned and before
//! `R_ToplevelContext` was restored; R's `findProfContext` walked a null
//! context pointer and the session died with a segmentation fault inside the
//! entropy solver. That window is well under a microsecond, so the crash was
//! rare, but every solver iteration of every profiled fit was exposed to it.
//!
//! Windows keeps the `R_ToplevelExec` path, because there a user's interrupt
//! is registered by `R_ProcessEvents` from inside `R_CheckUserInterrupt`
//! rather than by a signal handler, and polling the flag alone would miss it.
//! That branch therefore keeps the exposure: Windows profiles on a thread that
//! suspends the R thread and walks its contexts too, so the race is narrowed
//! to one platform rather than removed.
//!
//! One thing the Unix path gives up: `R_CheckUserInterrupt` also ran
//! `R_ProcessEvents`, which services polled GUI events and `setTimeLimit()`
//! checks. Neither matters to a solve, but a time limit set around
//! `balance()` now fires when the call returns rather than at the next
//! iteration.

#[cfg(windows)]
use std::ffi::c_void;

// R names these in its own style, which is not Rust's convention for statics.
#[cfg(not(windows))]
#[allow(non_upper_case_globals)]
unsafe extern "C" {
    /// Nonzero once R's `SIGINT` handler has recorded an interrupt, until R
    /// delivers it.
    static mut R_interrupts_pending: std::ffi::c_int;
    /// Nonzero while R has interrupt delivery suspended.
    static mut R_interrupts_suspended: std::ffi::c_int;
}

#[cfg(windows)]
unsafe extern "C" {
    fn R_CheckUserInterrupt();
    fn R_ToplevelExec(fun: extern "C" fn(*mut c_void), data: *mut c_void) -> i32;
}

#[cfg(windows)]
extern "C" fn check(_data: *mut c_void) {
    unsafe { R_CheckUserInterrupt() };
}

/// Return `true` when a user interrupt is pending.
#[cfg(not(windows))]
pub fn pending() -> bool {
    // Volatile because both flags change outside anything the compiler can
    // see, and a poll sits in a loop it would otherwise be hoisted out of.
    unsafe {
        std::ptr::read_volatile(&raw const R_interrupts_suspended) == 0
            && std::ptr::read_volatile(&raw const R_interrupts_pending) != 0
    }
}

/// Return `true` when a user interrupt is pending.
#[cfg(windows)]
pub fn pending() -> bool {
    // `R_ToplevelExec` returns TRUE (nonzero) when `check` completed and FALSE
    // (zero) when it longjmped, which happens exactly when an interrupt fired.
    unsafe { R_ToplevelExec(check, std::ptr::null_mut()) == 0 }
}
