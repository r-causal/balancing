//! User-interrupt polling that never unwinds through Rust.
//!
//! `R_CheckUserInterrupt` longjmps when an interrupt is pending, which would
//! skip Rust destructors if it fired inside the solver. Running it under
//! `R_ToplevelExec` contains the longjmp: the wrapper returns `FALSE` when the
//! check jumped, which the solver reads as a request to stop and then unwinds
//! normally. The check runs only on the main thread, between solver iterations,
//! never inside a parallel region.

use std::ffi::c_void;

unsafe extern "C" {
    fn R_CheckUserInterrupt();
    fn R_ToplevelExec(fun: extern "C" fn(*mut c_void), data: *mut c_void) -> i32;
}

extern "C" fn check(_data: *mut c_void) {
    unsafe { R_CheckUserInterrupt() };
}

/// Return `true` when a user interrupt is pending.
pub fn pending() -> bool {
    // `R_ToplevelExec` returns TRUE (nonzero) when `check` completed and FALSE
    // (zero) when it longjmped, which happens exactly when an interrupt fired.
    unsafe { R_ToplevelExec(check, std::ptr::null_mut()) == 0 }
}
