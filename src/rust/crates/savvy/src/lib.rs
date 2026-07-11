//! savvy exports for the balancing package.
//!
//! Every function here is a thin adapter over [`balancing_core`]: it converts
//! R data to plain Rust, calls into the core, and packs the result back into an
//! SEXP. No numerical work happens in this crate.

use savvy::{OwnedIntegerSexp, OwnedListSexp, OwnedStringSexp, savvy};

/// Report the parallel resources the Rust core observes.
///
/// @returns A list with two elements: `available`, the integer thread count, and
///   `cap_source`, a string naming the constraint that set it (`"system"` or
///   `"OMP_THREAD_LIMIT"`).
/// @export
#[savvy]
fn thread_info() -> savvy::Result<savvy::Sexp> {
    let (available, cap_source) = balancing_core::available_threads();

    let mut available_sexp = OwnedIntegerSexp::new(1)?;
    available_sexp[0] = i32::try_from(available).unwrap_or(i32::MAX);

    let mut cap_source_sexp = OwnedStringSexp::new(1)?;
    cap_source_sexp.set_elt(0, cap_source)?;

    let mut out = OwnedListSexp::new(2, true)?;
    out.set_name_and_value(0, "available", available_sexp)?;
    out.set_name_and_value(1, "cap_source", cap_source_sexp)?;

    Ok(out.into())
}

#[cfg(test)]
mod tests {
    #[test]
    fn core_reports_a_positive_thread_count() {
        let (available, _cap_source) = balancing_core::available_threads();
        assert!(available >= 1);
    }
}
