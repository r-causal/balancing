//! savvy exports for the balancing package.
//!
//! Every function here is a thin adapter over [`balancing_core`]: it converts
//! R data to plain Rust, calls into the core, and packs the result back into an
//! SEXP. No numerical work happens in this crate.

mod convert;
mod interrupt;

use balancing_core::methods::entropy::{
    EntropyInputs, EntropyResult, scale_estimating_output, solve_continuous, solve_discrete,
};
use savvy::{
    IntegerSexp, ListSexp, NullSexp, OwnedIntegerSexp, OwnedListSexp, OwnedLogicalSexp,
    OwnedRealSexp, OwnedStringSexp, RealSexp, savvy,
};

use convert::{parse_entropy_options, real_matrix, real_vector};

/// Report the parallel resources the Rust core observes.
///
/// A list with two elements: `available`, the integer thread count, and
/// `cap_source`, a string naming the constraint that set it (`"system"` or
/// `"OMP_THREAD_LIMIT"`).
///
/// This is an internal solver entry point, called from the R layer rather than
/// by users, so it is not exported. `@noRd` keeps it out of the reference and
/// out of NAMESPACE; savvy copies these doc lines into the generated wrapper, so
/// the tag survives wrapper regeneration.
/// @noRd
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

fn scalar_logical(value: bool) -> savvy::Result<OwnedLogicalSexp> {
    let mut out = OwnedLogicalSexp::new(1)?;
    out.set_elt(0, value)?;
    Ok(out)
}

fn scalar_integer(value: usize) -> savvy::Result<OwnedIntegerSexp> {
    let mut out = OwnedIntegerSexp::new(1)?;
    out[0] = i32::try_from(value).unwrap_or(i32::MAX);
    Ok(out)
}

fn scalar_real(value: f64) -> savvy::Result<OwnedRealSexp> {
    let mut out = OwnedRealSexp::new(1)?;
    out[0] = value;
    Ok(out)
}

fn scalar_string(value: &str) -> savvy::Result<OwnedStringSexp> {
    let mut out = OwnedStringSexp::new(1)?;
    out.set_elt(0, value)?;
    Ok(out)
}

/// Assemble the leading fields common to both entropy result lists into `out`:
/// `weights`, `duals`, `converged`, `iterations`, and `grad_norm` at indices
/// 0 through 4.
fn set_common_fields(out: &mut OwnedListSexp, result: &EntropyResult) -> savvy::Result<()> {
    out.set_name_and_value(0, "weights", real_vector(&result.weights)?)?;
    out.set_name_and_value(1, "duals", real_vector(&result.duals)?)?;
    out.set_name_and_value(2, "converged", scalar_logical(result.converged)?)?;
    out.set_name_and_value(3, "iterations", scalar_integer(result.iterations)?)?;
    out.set_name_and_value(4, "grad_norm", scalar_real(result.grad_norm)?)?;
    Ok(())
}

/// Set an optional matrix field, writing R `NULL` when the value is absent.
fn set_optional_matrix(
    out: &mut OwnedListSexp,
    index: usize,
    name: &str,
    value: &Option<Vec<f64>>,
    nrow: usize,
    ncol: usize,
) -> savvy::Result<()> {
    match value {
        Some(data) => out.set_name_and_value(index, name, real_matrix(data, nrow, ncol)?),
        None => out.set_name_and_value(index, name, NullSexp),
    }
}

/// Solve a discrete (binary or categorical) entropy balancing problem.
///
/// Internal solver entry point, called from the R layer rather than by users, so
/// it is not exported. `@noRd` keeps it out of the reference and out of
/// NAMESPACE, and survives wrapper regeneration because savvy copies these doc
/// lines into the generated wrapper.
/// @noRd
// The argument list is the fixed savvy boundary signature; the balancing inputs
// are irreducibly numerous.
#[allow(clippy::too_many_arguments)]
#[savvy]
fn solve_entropy(
    covs: RealSexp,
    group_idx: IntegerSexp,
    targets: RealSexp,
    base_weights: RealSexp,
    s_weights: RealSexp,
    tols: RealSexp,
    n_eff: f64,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let p = targets.len();
    let opts = parse_entropy_options(options)?;

    if covs.len() != n * p {
        return Err(savvy::Error::new(format!(
            "covs has {} elements but n * p = {n} * {p} = {}",
            covs.len(),
            n * p
        )));
    }
    if base_weights.len() != n || group_idx.len() != n {
        return Err(savvy::Error::new(
            "base_weights and group_idx must have length n",
        ));
    }
    if tols.len() != p {
        return Err(savvy::Error::new("tols must have length p"));
    }

    let inputs = EntropyInputs {
        covs: covs.as_slice(),
        n,
        p,
        targets: targets.as_slice(),
        tols: tols.as_slice(),
        base: base_weights.as_slice(),
        s: s_weights.as_slice(),
        n_eff,
        threads: opts.threads,
        max_iter: opts.max_iter,
        tol: opts.tol,
        solver: opts.solver,
    };

    let mut result = solve_discrete(&inputs, group_idx.as_slice(), &interrupt::pending);
    if let Some(scale) = &opts.esteq_scale {
        scale_estimating_output(&mut result, p, scale).map_err(savvy::Error::new)?;
    }
    let total_params = result.duals.len();

    let mut out = OwnedListSexp::new(10, true)?;
    set_common_fields(&mut out, &result)?;
    out.set_name_and_value(5, "solver", scalar_string(result.solver)?)?;
    set_optional_matrix(&mut out, 6, "psi", &result.psi, n, total_params)?;
    set_optional_matrix(&mut out, 7, "jac", &result.jac, total_params, total_params)?;
    set_optional_matrix(&mut out, 8, "dw_dbeta", &result.dw_dbeta, n, total_params)?;
    out.set_name_and_value(9, "interrupted", scalar_logical(result.interrupted)?)?;
    Ok(out.into())
}

/// Solve a continuous-exposure entropy balancing problem over the whole sample.
///
/// Internal solver entry point, called from the R layer rather than by users, so
/// it is not exported. `@noRd` keeps it out of the reference and out of
/// NAMESPACE, and survives wrapper regeneration because savvy copies these doc
/// lines into the generated wrapper.
/// @noRd
// The argument list is the fixed savvy boundary signature; the balancing inputs
// are irreducibly numerous.
#[allow(clippy::too_many_arguments)]
#[savvy]
fn solve_entropy_cont(
    covs: RealSexp,
    targets: RealSexp,
    tols: RealSexp,
    dist_ind: IntegerSexp,
    base_weights: RealSexp,
    s_weights: RealSexp,
    n_eff: f64,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let p = targets.len();
    let opts = parse_entropy_options(options)?;

    if covs.len() != n * p {
        return Err(savvy::Error::new(format!(
            "covs has {} elements but n * p = {n} * {p} = {}",
            covs.len(),
            n * p
        )));
    }
    if base_weights.len() != n {
        return Err(savvy::Error::new("base_weights must have length n"));
    }
    if tols.len() != p || dist_ind.len() != p {
        return Err(savvy::Error::new("tols and dist_ind must have length p"));
    }

    let inputs = EntropyInputs {
        covs: covs.as_slice(),
        n,
        p,
        targets: targets.as_slice(),
        tols: tols.as_slice(),
        base: base_weights.as_slice(),
        s: s_weights.as_slice(),
        n_eff,
        threads: opts.threads,
        max_iter: opts.max_iter,
        tol: opts.tol,
        solver: opts.solver,
    };

    let mut result = solve_continuous(&inputs, dist_ind.as_slice(), &interrupt::pending);
    if let Some(scale) = &opts.esteq_scale {
        scale_estimating_output(&mut result, p, scale).map_err(savvy::Error::new)?;
    }
    let total_params = result.duals.len();

    // The continuous list carries dw_dbeta alongside psi and jac for the
    // estimating-equations container, plus the solver that actually ran and the
    // interrupt flag.
    let mut out = OwnedListSexp::new(10, true)?;
    set_common_fields(&mut out, &result)?;
    set_optional_matrix(&mut out, 5, "psi", &result.psi, n, total_params)?;
    set_optional_matrix(&mut out, 6, "jac", &result.jac, total_params, total_params)?;
    set_optional_matrix(&mut out, 7, "dw_dbeta", &result.dw_dbeta, n, total_params)?;
    out.set_name_and_value(8, "solver", scalar_string(result.solver)?)?;
    out.set_name_and_value(9, "interrupted", scalar_logical(result.interrupted)?)?;
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
