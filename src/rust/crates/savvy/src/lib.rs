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
use balancing_core::methods::ipt::{self, IptInputs, IptResult};
use savvy::{
    IntegerSexp, ListSexp, NullSexp, OwnedIntegerSexp, OwnedListSexp, OwnedLogicalSexp,
    OwnedRealSexp, OwnedStringSexp, RealSexp, savvy,
};

use convert::{
    parse_binary_estimand, parse_entropy_options, parse_ipt_options, parse_link,
    parse_multi_estimand, real_matrix, real_vector,
};

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

/// Pack an inverse probability tilting result into its R list.
///
/// The list is `weights`, `ps`, `coefs`, `converged`, `iterations`,
/// `grad_norm`, `psi`, `jac`, `dw_dbeta`, and the interrupt flag, in that order.
fn ipt_result_list(result: &IptResult, n: usize) -> savvy::Result<savvy::Sexp> {
    let total_params = result.coefs.len();
    let mut out = OwnedListSexp::new(10, true)?;
    out.set_name_and_value(0, "weights", real_vector(&result.weights)?)?;
    out.set_name_and_value(1, "ps", real_vector(&result.ps)?)?;
    out.set_name_and_value(2, "coefs", real_vector(&result.coefs)?)?;
    out.set_name_and_value(3, "converged", scalar_logical(result.converged)?)?;
    out.set_name_and_value(4, "iterations", scalar_integer(result.iterations)?)?;
    out.set_name_and_value(5, "grad_norm", scalar_real(result.grad_norm)?)?;
    out.set_name_and_value(6, "psi", real_matrix(&result.psi, n, total_params)?)?;
    out.set_name_and_value(
        7,
        "jac",
        real_matrix(&result.jac, total_params, total_params)?,
    )?;
    out.set_name_and_value(
        8,
        "dw_dbeta",
        real_matrix(&result.dw_dbeta, n, total_params)?,
    )?;
    out.set_name_and_value(9, "interrupted", scalar_logical(result.interrupted)?)?;
    Ok(out.into())
}

/// Solve a binary inverse probability tilting problem.
///
/// `treat` holds the zero/one treatment indicator; `estimand` is one of `ate`,
/// `att`, or `atc`, and `link` is the propensity link.
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
fn solve_ipt(
    covs: RealSexp,
    treat: IntegerSexp,
    s_weights: RealSexp,
    estimand: &str,
    link: &str,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let p = if n == 0 { 0 } else { covs.len() / n };
    let opts = parse_ipt_options(options)?;
    let link = parse_link(link)?;
    let estimand = parse_binary_estimand(estimand)?;

    if n == 0 || covs.len() != n * p {
        return Err(savvy::Error::new(format!(
            "covs has {} elements but is not a multiple of n = {n}",
            covs.len()
        )));
    }
    if treat.len() != n {
        return Err(savvy::Error::new("treat must have length n"));
    }

    let inputs = IptInputs {
        covs: covs.as_slice(),
        n,
        p,
        treat: treat.as_slice(),
        n_levels: 2,
        s: s_weights.as_slice(),
        link,
        estimand,
        threads: opts.threads,
        max_iter: opts.max_iter,
        tol: opts.tol,
    };
    let result = ipt::solve(&inputs, &interrupt::pending);
    ipt_result_list(&result, n)
}

/// Solve a categorical inverse probability tilting problem.
///
/// `treat_idx` holds the zero-based level of each unit; `focal` is the focal
/// level index used by `att` and `atc` and ignored by `ate`.
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
fn solve_ipt_multi(
    covs: RealSexp,
    treat_idx: IntegerSexp,
    focal: i32,
    s_weights: RealSexp,
    estimand: &str,
    link: &str,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let p = if n == 0 { 0 } else { covs.len() / n };
    let opts = parse_ipt_options(options)?;
    let link = parse_link(link)?;

    if n == 0 || covs.len() != n * p {
        return Err(savvy::Error::new(format!(
            "covs has {} elements but is not a multiple of n = {n}",
            covs.len()
        )));
    }
    if treat_idx.len() != n {
        return Err(savvy::Error::new("treat_idx must have length n"));
    }
    let treat = treat_idx.as_slice();
    if treat.iter().any(|&t| t < 0) {
        return Err(savvy::Error::new("treat_idx values must be non-negative"));
    }
    let n_levels = treat.iter().copied().max().map_or(0, |m| m + 1) as usize;
    if focal < 0 || focal as usize >= n_levels {
        return Err(savvy::Error::new(
            "focal must be a level present in treat_idx",
        ));
    }
    let estimand = parse_multi_estimand(estimand, focal as usize)?;

    let inputs = IptInputs {
        covs: covs.as_slice(),
        n,
        p,
        treat,
        n_levels,
        s: s_weights.as_slice(),
        link,
        estimand,
        threads: opts.threads,
        max_iter: opts.max_iter,
        tol: opts.tol,
    };
    let result = ipt::solve(&inputs, &interrupt::pending);
    ipt_result_list(&result, n)
}

#[cfg(test)]
mod tests {
    #[test]
    fn core_reports_a_positive_thread_count() {
        let (available, _cap_source) = balancing_core::available_threads();
        assert!(available >= 1);
    }
}
