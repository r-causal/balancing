//! savvy exports for the balancing package.
//!
//! Every function here is a thin adapter over [`balancing_core`]: it converts
//! R data to plain Rust, calls into the core, and packs the result back into an
//! SEXP. No numerical work happens in this crate.

mod convert;
mod interrupt;

use balancing_core::dist::kernels::{Kernel, KernelParams, build_kernel};
use balancing_core::methods::cbps::{
    self, CbpsContInputs, CbpsInputs, CbpsMultiInputs, CbpsResult,
};
use balancing_core::methods::cfd::{self, CfdDiscreteInputs, CfdEstimand, CfdResult};
use balancing_core::methods::energy::{
    EnergyContInputs, EnergyDiscreteInputs, EnergyEstimand, EnergyResult,
};
use balancing_core::methods::entropy::{
    EntropyInputs, EntropyResult, scale_estimating_output, solve_continuous, solve_discrete,
};
use balancing_core::methods::ipt::{self, IptEstimand, IptInputs, IptResult};
use balancing_core::methods::sbw::{
    self, SbwContInputs, SbwDiscreteInputs, SbwEstimand, SbwResult,
};
use savvy::{
    IntegerSexp, ListSexp, LogicalSexp, NullSexp, OwnedIntegerSexp, OwnedListSexp,
    OwnedLogicalSexp, OwnedRealSexp, OwnedStringSexp, RealSexp, savvy,
};

use convert::{
    parse_binary_estimand, parse_cbps_estimand, parse_cbps_multi_estimand, parse_cfd_options,
    parse_distance, parse_entropy_options, parse_ipt_options, parse_kernel, parse_kernel_options,
    parse_link, parse_multi_estimand, parse_qp_options, parse_sbw_norm, parse_sbw_options,
    parse_smoothness, real_matrix, real_vector, require_binary_treat, require_dense_levels,
    require_finite,
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

/// Re-evaluate the discrete entropy estimating functions at a set of duals.
///
/// Given the solved duals in `coefs` (`p` per group, stacked in group order) and
/// the original solve inputs, returns the `n` by `P` per-unit estimating
/// functions at those parameters, with the same per-group renormalization
/// `esteq_scale` the solve output carried. It supports a finite-difference check
/// of the stored Jacobian without reimplementing the tilt math in R.
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
fn eval_psi_entropy(
    coefs: RealSexp,
    covs: RealSexp,
    group_idx: IntegerSexp,
    targets: RealSexp,
    base_weights: RealSexp,
    s_weights: RealSexp,
    n_eff: f64,
    esteq_scale: RealSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let p = targets.len();

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
    let n_groups = esteq_scale.len();
    if p == 0 || coefs.len() != n_groups * p {
        return Err(savvy::Error::new(
            "coefs must have length p times the number of groups",
        ));
    }

    let inputs = EntropyInputs {
        covs: covs.as_slice(),
        n,
        p,
        targets: targets.as_slice(),
        tols: &vec![0.0; p],
        base: base_weights.as_slice(),
        s: s_weights.as_slice(),
        n_eff,
        threads: 1,
        max_iter: 0,
        tol: 0.0,
        solver: balancing_core::methods::entropy::EntropySolver::Newton,
    };
    let psi = balancing_core::methods::entropy::eval_psi_discrete(
        &inputs,
        group_idx.as_slice(),
        coefs.as_slice(),
        esteq_scale.as_slice(),
    )
    .map_err(savvy::Error::new)?;
    Ok(real_matrix(&psi, n, n_groups * p)?.into())
}

/// Re-evaluate the discrete entropy balancing weights at a set of duals.
///
/// Given the solved duals in `coefs` (`p` per group, stacked in group order) and
/// the original solve inputs, returns the length-`n` weight vector at those
/// parameters, with the same per-group renormalization `esteq_scale` the solve
/// output carried. Units the solve leaves out, marked by a negative
/// `group_idx`, carry zero. It supports a sandwich variance that treats the
/// weights as a function of the duals without reimplementing the tilt math in R.
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
fn eval_weights_entropy(
    coefs: RealSexp,
    covs: RealSexp,
    group_idx: IntegerSexp,
    targets: RealSexp,
    base_weights: RealSexp,
    s_weights: RealSexp,
    n_eff: f64,
    esteq_scale: RealSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let p = targets.len();

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
    let n_groups = esteq_scale.len();
    if p == 0 || coefs.len() != n_groups * p {
        return Err(savvy::Error::new(
            "coefs must have length p times the number of groups",
        ));
    }

    let inputs = EntropyInputs {
        covs: covs.as_slice(),
        n,
        p,
        targets: targets.as_slice(),
        tols: &vec![0.0; p],
        base: base_weights.as_slice(),
        s: s_weights.as_slice(),
        n_eff,
        threads: 1,
        max_iter: 0,
        tol: 0.0,
        solver: balancing_core::methods::entropy::EntropySolver::Newton,
    };
    let weights = balancing_core::methods::entropy::eval_weights_discrete(
        &inputs,
        group_idx.as_slice(),
        coefs.as_slice(),
        esteq_scale.as_slice(),
    )
    .map_err(savvy::Error::new)?;
    Ok(real_vector(&weights)?.into())
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
    require_finite(covs.as_slice(), "covs")?;
    require_finite(s_weights.as_slice(), "s_weights")?;
    require_binary_treat(treat.as_slice())?;

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
    require_finite(covs.as_slice(), "covs")?;
    require_finite(s_weights.as_slice(), "s_weights")?;
    require_dense_levels(treat, n_levels)?;
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

/// Re-evaluate the inverse probability tilting estimating functions at a set of
/// coefficients.
///
/// Given the solved coefficients in `coefs` (`p` per block, stacked in block
/// order) and the original solve inputs, returns the `n` by `P` per-unit
/// estimating functions at those parameters. `treat_idx` holds the zero-based
/// level of each unit and `focal` the focal level index the focal estimands use.
/// The binary and categorical fits share this entrypoint. It supports a
/// finite-difference check of the stored Jacobian without reimplementing the
/// tilt math in R.
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
fn eval_psi_ipt(
    coefs: RealSexp,
    covs: RealSexp,
    treat_idx: IntegerSexp,
    focal: i32,
    s_weights: RealSexp,
    estimand: &str,
    link: &str,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let p = if n == 0 { 0 } else { covs.len() / n };
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
    require_finite(covs.as_slice(), "covs")?;
    require_finite(s_weights.as_slice(), "s_weights")?;
    require_dense_levels(treat, n_levels)?;
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
        threads: 1,
        max_iter: 0,
        tol: 0.0,
    };
    let n_blocks = match estimand {
        IptEstimand::Ate => n_levels,
        IptEstimand::Focal(_) => n_levels - 1,
    };
    let psi = ipt::eval_psi(&inputs, coefs.as_slice()).map_err(savvy::Error::new)?;
    Ok(real_matrix(&psi, n, n_blocks * p)?.into())
}

/// Re-evaluate the inverse probability tilting weights at a set of
/// coefficients.
///
/// Given the solved coefficients in `coefs` (`p` per block, stacked in block
/// order) and the original solve inputs, returns the length-`n` weight vector at
/// those parameters. `treat_idx` holds the zero-based level of each unit and
/// `focal` the focal level index the focal estimands use, whose units carry
/// weight one. The binary and categorical fits share this entrypoint. It
/// supports a sandwich variance that treats the weights as a function of the
/// coefficients without reimplementing the tilt math in R.
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
fn eval_weights_ipt(
    coefs: RealSexp,
    covs: RealSexp,
    treat_idx: IntegerSexp,
    focal: i32,
    s_weights: RealSexp,
    estimand: &str,
    link: &str,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let p = if n == 0 { 0 } else { covs.len() / n };
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
    require_finite(covs.as_slice(), "covs")?;
    require_finite(s_weights.as_slice(), "s_weights")?;
    require_dense_levels(treat, n_levels)?;
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
        threads: 1,
        max_iter: 0,
        tol: 0.0,
    };
    let weights = ipt::eval_weights(&inputs, coefs.as_slice()).map_err(savvy::Error::new)?;
    Ok(real_vector(&weights)?.into())
}

/// Pack a covariate balancing propensity score result into its R list.
///
/// The list is `weights`, `ps`, `coefs`, `converged`, `iterations`,
/// `obj_value`, then the estimating-equations matrices `psi`, `jac`, and
/// `dw_dbeta` (each `NULL` for the over-identified and continuous forms), the
/// `gmm_obj` criterion (`NULL` except for the over-identified form), and the
/// interrupt flag, in that order.
fn cbps_result_list(result: &CbpsResult, n: usize) -> savvy::Result<savvy::Sexp> {
    let total_params = result.coefs.len();
    let mut out = OwnedListSexp::new(11, true)?;
    out.set_name_and_value(0, "weights", real_vector(&result.weights)?)?;
    out.set_name_and_value(1, "ps", real_vector(&result.ps)?)?;
    out.set_name_and_value(2, "coefs", real_vector(&result.coefs)?)?;
    out.set_name_and_value(3, "converged", scalar_logical(result.converged)?)?;
    out.set_name_and_value(4, "iterations", scalar_integer(result.iterations)?)?;
    out.set_name_and_value(5, "obj_value", scalar_real(result.obj_value)?)?;
    set_optional_matrix(&mut out, 6, "psi", &result.psi, n, total_params)?;
    set_optional_matrix(&mut out, 7, "jac", &result.jac, total_params, total_params)?;
    set_optional_matrix(&mut out, 8, "dw_dbeta", &result.dw_dbeta, n, total_params)?;
    match result.gmm_obj {
        Some(value) => out.set_name_and_value(9, "gmm_obj", scalar_real(value)?)?,
        None => out.set_name_and_value(9, "gmm_obj", NullSexp)?,
    }
    out.set_name_and_value(10, "interrupted", scalar_logical(result.interrupted)?)?;
    Ok(out.into())
}

/// Solve a binary covariate balancing propensity score problem.
///
/// `covs_mod` is the propensity-model design and `covs_bal` the balance design;
/// for the just-identified form they must have the same shape. `treat` holds the
/// zero/one treatment indicator; `estimand` is one of `ate`, `att`, `atc`, or
/// `ato`; `link` is the propensity link. `over` selects the over-identified GMM
/// criterion, and `twostep` its two-step weighting matrix.
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
fn solve_cbps(
    covs_mod: RealSexp,
    covs_bal: RealSexp,
    treat: IntegerSexp,
    s_weights: RealSexp,
    estimand: &str,
    link: &str,
    over: bool,
    twostep: bool,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let p_mod = if n == 0 { 0 } else { covs_mod.len() / n };
    let p_bal = if n == 0 { 0 } else { covs_bal.len() / n };
    let opts = parse_ipt_options(options)?;
    let link = parse_link(link)?;
    let estimand = parse_cbps_estimand(estimand)?;

    if n == 0 || covs_mod.len() != n * p_mod {
        return Err(savvy::Error::new(format!(
            "covs_mod has {} elements but is not a multiple of n = {n}",
            covs_mod.len()
        )));
    }
    if covs_bal.len() != n * p_bal {
        return Err(savvy::Error::new(format!(
            "covs_bal has {} elements but is not a multiple of n = {n}",
            covs_bal.len()
        )));
    }
    if treat.len() != n {
        return Err(savvy::Error::new("treat must have length n"));
    }
    if !over && p_mod != p_bal {
        return Err(savvy::Error::new(
            "the just-identified fit requires covs_mod and covs_bal to have the same number of columns",
        ));
    }

    let inputs = CbpsInputs {
        covs_mod: covs_mod.as_slice(),
        covs_bal: covs_bal.as_slice(),
        n,
        p_mod,
        p_bal,
        treat: treat.as_slice(),
        s: s_weights.as_slice(),
        link,
        estimand,
        over,
        twostep,
        threads: opts.threads,
        max_iter: opts.max_iter,
        tol: opts.tol,
    };
    let result = cbps::solve(&inputs, &interrupt::pending);
    cbps_result_list(&result, n)
}

/// Re-evaluate the binary just-identified covariate balancing propensity score
/// estimating functions at a set of coefficients.
///
/// Given the solved coefficients in `coefs` and the original solve inputs,
/// returns the `n` by `p` per-unit estimating functions at those parameters:
/// column `j` is `s_i c_i(beta) x_ij`, the balancing factor times the model
/// covariate. It supports a finite-difference check of the stored Jacobian
/// without reimplementing the balancing-factor math in R.
///
/// Internal solver entry point, called from the R layer rather than by users, so
/// it is not exported. `@noRd` keeps it out of the reference and out of
/// NAMESPACE, and survives wrapper regeneration because savvy copies these doc
/// lines into the generated wrapper.
/// @noRd
#[savvy]
fn eval_psi_cbps(
    coefs: RealSexp,
    covs: RealSexp,
    treat: IntegerSexp,
    s_weights: RealSexp,
    estimand: &str,
    link: &str,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let p = if n == 0 { 0 } else { covs.len() / n };
    let link = parse_link(link)?;
    let estimand = parse_cbps_estimand(estimand)?;

    if n == 0 || covs.len() != n * p {
        return Err(savvy::Error::new(format!(
            "covs has {} elements but is not a multiple of n = {n}",
            covs.len()
        )));
    }
    if treat.len() != n {
        return Err(savvy::Error::new("treat must have length n"));
    }
    if coefs.len() != p {
        return Err(savvy::Error::new("coefs must have length p"));
    }

    let inputs = CbpsInputs {
        covs_mod: covs.as_slice(),
        covs_bal: covs.as_slice(),
        n,
        p_mod: p,
        p_bal: p,
        treat: treat.as_slice(),
        s: s_weights.as_slice(),
        link,
        estimand,
        over: false,
        twostep: false,
        threads: 1,
        max_iter: 0,
        tol: 0.0,
    };
    let psi = cbps::eval_psi_binary_just(&inputs, coefs.as_slice());
    Ok(real_matrix(&psi, n, p)?.into())
}

/// Re-evaluate the binary just-identified covariate balancing propensity score
/// weights at a set of coefficients.
///
/// Given the solved coefficients in `coefs` and the original solve inputs,
/// returns the length-`n` weight vector at those parameters: the estimand's
/// weight function evaluated at the modeled propensity and the unit's treatment
/// indicator. It supports a sandwich variance that treats the weights as a
/// function of the coefficients without reimplementing the propensity math in R.
///
/// Internal solver entry point, called from the R layer rather than by users, so
/// it is not exported. `@noRd` keeps it out of the reference and out of
/// NAMESPACE, and survives wrapper regeneration because savvy copies these doc
/// lines into the generated wrapper.
/// @noRd
#[savvy]
fn eval_weights_cbps(
    coefs: RealSexp,
    covs: RealSexp,
    treat: IntegerSexp,
    s_weights: RealSexp,
    estimand: &str,
    link: &str,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let p = if n == 0 { 0 } else { covs.len() / n };
    let link = parse_link(link)?;
    let estimand = parse_cbps_estimand(estimand)?;

    if n == 0 || covs.len() != n * p {
        return Err(savvy::Error::new(format!(
            "covs has {} elements but is not a multiple of n = {n}",
            covs.len()
        )));
    }
    if treat.len() != n {
        return Err(savvy::Error::new("treat must have length n"));
    }
    if coefs.len() != p {
        return Err(savvy::Error::new("coefs must have length p"));
    }

    let inputs = CbpsInputs {
        covs_mod: covs.as_slice(),
        covs_bal: covs.as_slice(),
        n,
        p_mod: p,
        p_bal: p,
        treat: treat.as_slice(),
        s: s_weights.as_slice(),
        link,
        estimand,
        over: false,
        twostep: false,
        threads: 1,
        max_iter: 0,
        tol: 0.0,
    };
    let weights = cbps::eval_weights_binary_just(&inputs, coefs.as_slice());
    Ok(real_vector(&weights)?.into())
}

/// Solve a categorical covariate balancing propensity score problem.
///
/// `treat_idx` holds the zero-based level of each unit; `focal` is the focal
/// level index used by `att` and ignored by `ate`. The categorical form is
/// always just-identified.
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
fn solve_cbps_multi(
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
    let estimand = parse_cbps_multi_estimand(estimand)?;

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

    let inputs = CbpsMultiInputs {
        covs: covs.as_slice(),
        n,
        p,
        treat,
        n_levels,
        focal: focal as usize,
        s: s_weights.as_slice(),
        link,
        estimand,
        threads: opts.threads,
        max_iter: opts.max_iter,
        tol: opts.tol,
    };
    let result = cbps::solve_multi(&inputs, &interrupt::pending);
    cbps_result_list(&result, n)
}

/// Solve a continuous-exposure covariate balancing propensity score problem.
///
/// `expo` holds the continuous exposure of each unit. The continuous form
/// targets the average treatment effect and supplies no estimating equations.
///
/// Internal solver entry point, called from the R layer rather than by users, so
/// it is not exported. `@noRd` keeps it out of the reference and out of
/// NAMESPACE, and survives wrapper regeneration because savvy copies these doc
/// lines into the generated wrapper.
/// @noRd
#[savvy]
fn solve_cbps_cont(
    covs: RealSexp,
    expo: RealSexp,
    s_weights: RealSexp,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let p = if n == 0 { 0 } else { covs.len() / n };
    let opts = parse_ipt_options(options)?;

    if n == 0 || covs.len() != n * p {
        return Err(savvy::Error::new(format!(
            "covs has {} elements but is not a multiple of n = {n}",
            covs.len()
        )));
    }
    if expo.len() != n {
        return Err(savvy::Error::new("expo must have length n"));
    }

    let inputs = CbpsContInputs {
        covs: covs.as_slice(),
        n,
        p,
        expo: expo.as_slice(),
        s: s_weights.as_slice(),
        threads: opts.threads,
        max_iter: opts.max_iter,
        tol: opts.tol,
    };
    let result = cbps::solve_cont(&inputs, &interrupt::pending);
    cbps_result_list(&result, n)
}

/// Pack an energy balancing result into its R list.
///
/// The list is `weights`, `duals`, `converged`, `iterations`, `objective`,
/// `solver_status` (the backend identity), `status` (the solver's terminal
/// status name), `pri_res`, `dua_res`, and the interrupt flag, in that order. The
/// quadratic-program family has no estimating equations, so none appear.
fn energy_result_list(result: &EnergyResult) -> savvy::Result<savvy::Sexp> {
    let mut out = OwnedListSexp::new(10, true)?;
    out.set_name_and_value(0, "weights", real_vector(&result.weights)?)?;
    out.set_name_and_value(1, "duals", real_vector(&result.duals)?)?;
    out.set_name_and_value(2, "converged", scalar_logical(result.converged)?)?;
    out.set_name_and_value(3, "iterations", scalar_integer(result.iterations)?)?;
    out.set_name_and_value(4, "objective", scalar_real(result.objective)?)?;
    out.set_name_and_value(5, "solver_status", scalar_string(result.backend)?)?;
    out.set_name_and_value(6, "status", scalar_string(result.status)?)?;
    out.set_name_and_value(7, "pri_res", scalar_real(result.pri_res)?)?;
    out.set_name_and_value(8, "dua_res", scalar_real(result.dua_res)?)?;
    out.set_name_and_value(9, "interrupted", scalar_logical(result.interrupted)?)?;
    Ok(out.into())
}

/// Validate the shared shape of the energy inputs, returning the covariate column
/// count `p` and the moment-constraint column count `q`.
fn energy_dims(
    n: usize,
    covs: &RealSexp,
    s_weights: &RealSexp,
    moment_covs: &RealSexp,
    targets: &RealSexp,
    tols: &RealSexp,
) -> savvy::Result<(usize, usize)> {
    if n == 0 || covs.len() % n != 0 {
        return Err(savvy::Error::new(format!(
            "covs has {} elements but is not a multiple of n = {n}",
            covs.len()
        )));
    }
    let p = covs.len() / n;
    if s_weights.len() != n {
        return Err(savvy::Error::new("s_weights must have length n"));
    }
    let q = targets.len();
    if tols.len() != q || moment_covs.len() != n * q {
        return Err(savvy::Error::new(
            "moment_covs must be n by length(targets), and tols must match targets",
        ));
    }
    Ok((p, q))
}

/// Solve a binary-exposure energy balancing problem.
///
/// `treat` holds the zero/one exposure indicator; `distance` names the covariate
/// distance definition; `estimand` is one of `ate`, `att`, or `atc`; `improved`
/// selects the between-group term of the improved average-treatment-effect
/// variant. `moment_covs` are the standardized moment-constraint columns with
/// `targets` and `tols`, both empty when no moment constraints are requested.
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
fn solve_energy(
    covs: RealSexp,
    treat: IntegerSexp,
    s_weights: RealSexp,
    distance: &str,
    estimand: &str,
    improved: bool,
    moment_covs: RealSexp,
    targets: RealSexp,
    tols: RealSexp,
    min_weight: f64,
    weight_penalty: f64,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let (p, q) = energy_dims(n, &covs, &s_weights, &moment_covs, &targets, &tols)?;
    if treat.len() != n {
        return Err(savvy::Error::new("treat must have length n"));
    }
    let treat_slice = treat.as_slice();
    if treat_slice.iter().any(|&t| t != 0 && t != 1) {
        return Err(savvy::Error::new("treat must hold only zero and one"));
    }
    let opts = parse_qp_options(options)?;
    let dist = parse_distance(distance)?;
    let energy_estimand = match estimand {
        "ate" => EnergyEstimand::Ate { improved },
        "att" => EnergyEstimand::Focal { focal: 1 },
        "atc" => EnergyEstimand::Focal { focal: 0 },
        other => {
            return Err(savvy::Error::new(format!(
                "unknown estimand `{other}`; expected ate, att, or atc"
            )));
        }
    };

    let inputs = EnergyDiscreteInputs {
        covs: covs.as_slice(),
        n,
        p,
        distance: dist,
        levels: treat_slice,
        n_levels: 2,
        estimand: energy_estimand,
        s: s_weights.as_slice(),
        min_weight,
        lambda: weight_penalty,
        moment_covs: moment_covs.as_slice(),
        n_moments: q,
        targets: targets.as_slice(),
        tols: tols.as_slice(),
        threads: opts.threads,
        qp: opts.qp,
    };
    let result = balancing_core::methods::energy::solve_discrete(&inputs, &interrupt::pending);
    energy_result_list(&result)
}

/// Solve a multi-category-exposure energy balancing problem.
///
/// `treat_idx` holds the zero-based level of each unit; `focal` is the focal
/// level index used by `att` and ignored by `ate`; `estimand` is `ate` or `att`.
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
fn solve_energy_multi(
    covs: RealSexp,
    treat_idx: IntegerSexp,
    focal: i32,
    s_weights: RealSexp,
    distance: &str,
    estimand: &str,
    improved: bool,
    moment_covs: RealSexp,
    targets: RealSexp,
    tols: RealSexp,
    min_weight: f64,
    weight_penalty: f64,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let (p, q) = energy_dims(n, &covs, &s_weights, &moment_covs, &targets, &tols)?;
    if treat_idx.len() != n {
        return Err(savvy::Error::new("treat_idx must have length n"));
    }
    let levels = treat_idx.as_slice();
    if levels.iter().any(|&t| t < 0) {
        return Err(savvy::Error::new("treat_idx values must be non-negative"));
    }
    let n_levels = levels.iter().copied().max().map_or(0, |m| m + 1) as usize;
    if focal < 0 || focal as usize >= n_levels {
        return Err(savvy::Error::new(
            "focal must be a level present in treat_idx",
        ));
    }
    let opts = parse_qp_options(options)?;
    let dist = parse_distance(distance)?;
    let energy_estimand = match estimand {
        "ate" => EnergyEstimand::Ate { improved },
        "att" | "atc" => EnergyEstimand::Focal {
            focal: focal as usize,
        },
        other => {
            return Err(savvy::Error::new(format!(
                "unknown estimand `{other}`; expected ate or att"
            )));
        }
    };

    let inputs = EnergyDiscreteInputs {
        covs: covs.as_slice(),
        n,
        p,
        distance: dist,
        levels,
        n_levels,
        estimand: energy_estimand,
        s: s_weights.as_slice(),
        min_weight,
        lambda: weight_penalty,
        moment_covs: moment_covs.as_slice(),
        n_moments: q,
        targets: targets.as_slice(),
        tols: tols.as_slice(),
        threads: opts.threads,
        qp: opts.qp,
    };
    let result = balancing_core::methods::energy::solve_discrete(&inputs, &interrupt::pending);
    energy_result_list(&result)
}

/// Solve a continuous-exposure energy balancing problem.
///
/// `treat` holds the continuous exposure. `d_covs` and `d_treat` are the
/// distribution-moment columns for the covariates and the exposure, centered and
/// scaled on the R side, and each is held exactly at a weighted mean of zero.
/// The R side centers them on the base measure, so pinning the rows at zero
/// holds every marginal at its base-measure sample value; `bal_covs` and
/// `bal_tols` are the correlation-constraint covariates and their tolerances.
/// `dimension_adj` weights the covariate energy distance by the dimensionality
/// adjustment.
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
fn solve_energy_cont(
    covs: RealSexp,
    treat: RealSexp,
    s_weights: RealSexp,
    distance: &str,
    dimension_adj: bool,
    min_weight: f64,
    weight_penalty: f64,
    d_covs: RealSexp,
    d_treat: RealSexp,
    bal_covs: RealSexp,
    bal_tols: RealSexp,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    if n == 0 || covs.len() % n != 0 {
        return Err(savvy::Error::new(format!(
            "covs has {} elements but is not a multiple of n = {n}",
            covs.len()
        )));
    }
    let p = covs.len() / n;
    if treat.len() != n {
        return Err(savvy::Error::new("treat must have length n"));
    }
    let n_d_covs = if d_covs.is_empty() {
        0
    } else {
        d_covs.len() / n
    };
    let n_d_treat = if d_treat.is_empty() {
        0
    } else {
        d_treat.len() / n
    };
    let n_bal = bal_tols.len();
    if d_covs.len() != n * n_d_covs || d_treat.len() != n * n_d_treat || bal_covs.len() != n * n_bal
    {
        return Err(savvy::Error::new(
            "d_covs, d_treat, and bal_covs must each have n rows",
        ));
    }
    let opts = parse_qp_options(options)?;
    let dist = parse_distance(distance)?;

    let inputs = EnergyContInputs {
        covs: covs.as_slice(),
        n,
        p,
        treat: treat.as_slice(),
        distance: dist,
        s: s_weights.as_slice(),
        min_weight,
        lambda: weight_penalty,
        dimension_adj,
        d_covs: d_covs.as_slice(),
        n_d_covs,
        d_treat: d_treat.as_slice(),
        n_d_treat,
        bal_covs: bal_covs.as_slice(),
        n_bal,
        bal_tols: bal_tols.as_slice(),
        threads: opts.threads,
        qp: opts.qp,
    };
    let result = balancing_core::methods::energy::solve_cont(&inputs, &interrupt::pending);
    energy_result_list(&result)
}

/// Pack a stable balancing result into its R list.
///
/// The list is `weights`, `duals`, `converged`, `iterations`, `objective`,
/// `solver_status` (the backend identity), `status` (the solver's terminal
/// status name), `pri_res`, `dua_res`, and the interrupt flag, in that order. The
/// quadratic-program family has no estimating equations, so none appear.
fn sbw_result_list(result: &SbwResult) -> savvy::Result<savvy::Sexp> {
    let mut out = OwnedListSexp::new(11, true)?;
    out.set_name_and_value(0, "weights", real_vector(&result.weights)?)?;
    out.set_name_and_value(1, "duals", real_vector(&result.duals)?)?;
    out.set_name_and_value(2, "converged", scalar_logical(result.converged)?)?;
    out.set_name_and_value(3, "iterations", scalar_integer(result.iterations)?)?;
    out.set_name_and_value(4, "objective", scalar_real(result.objective)?)?;
    out.set_name_and_value(5, "solver_status", scalar_string(result.backend)?)?;
    out.set_name_and_value(6, "status", scalar_string(result.status)?)?;
    out.set_name_and_value(7, "pri_res", scalar_real(result.pri_res)?)?;
    out.set_name_and_value(8, "dua_res", scalar_real(result.dua_res)?)?;
    out.set_name_and_value(9, "interrupted", scalar_logical(result.interrupted)?)?;
    out.set_name_and_value(10, "fell_back", scalar_logical(result.fell_back)?)?;
    Ok(out.into())
}

/// Validate the shared shape of the discrete stable balancing moment inputs,
/// returning the moment-constraint column count `q`.
fn sbw_moment_dims(
    n: usize,
    s_weights: &RealSexp,
    moment_covs: &RealSexp,
    targets: &RealSexp,
    tols: &RealSexp,
) -> savvy::Result<usize> {
    if s_weights.len() != n {
        return Err(savvy::Error::new("s_weights must have length n"));
    }
    let q = targets.len();
    if tols.len() != q || moment_covs.len() != n * q {
        return Err(savvy::Error::new(
            "moment_covs must be n by length(targets), and tols must match targets",
        ));
    }
    Ok(q)
}

/// Solve a binary-exposure stable balancing problem.
///
/// `treat` holds the zero/one exposure indicator; `estimand` is one of `ate`,
/// `att`, or `atc`; `norm` names the dispersion norm to minimize. `moment_covs`
/// are the standardized balance columns with their `targets` and `tols`, all
/// empty when no balance constraints are requested.
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
fn solve_sbw(
    treat: IntegerSexp,
    s_weights: RealSexp,
    estimand: &str,
    norm: &str,
    moment_covs: RealSexp,
    targets: RealSexp,
    tols: RealSexp,
    min_weight: f64,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let q = sbw_moment_dims(n, &s_weights, &moment_covs, &targets, &tols)?;
    if treat.len() != n {
        return Err(savvy::Error::new("treat must have length n"));
    }
    let treat_slice = treat.as_slice();
    if treat_slice.iter().any(|&t| t != 0 && t != 1) {
        return Err(savvy::Error::new("treat must hold only zero and one"));
    }
    let qp = parse_sbw_options(options)?;
    let norm = parse_sbw_norm(norm)?;
    let sbw_estimand = match estimand {
        "ate" => SbwEstimand::Ate,
        "att" => SbwEstimand::Focal { focal: 1 },
        "atc" => SbwEstimand::Focal { focal: 0 },
        other => {
            return Err(savvy::Error::new(format!(
                "unknown estimand `{other}`; expected ate, att, or atc"
            )));
        }
    };

    let inputs = SbwDiscreteInputs {
        n,
        levels: treat_slice,
        n_levels: 2,
        estimand: sbw_estimand,
        s: s_weights.as_slice(),
        min_weight,
        norm,
        moment_covs: moment_covs.as_slice(),
        n_moments: q,
        targets: targets.as_slice(),
        tols: tols.as_slice(),
        qp,
    };
    let result = sbw::solve_discrete(&inputs, &interrupt::pending).map_err(savvy::Error::new)?;
    sbw_result_list(&result)
}

/// Solve a multi-category-exposure stable balancing problem.
///
/// `treat_idx` holds the zero-based level of each unit; `focal` is the focal
/// level index used by `att` and ignored by `ate`; `estimand` is `ate` or `att`.
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
fn solve_sbw_multi(
    treat_idx: IntegerSexp,
    focal: i32,
    s_weights: RealSexp,
    estimand: &str,
    norm: &str,
    moment_covs: RealSexp,
    targets: RealSexp,
    tols: RealSexp,
    min_weight: f64,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let q = sbw_moment_dims(n, &s_weights, &moment_covs, &targets, &tols)?;
    if treat_idx.len() != n {
        return Err(savvy::Error::new("treat_idx must have length n"));
    }
    let levels = treat_idx.as_slice();
    if levels.iter().any(|&t| t < 0) {
        return Err(savvy::Error::new("treat_idx values must be non-negative"));
    }
    let n_levels = levels.iter().copied().max().map_or(0, |m| m + 1) as usize;
    if focal < 0 || focal as usize >= n_levels {
        return Err(savvy::Error::new(
            "focal must be a level present in treat_idx",
        ));
    }
    let qp = parse_sbw_options(options)?;
    let norm = parse_sbw_norm(norm)?;
    let sbw_estimand = match estimand {
        "ate" => SbwEstimand::Ate,
        "att" | "atc" => SbwEstimand::Focal {
            focal: focal as usize,
        },
        other => {
            return Err(savvy::Error::new(format!(
                "unknown estimand `{other}`; expected ate or att"
            )));
        }
    };

    let inputs = SbwDiscreteInputs {
        n,
        levels,
        n_levels,
        estimand: sbw_estimand,
        s: s_weights.as_slice(),
        min_weight,
        norm,
        moment_covs: moment_covs.as_slice(),
        n_moments: q,
        targets: targets.as_slice(),
        tols: tols.as_slice(),
        qp,
    };
    let result = sbw::solve_discrete(&inputs, &interrupt::pending).map_err(savvy::Error::new)?;
    sbw_result_list(&result)
}

/// Solve a continuous-exposure stable balancing problem.
///
/// `treat` holds the continuous exposure; `covs` are the standardized covariate
/// columns held in weighted correlation with the exposure within `tols`; `norm`
/// names the dispersion norm to minimize.
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
fn solve_sbw_cont(
    treat: RealSexp,
    covs: RealSexp,
    s_weights: RealSexp,
    norm: &str,
    tols: RealSexp,
    min_weight: f64,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    if treat.len() != n {
        return Err(savvy::Error::new("treat must have length n"));
    }
    let n_covs = tols.len();
    if covs.len() != n * n_covs {
        return Err(savvy::Error::new("covs must be n by length(tols)"));
    }
    let qp = parse_sbw_options(options)?;
    let norm = parse_sbw_norm(norm)?;

    let inputs = SbwContInputs {
        n,
        treat: treat.as_slice(),
        covs: covs.as_slice(),
        n_covs,
        s: s_weights.as_slice(),
        min_weight,
        norm,
        tols: tols.as_slice(),
        qp,
    };
    let result = sbw::solve_cont(&inputs, &interrupt::pending).map_err(savvy::Error::new)?;
    sbw_result_list(&result)
}

/// Pack a characteristic function distance result into its R list.
///
/// The list is `weights`, `duals`, `converged`, `iterations`, `objective`,
/// `solver_status` (the backend identity), `status` (the solver's terminal
/// status name), `pri_res`, `dua_res`, the interrupt flag, and the fallback flag,
/// in that order. The quadratic-program family has no estimating equations, so
/// none appear.
fn cfd_result_list(result: &CfdResult) -> savvy::Result<savvy::Sexp> {
    let mut out = OwnedListSexp::new(11, true)?;
    out.set_name_and_value(0, "weights", real_vector(&result.weights)?)?;
    out.set_name_and_value(1, "duals", real_vector(&result.duals)?)?;
    out.set_name_and_value(2, "converged", scalar_logical(result.converged)?)?;
    out.set_name_and_value(3, "iterations", scalar_integer(result.iterations)?)?;
    out.set_name_and_value(4, "objective", scalar_real(result.objective)?)?;
    out.set_name_and_value(5, "solver_status", scalar_string(result.backend)?)?;
    out.set_name_and_value(6, "status", scalar_string(result.status)?)?;
    out.set_name_and_value(7, "pri_res", scalar_real(result.pri_res)?)?;
    out.set_name_and_value(8, "dua_res", scalar_real(result.dua_res)?)?;
    out.set_name_and_value(9, "interrupted", scalar_logical(result.interrupted)?)?;
    out.set_name_and_value(10, "fell_back", scalar_logical(result.fell_back)?)?;
    Ok(out.into())
}

/// Resolve the Monte Carlo draw count for a kernel and validate the projection
/// matrix. The projections carry only for the t kernel, where they are a column-
/// major `p` by `n_draws` matrix drawn on the R side; the other kernels take an
/// empty matrix and no draws.
fn kernel_draw_count(kernel: Kernel, t_proj: &RealSexp, p: usize) -> savvy::Result<usize> {
    if kernel != Kernel::T {
        return Ok(0);
    }
    if p == 0 {
        return Err(savvy::Error::new(
            "the t kernel requires at least one covariate column",
        ));
    }
    if t_proj.len() % p != 0 {
        return Err(savvy::Error::new(
            "t_proj must be a matrix with one row per covariate column",
        ));
    }
    let n_draws = t_proj.len() / p;
    if n_draws == 0 {
        return Err(savvy::Error::new(
            "the t kernel requires at least one Monte Carlo draw",
        ));
    }
    Ok(n_draws)
}

/// Solve a binary-exposure characteristic function distance balancing problem.
///
/// `treat` holds the zero/one exposure indicator; `kernel` names the kernel;
/// `bw_scale` scales the median bandwidth; `smoothness` is the Matern smoothness;
/// `t_proj` is the column-major `p` by `n_draws` t-kernel projection matrix, empty
/// for the other kernels; `improved` selects the between-group term of the improved
/// average-treatment-effect variant; `estimand` is one of `ate`, `att`, or `atc`.
/// `moment_covs` are the standardized moment-constraint columns with `targets` and
/// `tols`, both empty when no moment constraints are requested.
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
fn solve_cfd(
    covs: RealSexp,
    treat: IntegerSexp,
    s_weights: RealSexp,
    kernel: &str,
    bw_scale: f64,
    smoothness: f64,
    t_proj: RealSexp,
    improved: bool,
    estimand: &str,
    moment_covs: RealSexp,
    targets: RealSexp,
    tols: RealSexp,
    min_weight: f64,
    weight_penalty: f64,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let (p, q) = energy_dims(n, &covs, &s_weights, &moment_covs, &targets, &tols)?;
    if treat.len() != n {
        return Err(savvy::Error::new("treat must have length n"));
    }
    let treat_slice = treat.as_slice();
    if treat_slice.iter().any(|&t| t != 0 && t != 1) {
        return Err(savvy::Error::new("treat must hold only zero and one"));
    }
    let opts = parse_cfd_options(options)?;
    let kernel = parse_kernel(kernel)?;
    let matern_nu = parse_smoothness(smoothness)?;
    let n_draws = kernel_draw_count(kernel, &t_proj, p)?;
    let cfd_estimand = match estimand {
        "ate" => CfdEstimand::Ate { improved },
        "att" => CfdEstimand::Focal { focal: 1 },
        "atc" => CfdEstimand::Focal { focal: 0 },
        other => {
            return Err(savvy::Error::new(format!(
                "unknown estimand `{other}`; expected ate, att, or atc"
            )));
        }
    };

    let inputs = CfdDiscreteInputs {
        covs: covs.as_slice(),
        n,
        p,
        kernel: KernelParams {
            kernel,
            bw_scale,
            matern_nu,
            t_proj: t_proj.as_slice(),
            n_draws,
        },
        levels: treat_slice,
        n_levels: 2,
        estimand: cfd_estimand,
        s: s_weights.as_slice(),
        min_weight,
        lambda: weight_penalty,
        moment_covs: moment_covs.as_slice(),
        n_moments: q,
        targets: targets.as_slice(),
        tols: tols.as_slice(),
        threads: opts.threads,
        qp: opts.qp,
    };
    let result = cfd::solve_discrete(&inputs, &interrupt::pending).map_err(savvy::Error::new)?;
    cfd_result_list(&result)
}

/// Solve a multi-category-exposure characteristic function distance balancing
/// problem.
///
/// `treat_idx` holds the zero-based level of each unit; `focal` is the focal level
/// index used by `att` and ignored by `ate`; `estimand` is `ate` or `att`. The
/// kernel arguments match [`solve_cfd`].
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
fn solve_cfd_multi(
    covs: RealSexp,
    treat_idx: IntegerSexp,
    focal: i32,
    s_weights: RealSexp,
    kernel: &str,
    bw_scale: f64,
    smoothness: f64,
    t_proj: RealSexp,
    improved: bool,
    estimand: &str,
    moment_covs: RealSexp,
    targets: RealSexp,
    tols: RealSexp,
    min_weight: f64,
    weight_penalty: f64,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    let (p, q) = energy_dims(n, &covs, &s_weights, &moment_covs, &targets, &tols)?;
    if treat_idx.len() != n {
        return Err(savvy::Error::new("treat_idx must have length n"));
    }
    let levels = treat_idx.as_slice();
    if levels.iter().any(|&t| t < 0) {
        return Err(savvy::Error::new("treat_idx values must be non-negative"));
    }
    let n_levels = levels.iter().copied().max().map_or(0, |m| m + 1) as usize;
    if focal < 0 || focal as usize >= n_levels {
        return Err(savvy::Error::new(
            "focal must be a level present in treat_idx",
        ));
    }
    let opts = parse_cfd_options(options)?;
    let kernel = parse_kernel(kernel)?;
    let matern_nu = parse_smoothness(smoothness)?;
    let n_draws = kernel_draw_count(kernel, &t_proj, p)?;
    let cfd_estimand = match estimand {
        "ate" => CfdEstimand::Ate { improved },
        "att" | "atc" => CfdEstimand::Focal {
            focal: focal as usize,
        },
        other => {
            return Err(savvy::Error::new(format!(
                "unknown estimand `{other}`; expected ate or att"
            )));
        }
    };

    let inputs = CfdDiscreteInputs {
        covs: covs.as_slice(),
        n,
        p,
        kernel: KernelParams {
            kernel,
            bw_scale,
            matern_nu,
            t_proj: t_proj.as_slice(),
            n_draws,
        },
        levels,
        n_levels,
        estimand: cfd_estimand,
        s: s_weights.as_slice(),
        min_weight,
        lambda: weight_penalty,
        moment_covs: moment_covs.as_slice(),
        n_moments: q,
        targets: targets.as_slice(),
        tols: tols.as_slice(),
        threads: opts.threads,
        qp: opts.qp,
    };
    let result = cfd::solve_discrete(&inputs, &interrupt::pending).map_err(savvy::Error::new)?;
    cfd_result_list(&result)
}

/// Build a kernel matrix for the covariates.
///
/// `kernel` names the kernel; `bw_scale` scales the median bandwidth;
/// `smoothness` is the Matern smoothness; `t_proj` is the column-major `p` by
/// `n_draws` t-kernel projection matrix, empty for the other kernels; `s_weights`
/// standardize the covariates; `discarded` marks units excluded from the bandwidth
/// median, empty to discard none. The result is a column-major `n` by `n` symmetric
/// matrix. This entry point is mandated by the boundary contract and is available
/// for building a kernel matrix directly, for diagnostics. The fit path does not
/// call it: the solve entry points build the kernel internally, so the `n` by `n`
/// matrix never crosses the boundary during a fit.
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
fn kernel_matrix(
    covs: RealSexp,
    kernel: &str,
    bw_scale: f64,
    smoothness: f64,
    t_proj: RealSexp,
    s_weights: RealSexp,
    discarded: LogicalSexp,
    options: ListSexp,
) -> savvy::Result<savvy::Sexp> {
    let n = s_weights.len();
    if n == 0 || covs.len() % n != 0 {
        return Err(savvy::Error::new(format!(
            "covs has {} elements but is not a multiple of n = {n}",
            covs.len()
        )));
    }
    let p = covs.len() / n;
    if !discarded.is_empty() && discarded.len() != n {
        return Err(savvy::Error::new(
            "discarded must be empty or have length n",
        ));
    }
    let threads = parse_kernel_options(options)?;
    let kernel = parse_kernel(kernel)?;
    let matern_nu = parse_smoothness(smoothness)?;
    let n_draws = kernel_draw_count(kernel, &t_proj, p)?;
    let discarded_mask: Vec<bool> = discarded.iter().collect();

    let params = KernelParams {
        kernel,
        bw_scale,
        matern_nu,
        t_proj: t_proj.as_slice(),
        n_draws,
    };
    let matrix = build_kernel(
        covs.as_slice(),
        n,
        p,
        s_weights.as_slice(),
        &params,
        &discarded_mask,
        threads,
    );
    Ok(real_matrix(&matrix, n, n)?.into())
}

#[cfg(test)]
mod tests {
    #[test]
    fn core_reports_a_positive_thread_count() {
        let (available, _cap_source) = balancing_core::available_threads();
        assert!(available >= 1);
    }
}
