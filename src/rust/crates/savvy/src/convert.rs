//! Conversions between R data and the core's plain Rust types.
//!
//! Inputs are read through borrowed slices with no copy; outputs are allocated
//! once as owned SEXPs and filled directly. Option lists are parsed strictly:
//! an unrecognized option name is a contract violation and becomes an error.

use balancing_core::dist::Distance;
use balancing_core::dist::kernels::{Kernel, MaternNu};
use balancing_core::links::Link;
use balancing_core::methods::cbps::CbpsEstimand;
use balancing_core::methods::entropy::EntropySolver;
use balancing_core::methods::ipt::IptEstimand;
use balancing_core::methods::sbw::SbwNorm;
use balancing_core::qp::{QpBackendChoice, QpOptions};
use savvy::{ListSexp, OwnedRealSexp, RealSexp, Sexp};

/// Resolved solver options shared by the entropy entrypoints.
pub struct EntropyOptions {
    pub threads: usize,
    pub solver: EntropySolver,
    pub max_iter: usize,
    pub tol: f64,
    /// Per-parameter-block scalar applied to the estimating-equation output
    /// before it crosses back to R, one value per solved group. Absent when the
    /// caller does no renormalization.
    pub esteq_scale: Option<Vec<f64>>,
}

/// The entropy defaults the boundary resolves when the option list omits a value.
///
/// The two solver fields are reached by different routes. `max_iterations` has no
/// R-side default, so it is absent from the option list on every fit whose caller
/// has not named one, which makes the cap here the budget every default entropy
/// fit runs under. `convergence_tolerance` does carry an R-side default and is
/// forwarded whenever it is non-NULL, so the tolerance here applies only to a
/// caller who passes NULL explicitly.
///
/// Both are deliberately separate from `SolveOptions::default()`, which serves
/// direct core callers and is stricter on the gradient while allowing fewer
/// iterations to reach it.
impl Default for EntropyOptions {
    fn default() -> Self {
        Self {
            threads: balancing_core::available_threads().0,
            solver: EntropySolver::Newton,
            max_iter: 1000,
            tol: 1e-10,
            esteq_scale: None,
        }
    }
}

/// Read a scalar option as an `f64`, accepting either an integer or a double.
fn option_f64(value: Sexp, name: &str) -> savvy::Result<f64> {
    if value.is_integer() {
        Ok(f64::from(i32::try_from(value)?))
    } else if value.is_real() {
        f64::try_from(value)
    } else {
        Err(savvy::Error::new(format!(
            "option `{name}` must be a numeric scalar"
        )))
    }
}

/// Read a scalar option as a `usize`, accepting either an integer or a double.
fn option_usize(value: Sexp, name: &str) -> savvy::Result<usize> {
    let x = option_f64(value, name)?;
    if x < 0.0 || x.fract() != 0.0 {
        return Err(savvy::Error::new(format!(
            "option `{name}` must be a non-negative whole number"
        )));
    }
    Ok(x as usize)
}

/// Parse the option list for an entropy solve, rejecting unknown names.
pub fn parse_entropy_options(options: ListSexp) -> savvy::Result<EntropyOptions> {
    const ALLOWED: [&str; 5] = [
        "threads",
        "solver",
        "max_iterations",
        "convergence_tolerance",
        "esteq_scale",
    ];

    for name in options.names_iter() {
        if !ALLOWED.contains(&name) {
            return Err(savvy::Error::new(format!(
                "unknown option `{name}`; allowed options are {}",
                ALLOWED.join(", ")
            )));
        }
    }

    let mut resolved = EntropyOptions::default();

    if let Some(value) = options.get("threads") {
        resolved.threads = option_usize(value, "threads")?.max(1);
    }
    if let Some(value) = options.get("max_iterations") {
        resolved.max_iter = option_usize(value, "max_iterations")?.max(1);
    }
    if let Some(value) = options.get("convergence_tolerance") {
        resolved.tol = option_f64(value, "convergence_tolerance")?;
    }
    if let Some(value) = options.get("solver") {
        let name = <&str>::try_from(value)
            .map_err(|_| savvy::Error::new("option `solver` must be a string"))?;
        resolved.solver = match name {
            "newton" => EntropySolver::Newton,
            "lbfgs" => EntropySolver::Lbfgs,
            "lbfgs_then_newton" => EntropySolver::LbfgsThenNewton,
            other => {
                return Err(savvy::Error::new(format!(
                    "unknown solver `{other}`; expected newton, lbfgs, or lbfgs_then_newton"
                )));
            }
        };
    }
    if let Some(value) = options.get("esteq_scale") {
        let scale = RealSexp::try_from(value)
            .map_err(|_| savvy::Error::new("option `esteq_scale` must be a numeric vector"))?;
        resolved.esteq_scale = Some(scale.as_slice().to_vec());
    }

    Ok(resolved)
}

/// Resolved solver options for the inverse probability tilting entrypoints.
pub struct IptOptions {
    pub threads: usize,
    pub max_iter: usize,
    pub tol: f64,
}

/// The defaults the boundary resolves when the option list omits a value, reached
/// by the same two routes as the entropy defaults above.
///
/// This parser serves the covariate balancing propensity score entrypoints as well
/// as the tilting ones, so the cap here is the budget for both families, and both
/// R constructors leave `max_iterations` unset by default. It matches the entropy
/// cap so that methods solving the same fit at the same tolerance are given the
/// same budget to reach it.
impl Default for IptOptions {
    fn default() -> Self {
        Self {
            threads: balancing_core::available_threads().0,
            max_iter: 1000,
            tol: 1e-10,
        }
    }
}

/// Parse the option list for an inverse probability tilting solve, rejecting
/// unknown names. The link and estimand cross the boundary as their own
/// arguments, so the option list carries only the solver tuning.
pub fn parse_ipt_options(options: ListSexp) -> savvy::Result<IptOptions> {
    const ALLOWED: [&str; 3] = ["threads", "max_iterations", "convergence_tolerance"];

    for name in options.names_iter() {
        if !ALLOWED.contains(&name) {
            return Err(savvy::Error::new(format!(
                "unknown option `{name}`; allowed options are {}",
                ALLOWED.join(", ")
            )));
        }
    }

    let mut resolved = IptOptions::default();

    if let Some(value) = options.get("threads") {
        resolved.threads = option_usize(value, "threads")?.max(1);
    }
    if let Some(value) = options.get("max_iterations") {
        resolved.max_iter = option_usize(value, "max_iterations")?.max(1);
    }
    if let Some(value) = options.get("convergence_tolerance") {
        resolved.tol = option_f64(value, "convergence_tolerance")?;
    }

    Ok(resolved)
}

/// Resolved options for an energy balancing solve: the worker-thread count and
/// the quadratic-program tuning.
pub struct EnergyOptions {
    pub threads: usize,
    pub qp: QpOptions,
}

/// Parse the option list for an energy balancing solve, rejecting unknown names.
///
/// `convergence_tolerance` sets both the absolute and relative solver tolerances;
/// `backend` is accepted for forward compatibility but energy balancing always
/// solves through the ADMM backend, so a value other than `osqp` or `auto` is an
/// error rather than a silent override.
pub fn parse_qp_options(options: ListSexp) -> savvy::Result<EnergyOptions> {
    const ALLOWED: [&str; 5] = [
        "threads",
        "convergence_tolerance",
        "max_iterations",
        "polish",
        "backend",
    ];

    for name in options.names_iter() {
        if !ALLOWED.contains(&name) {
            return Err(savvy::Error::new(format!(
                "unknown option `{name}`; allowed options are {}",
                ALLOWED.join(", ")
            )));
        }
    }

    let mut qp = QpOptions::default();
    let mut threads = balancing_core::available_threads().0;

    if let Some(value) = options.get("threads") {
        threads = option_usize(value, "threads")?.max(1);
    }
    if let Some(value) = options.get("convergence_tolerance") {
        let tol = option_f64(value, "convergence_tolerance")?;
        qp.eps_abs = tol;
        qp.eps_rel = tol;
    }
    if let Some(value) = options.get("max_iterations") {
        qp.max_iter = option_usize(value, "max_iterations")?.max(1);
    }
    if let Some(value) = options.get("polish") {
        qp.polish = bool::try_from(value)
            .map_err(|_| savvy::Error::new("option `polish` must be a logical scalar"))?;
    }
    if let Some(value) = options.get("backend") {
        let name = <&str>::try_from(value)
            .map_err(|_| savvy::Error::new("option `backend` must be a string"))?;
        if name != "osqp" && name != "auto" {
            return Err(savvy::Error::new(format!(
                "unknown backend `{name}`; energy balancing solves through `osqp`"
            )));
        }
    }

    Ok(EnergyOptions { threads, qp })
}

/// Parse the option list for a stable balancing solve, rejecting unknown names.
///
/// The tuning matches the other quadratic-program methods: `convergence_tolerance`
/// sets both the absolute and relative solver tolerances, `max_iterations` caps
/// the iterations, and `polish` toggles the solution polish. `backend` is one of
/// `auto`, `osqp`, or `clarabel`; `auto` (the default) solves with osqp first and
/// re-solves with clarabel on an osqp primal-infeasibility certificate, and any
/// other value is an error rather than a silent override.
pub fn parse_sbw_options(options: ListSexp) -> savvy::Result<QpOptions> {
    const ALLOWED: [&str; 5] = [
        "threads",
        "convergence_tolerance",
        "max_iterations",
        "polish",
        "backend",
    ];

    for name in options.names_iter() {
        if !ALLOWED.contains(&name) {
            return Err(savvy::Error::new(format!(
                "unknown option `{name}`; allowed options are {}",
                ALLOWED.join(", ")
            )));
        }
    }

    let mut qp = QpOptions::default();

    if let Some(value) = options.get("convergence_tolerance") {
        let tol = option_f64(value, "convergence_tolerance")?;
        qp.eps_abs = tol;
        qp.eps_rel = tol;
    }
    if let Some(value) = options.get("max_iterations") {
        qp.max_iter = option_usize(value, "max_iterations")?.max(1);
    }
    if let Some(value) = options.get("polish") {
        qp.polish = bool::try_from(value)
            .map_err(|_| savvy::Error::new("option `polish` must be a logical scalar"))?;
    }
    if let Some(value) = options.get("backend") {
        let name = <&str>::try_from(value)
            .map_err(|_| savvy::Error::new("option `backend` must be a string"))?;
        qp.backend = match name {
            "auto" => QpBackendChoice::Auto,
            "osqp" => QpBackendChoice::Osqp,
            "clarabel" => QpBackendChoice::Clarabel,
            other => {
                return Err(savvy::Error::new(format!(
                    "unknown backend `{other}`; expected auto, osqp, or clarabel"
                )));
            }
        };
    }

    Ok(qp)
}

/// Resolve the dispersion norm named by the R layer.
pub fn parse_sbw_norm(norm: &str) -> savvy::Result<SbwNorm> {
    match norm {
        "l2" => Ok(SbwNorm::L2),
        "l1" => Ok(SbwNorm::L1),
        "linf" => Ok(SbwNorm::Linf),
        other => Err(savvy::Error::new(format!(
            "unknown norm `{other}`; expected l2, l1, or linf"
        ))),
    }
}

/// Parse the option list for a characteristic function distance solve, rejecting
/// unknown names.
///
/// The tuning matches the other quadratic-program methods, and the kernels are
/// positive semidefinite by construction (except the energy kernel), so the full
/// backend set applies: `backend` is one of `auto`, `osqp`, or `clarabel`, `auto`
/// (the default) solving with osqp first and re-solving with clarabel on an osqp
/// primal-infeasibility certificate. Any other value is an error rather than a
/// silent override.
pub fn parse_cfd_options(options: ListSexp) -> savvy::Result<EnergyOptions> {
    const ALLOWED: [&str; 5] = [
        "threads",
        "convergence_tolerance",
        "max_iterations",
        "polish",
        "backend",
    ];

    for name in options.names_iter() {
        if !ALLOWED.contains(&name) {
            return Err(savvy::Error::new(format!(
                "unknown option `{name}`; allowed options are {}",
                ALLOWED.join(", ")
            )));
        }
    }

    let mut qp = QpOptions::default();
    let mut threads = balancing_core::available_threads().0;

    if let Some(value) = options.get("threads") {
        threads = option_usize(value, "threads")?.max(1);
    }
    if let Some(value) = options.get("convergence_tolerance") {
        let tol = option_f64(value, "convergence_tolerance")?;
        qp.eps_abs = tol;
        qp.eps_rel = tol;
    }
    if let Some(value) = options.get("max_iterations") {
        qp.max_iter = option_usize(value, "max_iterations")?.max(1);
    }
    if let Some(value) = options.get("polish") {
        qp.polish = bool::try_from(value)
            .map_err(|_| savvy::Error::new("option `polish` must be a logical scalar"))?;
    }
    if let Some(value) = options.get("backend") {
        let name = <&str>::try_from(value)
            .map_err(|_| savvy::Error::new("option `backend` must be a string"))?;
        qp.backend = match name {
            "auto" => QpBackendChoice::Auto,
            "osqp" => QpBackendChoice::Osqp,
            "clarabel" => QpBackendChoice::Clarabel,
            other => {
                return Err(savvy::Error::new(format!(
                    "unknown backend `{other}`; expected auto, osqp, or clarabel"
                )));
            }
        };
    }

    Ok(EnergyOptions { threads, qp })
}

/// Parse the option list for a standalone kernel-matrix build, rejecting unknown
/// names. The only tuning is the worker-thread count.
pub fn parse_kernel_options(options: ListSexp) -> savvy::Result<usize> {
    const ALLOWED: [&str; 1] = ["threads"];

    for name in options.names_iter() {
        if !ALLOWED.contains(&name) {
            return Err(savvy::Error::new(format!(
                "unknown option `{name}`; allowed options are {}",
                ALLOWED.join(", ")
            )));
        }
    }

    let mut threads = balancing_core::available_threads().0;
    if let Some(value) = options.get("threads") {
        threads = option_usize(value, "threads")?.max(1);
    }
    Ok(threads)
}

/// Resolve the kernel named by the R layer.
pub fn parse_kernel(kernel: &str) -> savvy::Result<Kernel> {
    Kernel::from_name(kernel).ok_or_else(|| {
        savvy::Error::new(format!(
            "unknown kernel `{kernel}`; expected energy, gaussian, laplace, matern, or t"
        ))
    })
}

/// Resolve the Matern smoothness passed by the R layer to a supported
/// half-integer.
///
/// This resolver admits a value within 1e-9 of a closed form and refuses every
/// other one, and the R constructor is written to that rule: it snaps a
/// near-canonical request onto the exact order and its validator then admits only
/// 0.5, 1.5, and 2.5, so an unsupported value cannot arrive through the public
/// API. The two rules have to be kept in step. An R-side tolerance wider than the
/// one here, which is what comparing with `all.equal()` alone gave, passes
/// construction and fails at this boundary instead. A direct core caller is bound
/// by neither rule, which is why resolving reports the value rather than silently
/// substituting a default, matching how `parse_kernel` treats an unrecognized
/// name. The value is ignored for the non-Matern kernels, which pass their default
/// smoothness.
pub fn parse_smoothness(smoothness: f64) -> savvy::Result<MaternNu> {
    MaternNu::from_smoothness(smoothness).ok_or_else(|| {
        savvy::Error::new(format!(
            "unsupported Matern smoothness `{smoothness}`; expected 0.5, 1.5, or 2.5"
        ))
    })
}

/// Resolve the distance definition named by the R layer.
pub fn parse_distance(distance: &str) -> savvy::Result<Distance> {
    Distance::from_name(distance).ok_or_else(|| {
        savvy::Error::new(format!(
            "unknown distance `{distance}`; expected scaled_euclidean, mahalanobis, or euclidean"
        ))
    })
}

/// Resolve the propensity link named by the R layer.
pub fn parse_link(link: &str) -> savvy::Result<Link> {
    Link::from_name(link).ok_or_else(|| {
        savvy::Error::new(format!(
            "unknown link `{link}`; expected logit, probit, or cloglog"
        ))
    })
}

/// Resolve a binary-exposure estimand to its tilting target. The treated level
/// is `1` and the control level is `0`, so `att` tilts toward the treated and
/// `atc` toward the controls.
///
/// The R fit path re-encodes a focal estimand as `att` with the focal exposure
/// level mapped to the treated encoding (`1`), so only `ate` and `att` arrive
/// from the package; the `atc` arm is unreachable from the package and is kept
/// for boundary symmetry with the categorical and quadratic-program parsers.
/// The R layer's own canonical untreated target is `atu`, the propensity synonym
/// for `atc`, which never reaches this parser because of that re-encoding.
pub fn parse_binary_estimand(estimand: &str) -> savvy::Result<IptEstimand> {
    match estimand {
        "ate" => Ok(IptEstimand::Ate),
        "att" => Ok(IptEstimand::Focal(1)),
        "atc" => Ok(IptEstimand::Focal(0)),
        other => Err(savvy::Error::new(format!(
            "unknown estimand `{other}`; expected ate, att, or atc"
        ))),
    }
}

/// Resolve a categorical-exposure estimand to its tilting target, given the
/// focal level index the R layer passes for the focal estimands.
pub fn parse_multi_estimand(estimand: &str, focal: usize) -> savvy::Result<IptEstimand> {
    match estimand {
        "ate" => Ok(IptEstimand::Ate),
        "att" | "atc" => Ok(IptEstimand::Focal(focal)),
        other => Err(savvy::Error::new(format!(
            "unknown estimand `{other}`; expected ate, att, or atc"
        ))),
    }
}

/// Resolve a binary covariate balancing propensity score estimand. The overlap
/// estimand is legal only for a binary exposure, so it is accepted here.
pub fn parse_cbps_estimand(estimand: &str) -> savvy::Result<CbpsEstimand> {
    match estimand {
        "ate" => Ok(CbpsEstimand::Ate),
        "att" => Ok(CbpsEstimand::Att),
        "atc" => Ok(CbpsEstimand::Atc),
        "ato" => Ok(CbpsEstimand::Ato),
        other => Err(savvy::Error::new(format!(
            "unknown estimand `{other}`; expected ate, att, atc, or ato"
        ))),
    }
}

/// Resolve a categorical covariate balancing propensity score estimand, which
/// admits only the average treatment effect and the effect on the treated.
pub fn parse_cbps_multi_estimand(estimand: &str) -> savvy::Result<CbpsEstimand> {
    match estimand {
        "ate" => Ok(CbpsEstimand::Ate),
        "att" => Ok(CbpsEstimand::Att),
        other => Err(savvy::Error::new(format!(
            "unknown estimand `{other}`; expected ate or att"
        ))),
    }
}

/// Reject a non-finite value in a boundary numeric input.
///
/// The R layer validates missing values (`NA`, which includes `NaN`) but not an
/// infinity, which slips through `anyNA` and standardizes to another infinity, so
/// this catches both as defense in depth before the value reaches a solver where
/// it would silently poison the arithmetic. One pass over the slice, cheap
/// relative to the solve that reads it repeatedly.
pub fn require_finite(values: &[f64], name: &str) -> savvy::Result<()> {
    if let Some(i) = values.iter().position(|x| !x.is_finite()) {
        return Err(savvy::Error::new(format!(
            "`{name}` contains a non-finite value at position {}",
            i + 1
        )));
    }
    Ok(())
}

/// Validate a binary treatment vector: every value is `0` or `1`, and both
/// levels are present.
///
/// The R fit path always encodes the treatment as `0`/`1`, so this makes the
/// binary boundary symmetric with the categorical entry (which already checks
/// its encoding) rather than leaving an out-of-range value to produce a unit
/// weight and a `NaN` propensity, or an all-one-level vector to drive a solve on
/// an empty group.
pub fn require_binary_treat(treat: &[i32]) -> savvy::Result<()> {
    let mut seen_zero = false;
    let mut seen_one = false;
    for &t in treat {
        match t {
            0 => seen_zero = true,
            1 => seen_one = true,
            other => {
                return Err(savvy::Error::new(format!(
                    "treat must be 0 or 1; found {other}"
                )));
            }
        }
    }
    if !(seen_zero && seen_one) {
        return Err(savvy::Error::new("treat must contain both levels 0 and 1"));
    }
    Ok(())
}

/// Validate that a categorical treatment vector fills every level `0..n_levels`
/// with at least one unit, so no solved block is left with an empty group.
///
/// The categorical entries derive `n_levels` from the largest level index, so a
/// gap (an unused interior exposure level) or an all-one-level vector would leave
/// a block with no units, whose target and Newton step are degenerate. Callers
/// have already checked non-negativity and the focal range.
pub fn require_dense_levels(treat: &[i32], n_levels: usize) -> savvy::Result<()> {
    let mut present = vec![false; n_levels];
    for &t in treat {
        if let Ok(level) = usize::try_from(t)
            && level < n_levels
        {
            present[level] = true;
        }
    }
    if let Some(level) = present.iter().position(|&seen| !seen) {
        return Err(savvy::Error::new(format!(
            "treat leaves level {level} empty; every level 0..{n_levels} must be present"
        )));
    }
    Ok(())
}

/// Build an R matrix SEXP from column-major data.
pub fn real_matrix(data: &[f64], nrow: usize, ncol: usize) -> savvy::Result<OwnedRealSexp> {
    debug_assert_eq!(data.len(), nrow * ncol);
    let mut out = OwnedRealSexp::new(data.len())?;
    out.as_mut_slice().copy_from_slice(data);
    out.set_dim(&[nrow as i32, ncol as i32])?;
    Ok(out)
}

/// Build an R numeric vector SEXP from a slice.
pub fn real_vector(data: &[f64]) -> savvy::Result<OwnedRealSexp> {
    let mut out = OwnedRealSexp::new(data.len())?;
    out.as_mut_slice().copy_from_slice(data);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    // The iteration cap the boundary resolves is a portability contract rather
    // than a tuning preference, so it is pinned here. An option list arriving
    // without `max_iterations` leaves the parser at these values, and an entropy
    // or tilting fit that needs more iterations than the cap allows stops short
    // of its tolerance and warns. The solvers walk a different floating-point
    // path on each platform, so a cap tight enough that a fit converges just
    // under it on one platform leaves the same fit warning on another. The cap
    // has to be wide enough that the platform spread sits well inside it. The
    // tolerance is pinned alongside it because the two halves only mean
    // something together: a budget is generous or tight only relative to the
    // convergence target it is spent reaching, and the documentation states
    // both, so a change to either has to be a deliberate one.
    #[test]
    fn the_entropy_solver_defaults_pin_the_iteration_cap_and_tolerance() {
        assert_eq!(EntropyOptions::default().max_iter, 1000);
        assert_eq!(EntropyOptions::default().tol, 1e-10);
    }

    #[test]
    fn the_tilting_solver_defaults_pin_the_iteration_cap_and_tolerance() {
        assert_eq!(IptOptions::default().max_iter, 1000);
        assert_eq!(IptOptions::default().tol, 1e-10);
    }

    // The quadratic-program defaults reach the boundary by a different route:
    // the parsers do not carry their own copies, they start from
    // `QpOptions::default()` and overwrite only what the option list names. That
    // makes the core's defaults the ones the sbw, energy and cfd fits run under,
    // and the ones the documentation states, so they are pinned here too.
    #[test]
    fn the_quadratic_program_defaults_pin_the_tolerances_and_iteration_cap() {
        let qp = QpOptions::default();
        assert_eq!(qp.eps_abs, 1e-8);
        assert_eq!(qp.eps_rel, 1e-8);
        assert_eq!(qp.max_iter, 200_000);
    }

    #[test]
    fn parse_smoothness_resolves_the_supported_half_integers() {
        assert_eq!(parse_smoothness(0.5).unwrap(), MaternNu::Half);
        assert_eq!(parse_smoothness(1.5).unwrap(), MaternNu::ThreeHalves);
        assert_eq!(parse_smoothness(2.5).unwrap(), MaternNu::FiveHalves);
    }

    #[test]
    fn parse_smoothness_errors_on_an_unsupported_value() {
        // Unreachable through the S7 validator, but the boundary reports rather
        // than silently substituting a default.
        assert!(parse_smoothness(2.0).is_err());
        assert!(parse_smoothness(3.5).is_err());
    }

    #[test]
    fn require_finite_accepts_finite_and_rejects_nan_and_infinity() {
        assert!(require_finite(&[0.0, -1.5, 3.0], "covs").is_ok());
        let nan = require_finite(&[1.0, f64::NAN, 2.0], "covs").unwrap_err();
        assert!(nan.to_string().contains("position 2"), "{nan}");
        let inf = require_finite(&[1.0, 2.0, f64::INFINITY], "s_weights").unwrap_err();
        assert!(inf.to_string().contains("s_weights"), "{inf}");
    }

    #[test]
    fn require_binary_treat_accepts_both_levels() {
        assert!(require_binary_treat(&[0, 1, 1, 0]).is_ok());
    }

    #[test]
    fn require_binary_treat_rejects_out_of_range_and_single_level() {
        assert!(require_binary_treat(&[0, 2, 1]).is_err());
        assert!(require_binary_treat(&[-1, 0, 1]).is_err());
        assert!(require_binary_treat(&[1, 1, 1]).is_err());
        assert!(require_binary_treat(&[0, 0, 0]).is_err());
    }

    #[test]
    fn require_dense_levels_accepts_a_full_encoding() {
        assert!(require_dense_levels(&[0, 1, 2, 1, 0], 3).is_ok());
    }

    #[test]
    fn require_dense_levels_rejects_a_gap() {
        // Level 1 is unused, so its block would solve on an empty group.
        let err = require_dense_levels(&[0, 2, 2, 0], 3).unwrap_err();
        assert!(err.to_string().contains("level 1"), "{err}");
    }
}
