//! Conversions between R data and the core's plain Rust types.
//!
//! Inputs are read through borrowed slices with no copy; outputs are allocated
//! once as owned SEXPs and filled directly. Option lists are parsed strictly:
//! an unrecognized option name is a contract violation and becomes an error.

use balancing_core::dist::Distance;
use balancing_core::dist::kernels::Kernel;
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

    let mut resolved = EntropyOptions {
        threads: balancing_core::available_threads().0,
        solver: EntropySolver::Newton,
        max_iter: 200,
        tol: 1e-10,
        esteq_scale: None,
    };

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

    let mut resolved = IptOptions {
        threads: balancing_core::available_threads().0,
        max_iter: 200,
        tol: 1e-10,
    };

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
