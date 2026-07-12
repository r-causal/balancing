//! Conversions between R data and the core's plain Rust types.
//!
//! Inputs are read through borrowed slices with no copy; outputs are allocated
//! once as owned SEXPs and filled directly. Option lists are parsed strictly:
//! an unrecognized option name is a contract violation and becomes an error.

use balancing_core::methods::entropy::EntropySolver;
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
