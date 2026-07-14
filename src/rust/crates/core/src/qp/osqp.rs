//! The OSQP backend: an ADMM solver that is the default for every method.
//!
//! ADMM regularizes the KKT system it factors, so it solves the indefinite
//! quadratic form of energy balancing where an interior-point method fails.
//! OSQP has no mid-solve callback, so the interrupt is honored by solving in
//! bounded iteration chunks and warm starting between them: after each chunk the
//! interrupt closure is polled, and the primal and dual iterates carry into the
//! next chunk so the split does not change the trajectory.

use osqp::{CscMatrix, Problem, Settings, Status};

use super::{
    QpBackend, QpError, QpOptions, QpSolution, QpSpec, QpStatus, objective, upper_triangular_csc,
};

/// Iterations between rho adaptations, fixed so a solve is reproducible across
/// runs. OSQP otherwise adapts rho on elapsed time, which varies run to run even
/// with identical inputs.
const ADAPTIVE_RHO_INTERVAL: u32 = 10;

/// The magnitude OSQP treats as infinity. A bound at or beyond this is an open
/// side; `f64::INFINITY` must be brought to this finite sentinel because the C
/// core reads a literal infinity as invalid rather than as an open bound.
const OSQP_INFTY: f64 = 1e30;

/// Clamp a bound vector so open sides reach the finite sentinel OSQP expects.
fn clamp_bounds(bounds: &[f64]) -> Vec<f64> {
    bounds
        .iter()
        .map(|&b| b.clamp(-OSQP_INFTY, OSQP_INFTY))
        .collect()
}

/// The OSQP backend.
pub struct Osqp;

impl QpBackend for Osqp {
    fn name(&self) -> &'static str {
        "osqp"
    }

    fn solve(
        &self,
        spec: &QpSpec,
        opts: &QpOptions,
        interrupt: &dyn Fn() -> bool,
    ) -> Result<QpSolution, QpError> {
        let n = spec.n;
        let (p_indptr, p_indices, p_values) = upper_triangular_csc(&spec.p, n);
        let p_csc = CscMatrix {
            nrows: n,
            ncols: n,
            indptr: p_indptr.into(),
            indices: p_indices.into(),
            data: p_values.into(),
        };
        let a_csc = CscMatrix {
            nrows: spec.m,
            ncols: n,
            indptr: spec.a_indptr.clone().into(),
            indices: spec.a_indices.clone().into(),
            data: spec.a_values.clone().into(),
        };

        // The chunk cap bounds the work between interrupt checks; OSQP reports
        // reaching it as MaxIterationsReached, which the loop distinguishes from
        // a genuine cap by tracking the total against `opts.max_iter`.
        let chunk = opts.chunk_iters.clamp(1, opts.max_iter.max(1)) as u32;
        let settings = Settings::default()
            .eps_abs(opts.eps_abs)
            .eps_rel(opts.eps_rel)
            .max_iter(chunk)
            .polishing(opts.polish)
            .adaptive_rho(opts.adaptive_rho)
            .adaptive_rho_interval(Some(ADAPTIVE_RHO_INTERVAL))
            .warm_starting(true)
            .verbose(false);

        let l = clamp_bounds(&spec.l);
        let u = clamp_bounds(&spec.u);
        let mut problem = Problem::new(p_csc, &spec.q, a_csc, &l, &u, &settings)
            .map_err(|e| QpError::Setup(format!("{e:?}")))?;

        let mut total_iters = 0usize;
        loop {
            let status = problem.solve();
            total_iters += status.iter() as usize;

            let terminal = matches!(
                status,
                Status::Solved(_)
                    | Status::SolvedInaccurate(_)
                    | Status::PrimalInfeasible(_)
                    | Status::PrimalInfeasibleInaccurate(_)
                    | Status::DualInfeasible(_)
                    | Status::DualInfeasibleInaccurate(_)
                    | Status::NonConvex(_)
                    | Status::TimeLimitReached(_)
            );

            if terminal || total_iters >= opts.max_iter {
                return Ok(finish(spec, &status, total_iters, false));
            }

            // Not terminal and under the cap: this was a chunk boundary. Poll the
            // interrupt, then warm start the next chunk from the current iterate.
            if interrupt() {
                return Ok(finish(spec, &status, total_iters, true));
            }
            let (x, y) = match status.solution() {
                Some(sol) => (sol.x().to_vec(), sol.y().to_vec()),
                None => return Ok(finish(spec, &status, total_iters, false)),
            };
            problem.warm_start(&x, &y);
        }
    }
}

/// Build a [`QpSolution`] from an OSQP status. When the solve was interrupted or
/// produced no usable iterate the primal is reported as zeros so the caller can
/// react to the status rather than to a partial vector.
fn finish(spec: &QpSpec, status: &Status<'_>, iterations: usize, interrupted: bool) -> QpSolution {
    let n = spec.n;
    let m = spec.m;
    let (x, duals) = match status.solution() {
        Some(sol) => (sol.x().to_vec(), sol.y().to_vec()),
        None => (vec![0.0; n], vec![0.0; m]),
    };
    let mapped = if interrupted {
        QpStatus::Interrupted
    } else {
        map_status(status)
    };
    let obj = objective(spec, &x);
    QpSolution {
        x,
        duals,
        status: mapped,
        iterations,
        obj,
        pri_res: 0.0,
        dua_res: 0.0,
        interrupted,
    }
}

/// Translate an OSQP status into the backend-agnostic [`QpStatus`].
fn map_status(status: &Status<'_>) -> QpStatus {
    match status {
        Status::Solved(_) => QpStatus::Solved,
        Status::SolvedInaccurate(_) => QpStatus::SolvedInaccurate,
        Status::MaxIterationsReached(_) => QpStatus::MaxIter,
        Status::PrimalInfeasible(_) | Status::PrimalInfeasibleInaccurate(_) => {
            QpStatus::PrimalInfeasible
        }
        Status::DualInfeasible(_) | Status::DualInfeasibleInaccurate(_) => QpStatus::DualInfeasible,
        Status::NonConvex(_) => QpStatus::NonConvex,
        Status::TimeLimitReached(_) => QpStatus::TimeLimit,
        _ => QpStatus::MaxIter,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qp::{Convexity, PMat};

    // Minimize 0.5 (x0^2 + x1^2) subject to x0 + x1 = 1 and x >= 0. The optimum
    // is x = (0.5, 0.5). The spec stores the doubled quadratic term diag(1, 1)
    // so 0.5 x' P x equals 0.5 (x0^2 + x1^2).
    fn simplex_spec() -> QpSpec {
        // Constraints as CSC of the 3 by 2 matrix
        //   row 0: [1, 0]  (x0 box)
        //   row 1: [0, 1]  (x1 box)
        //   row 2: [1, 1]  (sum)
        // Column 0 has rows {0, 2}, column 1 has rows {1, 2}.
        QpSpec {
            n: 2,
            m: 3,
            p: PMat::Diagonal(vec![1.0, 1.0]),
            q: vec![0.0, 0.0],
            a_indptr: vec![0, 2, 4],
            a_indices: vec![0, 2, 1, 2],
            a_values: vec![1.0, 1.0, 1.0, 1.0],
            l: vec![0.0, 0.0, 1.0],
            u: vec![f64::INFINITY, f64::INFINITY, 1.0],
            convexity: Convexity::Psd,
        }
    }

    #[test]
    fn solves_a_simplex_projection() {
        let spec = simplex_spec();
        let sol = Osqp
            .solve(&spec, &QpOptions::default(), &|| false)
            .expect("setup succeeds");
        assert!(sol.status.is_solved(), "status {:?}", sol.status);
        assert!((sol.x[0] - 0.5).abs() < 1e-5, "x0 = {}", sol.x[0]);
        assert!((sol.x[1] - 0.5).abs() < 1e-5, "x1 = {}", sol.x[1]);
        assert_eq!(sol.duals.len(), 3);
    }

    #[test]
    fn a_pending_interrupt_stops_the_solve() {
        let spec = simplex_spec();
        // A tiny chunk forces a boundary before convergence; the always-pending
        // interrupt then stops at the first boundary.
        let opts = QpOptions {
            chunk_iters: 1,
            ..QpOptions::default()
        };
        let sol = Osqp.solve(&spec, &opts, &|| true).expect("setup succeeds");
        assert!(sol.interrupted);
        assert_eq!(sol.status, QpStatus::Interrupted);
    }

    #[test]
    fn open_upper_bounds_expressed_as_infinity_are_accepted() {
        // The box upper bound is a literal infinity; the backend must bring it to
        // the finite sentinel so OSQP reads it as an open side rather than
        // rejecting the data.
        let spec = simplex_spec();
        assert!(spec.u.iter().any(|b| b.is_infinite()));
        let sol = Osqp
            .solve(&spec, &QpOptions::default(), &|| false)
            .expect("setup succeeds");
        assert!(sol.status.is_solved());
    }
}
