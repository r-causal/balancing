//! The OSQP backend: an ADMM solver that is the default for every method.
//!
//! ADMM regularizes the KKT system it factors, so it solves the indefinite
//! quadratic form of energy balancing where an interior-point method fails.
//! OSQP has no mid-solve callback, so the interrupt is honored by solving in
//! bounded iteration chunks and warm starting between them: after each chunk the
//! interrupt closure is polled, and the primal and dual iterates carry into the
//! next chunk so the split does not change the trajectory.

use osqp::{CscMatrix, Problem, Settings, Solution, Status};

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
        // reaching it as MaxIterationsReached or, when the looser tolerances are
        // already met, SolvedInaccurate. Both are chunk-boundary outcomes, not
        // genuine terminations, so the loop distinguishes them from the true cap
        // by tracking the total against the budget, and trims the last chunk to
        // what is left of it so a cap that is not a whole number of chunks is
        // still honored exactly. OSQP will not run fewer than one iteration, so a
        // zero cap is served as one.
        let budget = opts.max_iter.max(1);
        let chunk = opts.chunk_iters.clamp(1, budget) as u32;
        let mut settings = Settings::default()
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
        let mut cap = chunk;
        loop {
            let status = problem.solve();
            total_iters += status.iter() as usize;

            // A chunk-boundary outcome (the iteration cap reached, whether the
            // looser tolerances were met or not) is non-terminal: the solve carries
            // on into the next chunk until the total reaches the budget. Every
            // other status is a genuine termination.
            let terminal = !matches!(
                status,
                Status::MaxIterationsReached(_) | Status::SolvedInaccurate(_)
            );

            if terminal || total_iters >= budget {
                return Ok(finish(spec, &status, total_iters, false));
            }

            // A chunk boundary under the cap: poll the interrupt, then warm start
            // the next chunk from the current iterate.
            if interrupt() {
                return Ok(finish(spec, &status, total_iters, true));
            }
            let (x, y) = match extract_solution(&status) {
                Some(sol) => sol,
                None => return Ok(finish(spec, &status, total_iters, false)),
            };
            problem.warm_start(&x, &y);

            // Trim the next chunk to the remaining budget. The comparison keeps
            // the settings update to the one chunk that needs it, since OSQP
            // revalidates the whole settings block on every call.
            let remaining = u32::try_from(budget - total_iters).unwrap_or(u32::MAX);
            let next_cap = chunk.min(remaining);
            if next_cap != cap {
                settings = settings.max_iter(next_cap);
                problem.update_settings(&settings);
                cap = next_cap;
            }
        }
    }
}

/// The OSQP solution a status carries, when it carries one.
///
/// OSQP attaches a `Solution` to `Solved`, `SolvedInaccurate`,
/// `MaxIterationsReached`, and `TimeLimitReached`; the infeasibility and
/// non-convex certificates carry none. The crate's own `Status::solution()` is
/// `Some` only for `Solved`, so a solve that stopped at the iteration cap with a
/// usable iterate would otherwise read as all zeros.
fn solution_of<'a>(status: &Status<'a>) -> Option<Solution<'a>> {
    match status {
        Status::Solved(s)
        | Status::SolvedInaccurate(s)
        | Status::MaxIterationsReached(s)
        | Status::TimeLimitReached(s) => Some(s.clone()),
        _ => None,
    }
}

/// Copy the primal and dual iterates out of a status that carries a solution.
/// Copying releases the borrow of the problem, which the warm start of the next
/// chunk needs mutably.
fn extract_solution(status: &Status<'_>) -> Option<(Vec<f64>, Vec<f64>)> {
    solution_of(status).map(|s| (s.x().to_vec(), s.y().to_vec()))
}

/// Build a [`QpSolution`] from an OSQP status.
///
/// The residuals are the ones OSQP measured at the iterate it returns, read off
/// the solution it attaches to every status that carries one, which includes the
/// iteration-cap and time-limit outcomes an interrupted or capped solve ends on.
/// A status that carries no iterate at all (an infeasibility or non-convex
/// certificate) has no residual to report either, so the primal, the duals, and
/// both residuals are zeros together and the caller reacts to the status rather
/// than to a number describing nothing.
fn finish(spec: &QpSpec, status: &Status<'_>, iterations: usize, interrupted: bool) -> QpSolution {
    let n = spec.n;
    let m = spec.m;
    let (x, duals, pri_res, dua_res) = match solution_of(status) {
        Some(s) => (s.x().to_vec(), s.y().to_vec(), s.pri_res(), s.dua_res()),
        None => (vec![0.0; n], vec![0.0; m], 0.0, 0.0),
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
        pri_res,
        dua_res,
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
    fn a_multi_chunk_solve_matches_a_single_chunk_solve() {
        // A chunk smaller than the iterations the problem needs forces the loop to
        // warm start across several chunks. The result must match the single-chunk
        // solve, proving the warm start carries the iterate forward rather than
        // restarting or bailing out with zeros at the first cap. The chunk stays at
        // or above OSQP's rho-adaptation interval so the step size still adapts
        // within each chunk, the regime real solves run in.
        let spec = simplex_spec();
        let chunked = Osqp
            .solve(
                &spec,
                &QpOptions {
                    chunk_iters: 25,
                    ..QpOptions::default()
                },
                &|| false,
            )
            .expect("setup succeeds");
        let single = Osqp
            .solve(&spec, &QpOptions::default(), &|| false)
            .expect("setup succeeds");
        assert!(chunked.status.is_solved(), "status {:?}", chunked.status);
        assert!(
            chunked.iterations > 25,
            "solve took {} iterations, not a multi-chunk run",
            chunked.iterations
        );
        assert!((chunked.x[0] - single.x[0]).abs() < 1e-6, "x0 diverged");
        assert!((chunked.x[1] - single.x[1]).abs() < 1e-6, "x1 diverged");
        assert!((chunked.x[0] - 0.5).abs() < 1e-5, "x0 = {}", chunked.x[0]);
    }

    /// The sup norm of the constraint violation of `x`, the distance from `A x`
    /// to `[l, u]` row by row. OSQP measures its primal residual against a slack
    /// that lies inside the bounds, so the residual it reports is never smaller
    /// than this.
    fn sup_norm_violation(spec: &QpSpec, x: &[f64]) -> f64 {
        let mut ax = vec![0.0; spec.m];
        for (j, &xj) in x.iter().enumerate().take(spec.n) {
            for k in spec.a_indptr[j]..spec.a_indptr[j + 1] {
                ax[spec.a_indices[k]] += spec.a_values[k] * xj;
            }
        }
        (0..spec.m)
            .map(|r| (spec.l[r] - ax[r]).max(ax[r] - spec.u[r]).max(0.0))
            .fold(0.0_f64, f64::max)
    }

    #[test]
    fn a_capped_solve_reports_the_backend_residuals() {
        // The solution contract documents both residuals as backend-reported. A
        // single ADMM iteration from the origin is nowhere near the simplex, so a
        // one-iteration cap must report a primal residual that is genuinely
        // nonzero and no smaller than the violation its own iterate carries.
        let spec = simplex_spec();
        let opts = QpOptions {
            max_iter: 1,
            chunk_iters: 1,
            polish: false,
            ..QpOptions::default()
        };
        let capped = Osqp.solve(&spec, &opts, &|| false).expect("setup succeeds");
        assert_eq!(capped.status, QpStatus::MaxIter);
        assert!(
            capped.pri_res.is_finite() && capped.pri_res > 0.0,
            "pri_res {} at the iteration cap",
            capped.pri_res
        );
        assert!(capped.dua_res.is_finite(), "dua_res {}", capped.dua_res);
        let violation = sup_norm_violation(&spec, &capped.x);
        assert!(
            capped.pri_res >= violation - 1e-12,
            "pri_res {} is below the iterate's own violation {violation}",
            capped.pri_res
        );

        // The other side of the same evidence: a converged solve reports a
        // residual at its tolerance rather than a constant.
        let solved = Osqp
            .solve(&spec, &QpOptions::default(), &|| false)
            .expect("setup succeeds");
        assert!(solved.status.is_solved(), "status {:?}", solved.status);
        assert!(
            solved.pri_res.is_finite() && solved.pri_res < 1e-6,
            "pri_res {} at convergence",
            solved.pri_res
        );
    }

    #[test]
    fn the_iteration_total_never_overruns_the_cap() {
        // A cap that is not a multiple of the chunk. Untrimmed chunks of three
        // would overshoot a cap of seven by two iterations, so the final chunk
        // must be trimmed to the remaining budget. The chunk sits below OSQP's
        // termination-check interval, so no chunk can report convergence and the
        // loop runs the cap out.
        let spec = simplex_spec();
        let opts = QpOptions {
            max_iter: 7,
            chunk_iters: 3,
            polish: false,
            ..QpOptions::default()
        };
        let sol = Osqp.solve(&spec, &opts, &|| false).expect("setup succeeds");
        assert_eq!(sol.status, QpStatus::MaxIter, "status {:?}", sol.status);
        assert!(
            sol.iterations <= opts.max_iter,
            "{} iterations against a cap of {}",
            sol.iterations,
            opts.max_iter
        );
        assert_eq!(
            sol.iterations, opts.max_iter,
            "an exhausted cap should report exactly the cap"
        );
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
