//! Estimating-equation solvers.
//!
//! A problem whose solution is a root of a system of estimating equations, or
//! equivalently the minimizer of a smooth convex objective, implements
//! [`EsteqProblem`]. The damped Newton method in [`newton`] is the default and
//! the fastest path for the strictly convex members of the family; [`lbfgs`]
//! offers a quasi-Newton alternative behind the `solver-basin` feature, and
//! [`fista`] handles the composite (smooth plus L1) problems the Newton method
//! cannot address directly.

pub mod fista;
#[cfg(feature = "solver-basin")]
pub mod lbfgs;
pub mod newton;

use faer::MatMut;

/// A smooth estimating-equation problem in `n_params` parameters over `n_units`
/// observations.
///
/// The four evaluation methods match the estimating-equation contract exactly.
/// [`EsteqProblem::value_grad_hess`] is the fused hot path the Newton method
/// calls: it computes the per-unit tilt once and accumulates the gradient and
/// Hessian in a single pass. Its default delegates to the individual methods, so
/// a problem may implement only those and still solve, at the cost of the extra
/// pass.
pub trait EsteqProblem {
    /// Number of parameters, the dimension of `beta`.
    fn n_params(&self) -> usize;

    /// Number of observations contributing to the estimating equations.
    fn n_units(&self) -> usize;

    /// The smooth objective value when one exists (the entropy dual, a GMM
    /// norm). `None` marks a pure root-finding problem, for which the Newton
    /// method minimizes the merit `0.5 * ||g||^2` instead.
    fn value(&self, beta: &[f64]) -> Option<f64>;

    /// Write the gradient of the objective (equivalently the summed estimating
    /// function) at `beta` into `g`.
    fn gradient(&self, beta: &[f64], g: &mut [f64]);

    /// Write the Hessian (the Jacobian of the summed estimating function) at
    /// `beta` into the `n_params` by `n_params` matrix `h`.
    fn hessian(&self, beta: &[f64], h: MatMut<'_, f64>);

    /// Write the `n_units` by `n_params` matrix of per-unit estimating
    /// functions at `beta` into `out`.
    fn psi(&self, beta: &[f64], out: MatMut<'_, f64>);

    /// Fused evaluation of the value, gradient, and Hessian at `beta`.
    ///
    /// Returns the same value as [`EsteqProblem::value`]. The default calls the
    /// individual methods; problems override it to compute the shared per-unit
    /// quantities once.
    fn value_grad_hess(&self, beta: &[f64], g: &mut [f64], h: MatMut<'_, f64>) -> Option<f64> {
        self.gradient(beta, g);
        self.hessian(beta, h);
        self.value(beta)
    }
}

/// Which solver drives the estimating-equation iteration.
///
/// The choice is a runtime option so quasi-Newton and Newton can be compared
/// without rebuilding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Solver {
    /// Damped Newton with ridge escalation and Armijo backtracking.
    Newton,
    /// Limited-memory BFGS through the basin backend.
    Lbfgs,
    /// L-BFGS to a loose tolerance, then Newton to polish.
    LbfgsThenNewton,
}

/// Loosest gradient sup norm the hybrid's L-BFGS warm start is asked to reach.
/// The warm start hands off to the Newton polish here rather than pursuing the
/// full tolerance itself, where the Hessian-based step is far more efficient.
const HYBRID_WARM_GRAD_TOL: f64 = 1e-6;

/// Tuning shared across the estimating-equation solvers.
#[derive(Debug, Clone, Copy)]
pub struct SolveOptions {
    /// Maximum outer iterations.
    pub max_iter: usize,
    /// Convergence threshold on the sup norm of the gradient.
    pub grad_tol: f64,
    /// Relative change in the full loss at which FISTA stops.
    pub fista_rel_tol: f64,
}

impl Default for SolveOptions {
    fn default() -> Self {
        Self {
            max_iter: 100,
            grad_tol: 1e-12,
            fista_rel_tol: 1e-10,
        }
    }
}

/// Outcome of a solve.
#[derive(Debug, Clone, Copy)]
pub struct SolveReport {
    /// Whether the gradient tolerance was met.
    pub converged: bool,
    /// Whether the solve stopped because a user interrupt was pending.
    pub interrupted: bool,
    /// Outer iterations performed.
    pub iterations: usize,
    /// Sup norm of the gradient at the returned parameters.
    pub grad_norm: f64,
    /// Objective value at the returned parameters, or the merit for a
    /// root-finding problem.
    pub final_value: f64,
    /// Largest ridge added to the Hessian across all iterations.
    pub ridge_max: f64,
    /// The solver that produced the result.
    pub solver_used: Solver,
}

/// Solve `problem` in place, starting from `beta`, with the requested solver.
///
/// `interrupt` is polled between iterations; returning `true` stops the solve
/// and reports `converged = false`.
pub fn solve<P: EsteqProblem>(
    problem: &P,
    beta: &mut [f64],
    solver: Solver,
    opts: &SolveOptions,
    interrupt: &dyn Fn() -> bool,
) -> SolveReport {
    match solver {
        Solver::Newton => newton::solve(problem, beta, opts, interrupt),
        Solver::Lbfgs => lbfgs_solve(problem, beta, opts, interrupt),
        Solver::LbfgsThenNewton => {
            // Warm-start with L-BFGS to a loose tolerance. Its gradient-only
            // iterations avoid the O(n p^2) Hessian, so a coarse solution is
            // cheap; the square root of the target tolerance, floored so the
            // warm start never chases machine precision itself, hands a good
            // point to the polish.
            let warm_opts = SolveOptions {
                grad_tol: opts.grad_tol.sqrt().max(HYBRID_WARM_GRAD_TOL),
                ..*opts
            };
            let warm = lbfgs_solve(problem, beta, &warm_opts, interrupt);
            if warm.interrupted || interrupt() {
                // The polish never runs, so the estimate sits at the loose warm
                // tolerance. Report the interrupt so the caller re-signals it and
                // do not claim convergence at that point, regardless of what the
                // warm start recorded before the interrupt was observed.
                return SolveReport {
                    solver_used: Solver::LbfgsThenNewton,
                    interrupted: true,
                    converged: false,
                    ..warm
                };
            }
            // Polish with Newton to the full tolerance, taking at least one
            // Newton step so the estimating-equation output is evaluated at a
            // machine-precision solution rather than the loose warm-start point.
            let polished = newton::solve_polish(problem, beta, opts, interrupt);
            SolveReport {
                solver_used: Solver::LbfgsThenNewton,
                iterations: warm.iterations + polished.iterations,
                ..polished
            }
        }
    }
}

#[cfg(feature = "solver-basin")]
fn lbfgs_solve<P: EsteqProblem>(
    problem: &P,
    beta: &mut [f64],
    opts: &SolveOptions,
    interrupt: &dyn Fn() -> bool,
) -> SolveReport {
    lbfgs::solve(problem, beta, opts, interrupt)
}

// Without the basin backend the L-BFGS path falls back to Newton, so a runtime
// request for it still solves rather than failing.
#[cfg(not(feature = "solver-basin"))]
fn lbfgs_solve<P: EsteqProblem>(
    problem: &P,
    beta: &mut [f64],
    opts: &SolveOptions,
    interrupt: &dyn Fn() -> bool,
) -> SolveReport {
    newton::solve(problem, beta, opts, interrupt)
}
