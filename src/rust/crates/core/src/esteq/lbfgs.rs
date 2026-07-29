//! Limited-memory BFGS adapter over the basin backend.
//!
//! This exposes the same call shape as [`super::newton::solve`] so the solver is
//! a runtime choice: quasi-Newton and Newton can be compared on the same problem
//! without a rebuild. L-BFGS uses only the objective and its gradient, so it
//! applies to the smooth members of the estimating-equation family.

use basin::core::problem::{CostFunction, Gradient};
use basin::solver::lbfgs::Unbounded;
use basin::{Executor, GradientTolerance, Lbfgs, LbfgsState, MaxIter, State};

use super::{EsteqProblem, SolveOptions, SolveReport, Solver};

/// History length for the L-BFGS two-loop recursion. Within the recommended
/// range for problems whose parameter count is small.
const HISTORY: usize = 10;

/// Wraps an [`EsteqProblem`] as a basin cost-and-gradient problem.
struct Adapter<'a, P: EsteqProblem> {
    inner: &'a P,
}

impl<P: EsteqProblem> CostFunction for Adapter<'_, P> {
    type Param = Vec<f64>;
    type Output = f64;
    type Error = std::convert::Infallible;

    fn cost(&self, param: &Vec<f64>) -> Result<f64, Self::Error> {
        Ok(self
            .inner
            .value(param)
            .expect("L-BFGS requires a smooth objective"))
    }
}

impl<P: EsteqProblem> Gradient for Adapter<'_, P> {
    type Gradient = Vec<f64>;

    fn gradient(&self, param: &Vec<f64>) -> Result<Vec<f64>, Self::Error> {
        let mut g = vec![0.0; self.inner.n_params()];
        self.inner.gradient(param, &mut g);
        Ok(g)
    }
}

/// Minimize a smooth estimating-equation problem by L-BFGS, writing the solution
/// back into `beta`.
pub fn solve<P: EsteqProblem>(
    problem: &P,
    beta: &mut [f64],
    opts: &SolveOptions,
    interrupt: &dyn Fn() -> bool,
) -> SolveReport {
    let adapter = Adapter { inner: problem };
    let state = LbfgsState::new(beta.to_vec(), HISTORY);
    let solver = Lbfgs::<Unbounded>::new();
    // The tolerance is read at the scale the problem's gradient carries, the
    // same convention the Newton method applies, so the solver choice never
    // moves the verdict.
    let grad_tol = opts.grad_tol * problem.residual_scale();

    let result = Executor::new(adapter, solver, state)
        .terminate_on(MaxIter(opts.max_iter as u64))
        .terminate_on(GradientTolerance(grad_tol))
        .run()
        .expect("basin L-BFGS is infallible for an infallible problem");

    let solution = result.param();
    beta.copy_from_slice(solution);
    let iterations = result.state.iter() as usize;

    let mut g = vec![0.0; problem.n_params()];
    problem.gradient(beta, &mut g);
    let grad_norm = g.iter().fold(0.0_f64, |m, gi| m.max(gi.abs()));
    let final_value = problem
        .value(beta)
        .expect("L-BFGS requires a smooth objective");

    // basin runs to completion without polling, so a pending interrupt is
    // observed here, once control returns to this thread.
    let interrupted = interrupt();

    // Convergence is judged on the recomputed sup norm against the tolerance, not
    // on basin's exit reason. basin terminates on the L2 gradient norm, which is
    // never smaller than the sup norm, so a run can hit the iteration cap with its
    // L2 norm above the tolerance while the sup norm is already below it; that run
    // has converged by this criterion.
    SolveReport {
        converged: grad_norm <= grad_tol,
        interrupted,
        iterations,
        grad_norm,
        final_value,
        ridge_max: 0.0,
        solver_used: Solver::Lbfgs,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::MatMut;
    use faer::prelude::ReborrowMut;

    // A separable quadratic with a known minimizer, wrapped as an EsteqProblem.
    struct Quad {
        center: Vec<f64>,
        scale: Vec<f64>,
    }

    impl EsteqProblem for Quad {
        fn n_params(&self) -> usize {
            self.center.len()
        }
        fn n_units(&self) -> usize {
            self.center.len()
        }
        fn value(&self, beta: &[f64]) -> Option<f64> {
            let mut v = 0.0;
            for (j, &b) in beta.iter().enumerate() {
                v += 0.5 * self.scale[j] * (b - self.center[j]).powi(2);
            }
            Some(v)
        }
        fn gradient(&self, beta: &[f64], g: &mut [f64]) {
            for (j, gj) in g.iter_mut().enumerate() {
                *gj = self.scale[j] * (beta[j] - self.center[j]);
            }
        }
        fn hessian(&self, _beta: &[f64], mut h: MatMut<'_, f64>) {
            let p = self.center.len();
            for i in 0..p {
                for j in 0..p {
                    *h.rb_mut().get_mut(i, j) = if i == j { self.scale[j] } else { 0.0 };
                }
            }
        }
        fn psi(&self, _beta: &[f64], _out: MatMut<'_, f64>) {}
    }

    #[test]
    fn lbfgs_finds_the_quadratic_minimizer() {
        let problem = Quad {
            center: vec![1.5, -2.0, 3.25],
            scale: vec![1.0, 4.0, 0.5],
        };
        let mut beta = vec![0.0, 0.0, 0.0];
        let report = solve(
            &problem,
            &mut beta,
            &SolveOptions {
                max_iter: 500,
                grad_tol: 1e-10,
                fista_rel_tol: 1e-10,
            },
            &|| false,
        );
        assert!(report.converged);
        for (b, c) in beta.iter().zip(&problem.center) {
            assert!((b - c).abs() < 1e-6);
        }
    }

    // Regression: convergence is judged on the final gradient sup norm against the
    // tolerance, independent of how basin's loop exited. basin's own gradient
    // termination uses the L2 norm, so a run can reach the iteration cap with its
    // L2 norm still above the tolerance while the sup norm is already below it; the
    // adapter must report that run converged. The problem below places the start
    // point at a gradient of 0.6e-8 in every one of four components: the sup norm
    // is 0.6e-8 (below the 1e-8 tolerance) while the L2 norm is 1.2e-8 (above it),
    // so basin exits on the iteration cap rather than its gradient criterion.
    #[test]
    fn a_gradient_below_tolerance_at_the_iteration_cap_reports_converged() {
        let grad_tol = 1e-8;
        let component = 0.6e-8;
        let problem = Quad {
            // gradient(0) = scale * (0 - center) = -center, so center = -component
            // puts every start-point gradient component at +component.
            center: vec![-component; 4],
            scale: vec![1.0; 4],
        };
        let mut beta = vec![0.0; 4];
        let report = solve(
            &problem,
            &mut beta,
            &SolveOptions {
                // A zero cap forces the run to exit on the iteration limit at the
                // start point, the deterministic form of "cap reached".
                max_iter: 0,
                grad_tol,
                fista_rel_tol: 1e-10,
            },
            &|| false,
        );

        // The exit was the cap, not basin's gradient criterion: the L2 norm at the
        // start point exceeds the tolerance.
        let l2 = (4.0 * component * component).sqrt();
        assert!(
            l2 > grad_tol,
            "the scenario needs the L2 norm above tolerance"
        );
        assert_eq!(report.iterations, 0, "the iteration cap should bind");
        assert!(
            report.grad_norm <= grad_tol,
            "sup norm {} should be below tolerance",
            report.grad_norm
        );
        assert!(
            report.converged,
            "a sup-norm gradient below tolerance is convergence regardless of the exit reason"
        );
    }
}
