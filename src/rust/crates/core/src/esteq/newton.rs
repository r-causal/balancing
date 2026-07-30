//! Damped Newton method with ridge escalation and Armijo backtracking.

use faer::MatMut;

use super::{EsteqProblem, SolveOptions, SolveReport, Solver};
use crate::linalg::solve_symmetric_ridge;

/// Armijo sufficient-decrease constant.
const ARMIJO_C: f64 = 1e-4;
/// Maximum step halvings per iteration.
const MAX_BACKTRACKS: usize = 30;
/// Maximum ridge escalations before an iteration is abandoned.
const MAX_ESCALATIONS: usize = 12;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn sup_norm(a: &[f64]) -> f64 {
    a.iter().fold(0.0_f64, |m, x| m.max(x.abs()))
}

/// Symmetric matrix-vector product for a column-major `p` by `p` slice.
fn sym_matvec(h: &[f64], p: usize, x: &[f64], out: &mut [f64]) {
    for i in 0..p {
        let mut acc = 0.0;
        for j in 0..p {
            acc += h[j * p + i] * x[j];
        }
        out[i] = acc;
    }
}

/// A positive ridge seed scaled to the Hessian magnitude, used the first time a
/// zero ridge fails to yield a usable direction.
fn ridge_seed(h: &[f64], p: usize) -> f64 {
    let mut trace = 0.0;
    for i in 0..p {
        trace += h[i * p + i].abs();
    }
    let scale = (trace / p as f64).max(1.0);
    1e-10 * scale
}

/// Solve a smooth (or root-finding) estimating-equation problem by damped Newton.
pub fn solve<P: EsteqProblem>(
    problem: &P,
    beta: &mut [f64],
    opts: &SolveOptions,
    interrupt: &dyn Fn() -> bool,
) -> SolveReport {
    solve_with_min_iter(problem, beta, opts, 0, interrupt)
}

/// Polish an estimate that is already near the solution, guaranteeing at least
/// one full Newton step.
///
/// The hybrid solver warm-starts with L-BFGS to a loose tolerance, which may
/// leave the estimate already inside the gradient tolerance. Forcing a Newton
/// step drives the estimating-equation output to a machine-precision solution,
/// which is what the second-order method uniquely provides and what linearized
/// inference relies on.
pub fn solve_polish<P: EsteqProblem>(
    problem: &P,
    beta: &mut [f64],
    opts: &SolveOptions,
    interrupt: &dyn Fn() -> bool,
) -> SolveReport {
    solve_with_min_iter(problem, beta, opts, 1, interrupt)
}

/// Damped Newton that performs at least `min_iter` steps before the convergence
/// check can stop it.
fn solve_with_min_iter<P: EsteqProblem>(
    problem: &P,
    beta: &mut [f64],
    opts: &SolveOptions,
    min_iter: usize,
    interrupt: &dyn Fn() -> bool,
) -> SolveReport {
    let p = problem.n_params();
    let smooth = problem.value(beta).is_some();
    // The tolerance is read at the scale the problem's gradient carries, so a
    // problem whose estimating function is a sampling-weight total gets the same
    // verdict whatever units those weights are expressed in. The scale is a
    // constant of the problem, so it is resolved once here rather than at every
    // convergence check.
    let grad_tol = opts.grad_tol * problem.residual_scale();

    let mut g = vec![0.0; p];
    let mut h = vec![0.0; p * p];
    let mut dir = vec![0.0; p];
    let mut beta_trial = vec![0.0; p];
    let mut merit_grad = vec![0.0; p];
    let mut g_trial = vec![0.0; p];

    let mut ridge_max = 0.0_f64;
    let mut iterations = 0;
    let mut converged = false;
    let mut interrupted = false;

    for iter in 0..opts.max_iter {
        if interrupt() {
            interrupted = true;
            break;
        }

        let value = {
            let h_mat = MatMut::from_column_major_slice_mut(&mut h, p, p);
            problem.value_grad_hess(beta, &mut g, h_mat)
        };
        let grad_norm = sup_norm(&g);
        if iter >= min_iter && grad_norm <= grad_tol {
            converged = true;
            break;
        }

        // Newton direction: solve (H + ridge I) dir = -g, growing the ridge from
        // zero until the factorization succeeds and the direction descends.
        let mut ridge = 0.0;
        let mut have_dir = false;
        for escalation in 0..=MAX_ESCALATIONS {
            for k in 0..p {
                dir[k] = -g[k];
            }
            if solve_symmetric_ridge(&h, p, ridge, &mut dir) && dot(&g, &dir) < 0.0 {
                have_dir = true;
                break;
            }
            ridge = if ridge == 0.0 {
                ridge_seed(&h, p)
            } else {
                ridge * 10.0
            };
            ridge_max = ridge_max.max(ridge);
            if escalation == MAX_ESCALATIONS {
                break;
            }
        }
        if !have_dir {
            break;
        }

        // Directional derivative of the merit being reduced. For a smooth
        // objective it is g . dir; for a root-finding merit 0.5||g||^2 it is
        // (H g) . dir.
        let base = value.unwrap_or_else(|| 0.5 * dot(&g, &g));
        let dderiv = if smooth {
            dot(&g, &dir)
        } else {
            sym_matvec(&h, p, &g, &mut merit_grad);
            dot(&merit_grad, &dir)
        };

        let mut step = 1.0;
        let mut accepted = false;
        for _ in 0..MAX_BACKTRACKS {
            for k in 0..p {
                beta_trial[k] = beta[k] + step * dir[k];
            }
            let trial = if smooth {
                problem
                    .value(&beta_trial)
                    .expect("smooth problem returns a value")
            } else {
                problem.gradient(&beta_trial, &mut g_trial);
                0.5 * dot(&g_trial, &g_trial)
            };
            if trial <= base + ARMIJO_C * step * dderiv {
                accepted = true;
                break;
            }
            step *= 0.5;
        }
        // A descent step that no longer reduces the objective means the
        // predicted decrease has fallen below the objective's floating-point
        // resolution: the iterate is at the numerical optimum, so stop rather
        // than crawl in ever-smaller steps.
        if !accepted {
            break;
        }
        for k in 0..p {
            beta[k] += step * dir[k];
        }
        // The count follows the accepted steps, so it states how many the returned
        // parameters embody. Recording the pass index at the top of the loop
        // instead would report one fewer than were taken when the cap is
        // exhausted, and `iterations == max_iter` is how a caller recognizes a fit
        // that ran out of its budget. Every other exit leaves the loop before
        // moving the iterate, so those counts are unchanged.
        iterations = iter + 1;
    }

    // Report the gradient and objective at the returned parameters. The reported
    // norm stays on the problem's own scale, the scale its estimating functions
    // are stored at; only the tolerance it is judged against moves.
    problem.gradient(beta, &mut g);
    let grad_norm = sup_norm(&g);
    if grad_norm <= grad_tol {
        converged = true;
    }
    let final_value = problem.value(beta).unwrap_or_else(|| 0.5 * dot(&g, &g));

    SolveReport {
        converged,
        interrupted,
        iterations,
        grad_norm,
        final_value,
        ridge_max,
        solver_used: Solver::Newton,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::esteq::EsteqProblem;
    use faer::prelude::ReborrowMut;

    /// `0.5 (beta - center)^2` per coordinate: one Newton step reaches the
    /// minimizer exactly, so the number of accepted steps is known in advance.
    struct Quadratic {
        center: Vec<f64>,
    }

    impl EsteqProblem for Quadratic {
        fn n_params(&self) -> usize {
            self.center.len()
        }
        fn n_units(&self) -> usize {
            self.center.len()
        }
        fn value(&self, beta: &[f64]) -> Option<f64> {
            Some(
                beta.iter()
                    .zip(&self.center)
                    .map(|(b, c)| 0.5 * (b - c) * (b - c))
                    .sum(),
            )
        }
        fn gradient(&self, beta: &[f64], g: &mut [f64]) {
            for (j, gj) in g.iter_mut().enumerate() {
                *gj = beta[j] - self.center[j];
            }
        }
        fn hessian(&self, _beta: &[f64], mut h: MatMut<'_, f64>) {
            let p = self.center.len();
            for i in 0..p {
                for j in 0..p {
                    *h.rb_mut().get_mut(i, j) = f64::from(i == j);
                }
            }
        }
        fn psi(&self, _beta: &[f64], _out: MatMut<'_, f64>) {}
    }

    /// `0.25 beta^4`, whose Newton step is `beta -> (2/3) beta`. The iteration
    /// converges only linearly, so a small cap is always exhausted.
    struct Quartic;

    impl EsteqProblem for Quartic {
        fn n_params(&self) -> usize {
            1
        }
        fn n_units(&self) -> usize {
            1
        }
        fn value(&self, beta: &[f64]) -> Option<f64> {
            Some(0.25 * beta[0].powi(4))
        }
        fn gradient(&self, beta: &[f64], g: &mut [f64]) {
            g[0] = beta[0].powi(3);
        }
        fn hessian(&self, beta: &[f64], mut h: MatMut<'_, f64>) {
            *h.rb_mut().get_mut(0, 0) = 3.0 * beta[0] * beta[0];
        }
        fn psi(&self, _beta: &[f64], _out: MatMut<'_, f64>) {}
    }

    fn opts(max_iter: usize) -> SolveOptions {
        SolveOptions {
            max_iter,
            grad_tol: 1e-12,
            fista_rel_tol: 1e-12,
        }
    }

    #[test]
    fn an_exhausted_cap_reports_every_accepted_step() {
        // Three passes of the loop take three accepted steps and the quartic is
        // nowhere near the gradient tolerance afterwards, so the report must read
        // three: `iterations == max_iter` is how a caller recognizes a fit that
        // ran out of budget.
        let mut beta = vec![1.0];
        let report = solve(&Quartic, &mut beta, &opts(3), &|| false);
        assert!(
            !report.converged,
            "the quartic cannot converge in three steps"
        );
        assert_eq!(report.iterations, 3);
        // Each step multiplies the iterate by 2/3, so three steps land on (2/3)^3.
        assert!(
            (beta[0] - (2.0_f64 / 3.0).powi(3)).abs() < 1e-12,
            "beta {} after three Newton steps",
            beta[0]
        );
    }

    #[test]
    fn a_converged_solve_counts_only_the_steps_it_took() {
        // One Newton step solves a quadratic exactly and the next pass sees a zero
        // gradient, so the count is one. The convergence path is exact today and
        // must stay exact.
        let problem = Quadratic {
            center: vec![1.5, -2.0],
        };
        let mut beta = vec![0.0, 0.0];
        let report = solve(&problem, &mut beta, &opts(50), &|| false);
        assert!(report.converged);
        assert_eq!(report.iterations, 1);
    }

    #[test]
    fn an_interrupt_before_the_first_step_reports_no_iterations() {
        let problem = Quadratic {
            center: vec![1.5, -2.0],
        };
        let mut beta = vec![0.0, 0.0];
        let report = solve(&problem, &mut beta, &opts(50), &|| true);
        assert!(report.interrupted);
        assert_eq!(report.iterations, 0);
        assert_eq!(beta, vec![0.0, 0.0]);
    }
}
