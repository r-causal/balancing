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

/// The smallest change in an objective of magnitude `value` that floating-point
/// arithmetic resolves: a predicted decrease at or below it cannot be
/// distinguished from no decrease at all.
///
/// Objectives smaller than one in magnitude share the resolution of one, the
/// same floor the FISTA relative-loss rule applies, so an objective that happens
/// to pass near zero does not report a resolution of zero. The floor also makes
/// this a lower bound on the true resolution of an accumulated objective, whose
/// rounding grows with the number of terms summed, which is the conservative
/// direction for anything that reads it as a stationarity certificate.
///
/// An objective that is not finite carries no resolution to speak of, and the
/// floor above would otherwise read a `NaN` as an objective of magnitude one:
/// `f64::max` returns its finite argument when the other is `NaN`. Reporting
/// zero instead refuses the certificate, which is the direction a caller reading
/// this as evidence of stationarity needs, and leaves the stall to be reported
/// as the non-convergence it is.
fn value_resolution(value: f64) -> f64 {
    if !value.is_finite() {
        return 0.0;
    }
    f64::EPSILON * value.abs().max(1.0)
}

/// A positive ridge seed scaled to the Hessian magnitude, used the first time a
/// zero ridge fails to yield a usable direction.
///
/// The seed is a small fraction of the Hessian's average diagonal magnitude, so
/// it is a relative perturbation on a well-scaled problem and grows with the
/// Hessian rather than against it. The floor of one on that magnitude is the
/// calibration assumption: for a Hessian whose diagonal is far below one, which a
/// sampling-weight scale far below one would produce, the seed stops tracking the
/// matrix and becomes the absolute `1e-10`, large enough to dominate the
/// curvature. The floor stands because the seed is only reachable when the
/// factorization of the unridged Hessian fails or its solution does not descend,
/// and the positive definite Jacobians the members of this family present do not
/// reach that branch; where they do, through an exactly singular direction, the
/// escalation from this seed is what recovers a usable direction rather than a
/// quantity any estimate depends on.
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
        //
        // Whether that numerical optimum also satisfies the gradient tolerance is
        // a separate question, and on a singular Hessian the answer can be no
        // however well the problem is solved. A dual with an exactly flat
        // direction, which the level indicators of a factor create by summing to
        // the constant function, loses the curvature that drives the gradient to
        // zero, and the gradient the solve can then reach is set by the
        // accumulated rounding of the weighted means rather than by the tolerance.
        // The decrease the full Newton step predicts, `-g . dir`, is the measure
        // of how far the objective still has to fall, and it is invariant to
        // reparameterization where the gradient sup norm is not. A predicted
        // decrease at or below the objective's own resolution therefore certifies
        // that the iterate minimizes the objective to the precision the
        // arithmetic carries, and convergence is claimed on that certificate
        // rather than on the stall alone: a solve that stalls with real progress
        // still available reports the non-convergence it should, as does one that
        // exhausts its iteration cap. A root-finding problem has no objective of
        // its own, only the merit built from the gradient it is trying to zero, so
        // there is no independent resolution to certify against and its stall
        // keeps the plain verdict.
        if !accepted {
            if smooth && -dderiv <= value_resolution(base) {
                converged = true;
            }
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

    // A stall is certified as convergence only when the predicted decrease sits
    // below the objective's resolution, so an objective that cannot be read must
    // report no resolution at all. Both non-finite readings would otherwise pass
    // a stall off as an optimum: `f64::max` hands a `NaN` back the finite floor
    // of one, and an infinite objective would report an infinite resolution that
    // every predicted decrease falls below.
    #[test]
    fn a_non_finite_objective_has_no_resolution() {
        assert_eq!(value_resolution(f64::NAN), 0.0);
        assert_eq!(value_resolution(f64::INFINITY), 0.0);
        assert_eq!(value_resolution(f64::NEG_INFINITY), 0.0);
        assert_eq!(value_resolution(0.0), f64::EPSILON);
        assert_eq!(value_resolution(-8.0), 8.0 * f64::EPSILON);
    }

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

    /// `0.5 beta^2` accumulated against a large offset, so the objective's
    /// floating-point resolution near the minimizer is `eps * offset` rather than
    /// `eps * value`: the decrease a Newton step predicts vanishes into the
    /// rounding of the sum while the gradient is still far above a tight
    /// tolerance. The gradient and Hessian are exact.
    ///
    /// This is the shape of the entropy dual at a singular Hessian. Its value is
    /// a log-sum over the sample, and once the exactly flat direction a factor's
    /// level indicators create has removed the curvature that would drive the
    /// gradient to zero, the remaining decrease sits orders of magnitude below
    /// that sum's rounding.
    struct OffsetQuadratic {
        offset: f64,
    }

    impl EsteqProblem for OffsetQuadratic {
        fn n_params(&self) -> usize {
            1
        }
        fn n_units(&self) -> usize {
            1
        }
        fn value(&self, beta: &[f64]) -> Option<f64> {
            Some((self.offset + 0.5 * beta[0] * beta[0]) - self.offset)
        }
        fn gradient(&self, beta: &[f64], g: &mut [f64]) {
            g[0] = beta[0];
        }
        fn hessian(&self, _beta: &[f64], mut h: MatMut<'_, f64>) {
            *h.rb_mut().get_mut(0, 0) = 1.0;
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
    fn a_stall_at_the_objectives_resolution_reports_converged() {
        // The start point sits where the predicted decrease, 1e-18, is below the
        // objective's floating-point resolution of one machine epsilon: every
        // trial step returns the same value, so the line search cannot accept
        // one. The iterate is at the numerical optimum, which is convergence,
        // even though the gradient of 1e-9 never reaches the 1e-12 tolerance.
        let problem = OffsetQuadratic { offset: 1.0 };
        let mut beta = vec![1e-9];
        let report = solve(&problem, &mut beta, &opts(50), &|| false);
        assert!(
            report.grad_norm > opts(50).grad_tol,
            "the scenario needs a gradient above the tolerance, got {}",
            report.grad_norm
        );
        assert_eq!(report.iterations, 0, "no step can be accepted");
        assert!(
            report.converged,
            "a stall with the predicted decrease below the objective's \
             resolution is the numerical optimum"
        );
    }

    #[test]
    fn a_stall_with_progress_still_available_does_not_claim_convergence() {
        // A far larger offset coarsens the objective's resolution to 1e-4, so the
        // line search stalls at a gradient of 1e-3 with a predicted decrease of
        // 1e-6. The certificate reads the resolution at the value's own
        // magnitude, so it does not fire, and the solve reports the
        // non-convergence it should: the stall alone is not evidence of an
        // optimum.
        let problem = OffsetQuadratic { offset: 1e12 };
        let mut beta = vec![1e-3];
        let report = solve(&problem, &mut beta, &opts(50), &|| false);
        assert_eq!(report.iterations, 0, "no step can be accepted");
        assert!(
            !report.converged,
            "a stall short of the resolution must not claim convergence"
        );
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
