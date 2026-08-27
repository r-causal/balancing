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
/// How many multiples of the objective's resolution an accepted step may buy and
/// still count as having bought nothing. See [`solve_with_min_iter`] for why the
/// factor is loose and why that costs nothing.
const UNPRODUCTIVE_DECREASE_FACTOR: f64 = 64.0;
/// How many consecutive unproductive accepted steps certify the iterate.
const UNPRODUCTIVE_STEPS: usize = 2;

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
    let mut unproductive = 0_usize;

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

        // The decrease the full step predicts is the measure of how far the
        // objective still has to fall. Once it sits at or below the objective's own
        // floating-point resolution, the objective has stopped being able to judge
        // a step at all: every trial the line search could offer returns a value
        // equal to the current one to within rounding, and the sufficient-decrease
        // test then decides between them by ulps. It asks for a ten thousandth of a
        // predicted decrease that is already beneath the resolution, so a trial one
        // ulp above the current value is refused, halved until some halving happens
        // to read at or below it, and that halving is accepted as a step which
        // moves nothing. Every pass repeats it and the solve spends its budget
        // standing still.
        //
        // What to do instead turns on how much accuracy the direction still has to
        // offer, and the two answers below are the same optimality reading applied
        // to problems that differ in that one respect.
        //
        // A direction built from a surrogate Hessian has none left: the surrogate's
        // progress is bounded by the same resolution, so a predicted decrease
        // beneath it means the iterate minimizes the objective to the precision the
        // arithmetic carries, and the iterations that follow only shuffle the
        // parameters inside the flat region. A generalized-method-of-moments
        // criterion left to run two thousand such steps moves in no printed digit.
        // That iterate is certified here, one step ahead of the stalled line search
        // below, and reading it earlier is what makes the certificate reachable: a
        // line search need not stall to be finished.
        //
        // An exact Hessian is the other case, and it used to be excluded from this
        // block outright, which left it the crawl rather than a verdict. Its steps
        // are superlinear and go on sharpening the parameters well after the value
        // has stopped registering the improvement, which is the accuracy the polish
        // exists to collect, so the step is worth taking. The objective cannot
        // referee it, so the line search is skipped and the full step is taken
        // unguarded, then judged by the one quantity that still carries information
        // at this scale: the gradient the step produced. A gradient that fell is
        // progress the objective could not see, so the step is kept and counted as
        // the accepted iteration it is, and the next pass's gradient test usually
        // certifies outright. A gradient that did not fall says the arithmetic has
        // nothing further to give, so the step is handed back, nothing is counted,
        // and the iterate is certified on the same optimality reading the surrogate
        // path uses. Judging by the gradient rather than by the value is what keeps
        // the exact carve-out from being an unguarded crawl: a rounding floor
        // cannot pass that test twice, so at most one probe is spent discovering
        // it.
        //
        // A gradient that is not finite is no evidence of progress however small
        // its sup norm reads, because `f64::max` hands back the finite argument of
        // a pair and would report a norm of zero for a direction that destroyed the
        // iterate. Such a step is handed back like any other that bought nothing.
        //
        // Two exclusions remain. A predicted decrease the objective can still
        // resolve goes to the line search as before, so genuine progress is
        // untouched. A root-finding problem has no objective of its own, only the
        // merit built from the gradient it is zeroing, so there is no independent
        // resolution to certify against. The minimum-iteration contract is honored
        // the way the gradient test honors it, so a polish still takes the step it
        // promises, through the line search, before any verdict.
        if smooth && iter >= min_iter && -dderiv <= value_resolution(base) {
            if !problem.hessian_is_exact() {
                converged = true;
                break;
            }
            for k in 0..p {
                beta_trial[k] = beta[k] + dir[k];
            }
            problem.gradient(&beta_trial, &mut g_trial);
            let probe_norm = sup_norm(&g_trial);
            if g_trial.iter().all(|x| x.is_finite()) && probe_norm < grad_norm {
                beta.copy_from_slice(&beta_trial);
                iterations = iter + 1;
                continue;
            }
            converged = true;
            break;
        }

        let mut step = 1.0;
        let mut accepted = false;
        let mut accepted_value = base;
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
                accepted_value = trial;
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
        //
        // For a smooth problem past its minimum-iteration count the block above now
        // takes the unresolvable case, by either route, so what remains here is the
        // stall a promised polish step runs into and the stall that arrives with
        // real progress still on the table.
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

        // What the step actually bought, read after the fact, as against what the
        // model predicted it would buy before it was taken. Consecutive steps that
        // buy nothing the objective can tell apart from its own rounding certify
        // the iterate: the solve is grinding at the resolution floor, and the
        // passes that follow only shuffle the parameters inside it.
        //
        // This route exists because the predicted decrease cannot decide the
        // question on its own. At a numerical minimizer the predicted decrease and
        // the resolution it is compared against are the same order of magnitude,
        // so a criterion sitting a hair above the bar is refused and the identical
        // fit on another machine, or from another draw of the same design, is
        // certified. What a step achieved is the sturdier reading, and it is read
        // over consecutive steps because one unproductive step is a coincidence a
        // wandering iterate can produce while real progress is still available.
        // Two in a row is not: the measured stalls run four steps and longer, so
        // waiting for the second costs one pass and rules the coincidence out.
        //
        // The factor is loose on purpose, and it costs nothing because the two
        // populations are nine orders of magnitude apart. An unproductive step in
        // the runs measured buys between a half and fourteen times the resolution;
        // a productive one, in every fixture here, buys upwards of two billion
        // times it. Sixty-four sits between them, though nowhere near the middle:
        // it clears the unproductive steps by about a factor of four and stays some
        // seven decades below the productive ones. The narrow side is the one to
        // respect. Lowering the factor past fourteen would stop the measured stalls
        // counting as unproductive and lose the verdict this route exists to reach,
        // where raising it has decades of slack before it could swallow real
        // progress. The exact bar of the predicted-decrease route, by contrast,
        // sits at the top edge of the first population with no room at all.
        //
        // Only accepted steps are evidence, which is what keeps this away from the
        // stall. A rejected trial step says the objective declined to move in the
        // direction offered, and an objective too coarse to register its own
        // improvement declines in exactly the same way as one already at its
        // optimum. A stall therefore keeps the predicted-decrease reading it always
        // had, and a solve that stalls before taking a single step, which is what
        // an unreadable objective does, never reaches this route at all. The
        // exclusions the predicted-decrease certificate carries apply here for the
        // same reasons: a root-finding merit has no objective of its own to resolve
        // against, an exact Hessian goes on sharpening the parameters after the
        // value has stopped moving, and the polish takes the step it promises
        // before any verdict.
        if smooth && iter >= min_iter && !problem.hessian_is_exact() {
            let bought = (base - accepted_value).abs();
            if bought <= UNPRODUCTIVE_DECREASE_FACTOR * value_resolution(base) {
                unproductive += 1;
            } else {
                unproductive = 0;
            }
            if unproductive >= UNPRODUCTIVE_STEPS {
                converged = true;
                break;
            }
        }
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
    ///
    /// `surrogate` changes only the answer given to
    /// [`EsteqProblem::hessian_is_exact`], the declaration the solver reads, so
    /// the two configurations isolate that branch and nothing else.
    struct OffsetQuadratic {
        offset: f64,
        surrogate: bool,
    }

    impl EsteqProblem for OffsetQuadratic {
        fn n_params(&self) -> usize {
            1
        }
        fn n_units(&self) -> usize {
            1
        }
        fn hessian_is_exact(&self) -> bool {
            !self.surrogate
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

    /// `offset + 0.5 beta^2` reported at its accumulated magnitude rather than
    /// relative to the offset, so near the minimizer the objective's resolution is
    /// `eps * offset`: the decrease a full step predicts disappears into the
    /// rounding of the sum while the gradient is still far above a tight tolerance,
    /// and the line search accepts the step anyway, because an objective that does
    /// not move satisfies a sufficient-decrease test asking for a ten thousandth of
    /// nothing. That combination is the micro-crawl, and it is what separates this
    /// from [`OffsetQuadratic`], whose line search stalls outright.
    ///
    /// The gradient and Hessian are exact in both configurations. `surrogate`
    /// changes only the answer given to `EsteqProblem::hessian_is_exact`, the
    /// declaration the solver reads, so the two configurations isolate that branch
    /// and nothing else about the problem.
    struct FlatQuadratic {
        offset: f64,
        surrogate: bool,
    }

    impl EsteqProblem for FlatQuadratic {
        fn n_params(&self) -> usize {
            1
        }
        fn n_units(&self) -> usize {
            1
        }
        fn hessian_is_exact(&self) -> bool {
            !self.surrogate
        }
        fn value(&self, beta: &[f64]) -> Option<f64> {
            Some(self.offset + 0.5 * beta[0] * beta[0])
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
        // objective's floating-point resolution of one machine epsilon, and the
        // direction comes from a surrogate whose progress that resolution bounds.
        // Neither route can move this iterate: the surrogate has no accuracy left
        // to collect, and the line search that would otherwise run cannot accept a
        // step either, since every trial returns the same value. The iterate is at
        // the numerical optimum, which is convergence, even though the gradient of
        // 1e-9 never reaches the 1e-12 tolerance.
        //
        // The exact-Hessian reading of the same objective is a different verdict
        // and belongs to
        // `an_exact_hessian_keeps_stepping_where_the_objective_has_flattened`.
        let problem = OffsetQuadratic {
            offset: 1.0,
            surrogate: true,
        };
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
        let problem = OffsetQuadratic {
            offset: 1e12,
            surrogate: false,
        };
        let mut beta = vec![1e-3];
        let report = solve(&problem, &mut beta, &opts(50), &|| false);
        assert_eq!(report.iterations, 0, "no step can be accepted");
        assert!(
            !report.converged,
            "a stall short of the resolution must not claim convergence"
        );
    }

    #[test]
    fn a_surrogate_hessian_certifies_an_unresolvable_predicted_decrease() {
        // The predicted decrease of 1e-18 is below the objective's resolution of
        // one machine epsilon, and the direction comes from a surrogate whose
        // progress that resolution bounds, so the iterate is the optimum the method
        // can reach. The verdict is read before the line search, which is the point:
        // this line search would have accepted the step and crawled to the cap.
        let problem = FlatQuadratic {
            offset: 1.0,
            surrogate: true,
        };
        let mut beta = vec![1e-9];
        let report = solve(&problem, &mut beta, &opts(50), &|| false);
        assert!(
            report.grad_norm > opts(50).grad_tol,
            "the scenario needs a gradient above the tolerance, got {}",
            report.grad_norm
        );
        assert_eq!(report.iterations, 0, "the certificate precedes any step");
        assert!(report.converged);
    }

    #[test]
    fn an_exact_hessian_keeps_stepping_where_the_objective_has_flattened() {
        // The identical objective and start point, differing only in the answer to
        // `hessian_is_exact`. A superlinear step goes on sharpening the parameters
        // after the value has stopped registering the improvement, so this one is
        // taken rather than certified away, and it lands on the minimizer exactly.
        let problem = FlatQuadratic {
            offset: 1.0,
            surrogate: false,
        };
        let mut beta = vec![1e-9];
        let report = solve(&problem, &mut beta, &opts(50), &|| false);
        assert!(report.converged);
        assert_eq!(report.iterations, 1, "the exact step is worth taking");
        assert_eq!(beta[0], 0.0, "and it reaches the minimizer");
    }

    #[test]
    fn a_resolvable_predicted_decrease_is_taken_not_certified() {
        // A surrogate Hessian on the same objective, started where the predicted
        // decrease of 1e-6 sits far above the resolution. The certificate must not
        // fire on progress the objective can still measure, so the step is taken
        // and the solve converges on the gradient as it always did.
        let problem = FlatQuadratic {
            offset: 1.0,
            surrogate: true,
        };
        let mut beta = vec![1e-3];
        let report = solve(&problem, &mut beta, &opts(50), &|| false);
        assert!(report.converged);
        assert_eq!(report.iterations, 1);
        assert_eq!(beta[0], 0.0);
    }

    #[test]
    fn a_polish_steps_before_the_surrogate_certificate_can_fire() {
        // The configuration `solve` certifies without moving, run through the
        // polish instead. The polish exists to collect the machine-precision
        // accuracy a warm start leaves on the table, so the step it promises has to
        // outrank every convergence verdict, the certificate included. Both runs
        // start from the same iterate and differ only in the minimum-iteration
        // contract, which is what makes the step count the whole difference.
        let problem = FlatQuadratic {
            offset: 1.0,
            surrogate: true,
        };
        let mut certified_beta = vec![1e-9];
        let certified = solve(&problem, &mut certified_beta, &opts(50), &|| false);
        assert_eq!(
            certified.iterations, 0,
            "the fixture must be one the certificate stops before any step"
        );

        let mut beta = vec![1e-9];
        let report = solve_polish(&problem, &mut beta, &opts(50), &|| false);
        assert_eq!(
            report.iterations, 1,
            "the promised step precedes the certificate"
        );
        assert_eq!(beta[0], 0.0, "and it reaches the minimizer");
        assert!(report.converged);
    }

    /// The distance between adjacent doubles at ten million, two to the minus
    /// twenty-ninth. It is the whole step available to a parameter of that
    /// magnitude, and both fixtures below are built around it.
    const DOUBLE_SPACING: f64 = 1.862645149230957e-9;

    /// A one-parameter least squares whose four responses sit ten million from the
    /// origin, so its minimizer does too. Two resolutions are in play and holding
    /// them apart is what the fixture exists for. The residuals at the minimizer
    /// are of order one, so the objective is near five and cannot tell apart two
    /// values closer than about `eps` times that; the parameter, at ten million,
    /// moves only in whole doubles of [`DOUBLE_SPACING`]. A Newton step of a
    /// couple of those predicts a decrease of order `1e-17`, well beneath anything
    /// the objective registers, while the gradient the step would remove is still
    /// orders of magnitude above a tolerance of `1e-10`.
    ///
    /// That is the shape of the entropy dual at an iterate whose value has
    /// flattened. The dual is a sum over the sample and stops moving in its last
    /// digit long before the estimating equations are solved, and the line search
    /// is then comparing two values that are equal to rounding.
    ///
    /// The gradient and Hessian are exact, and the default answer to
    /// [`EsteqProblem::hessian_is_exact`] is the one that matters here. Each
    /// residual is the difference of two nearby doubles and is computed without
    /// error, and the Hessian is the unit count.
    struct DistantLeastSquares {
        y: [f64; 4],
    }

    impl EsteqProblem for DistantLeastSquares {
        fn n_params(&self) -> usize {
            1
        }
        fn n_units(&self) -> usize {
            self.y.len()
        }
        fn value(&self, beta: &[f64]) -> Option<f64> {
            let mut acc = 0.0;
            for y in &self.y {
                let residual = beta[0] - y;
                acc += 0.5 * residual * residual;
            }
            Some(acc)
        }
        fn gradient(&self, beta: &[f64], g: &mut [f64]) {
            let mut acc = 0.0;
            for y in &self.y {
                acc += beta[0] - y;
            }
            g[0] = acc;
        }
        fn hessian(&self, _beta: &[f64], mut h: MatMut<'_, f64>) {
            *h.rb_mut().get_mut(0, 0) = self.y.len() as f64;
        }
        fn psi(&self, _beta: &[f64], _out: MatMut<'_, f64>) {}
    }

    /// The options both resolution fixtures run under: a tolerance loose enough
    /// that the gradients in play sit far above it, and a cap large enough that
    /// reaching it is unambiguously a crawl rather than a tight budget.
    fn resolution_opts() -> SolveOptions {
        SolveOptions {
            max_iter: 1000,
            grad_tol: 1e-10,
            fista_rel_tol: 1e-12,
        }
    }

    #[test]
    fn an_exact_step_is_taken_where_the_objective_cannot_resolve_the_decrease() {
        // The iterate is two doubles above the minimizer. Its gradient is 1.5e-8,
        // a hundred and fifty times the tolerance, and one Newton step removes it
        // exactly. What the objective sees is nothing: the step predicts a
        // decrease of 5.6e-17 against a resolution of 1.1e-15, so the
        // sufficient-decrease test compares two values that are equal to rounding
        // and decides between them by ulps. It refuses the full step, refuses the
        // half, and accepts the quarter, which is below half a double at this
        // magnitude and so moves the iterate nowhere at all. Every pass that
        // follows repeats it, and the solve spends its whole budget standing
        // still.
        //
        // The exact Hessian is what makes the step worth taking rather than the
        // iterate worth certifying. A superlinear step goes on sharpening the
        // parameter after the value has stopped registering the improvement, and
        // here it lands on the solution, so the objective's silence is no evidence
        // that the iteration has finished.
        let problem = DistantLeastSquares {
            y: [10000001.8, 9999997.7, 10000001.1, 9999999.4],
        };
        let opts = resolution_opts();
        let start = 1e7 + 2.0 * DOUBLE_SPACING;

        let mut g = vec![0.0];
        problem.gradient(&[start], &mut g);
        let base = problem.value(&[start]).expect("smooth objective");
        let predicted = g[0] * g[0] / problem.n_units() as f64;
        assert!(
            predicted <= value_resolution(base),
            "the scenario needs a predicted decrease the objective cannot resolve, \
             got {predicted:e} against a resolution of {:e}",
            value_resolution(base)
        );
        assert!(
            g[0].abs() > 100.0 * opts.grad_tol,
            "and a gradient far above the tolerance, got {:e}",
            g[0]
        );

        let mut beta = vec![start];
        let report = solve(&problem, &mut beta, &opts, &|| false);
        assert!(
            report.grad_norm <= opts.grad_tol,
            "the step the objective cannot resolve is the one that solves the \
             estimating equations, so the returned gradient must meet the \
             tolerance, got {:e}",
            report.grad_norm
        );
        assert!(report.converged);
        assert!(
            report.iterations <= 5,
            "a few steps rather than the iteration cap, got {}",
            report.iterations
        );
        assert_eq!(beta[0], 1e7, "and the step reaches the minimizer");
    }

    #[test]
    fn a_gradient_that_cannot_fall_is_certified_rather_than_crawled() {
        // The same objective with one response displaced by a single double, which
        // is enough that no representable parameter zeroes the four residuals. The
        // gradient reaches one double-spacing, 1.9e-9, and no less, however well
        // the problem is solved; the predicted decrease there is 8.7e-19 against a
        // resolution of 1.1e-15, and the full Newton step is a quarter of a double
        // and moves nothing.
        //
        // This is the case the exact step cannot rescue, and the solve has to say
        // so rather than spend its budget rediscovering it. Taking the step and
        // finding the gradient unmoved is the evidence: the iterate is as far as
        // the arithmetic reaches, so the step is given back and convergence is
        // reported on the resolution certificate, at whatever count the real steps
        // came to.
        let problem = DistantLeastSquares {
            y: [
                10000001.8,
                9999997.7,
                10000001.1,
                9999999.4 + DOUBLE_SPACING,
            ],
        };
        let opts = resolution_opts();
        let start = 1e7;

        let mut g = vec![0.0];
        problem.gradient(&[start], &mut g);
        let base = problem.value(&[start]).expect("smooth objective");
        let predicted = g[0] * g[0] / problem.n_units() as f64;
        assert!(
            predicted <= value_resolution(base),
            "the scenario needs a predicted decrease the objective cannot resolve, \
             got {predicted:e} against a resolution of {:e}",
            value_resolution(base)
        );
        // The floor is a property of the problem rather than of wherever the solve
        // happened to stop: every neighbouring double carries a larger gradient,
        // so no step reaches a smaller one.
        for offset in [-2.0_f64, -1.0, 1.0, 2.0] {
            let mut neighbour = vec![0.0];
            problem.gradient(&[start + offset * DOUBLE_SPACING], &mut neighbour);
            assert!(
                neighbour[0].abs() > g[0].abs(),
                "the double {offset} away carries a gradient of {:e}, not above \
                 the floor of {:e}",
                neighbour[0],
                g[0]
            );
        }

        let mut beta = vec![start];
        let report = solve(&problem, &mut beta, &opts, &|| false);
        assert_eq!(
            report.iterations, 0,
            "a reverted probe embodies nothing, so the count must stay at zero \
             rather than crawl toward the cap"
        );
        assert!(
            report.converged,
            "an iterate the arithmetic cannot improve on is the optimum it is"
        );
        assert_eq!(beta[0], start, "the probe step is given back");
        assert_eq!(
            report.grad_norm, DOUBLE_SPACING,
            "and the gradient stands at the floor the problem sets"
        );
        assert!(
            report.grad_norm > opts.grad_tol,
            "which is a floor the tolerance cannot reach"
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
