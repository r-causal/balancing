//! Inverse probability tilting: propensity weights whose tilted score equations
//! force the weighted covariate means to their estimand targets.
//!
//! For each treatment level a propensity model `p_i = G(x_i . beta_k)` is fit,
//! not by maximum likelihood but by the tilting moment condition
//!
//! ```text
//! sum_i s_i ( tau_i - m_ik w_i(beta_k) ) x_i = 0,
//! ```
//!
//! where `m_ik` marks the units in level `k`, `w_i` is the level's weight
//! function, and `tau_i` is the target population indicator. The moment forces
//! the level's weighted covariate total to the target total, so balance on the
//! constraint moments is exact by construction. The average-treatment-effect
//! estimand tilts every level to the whole sample (`tau_i = 1`, `w_i = 1 / p_i`);
//! a focal estimand tilts each non-focal level to the focal level
//! (`tau_i = 1{T_i = focal}`, `w_i = (1 - p_i) / p_i`) and leaves the focal
//! units at weight one.
//!
//! Each level's moment is an independent root-finding problem whose Jacobian
//! `sum_i s_i (mu_eta_i / p_i^2) x_i x_i'` is positive definite, so the damped
//! Newton method in merit mode converges from the GLM warm start. The blocks are
//! solved in level order and their coefficients, Jacobians, and weight
//! derivatives stack block-diagonally for the estimating-equations output.

use std::sync::Arc;

use faer::MatMut;
use faer::prelude::ReborrowMut;
use rayon::ThreadPool;
use rayon::iter::{IndexedParallelIterator, ParallelIterator};
use rayon::slice::ParallelSliceMut;

use crate::esteq::{self, EsteqProblem, SolveOptions, Solver};
use crate::glm;
use crate::links::Link;
use crate::threads::{deterministic_map_reduce, get_pool};
use crate::weights::{WeightForm, weight_deriv_eta};

/// Which population the tilt targets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IptEstimand {
    /// Average treatment effect: every level is tilted to the whole sample.
    Ate,
    /// A focal estimand: every non-focal level is tilted to the focal level,
    /// whose units keep weight one. The treatment effect on the treated is the
    /// focal-equals-treated case; on the controls, the focal-equals-control case.
    Focal(usize),
}

/// Inputs for an inverse probability tilting solve.
pub struct IptInputs<'a> {
    /// Column-major `n` by `p` design matrix, including any intercept column.
    pub covs: &'a [f64],
    /// Number of units.
    pub n: usize,
    /// Number of model terms.
    pub p: usize,
    /// Treatment level of each unit, values `0..n_levels`.
    pub treat: &'a [i32],
    /// Number of distinct treatment levels.
    pub n_levels: usize,
    /// Sampling weights `s_i`.
    pub s: &'a [f64],
    /// Propensity link.
    pub link: Link,
    /// Target population.
    pub estimand: IptEstimand,
    /// Worker threads for the reductions.
    pub threads: usize,
    /// Maximum Newton iterations per block.
    pub max_iter: usize,
    /// Convergence threshold on the moment sup norm.
    pub tol: f64,
}

/// Solution of an inverse probability tilting problem.
pub struct IptResult {
    /// Balancing weights, one per unit, in the original unit order.
    pub weights: Vec<f64>,
    /// Modeled propensity of each unit under its level's block. Units in a level
    /// with no block (the focal level of a focal estimand) carry `NaN`, having no
    /// fitted model.
    pub ps: Vec<f64>,
    /// Link coefficients, `p` per block, stacked in block order.
    pub coefs: Vec<f64>,
    /// Whether every block met its convergence criterion.
    pub converged: bool,
    /// Whether any block's solve stopped because a user interrupt was pending.
    pub interrupted: bool,
    /// Largest iteration count across blocks.
    pub iterations: usize,
    /// Largest moment sup norm across blocks.
    pub grad_norm: f64,
    /// Column-major `n` by `P` per-unit estimating functions, `P = p` times the
    /// number of blocks.
    pub psi: Vec<f64>,
    /// Column-major `P` by `P` block-diagonal Jacobian of the summed estimating
    /// functions.
    pub jac: Vec<f64>,
    /// Column-major `n` by `P` derivative of each unit's weight in the parameters.
    pub dw_dbeta: Vec<f64>,
    /// Treatment level each parameter block belongs to, in block order.
    pub block_levels: Vec<usize>,
}

/// The tilting moment problem for a single treatment level.
///
/// `idx` selects the level's units, whose weights carry the parameter
/// dependence; `target` is the fixed population total the moment balances to.
struct IptBlock<'a> {
    covs: &'a [f64],
    n: usize,
    p: usize,
    idx: &'a [usize],
    s: &'a [f64],
    target: &'a [f64],
    link: Link,
    form: WeightForm,
    pool: Arc<ThreadPool>,
}

/// Per-evaluation sufficient statistics: the achieved weighted covariate total
/// and the moment Jacobian.
struct BlockAccum {
    a: Vec<f64>,
    jac: Vec<f64>,
    row: Vec<f64>,
}

impl BlockAccum {
    fn zeros(p: usize, want_jac: bool) -> Self {
        Self {
            a: vec![0.0; p],
            jac: if want_jac {
                vec![0.0; p * p]
            } else {
                Vec::new()
            },
            row: vec![0.0; p],
        }
    }
}

impl IptBlock<'_> {
    /// Linear predictor `x_i . beta` for the level unit at global index `gi`.
    fn lin(&self, gi: usize, beta: &[f64]) -> f64 {
        let mut acc = 0.0;
        for (j, &b) in beta.iter().enumerate() {
            acc += self.covs[j * self.n + gi] * b;
        }
        acc
    }

    /// Fold the level's tilted contributions: the achieved total `a` always, and
    /// the moment Jacobian when `want_jac`.
    fn accumulate(&self, beta: &[f64], want_jac: bool) -> BlockAccum {
        let p = self.p;
        deterministic_map_reduce(
            &self.pool,
            self.idx.len(),
            || BlockAccum::zeros(p, want_jac),
            |acc, local| {
                let gi = self.idx[local];
                for j in 0..p {
                    acc.row[j] = self.covs[j * self.n + gi];
                }
                let eta: f64 = acc.row.iter().zip(beta).map(|(c, b)| c * b).sum();
                let prob = self.link.linkinv(eta);
                let w = self.form.weight(prob);
                let sw = self.s[gi] * w;
                for j in 0..p {
                    acc.a[j] += sw * acc.row[j];
                }
                if want_jac {
                    // dg/dbeta = -dA/dbeta = sum_i s_i (mu_eta / p^2) x_i x_i',
                    // which is positive definite and drives the Newton step.
                    let mu_eta = self.link.mu_eta(eta);
                    let curv = self.s[gi] * mu_eta / (prob * prob);
                    for j in 0..p {
                        let cj = curv * acc.row[j];
                        for k in 0..p {
                            acc.jac[j * p + k] += cj * acc.row[k];
                        }
                    }
                }
            },
            |acc, other| {
                for j in 0..p {
                    acc.a[j] += other.a[j];
                }
                if want_jac {
                    for jk in 0..p * p {
                        acc.jac[jk] += other.jac[jk];
                    }
                }
            },
        )
    }
}

impl EsteqProblem for IptBlock<'_> {
    fn n_params(&self) -> usize {
        self.p
    }

    fn n_units(&self) -> usize {
        self.idx.len()
    }

    fn value(&self, _beta: &[f64]) -> Option<f64> {
        // A pure root-finding problem: the Newton method minimizes the merit
        // 0.5 ||g||^2 rather than an objective.
        None
    }

    fn gradient(&self, beta: &[f64], g: &mut [f64]) {
        let acc = self.accumulate(beta, false);
        for (j, gj) in g.iter_mut().enumerate() {
            *gj = self.target[j] - acc.a[j];
        }
    }

    fn hessian(&self, beta: &[f64], mut h: MatMut<'_, f64>) {
        let acc = self.accumulate(beta, true);
        for j in 0..self.p {
            for k in 0..self.p {
                *h.rb_mut().get_mut(j, k) = acc.jac[j * self.p + k];
            }
        }
    }

    fn value_grad_hess(&self, beta: &[f64], g: &mut [f64], mut h: MatMut<'_, f64>) -> Option<f64> {
        let acc = self.accumulate(beta, true);
        for (j, gj) in g.iter_mut().enumerate() {
            *gj = self.target[j] - acc.a[j];
            for k in 0..self.p {
                *h.rb_mut().get_mut(j, k) = acc.jac[j * self.p + k];
            }
        }
        None
    }

    fn psi(&self, beta: &[f64], mut out: MatMut<'_, f64>) {
        // The level's own-unit contributions; the full per-unit matrix, which
        // also carries every unit's target term, is assembled in `solve`. Not on
        // the solver's path, this exists to satisfy the trait.
        for (local, &gi) in self.idx.iter().enumerate() {
            let prob = self.link.linkinv(self.lin(gi, beta));
            let w = self.form.weight(prob);
            for j in 0..self.p {
                *out.rb_mut().get_mut(local, j) = -self.s[gi] * w * self.covs[j * self.n + gi];
            }
        }
    }
}

/// The blocks to solve, as `(level, weight form)` pairs, and the target
/// population indicator `tau_i` shared across them.
fn plan(inputs: &IptInputs<'_>) -> (Vec<(usize, WeightForm)>, Vec<f64>) {
    match inputs.estimand {
        IptEstimand::Ate => {
            let blocks = (0..inputs.n_levels)
                .map(|k| (k, WeightForm::Inverse))
                .collect();
            let tau = vec![1.0; inputs.n];
            (blocks, tau)
        }
        IptEstimand::Focal(focal) => {
            let blocks = (0..inputs.n_levels)
                .filter(|&k| k != focal)
                .map(|k| (k, WeightForm::InverseComplement))
                .collect();
            let tau = inputs
                .treat
                .iter()
                .map(|&t| if t as usize == focal { 1.0 } else { 0.0 })
                .collect();
            (blocks, tau)
        }
    }
}

/// Solve an inverse probability tilting problem across its treatment levels.
pub fn solve(inputs: &IptInputs<'_>, interrupt: &dyn Fn() -> bool) -> IptResult {
    let n = inputs.n;
    let p = inputs.p;
    let pool = get_pool(inputs.threads);
    let (blocks, tau) = plan(inputs);
    let n_blocks = blocks.len();
    let total_params = n_blocks * p;

    // The target total sum_i s_i tau_i x_i is shared by every block.
    let mut target = vec![0.0; p];
    for (i, &tau_i) in tau.iter().enumerate() {
        let stau = inputs.s[i] * tau_i;
        if stau != 0.0 {
            for (j, t) in target.iter_mut().enumerate() {
                *t += stau * inputs.covs[j * n + i];
            }
        }
    }

    // A level's units carry weight one unless a block overrides them, and their
    // propensity is undefined without a block.
    let mut weights = vec![1.0; n];
    let mut ps = vec![f64::NAN; n];
    let mut coefs = vec![0.0; total_params];
    let mut psi = vec![0.0; n * total_params];
    let mut jac = vec![0.0; total_params * total_params];
    let mut dw = vec![0.0; n * total_params];
    let mut block_levels = Vec::with_capacity(n_blocks);

    let mut converged = true;
    let mut interrupted = false;
    let mut iterations = 0;
    let mut grad_norm = 0.0_f64;

    let solve_opts = SolveOptions {
        max_iter: inputs.max_iter,
        grad_tol: inputs.tol,
        fista_rel_tol: inputs.tol,
    };

    for (b, &(level, form)) in blocks.iter().enumerate() {
        block_levels.push(level);
        let idx: Vec<usize> = inputs
            .treat
            .iter()
            .enumerate()
            .filter_map(|(i, &t)| (t as usize == level).then_some(i))
            .collect();
        if idx.is_empty() {
            continue;
        }

        let mut beta = warm_start(inputs, level, form);

        let block = IptBlock {
            covs: inputs.covs,
            n,
            p,
            idx: &idx,
            s: inputs.s,
            target: &target,
            link: inputs.link,
            form,
            pool: Arc::clone(&pool),
        };
        let report = esteq::solve(&block, &mut beta, Solver::Newton, &solve_opts, interrupt);
        converged &= report.converged;
        interrupted |= report.interrupted;
        iterations = iterations.max(report.iterations);
        grad_norm = grad_norm.max(report.grad_norm);

        coefs[b * p..b * p + p].copy_from_slice(&beta);
        fill_block_output(
            &pool,
            inputs,
            &block,
            &beta,
            &idx,
            &tau,
            b,
            total_params,
            &mut weights,
            &mut ps,
            &mut psi,
            &mut jac,
            &mut dw,
        );
    }

    IptResult {
        weights,
        ps,
        coefs,
        converged: converged && !interrupted,
        interrupted,
        iterations,
        grad_norm,
        psi,
        jac,
        dw_dbeta: dw,
        block_levels,
    }
}

/// GLM warm start for a block. The average-treatment-effect blocks fit each
/// level against the whole sample; the focal blocks fit each non-focal level
/// against the focal level only, matching the two-class comparison the tilt
/// weight expresses.
fn warm_start(inputs: &IptInputs<'_>, level: usize, form: WeightForm) -> Vec<f64> {
    let (fit_idx, response): (Vec<usize>, Vec<f64>) = match form {
        WeightForm::Inverse => {
            let response = inputs
                .treat
                .iter()
                .map(|&t| f64::from(t as usize == level))
                .collect();
            ((0..inputs.n).collect(), response)
        }
        WeightForm::InverseComplement => {
            let IptEstimand::Focal(focal) = inputs.estimand else {
                // A complement weight is only ever planned for a focal estimand.
                return vec![0.0; inputs.p];
            };
            let mut fit_idx = Vec::new();
            let mut response = Vec::new();
            for (i, &t) in inputs.treat.iter().enumerate() {
                let t = t as usize;
                if t == level || t == focal {
                    fit_idx.push(i);
                    response.push(f64::from(t == level));
                }
            }
            (fit_idx, response)
        }
    };
    glm::irls(
        inputs.covs,
        inputs.n,
        inputs.p,
        &fit_idx,
        &response,
        inputs.s,
        inputs.link,
        50,
    )
}

/// Fill one block's weights, propensities, and estimating-equation blocks.
#[allow(clippy::too_many_arguments)]
fn fill_block_output(
    pool: &ThreadPool,
    inputs: &IptInputs<'_>,
    block: &IptBlock<'_>,
    beta: &[f64],
    idx: &[usize],
    tau: &[f64],
    b: usize,
    total_params: usize,
    weights: &mut [f64],
    ps: &mut [f64],
    psi: &mut [f64],
    jac: &mut [f64],
    dw: &mut [f64],
) {
    let n = inputs.n;
    let p = inputs.p;
    let offset = b * p;

    // Per-level-unit propensity, weight, and weight derivative at the solution,
    // computed once so the column fills below only gather covariates.
    let mut wvec = vec![0.0; idx.len()];
    let mut dwe = vec![0.0; idx.len()];
    for (local, &gi) in idx.iter().enumerate() {
        let eta = block.lin(gi, beta);
        let prob = inputs.link.linkinv(eta);
        let mu_eta = inputs.link.mu_eta(eta);
        let w = block.form.weight(prob);
        weights[gi] = w;
        ps[gi] = prob;
        wvec[local] = w;
        dwe[local] = weight_deriv_eta(prob, mu_eta);
    }

    // psi column j of this block: s_i tau_i x_ij for every unit, less the level's
    // own weighted term. dw column j: the weight derivative on the level's units.
    // Both are elementwise per column, so they fill in parallel over columns.
    let psi_block = &mut psi[offset * n..(offset + p) * n];
    let dw_block = &mut dw[offset * n..(offset + p) * n];
    pool.install(|| {
        psi_block
            .par_chunks_mut(n)
            .zip(dw_block.par_chunks_mut(n))
            .enumerate()
            .for_each(|(j, (psi_col, dw_col))| {
                let covs_col = &inputs.covs[j * n..(j + 1) * n];
                for i in 0..n {
                    psi_col[i] = inputs.s[i] * tau[i] * covs_col[i];
                }
                for (local, &gi) in idx.iter().enumerate() {
                    psi_col[gi] -= inputs.s[gi] * wvec[local] * covs_col[gi];
                    dw_col[gi] = dwe[local] * covs_col[gi];
                }
            });
    });

    // The Jacobian block equals the moment Jacobian, positive definite by
    // construction, placed on the diagonal block for this level.
    let acc = block.accumulate(beta, true);
    for j in 0..p {
        let row = offset + j;
        for k in 0..p {
            let col = offset + k;
            jac[col * total_params + row] = acc.jac[j * p + k];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn no_interrupt() -> impl Fn() -> bool {
        || false
    }

    /// A saturated two-cell binary design: covariates are an intercept and a 0/1
    /// cell indicator, so the tilt recovers exact empirical cell propensities.
    ///
    /// Cell z = 0: 4 units, 2 treated. Cell z = 1: 6 units, 2 treated.
    fn saturated_design() -> (Vec<f64>, Vec<i32>, usize) {
        let z = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        let treat = vec![1, 1, 0, 0, 1, 1, 0, 0, 0, 0];
        let n = z.len();
        let mut covs = vec![1.0; n];
        covs.extend_from_slice(&z);
        (covs, treat, n)
    }

    fn run(estimand: IptEstimand, threads: usize) -> IptResult {
        let (covs, treat, n) = saturated_design();
        let s = vec![1.0; n];
        let inputs = IptInputs {
            covs: &covs,
            n,
            p: 2,
            treat: &treat,
            n_levels: 2,
            s: &s,
            link: Link::Logit,
            estimand,
            threads,
            max_iter: 200,
            tol: 1e-12,
        };
        solve(&inputs, &no_interrupt())
    }

    // On a saturated design the average-treatment-effect weights are the inverse
    // empirical cell propensities: treated 1/e, control 1/(1 - e).
    #[test]
    fn ate_recovers_inverse_cell_propensities() {
        let result = run(IptEstimand::Ate, 1);
        assert!(result.converged);
        // Cell 0: e = 0.5 -> treated 2, control 2. Cell 1: e = 1/3 -> treated 3,
        // control 1.5.
        let expected = [2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 1.5, 1.5, 1.5, 1.5];
        for (i, &want) in expected.iter().enumerate() {
            assert!(
                (result.weights[i] - want).abs() < 1e-8,
                "weight[{i}] = {} expected {want}",
                result.weights[i]
            );
        }
        assert!(result.weights.iter().all(|&w| w >= 0.0));
    }

    // Focal on the treated: treated units keep weight one; each control's weight
    // is the treated-to-control count ratio in its cell.
    #[test]
    fn att_keeps_treated_at_one_and_tilts_controls() {
        let result = run(IptEstimand::Focal(1), 1);
        assert!(result.converged);
        // Cell 0: 2 treated / 2 control -> control weight 1. Cell 1: 2 treated /
        // 4 control -> control weight 0.5. Treated stay at 1.
        let expected = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5];
        for (i, &want) in expected.iter().enumerate() {
            assert!(
                (result.weights[i] - want).abs() < 1e-8,
                "weight[{i}] = {} expected {want}",
                result.weights[i]
            );
        }
    }

    // The estimating functions sum to zero at the solution, column by column.
    #[test]
    fn psi_columns_sum_to_zero() {
        for estimand in [IptEstimand::Ate, IptEstimand::Focal(1)] {
            let result = run(estimand, 1);
            let n = result.weights.len();
            let total = result.coefs.len();
            for col in 0..total {
                let s: f64 = (0..n).map(|i| result.psi[col * n + i]).sum();
                assert!(s.abs() < 1e-8, "column {col} sums to {s}");
            }
        }
    }

    // The weight derivative matches a central difference of the raw weight in the
    // block's coefficients, for an interior unit.
    #[test]
    fn weight_derivative_matches_finite_difference() {
        let (covs, treat, n) = saturated_design();
        let s = vec![1.0; n];
        let mk = |covs: &[f64]| {
            let inputs = IptInputs {
                covs,
                n,
                p: 2,
                treat: &treat,
                n_levels: 2,
                s: &s,
                link: Link::Logit,
                estimand: IptEstimand::Ate,
                threads: 1,
                max_iter: 200,
                tol: 1e-12,
            };
            solve(&inputs, &no_interrupt())
        };
        let result = mk(&covs);
        // Unit 0 is in level 1, whose block occupies columns [p .. 2p). Perturb
        // its own coefficients and read the weight change; here we instead check
        // the analytic derivative against a direct link computation.
        let unit = 0usize;
        let block = 1usize; // level 1 is the second block in ATE order
        let p = 2usize;
        let beta = &result.coefs[block * p..block * p + p];
        let link = Link::Logit;
        let eta: f64 = (0..p).map(|j| covs[j * n + unit] * beta[j]).sum();
        let h = 1e-6;
        for j in 0..p {
            let x = covs[j * n + unit];
            let w = |d: f64| 1.0 / link.linkinv(eta + d * x);
            let fd = (w(h) - w(-h)) / (2.0 * h);
            let analytic = result.dw_dbeta[(block * p + j) * n + unit];
            assert!(
                (fd - analytic).abs() < 1e-5,
                "col {j}: fd {fd} analytic {analytic}"
            );
        }
    }

    // Same binary inputs, any thread count: bit-identical weights and coefficients.
    #[test]
    fn solve_is_deterministic_across_thread_counts() {
        let one = run(IptEstimand::Ate, 1);
        for threads in [2, 4, 8] {
            let many = run(IptEstimand::Ate, threads);
            for i in 0..one.weights.len() {
                assert_eq!(one.weights[i].to_bits(), many.weights[i].to_bits());
            }
            for i in 0..one.coefs.len() {
                assert_eq!(one.coefs[i].to_bits(), many.coefs[i].to_bits());
            }
        }
    }

    // A pending interrupt stops the solve and is reported rather than swallowed.
    #[test]
    fn interrupt_is_surfaced() {
        let (covs, treat, n) = saturated_design();
        let s = vec![1.0; n];
        let inputs = IptInputs {
            covs: &covs,
            n,
            p: 2,
            treat: &treat,
            n_levels: 2,
            s: &s,
            link: Link::Logit,
            estimand: IptEstimand::Ate,
            threads: 1,
            max_iter: 200,
            tol: 1e-12,
        };
        let result = solve(&inputs, &|| true);
        assert!(result.interrupted);
        assert!(!result.converged);
    }

    // Every link fits the saturated design and returns non-negative weights.
    #[test]
    fn every_link_balances_the_saturated_design() {
        let (covs, treat, n) = saturated_design();
        let s = vec![1.0; n];
        for link in [Link::Logit, Link::Probit, Link::Cloglog] {
            let inputs = IptInputs {
                covs: &covs,
                n,
                p: 2,
                treat: &treat,
                n_levels: 2,
                s: &s,
                link,
                estimand: IptEstimand::Ate,
                threads: 1,
                max_iter: 200,
                tol: 1e-12,
            };
            let result = solve(&inputs, &no_interrupt());
            assert!(result.converged, "link {link:?} did not converge");
            // The saturated cell propensities are identical across links, so the
            // weights match the inverse empirical propensities regardless.
            let expected = [2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 1.5, 1.5, 1.5, 1.5];
            for (i, &want) in expected.iter().enumerate() {
                assert!(
                    (result.weights[i] - want).abs() < 1e-7,
                    "link {link:?} weight[{i}] = {} expected {want}",
                    result.weights[i]
                );
            }
        }
    }
}
