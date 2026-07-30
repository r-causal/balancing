//! Covariate balancing propensity score: a propensity model whose coefficients
//! are chosen so the weighted covariates balance.
//!
//! The just-identified estimator replaces the maximum-likelihood score
//! equations with the covariate balancing moment conditions, so the number of
//! conditions equals the number of parameters and balance holds exactly at the
//! solution. For a binary treatment with a single propensity model
//! `p_i = G(x_i . beta)` the balancing condition is
//!
//! ```text
//! sum_i s_i c_i(beta) x_i = 0,   c_i = (p_i - T_i) / d_i,
//! ```
//!
//! where the denominator `d_i` selects the estimand: `p_i (1 - p_i)` for the
//! average treatment effect, `1 - p_i` for the effect on the treated, `p_i` for
//! the effect on the controls, and one for the overlap-weighted estimand. Each
//! condition equals the requirement that the weighted covariate totals of the
//! two treatment arms agree, so the inverse-propensity weights
//! `w_i = T_i / p_i + (1 - T_i) / (1 - p_i)` (and their focal and overlap
//! analogues) balance the covariates by construction. The moment Jacobian
//! `sum_i s_i (dc_i / dp) mu_eta_i x_i x_i'` is positive semidefinite, so the
//! damped Newton method in merit mode converges from the maximum-likelihood
//! warm start.
//!
//! The over-identified estimator stacks the response-residual score conditions
//! `sum_i s_i (T_i - p_i) x_i` onto the balancing conditions and minimizes the
//! generalized-method-of-moments criterion `m(beta)' W m(beta)`. The weighting
//! matrix `W` is the pseudo-inverse of the moment covariance, held fixed at a
//! preliminary estimate for the two-step variant or recomputed at each iterate
//! for the continuously-updated variant. The criterion is minimized by the
//! limited-memory BFGS solver followed by a damped Newton polish: the quasi-Newton
//! pass needs only the objective and its gradient and does the work of arriving,
//! and the polish reads the stationarity certificate that recognizes a criterion
//! flat to the resolution of its own arithmetic. Both stages run against the
//! Gauss-Newton Hessian `2 G' W G`, and where the quasi-Newton backend is not
//! compiled in the Newton method takes over the first stage too.
//!
//! A continuous exposure balances the weighted exposure-covariate covariance.
//! The covariate balancing conditions ask the weighted exposure mean to match
//! the sample mean and the weighted covariance between the exposure and every
//! covariate to vanish; the minimum-divergence positive weights that meet them
//! are the exponential tilt whose parameters solve a convex dual, the same
//! Newton machinery the discrete estimators use.

use std::sync::Arc;

use faer::MatMut;
use faer::prelude::ReborrowMut;
use rayon::ThreadPool;
use rayon::iter::{IndexedParallelIterator, ParallelIterator};
use rayon::slice::ParallelSliceMut;

use crate::esteq::{self, EsteqProblem, SolveOptions, Solver, sampling_weight_scale};
use crate::glm;
use crate::linalg::pseudo_inverse_symmetric;
use crate::links::Link;
use crate::methods::ipt::{self, IptEstimand, IptInputs};
use crate::threads::{deterministic_map_reduce, get_pool};

/// Propensities are held away from zero and one inside the generalized-method-of
/// -moments objective. The average-treatment-effect balancing factor divides by
/// `p (1 - p)`, so an unclamped line-search excursion toward the boundary would
/// overflow the criterion; the interior optimum is unaffected by the clamp.
const GMM_PROB_FLOOR: f64 = 1e-8;

/// Clamp a propensity into the open unit interval for the GMM objective.
fn clamp_prob(p: f64) -> f64 {
    p.clamp(GMM_PROB_FLOOR, 1.0 - GMM_PROB_FLOOR)
}

/// Which population the covariate balancing propensity score targets, for a
/// binary treatment coded `1` for treated and `0` for control.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CbpsEstimand {
    /// Average treatment effect: both arms weighted to the whole sample.
    Ate,
    /// Effect on the treated: treated units keep weight one, controls tilted.
    Att,
    /// Effect on the controls: control units keep weight one, treated tilted.
    Atc,
    /// Overlap-weighted estimand: weights proportional to `p (1 - p)`.
    Ato,
}

impl CbpsEstimand {
    /// The balancing weight a unit carries at propensity `p`, given its
    /// treatment indicator `t` in `{0, 1}`.
    fn weight(self, p: f64, t: f64) -> f64 {
        match self {
            CbpsEstimand::Ate => t / p + (1.0 - t) / (1.0 - p),
            CbpsEstimand::Att => t + (1.0 - t) * p / (1.0 - p),
            CbpsEstimand::Atc => t * (1.0 - p) / p + (1.0 - t),
            CbpsEstimand::Ato => t * (1.0 - p) + (1.0 - t) * p,
        }
    }

    /// The derivative of the weight in the propensity, `dw / dp`.
    fn weight_deriv_p(self, p: f64, t: f64) -> f64 {
        match self {
            CbpsEstimand::Ate => -t / (p * p) + (1.0 - t) / ((1.0 - p) * (1.0 - p)),
            CbpsEstimand::Att => (1.0 - t) / ((1.0 - p) * (1.0 - p)),
            CbpsEstimand::Atc => -t / (p * p),
            CbpsEstimand::Ato => 1.0 - 2.0 * t,
        }
    }

    /// The scalar balancing factor `c_i = (p - t) / d`, whose weighted covariate
    /// total is the moment condition that is driven to zero.
    fn bal_factor(self, p: f64, t: f64) -> f64 {
        match self {
            CbpsEstimand::Ate => (p - t) / (p * (1.0 - p)),
            CbpsEstimand::Att => (p - t) / (1.0 - p),
            CbpsEstimand::Atc => (p - t) / p,
            CbpsEstimand::Ato => p - t,
        }
    }

    /// The derivative of the balancing factor in the propensity, `dc / dp`. It
    /// is non-negative for every estimand, so the moment Jacobian is positive
    /// semidefinite.
    fn bal_factor_deriv_p(self, p: f64, t: f64) -> f64 {
        match self {
            CbpsEstimand::Ate => t / (p * p) + (1.0 - t) / ((1.0 - p) * (1.0 - p)),
            CbpsEstimand::Att => (1.0 - t) / ((1.0 - p) * (1.0 - p)),
            CbpsEstimand::Atc => t / (p * p),
            CbpsEstimand::Ato => 1.0,
        }
    }
}

/// Inputs for a binary covariate balancing propensity score solve.
pub struct CbpsInputs<'a> {
    /// Column-major `n` by `p_mod` propensity-model design, including intercept.
    pub covs_mod: &'a [f64],
    /// Column-major `n` by `p_bal` balance design. For the just-identified form
    /// it must have the same shape as `covs_mod`, which is the standard case.
    pub covs_bal: &'a [f64],
    /// Number of units.
    pub n: usize,
    /// Number of propensity-model terms, the parameter count.
    pub p_mod: usize,
    /// Number of balance terms.
    pub p_bal: usize,
    /// Treatment indicator, `1` treated and `0` control.
    pub treat: &'a [i32],
    /// Sampling weights `s_i`.
    pub s: &'a [f64],
    /// Propensity link.
    pub link: Link,
    /// Target population.
    pub estimand: CbpsEstimand,
    /// Whether to stack the score conditions and minimize the GMM criterion.
    pub over: bool,
    /// Whether the over-identified variant uses the two-step weighting matrix
    /// rather than continuous updating. Ignored when `over` is false.
    pub twostep: bool,
    /// Worker threads for the reductions.
    pub threads: usize,
    /// Maximum solver iterations.
    pub max_iter: usize,
    /// Convergence threshold on the moment sup norm or the objective gradient.
    pub tol: f64,
}

/// Solution of a covariate balancing propensity score problem.
pub struct CbpsResult {
    /// Balancing weights, one per unit, in the original unit order.
    pub weights: Vec<f64>,
    /// Fitted propensity of each unit, `NaN` for a continuous exposure.
    pub ps: Vec<f64>,
    /// Link coefficients.
    pub coefs: Vec<f64>,
    /// Whether the solve met its convergence criterion.
    pub converged: bool,
    /// Whether the solve stopped because a user interrupt was pending.
    pub interrupted: bool,
    /// Iterations performed.
    pub iterations: usize,
    /// Moment sup norm (just-identified) or GMM criterion (over-identified).
    pub obj_value: f64,
    /// Column-major `n` by `P` per-unit estimating functions, `None` for the
    /// over-identified and continuous forms.
    pub psi: Option<Vec<f64>>,
    /// Column-major `P` by `P` Jacobian of the summed estimating functions,
    /// `None` for the over-identified and continuous forms.
    pub jac: Option<Vec<f64>>,
    /// Column-major `n` by `P` weight derivatives, `None` for the over-identified
    /// and continuous forms.
    pub dw_dbeta: Option<Vec<f64>>,
    /// The generalized-method-of-moments criterion, `Some` only for the
    /// over-identified form.
    pub gmm_obj: Option<f64>,
}

// ---- Binary just-identified -------------------------------------------------

/// The just-identified balancing moment problem for a binary treatment.
///
/// A single propensity model over all units carries the parameter dependence;
/// the moment `sum_i s_i c_i(beta) x_i` balances the covariates and its Jacobian
/// is positive semidefinite, which drives the merit Newton step.
struct JustBlock<'a> {
    covs: &'a [f64],
    n: usize,
    p: usize,
    treat: &'a [i32],
    s: &'a [f64],
    link: Link,
    estimand: CbpsEstimand,
    residual_scale: f64,
    pool: Arc<ThreadPool>,
}

/// The accumulated moment vector and its Jacobian at one parameter value.
struct JustAccum {
    g: Vec<f64>,
    jac: Vec<f64>,
    row: Vec<f64>,
}

impl JustAccum {
    fn zeros(p: usize, want_jac: bool) -> Self {
        Self {
            g: vec![0.0; p],
            jac: if want_jac {
                vec![0.0; p * p]
            } else {
                Vec::new()
            },
            row: vec![0.0; p],
        }
    }
}

impl JustBlock<'_> {
    /// Linear predictor `x_i . beta` for unit `i`.
    fn lin(&self, i: usize, beta: &[f64]) -> f64 {
        let mut acc = 0.0;
        for (j, &b) in beta.iter().enumerate() {
            acc += self.covs[j * self.n + i] * b;
        }
        acc
    }

    /// Fold the moment `g` and, when requested, its Jacobian.
    fn accumulate(&self, beta: &[f64], want_jac: bool) -> JustAccum {
        let p = self.p;
        deterministic_map_reduce(
            &self.pool,
            self.n,
            || JustAccum::zeros(p, want_jac),
            |acc, i| {
                for j in 0..p {
                    acc.row[j] = self.covs[j * self.n + i];
                }
                let eta: f64 = acc.row.iter().zip(beta).map(|(c, b)| c * b).sum();
                let prob = self.link.linkinv(eta);
                let t = f64::from(self.treat[i]);
                let c = self.estimand.bal_factor(prob, t);
                let sc = self.s[i] * c;
                for j in 0..p {
                    acc.g[j] += sc * acc.row[j];
                }
                if want_jac {
                    let mu_eta = self.link.mu_eta(eta);
                    let curv = self.s[i] * self.estimand.bal_factor_deriv_p(prob, t) * mu_eta;
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
                    acc.g[j] += other.g[j];
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

impl EsteqProblem for JustBlock<'_> {
    fn n_params(&self) -> usize {
        self.p
    }

    fn n_units(&self) -> usize {
        self.n
    }

    fn value(&self, _beta: &[f64]) -> Option<f64> {
        // A root-finding problem: the Newton method minimizes the merit
        // 0.5 ||g||^2 rather than an objective.
        None
    }

    fn residual_scale(&self) -> f64 {
        self.residual_scale
    }

    fn gradient(&self, beta: &[f64], g: &mut [f64]) {
        let acc = self.accumulate(beta, false);
        g.copy_from_slice(&acc.g);
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
        g.copy_from_slice(&acc.g);
        for j in 0..self.p {
            for k in 0..self.p {
                *h.rb_mut().get_mut(j, k) = acc.jac[j * self.p + k];
            }
        }
        None
    }

    fn psi(&self, beta: &[f64], mut out: MatMut<'_, f64>) {
        for i in 0..self.n {
            let prob = self.link.linkinv(self.lin(i, beta));
            let t = f64::from(self.treat[i]);
            let c = self.estimand.bal_factor(prob, t);
            for j in 0..self.p {
                *out.rb_mut().get_mut(i, j) = self.s[i] * c * self.covs[j * self.n + i];
            }
        }
    }
}

/// Solve the binary just-identified problem and fill the estimating-equations
/// output.
fn solve_binary_just(inputs: &CbpsInputs<'_>, interrupt: &dyn Fn() -> bool) -> CbpsResult {
    let n = inputs.n;
    let p = inputs.p_mod;
    let pool = get_pool(inputs.threads);

    // Maximum-likelihood warm start: a binomial fit of the treatment on the
    // model covariates gives the Newton method a good starting point.
    let response: Vec<f64> = inputs.treat.iter().map(|&t| f64::from(t)).collect();
    let all: Vec<usize> = (0..n).collect();
    let mut beta = glm::irls(
        inputs.covs_mod,
        n,
        p,
        &all,
        &response,
        inputs.s,
        inputs.link,
        50,
    );

    let block = JustBlock {
        covs: inputs.covs_mod,
        n,
        p,
        treat: inputs.treat,
        s: inputs.s,
        link: inputs.link,
        estimand: inputs.estimand,
        residual_scale: sampling_weight_scale(inputs.s),
        pool: Arc::clone(&pool),
    };

    let opts = SolveOptions {
        max_iter: inputs.max_iter,
        grad_tol: inputs.tol,
        fista_rel_tol: inputs.tol,
    };
    let report = esteq::solve(&block, &mut beta, Solver::Newton, &opts, interrupt);

    // Per-unit propensity, weight, and derivatives at the solution.
    let mut weights = vec![0.0; n];
    let mut ps = vec![0.0; n];
    let mut dw_scale = vec![0.0; n];
    for i in 0..n {
        let eta = block.lin(i, &beta);
        let prob = inputs.link.linkinv(eta);
        let mu_eta = inputs.link.mu_eta(eta);
        let t = f64::from(inputs.treat[i]);
        ps[i] = prob;
        weights[i] = inputs.estimand.weight(prob, t);
        dw_scale[i] = inputs.estimand.weight_deriv_p(prob, t) * mu_eta;
    }

    // psi column j: s_i c_i x_ij. dw column j: (dw / dp) mu_eta x_ij. Both are
    // elementwise per column, so they fill in parallel over columns.
    let mut psi = vec![0.0; n * p];
    let mut dw = vec![0.0; n * p];
    let c_vec: Vec<f64> = (0..n)
        .map(|i| {
            let prob = ps[i];
            let t = f64::from(inputs.treat[i]);
            inputs.s[i] * inputs.estimand.bal_factor(prob, t)
        })
        .collect();
    pool.install(|| {
        psi.par_chunks_mut(n)
            .zip(dw.par_chunks_mut(n))
            .enumerate()
            .for_each(|(j, (psi_col, dw_col))| {
                let covs_col = &inputs.covs_mod[j * n..(j + 1) * n];
                for i in 0..n {
                    psi_col[i] = c_vec[i] * covs_col[i];
                    dw_col[i] = dw_scale[i] * covs_col[i];
                }
            });
    });

    // The Jacobian of the summed estimating functions, positive semidefinite.
    let acc = block.accumulate(&beta, true);

    CbpsResult {
        weights,
        ps,
        coefs: beta,
        converged: report.converged && !report.interrupted,
        interrupted: report.interrupted,
        iterations: report.iterations,
        obj_value: report.grad_norm,
        psi: Some(psi),
        jac: Some(acc.jac),
        dw_dbeta: Some(dw),
        gmm_obj: None,
    }
}

/// The modeled propensity of each unit at `beta`, read off the
/// propensity-model design.
///
/// Both re-evaluation entrypoints start here: the estimating functions turn the
/// propensity into a balancing factor, the weights into the estimand's weight.
fn model_propensities(inputs: &CbpsInputs<'_>, beta: &[f64]) -> Vec<f64> {
    let n = inputs.n;
    let p = inputs.p_mod;
    // The linear predictor sums over the model columns, not over whatever the
    // caller supplied, so a short coefficient vector is a broken call rather than
    // a model on fewer columns. Iterating the slice instead would build the
    // predictor from a prefix of the design and return propensities that look
    // ordinary. The public re-evaluation entrypoints refuse the mismatch before
    // reaching here, so this asserts the invariant they establish.
    assert!(
        beta.len() == p,
        "beta has length {} but the propensity model has {p} coefficient(s)",
        beta.len()
    );
    (0..n)
        .map(|i| {
            let eta: f64 = (0..p).map(|j| inputs.covs_mod[j * n + i] * beta[j]).sum();
            inputs.link.linkinv(eta)
        })
        .collect()
}

/// Check that `beta` carries one coefficient per model covariate, the shared
/// precondition of the re-evaluation entrypoints.
///
/// A short vector would otherwise build the linear predictor from a prefix of the
/// design and a long one would ignore its tail, so the mismatch is returned for
/// the boundary layer to surface rather than left to produce a plausible-looking
/// answer from the wrong model.
fn check_beta(p: usize, beta: &[f64]) -> Result<(), String> {
    if p == 0 || beta.len() != p {
        return Err(format!(
            "beta has length {} but the propensity model has {p} coefficient(s)",
            beta.len()
        ));
    }
    Ok(())
}

/// Re-evaluate the binary just-identified estimating functions at a supplied set
/// of coefficients.
///
/// The estimating-equations container stores `psi` at the solution; a sandwich
/// variance that needs a finite difference must re-evaluate it at perturbed
/// parameters. This recomputes the `n` by `p` matrix from `beta` without
/// solving, matching [`solve_binary_just`]: column `j` is
/// `s_i c_i(beta) x_ij`, the balancing factor times the model covariate.
pub fn eval_psi_binary_just(inputs: &CbpsInputs<'_>, beta: &[f64]) -> Result<Vec<f64>, String> {
    let n = inputs.n;
    let p = inputs.p_mod;
    check_beta(p, beta)?;
    let ps = model_propensities(inputs, beta);
    let mut psi = vec![0.0; n * p];
    for (i, &prob) in ps.iter().enumerate() {
        let t = f64::from(inputs.treat[i]);
        let sc = inputs.s[i] * inputs.estimand.bal_factor(prob, t);
        for j in 0..p {
            psi[j * n + i] = sc * inputs.covs_mod[j * n + i];
        }
    }
    Ok(psi)
}

/// Re-evaluate the binary just-identified balancing weights at a supplied set
/// of coefficients.
///
/// A sandwich variance that treats the weights as a function of the propensity
/// coefficients needs the weight map itself, not only its derivative at the
/// solution. This recomputes the length-`n` weight vector from `beta` without
/// solving, at the same scale [`solve_binary_just`] reports: the estimand's
/// weight function evaluated at the modeled propensity and the unit's treatment
/// indicator.
pub fn eval_weights_binary_just(inputs: &CbpsInputs<'_>, beta: &[f64]) -> Result<Vec<f64>, String> {
    check_beta(inputs.p_mod, beta)?;
    Ok(model_propensities(inputs, beta)
        .into_iter()
        .enumerate()
        .map(|(i, prob)| inputs.estimand.weight(prob, f64::from(inputs.treat[i])))
        .collect())
}

// ---- Binary over-identified GMM --------------------------------------------

/// Weighting-matrix policy for the over-identified criterion.
enum GmmWeighting {
    /// A fixed weighting matrix, the pseudo-inverse of the moment covariance at
    /// a preliminary estimate.
    TwoStep(Vec<f64>),
    /// Continuous updating: the weighting matrix is recomputed at each iterate.
    Continuous,
}

/// The over-identified GMM problem: score conditions stacked on balancing
/// conditions, minimized as `m(beta)' W m(beta)`.
struct GmmProblem<'a> {
    inputs: &'a CbpsInputs<'a>,
    /// Total moments, `p_mod` score plus `p_bal` balance.
    m_total: usize,
    weighting: GmmWeighting,
    pool: Arc<ThreadPool>,
}

/// Per-evaluation GMM statistics: the mean moment, its parameter Jacobian, and
/// the moment covariance for the continuously-updated weighting.
struct GmmAccum {
    m: Vec<f64>,
    grad: Vec<f64>,
    cov: Vec<f64>,
    g_row: Vec<f64>,
    mod_row: Vec<f64>,
    bal_row: Vec<f64>,
}

impl GmmAccum {
    fn zeros(m_total: usize, p: usize, want_cov: bool) -> Self {
        Self {
            m: vec![0.0; m_total],
            grad: vec![0.0; m_total * p],
            cov: if want_cov {
                vec![0.0; m_total * m_total]
            } else {
                Vec::new()
            },
            g_row: vec![0.0; m_total],
            mod_row: vec![0.0; p],
            bal_row: vec![0.0; m_total - p],
        }
    }
}

impl GmmProblem<'_> {
    /// The divisor that turns the per-unit sums below into the normalized
    /// quantities the criterion is built from: the sample size times the average
    /// sampling weight.
    ///
    /// Every accumulated sum here carries exactly one sampling weight per unit,
    /// so dividing by the average weight as well as by the sample size makes the
    /// mean moment, its parameter Jacobian, and the moment covariance all degree
    /// zero in the sampling weights. Expressing the same design in survey-
    /// expansion units then leaves each of them, and so the criterion `m' W m`
    /// and its gradient, at exactly the numbers they held before, which is what
    /// lets a fixed convergence tolerance read the fit rather than the units the
    /// weights happen to be written in. Unit sampling weights average to one
    /// exactly, so the arithmetic there is untouched.
    fn normalizer(&self) -> f64 {
        1.0 / (self.inputs.n as f64 * sampling_weight_scale(self.inputs.s))
    }

    /// Fold the mean moment `m`, the moment-parameter Jacobian `grad` (column
    /// major `m_total` by `p`), and, when requested, the moment covariance, each
    /// normalized by [`GmmProblem::normalizer`].
    fn accumulate(&self, beta: &[f64], want_cov: bool) -> GmmAccum {
        let inputs = self.inputs;
        let n = inputs.n;
        let p = inputs.p_mod;
        let p_bal = inputs.p_bal;
        let m_total = self.m_total;
        let normalizer = self.normalizer();
        let mut acc = deterministic_map_reduce(
            &self.pool,
            n,
            || GmmAccum::zeros(m_total, p, want_cov),
            |acc, i| {
                for j in 0..p {
                    acc.mod_row[j] = inputs.covs_mod[j * n + i];
                }
                for j in 0..p_bal {
                    acc.bal_row[j] = inputs.covs_bal[j * n + i];
                }
                let eta: f64 = acc.mod_row.iter().zip(beta).map(|(c, b)| c * b).sum();
                let prob = clamp_prob(inputs.link.linkinv(eta));
                let mu_eta = inputs.link.mu_eta(eta);
                let t = f64::from(inputs.treat[i]);
                let s = inputs.s[i];

                // Score residual moment on the model covariates, then the
                // balancing moment on the balance covariates.
                let score = t - prob;
                let c = inputs.estimand.bal_factor(prob, t);
                for j in 0..p {
                    acc.g_row[j] = score * acc.mod_row[j];
                }
                let dc = inputs.estimand.bal_factor_deriv_p(prob, t);
                for j in 0..p_bal {
                    acc.g_row[p + j] = c * acc.bal_row[j];
                }
                for a in 0..m_total {
                    acc.m[a] += s * acc.g_row[a];
                }
                // The per-unit moment Jacobian in beta. The score block is
                // -mu_eta x_mod x_mod'; the balance block is dc mu_eta x_bal x_mod'.
                for k in 0..p {
                    let base = k * m_total;
                    let xk = acc.mod_row[k];
                    let score_scale = -mu_eta * xk * s;
                    for j in 0..p {
                        acc.grad[base + j] += score_scale * acc.mod_row[j];
                    }
                    let bal_scale = dc * mu_eta * xk * s;
                    for j in 0..p_bal {
                        acc.grad[base + p + j] += bal_scale * acc.bal_row[j];
                    }
                }
                // The moment covariance accumulates a single sampling weight per
                // unit rather than its square. Any positive semidefinite
                // weighting matrix leaves the generalized-method-of-moments
                // estimator consistent, so this is a choice of efficiency and not
                // of identification: the sampling weights enter as a design
                // measure here, weighting each unit's outer product by the mass
                // it represents, where squaring them would instead treat them as
                // the variance of a heteroskedastic error and change the
                // over-identified point estimates. That choice also gives the
                // covariance the same single power of the sampling weights the
                // moment itself carries, so the normalization applied below
                // clears both at once.
                if want_cov {
                    for a in 0..m_total {
                        let ga = s * acc.g_row[a];
                        for b in 0..m_total {
                            acc.cov[a * m_total + b] += ga * acc.g_row[b];
                        }
                    }
                }
            },
            |acc, other| {
                for a in 0..m_total {
                    acc.m[a] += other.m[a];
                }
                for e in 0..m_total * p {
                    acc.grad[e] += other.grad[e];
                }
                if want_cov {
                    for e in 0..m_total * m_total {
                        acc.cov[e] += other.cov[e];
                    }
                }
            },
        );
        for value in &mut acc.m {
            *value *= normalizer;
        }
        for value in &mut acc.grad {
            *value *= normalizer;
        }
        if want_cov {
            for value in &mut acc.cov {
                *value *= normalizer;
            }
        }
        acc
    }

    /// The weighting matrix at `beta`: fixed for two-step, recomputed for
    /// continuous updating from the supplied covariance. The continuously-updated
    /// covariance is validated at the maximum-likelihood anchor in `solve_binary_
    /// over`, so a decomposition failure at an interior iterate is a genuine bug
    /// rather than a data condition and is surfaced rather than silently zeroed.
    fn weighting_matrix(&self, cov: &[f64]) -> Vec<f64> {
        match &self.weighting {
            GmmWeighting::TwoStep(w) => w.clone(),
            GmmWeighting::Continuous => pseudo_inverse_symmetric(cov, self.m_total)
                .expect("moment covariance decomposition failed at an interior iterate"),
        }
    }
}

impl EsteqProblem for GmmProblem<'_> {
    fn n_params(&self) -> usize {
        self.inputs.p_mod
    }

    fn n_units(&self) -> usize {
        self.inputs.n
    }

    fn value(&self, beta: &[f64]) -> Option<f64> {
        let want_cov = matches!(self.weighting, GmmWeighting::Continuous);
        let acc = self.accumulate(beta, want_cov);
        let w = self.weighting_matrix(&acc.cov);
        Some(quad_form(&w, &acc.m, self.m_total))
    }

    // The criterion is already degree zero in the sampling weights, normalized
    // there by `accumulate`, so it keeps the trait's default residual scale of
    // one and the tolerance is read exactly as written.

    fn hessian_is_exact(&self) -> bool {
        // The Gauss-Newton surrogate below drops the moment curvature, so the
        // steps it generates descend at a linear rate rather than a superlinear
        // one. What the criterion can still be reduced by is therefore bounded by
        // what it can resolve, which is what lets the solver read an unresolvable
        // predicted decrease as the optimum.
        false
    }

    fn gradient(&self, beta: &[f64], grad: &mut [f64]) {
        let p = self.inputs.p_mod;
        let m_total = self.m_total;
        let continuous = matches!(self.weighting, GmmWeighting::Continuous);
        let acc = self.accumulate(beta, continuous);
        let w = self.weighting_matrix(&acc.cov);
        // u = W m, so the base gradient is 2 G' u.
        let u = mat_vec(&w, &acc.m, m_total, m_total);
        for (k, grad_k) in grad.iter_mut().enumerate().take(p) {
            let mut gk = 0.0;
            for (a, &ua) in u.iter().enumerate() {
                gk += acc.grad[k * m_total + a] * ua;
            }
            *grad_k = 2.0 * gk;
        }
        if continuous {
            // The continuously-updated criterion adds the derivative of the
            // weighting matrix: -2/N sum_i (g_i . u)(G_i' u).
            let extra = self.continuous_extra(beta, &u);
            for k in 0..p {
                grad[k] += extra[k];
            }
        }
    }

    /// The Gauss-Newton Hessian `2 G' W G`, the criterion's second derivative with
    /// the moment curvature dropped.
    ///
    /// The exact Hessian of `m' W m` carries the second derivative of the moments
    /// in the parameters, and for continuous updating the second derivative of the
    /// weighting matrix as well. Neither is needed to generate a search direction:
    /// `W` is a pseudo-inverse of a covariance and so positive semidefinite, which
    /// makes this approximation positive semidefinite too and a descent direction
    /// whenever it is invertible, with the solver's ridge escalation covering the
    /// singular case. The line search and the convergence test both run on the
    /// exact criterion and its exact gradient, so the stationary point the solve
    /// reaches is the true minimizer regardless of the approximation. At a
    /// saturated design, where every moment vanishes at the solution, dropping the
    /// curvature term costs nothing at all.
    fn hessian(&self, beta: &[f64], mut h: MatMut<'_, f64>) {
        let p = self.inputs.p_mod;
        let m_total = self.m_total;
        let continuous = matches!(self.weighting, GmmWeighting::Continuous);
        let acc = self.accumulate(beta, continuous);
        let w = self.weighting_matrix(&acc.cov);

        // W G, one column of the moment-parameter Jacobian at a time.
        let mut wg = vec![0.0; m_total * p];
        for k in 0..p {
            let g_col = &acc.grad[k * m_total..(k + 1) * m_total];
            wg[k * m_total..(k + 1) * m_total]
                .copy_from_slice(&mat_vec(&w, g_col, m_total, m_total));
        }
        for j in 0..p {
            for k in 0..p {
                let mut entry = 0.0;
                for a in 0..m_total {
                    entry += acc.grad[j * m_total + a] * wg[k * m_total + a];
                }
                *h.rb_mut().get_mut(j, k) = 2.0 * entry;
            }
        }
    }

    fn psi(&self, _beta: &[f64], _out: MatMut<'_, f64>) {
        // The over-identified form supplies no estimating equations.
    }
}

impl GmmProblem<'_> {
    /// The continuously-updated gradient correction `-2 (g_i . u)(G_i' u)`
    /// summed over units and normalized by [`GmmProblem::normalizer`], for the
    /// fixed contraction `u = W m`.
    ///
    /// The correction differentiates the weighting matrix, and the covariance it
    /// inverts is normalized the same way the moment is, so this term takes the
    /// same divisor as the rest of the gradient.
    fn continuous_extra(&self, beta: &[f64], u: &[f64]) -> Vec<f64> {
        let inputs = self.inputs;
        let n = inputs.n;
        let p = inputs.p_mod;
        let p_bal = inputs.p_bal;
        let normalizer = self.normalizer();
        let mut extra = deterministic_map_reduce(
            &self.pool,
            n,
            || vec![0.0; p],
            |acc, i| {
                let eta: f64 = (0..p).map(|j| inputs.covs_mod[j * n + i] * beta[j]).sum();
                let prob = clamp_prob(inputs.link.linkinv(eta));
                let mu_eta = inputs.link.mu_eta(eta);
                let t = f64::from(inputs.treat[i]);
                let s = inputs.s[i];
                let score = t - prob;
                let c = inputs.estimand.bal_factor(prob, t);
                let dc = inputs.estimand.bal_factor_deriv_p(prob, t);

                // The exposure-model part of u contracted with the model row,
                // shared across every parameter of the score block.
                let mut mod_dot = 0.0;
                for (j, &uj) in u.iter().take(p).enumerate() {
                    mod_dot += inputs.covs_mod[j * n + i] * uj;
                }
                let mut bal_dot = 0.0;
                for j in 0..p_bal {
                    bal_dot += inputs.covs_bal[j * n + i] * u[p + j];
                }
                // q = s (g_i . u), a scalar. The score block contributes
                // score * (x_mod . u) and the balance block c * (x_bal . u_bal).
                let q = s * (score * mod_dot + c * bal_dot);
                // a = G_i' u (no sampling weight; q already carries the single s
                // that matches the moment covariance): the score block gives
                // -mu_eta x_k (x_mod . u), the balance block dc mu_eta x_k
                // (x_bal . u_bal).
                for (k, ak) in acc.iter_mut().enumerate() {
                    let xk = inputs.covs_mod[k * n + i];
                    let val = mu_eta * xk * (-mod_dot + dc * bal_dot);
                    *ak += q * val;
                }
            },
            |acc, other| {
                for k in 0..p {
                    acc[k] += other[k];
                }
            },
        );
        for value in &mut extra {
            *value *= -2.0 * normalizer;
        }
        extra
    }
}

/// Quadratic form `v' A v` for a column-major `m` by `m` matrix.
fn quad_form(a: &[f64], v: &[f64], m: usize) -> f64 {
    let mut acc = 0.0;
    for j in 0..m {
        let vj = v[j];
        for i in 0..m {
            acc += vj * a[j * m + i] * v[i];
        }
    }
    acc
}

/// Matrix-vector product `A v` for a column-major `rows` by `cols` matrix.
fn mat_vec(a: &[f64], v: &[f64], rows: usize, cols: usize) -> Vec<f64> {
    let mut out = vec![0.0; rows];
    for j in 0..cols {
        let vj = v[j];
        for (i, oi) in out.iter_mut().enumerate() {
            *oi += a[j * rows + i] * vj;
        }
    }
    out
}

/// The per-unit propensity and estimand weight at a coefficient vector, shared by
/// the over-identified solve and its non-converged early return.
fn binary_weights_ps(inputs: &CbpsInputs<'_>, beta: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = inputs.n;
    let mut weights = vec![0.0; n];
    let mut ps = vec![0.0; n];
    for i in 0..n {
        let mut eta = 0.0;
        for (j, &b) in beta.iter().enumerate() {
            eta += inputs.covs_mod[j * n + i] * b;
        }
        let prob = clamp_prob(inputs.link.linkinv(eta));
        let t = f64::from(inputs.treat[i]);
        ps[i] = prob;
        weights[i] = inputs.estimand.weight(prob, t);
    }
    (weights, ps)
}

/// Solve the binary over-identified GMM criterion.
fn solve_binary_over(inputs: &CbpsInputs<'_>, interrupt: &dyn Fn() -> bool) -> CbpsResult {
    let n = inputs.n;
    let p = inputs.p_mod;
    let m_total = inputs.p_mod + inputs.p_bal;
    let pool = get_pool(inputs.threads);

    // Preliminary estimate: the maximum-likelihood fit, which also anchors the
    // two-step weighting matrix.
    let response: Vec<f64> = inputs.treat.iter().map(|&t| f64::from(t)).collect();
    let all: Vec<usize> = (0..n).collect();
    let mut beta = glm::irls(
        inputs.covs_mod,
        n,
        p,
        &all,
        &response,
        inputs.s,
        inputs.link,
        50,
    );

    // Validate the moment covariance at the maximum-likelihood anchor, which the
    // two-step policy also uses as its fixed weighting. A failed eigendecomposition
    // leaves the criterion without a usable weighting, so surface it as a
    // non-converged fit rather than letting a zero weighting report an instant
    // spurious optimum at the start.
    let probe = GmmProblem {
        inputs,
        m_total,
        weighting: GmmWeighting::Continuous,
        pool: Arc::clone(&pool),
    };
    let anchor_cov = probe.accumulate(&beta, true).cov;
    let anchor_weighting = match pseudo_inverse_symmetric(&anchor_cov, m_total) {
        Some(w) => w,
        None => {
            let (weights, ps) = binary_weights_ps(inputs, &beta);
            return CbpsResult {
                weights,
                ps,
                coefs: beta,
                converged: false,
                interrupted: false,
                iterations: 0,
                obj_value: f64::NAN,
                psi: None,
                jac: None,
                dw_dbeta: None,
                gmm_obj: Some(f64::NAN),
            };
        }
    };

    let weighting = if inputs.twostep {
        GmmWeighting::TwoStep(anchor_weighting)
    } else {
        GmmWeighting::Continuous
    };
    let problem = GmmProblem {
        inputs,
        m_total,
        weighting,
        pool: Arc::clone(&pool),
    };

    let opts = SolveOptions {
        max_iter: inputs.max_iter,
        grad_tol: inputs.tol,
        fista_rel_tol: inputs.tol,
    };
    // The hybrid rather than the quasi-Newton solver alone. The criterion is a
    // squared quantity, so its gradient reaches the square root of the arithmetic's
    // precision and no further; a run judged on the gradient alone spends its whole
    // budget inside that flat region. The Newton polish carries the stationarity
    // certificate that recognizes the region for what it is, and the L-BFGS warm
    // start still does the work of getting there without the Hessian.
    let report = esteq::solve(
        &problem,
        &mut beta,
        Solver::LbfgsThenNewton,
        &opts,
        interrupt,
    );

    let (weights, ps) = binary_weights_ps(inputs, &beta);
    let criterion = problem.value(&beta).unwrap_or(f64::NAN);

    CbpsResult {
        weights,
        ps,
        coefs: beta,
        converged: report.converged && !report.interrupted,
        interrupted: report.interrupted,
        iterations: report.iterations,
        obj_value: criterion,
        psi: None,
        jac: None,
        dw_dbeta: None,
        gmm_obj: Some(criterion),
    }
}

/// Solve a binary covariate balancing propensity score problem.
pub fn solve(inputs: &CbpsInputs<'_>, interrupt: &dyn Fn() -> bool) -> CbpsResult {
    if inputs.over {
        solve_binary_over(inputs, interrupt)
    } else {
        solve_binary_just(inputs, interrupt)
    }
}

// ---- Categorical just-identified -------------------------------------------

/// Inputs for a categorical covariate balancing propensity score solve.
pub struct CbpsMultiInputs<'a> {
    /// Column-major `n` by `p` propensity-model design.
    pub covs: &'a [f64],
    /// Number of units.
    pub n: usize,
    /// Number of model terms.
    pub p: usize,
    /// Treatment level of each unit, values `0..n_levels`.
    pub treat: &'a [i32],
    /// Number of distinct treatment levels.
    pub n_levels: usize,
    /// Focal level for the effect on the treated; ignored by the average
    /// treatment effect.
    pub focal: usize,
    /// Sampling weights `s_i`.
    pub s: &'a [f64],
    /// Propensity link.
    pub link: Link,
    /// Target population, `Ate` or `Att`.
    pub estimand: CbpsEstimand,
    /// Worker threads for the reductions.
    pub threads: usize,
    /// Maximum solver iterations.
    pub max_iter: usize,
    /// Convergence threshold on the moment sup norm.
    pub tol: f64,
}

/// Solve a categorical covariate balancing propensity score problem.
///
/// For a single covariate model the just-identified categorical estimator
/// coincides with per-level covariate balancing: each level is tilted so its
/// weighted covariate total matches the target population, exactly the
/// construction the tilting solver provides. The estimating-equations output
/// therefore stacks block-diagonally across levels.
pub fn solve_multi(inputs: &CbpsMultiInputs<'_>, interrupt: &dyn Fn() -> bool) -> CbpsResult {
    let estimand = match inputs.estimand {
        CbpsEstimand::Ate => IptEstimand::Ate,
        CbpsEstimand::Att | CbpsEstimand::Atc => IptEstimand::Focal(inputs.focal),
        // The overlap estimand is binary-only; the savvy layer rejects it for a
        // categorical exposure, so it never reaches the categorical solver.
        CbpsEstimand::Ato => unreachable!("the overlap estimand is binary-only"),
    };
    let ipt_inputs = IptInputs {
        covs: inputs.covs,
        n: inputs.n,
        p: inputs.p,
        treat: inputs.treat,
        n_levels: inputs.n_levels,
        s: inputs.s,
        link: inputs.link,
        estimand,
        threads: inputs.threads,
        max_iter: inputs.max_iter,
        tol: inputs.tol,
    };
    let result = ipt::solve(&ipt_inputs, interrupt);
    CbpsResult {
        weights: result.weights,
        ps: result.ps,
        coefs: result.coefs,
        converged: result.converged,
        interrupted: result.interrupted,
        iterations: result.iterations,
        obj_value: result.grad_norm,
        psi: Some(result.psi),
        jac: Some(result.jac),
        dw_dbeta: Some(result.dw_dbeta),
        gmm_obj: None,
    }
}

// ---- Continuous covariate balancing ----------------------------------------

/// Inputs for a continuous-exposure covariate balancing propensity score solve.
pub struct CbpsContInputs<'a> {
    /// Column-major `n` by `p` design, including an intercept column.
    pub covs: &'a [f64],
    /// Number of units.
    pub n: usize,
    /// Number of model terms.
    pub p: usize,
    /// Continuous exposure value of each unit.
    pub expo: &'a [f64],
    /// Sampling weights `s_i`.
    pub s: &'a [f64],
    /// Worker threads for the reductions.
    pub threads: usize,
    /// Maximum solver iterations.
    pub max_iter: usize,
    /// Convergence threshold on the covariance-condition sup norm.
    pub tol: f64,
}

/// The continuous covariance-balancing problem.
///
/// The covariate balancing conditions for a continuous exposure require the
/// weighted exposure mean to match the sample mean and the weighted covariance
/// between the exposure and every covariate to vanish. Positive weights that
/// satisfy these conditions with the least departure from uniformity are the
/// exponential tilt `w_i = exp(-f_i . gamma)`, whose tilt parameters solve the
/// convex dual below. This is the non-parametric covariate balancing solution:
/// the conditions define the estimator, and the tilt is the minimum-divergence
/// weighting that meets them.
struct ContProblem<'a> {
    /// Column-major `n` by `p` feature matrix whose weighted means are the
    /// balancing conditions: the centered exposure for the exposure-mean
    /// condition and the centered-covariate times centered-exposure product for
    /// each covariance condition.
    feat: &'a [f64],
    n: usize,
    p: usize,
    s: &'a [f64],
    pool: Arc<ThreadPool>,
}

/// Accumulated dual sufficient statistics for one tilt value.
struct ContAccum {
    z: f64,
    m: Vec<f64>,
    smat: Vec<f64>,
    row: Vec<f64>,
}

impl ContAccum {
    fn zeros(p: usize, want_smat: bool) -> Self {
        Self {
            z: 0.0,
            m: vec![0.0; p],
            smat: if want_smat {
                vec![0.0; p * p]
            } else {
                Vec::new()
            },
            row: vec![0.0; p],
        }
    }
}

impl ContProblem<'_> {
    /// Fold the tilted feature sums into the dual statistics, filling the second
    /// moment matrix only when the Hessian is needed.
    fn accumulate(&self, gamma: &[f64], want_smat: bool) -> ContAccum {
        let n = self.n;
        let p = self.p;
        deterministic_map_reduce(
            &self.pool,
            n,
            || ContAccum::zeros(p, want_smat),
            |acc, i| {
                for j in 0..p {
                    acc.row[j] = self.feat[j * n + i];
                }
                let lin: f64 = acc.row.iter().zip(gamma).map(|(f, g)| f * g).sum();
                let e = self.s[i] * (-lin).exp();
                acc.z += e;
                for j in 0..p {
                    let efj = e * acc.row[j];
                    acc.m[j] += efj;
                    if want_smat {
                        for k in 0..p {
                            acc.smat[j * p + k] += efj * acc.row[k];
                        }
                    }
                }
            },
            |acc, other| {
                acc.z += other.z;
                for j in 0..p {
                    acc.m[j] += other.m[j];
                }
                if want_smat {
                    for jk in 0..p * p {
                        acc.smat[jk] += other.smat[jk];
                    }
                }
            },
        )
    }
}

impl EsteqProblem for ContProblem<'_> {
    fn n_params(&self) -> usize {
        self.p
    }

    fn n_units(&self) -> usize {
        self.n
    }

    fn value(&self, gamma: &[f64]) -> Option<f64> {
        // The targets are zero, so the dual is just the log partition.
        Some(self.accumulate(gamma, false).z.ln())
    }

    fn gradient(&self, gamma: &[f64], g: &mut [f64]) {
        let acc = self.accumulate(gamma, false);
        for (j, gj) in g.iter_mut().enumerate() {
            *gj = -acc.m[j] / acc.z;
        }
    }

    fn hessian(&self, gamma: &[f64], mut h: MatMut<'_, f64>) {
        let acc = self.accumulate(gamma, true);
        let z = acc.z;
        for j in 0..self.p {
            let mbar_j = acc.m[j] / z;
            for k in 0..self.p {
                let mbar_k = acc.m[k] / z;
                *h.rb_mut().get_mut(j, k) = acc.smat[j * self.p + k] / z - mbar_j * mbar_k;
            }
        }
    }

    fn value_grad_hess(&self, gamma: &[f64], g: &mut [f64], mut h: MatMut<'_, f64>) -> Option<f64> {
        let acc = self.accumulate(gamma, true);
        let z = acc.z;
        for (j, gj) in g.iter_mut().enumerate() {
            let mbar_j = acc.m[j] / z;
            *gj = -mbar_j;
            for k in 0..self.p {
                let mbar_k = acc.m[k] / z;
                *h.rb_mut().get_mut(j, k) = acc.smat[j * self.p + k] / z - mbar_j * mbar_k;
            }
        }
        Some(z.ln())
    }

    fn psi(&self, _gamma: &[f64], _out: MatMut<'_, f64>) {
        // A continuous exposure supplies no estimating equations.
    }
}

/// Solve a continuous-exposure covariate balancing propensity score problem.
pub fn solve_cont(inputs: &CbpsContInputs<'_>, interrupt: &dyn Fn() -> bool) -> CbpsResult {
    let n = inputs.n;
    let p = inputs.p;
    let pool = get_pool(inputs.threads);

    // Weighted exposure and covariate means for centering.
    let sbar: f64 = inputs.s.iter().sum();
    let m_expo: f64 = inputs
        .s
        .iter()
        .zip(inputs.expo)
        .map(|(s, t)| s * t)
        .sum::<f64>()
        / sbar;
    let mut xbar = vec![0.0; p];
    for (j, xbar_j) in xbar.iter_mut().enumerate() {
        let mut acc = 0.0;
        for i in 0..n {
            acc += inputs.s[i] * inputs.covs[j * n + i];
        }
        *xbar_j = acc / sbar;
    }

    // The balancing features. The constant (intercept) column carries the
    // centered exposure, giving the exposure-mean condition; every other column
    // carries the centered covariate times the centered exposure, giving the
    // covariance condition. The features are scaled by the exposure spread so the
    // tilt parameters are commensurate and the dual is well conditioned.
    let expo_sd = {
        let var: f64 = (0..n)
            .map(|i| {
                let de = inputs.expo[i] - m_expo;
                inputs.s[i] * de * de
            })
            .sum::<f64>()
            / sbar;
        var.sqrt().max(1e-8)
    };
    let intercept = (0..p).find(|&j| {
        let first = inputs.covs[j * n];
        first != 0.0 && (0..n).all(|i| (inputs.covs[j * n + i] - first).abs() < 1e-12)
    });
    let mut feat = vec![0.0; n * p];
    for j in 0..p {
        for i in 0..n {
            let de = (inputs.expo[i] - m_expo) / expo_sd;
            feat[j * n + i] = if Some(j) == intercept {
                de
            } else {
                (inputs.covs[j * n + i] - xbar[j]) * de
            };
        }
    }

    let problem = ContProblem {
        feat: &feat,
        n,
        p,
        s: inputs.s,
        pool: Arc::clone(&pool),
    };

    let opts = SolveOptions {
        max_iter: inputs.max_iter,
        grad_tol: inputs.tol,
        fista_rel_tol: inputs.tol,
    };
    let mut gamma = vec![0.0; p];
    let report = esteq::solve(&problem, &mut gamma, Solver::Newton, &opts, interrupt);

    // The tilt weights, normalized so the weighted sample size is preserved.
    let mut weights = vec![0.0; n];
    for (i, wi) in weights.iter_mut().enumerate() {
        let mut lin = 0.0;
        for j in 0..p {
            lin += feat[j * n + i] * gamma[j];
        }
        *wi = (-lin).exp();
    }
    let wsum: f64 = inputs.s.iter().zip(&weights).map(|(s, w)| s * w).sum();
    let scale = sbar / wsum;
    for w in &mut weights {
        *w *= scale;
    }

    CbpsResult {
        weights,
        ps: vec![f64::NAN; n],
        coefs: gamma,
        converged: report.converged && !report.interrupted,
        interrupted: report.interrupted,
        iterations: report.iterations,
        obj_value: report.grad_norm,
        psi: None,
        jac: None,
        dw_dbeta: None,
        gmm_obj: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn no_interrupt() -> impl Fn() -> bool {
        || false
    }

    /// A saturated two-cell binary design: an intercept and a 0/1 cell
    /// indicator, so the fit recovers exact empirical cell propensities. Cell
    /// z = 0 has 4 units, 2 treated (p = 1/2); cell z = 1 has 6 units, 2 treated
    /// (p = 1/3).
    fn saturated_design() -> (Vec<f64>, Vec<i32>, usize) {
        let z = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        let treat = vec![1, 1, 0, 0, 1, 1, 0, 0, 0, 0];
        let n = z.len();
        let mut covs = vec![1.0; n];
        covs.extend_from_slice(&z);
        (covs, treat, n)
    }

    fn run_binary(estimand: CbpsEstimand, threads: usize) -> CbpsResult {
        let (covs, treat, n) = saturated_design();
        let s = vec![1.0; n];
        let inputs = CbpsInputs {
            covs_mod: &covs,
            covs_bal: &covs,
            n,
            p_mod: 2,
            p_bal: 2,
            treat: &treat,
            s: &s,
            link: Link::Logit,
            estimand,
            over: false,
            twostep: true,
            threads,
            max_iter: 200,
            tol: 1e-12,
        };
        solve(&inputs, &no_interrupt())
    }

    /// The largest absolute balancing moment at the returned weights, computed
    /// directly from the estimand's weight form: the weighted covariate total of
    /// the two arms must agree (or match the target arm).
    fn balance_gap(result: &CbpsResult, covs: &[f64], treat: &[i32], n: usize, p: usize) -> f64 {
        let mut gap = 0.0_f64;
        for j in 0..p {
            let mut treated = 0.0;
            let mut control = 0.0;
            for i in 0..n {
                let x = covs[j * n + i];
                if treat[i] == 1 {
                    treated += result.weights[i] * x;
                } else {
                    control += result.weights[i] * x;
                }
            }
            gap = gap.max((treated - control).abs());
        }
        gap
    }

    #[test]
    fn ate_recovers_inverse_cell_propensities() {
        let result = run_binary(CbpsEstimand::Ate, 1);
        assert!(result.converged);
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

    #[test]
    fn ate_balances_the_covariates_exactly() {
        let (covs, treat, n) = saturated_design();
        let result = run_binary(CbpsEstimand::Ate, 1);
        assert!(balance_gap(&result, &covs, &treat, n, 2) < 1e-7);
    }

    #[test]
    fn att_keeps_treated_at_one_and_tilts_controls() {
        let result = run_binary(CbpsEstimand::Att, 1);
        assert!(result.converged);
        // Treated stay at one; each control's weight is the cell odds of
        // treatment: cell 0 gives 1, cell 1 gives (1/3)/(2/3) = 1/2.
        let expected = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5];
        for (i, &want) in expected.iter().enumerate() {
            assert!(
                (result.weights[i] - want).abs() < 1e-8,
                "weight[{i}] = {} expected {want}",
                result.weights[i]
            );
        }
    }

    #[test]
    fn ato_matches_the_overlap_tilt_closed_form() {
        let result = run_binary(CbpsEstimand::Ato, 1);
        assert!(result.converged);
        // Overlap weights: treated by (1 - p), control by p, at the empirical
        // cell propensity. Cell 0 (p = 1/2): both 1/2. Cell 1 (p = 1/3): treated
        // 2/3, control 1/3.
        let expected = [
            0.5,
            0.5,
            0.5,
            0.5,
            2.0 / 3.0,
            2.0 / 3.0,
            1.0 / 3.0,
            1.0 / 3.0,
            1.0 / 3.0,
            1.0 / 3.0,
        ];
        for (i, &want) in expected.iter().enumerate() {
            assert!(
                (result.weights[i] - want).abs() < 1e-8,
                "weight[{i}] = {} expected {want}",
                result.weights[i]
            );
        }
        let (covs, treat, n) = saturated_design();
        assert!(balance_gap(&result, &covs, &treat, n, 2) < 1e-7);
    }

    #[test]
    fn atc_keeps_controls_at_one() {
        let result = run_binary(CbpsEstimand::Atc, 1);
        assert!(result.converged);
        // Controls stay at one; treated weight is the cell odds against
        // treatment (1 - p) / p. Cell 0: 1. Cell 1: (2/3)/(1/3) = 2.
        let expected = [1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0];
        for (i, &want) in expected.iter().enumerate() {
            assert!(
                (result.weights[i] - want).abs() < 1e-8,
                "weight[{i}] = {} expected {want}",
                result.weights[i]
            );
        }
    }

    #[test]
    fn psi_columns_sum_to_zero() {
        for estimand in [
            CbpsEstimand::Ate,
            CbpsEstimand::Att,
            CbpsEstimand::Atc,
            CbpsEstimand::Ato,
        ] {
            let result = run_binary(estimand, 1);
            let psi = result.psi.as_ref().expect("just-identified has psi");
            let n = result.weights.len();
            for col in 0..result.coefs.len() {
                let s: f64 = (0..n).map(|i| psi[col * n + i]).sum();
                assert!(
                    s.abs() < 1e-7,
                    "estimand {estimand:?} column {col} sums to {s}"
                );
            }
        }
    }

    #[test]
    fn weight_derivative_matches_finite_difference() {
        // Perturb one coefficient and read the ATE weight of an interior unit
        // against the stored analytic derivative.
        let (covs, treat, n) = saturated_design();
        let result = run_binary(CbpsEstimand::Ate, 1);
        let beta = &result.coefs;
        let dw = result.dw_dbeta.as_ref().expect("just-identified has dw");
        let unit = 4usize; // a treated unit in cell 1
        let link = Link::Logit;
        let eta: f64 = (0..2).map(|j| covs[j * n + unit] * beta[j]).sum();
        let t = f64::from(treat[unit]);
        let h = 1e-6;
        for j in 0..2 {
            let x = covs[j * n + unit];
            let w = |d: f64| CbpsEstimand::Ate.weight(link.linkinv(eta + d * x), t);
            let fd = (w(h) - w(-h)) / (2.0 * h);
            let analytic = dw[j * n + unit];
            assert!(
                (fd - analytic).abs() < 1e-5,
                "col {j}: fd {fd} analytic {analytic}"
            );
        }
    }

    // Re-evaluating psi at the solved coefficients reproduces the solve's own
    // psi, and a central finite difference of the column sums reproduces the
    // analytic Jacobian, across every binary estimand.
    #[test]
    fn eval_psi_matches_solve_and_jacobian() {
        let (covs, treat, n) = saturated_design();
        let s = vec![1.0; n];
        let p = 2;
        for estimand in [
            CbpsEstimand::Ate,
            CbpsEstimand::Att,
            CbpsEstimand::Atc,
            CbpsEstimand::Ato,
        ] {
            let inputs = CbpsInputs {
                covs_mod: &covs,
                covs_bal: &covs,
                n,
                p_mod: p,
                p_bal: p,
                treat: &treat,
                s: &s,
                link: Link::Logit,
                estimand,
                over: false,
                twostep: true,
                threads: 1,
                max_iter: 200,
                tol: 1e-12,
            };
            let result = solve(&inputs, &no_interrupt());
            let stored = result.psi.as_ref().expect("just-identified has psi");
            let jac = result.jac.as_ref().expect("just-identified has jac");

            let recomputed =
                eval_psi_binary_just(&inputs, &result.coefs).expect("beta length matches");
            for (a, b) in stored.iter().zip(&recomputed) {
                assert!((a - b).abs() < 1e-10, "psi mismatch: {a} vs {b}");
            }

            let eps = 1e-6;
            for col in 0..p {
                let mut up = result.coefs.clone();
                let mut down = result.coefs.clone();
                up[col] += eps;
                down[col] -= eps;
                let psi_up = eval_psi_binary_just(&inputs, &up).expect("beta length matches");
                let psi_down = eval_psi_binary_just(&inputs, &down).expect("beta length matches");
                for row in 0..p {
                    let cs_up: f64 = (0..n).map(|i| psi_up[row * n + i]).sum();
                    let cs_down: f64 = (0..n).map(|i| psi_down[row * n + i]).sum();
                    let fd = (cs_up - cs_down) / (2.0 * eps);
                    let analytic = jac[col * p + row];
                    assert!(
                        (fd - analytic).abs() < 1e-4,
                        "estimand {estimand:?} jac[{row},{col}]: fd {fd} analytic {analytic}"
                    );
                }
            }
        }
    }

    // Re-evaluating the weights at the solved coefficients reproduces the
    // solve's own weight vector, and a central finite difference of that vector
    // reproduces the analytic weight derivative, across every binary estimand.
    #[test]
    fn eval_weights_binary_just_matches_solve_and_weight_jacobian() {
        let (covs, treat, n) = saturated_design();
        let s = vec![1.0; n];
        let p = 2;
        for estimand in [
            CbpsEstimand::Ate,
            CbpsEstimand::Att,
            CbpsEstimand::Atc,
            CbpsEstimand::Ato,
        ] {
            let inputs = CbpsInputs {
                covs_mod: &covs,
                covs_bal: &covs,
                n,
                p_mod: p,
                p_bal: p,
                treat: &treat,
                s: &s,
                link: Link::Logit,
                estimand,
                over: false,
                twostep: true,
                threads: 1,
                max_iter: 200,
                tol: 1e-12,
            };
            let result = solve(&inputs, &no_interrupt());
            assert!(result.converged, "estimand {estimand:?} did not converge");
            let dw = result.dw_dbeta.as_ref().expect("just-identified has dw");

            let recomputed =
                eval_weights_binary_just(&inputs, &result.coefs).expect("beta length matches");
            assert_eq!(recomputed.len(), n);
            for (i, (a, b)) in result.weights.iter().zip(&recomputed).enumerate() {
                assert!(
                    (a - b).abs() < 1e-12,
                    "estimand {estimand:?} weight[{i}]: solve {a} versus eval {b}"
                );
            }

            let eps = 1e-6;
            for col in 0..p {
                let mut up = result.coefs.clone();
                let mut down = result.coefs.clone();
                up[col] += eps;
                down[col] -= eps;
                let w_up = eval_weights_binary_just(&inputs, &up).expect("beta length matches");
                let w_down = eval_weights_binary_just(&inputs, &down).expect("beta length matches");
                for i in 0..n {
                    let fd = (w_up[i] - w_down[i]) / (2.0 * eps);
                    let analytic = dw[col * n + i];
                    assert!(
                        (fd - analytic).abs() < 1e-5,
                        "estimand {estimand:?} dw[{i},{col}]: fd {fd} analytic {analytic}"
                    );
                }
            }
        }
    }

    // Non-unit sampling weights move the solved coefficients, so the weight
    // re-evaluation and its finite difference are checked again on the
    // sampling-weighted fixture. The weight map itself does not read the
    // sampling weights, which is exactly why the check belongs here: an
    // implementation that folded them in would agree with the solve output
    // under unit weights and disagree here.
    #[test]
    fn eval_weights_binary_just_under_sampling_weights() {
        let covs = vec![
            // intercept column
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, //
            // covariate column
            -1.2, 0.4, 0.9, -0.5, 1.5, -0.8, 0.2, 1.1,
        ];
        let treat = vec![1, 0, 1, 0, 1, 0, 1, 0];
        let s = vec![0.7, 1.3, 0.5, 1.8, 1.1, 0.9, 1.4, 0.6];
        let n = 8;
        let p = 2;
        let inputs = CbpsInputs {
            covs_mod: &covs,
            covs_bal: &covs,
            n,
            p_mod: p,
            p_bal: p,
            treat: &treat,
            s: &s,
            link: Link::Logit,
            estimand: CbpsEstimand::Ate,
            over: false,
            twostep: true,
            threads: 1,
            max_iter: 200,
            tol: 1e-12,
        };
        let result = solve(&inputs, &no_interrupt());
        assert!(result.converged, "moment norm {}", result.obj_value);
        let dw = result.dw_dbeta.as_ref().expect("just-identified has dw");

        let recomputed =
            eval_weights_binary_just(&inputs, &result.coefs).expect("beta length matches");
        for (i, (a, b)) in result.weights.iter().zip(&recomputed).enumerate() {
            assert!(
                (a - b).abs() < 1e-12,
                "weight[{i}]: solve {a} versus eval {b}"
            );
        }

        let eps = 1e-6;
        for col in 0..p {
            let mut up = result.coefs.clone();
            let mut down = result.coefs.clone();
            up[col] += eps;
            down[col] -= eps;
            let w_up = eval_weights_binary_just(&inputs, &up).expect("beta length matches");
            let w_down = eval_weights_binary_just(&inputs, &down).expect("beta length matches");
            for i in 0..n {
                let fd = (w_up[i] - w_down[i]) / (2.0 * eps);
                let analytic = dw[col * n + i];
                assert!(
                    (fd - analytic).abs() < 1e-5,
                    "dw[{i},{col}]: fd {fd} analytic {analytic}"
                );
            }
        }
    }

    // A coefficient vector shorter than the model design used to build the linear
    // predictor from however many coefficients it held, so a caller that passed
    // the wrong length got propensities from a subset of the columns rather than a
    // complaint. Both re-evaluation entrypoints report the mismatch, matching how
    // inverse probability tilting guards its own.
    #[test]
    fn eval_psi_rejects_wrong_length_beta() {
        let (covs, treat, n) = saturated_design();
        let s = vec![1.0; n];
        let inputs = CbpsInputs {
            covs_mod: &covs,
            covs_bal: &covs,
            n,
            p_mod: 2,
            p_bal: 2,
            treat: &treat,
            s: &s,
            link: Link::Logit,
            estimand: CbpsEstimand::Ate,
            over: false,
            twostep: false,
            threads: 1,
            max_iter: 0,
            tol: 0.0,
        };
        let err = eval_psi_binary_just(&inputs, &[0.0]).unwrap_err();
        assert!(err.contains("beta has length 1"), "message was: {err}");
        let err = eval_psi_binary_just(&inputs, &[0.0, 0.0, 0.0]).unwrap_err();
        assert!(err.contains("beta has length 3"), "message was: {err}");
    }

    #[test]
    fn eval_weights_rejects_wrong_length_beta() {
        let (covs, treat, n) = saturated_design();
        let s = vec![1.0; n];
        let inputs = CbpsInputs {
            covs_mod: &covs,
            covs_bal: &covs,
            n,
            p_mod: 2,
            p_bal: 2,
            treat: &treat,
            s: &s,
            link: Link::Logit,
            estimand: CbpsEstimand::Ate,
            over: false,
            twostep: false,
            threads: 1,
            max_iter: 0,
            tol: 0.0,
        };
        let err = eval_weights_binary_just(&inputs, &[0.0]).unwrap_err();
        assert!(err.contains("beta has length 1"), "message was: {err}");
    }

    // The propensity map is the shared inner step, and it used to truncate the
    // linear predictor silently. The public entrypoints refuse a wrong length
    // first, so this pins that the inner step fails loudly rather than returning a
    // partial model if a future caller reaches it directly.
    #[test]
    #[should_panic(expected = "beta has length 1")]
    fn model_propensities_fails_loudly_on_short_beta() {
        let (covs, treat, n) = saturated_design();
        let s = vec![1.0; n];
        let inputs = CbpsInputs {
            covs_mod: &covs,
            covs_bal: &covs,
            n,
            p_mod: 2,
            p_bal: 2,
            treat: &treat,
            s: &s,
            link: Link::Logit,
            estimand: CbpsEstimand::Ate,
            over: false,
            twostep: false,
            threads: 1,
            max_iter: 0,
            tol: 0.0,
        };
        let _ = model_propensities(&inputs, &[0.0]);
    }

    #[test]
    fn solve_is_deterministic_across_thread_counts() {
        let one = run_binary(CbpsEstimand::Ate, 1);
        for threads in [2, 4, 8] {
            let many = run_binary(CbpsEstimand::Ate, threads);
            for i in 0..one.weights.len() {
                assert_eq!(one.weights[i].to_bits(), many.weights[i].to_bits());
            }
            for i in 0..one.coefs.len() {
                assert_eq!(one.coefs[i].to_bits(), many.coefs[i].to_bits());
            }
        }
    }

    // The balancing moment is a sampling-weighted total, so expressing the same
    // design in survey-expansion units multiplies the residual by the expansion
    // factor while leaving the solution untouched. The verdict must follow the
    // solution rather than the units.
    #[test]
    fn convergence_verdict_is_invariant_to_the_sampling_weight_scale() {
        let (covs, treat, n) = saturated_design();
        let solve_at = |scale: f64| {
            let s = vec![scale; n];
            let inputs = CbpsInputs {
                covs_mod: &covs,
                covs_bal: &covs,
                n,
                p_mod: 2,
                p_bal: 2,
                treat: &treat,
                s: &s,
                link: Link::Logit,
                estimand: CbpsEstimand::Ate,
                over: false,
                twostep: true,
                threads: 1,
                max_iter: 200,
                tol: 1e-10,
            };
            solve(&inputs, &no_interrupt())
        };
        let base = solve_at(1.0);
        let expanded = solve_at(1e6);

        for i in 0..n {
            assert!(
                (base.weights[i] - expanded.weights[i]).abs() < 1e-9,
                "weight[{i}]: unscaled {} versus expanded {}",
                base.weights[i],
                expanded.weights[i]
            );
        }
        assert!(base.converged, "the unscaled fit must converge");
        assert!(
            expanded.converged,
            "the expanded fit solved to the same weights but reported \
             a moment sup norm of {} against tolerance 1e-10",
            expanded.obj_value
        );
    }

    #[test]
    fn interrupt_is_surfaced() {
        let (covs, treat, n) = saturated_design();
        let s = vec![1.0; n];
        let inputs = CbpsInputs {
            covs_mod: &covs,
            covs_bal: &covs,
            n,
            p_mod: 2,
            p_bal: 2,
            treat: &treat,
            s: &s,
            link: Link::Logit,
            estimand: CbpsEstimand::Ate,
            over: false,
            twostep: true,
            threads: 1,
            max_iter: 200,
            tol: 1e-12,
        };
        let result = solve(&inputs, &|| true);
        assert!(result.interrupted);
        assert!(!result.converged);
    }

    #[test]
    fn every_link_recovers_the_saturated_propensities() {
        let (covs, treat, n) = saturated_design();
        let s = vec![1.0; n];
        for link in [Link::Logit, Link::Probit, Link::Cloglog] {
            let inputs = CbpsInputs {
                covs_mod: &covs,
                covs_bal: &covs,
                n,
                p_mod: 2,
                p_bal: 2,
                treat: &treat,
                s: &s,
                link,
                estimand: CbpsEstimand::Ate,
                over: false,
                twostep: true,
                threads: 1,
                max_iter: 200,
                tol: 1e-12,
            };
            let result = solve(&inputs, &no_interrupt());
            assert!(result.converged, "link {link:?} did not converge");
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

    // ---- Over-identified GMM ------------------------------------------------

    fn run_over(
        twostep: bool,
        threads: usize,
        covs: &[f64],
        treat: &[i32],
        n: usize,
        p: usize,
    ) -> CbpsResult {
        let s = vec![1.0; n];
        let inputs = CbpsInputs {
            covs_mod: covs,
            covs_bal: covs,
            n,
            p_mod: p,
            p_bal: p,
            treat,
            s: &s,
            link: Link::Logit,
            estimand: CbpsEstimand::Ate,
            over: true,
            twostep,
            threads,
            max_iter: 500,
            tol: 1e-10,
        };
        solve(&inputs, &no_interrupt())
    }

    #[test]
    fn over_identified_criterion_vanishes_on_a_saturated_design() {
        // On a saturated design the maximum-likelihood fit satisfies both the
        // score and the balancing conditions, so every moment is zero and the
        // two-step criterion attains its floor of zero.
        let (covs, treat, n) = saturated_design();
        let result = run_over(true, 1, &covs, &treat, n, 2);
        let obj = result
            .gmm_obj
            .expect("over-identified reports the criterion");
        assert!(obj >= 0.0);
        assert!(obj < 1e-10, "criterion {obj} should vanish");
        assert!(result.weights.iter().all(|&w| w >= 0.0));
        assert!(result.psi.is_none());
    }

    #[test]
    fn the_gmm_criterion_minimizes_under_the_newton_solver() {
        // Every build reaches the Newton method on this criterion: it polishes the
        // quasi-Newton pass where the basin backend is compiled in and drives both
        // stages where it is not, so the problem has to supply a Hessian the Newton
        // step can use. Driving Newton directly exercises that path in either
        // build: on the saturated design every moment can be driven to zero, so the
        // criterion reaches its floor.
        let (covs, treat, n) = saturated_design();
        let s = vec![1.0; n];
        let inputs = CbpsInputs {
            covs_mod: &covs,
            covs_bal: &covs,
            n,
            p_mod: 2,
            p_bal: 2,
            treat: &treat,
            s: &s,
            link: Link::Logit,
            estimand: CbpsEstimand::Ate,
            over: true,
            twostep: true,
            threads: 1,
            max_iter: 200,
            tol: 1e-10,
        };
        let pool = get_pool(1);
        let m_total = inputs.p_mod + inputs.p_bal;
        let mut beta = vec![0.0; inputs.p_mod];

        // The two-step weighting anchored at the start point, the same
        // construction `solve_binary_over` uses.
        let probe = GmmProblem {
            inputs: &inputs,
            m_total,
            weighting: GmmWeighting::Continuous,
            pool: Arc::clone(&pool),
        };
        let anchor = pseudo_inverse_symmetric(&probe.accumulate(&beta, true).cov, m_total)
            .expect("the moment covariance at the start point decomposes");
        let problem = GmmProblem {
            inputs: &inputs,
            m_total,
            weighting: GmmWeighting::TwoStep(anchor),
            pool,
        };

        let opts = SolveOptions {
            max_iter: 200,
            grad_tol: 1e-10,
            fista_rel_tol: 1e-10,
        };
        let report = esteq::solve(&problem, &mut beta, Solver::Newton, &opts, &no_interrupt());
        let criterion = problem.value(&beta).expect("the criterion has a value");
        assert!(
            criterion < 1e-10,
            "criterion {criterion} after {} Newton iterations",
            report.iterations
        );
    }

    // A non-saturated design with a continuous covariate whose treated and
    // untreated supports interleave: relative to `separated_design`, the units at
    // x = -0.2 and x = 0.4 trade assignments, so neither group's covariate values
    // sit entirely above or entirely below the other's. The score and
    // balancing conditions still cannot both hold, so the weighting matrix
    // matters, but the logistic maximum-likelihood fit that anchors the solve now
    // exists and both weighting policies reach a minimizer.
    //
    // This is the general-purpose over-identified fixture, and the overlap is what
    // qualifies it: a test written on it reads the behavior it means to read
    // rather than a failure to solve. Where the failure itself is the subject, use
    // `separated_design`.
    fn overlapping_cov_design() -> (Vec<f64>, Vec<i32>, usize) {
        let x = [
            -1.3, -0.7, -0.2, 0.4, 0.9, 1.5, -1.1, -0.4, 0.1, 0.6, 1.2, 1.8, -0.9, 0.3, 1.0, -0.5,
        ];
        let treat = vec![0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0];
        let n = x.len();
        let mut covs = vec![1.0; n];
        covs.extend_from_slice(&x);
        (covs, treat, n)
    }

    // The same covariate values under a perfectly separated assignment: every
    // treated unit sits at a positive covariate value and every untreated unit at
    // a negative one, so no logistic propensity model fits it. The
    // maximum-likelihood coefficient on the covariate grows without bound, the
    // warm-start iteration leaves the finite range within ten passes, and the
    // over-identified solve is left anchored at the all-zero fallback, a point
    // that carries no information about the design. The continuously-updated arm
    // barely moves from there and stops with a gradient near a tenth whatever
    // iteration budget it is given.
    //
    // That is what the fixture is for. A design the criterion cannot be solved on
    // is the one place to pin that nothing in the convergence machinery, neither
    // the stalled-iterate certificate nor a tolerance read at the problem's own
    // scale, ever reports a non-stationary stop as convergence.
    fn separated_design() -> (Vec<f64>, Vec<i32>, usize) {
        let x = [
            -1.3, -0.7, -0.2, 0.4, 0.9, 1.5, -1.1, -0.4, 0.1, 0.6, 1.2, 1.8, -0.9, 0.3, 1.0, -0.5,
        ];
        let treat = vec![0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0];
        let n = x.len();
        let mut covs = vec![1.0; n];
        covs.extend_from_slice(&x);
        (covs, treat, n)
    }

    // A stop short of a stationary point is not convergence, however the verdict
    // is reached. Both halves are asserted here: the separation that makes the
    // stopping point non-stationary, and the verdict that has to follow from it.
    #[test]
    fn a_stall_on_a_separated_design_is_not_reported_as_convergence() {
        let (covs, treat, n) = separated_design();
        for i in 0..n {
            assert_eq!(
                treat[i] == 1,
                covs[n + i] > 0.0,
                "unit {i} breaks the separation the fixture is built on"
            );
        }

        let s = vec![1.0; n];
        let inputs = CbpsInputs {
            covs_mod: &covs,
            covs_bal: &covs,
            n,
            p_mod: 2,
            p_bal: 2,
            treat: &treat,
            s: &s,
            link: Link::Logit,
            estimand: CbpsEstimand::Ate,
            over: true,
            twostep: false,
            threads: 1,
            max_iter: 500,
            tol: 1e-10,
        };
        let result = solve(&inputs, &no_interrupt());

        // The criterion rebuilt at the returned coefficients under the same
        // continuously-updated weighting the solve ran on, so the gradient read
        // here is the one the solve was asked to drive to zero.
        let problem = GmmProblem {
            inputs: &inputs,
            m_total: 4,
            weighting: GmmWeighting::Continuous,
            pool: get_pool(1),
        };
        let mut grad = vec![0.0; 2];
        problem.gradient(&result.coefs, &mut grad);
        let grad_norm = grad.iter().fold(0.0_f64, |acc, g| acc.max(g.abs()));
        assert!(
            grad_norm > 1e-3,
            "the fixture must leave the solve away from a stationary point, \
             gradient sup norm {grad_norm}"
        );
        assert!(
            !result.converged,
            "a stop at a gradient of {grad_norm} is not convergence"
        );
    }

    #[test]
    fn two_step_and_continuous_updating_differ() {
        let (covs, treat, n) = overlapping_cov_design();
        let two = run_over(true, 1, &covs, &treat, n, 2);
        let cue = run_over(false, 1, &covs, &treat, n, 2);
        assert!(two.gmm_obj.unwrap().is_finite());
        assert!(cue.gmm_obj.unwrap().is_finite());
        let diff = (0..n)
            .map(|i| (two.weights[i] - cue.weights[i]).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            diff > 1e-8,
            "two-step and continuous updating gave the same weights"
        );
    }

    #[test]
    fn over_identified_is_deterministic_across_thread_counts() {
        let (covs, treat, n) = overlapping_cov_design();
        let one = run_over(true, 1, &covs, &treat, n, 2);
        for threads in [2, 4] {
            let many = run_over(true, threads, &covs, &treat, n, 2);
            for i in 0..n {
                assert_eq!(one.weights[i].to_bits(), many.weights[i].to_bits());
            }
        }
    }

    // The over-identified criterion is degree zero in the sampling weights: the
    // mean moment, its parameter Jacobian, and the moment covariance the
    // weighting matrix inverts are each divided by the average sampling weight,
    // which leaves every one of them unchanged when the same design is expressed
    // in survey-expansion units. The criterion and its gradient are then the same
    // numbers at any expansion, not merely proportional to each other, so a fixed
    // convergence tolerance reads the fit rather than the units the weights are
    // written in. This pins the structure directly, at a fixed parameter vector
    // where no solver trajectory intervenes. Both weighting policies share it, so
    // both are checked.
    #[test]
    fn the_gmm_criterion_is_degree_zero_in_the_sampling_weights() {
        let (covs, treat, n) = gmm_gradient_design();
        let base_s = gmm_gradient_sampling_weights();
        let factor = 1e6;
        let beta = [0.15, -0.3];
        let p = 2;
        let m_total = 2 * p;

        for twostep in [true, false] {
            let mut criteria = Vec::new();
            let mut gradients = Vec::new();
            for expansion in [1.0, factor] {
                let s: Vec<f64> = base_s.iter().map(|si| si * expansion).collect();
                let inputs = CbpsInputs {
                    covs_mod: &covs,
                    covs_bal: &covs,
                    n,
                    p_mod: p,
                    p_bal: p,
                    treat: &treat,
                    s: &s,
                    link: Link::Logit,
                    estimand: CbpsEstimand::Ate,
                    over: true,
                    twostep,
                    threads: 1,
                    max_iter: 1,
                    tol: 1e-12,
                };
                let pool = get_pool(1);
                let weighting = if twostep {
                    let probe = GmmProblem {
                        inputs: &inputs,
                        m_total,
                        weighting: GmmWeighting::Continuous,
                        pool: Arc::clone(&pool),
                    };
                    let cov = probe.accumulate(&beta, true).cov;
                    GmmWeighting::TwoStep(
                        pseudo_inverse_symmetric(&cov, m_total).expect("covariance decomposes"),
                    )
                } else {
                    GmmWeighting::Continuous
                };
                let problem = GmmProblem {
                    inputs: &inputs,
                    m_total,
                    weighting,
                    pool,
                };
                let mut grad = vec![0.0; p];
                problem.gradient(&beta, &mut grad);
                criteria.push(problem.value(&beta).expect("the criterion has a value"));
                gradients.push(grad);
            }

            // The bound is round-off on a million-fold expansion, far below the
            // factor-of-a-million discrepancy a degree-one criterion would show.
            assert!(
                (criteria[1] - criteria[0]).abs() <= 1e-9 * criteria[0].abs(),
                "two_step = {twostep} criterion moved from {} to {} under a \
                 {factor}-fold expansion",
                criteria[0],
                criteria[1]
            );
            for (k, unscaled) in gradients[0].iter().enumerate() {
                assert!(
                    (gradients[1][k] - unscaled).abs() <= 1e-9 * unscaled.abs(),
                    "two_step = {twostep} gradient[{k}] moved from {unscaled} to {} \
                     under a {factor}-fold expansion",
                    gradients[1][k]
                );
            }
        }
    }

    // The verdict the criterion above earns, end to end through the solve. The
    // two-step policy is the one this design holds the tolerance on at both
    // expansions, so it is the policy that can pin a verdict here: an arm that
    // misses the tolerance at one of the two readings would pin the backend's
    // difficulty rather than the scale. The continuously-updated criterion's own
    // degree in the sampling weights is pinned above, and both policies are
    // exercised end to end by the R suite.
    #[test]
    fn the_over_identified_verdict_is_invariant_to_the_sampling_weight_scale() {
        let (covs, treat, n) = overlapping_cov_design();
        let solve_at = |scale: f64| {
            let s = vec![scale; n];
            let inputs = CbpsInputs {
                covs_mod: &covs,
                covs_bal: &covs,
                n,
                p_mod: 2,
                p_bal: 2,
                treat: &treat,
                s: &s,
                link: Link::Logit,
                estimand: CbpsEstimand::Ate,
                over: true,
                twostep: true,
                threads: 1,
                max_iter: 500,
                tol: 1e-10,
            };
            solve(&inputs, &no_interrupt())
        };
        let base = solve_at(1.0);
        let expanded = solve_at(1e6);

        assert!(base.converged, "the unscaled fit must converge");
        assert!(
            expanded.converged,
            "the expanded fit reported a criterion of {} against tolerance 1e-10",
            expanded.gmm_obj.unwrap_or(f64::NAN)
        );
        // The reported criterion is the same number at either expansion, not the
        // same number times the expansion factor.
        let base_obj = base.gmm_obj.expect("the criterion is reported");
        let expanded_obj = expanded.gmm_obj.expect("the criterion is reported");
        assert!(
            (expanded_obj - base_obj).abs() <= 1e-6 * base_obj.abs().max(1e-12),
            "criterion {base_obj} unscaled against {expanded_obj} expanded"
        );
        for i in 0..n {
            assert!(
                (base.weights[i] - expanded.weights[i]).abs() < 1e-6,
                "weight[{i}]: unscaled {} versus expanded {}",
                base.weights[i],
                expanded.weights[i]
            );
        }
    }

    #[test]
    fn over_identified_interrupt_is_surfaced() {
        let (covs, treat, n) = overlapping_cov_design();
        let s = vec![1.0; n];
        let inputs = CbpsInputs {
            covs_mod: &covs,
            covs_bal: &covs,
            n,
            p_mod: 2,
            p_bal: 2,
            treat: &treat,
            s: &s,
            link: Link::Logit,
            estimand: CbpsEstimand::Ate,
            over: true,
            twostep: true,
            threads: 1,
            max_iter: 500,
            tol: 1e-10,
        };
        let result = solve(&inputs, &|| true);
        assert!(result.interrupted);
        assert!(!result.converged);
    }

    /// A confounded binary design of `n` units in the shape the R suite's
    /// `sim_binary` produces: two standard normal covariates, an exposure drawn
    /// from a logit model on them, and sampling weights spread over the interval
    /// from a half to three halves. A small linear congruential generator with a
    /// Box-Muller transform supplies the draws, so the design is byte identical on
    /// every run and needs no dependency.
    ///
    /// Size is the point. At a few hundred units the over-identified criterion
    /// stops being resolvable well before its gradient reaches the default
    /// tolerance: the smallest parameter displacement the criterion can tell apart
    /// from no displacement already leaves a gradient near a billionth, so a solve
    /// that has arrived at the numerical minimizer still reads a gradient far above
    /// the tolerance. The sixteen-unit fixtures above are too small to show it.
    fn confounded_binary_design(n: usize) -> (Vec<f64>, Vec<i32>, Vec<f64>) {
        let mut state: u64 = 0x2545F4914F6CDD1D;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            // Shifted off zero so the logarithm in the Box-Muller transform stays
            // finite.
            ((state >> 11) as f64 + 0.5) / ((1u64 << 53) as f64)
        };
        let mut normal = || {
            let u1: f64 = next();
            let u2: f64 = next();
            (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
        };
        let x1: Vec<f64> = (0..n).map(|_| normal()).collect();
        let x2: Vec<f64> = (0..n).map(|_| normal()).collect();
        let mut covs = vec![1.0; n];
        covs.extend_from_slice(&x1);
        covs.extend_from_slice(&x2);
        let mut treat = vec![0_i32; n];
        let mut s = vec![0.0; n];
        for i in 0..n {
            let eta = 0.8 * x1[i] - 0.6 * x2[i];
            treat[i] = i32::from(next() < 1.0 / (1.0 + (-eta).exp()));
            s[i] = 0.5 + next();
        }
        (covs, treat, s)
    }

    // A solve that has reached the criterion's numerical minimizer has converged,
    // whatever the gradient sup norm there happens to read. The over-identified
    // criterion is a squared quantity, so the smallest parameter displacement its
    // own arithmetic resolves leaves a gradient near the square root of machine
    // epsilon times the criterion's scale, orders of magnitude above the default
    // tolerance. Without a stationarity certificate the verdict on this design
    // records how far the backend happened to crawl before its budget ran out
    // rather than where it arrived, and the two weighting policies disagree about
    // a fit they both reach to nine digits.
    //
    // The tolerance is the default, not loosened: the certificate is what makes
    // the verdict readable, and it is asserted for both policies.
    #[test]
    fn the_over_identified_verdict_certifies_the_numerical_minimizer() {
        let n = 500;
        let (covs, treat, s) = confounded_binary_design(n);
        for twostep in [true, false] {
            let inputs = CbpsInputs {
                covs_mod: &covs,
                covs_bal: &covs,
                n,
                p_mod: 3,
                p_bal: 3,
                treat: &treat,
                s: &s,
                link: Link::Logit,
                estimand: CbpsEstimand::Ate,
                over: true,
                twostep,
                threads: 1,
                max_iter: 200,
                tol: 1e-10,
            };
            let result = solve(&inputs, &no_interrupt());
            assert!(
                result.converged,
                "two_step = {twostep} stopped short of a verdict after {} \
                 iterations at a criterion of {}",
                result.iterations,
                result.gmm_obj.unwrap_or(f64::NAN)
            );
            assert!(
                result.iterations < 200,
                "two_step = {twostep} burnt its cap"
            );
        }
    }

    // ---- Categorical --------------------------------------------------------

    fn categorical_design() -> (Vec<f64>, Vec<i32>, usize) {
        // Three cells crossed with three treatment levels, saturated so balance
        // holds exactly. Intercept plus a cell covariate.
        let z = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 0.0, 1.0, 2.0];
        let treat = vec![0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2];
        let n = z.len();
        let mut covs = vec![1.0; n];
        covs.extend_from_slice(&z);
        (covs, treat, n)
    }

    #[test]
    fn categorical_ate_balances_and_sums_psi_to_zero() {
        let (covs, treat, n) = categorical_design();
        let s = vec![1.0; n];
        let inputs = CbpsMultiInputs {
            covs: &covs,
            n,
            p: 2,
            treat: &treat,
            n_levels: 3,
            focal: 0,
            s: &s,
            link: Link::Logit,
            estimand: CbpsEstimand::Ate,
            threads: 1,
            max_iter: 200,
            tol: 1e-12,
        };
        let result = solve_multi(&inputs, &no_interrupt());
        assert!(result.converged);
        assert!(result.weights.iter().all(|&w| w >= 0.0));
        // Every level's weighted covariate mean matches the whole sample.
        let sample_mean = |j: usize| (0..n).map(|i| covs[j * n + i]).sum::<f64>() / n as f64;
        for level in 0..3 {
            for j in 0..2 {
                let mut wsum = 0.0;
                let mut wx = 0.0;
                for i in 0..n {
                    if treat[i] as usize == level {
                        wsum += result.weights[i];
                        wx += result.weights[i] * covs[j * n + i];
                    }
                }
                assert!(
                    (wx / wsum - sample_mean(j)).abs() < 1e-6,
                    "level {level} covariate {j} not balanced"
                );
            }
        }
        let psi = result.psi.as_ref().unwrap();
        for col in 0..result.coefs.len() {
            let sum: f64 = (0..n).map(|i| psi[col * n + i]).sum();
            assert!(sum.abs() < 1e-6, "psi column {col} sums to {sum}");
        }
    }

    // ---- Continuous ---------------------------------------------------------

    fn continuous_exposure_design() -> (Vec<f64>, Vec<f64>, usize) {
        // Exposure correlated with a covariate through a linear term, plus noise
        // that is independent of the covariate, so the unweighted covariance is
        // non-trivial yet a moderate reweighting can decorrelate the two.
        let x = [
            -1.5, -1.1, -0.8, -0.3, 0.1, 0.5, 0.9, 1.2, 1.6, -0.6, 0.2, 0.7, -1.2, 0.4, 1.0, -0.9,
        ];
        let noise = [
            0.3, -0.5, 0.7, -0.2, 0.4, -0.6, 0.1, 0.5, -0.4, 0.2, -0.3, 0.6, -0.7, 0.35, -0.15,
            0.25,
        ];
        let expo: Vec<f64> = x
            .iter()
            .zip(noise.iter())
            .map(|(xi, ni)| 0.6 * xi + ni)
            .collect();
        let n = x.len();
        let mut covs = vec![1.0; n];
        covs.extend_from_slice(&x);
        (covs, expo, n)
    }

    /// Largest absolute weighted covariance between the exposure and a covariate.
    fn weighted_cov_gap(covs: &[f64], expo: &[f64], weights: &[f64], n: usize, p: usize) -> f64 {
        let wsum: f64 = weights.iter().sum();
        let ebar: f64 = weights.iter().zip(expo).map(|(w, e)| w * e).sum::<f64>() / wsum;
        let mut gap = 0.0_f64;
        for j in 1..p {
            let xbar: f64 = weights
                .iter()
                .enumerate()
                .map(|(i, w)| w * covs[j * n + i])
                .sum::<f64>()
                / wsum;
            let cov: f64 = (0..n)
                .map(|i| weights[i] * (covs[j * n + i] - xbar) * (expo[i] - ebar))
                .sum::<f64>()
                / wsum;
            gap = gap.max(cov.abs());
        }
        gap
    }

    #[test]
    fn continuous_drives_the_weighted_covariance_to_zero() {
        let (covs, expo, n) = continuous_exposure_design();
        let s = vec![1.0; n];
        // The unweighted covariance is well away from zero to start.
        let uniform = vec![1.0; n];
        assert!(weighted_cov_gap(&covs, &expo, &uniform, n, 2) > 0.1);

        let inputs = CbpsContInputs {
            covs: &covs,
            n,
            p: 2,
            expo: &expo,
            s: &s,
            threads: 1,
            max_iter: 500,
            tol: 1e-12,
        };
        let result = solve_cont(&inputs, &no_interrupt());
        assert!(result.converged, "moment norm {}", result.obj_value);
        assert!(result.weights.iter().all(|&w| w > 0.0));
        assert!(weighted_cov_gap(&covs, &expo, &result.weights, n, 2) < 1e-6);
    }

    #[test]
    fn continuous_is_deterministic_across_thread_counts() {
        let (covs, expo, n) = continuous_exposure_design();
        let s = vec![1.0; n];
        let run = |threads: usize| {
            let inputs = CbpsContInputs {
                covs: &covs,
                n,
                p: 2,
                expo: &expo,
                s: &s,
                threads,
                max_iter: 500,
                tol: 1e-12,
            };
            solve_cont(&inputs, &no_interrupt())
        };
        let one = run(1);
        for threads in [2, 4] {
            let many = run(threads);
            for i in 0..n {
                assert_eq!(one.weights[i].to_bits(), many.weights[i].to_bits());
            }
        }
    }

    /// An eight-unit binary design with an intercept and one covariate, used for
    /// the criterion's analytic checks at a fixed parameter vector rather than for
    /// a solve.
    fn gmm_gradient_design() -> (Vec<f64>, Vec<i32>, usize) {
        let covs = vec![
            // intercept column
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, //
            // covariate column
            -1.2, 0.4, 0.9, -0.5, 1.5, -0.8, 0.2, 1.1,
        ];
        let treat = vec![1, 0, 1, 0, 1, 0, 1, 0];
        (covs, treat, 8)
    }

    /// Non-unit sampling weights for that design, so a term that reads them the
    /// wrong number of times shows up where unit weights would hide it.
    fn gmm_gradient_sampling_weights() -> Vec<f64> {
        vec![0.7, 1.3, 0.5, 1.8, 1.1, 0.9, 1.4, 0.6]
    }

    /// The generalized-method-of-moments gradient must match a central finite
    /// difference of the criterion for both weighting policies, under non-unit
    /// sampling weights. The two-step criterion holds its weighting fixed, so the
    /// gradient is the plain moment-Jacobian contraction; the continuously-updated
    /// criterion differentiates the weighting too, and a sampling-weight
    /// double-count in that correction shows up here but not under unit weights.
    #[test]
    fn gmm_gradient_matches_finite_differences() {
        let (covs, treat, n) = gmm_gradient_design();
        let s = gmm_gradient_sampling_weights();
        let p = 2;
        let m_total = 2 * p;
        let beta = [0.15, -0.3];

        for &twostep in &[true, false] {
            let inputs = CbpsInputs {
                covs_mod: &covs,
                covs_bal: &covs,
                n,
                p_mod: p,
                p_bal: p,
                treat: &treat,
                s: &s,
                link: Link::Logit,
                estimand: CbpsEstimand::Ate,
                over: true,
                twostep,
                threads: 1,
                max_iter: 1,
                tol: 1e-12,
            };
            let pool = get_pool(1);
            let weighting = if twostep {
                // The two-step weighting is fixed at the evaluation point, so the
                // finite difference sees a constant matrix.
                let probe = GmmProblem {
                    inputs: &inputs,
                    m_total,
                    weighting: GmmWeighting::Continuous,
                    pool: Arc::clone(&pool),
                };
                let cov = probe.accumulate(&beta, true).cov;
                GmmWeighting::TwoStep(
                    pseudo_inverse_symmetric(&cov, m_total).expect("covariance decomposes"),
                )
            } else {
                GmmWeighting::Continuous
            };
            let problem = GmmProblem {
                inputs: &inputs,
                m_total,
                weighting,
                pool: Arc::clone(&pool),
            };

            let mut analytic = vec![0.0; p];
            problem.gradient(&beta, &mut analytic);

            let h = 1e-6;
            for k in 0..p {
                let mut bp = beta;
                let mut bm = beta;
                bp[k] += h;
                bm[k] -= h;
                let fp = problem.value(&bp).expect("criterion is defined");
                let fm = problem.value(&bm).expect("criterion is defined");
                let fd = (fp - fm) / (2.0 * h);
                // The bound sits far below the sampling-weight double-count this
                // guards against (order 1e-2) and well above central-difference
                // noise (order 1e-6).
                assert!(
                    (analytic[k] - fd).abs() < 1e-4,
                    "twostep={twostep} k={k}: analytic {} versus finite difference {}",
                    analytic[k],
                    fd,
                );
            }
        }
    }
}
