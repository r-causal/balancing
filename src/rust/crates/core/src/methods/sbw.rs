//! Stable balancing weights: minimum-dispersion weights under approximate
//! covariate balance.
//!
//! Among all weightings that hold each reweighted group's covariate means inside
//! a tolerance band, stable balancing weights select the one of least dispersion.
//! The default objective minimizes the sum of squared weights, so for a fixed
//! per-group total it minimizes the weight variance, following Zubizarreta. The
//! problem is a strictly convex quadratic program: a positive-definite diagonal
//! objective over a simplex-type constraint set, so both quadratic-program
//! backends accept it, unlike the indefinite energy form.
//!
//! The constraint set reuses the shared assembler: a box on each weight with a
//! minimum-weight floor and pinned zero-sampling-weight units, one group-sum row
//! per reweighted group, and one moment row per covariate holding the group's
//! weighted mean within `target +/- tolerance`. The tolerance is applied at its
//! full width against the target for every estimand, following Zubizarreta's
//! formulation, so the achieved arm-to-target standardized mean difference matches
//! the tolerance the method is asked for. For a focal estimand this is the
//! external reference's own single held-fixed band; for the average treatment
//! effect it is a strict superset of the reference's set, which additionally pins
//! the pair average and so binds each arm at half width. In both cases the
//! reference weights are feasible in this band, so their dispersion bounds the
//! minimum-dispersion solution from above, the comparison the objective-level
//! parity check relies on.
//!
//! For a focal estimand the focal group is held at unit weight and the non-focal
//! groups are pulled to the focal group's means; for the average treatment effect
//! every group is pulled to the shared target. A continuous exposure keeps the
//! same minimum-dispersion objective under a single total-sum constraint and one
//! bounded weighted-correlation row per covariate.

use crate::qp::{
    Convexity, PMat, QpError, QpOptions, QpSolution, QpSpec, QpStatus, RoutedSolution, solve_psd,
};

use super::qp_balance::{ConstraintBuilder, ZERO_SW, expand_and_floor, group_normalized};

/// The dispersion norm the objective minimizes.
///
/// Each norm measures the spread of the group-normalized weights around the
/// uniform baseline of one. `L2` minimizes the sum of squared deviations through
/// a diagonal quadratic program; `L1` and `Linf` minimize the summed and the
/// largest absolute deviation through a linear program in auxiliary variables.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SbwNorm {
    /// Minimize the sum of squared weights (minimum variance).
    L2,
    /// Minimize the sum of absolute weight deviations.
    L1,
    /// Minimize the largest weight deviation.
    Linf,
}

/// The estimand a stable balancing solve targets.
#[derive(Debug, Clone, Copy)]
pub enum SbwEstimand {
    /// Average treatment effect: every group is pulled to the shared target.
    Ate,
    /// A focal-group effect: the non-focal groups are pulled to the focal group,
    /// whose units keep unit weight.
    Focal { focal: usize },
}

/// A reason a stable balancing solve cannot proceed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SbwError {
    /// An explicitly requested backend is not compiled into this build; the field
    /// names the missing backend.
    BackendUnavailable(&'static str),
}

impl std::fmt::Display for SbwError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SbwError::BackendUnavailable(name) => write!(
                f,
                "the `{name}` backend is not compiled into this build; rebuild with the `qp-clarabel` feature or choose another backend"
            ),
        }
    }
}

/// Outcome of a stable balancing solve, carrying the same diagnostic fields as
/// the other quadratic-program methods.
#[derive(Debug, Clone)]
pub struct SbwResult {
    /// Balancing weights, one per unit, in the original unit order.
    pub weights: Vec<f64>,
    /// Constraint dual variables from the quadratic-program solve.
    pub duals: Vec<f64>,
    /// Whether the solve reached a solved status.
    pub converged: bool,
    /// Whether the solve stopped because a user interrupt was pending.
    pub interrupted: bool,
    /// Iterations performed by the backend.
    pub iterations: usize,
    /// Objective value at the solution.
    pub objective: f64,
    /// The backend that produced the result, `"osqp"` or `"clarabel"`, carried to
    /// R as the solver identity.
    pub backend: &'static str,
    /// Whether the automatic clarabel fallback engaged after an osqp
    /// primal-infeasibility certificate.
    pub fell_back: bool,
    /// The backend's terminal status name.
    pub status: &'static str,
    /// Primal residual from the backend.
    pub pri_res: f64,
    /// Dual residual from the backend.
    pub dua_res: f64,
}

/// Inputs for a binary or multi-category stable balancing solve.
pub struct SbwDiscreteInputs<'a> {
    /// Number of units.
    pub n: usize,
    /// Zero-based exposure level of each unit; a negative level excludes the unit.
    pub levels: &'a [i32],
    /// Number of exposure levels.
    pub n_levels: usize,
    /// The estimand.
    pub estimand: SbwEstimand,
    /// Sampling weights.
    pub s: &'a [f64],
    /// Minimum allowable weight.
    pub min_weight: f64,
    /// The dispersion norm to minimize.
    pub norm: SbwNorm,
    /// Column-major `n` by `q` moment-constraint covariates, standardized on the
    /// R side; empty for no moment constraints.
    pub moment_covs: &'a [f64],
    /// Number of moment-constraint columns `q`.
    pub n_moments: usize,
    /// Target weighted mean for each moment column.
    pub targets: &'a [f64],
    /// Tolerance band width for each moment column, applied at full width against
    /// the target on each side.
    pub tols: &'a [f64],
    /// Quadratic-program tuning.
    pub qp: QpOptions,
}

/// Assemble the doubled diagonal quadratic term for the `L2` objective over the
/// active variables. The loss is the sum of squared weights, so the loss matrix
/// is the identity and the doubled term the spec stores is `2 I`.
fn l2_pmat(nvar: usize) -> PMat {
    PMat::Diagonal(vec![2.0; nvar])
}

/// The uniform baseline the absolute-deviation norms measure spread against. Each
/// reweighted group is normalized so its group-normalized weights average one, so
/// a weight of one is the no-reweighting reference the `L1` and `Linf` objectives
/// linearize `|w - reference|` around, matching the reference implementation.
const SBW_REFERENCE: f64 = 1.0;

/// A small quadratic weight on the balancing variables of the absolute-deviation
/// norms. The pure L1 and Linf objectives are linear, and the alternating-direction
/// backend can stall on a linear program without strong convexity; this ridge
/// restores it. It is small enough that the minimizer stays a minimum-deviation
/// weighting to well within the objective-parity tolerance, and it is carried into
/// the reported objective so the golden comparison remains exact. The value doubles
/// into the spec's quadratic block by the doubling convention, so each weight
/// variable's diagonal entry is `2 * SBW_LP_RIDGE`.
const SBW_LP_RIDGE: f64 = 1e-4;

/// Add the box row for each weight variable, floored at `min_weight` or pinned to
/// one where a zero sampling weight excludes the unit from reweighting. This
/// mirrors [`ConstraintBuilder::add_box`] but writes only the first `nvar`
/// variables, leaving the trailing auxiliary variables to their own bounds.
fn add_weight_box(builder: &mut ConstraintBuilder, min_weight: f64, pinned: &[bool]) {
    for (i, &is_pinned) in pinned.iter().enumerate() {
        let (l, u) = if is_pinned {
            (1.0, 1.0)
        } else {
            (min_weight, f64::INFINITY)
        };
        builder.add_sparse_row(&[(i, 1.0)], l, u);
    }
}

/// Append the group-sum and moment rows shared by every norm. Each row's
/// coefficients span only the `nvar` weight variables; the auxiliary columns of
/// the absolute-deviation norms carry no balance coefficients, so a dense row of
/// length `nvar` addresses exactly the weight block whatever the total variable
/// count.
fn add_balance_rows(
    builder: &mut ConstraintBuilder,
    inputs: &SbwDiscreteInputs<'_>,
    active: &[usize],
    group_levels: &[usize],
    swnt: &[f64],
) {
    let n = inputs.n;
    // One group-sum row per reweighted group, fixing the group-normalized total
    // to one so each group's mean weight is one.
    for &t in group_levels {
        let coeffs: Vec<f64> = active
            .iter()
            .map(|&i| {
                if inputs.levels[i] as usize == t {
                    swnt[i]
                } else {
                    0.0
                }
            })
            .collect();
        builder.add_dense_row(&coeffs, 1.0, 1.0);
    }

    // One moment row per group per covariate, holding the group's weighted mean
    // within the full tolerance band around the target. Applying the tolerance at
    // its full width against the target keeps the reference weights feasible in
    // this band, so their dispersion bounds the minimum-dispersion solution the
    // comparison is measured against.
    for &t in group_levels {
        for c in 0..inputs.n_moments {
            let coeffs: Vec<f64> = active
                .iter()
                .map(|&i| {
                    if inputs.levels[i] as usize == t {
                        inputs.moment_covs[c * n + i] * swnt[i]
                    } else {
                        0.0
                    }
                })
                .collect();
            let band = inputs.tols[c].abs();
            builder.add_dense_row(&coeffs, inputs.targets[c] - band, inputs.targets[c] + band);
        }
    }
}

/// Append the auxiliary box rows and the conversion rows that pin each auxiliary
/// deviation to the absolute weight departure from the reference.
///
/// For `L1` there is one auxiliary `t_j` per active weight, bounding `|w_j - 1|`
/// from below, and the objective sums them. For `Linf` a single auxiliary bounds
/// the largest departure across every active weight, and the objective is that
/// auxiliary alone. In both cases the pair of conversion rows
/// `w_j + a >= 1` and `w_j - a <= 1` forces `a >= |w_j - 1|`, and minimizing the
/// auxiliary drives it to equality.
fn add_deviation_rows(builder: &mut ConstraintBuilder, nvar: usize, per_weight_aux: bool) {
    let aux = |j: usize| if per_weight_aux { nvar + j } else { nvar };
    // Each auxiliary is non-negative. A shared auxiliary is bounded once.
    let n_aux = if per_weight_aux { nvar } else { 1 };
    for k in 0..n_aux {
        builder.add_sparse_row(&[(nvar + k, 1.0)], 0.0, f64::INFINITY);
    }
    for j in 0..nvar {
        let a = aux(j);
        builder.add_sparse_row(&[(j, 1.0), (a, 1.0)], SBW_REFERENCE, f64::INFINITY);
        builder.add_sparse_row(&[(j, 1.0), (a, -1.0)], f64::NEG_INFINITY, SBW_REFERENCE);
    }
}

/// Assemble the quadratic-program specification for a discrete solve under the
/// requested norm. The `L2` form is a diagonal quadratic program over the weight
/// variables alone; the `L1` and `Linf` forms are linear programs that add
/// auxiliary deviation variables with a zero quadratic block, which stays
/// positive semidefinite so the routed backend policy applies unchanged.
fn assemble_discrete_spec(
    inputs: &SbwDiscreteInputs<'_>,
    active: &[usize],
    group_levels: &[usize],
    swnt: &[f64],
    pinned: &[bool],
) -> QpSpec {
    let nvar = active.len();
    match inputs.norm {
        SbwNorm::L2 => {
            let mut builder = ConstraintBuilder::new(nvar);
            builder.add_box(inputs.min_weight, pinned);
            add_balance_rows(&mut builder, inputs, active, group_levels, swnt);
            let (m, indptr, indices, values, l, u) = builder.finish();
            QpSpec {
                n: nvar,
                m,
                p: l2_pmat(nvar),
                q: vec![0.0; nvar],
                a_indptr: indptr,
                a_indices: indices,
                a_values: values,
                l,
                u,
                convexity: Convexity::Psd,
            }
        }
        SbwNorm::L1 | SbwNorm::Linf => {
            let per_weight_aux = matches!(inputs.norm, SbwNorm::L1);
            let n_aux = if per_weight_aux { nvar } else { 1 };
            let total = nvar + n_aux;
            let mut builder = ConstraintBuilder::new(total);
            add_weight_box(&mut builder, inputs.min_weight, pinned);
            add_balance_rows(&mut builder, inputs, active, group_levels, swnt);
            add_deviation_rows(&mut builder, nvar, per_weight_aux);
            let (m, indptr, indices, values, l, u) = builder.finish();
            // The objective is linear: minimize the sum of the auxiliary
            // deviations (L1) or the single shared deviation (Linf). A small ridge
            // on the weight block restores the strong convexity the backend needs;
            // the auxiliary block stays at zero. The whole term is positive
            // semidefinite by construction.
            let mut q = vec![0.0; total];
            for qi in q.iter_mut().skip(nvar) {
                *qi = 1.0;
            }
            let mut diag = vec![0.0; total];
            for di in diag.iter_mut().take(nvar) {
                *di = 2.0 * SBW_LP_RIDGE;
            }
            QpSpec {
                n: total,
                m,
                p: PMat::Diagonal(diag),
                q,
                a_indptr: indptr,
                a_indices: indices,
                a_values: values,
                l,
                u,
                convexity: Convexity::Psd,
            }
        }
    }
}

/// The active variables, the group rows that carry constraints, and the focal
/// level for a focal estimand. For the average treatment effect every present
/// unit is a variable and every level is constrained; for a focal estimand only
/// the non-focal units are variables and only the non-focal levels are
/// constrained.
fn active_layout(
    n: usize,
    levels: &[i32],
    n_levels: usize,
    estimand: SbwEstimand,
) -> (Vec<usize>, Vec<usize>) {
    match estimand {
        SbwEstimand::Ate => {
            let active: Vec<usize> = (0..n).filter(|&i| levels[i] >= 0).collect();
            let groups: Vec<usize> = (0..n_levels).collect();
            (active, groups)
        }
        SbwEstimand::Focal { focal } => {
            let active: Vec<usize> = (0..n)
                .filter(|&i| levels[i] >= 0 && levels[i] as usize != focal)
                .collect();
            let groups: Vec<usize> = (0..n_levels).filter(|&t| t != focal).collect();
            (active, groups)
        }
    }
}

/// Solve a binary or multi-category stable balancing problem.
pub fn solve_discrete(
    inputs: &SbwDiscreteInputs<'_>,
    interrupt: &dyn Fn() -> bool,
) -> Result<SbwResult, SbwError> {
    let n = inputs.n;
    let s_norm = group_normalized(inputs.s, inputs.levels, inputs.n_levels);

    // Level sizes and each unit's own group-normalization value, so the group-sum
    // row reads as a mean of one and the moment row as a weighted mean.
    let mut n_t = vec![0usize; inputs.n_levels];
    for &g in inputs.levels {
        if g >= 0 {
            n_t[g as usize] += 1;
        }
    }
    let swnt: Vec<f64> = inputs
        .levels
        .iter()
        .enumerate()
        .map(|(i, &g)| {
            if g >= 0 && n_t[g as usize] > 0 {
                s_norm[i] / n_t[g as usize] as f64
            } else {
                0.0
            }
        })
        .collect();

    let (active, group_levels) = active_layout(n, inputs.levels, inputs.n_levels, inputs.estimand);

    let pinned: Vec<bool> = active
        .iter()
        .map(|&i| inputs.s[i].abs() < ZERO_SW)
        .collect();

    let spec = assemble_discrete_spec(inputs, &active, &group_levels, &swnt, &pinned);

    let routed = route(&spec, &inputs.qp, interrupt)?;
    // The weight variables lead the decision vector; any auxiliary deviations
    // trail them and are dropped by taking the first entry per active unit.
    let weights = expand_and_floor(n, &active, &routed.solution.x, inputs.min_weight);
    Ok(package(weights, &routed))
}

/// Inputs for a continuous-exposure stable balancing solve.
pub struct SbwContInputs<'a> {
    /// Number of units.
    pub n: usize,
    /// Continuous exposure of each unit.
    pub treat: &'a [f64],
    /// Column-major `n` by `q` covariate columns, standardized on the R side, held
    /// in weighted correlation with the exposure within a tolerance.
    pub covs: &'a [f64],
    /// Number of covariate columns `q`.
    pub n_covs: usize,
    /// Sampling weights.
    pub s: &'a [f64],
    /// Minimum allowable weight.
    pub min_weight: f64,
    /// The dispersion norm to minimize.
    pub norm: SbwNorm,
    /// Tolerance on the weighted exposure-covariate correlation for each column.
    pub tols: &'a [f64],
    /// Quadratic-program tuning.
    pub qp: QpOptions,
}

/// Reliability-weighted variance of a vector, matching the denominator the
/// distance transforms use, so a bounded weighted product of standardized columns
/// reads as a bounded correlation.
fn weighted_variance(x: &[f64], w: &[f64]) -> f64 {
    let mut sw = 0.0;
    let mut sw2 = 0.0;
    let mut swx = 0.0;
    let mut swxx = 0.0;
    for (&xi, &wi) in x.iter().zip(w) {
        sw += wi;
        sw2 += wi * wi;
        swx += wi * xi;
        swxx += wi * xi * xi;
    }
    if sw <= 0.0 {
        return 0.0;
    }
    let mean = swx / sw;
    let denom = 1.0 - sw2 / (sw * sw);
    if denom > 0.0 {
        ((swxx / sw - mean * mean) / denom).max(0.0)
    } else {
        0.0
    }
}

/// Solve a continuous-exposure stable balancing problem.
pub fn solve_cont(
    inputs: &SbwContInputs<'_>,
    interrupt: &dyn Fn() -> bool,
) -> Result<SbwResult, SbwError> {
    let n = inputs.n;
    let nf = n as f64;

    // Sampling weights scaled to sum to n, the normalization the constraints
    // assume, so the total-sum row reads as a mean weight of one.
    let s_sum: f64 = inputs.s.iter().sum();
    let s: Vec<f64> = if s_sum > 0.0 {
        inputs.s.iter().map(|&si| si * nf / s_sum).collect()
    } else {
        vec![1.0; n]
    };

    // The exposure centered and scaled to unit weighted variance, so a bounded
    // weighted mean of the exposure-covariate product is a bounded correlation.
    let a_mean = inputs
        .s
        .iter()
        .zip(inputs.treat)
        .map(|(&wi, &ti)| wi * ti)
        .sum::<f64>()
        / s_sum.max(f64::MIN_POSITIVE);
    let a_var = weighted_variance(inputs.treat, inputs.s);
    let a_sd = if a_var > 0.0 { a_var.sqrt() } else { 1.0 };
    let treat_std: Vec<f64> = inputs.treat.iter().map(|&t| (t - a_mean) / a_sd).collect();

    let pinned: Vec<bool> = inputs.s.iter().map(|&si| si.abs() < ZERO_SW).collect();

    // The auxiliary layout mirrors the discrete solve: L2 keeps only the n weight
    // variables under a diagonal quadratic term, while L1 and Linf append their
    // deviation variables and switch to the linear objective.
    let (per_weight_aux, n_aux) = match inputs.norm {
        SbwNorm::L2 => (false, 0usize),
        SbwNorm::L1 => (true, n),
        SbwNorm::Linf => (false, 1),
    };
    let total = n + n_aux;

    let mut builder = ConstraintBuilder::new(total);
    match inputs.norm {
        SbwNorm::L2 => builder.add_box(inputs.min_weight, &pinned),
        SbwNorm::L1 | SbwNorm::Linf => add_weight_box(&mut builder, inputs.min_weight, &pinned),
    }

    // The single total-sum row fixes the weighted total to n. The auxiliary
    // columns carry no total-sum coefficient, so a length-n row addresses the
    // weight block alone.
    builder.add_dense_row(&s, nf, nf);

    // One bounded weighted-correlation row per covariate. The average divides by n,
    // matching the reference's denominator so a boundary-feasible instance is not
    // rejected by a spuriously tighter (n - 1)/n band.
    let denom = nf.max(1.0);
    for c in 0..inputs.n_covs {
        let coeffs: Vec<f64> = (0..n)
            .map(|i| inputs.covs[c * n + i] * treat_std[i] * s[i] / denom)
            .collect();
        let tol = inputs.tols[c].abs();
        builder.add_dense_row(&coeffs, -tol, tol);
    }

    let (p, q) = match inputs.norm {
        SbwNorm::L2 => (l2_pmat(n), vec![0.0; n]),
        SbwNorm::L1 | SbwNorm::Linf => {
            add_deviation_rows(&mut builder, n, per_weight_aux);
            let mut q = vec![0.0; total];
            for qi in q.iter_mut().skip(n) {
                *qi = 1.0;
            }
            let mut diag = vec![0.0; total];
            for di in diag.iter_mut().take(n) {
                *di = 2.0 * SBW_LP_RIDGE;
            }
            (PMat::Diagonal(diag), q)
        }
    };

    let (m, indptr, indices, values, l, u) = builder.finish();
    let spec = QpSpec {
        n: total,
        m,
        p,
        q,
        a_indptr: indptr,
        a_indices: indices,
        a_values: values,
        l,
        u,
        convexity: Convexity::Psd,
    };

    let routed = route(&spec, &inputs.qp, interrupt)?;
    let active: Vec<usize> = (0..n).collect();
    let weights = expand_and_floor(n, &active, &routed.solution.x, inputs.min_weight);
    Ok(package(weights, &routed))
}

/// Solve the spec through the routed backend policy. An explicitly requested
/// backend that is not compiled in is a deliberate user choice, so it propagates
/// as an error naming the missing feature; any other setup failure falls back to a
/// degenerate solution so the R layer raises the convergence condition rather than
/// propagating a bug as a panic.
fn route(
    spec: &QpSpec,
    opts: &QpOptions,
    interrupt: &dyn Fn() -> bool,
) -> Result<RoutedSolution, SbwError> {
    match solve_psd(spec, opts, interrupt) {
        Ok(routed) => Ok(routed),
        Err(QpError::BackendUnavailable(name)) => Err(SbwError::BackendUnavailable(name)),
        Err(_) => Ok(RoutedSolution {
            solution: degenerate(spec),
            backend: "osqp",
            fell_back: false,
        }),
    }
}

/// Pack floored weights and a routed solution into the result record.
fn package(weights: Vec<f64>, routed: &RoutedSolution) -> SbwResult {
    let solution = &routed.solution;
    SbwResult {
        weights,
        duals: solution.duals.clone(),
        converged: solution.status.is_solved(),
        interrupted: solution.interrupted,
        iterations: solution.iterations,
        objective: solution.obj,
        backend: routed.backend,
        fell_back: routed.fell_back,
        status: solution.status.as_str(),
        pri_res: solution.pri_res,
        dua_res: solution.dua_res,
    }
}

/// A degenerate solution for a spec a backend refused to set up: zero weights
/// carrying the failure status so the R layer raises the convergence condition.
fn degenerate(spec: &QpSpec) -> QpSolution {
    QpSolution {
        x: vec![0.0; spec.n],
        duals: vec![0.0; spec.m],
        status: QpStatus::NonConvex,
        iterations: 0,
        obj: 0.0,
        pri_res: 0.0,
        dua_res: 0.0,
        interrupted: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn discrete_inputs<'a>(
        levels: &'a [i32],
        n_levels: usize,
        s: &'a [f64],
        estimand: SbwEstimand,
        moment_covs: &'a [f64],
        targets: &'a [f64],
        tols: &'a [f64],
    ) -> SbwDiscreteInputs<'a> {
        SbwDiscreteInputs {
            n: levels.len(),
            levels,
            n_levels,
            estimand,
            s,
            min_weight: 1e-8,
            norm: SbwNorm::L2,
            moment_covs,
            n_moments: targets.len(),
            targets,
            tols,
            qp: QpOptions::default(),
        }
    }

    /// The group-normalized weighted mean of a moment column within one level.
    fn group_mean(levels: &[i32], w: &[f64], z: &[f64], level: i32) -> f64 {
        let mut num = 0.0;
        let mut den = 0.0;
        for i in 0..levels.len() {
            if levels[i] == level {
                num += w[i] * z[i];
                den += w[i];
            }
        }
        num / den
    }

    #[test]
    fn a_slack_tolerance_keeps_uniform_weights() {
        // A moment column already balanced across the two groups at uniform
        // weights, so the minimum-dispersion solution is the uniform weighting.
        let levels = [0, 0, 1, 1];
        let s = vec![1.0; 4];
        // Each group has mean zero, so uniform weights meet any target-zero band.
        let z = vec![1.0, -1.0, 1.0, -1.0];
        let inputs = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.1]);
        let result = solve_discrete(&inputs, &|| false).unwrap();
        assert!(result.converged, "status {}", result.status);
        for w in &result.weights {
            assert!((w - 1.0).abs() < 1e-4, "weight {w} is far from one");
        }
    }

    #[test]
    fn a_tight_moment_binds_and_is_hand_checkable() {
        // Group zero carries an imbalanced column z = (2, -1); the exact target of
        // zero forces w0 + w1 = 2 and 2 w0 - w1 = 0, whose minimum-norm solution
        // is w0 = 2/3, w1 = 4/3. Group one is balanced at uniform weights.
        let levels = [0, 0, 1, 1];
        let s = vec![1.0; 4];
        let z = vec![2.0, -1.0, 1.0, -1.0];
        let inputs = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.0]);
        let result = solve_discrete(&inputs, &|| false).unwrap();
        assert!(result.converged, "status {}", result.status);
        let w = &result.weights;
        assert!((w[0] - 2.0 / 3.0).abs() < 1e-4, "w0 = {}", w[0]);
        assert!((w[1] - 4.0 / 3.0).abs() < 1e-4, "w1 = {}", w[1]);
        // The weighted mean of the constrained group returns to the target.
        assert!(group_mean(&levels, w, &z, 0).abs() < 1e-6);
        // Both group totals equal their sizes.
        assert!((w[0] + w[1] - 2.0).abs() < 1e-4);
        assert!((w[2] + w[3] - 2.0).abs() < 1e-4);
    }

    #[test]
    fn an_att_holds_the_focal_group_at_unit_weight() {
        // Focal on level one: level-one units keep unit weight, level-zero units
        // are pulled to the level-one column mean, and both group totals equal the
        // focal count.
        let levels = [1, 1, 1, 0, 0, 0];
        let s = vec![1.0; 6];
        // Level one mean is 1; level zero straddles it so the pull is feasible.
        let z = vec![1.0, 1.0, 1.0, 3.0, 1.0, -1.0];
        let target = 1.0;
        let targets = [target];
        let inputs = discrete_inputs(
            &levels,
            2,
            &s,
            SbwEstimand::Focal { focal: 1 },
            &z,
            &targets,
            &[0.0],
        );
        let result = solve_discrete(&inputs, &|| false).unwrap();
        assert!(result.converged, "status {}", result.status);
        let w = &result.weights;
        for wi in w.iter().take(3) {
            assert!((wi - 1.0).abs() < 1e-9, "focal weight {wi}");
        }
        let n_focal = 3.0;
        let sum_focal: f64 = (0..6).filter(|&i| levels[i] == 1).map(|i| w[i]).sum();
        let sum_other: f64 = (0..6).filter(|&i| levels[i] == 0).map(|i| w[i]).sum();
        assert!((sum_focal - n_focal).abs() < 1e-9, "focal sum {sum_focal}");
        assert!((sum_other - n_focal).abs() < 1e-4, "other sum {sum_other}");
        // The reweighted level-zero group reaches the focal mean.
        assert!((group_mean(&levels, w, &z, 0) - target).abs() < 1e-6);
    }

    #[test]
    fn the_minimum_weight_floor_binds_the_box() {
        // A minimum weight equal to the group mean forces every weight in a group
        // that sums to its size onto the floor: three units summing to three with
        // each at least one are all exactly one.
        let levels = [0, 0, 0, 1, 1, 1];
        let s = vec![1.0; 6];
        let z = vec![0.5, 0.0, -0.5, 0.5, 0.0, -0.5];
        let mut inputs = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.5]);
        inputs.min_weight = 1.0;
        let result = solve_discrete(&inputs, &|| false).unwrap();
        assert!(result.converged, "status {}", result.status);
        for w in &result.weights {
            assert!(*w >= 1.0 - 1e-9, "weight {w} below the floor");
            assert!((w - 1.0).abs() < 1e-6, "weight {w} not at the floor");
        }
    }

    #[test]
    fn tightening_the_tolerance_cannot_lower_the_dispersion() {
        // A smaller tolerance is a smaller feasible set, so the minimized weight
        // dispersion cannot fall. The objective the solver reports is the sum of
        // squared weights, monotone in the tolerance.
        let levels = [0, 0, 0, 1, 1, 1];
        let s = vec![1.0; 6];
        let z = vec![1.5, -0.3, -1.2, 0.4, -0.1, -0.3];
        let solve_at = |tol: f64| {
            let tols = [tol];
            let inputs = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &tols);
            solve_discrete(&inputs, &|| false).unwrap()
        };
        let tight = solve_at(0.05);
        let loose = solve_at(0.4);
        assert!(tight.converged && loose.converged);
        assert!(
            tight.objective >= loose.objective - 1e-6,
            "tight {} !>= loose {}",
            tight.objective,
            loose.objective
        );
    }

    #[test]
    fn an_infeasible_band_reports_primal_infeasibility() {
        // A column that perfectly separates the groups cannot be balanced: every
        // level-zero unit shares one value and every level-one unit another, so no
        // reweighting moves the group means and a tight band is infeasible.
        let levels = [0, 0, 0, 1, 1, 1];
        let s = vec![1.0; 6];
        let z = vec![-1.0, -1.0, -1.0, 1.0, 1.0, 1.0];
        let inputs = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.01]);
        let result = solve_discrete(&inputs, &|| false).unwrap();
        assert!(!result.converged);
        assert_eq!(result.status, "primal_infeasible");
    }

    #[test]
    fn an_l1_tight_moment_recovers_the_determinate_weights() {
        // The determinate tight-moment instance has a unique feasible weighting, so
        // every dispersion norm returns it. Under L1 the objective is the summed
        // absolute departure from one: |2/3 - 1| + |4/3 - 1| in the constrained
        // group and zero in the balanced group, which is 2/3.
        let levels = [0, 0, 1, 1];
        let s = vec![1.0; 4];
        let z = vec![2.0, -1.0, 1.0, -1.0];
        let mut inputs = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.0]);
        inputs.norm = SbwNorm::L1;
        let result = solve_discrete(&inputs, &|| false).unwrap();
        assert!(result.converged, "status {}", result.status);
        let w = &result.weights;
        assert!((w[0] - 2.0 / 3.0).abs() < 1e-4, "w0 = {}", w[0]);
        assert!((w[1] - 4.0 / 3.0).abs() < 1e-4, "w1 = {}", w[1]);
        assert!(
            (result.objective - 2.0 / 3.0).abs() < 1e-3,
            "obj {}",
            result.objective
        );
    }

    #[test]
    fn an_l1_slack_tolerance_keeps_uniform_weights() {
        // A column already balanced at uniform weights leaves the summed absolute
        // deviation at its floor of zero, so the L1 solution is the uniform
        // weighting just as the L2 solution is.
        let levels = [0, 0, 1, 1];
        let s = vec![1.0; 4];
        let z = vec![1.0, -1.0, 1.0, -1.0];
        let mut inputs = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.1]);
        inputs.norm = SbwNorm::L1;
        let result = solve_discrete(&inputs, &|| false).unwrap();
        assert!(result.converged, "status {}", result.status);
        for w in &result.weights {
            assert!((w - 1.0).abs() < 1e-3, "weight {w} is far from one");
        }
        // The uniform weighting incurs no absolute deviation, so the objective is
        // the stabilizing ridge alone.
        assert!(result.objective < 1e-3, "obj {}", result.objective);
    }

    #[test]
    fn the_l1_and_linf_objectives_are_hand_checkable_and_differ() {
        // Each group carries the column (1, 0, 0) held exactly at a target of 0.5.
        // The first unit's group-normalized mean is w0 / 3, so the exact target
        // pins w0 = 1.5 in each group and leaves the remaining pair summing to 1.5.
        // Under L1 the per-group objective is |1.5 - 1| plus the pair's summed
        // shortfall of (1 - w1) + (1 - w2) = 0.5, so 1.0 per group and 2.0 across
        // both. Under Linf the single shared deviation is the largest departure,
        // 0.5, dominated by the pinned first units.
        let levels = [0, 0, 0, 1, 1, 1];
        let s = vec![1.0; 6];
        let z = vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let targets = [0.5];
        let tols = [0.0];

        let mut l1 = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &targets, &tols);
        l1.norm = SbwNorm::L1;
        let l1_result = solve_discrete(&l1, &|| false).unwrap();
        assert!(l1_result.converged, "l1 status {}", l1_result.status);
        assert!(
            (l1_result.weights[0] - 1.5).abs() < 1e-4,
            "w0 = {}",
            l1_result.weights[0]
        );
        assert!(
            (l1_result.weights[3] - 1.5).abs() < 1e-4,
            "w3 = {}",
            l1_result.weights[3]
        );
        assert!(
            (l1_result.objective - 2.0).abs() < 1e-3,
            "l1 obj {}",
            l1_result.objective
        );

        let mut linf = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &targets, &tols);
        linf.norm = SbwNorm::Linf;
        let linf_result = solve_discrete(&linf, &|| false).unwrap();
        assert!(linf_result.converged, "linf status {}", linf_result.status);
        assert!(
            (linf_result.weights[0] - 1.5).abs() < 1e-4,
            "w0 = {}",
            linf_result.weights[0]
        );
        // The reported objective adds the stabilizing ridge to the largest
        // deviation, so it sits just above 0.5; the deviation itself is read from
        // the weights.
        let linf_max = linf_result
            .weights
            .iter()
            .map(|w| (w - 1.0).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            (linf_max - 0.5).abs() < 1e-3,
            "linf max deviation {linf_max}"
        );
        // The largest departure under Linf is no larger than the one the L1
        // solution incurs, the defining property of the supremum norm.
        let l1_max = l1_result
            .weights
            .iter()
            .map(|w| (w - 1.0).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            linf_max <= l1_max + 1e-6,
            "linf max {linf_max} > l1 max {l1_max}"
        );
    }

    #[test]
    fn the_min_weight_floor_binds_under_l1() {
        // A minimum weight equal to the group mean forces every weight in a group
        // that sums to its size onto the floor, whatever the norm.
        let levels = [0, 0, 0, 1, 1, 1];
        let s = vec![1.0; 6];
        let z = vec![0.5, 0.0, -0.5, 0.5, 0.0, -0.5];
        let mut inputs = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.5]);
        inputs.norm = SbwNorm::L1;
        inputs.min_weight = 1.0;
        let result = solve_discrete(&inputs, &|| false).unwrap();
        assert!(result.converged, "status {}", result.status);
        for w in &result.weights {
            assert!(*w >= 1.0 - 1e-9, "weight {w} below the floor");
            assert!((w - 1.0).abs() < 1e-6, "weight {w} not at the floor");
        }
    }

    #[test]
    fn an_infeasible_band_reports_primal_infeasibility_under_linf() {
        // A perfectly separating column cannot be balanced under any norm, so a
        // tight band certifies primal infeasibility just as it does for L2.
        let levels = [0, 0, 0, 1, 1, 1];
        let s = vec![1.0; 6];
        let z = vec![-1.0, -1.0, -1.0, 1.0, 1.0, 1.0];
        let mut inputs = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.01]);
        inputs.norm = SbwNorm::Linf;
        let result = solve_discrete(&inputs, &|| false).unwrap();
        assert!(!result.converged);
        assert_eq!(result.status, "primal_infeasible");
    }

    #[test]
    fn the_l1_solve_is_deterministic() {
        let levels = [0, 0, 0, 1, 1, 1];
        let s = vec![1.0; 6];
        let z = vec![1.5, -0.3, -1.2, 0.4, -0.1, -0.3];
        let mut inputs = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.05]);
        inputs.norm = SbwNorm::L1;
        let a = solve_discrete(&inputs, &|| false).unwrap();
        let b = solve_discrete(&inputs, &|| false).unwrap();
        for (x, y) in a.weights.iter().zip(&b.weights) {
            assert_eq!(x.to_bits(), y.to_bits());
        }
    }

    #[test]
    fn a_continuous_l1_solve_meets_the_correlation_tolerance() {
        // The continuous L1 solve bounds the same weighted correlation the L2 solve
        // does; only the dispersion objective changes. The weights stay
        // non-negative and the bounded row sits inside the tolerance.
        let n = 40;
        let treat: Vec<f64> = (0..n)
            .map(|i| ((i as f64) * 0.13).sin() * 2.0 + (i as f64) * 0.02)
            .collect();
        let raw: Vec<f64> = (0..n)
            .map(|i| treat[i] * 0.7 + ((i as f64) * 0.37).cos())
            .collect();
        let mean: f64 = raw.iter().sum::<f64>() / n as f64;
        let sd = (raw.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0)).sqrt();
        let covs: Vec<f64> = raw.iter().map(|v| (v - mean) / sd).collect();
        let s = vec![1.0; n];
        let tol = 0.1;
        let inputs = SbwContInputs {
            n,
            treat: &treat,
            covs: &covs,
            n_covs: 1,
            s: &s,
            min_weight: 1e-8,
            norm: SbwNorm::L1,
            tols: &[tol],
            qp: QpOptions::default(),
        };
        let result = solve_cont(&inputs, &|| false).unwrap();
        assert!(result.converged, "status {}", result.status);
        assert!(result.weights.iter().all(|&w| w >= 0.0));

        let a_mean: f64 = treat.iter().sum::<f64>() / n as f64;
        let a_sd = weighted_variance(&treat, &s).sqrt();
        let constrained = result
            .weights
            .iter()
            .enumerate()
            .map(|(i, &w)| w * ((treat[i] - a_mean) / a_sd) * covs[i])
            .sum::<f64>()
            / n as f64;
        assert!(
            constrained.abs() <= tol + 1e-6,
            "constrained correlation {constrained}"
        );
    }

    #[test]
    fn a_pending_interrupt_surfaces_on_the_result() {
        let levels = [0, 0, 0, 1, 1, 1];
        let s = vec![1.0; 6];
        let z = vec![1.5, -0.3, -1.2, 0.4, -0.1, -0.3];
        let mut inputs = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.02]);
        inputs.qp = QpOptions {
            chunk_iters: 1,
            ..QpOptions::default()
        };
        let result = solve_discrete(&inputs, &|| true).unwrap();
        assert!(result.interrupted);
        assert_eq!(result.status, "interrupted");
        assert!(!result.converged);
    }

    #[test]
    #[cfg(feature = "qp-clarabel")]
    fn an_empty_active_set_fails_the_same_way_under_either_backend() {
        // Every unit absent from every level leaves no weight variable to solve for
        // while the group and moment rows survive, so there is no problem to hand a
        // backend. Both refuse it, and the refusal becomes the same failure status
        // whichever backend was asked, so the R layer raises one condition rather
        // than two.
        use crate::qp::QpBackendChoice;
        let levels = [-1, -1, -1, -1];
        let s = vec![1.0; 4];
        let z = vec![0.5, -0.5, 0.5, -0.5];
        let mut osqp = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.1]);
        osqp.qp.backend = QpBackendChoice::Osqp;
        let osqp_result = solve_discrete(&osqp, &|| false).unwrap();
        let mut clarabel = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.1]);
        clarabel.qp.backend = QpBackendChoice::Clarabel;
        let clarabel_result = solve_discrete(&clarabel, &|| false).unwrap();
        assert!(!osqp_result.converged, "osqp status {}", osqp_result.status);
        assert_eq!(
            clarabel_result.status, osqp_result.status,
            "clarabel reported {} against osqp {}",
            clarabel_result.status, osqp_result.status
        );
    }

    #[test]
    #[cfg(feature = "qp-clarabel")]
    fn a_pending_interrupt_surfaces_from_the_clarabel_backend() {
        // The interior-point backend has no iteration chunking to stop at, so it
        // stops through its own termination callback. The result reads the same as
        // the ADMM backend's.
        use crate::qp::QpBackendChoice;
        let levels = [0, 0, 0, 1, 1, 1];
        let s = vec![1.0; 6];
        let z = vec![1.5, -0.3, -1.2, 0.4, -0.1, -0.3];
        let mut inputs = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.02]);
        inputs.qp.backend = QpBackendChoice::Clarabel;
        let result = solve_discrete(&inputs, &|| true).unwrap();
        assert!(result.interrupted);
        assert_eq!(result.status, "interrupted");
        assert!(!result.converged);
        assert_eq!(result.backend, "clarabel");
    }

    #[test]
    fn a_categorical_ate_normalizes_each_group_to_its_size() {
        // Three levels, each pulled to a shared target of zero; every group's total
        // returns to its own size.
        let levels = [0, 0, 0, 1, 1, 1, 2, 2, 2];
        let s = vec![1.0; 9];
        let z = vec![
            0.6, -0.2, -0.4, 0.5, -0.1, -0.4, 0.3, 0.0, -0.3, // one moment column
        ];
        let inputs = discrete_inputs(&levels, 3, &s, SbwEstimand::Ate, &z, &[0.0], &[0.05]);
        let result = solve_discrete(&inputs, &|| false).unwrap();
        assert!(result.converged, "status {}", result.status);
        for level in 0..3 {
            let sum: f64 = (0..9)
                .filter(|&i| levels[i] == level)
                .map(|i| result.weights[i])
                .sum();
            assert!((sum - 3.0).abs() < 1e-4, "level {level} sum {sum}");
        }
    }

    #[test]
    fn a_continuous_solve_meets_the_correlation_tolerance() {
        // An exposure correlated with the covariate: uniform weights leave a large
        // weighted correlation, and the bounded-correlation row pulls it within the
        // tolerance while the dispersion objective keeps the weights close to one.
        let n = 40;
        let treat: Vec<f64> = (0..n)
            .map(|i| ((i as f64) * 0.13).sin() * 2.0 + (i as f64) * 0.02)
            .collect();
        // A covariate standardized to roughly unit variance, correlated with the
        // exposure but not collinear, so the correlation constraint is well posed.
        let raw: Vec<f64> = (0..n)
            .map(|i| treat[i] * 0.7 + ((i as f64) * 0.37).cos())
            .collect();
        let mean: f64 = raw.iter().sum::<f64>() / n as f64;
        let sd = (raw.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0)).sqrt();
        let covs: Vec<f64> = raw.iter().map(|v| (v - mean) / sd).collect();
        let s = vec![1.0; n];
        let tol = 0.1;
        let inputs = SbwContInputs {
            n,
            treat: &treat,
            covs: &covs,
            n_covs: 1,
            s: &s,
            min_weight: 1e-8,
            norm: SbwNorm::L2,
            tols: &[tol],
            qp: QpOptions::default(),
        };
        let result = solve_cont(&inputs, &|| false).unwrap();
        assert!(result.converged, "status {}", result.status);
        assert!(result.weights.iter().all(|&w| w >= 0.0));

        // The row the solver enforces bounds the weighted product of the covariate
        // and the exposure standardized to the sampling-weight scale. Checking that
        // same quantity confirms the constraint the quadratic program imposes; the
        // exposure standardization is fixed from the sampling weights, matching the
        // constraint assembly rather than the reweighted sample.
        let a_mean: f64 = treat.iter().sum::<f64>() / n as f64;
        let a_sd = weighted_variance(&treat, &s).sqrt();
        let constrained = result
            .weights
            .iter()
            .enumerate()
            .map(|(i, &w)| w * ((treat[i] - a_mean) / a_sd) * covs[i])
            .sum::<f64>()
            / n as f64;
        assert!(
            constrained.abs() <= tol + 1e-6,
            "constrained correlation {constrained}"
        );
        // The pull reduces the association below the large uniform-weight value.
        let uniform: f64 = (0..n)
            .map(|i| ((treat[i] - a_mean) / a_sd) * covs[i])
            .sum::<f64>()
            / n as f64;
        assert!(
            uniform.abs() > tol,
            "uniform correlation {uniform} should exceed tol"
        );
    }

    #[test]
    fn the_solve_is_deterministic() {
        let levels = [0, 0, 0, 1, 1, 1];
        let s = vec![1.0; 6];
        let z = vec![1.5, -0.3, -1.2, 0.4, -0.1, -0.3];
        let inputs = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.05]);
        let a = solve_discrete(&inputs, &|| false).unwrap();
        let b = solve_discrete(&inputs, &|| false).unwrap();
        for (x, y) in a.weights.iter().zip(&b.weights) {
            assert_eq!(x.to_bits(), y.to_bits());
        }
    }

    /// SplitMix64, replicating the benchmark generator so the ill-scaled instance
    /// is bit-identical to the one the backend comparison certified.
    struct SplitMix64(u64);

    impl SplitMix64 {
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }

        fn unit(&mut self) -> f64 {
            let bits = self.next_u64() >> 11;
            (bits as f64 + 0.5) * (1.0 / (1u64 << 53) as f64)
        }

        fn normal(&mut self) -> f64 {
            let u1 = self.unit();
            let u2 = self.unit();
            (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
        }
    }

    /// The ill-scaled binary average-treatment-effect instance the backend
    /// benchmark found osqp falsely certifies primal infeasible: four covariate
    /// columns spanning 1e-3 to 1e3, left unstandardized so the constraint matrix
    /// is badly conditioned, with per-column bands proportional to each column's
    /// own spread. This is the small-n form of the same misfire the n=20000
    /// average-treatment-effect instance shows. Returns levels, covariates,
    /// targets, and tolerances.
    fn ill_scaled_instance(n: usize) -> (Vec<i32>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let p = 4usize;
        let mut rng = SplitMix64(0x5B_15CA);
        let mut covs = vec![0.0; n * p];
        for j in 0..p {
            let scale = 10f64.powi(2 * j as i32 - 3);
            for i in 0..n {
                covs[j * n + i] = rng.normal() * scale;
            }
        }
        let mut levels = vec![0i32; n];
        for (i, level) in levels.iter_mut().enumerate() {
            let lp = 0.5 * (covs[i] / 1e-3) + 0.4 * (covs[n + i] / 0.1);
            let pr = 1.0 / (1.0 + (-lp).exp());
            *level = if rng.unit() < pr { 1 } else { 0 };
        }
        let targets: Vec<f64> = (0..p)
            .map(|j| covs[j * n..(j + 1) * n].iter().sum::<f64>() / n as f64)
            .collect();
        let tols: Vec<f64> = (0..p)
            .map(|j| {
                let col = &covs[j * n..(j + 1) * n];
                let mean = col.iter().sum::<f64>() / n as f64;
                let sd =
                    (col.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0)).sqrt();
                0.05 * sd
            })
            .collect();
        (levels, covs, targets, tols)
    }

    #[test]
    #[cfg(feature = "qp-clarabel")]
    fn an_ill_scaled_instance_falls_back_to_clarabel() {
        use crate::qp::QpBackendChoice;
        let n = 2000usize;
        let (levels, covs, targets, tols) = ill_scaled_instance(n);
        let s = vec![1.0; n];

        // osqp alone falsely certifies this feasible instance primal infeasible.
        let mut osqp = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &covs, &targets, &tols);
        osqp.qp.backend = QpBackendChoice::Osqp;
        let osqp_result = solve_discrete(&osqp, &|| false).unwrap();
        assert_eq!(
            osqp_result.status, "primal_infeasible",
            "osqp did not reproduce the misfire; status {}",
            osqp_result.status
        );
        assert_eq!(osqp_result.backend, "osqp");
        assert!(!osqp_result.fell_back);

        // The default auto routing detects the certificate and re-solves with
        // clarabel, which finds the feasible minimum.
        let auto = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &covs, &targets, &tols);
        assert_eq!(auto.qp.backend, QpBackendChoice::Auto);
        let auto_result = solve_discrete(&auto, &|| false).unwrap();
        assert!(
            auto_result.converged,
            "fallback did not solve; status {}",
            auto_result.status
        );
        assert_eq!(auto_result.backend, "clarabel");
        assert!(auto_result.fell_back);
        assert!(
            auto_result
                .weights
                .iter()
                .all(|&w| w.is_finite() && w >= 0.0),
            "fallback weights are not finite and non-negative"
        );
    }

    #[test]
    #[cfg(not(feature = "qp-clarabel"))]
    fn an_explicit_clarabel_choice_without_the_feature_errors() {
        // Without the clarabel feature an explicit clarabel request is refused with
        // a named error rather than silently substituting osqp, so a caller that
        // asked for the interior-point backend by name learns it is unavailable.
        // The shipped R package always compiles the feature in, so this path is not
        // reachable from R; it is covered here at the Rust boundary instead.
        use crate::qp::QpBackendChoice;
        let levels = [0, 0, 1, 1];
        let s = vec![1.0; 4];
        let z = vec![1.0, -1.0, 1.0, -1.0];
        let mut inputs = discrete_inputs(&levels, 2, &s, SbwEstimand::Ate, &z, &[0.0], &[0.1]);
        inputs.qp.backend = QpBackendChoice::Clarabel;
        let err = solve_discrete(&inputs, &|| false).unwrap_err();
        assert_eq!(err, SbwError::BackendUnavailable("clarabel"));
    }
}
