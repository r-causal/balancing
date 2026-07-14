//! Energy balancing: weights that minimize an energy statistic of covariate
//! balance under a simplex and optional moment constraints.
//!
//! For a binary or multi-category exposure the objective is the energy distance
//! between each reweighted group and a target sample; the improved variant for
//! the average treatment effect adds the energy distance between each pair of
//! groups. Both reduce to a quadratic form in the weights built from the pairwise
//! covariate distance matrix: the quadratic term is the negated distance matrix
//! scaled by the group-normalization outer product, and the linear term is the
//! group's cross energy with the full sample.
//!
//! For a continuous exposure the objective is the weighted distance covariance
//! between the exposure and the covariates, plus the energy distances between the
//! weighted and original exposure and covariate distributions, following Huling,
//! Greifer, and Chen. The distance covariance enters through the double-centered
//! product of the exposure and covariate distance matrices.
//!
//! Every form is a quadratic program with a simplex-type constraint set. The
//! quadratic term is indefinite by construction, so the solve routes to the ADMM
//! backend; the interior-point backend rejects it.

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::dist::{Distance, distance_matrix, pairwise};
use crate::qp::osqp::Osqp;
use crate::qp::{Convexity, QpBackend, QpOptions, QpSpec, QpStatus, objective};
use crate::threads::get_pool;

use super::qp_balance::{
    ConstraintBuilder, ZERO_SW, add_diagonal_penalty, doubled_dense, expand_and_floor,
    group_normalized,
};

/// The estimand an energy balancing solve targets.
#[derive(Debug, Clone, Copy)]
pub enum EnergyEstimand {
    /// Average treatment effect: every group is reweighted toward the full
    /// sample, optionally adding the between-group term of the improved variant.
    Ate { improved: bool },
    /// A focal-group effect (treated or control): the non-focal groups are
    /// reweighted toward the focal group, whose units keep unit weight.
    Focal { focal: usize },
}

/// Outcome of an energy balancing solve.
pub struct EnergyResult {
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
    /// The backend that produced the result, carried to R as the solver identity.
    pub backend: &'static str,
    /// The backend's terminal status name.
    pub status: &'static str,
    /// Primal residual from the backend.
    pub pri_res: f64,
    /// Dual residual from the backend.
    pub dua_res: f64,
}

/// Inputs for a binary or multi-category energy balancing solve.
pub struct EnergyDiscreteInputs<'a> {
    /// Column-major `n` by `p` covariates for the distance matrix.
    pub covs: &'a [f64],
    /// Number of units.
    pub n: usize,
    /// Number of covariate columns.
    pub p: usize,
    /// Distance definition the energy objective is built on.
    pub distance: Distance,
    /// Zero-based exposure level of each unit.
    pub levels: &'a [i32],
    /// Number of exposure levels.
    pub n_levels: usize,
    /// The estimand.
    pub estimand: EnergyEstimand,
    /// Sampling weights.
    pub s: &'a [f64],
    /// Minimum allowable weight.
    pub min_weight: f64,
    /// Weight penalty coefficient.
    pub lambda: f64,
    /// Column-major `n` by `q` moment-constraint covariates, standardized on the
    /// R side; empty for no moment constraints.
    pub moment_covs: &'a [f64],
    /// Number of moment-constraint columns `q`.
    pub n_moments: usize,
    /// Target weighted mean for each moment column.
    pub targets: &'a [f64],
    /// Tolerance band half-width source for each moment column.
    pub tols: &'a [f64],
    /// Worker threads for the distance assembly.
    pub threads: usize,
    /// Quadratic-program tuning.
    pub qp: QpOptions,
}

/// The interaction factor of the group-normalization outer product for two units
/// in levels `li` and `lj` under the estimand.
///
/// For the average treatment effect the standard variant pairs same-level units
/// with weight one; the improved variant adds the between-group term, which
/// raises the same-level weight to the number of levels and sets the
/// cross-level weight to minus one. For a focal estimand only same-level,
/// non-focal units interact.
fn nn_factor(estimand: EnergyEstimand, n_levels: usize, li: i32, lj: i32) -> f64 {
    match estimand {
        EnergyEstimand::Ate { improved: false } => {
            if li == lj {
                1.0
            } else {
                0.0
            }
        }
        EnergyEstimand::Ate { improved: true } => {
            if li == lj {
                n_levels as f64
            } else {
                -1.0
            }
        }
        EnergyEstimand::Focal { .. } => {
            if li == lj {
                1.0
            } else {
                0.0
            }
        }
    }
}

/// Solve a binary or multi-category energy balancing problem.
pub fn solve_discrete(
    inputs: &EnergyDiscreteInputs<'_>,
    interrupt: &dyn Fn() -> bool,
) -> EnergyResult {
    let n = inputs.n;
    let dist = distance_matrix(
        inputs.covs,
        n,
        inputs.p,
        inputs.distance,
        inputs.s,
        inputs.threads,
    );
    let s_norm = group_normalized(inputs.s, inputs.levels, inputs.n_levels);

    // Level sizes and each unit's own-column normalization value.
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

    // The active variables and the level rows that carry group and moment
    // constraints. For the average treatment effect every unit is a variable and
    // every level is constrained; for a focal estimand only the non-focal units
    // are variables and only the non-focal levels are constrained.
    let (active, group_levels, focal) = match inputs.estimand {
        EnergyEstimand::Ate { .. } => {
            let active: Vec<usize> = (0..n).filter(|&i| inputs.levels[i] >= 0).collect();
            let groups: Vec<usize> = (0..inputs.n_levels).collect();
            (active, groups, None)
        }
        EnergyEstimand::Focal { focal } => {
            let active: Vec<usize> = (0..n)
                .filter(|&i| inputs.levels[i] >= 0 && inputs.levels[i] as usize != focal)
                .collect();
            let groups: Vec<usize> = (0..inputs.n_levels).filter(|&t| t != focal).collect();
            (active, groups, Some(focal))
        }
    };
    let nvar = active.len();

    // Quadratic term P = -d * nn with the group-normalization outer product,
    // built over the active variables.
    let mut pmat = vec![0.0; nvar * nvar];
    for a in 0..nvar {
        let ia = active[a];
        for b in 0..nvar {
            let ib = active[b];
            let factor = nn_factor(
                inputs.estimand,
                inputs.n_levels,
                inputs.levels[ia],
                inputs.levels[ib],
            );
            if factor == 0.0 {
                continue;
            }
            let d = dist[ib * n + ia];
            pmat[b * nvar + a] = -d * swnt[ia] * swnt[ib] * factor;
        }
    }
    // The penalty scale is each active unit's own group-normalization value.
    let scale: Vec<f64> = active.iter().map(|&i| swnt[i]).collect();
    add_diagonal_penalty(&mut pmat, nvar, inputs.lambda, &scale);

    // Linear term q_a = cross_a * scale_a, where cross_a is the group's cross
    // energy with its target sample evaluated at active column a. The source runs
    // over the whole sample for the average treatment effect and over the focal
    // group for a focal estimand.
    let (src, mult): (Vec<f64>, f64) = match inputs.estimand {
        EnergyEstimand::Ate { .. } => (s_norm.clone(), 2.0 / n as f64),
        EnergyEstimand::Focal { .. } => {
            let f = focal.expect("focal estimand carries a focal level");
            let nf = n_t[f].max(1) as f64;
            let src: Vec<f64> = (0..n)
                .map(|i| {
                    if inputs.levels[i] >= 0 && inputs.levels[i] as usize == f {
                        s_norm[i] / nf
                    } else {
                        0.0
                    }
                })
                .collect();
            (src, 2.0)
        }
    };
    let cross = cross_energy(&dist, n, &active, &src, mult, inputs.threads);
    let q: Vec<f64> = (0..nvar).map(|a| cross[a] * scale[a]).collect();

    // Constraints: the box, the group sums, and the moment rows.
    let mut builder = ConstraintBuilder::new(nvar);
    let pinned: Vec<bool> = active
        .iter()
        .map(|&i| inputs.s[i].abs() < ZERO_SW)
        .collect();
    builder.add_box(inputs.min_weight, &pinned);

    for &t in &group_levels {
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

    // The tolerance is split across the constrained groups for the average
    // treatment effect, which bounds each group's mean around the shared target
    // so the between-group difference stays within the tolerance; a focal
    // estimand constrains one side and uses the full tolerance.
    let tol_half = match inputs.estimand {
        EnergyEstimand::Ate { .. } => 0.5,
        EnergyEstimand::Focal { .. } => 1.0,
    };
    for &t in &group_levels {
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
            let band = inputs.tols[c] * tol_half;
            builder.add_dense_row(&coeffs, inputs.targets[c] - band, inputs.targets[c] + band);
        }
    }

    let (m, indptr, indices, values, l, u) = builder.finish();
    let spec = QpSpec {
        n: nvar,
        m,
        p: doubled_dense(pmat),
        q,
        a_indptr: indptr,
        a_indices: indices,
        a_values: values,
        l,
        u,
        convexity: Convexity::Indefinite,
    };

    let solution = Osqp
        .solve(&spec, &inputs.qp, interrupt)
        .unwrap_or_else(|e| degenerate(&spec, e));
    package_discrete(inputs, &active, &solution)
}

/// The source-weighted cross energy of each active column.
///
/// The value at active column `a` is `mult` times the source-weighted sum of
/// distances from every unit to the unit at that column. Work is split over the
/// active columns and is deterministic by construction because each column is an
/// independent sum in a fixed index order.
fn cross_energy(
    dist: &[f64],
    n: usize,
    active: &[usize],
    src: &[f64],
    mult: f64,
    threads: usize,
) -> Vec<f64> {
    let pool = get_pool(threads);
    pool.install(|| {
        (0..active.len())
            .into_par_iter()
            .map(|a| {
                let j = active[a];
                let mut acc = 0.0;
                for i in 0..n {
                    acc += src[i] * dist[j * n + i];
                }
                mult * acc
            })
            .collect()
    })
}

/// Expand the active solution to the full unit order and apply the reference
/// post-processing: floor at the minimum weight, and snap negligible weights to
/// zero when the floor is essentially zero.
fn package_discrete(
    inputs: &EnergyDiscreteInputs<'_>,
    active: &[usize],
    solution: &crate::qp::QpSolution,
) -> EnergyResult {
    let weights = expand_and_floor(inputs.n, active, &solution.x, inputs.min_weight);
    EnergyResult {
        weights,
        duals: solution.duals.clone(),
        converged: solution.status.is_solved(),
        interrupted: solution.interrupted,
        iterations: solution.iterations,
        objective: solution.obj,
        backend: "osqp",
        status: solution.status.as_str(),
        pri_res: solution.pri_res,
        dua_res: solution.dua_res,
    }
}

/// A degenerate solution for a spec the backend refused to set up: zero weights
/// carrying the failure status so the R layer raises the convergence condition.
fn degenerate(spec: &QpSpec, _err: crate::qp::QpError) -> crate::qp::QpSolution {
    crate::qp::QpSolution {
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

/// Inputs for a continuous-exposure energy balancing solve.
pub struct EnergyContInputs<'a> {
    /// Column-major `n` by `p` covariates for the covariate distance matrix.
    pub covs: &'a [f64],
    /// Number of units.
    pub n: usize,
    /// Number of covariate columns, the dimension used by the adjustment.
    pub p: usize,
    /// Continuous exposure of each unit.
    pub treat: &'a [f64],
    /// Distance definition for the covariate distance matrix.
    pub distance: Distance,
    /// Sampling weights.
    pub s: &'a [f64],
    /// Minimum allowable weight.
    pub min_weight: f64,
    /// Weight penalty coefficient.
    pub lambda: f64,
    /// Whether to weight the covariate energy distance by the dimension
    /// adjustment.
    pub dimension_adj: bool,
    /// Column-major `n` by `qd` distribution-moment covariate columns, centered
    /// and scaled on the R side, held exactly at the unweighted sample value.
    pub d_covs: &'a [f64],
    /// Number of distribution-moment covariate columns.
    pub n_d_covs: usize,
    /// Column-major `n` by `qt` distribution-moment exposure columns, centered
    /// and scaled on the R side.
    pub d_treat: &'a [f64],
    /// Number of distribution-moment exposure columns.
    pub n_d_treat: usize,
    /// Column-major `n` by `qb` correlation-constraint covariates, scaled on the
    /// R side.
    pub bal_covs: &'a [f64],
    /// Number of correlation-constraint columns.
    pub n_bal: usize,
    /// Tolerance for each correlation constraint.
    pub bal_tols: &'a [f64],
    /// Worker threads for the distance assembly.
    pub threads: usize,
    /// Quadratic-program tuning.
    pub qp: QpOptions,
}

/// Reliability-weighted variance of a single vector, the denominator matching the
/// per-column variance used by the distance transforms.
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

/// Double-center a symmetric distance matrix: `A_ij + grand - row_i - row_j`.
fn double_center(d: &[f64], n: usize) -> Vec<f64> {
    let mut row_means = vec![0.0; n];
    for j in 0..n {
        let mut acc = 0.0;
        for i in 0..n {
            acc += d[j * n + i];
        }
        row_means[j] = acc / n as f64;
    }
    let grand = row_means.iter().sum::<f64>() / n as f64;
    let mut out = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            out[j * n + i] = d[j * n + i] + grand - row_means[i] - row_means[j];
        }
    }
    out
}

/// Solve a continuous-exposure energy balancing problem.
pub fn solve_cont(inputs: &EnergyContInputs<'_>, interrupt: &dyn Fn() -> bool) -> EnergyResult {
    let n = inputs.n;
    let nf = n as f64;
    let xdist = distance_matrix(
        inputs.covs,
        n,
        inputs.p,
        inputs.distance,
        inputs.s,
        inputs.threads,
    );

    // Exposure distance on the exposure scaled by its weighted standard
    // deviation, so a one-dimensional Euclidean distance is the absolute
    // difference of the scaled exposure.
    let a_var = weighted_variance(inputs.treat, inputs.s);
    let a_sd = if a_var > 0.0 { a_var.sqrt() } else { 1.0 };
    let a_scaled: Vec<f64> = inputs.treat.iter().map(|&t| t / a_sd).collect();
    let adist = pairwise::euclidean(&a_scaled, n, 1, inputs.threads);

    // Sampling weights scaled to sum to n, the normalization the continuous
    // objective assumes.
    let s_sum: f64 = inputs.s.iter().sum();
    let s: Vec<f64> = if s_sum > 0.0 {
        inputs.s.iter().map(|&si| si * nf / s_sum).collect()
    } else {
        vec![1.0; n]
    };

    let xa = double_center(&xdist, n);
    let aa = double_center(&adist, n);

    let q_a = if inputs.dimension_adj {
        1.0 / (1.0 + (inputs.p as f64).sqrt())
    } else {
        0.5
    };
    let q_x = 1.0 - q_a;
    let n2 = nf * nf;

    // Quadratic term: the distance covariance kernel plus the marginal energy
    // terms, each scaled by the sampling-weight outer product.
    let mut pmat = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            let dcov = xa[j * n + i] * aa[j * n + i] / n2;
            let eb_a = -adist[j * n + i] / n2 * q_a;
            let eb_x = -xdist[j * n + i] / n2 * q_x;
            pmat[j * n + i] = (dcov + eb_a + eb_x) * s[i] * s[j];
        }
    }
    // Linear term: the cross marginal energy, scaled by the sampling weight.
    let mut q = vec![0.0; n];
    for j in 0..n {
        let mut ca = 0.0;
        let mut cx = 0.0;
        for i in 0..n {
            ca += s[i] * adist[j * n + i];
            cx += s[i] * xdist[j * n + i];
        }
        q[j] = (ca * 2.0 / n2 * q_a + cx * 2.0 / n2 * q_x) * s[j];
    }

    add_diagonal_penalty(&mut pmat, n, inputs.lambda, &s);

    // Constraints: the box, the single sum fixing the total to n, the exact
    // distribution-moment rows, and the correlation-constraint rows.
    let mut builder = ConstraintBuilder::new(n);
    let pinned: Vec<bool> = inputs.s.iter().map(|&si| si.abs() < ZERO_SW).collect();
    builder.add_box(inputs.min_weight, &pinned);

    let sum_coeffs = s.clone();
    builder.add_dense_row(&sum_coeffs, nf, nf);

    for c in 0..inputs.n_d_covs {
        let coeffs: Vec<f64> = (0..n)
            .map(|i| inputs.d_covs[c * n + i] * s[i] / nf)
            .collect();
        builder.add_dense_row(&coeffs, 0.0, 0.0);
    }
    for c in 0..inputs.n_d_treat {
        let coeffs: Vec<f64> = (0..n)
            .map(|i| inputs.d_treat[c * n + i] * s[i] / nf)
            .collect();
        builder.add_dense_row(&coeffs, 0.0, 0.0);
    }

    // Correlation constraints use the exposure centered and scaled to unit
    // weighted variance, so a bounded weighted mean of the product is a bounded
    // exposure-covariate correlation.
    let a_mean = {
        let sw: f64 = inputs.s.iter().sum();
        inputs
            .s
            .iter()
            .zip(inputs.treat)
            .map(|(&wi, &ti)| wi * ti)
            .sum::<f64>()
            / sw.max(f64::MIN_POSITIVE)
    };
    let treat_std: Vec<f64> = inputs.treat.iter().map(|&t| (t - a_mean) / a_sd).collect();
    let denom = (n as f64 - 1.0).max(1.0);
    for c in 0..inputs.n_bal {
        let coeffs: Vec<f64> = (0..n)
            .map(|i| inputs.bal_covs[c * n + i] * treat_std[i] * s[i] / denom)
            .collect();
        let tol = inputs.bal_tols[c].abs();
        builder.add_dense_row(&coeffs, -tol, tol);
    }

    let (m, indptr, indices, values, l, u) = builder.finish();
    let spec = QpSpec {
        n,
        m,
        p: doubled_dense(pmat),
        q,
        a_indptr: indptr,
        a_indices: indices,
        a_values: values,
        l,
        u,
        convexity: Convexity::Indefinite,
    };

    let solution = Osqp
        .solve(&spec, &inputs.qp, interrupt)
        .unwrap_or_else(|e| degenerate(&spec, e));

    // Every unit is a decision variable in the continuous solve, so the active
    // set is the full index range.
    let active: Vec<usize> = (0..n).collect();
    let weights = expand_and_floor(n, &active, &solution.x, inputs.min_weight);
    let obj = objective(&spec, &solution.x);
    EnergyResult {
        weights,
        duals: solution.duals.clone(),
        converged: solution.status.is_solved(),
        interrupted: solution.interrupted,
        iterations: solution.iterations,
        objective: obj,
        backend: "osqp",
        status: solution.status.as_str(),
        pri_res: solution.pri_res,
        dua_res: solution.dua_res,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The scaled-Euclidean distance matrix on standardized covariates, mirroring
    // the test oracle so the objective comparisons use the same metric.
    fn oracle_energy_ate(dist: &[f64], n: usize, treated: &[bool], w: &[f64]) -> f64 {
        let idx_t: Vec<usize> = (0..n).filter(|&i| treated[i]).collect();
        let idx_c: Vec<usize> = (0..n).filter(|&i| !treated[i]).collect();
        energy_to_uniform(dist, n, &idx_t, w)
            + energy_to_uniform(dist, n, &idx_c, w)
            + between(dist, n, &idx_t, &idx_c, w)
    }

    fn energy_to_uniform(dist: &[f64], n: usize, idx: &[usize], w: &[f64]) -> f64 {
        let wsum: f64 = idx.iter().map(|&i| w[i]).sum();
        let wn: Vec<f64> = idx.iter().map(|&i| w[i] / wsum).collect();
        let mut cross = 0.0;
        for (k, &i) in idx.iter().enumerate() {
            let mut row = 0.0;
            for j in 0..n {
                row += dist[j * n + i];
            }
            cross += wn[k] * row;
        }
        cross *= 2.0 / n as f64;
        let mut within = 0.0;
        for (a, &i) in idx.iter().enumerate() {
            for (b, &j) in idx.iter().enumerate() {
                within += wn[a] * dist[j * n + i] * wn[b];
            }
        }
        cross - within
    }

    fn between(dist: &[f64], n: usize, idx_a: &[usize], idx_b: &[usize], w: &[f64]) -> f64 {
        let sa: f64 = idx_a.iter().map(|&i| w[i]).sum();
        let sb: f64 = idx_b.iter().map(|&i| w[i]).sum();
        let wa: Vec<f64> = idx_a.iter().map(|&i| w[i] / sa).collect();
        let wb: Vec<f64> = idx_b.iter().map(|&i| w[i] / sb).collect();
        let mut cross = 0.0;
        for (p, &i) in idx_a.iter().enumerate() {
            for (qq, &j) in idx_b.iter().enumerate() {
                cross += wa[p] * dist[j * n + i] * wb[qq];
            }
        }
        cross *= 2.0;
        let mut wa2 = 0.0;
        for (p, &i) in idx_a.iter().enumerate() {
            for (r, &j) in idx_a.iter().enumerate() {
                wa2 += wa[p] * dist[j * n + i] * wa[r];
            }
        }
        let mut wb2 = 0.0;
        for (p, &i) in idx_b.iter().enumerate() {
            for (r, &j) in idx_b.iter().enumerate() {
                wb2 += wb[p] * dist[j * n + i] * wb[r];
            }
        }
        cross - wa2 - wb2
    }

    fn empty() -> (Vec<f64>, Vec<f64>) {
        (Vec::new(), Vec::new())
    }

    #[test]
    fn a_symmetric_problem_keeps_uniform_weights() {
        // Two identical groups: treated points equal control points, so the
        // groups already share a distribution and uniform weights are optimal.
        let covs = vec![0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0];
        let n = 4;
        let levels = [1, 1, 0, 0];
        let s = vec![1.0; n];
        let (targets, tols) = empty();
        let inputs = EnergyDiscreteInputs {
            covs: &covs,
            n,
            p: 2,
            distance: Distance::ScaledEuclidean,
            levels: &levels,
            n_levels: 2,
            estimand: EnergyEstimand::Ate { improved: true },
            s: &s,
            min_weight: 1e-8,
            lambda: 1e-4,
            moment_covs: &[],
            n_moments: 0,
            targets: &targets,
            tols: &tols,
            threads: 1,
            qp: QpOptions::default(),
        };
        let result = solve_discrete(&inputs, &|| false);
        assert!(result.converged, "status {}", result.status);
        for w in &result.weights {
            assert!((w - 1.0).abs() < 1e-3, "weight {w} is far from one");
        }
    }

    #[test]
    fn energy_balancing_reduces_the_ate_energy_distance() {
        // A confounded four-by-two design where the groups differ, so reweighting
        // strictly lowers the energy distance below the uniform value.
        let n = 8;
        // Column-major 8 by 2.
        let covs = vec![
            2.0, 2.2, 1.8, 2.1, -1.0, -1.2, -0.8, -0.9, // covariate 1
            1.0, 0.8, 1.2, 0.9, -2.0, -1.8, -2.2, -1.9, // covariate 2
        ];
        let levels = [1, 1, 1, 1, 0, 0, 0, 0];
        // Perturb so the groups are not perfectly symmetric.
        let covs = {
            let mut c = covs;
            c[0] += 0.5;
            c[4] -= 0.3;
            c
        };
        let s = vec![1.0; n];
        let (targets, tols) = empty();
        let inputs = EnergyDiscreteInputs {
            covs: &covs,
            n,
            p: 2,
            distance: Distance::ScaledEuclidean,
            levels: &levels,
            n_levels: 2,
            estimand: EnergyEstimand::Ate { improved: true },
            s: &s,
            min_weight: 1e-8,
            lambda: 1e-4,
            moment_covs: &[],
            n_moments: 0,
            targets: &targets,
            tols: &tols,
            threads: 1,
            qp: QpOptions::default(),
        };
        let result = solve_discrete(&inputs, &|| false);
        assert!(result.converged, "status {}", result.status);

        // The oracle is evaluated on the same scaled-Euclidean distance matrix
        // the solver used, so the inequality tests the weighting, not the metric.
        let dist = distance_matrix(&covs, n, 2, Distance::ScaledEuclidean, &s, 1);
        let treated: Vec<bool> = levels.iter().map(|&g| g == 1).collect();
        let weighted = oracle_energy_ate(&dist, n, &treated, &result.weights);
        let unit = vec![1.0; n];
        let unweighted = oracle_energy_ate(&dist, n, &treated, &unit);
        assert!(
            weighted < unweighted,
            "weighted {weighted} !< unweighted {unweighted}"
        );
        // Group sums equal the group sizes.
        let sum_t: f64 = (0..n)
            .filter(|&i| treated[i])
            .map(|i| result.weights[i])
            .sum();
        assert!((sum_t - 4.0).abs() < 1e-3, "treated sum {sum_t}");
    }

    #[test]
    fn a_moment_constraint_holds_the_group_means() {
        // With an exact moment constraint on a standardized covariate, each
        // group's weighted mean of that column equals the target. The covariates
        // driving the distance are confounded, so the weights are not uniform,
        // but the moment column overlaps between groups so balancing it is
        // feasible.
        let n = 8;
        let covs = vec![
            1.5, 1.2, 0.9, 1.1, -1.0, -1.3, -0.7, -1.1, // covariate 1
            0.5, 0.2, -0.1, 0.3, -0.6, -0.2, 0.1, -0.4, // covariate 2
        ];
        let levels = [1, 1, 1, 1, 0, 0, 0, 0];
        let s = vec![1.0; n];
        // A standardized moment column that straddles zero within each group, so
        // an exact constraint to the overall mean is feasible.
        let moment = {
            let raw = [0.9, -0.7, 0.5, -0.6, 0.8, -0.9, 0.4, -0.5];
            let mut col = raw.to_vec();
            let mean: f64 = col.iter().sum::<f64>() / n as f64;
            let sd =
                (col.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0)).sqrt();
            for v in col.iter_mut() {
                *v = (*v - mean) / sd;
            }
            col
        };
        let target: f64 = moment.iter().sum::<f64>() / n as f64;
        let inputs = EnergyDiscreteInputs {
            covs: &covs,
            n,
            p: 2,
            distance: Distance::ScaledEuclidean,
            levels: &levels,
            n_levels: 2,
            estimand: EnergyEstimand::Ate { improved: true },
            s: &s,
            min_weight: 1e-8,
            lambda: 1e-4,
            moment_covs: &moment,
            n_moments: 1,
            targets: &[target],
            tols: &[0.0],
            threads: 1,
            qp: QpOptions::default(),
        };
        let result = solve_discrete(&inputs, &|| false);
        assert!(result.converged, "status {}", result.status);
        // Weighted mean of the moment column in each group equals the target.
        for level in [0, 1] {
            let (mut num, mut den) = (0.0, 0.0);
            for i in 0..n {
                if levels[i] == level {
                    num += result.weights[i] * moment[i];
                    den += result.weights[i];
                }
            }
            assert!(
                (num / den - target).abs() < 1e-4,
                "group {level} mean {}",
                num / den
            );
        }
    }

    fn discrete_inputs<'a>(
        covs: &'a [f64],
        levels: &'a [i32],
        n_levels: usize,
        s: &'a [f64],
        estimand: EnergyEstimand,
    ) -> EnergyDiscreteInputs<'a> {
        EnergyDiscreteInputs {
            covs,
            n: levels.len(),
            p: covs.len() / levels.len(),
            distance: Distance::ScaledEuclidean,
            levels,
            n_levels,
            estimand,
            s,
            min_weight: 1e-8,
            lambda: 1e-4,
            moment_covs: &[],
            n_moments: 0,
            targets: &[],
            tols: &[],
            threads: 1,
            qp: QpOptions::default(),
        }
    }

    #[test]
    fn an_att_fit_pins_the_treated_and_targets_their_total() {
        // Focal on the treated: control units are reweighted toward the treated
        // group, treated units keep unit weight, and both group totals equal the
        // treated count.
        let n = 8;
        let covs = vec![
            1.5, 1.2, 0.9, 1.1, -1.0, -1.3, -0.7, -1.1, // covariate 1
            0.5, 0.2, -0.1, 0.3, -0.6, -0.2, 0.1, -0.4, // covariate 2
        ];
        let levels = [1, 1, 1, 1, 0, 0, 0, 0];
        let s = vec![1.0; n];
        let inputs = discrete_inputs(&covs, &levels, 2, &s, EnergyEstimand::Focal { focal: 1 });
        let result = solve_discrete(&inputs, &|| false);
        assert!(result.converged, "status {}", result.status);
        let n_treated = 4.0;
        let sum_t: f64 = (0..n)
            .filter(|&i| levels[i] == 1)
            .map(|i| result.weights[i])
            .sum();
        let sum_c: f64 = (0..n)
            .filter(|&i| levels[i] == 0)
            .map(|i| result.weights[i])
            .sum();
        assert!((sum_t - n_treated).abs() < 1e-6, "treated sum {sum_t}");
        assert!((sum_c - n_treated).abs() < 1e-3, "control sum {sum_c}");
        for (level, weight) in levels.iter().zip(&result.weights) {
            if *level == 1 {
                assert!((weight - 1.0).abs() < 1e-12, "treated weight {weight}");
            }
            assert!(*weight >= 1e-8 - 1e-12);
        }
    }

    #[test]
    fn a_categorical_ate_normalizes_each_group_to_its_size() {
        // Three levels: every group is reweighted to the full sample and its
        // total returns to its own size.
        let n = 9;
        let covs = vec![
            1.0, 1.2, 0.8, 0.1, -0.1, 0.2, -1.0, -1.2, -0.9, // covariate 1
            0.3, 0.1, 0.5, 1.1, 0.9, 1.3, -0.4, -0.2, -0.6, // covariate 2
        ];
        let levels = [0, 0, 0, 1, 1, 1, 2, 2, 2];
        let s = vec![1.0; n];
        let inputs = discrete_inputs(
            &covs,
            &levels,
            3,
            &s,
            EnergyEstimand::Ate { improved: true },
        );
        let result = solve_discrete(&inputs, &|| false);
        assert!(result.converged, "status {}", result.status);
        for level in 0..3 {
            let sum: f64 = (0..n)
                .filter(|&i| levels[i] == level)
                .map(|i| result.weights[i])
                .sum();
            assert!((sum - 3.0).abs() < 1e-3, "level {level} sum {sum}");
        }
        assert!(result.weights.iter().all(|&w| w >= 1e-8 - 1e-12));
    }

    #[test]
    fn the_improved_variant_differs_from_the_plain_variant() {
        let n = 8;
        let covs = vec![
            1.5, 1.2, 0.9, 1.1, -1.0, -1.3, -0.7, -1.1, // covariate 1
            0.5, 0.2, -0.1, 0.3, -0.6, -0.2, 0.1, -0.4, // covariate 2
        ];
        let levels = [1, 1, 1, 1, 0, 0, 0, 0];
        let s = vec![1.0; n];
        let improved = solve_discrete(
            &discrete_inputs(
                &covs,
                &levels,
                2,
                &s,
                EnergyEstimand::Ate { improved: true },
            ),
            &|| false,
        );
        let plain = solve_discrete(
            &discrete_inputs(
                &covs,
                &levels,
                2,
                &s,
                EnergyEstimand::Ate { improved: false },
            ),
            &|| false,
        );
        let differ = improved
            .weights
            .iter()
            .zip(&plain.weights)
            .any(|(a, b)| (a - b).abs() > 1e-4);
        assert!(differ, "improved and plain weights should differ");
    }

    #[test]
    fn a_pending_interrupt_surfaces_on_the_result() {
        let n = 8;
        let covs = vec![
            1.5, 1.2, 0.9, 1.1, -1.0, -1.3, -0.7, -1.1, // covariate 1
            0.5, 0.2, -0.1, 0.3, -0.6, -0.2, 0.1, -0.4, // covariate 2
        ];
        let levels = [1, 1, 1, 1, 0, 0, 0, 0];
        let s = vec![1.0; n];
        let mut inputs = discrete_inputs(
            &covs,
            &levels,
            2,
            &s,
            EnergyEstimand::Ate { improved: true },
        );
        inputs.qp = QpOptions {
            chunk_iters: 1,
            ..QpOptions::default()
        };
        let result = solve_discrete(&inputs, &|| true);
        assert!(result.interrupted);
        assert_eq!(result.status, "interrupted");
        assert!(!result.converged);
    }

    #[test]
    fn a_low_iteration_cap_exercises_the_non_solved_status_path() {
        // A tight tolerance under a small iteration cap makes the backend stop at a
        // chunk boundary before a full solve, exercising the warm-started
        // SolvedInaccurate / MaxIter path on the energy slice. The solve must
        // return a non-solved terminal status with finite weights rather than
        // corrupting the iterate across the warm-start chunks.
        let n = 8;
        let covs = vec![
            1.5, 1.2, 0.9, 1.1, -1.0, -1.3, -0.7, -1.1, // covariate 1
            0.5, 0.2, -0.1, 0.3, -0.6, -0.2, 0.1, -0.4, // covariate 2
        ];
        let levels = [1, 1, 1, 1, 0, 0, 0, 0];
        let s = vec![1.0; n];
        let mut inputs = discrete_inputs(
            &covs,
            &levels,
            2,
            &s,
            EnergyEstimand::Ate { improved: true },
        );
        inputs.qp = QpOptions {
            eps_abs: 1e-12,
            eps_rel: 1e-12,
            max_iter: 2,
            chunk_iters: 1,
            polish: false,
            ..QpOptions::default()
        };
        let result = solve_discrete(&inputs, &|| false);
        eprintln!(
            "[energy low-cap] status={} converged={} iterations recorded via status",
            result.status, result.converged
        );
        assert!(
            matches!(result.status, "max_iter" | "solved_inaccurate"),
            "expected a chunk-boundary status, got {}",
            result.status
        );
        assert_eq!(result.converged, result.status == "solved_inaccurate");
        assert!(
            result.weights.iter().all(|w| w.is_finite()),
            "capped solve returned non-finite weights"
        );
    }

    #[test]
    fn continuous_energy_reduces_the_distance_covariance() {
        // An exposure correlated with the covariates: reweighting lowers the
        // weighted distance covariance below the uniform value.
        let n = 10;
        let treat: Vec<f64> = (0..n).map(|i| i as f64 - 4.5).collect();
        let covs: Vec<f64> = (0..n)
            .map(|i| (i as f64 - 4.5) * 0.8 + ((i as f64) * 1.3).sin())
            .collect();
        let s = vec![1.0; n];
        let inputs = EnergyContInputs {
            covs: &covs,
            n,
            p: 1,
            treat: &treat,
            distance: Distance::ScaledEuclidean,
            s: &s,
            min_weight: 1e-8,
            lambda: 1e-4,
            dimension_adj: true,
            d_covs: &[],
            n_d_covs: 0,
            d_treat: &[],
            n_d_treat: 0,
            bal_covs: &[],
            n_bal: 0,
            bal_tols: &[],
            threads: 1,
            qp: QpOptions::default(),
        };
        let result = solve_cont(&inputs, &|| false);
        assert!(result.converged, "status {}", result.status);
        let weighted = oracle_dcov(&treat, &covs, n, &result.weights);
        let unit = vec![1.0; n];
        let unweighted = oracle_dcov(&treat, &covs, n, &unit);
        assert!(
            weighted < unweighted,
            "weighted {weighted} !< unweighted {unweighted}"
        );
    }

    // Weighted distance covariance in the V-statistic form of the test oracle.
    fn oracle_dcov(treat: &[f64], covs: &[f64], n: usize, w: &[f64]) -> f64 {
        let wsum: f64 = w.iter().sum();
        let p: Vec<f64> = w.iter().map(|&wi| wi / wsum).collect();
        let a = pairwise::euclidean(treat, n, 1, 1);
        let b = pairwise::euclidean(covs, n, 1, 1);
        let u: Vec<f64> = (0..n)
            .map(|i| (0..n).map(|j| a[j * n + i] * p[j]).sum())
            .collect();
        let v: Vec<f64> = (0..n)
            .map(|i| (0..n).map(|j| b[j * n + i] * p[j]).sum())
            .collect();
        let mut s1 = 0.0;
        for i in 0..n {
            for j in 0..n {
                s1 += p[i] * a[j * n + i] * b[j * n + i] * p[j];
            }
        }
        let pap: f64 = (0..n).map(|i| p[i] * u[i]).sum();
        let pbp: f64 = (0..n).map(|i| p[i] * v[i]).sum();
        let s2 = pap * pbp;
        let s3: f64 = (0..n).map(|i| p[i] * u[i] * v[i]).sum();
        s1 + s2 - 2.0 * s3
    }

    #[test]
    fn distance_assembly_is_deterministic_across_threads() {
        let n = 60;
        let covs: Vec<f64> = (0..n * 3).map(|k| ((k as f64) * 0.31).sin()).collect();
        let s = vec![1.0; n];
        let one = distance_matrix(&covs, n, 3, Distance::ScaledEuclidean, &s, 1);
        for threads in [2, 4] {
            let many = distance_matrix(&covs, n, 3, Distance::ScaledEuclidean, &s, threads);
            for (a, b) in one.iter().zip(&many) {
                assert_eq!(a.to_bits(), b.to_bits());
            }
        }
    }
}
