//! Entropy balancing: the maximum-entropy reweighting whose weighted covariate
//! moments match a set of targets.
//!
//! For a group of units with base weights `q_i` and sampling weights `s_i`, the
//! weights minimize the Kullback-Leibler divergence `sum_i q_i-relative entropy`
//! subject to `weighted mean of C_j = target_j`. The Lagrangian dual is the
//! smooth convex function
//!
//! ```text
//! F(beta) = log( sum_i s_i q_i exp(-C_i . beta) ) + target . beta
//! ```
//!
//! whose gradient is `target - weighted mean of C` and whose Hessian is the
//! weighted covariance of `C`. The minimizing weights are the exponential tilt
//! `w_i proportional to q_i exp(-C_i . beta)`, normalized so `sum_i s_i w_i`
//! equals the requested effective size. When any per-constraint tolerance is
//! positive the equality constraints become a box, whose dual adds a weighted L1
//! penalty on `beta`, and the problem is solved by FISTA instead of Newton.
//!
//! Discrete exposures are handled by looping independent per-group problems that
//! share one set of targets; a continuous exposure is the single-group case over
//! the whole sample.

use std::sync::Arc;

use faer::MatMut;
use faer::prelude::ReborrowMut;
use rayon::ThreadPool;
use rayon::iter::{IndexedParallelIterator, ParallelIterator};
use rayon::slice::ParallelSliceMut;

use crate::esteq::fista;
use crate::esteq::{self, EsteqProblem, SolveOptions, Solver};
use crate::threads::{deterministic_map_reduce, get_pool};

/// Dual problem for a single group of units.
///
/// `covs` is the full column-major `n` by `p` constraint matrix; `idx` selects
/// the group's rows. `base` and `s` are the length-`n` base and sampling
/// weights. Weights are scaled so `sum_i s_i w_i = n_eff` within the group.
struct EntropyProblem<'a> {
    covs: &'a [f64],
    n: usize,
    p: usize,
    idx: &'a [usize],
    targets: &'a [f64],
    base: &'a [f64],
    s: &'a [f64],
    n_eff: f64,
    pool: Arc<ThreadPool>,
}

/// Accumulated sufficient statistics for one dual evaluation.
struct Accum {
    z: f64,
    m: Vec<f64>,
    smat: Vec<f64>,
    /// Scratch for one unit's covariate row, reused across the fold so the
    /// column-major matrix is read down a unit only once per unit rather than
    /// once per Hessian entry.
    crow: Vec<f64>,
}

impl Accum {
    fn zeros(p: usize, want_smat: bool) -> Self {
        Self {
            z: 0.0,
            m: vec![0.0; p],
            smat: if want_smat {
                vec![0.0; p * p]
            } else {
                Vec::new()
            },
            crow: vec![0.0; p],
        }
    }
}

impl EntropyProblem<'_> {
    /// Linear predictor `C_i . beta` for global unit `gi`.
    fn lin(&self, gi: usize, beta: &[f64]) -> f64 {
        let mut acc = 0.0;
        for (j, &b) in beta.iter().enumerate() {
            acc += self.covs[j * self.n + gi] * b;
        }
        acc
    }

    /// Fold the group's tilted contributions into `Accum`. Fills the second
    /// moment matrix only when `want_smat`.
    fn accumulate(&self, beta: &[f64], want_smat: bool) -> Accum {
        let p = self.p;
        deterministic_map_reduce(
            &self.pool,
            self.idx.len(),
            || Accum::zeros(p, want_smat),
            |acc, local| {
                let gi = self.idx[local];
                // Gather this unit's covariate row once. The constraint matrix
                // is column-major, so reading down a unit strides by n; doing it
                // a single time turns the inner p-by-p Hessian accumulation into
                // unit-stride reads from the scratch row.
                for j in 0..p {
                    acc.crow[j] = self.covs[j * self.n + gi];
                }
                let lin: f64 = acc.crow.iter().zip(beta).map(|(c, b)| c * b).sum();
                let e = self.s[gi] * self.base[gi] * (-lin).exp();
                acc.z += e;
                for j in 0..p {
                    let ecj = e * acc.crow[j];
                    acc.m[j] += ecj;
                    if want_smat {
                        for k in 0..p {
                            acc.smat[j * p + k] += ecj * acc.crow[k];
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

    /// Per-unit weights and the achieved weighted means at `beta`.
    fn solution_parts(&self, beta: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let acc = self.accumulate(beta, false);
        let z = acc.z;
        let mbar: Vec<f64> = acc.m.iter().map(|mj| mj / z).collect();
        let weights: Vec<f64> = self
            .idx
            .iter()
            .map(|&gi| {
                let tilt = self.base[gi] * (-self.lin(gi, beta)).exp();
                self.n_eff * tilt / z
            })
            .collect();
        (weights, mbar)
    }
}

impl EsteqProblem for EntropyProblem<'_> {
    fn n_params(&self) -> usize {
        self.p
    }

    fn n_units(&self) -> usize {
        self.idx.len()
    }

    fn value(&self, beta: &[f64]) -> Option<f64> {
        let acc = self.accumulate(beta, false);
        let mut dot = 0.0;
        for (j, &b) in beta.iter().enumerate() {
            dot += self.targets[j] * b;
        }
        Some(acc.z.ln() + dot)
    }

    fn gradient(&self, beta: &[f64], g: &mut [f64]) {
        let acc = self.accumulate(beta, false);
        for (j, gj) in g.iter_mut().enumerate() {
            *gj = self.targets[j] - acc.m[j] / acc.z;
        }
    }

    fn hessian(&self, beta: &[f64], mut h: MatMut<'_, f64>) {
        let acc = self.accumulate(beta, true);
        let z = acc.z;
        for j in 0..self.p {
            let mbar_j = acc.m[j] / z;
            for k in 0..self.p {
                let mbar_k = acc.m[k] / z;
                *h.rb_mut().get_mut(j, k) = acc.smat[j * self.p + k] / z - mbar_j * mbar_k;
            }
        }
    }

    fn value_grad_hess(&self, beta: &[f64], g: &mut [f64], mut h: MatMut<'_, f64>) -> Option<f64> {
        let acc = self.accumulate(beta, true);
        let z = acc.z;
        let mut dot = 0.0;
        for j in 0..self.p {
            let mbar_j = acc.m[j] / z;
            g[j] = self.targets[j] - mbar_j;
            dot += self.targets[j] * beta[j];
            for k in 0..self.p {
                let mbar_k = acc.m[k] / z;
                *h.rb_mut().get_mut(j, k) = acc.smat[j * self.p + k] / z - mbar_j * mbar_k;
            }
        }
        Some(z.ln() + dot)
    }

    fn psi(&self, beta: &[f64], mut out: MatMut<'_, f64>) {
        let (weights, _mbar) = self.solution_parts(beta);
        for (local, &gi) in self.idx.iter().enumerate() {
            let sw = self.s[gi] * weights[local];
            for j in 0..self.p {
                *out.rb_mut().get_mut(local, j) =
                    sw * (self.covs[j * self.n + gi] - self.targets[j]);
            }
        }
    }
}

/// The solver used for a discrete or continuous entropy problem.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EntropySolver {
    /// Damped Newton (default) for the exact problem.
    Newton,
    /// L-BFGS through the basin backend.
    Lbfgs,
    /// L-BFGS warm start followed by a Newton polish.
    LbfgsThenNewton,
}

impl EntropySolver {
    fn as_esteq(self) -> Solver {
        match self {
            EntropySolver::Newton => Solver::Newton,
            EntropySolver::Lbfgs => Solver::Lbfgs,
            EntropySolver::LbfgsThenNewton => Solver::LbfgsThenNewton,
        }
    }
}

/// Solution of an entropy balancing problem across one or more groups.
pub struct EntropyResult {
    /// Balancing weights, one per unit, in the original unit order.
    pub weights: Vec<f64>,
    /// Dual variables, `p` per group, stacked in group order.
    pub duals: Vec<f64>,
    /// Whether every group met its convergence criterion.
    pub converged: bool,
    /// Whether any group's solve stopped because a user interrupt was pending.
    pub interrupted: bool,
    /// Largest iteration count across groups.
    pub iterations: usize,
    /// Largest gradient sup norm across groups.
    pub grad_norm: f64,
    /// Name of the solver that produced the result.
    pub solver: &'static str,
    /// Column-major `n` by `P` estimating functions, `None` for the inexact
    /// (tolerance) problem. `P` is `p` times the number of groups.
    pub psi: Option<Vec<f64>>,
    /// Column-major `P` by `P` block-diagonal Jacobian, `None` for the inexact
    /// problem.
    pub jac: Option<Vec<f64>>,
    /// Column-major `n` by `P` weight derivatives, `None` for the inexact
    /// problem.
    pub dw_dbeta: Option<Vec<f64>>,
}

/// Inputs shared by the discrete and continuous entrypoints.
pub struct EntropyInputs<'a> {
    /// Column-major `n` by `p` constraint matrix.
    pub covs: &'a [f64],
    /// Number of units.
    pub n: usize,
    /// Number of constraints.
    pub p: usize,
    /// Target weighted mean for each constraint.
    pub targets: &'a [f64],
    /// Per-constraint tolerance. Any positive entry selects the inexact
    /// problem; the tolerances become the L1 penalty weights.
    pub tols: &'a [f64],
    /// Base weights `q_i`.
    pub base: &'a [f64],
    /// Sampling weights `s_i`.
    pub s: &'a [f64],
    /// Target for `sum_i s_i w_i` within each group.
    pub n_eff: f64,
    /// Worker threads for the reductions.
    pub threads: usize,
    /// Maximum solver iterations.
    pub max_iter: usize,
    /// Convergence threshold on the gradient sup norm (Newton) and the relative
    /// loss change (FISTA).
    pub tol: f64,
    /// Solver for the exact problem.
    pub solver: EntropySolver,
}

/// The number of solved groups `group_idx` implies, one past its largest
/// non-negative entry.
fn group_count(group_idx: &[i32]) -> usize {
    group_idx.iter().copied().max().map_or(0, |g| g + 1).max(0) as usize
}

/// Partition the units into their groups, in unit order.
///
/// A negative index marks a unit that belongs to no solved group, the focal
/// level a focal estimand leaves out. Those units appear in no partition, so
/// every output the group loop fills keeps its initial value for them.
fn partition_groups(group_idx: &[i32], n_groups: usize) -> Vec<Vec<usize>> {
    let mut groups: Vec<Vec<usize>> = vec![Vec::new(); n_groups];
    for (i, &g) in group_idx.iter().enumerate() {
        if g >= 0 {
            groups[g as usize].push(i);
        }
    }
    groups
}

/// Solve a discrete entropy problem whose units are partitioned into groups by
/// `group_idx` (values `0..G`), each balanced to the shared `targets`.
pub fn solve_discrete(
    inputs: &EntropyInputs<'_>,
    group_idx: &[i32],
    interrupt: &dyn Fn() -> bool,
) -> EntropyResult {
    let groups = partition_groups(group_idx, group_count(group_idx));
    // Every constraint column can be relaxed to its tolerance; there are no
    // exact-only marginal columns in a discrete problem.
    solve_groups(inputs, &groups, inputs.tols, interrupt)
}

/// Solve a continuous entropy problem: a single group over all units, balanced
/// to `targets` under the supplied constraint matrix.
///
/// `dist_ind` marks the marginal-distribution constraint columns (a nonzero
/// entry). Those columns define the estimand's reference marginals and are held
/// exactly even when `tols` requests a relaxed fit, so their L1 penalty weight
/// is forced to zero regardless of the requested tolerance.
pub fn solve_continuous(
    inputs: &EntropyInputs<'_>,
    dist_ind: &[i32],
    interrupt: &dyn Fn() -> bool,
) -> EntropyResult {
    let l1: Vec<f64> = inputs
        .tols
        .iter()
        .zip(dist_ind)
        .map(|(&t, &d)| if d != 0 { 0.0 } else { t })
        .collect();
    solve_groups(inputs, &single_group(inputs.n), &l1, interrupt)
}

fn solve_groups(
    inputs: &EntropyInputs<'_>,
    groups: &[Vec<usize>],
    l1: &[f64],
    interrupt: &dyn Fn() -> bool,
) -> EntropyResult {
    let p = inputs.p;
    let n = inputs.n;
    let n_groups = groups.len();
    let total_params = n_groups * p;
    let inexact = l1.iter().any(|&t| t > 0.0);
    let pool = get_pool(inputs.threads);

    let mut weights = vec![0.0; n];
    let mut duals = vec![0.0; total_params];
    let mut converged = true;
    let mut interrupted = false;
    let mut iterations = 0;
    let mut grad_norm = 0.0_f64;

    let mut psi = (!inexact).then(|| vec![0.0; n * total_params]);
    let mut jac = (!inexact).then(|| vec![0.0; total_params * total_params]);
    let mut dw = (!inexact).then(|| vec![0.0; n * total_params]);

    let solve_opts = SolveOptions {
        max_iter: inputs.max_iter,
        grad_tol: inputs.tol,
        fista_rel_tol: inputs.tol,
    };

    let solver_name = if inexact {
        "fista"
    } else {
        match inputs.solver {
            EntropySolver::Newton => "newton",
            EntropySolver::Lbfgs => "lbfgs",
            EntropySolver::LbfgsThenNewton => "lbfgs_then_newton",
        }
    };

    for (g, idx) in groups.iter().enumerate() {
        if idx.is_empty() {
            continue;
        }
        let problem = EntropyProblem {
            covs: inputs.covs,
            n,
            p,
            idx,
            targets: inputs.targets,
            base: inputs.base,
            s: inputs.s,
            n_eff: inputs.n_eff,
            pool: Arc::clone(&pool),
        };
        let mut beta = vec![0.0; p];

        if inexact {
            let report = fista::minimize(
                &mut beta,
                l1,
                |b| problem.value(b).expect("entropy dual is smooth"),
                |b, gout| problem.gradient(b, gout),
                inputs.max_iter,
                inputs.tol,
                interrupt,
            );
            converged &= report.converged;
            interrupted |= report.interrupted;
            iterations = iterations.max(report.iterations);
            grad_norm = grad_norm.max(report.grad_norm);
        } else {
            let report = esteq::solve(
                &problem,
                &mut beta,
                inputs.solver.as_esteq(),
                &solve_opts,
                interrupt,
            );
            converged &= report.converged;
            interrupted |= report.interrupted;
            iterations = iterations.max(report.iterations);
            grad_norm = grad_norm.max(report.grad_norm);
        }

        duals[g * p..g * p + p].copy_from_slice(&beta);
        let (group_weights, mbar) = problem.solution_parts(&beta);
        for (local, &gi) in idx.iter().enumerate() {
            weights[gi] = group_weights[local];
        }

        if let (Some(psi), Some(jac), Some(dw)) = (psi.as_mut(), jac.as_mut(), dw.as_mut()) {
            fill_estimating_output(
                &pool,
                inputs,
                idx,
                g,
                &group_weights,
                &mbar,
                total_params,
                psi,
                jac,
                dw,
            );
        }
    }

    EntropyResult {
        weights,
        duals,
        converged: converged && !interrupted,
        interrupted,
        iterations,
        grad_norm,
        solver: solver_name,
        psi,
        jac,
        dw_dbeta: dw,
    }
}

/// Scratch and partial for the Jacobian block reduction.
struct JacAccum {
    mat: Vec<f64>,
    dev: Vec<f64>,
}

impl JacAccum {
    fn zeros(p: usize) -> Self {
        Self {
            mat: vec![0.0; p * p],
            dev: vec![0.0; p],
        }
    }
}

/// Fill the estimating-equation blocks for one group into the global arrays.
///
/// For unit `i` in group `g`, with parameter block starting at column `g * p`:
/// `psi_ij = s_i w_i (C_ij - target_j)`, `dw_ik = w_i (mbar_k - C_ik)`, and the
/// Jacobian block is `-sum_i s_i w_i (C_i - target)(C_i - target)'`.
///
/// The per-unit `psi` and weight-derivative entries are elementwise, so they
/// fill in parallel over the block's columns with no reduction. The Jacobian
/// block is a sum over units of a rank-one outer product and goes through
/// [`deterministic_map_reduce`] so it is bit-identical across thread counts.
#[allow(clippy::too_many_arguments)]
fn fill_estimating_output(
    pool: &ThreadPool,
    inputs: &EntropyInputs<'_>,
    idx: &[usize],
    g: usize,
    group_weights: &[f64],
    mbar: &[f64],
    total_params: usize,
    psi: &mut [f64],
    jac: &mut [f64],
    dw: &mut [f64],
) {
    let n = inputs.n;
    let p = inputs.p;
    let offset = g * p;

    // psi and dw are n by P, column-major; the block for this group occupies the
    // contiguous columns `offset .. offset + p`. Each output entry is written
    // exactly once and its value is independent of the order columns are filled,
    // so the block fills in parallel over columns without a reduction.
    let psi_block = &mut psi[offset * n..(offset + p) * n];
    let dw_block = &mut dw[offset * n..(offset + p) * n];
    pool.install(|| {
        psi_block
            .par_chunks_mut(n)
            .zip(dw_block.par_chunks_mut(n))
            .enumerate()
            .for_each(|(j, (psi_col, dw_col))| {
                let target_j = inputs.targets[j];
                let mbar_j = mbar[j];
                let covs_col = &inputs.covs[j * n..(j + 1) * n];
                for (local, &gi) in idx.iter().enumerate() {
                    let w = group_weights[local];
                    let c = covs_col[gi];
                    psi_col[gi] = inputs.s[gi] * w * (c - target_j);
                    dw_col[gi] = w * (mbar_j - c);
                }
            });
    });

    let block = deterministic_map_reduce(
        pool,
        idx.len(),
        || JacAccum::zeros(p),
        |acc, local| {
            let gi = idx[local];
            let sw = inputs.s[gi] * group_weights[local];
            for j in 0..p {
                acc.dev[j] = inputs.covs[j * n + gi] - inputs.targets[j];
            }
            for j in 0..p {
                let sdev = sw * acc.dev[j];
                for k in 0..p {
                    acc.mat[j * p + k] += sdev * acc.dev[k];
                }
            }
        },
        |acc, other| {
            for i in 0..p * p {
                acc.mat[i] += other.mat[i];
            }
        },
    );
    for j in 0..p {
        let row = offset + j;
        for k in 0..p {
            let col = offset + k;
            jac[col * total_params + row] -= block.mat[j * p + k];
        }
    }
}

/// Apply a per-group scalar to the estimating-equation output in place.
///
/// The R layer renormalizes each solved group's weights to the estimand's
/// target sum after the solve. That scalar multiplies the group's `psi`
/// columns, weight-derivative columns, and Jacobian block, which the core
/// computed at its own within-group normalization. Performing the multiply here,
/// before the matrices cross the FFI boundary, keeps the large `n` by `P`
/// matrices from being copied a second time on the R side. `scales` holds one
/// value per parameter block; a block scaled by exactly one is left untouched.
///
/// The inexact problem carries no estimating equations, so the call is a no-op.
/// Otherwise the length of `scales` must match the number of parameter blocks,
/// `duals.len() / p`; a mismatch would silently mis-scale the output and is
/// returned as an error for the boundary layer to surface. A solve with no
/// parameters at all, from an empty covariate matrix or from a design whose every
/// unit sits outside the solved groups, is rejected on the same path: the total
/// parameter count is also the divisor the row count is recovered from below.
pub fn scale_estimating_output(
    result: &mut EntropyResult,
    p: usize,
    scales: &[f64],
) -> Result<(), String> {
    if result.psi.is_none() {
        return Ok(());
    }
    let total_params = result.duals.len();
    if p == 0 || total_params == 0 || scales.len() * p != total_params {
        return Err(format!(
            "estimating-equation scale has length {} but the solve has {} parameter block(s) of size {p}",
            scales.len(),
            total_params / p.max(1)
        ));
    }
    let (Some(psi), Some(jac), Some(dw)) = (
        result.psi.as_mut(),
        result.jac.as_mut(),
        result.dw_dbeta.as_mut(),
    ) else {
        return Ok(());
    };
    let n = psi.len() / total_params;
    for (g, &scale) in scales.iter().enumerate() {
        if scale == 1.0 {
            continue;
        }
        let cols = (g * p)..(g * p + p);
        for col in cols.clone() {
            for row in 0..n {
                psi[col * n + row] *= scale;
                dw[col * n + row] *= scale;
            }
        }
        for row in cols.clone() {
            for col in cols.clone() {
                jac[col * total_params + row] *= scale;
            }
        }
    }
    Ok(())
}

/// Check that the inputs the re-evaluation entrypoints read are as long as the
/// dimensions they are read at. Every entrypoint indexes the constraint matrix
/// at `n` by `p` and the weight vectors at `n`, so a short input would read past
/// the end of a slice. It is returned for the boundary layer to surface rather
/// than left to panic.
fn check_eval_dimensions(inputs: &EntropyInputs<'_>) -> Result<(), String> {
    let (n, p) = (inputs.n, inputs.p);
    if p == 0 {
        return Err("a re-evaluation needs at least one constraint".to_string());
    }
    if inputs.covs.len() != n * p {
        return Err(format!(
            "covs has {} elements but n * p = {n} * {p} = {}",
            inputs.covs.len(),
            n * p
        ));
    }
    if inputs.targets.len() != p {
        return Err(format!(
            "targets has {} element(s) but the constraint matrix has {p} column(s)",
            inputs.targets.len()
        ));
    }
    if inputs.base.len() != n || inputs.s.len() != n {
        return Err(format!(
            "the base and sampling weight vectors must both have length {n} (got {} and {})",
            inputs.base.len(),
            inputs.s.len()
        ));
    }
    Ok(())
}

/// Check that `coefs` carries one `p`-vector per group and `scales` one value
/// per group, the shared precondition of the re-evaluation entrypoints. A
/// mismatch would slice past the end of the dual vector, so it is returned for
/// the boundary layer to surface rather than left to panic.
fn check_discrete_args(
    p: usize,
    n_groups: usize,
    coefs: &[f64],
    scales: &[f64],
) -> Result<(), String> {
    let total_params = n_groups * p;
    if coefs.len() != total_params || scales.len() < n_groups {
        return Err(format!(
            "group_idx implies {n_groups} group(s) of size {p}, so coefs must have length {total_params} (got {}) and scales at least {n_groups} (got {})",
            coefs.len(),
            scales.len()
        ));
    }
    Ok(())
}

/// The per-unit weights a set of duals implies, at the reported scale.
///
/// Each group's exponential tilt is normalized so its sampling-weighted total
/// is `n_eff`, then multiplied by that group's entry of `scales`, the same
/// per-group renormalization [`scale_estimating_output`] applies to the
/// estimating-equation blocks. A unit in no group keeps zero, the value
/// [`solve_groups`] leaves for the focal level a focal estimand omits. A
/// continuous problem is the single-group case, so it reads this through one
/// group holding every unit.
fn tilt_weights(
    inputs: &EntropyInputs<'_>,
    groups: &[Vec<usize>],
    coefs: &[f64],
    scales: &[f64],
) -> Vec<f64> {
    let p = inputs.p;
    let pool = get_pool(inputs.threads);
    let mut weights = vec![0.0; inputs.n];
    for (g, idx) in groups.iter().enumerate() {
        if idx.is_empty() {
            continue;
        }
        let beta = &coefs[g * p..g * p + p];
        let problem = EntropyProblem {
            covs: inputs.covs,
            n: inputs.n,
            p,
            idx,
            targets: inputs.targets,
            base: inputs.base,
            s: inputs.s,
            n_eff: inputs.n_eff,
            pool: Arc::clone(&pool),
        };
        let (group_weights, _mbar) = problem.solution_parts(beta);
        let scale = scales.get(g).copied().unwrap_or(1.0);
        for (&gi, w) in idx.iter().zip(group_weights) {
            weights[gi] = w * scale;
        }
    }
    weights
}

/// The per-unit estimating functions a set of duals implies, at the scale
/// [`fill_estimating_output`] and [`scale_estimating_output`] leave in the
/// solve's own output: `psi_ij = s_i w_i (C_ij - target_j)`, with `w` the tilt
/// this group's entry of `scales` has already been applied to. The block for
/// group `g` occupies columns `g * p .. g * p + p`, and a unit in no group
/// contributes a zero row.
fn tilt_psi(
    inputs: &EntropyInputs<'_>,
    groups: &[Vec<usize>],
    coefs: &[f64],
    scales: &[f64],
) -> Vec<f64> {
    let n = inputs.n;
    let p = inputs.p;

    // The weights already carry the sampling weight's partner scale, so each
    // unit's row is its sampling weight times its scaled weight times the
    // centered constraint row.
    let weights = tilt_weights(inputs, groups, coefs, scales);
    let mut psi = vec![0.0; n * groups.len() * p];
    for (g, idx) in groups.iter().enumerate() {
        let offset = g * p;
        for &gi in idx {
            let sw = inputs.s[gi] * weights[gi];
            for j in 0..p {
                psi[(offset + j) * n + gi] = sw * (inputs.covs[j * n + gi] - inputs.targets[j]);
            }
        }
    }
    psi
}

/// Re-evaluate the discrete estimating functions at a supplied set of duals.
///
/// The estimating-equations container stores `psi` at the solution; a sandwich
/// variance that needs a finite difference must re-evaluate it at perturbed
/// parameters. This recomputes the `n` by `P` matrix from `coefs` (`p` duals per
/// group, stacked in group order) without solving, so the link math is never
/// reimplemented in R. `group_idx` partitions the units exactly as
/// [`solve_discrete`] does, and `scales` applies the same per-group
/// renormalization the solve output carried, one value per group. Units with a
/// negative group index contribute a zero row, matching the block structure of
/// the stored matrix.
pub fn eval_psi_discrete(
    inputs: &EntropyInputs<'_>,
    group_idx: &[i32],
    coefs: &[f64],
    scales: &[f64],
) -> Result<Vec<f64>, String> {
    let n_groups = group_count(group_idx);
    check_eval_dimensions(inputs)?;
    check_discrete_args(inputs.p, n_groups, coefs, scales)?;
    let groups = partition_groups(group_idx, n_groups);
    Ok(tilt_psi(inputs, &groups, coefs, scales))
}

/// Re-evaluate the discrete balancing weights at a supplied set of duals.
///
/// A sandwich variance that treats the weights as a function of the duals needs
/// the weight map itself, not only its derivative at the solution. This
/// recomputes the length-`n` weight vector from `coefs` (`p` duals per group,
/// stacked in group order) without solving, at the scale
/// [`solve_discrete`] reports: each group's exponential tilt normalized so its
/// sampling-weighted total is `n_eff`, multiplied by that group's entry of
/// `scales`, the same per-group renormalization
/// [`scale_estimating_output`] applies to the estimating-equation blocks. Units
/// with a negative group index belong to no solved group and carry zero,
/// matching the solve output.
pub fn eval_weights_discrete(
    inputs: &EntropyInputs<'_>,
    group_idx: &[i32],
    coefs: &[f64],
    scales: &[f64],
) -> Result<Vec<f64>, String> {
    let n_groups = group_count(group_idx);
    check_eval_dimensions(inputs)?;
    check_discrete_args(inputs.p, n_groups, coefs, scales)?;
    let groups = partition_groups(group_idx, n_groups);
    Ok(tilt_weights(inputs, &groups, coefs, scales))
}

/// Check the precondition of the continuous re-evaluation entrypoints: the
/// inputs are as long as their dimensions, and `coefs` carries the single dual
/// block the continuous problem solves, one value per constraint. A mismatch
/// would read past the end of a slice, so it is returned for the boundary layer
/// to surface rather than left to panic.
fn check_continuous_args(inputs: &EntropyInputs<'_>, coefs: &[f64]) -> Result<(), String> {
    check_eval_dimensions(inputs)?;
    if coefs.len() != inputs.p {
        return Err(format!(
            "a continuous solve carries one dual per constraint, so coefs must have length {} (got {})",
            inputs.p,
            coefs.len()
        ));
    }
    Ok(())
}

/// The whole sample as a single group, the partition a continuous problem
/// solves over.
fn single_group(n: usize) -> Vec<Vec<usize>> {
    vec![(0..n).collect()]
}

/// Re-evaluate the continuous estimating functions at a supplied set of duals.
///
/// A continuous exposure is the single-group case of the same tilt: one dual
/// block of length `p` over the constraint matrix [`solve_continuous`] is given,
/// whose marginal and product columns and whose targets the R layer builds. The
/// re-evaluation solves nothing; it recomputes the `n` by `p` matrix at `coefs`
/// so a sandwich variance can finite-difference the stored Jacobian without the
/// tilt math being reimplemented in R.
///
/// Scale convention, the one the solve's own container reports: `psi` is at the
/// raw M-estimator scale, `psi_ij = s_i w_i (C_ij - target_j)`, where `s` holds
/// the sampling weights and `w` is the exponential tilt normalized so
/// `sum_i s_i w_i = n_eff`, anchored to the base measure `s_i q_i` the discrete
/// path composes as well. `scale` is then the single group's renormalization,
/// the ratio between the weights the fit reports and that normalization, applied
/// exactly as [`scale_estimating_output`] applies it to the stored blocks. A
/// continuous fit reports its weights at `n_eff`, so the R layer passes one; the
/// argument keeps the returned `psi` on the container's scale should that
/// reporting rule ever change.
///
/// The tolerances and the marginal-column indicator play no part: they shape the
/// penalty the solve minimizes, and the inexact problem carries no estimating
/// equations at all.
pub fn eval_psi_continuous(
    inputs: &EntropyInputs<'_>,
    coefs: &[f64],
    scale: f64,
) -> Result<Vec<f64>, String> {
    check_continuous_args(inputs, coefs)?;
    Ok(tilt_psi(inputs, &single_group(inputs.n), coefs, &[scale]))
}

/// Re-evaluate the continuous balancing weights at a supplied set of duals.
///
/// The weight map itself, not only its derivative at the solution, is what a
/// sandwich variance treating the weights as a function of the duals needs. This
/// recomputes the length-`n` vector at `coefs` without solving, at the scale
/// [`solve_continuous`] reports: the single group's exponential tilt normalized
/// so its sampling-weighted total is `n_eff`, multiplied by `scale`, the same
/// renormalization [`eval_psi_continuous`] applies to the estimating functions.
pub fn eval_weights_continuous(
    inputs: &EntropyInputs<'_>,
    coefs: &[f64],
    scale: f64,
) -> Result<Vec<f64>, String> {
    check_continuous_args(inputs, coefs)?;
    Ok(tilt_weights(
        inputs,
        &single_group(inputs.n),
        coefs,
        &[scale],
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn no_interrupt() -> impl Fn() -> bool {
        || false
    }

    // The exact dual solved by FISTA with no L1 penalty must reach the same
    // point as Newton: the tolerance-zero problem is the exact problem.
    #[test]
    fn fista_matches_newton_with_zero_penalty() {
        let covs = vec![1.0, 1.0, 1.0, 0.0, 0.0, 0.4, 0.1, 0.9, 0.2, 0.7];
        let n = 5;
        let p = 2;
        let targets = vec![0.7, 0.5];
        let base = vec![1.0; n];
        let s = vec![1.0; n];
        let idx: Vec<usize> = (0..n).collect();
        let pool = get_pool(1);
        let problem = EntropyProblem {
            covs: &covs,
            n,
            p,
            idx: &idx,
            targets: &targets,
            base: &base,
            s: &s,
            n_eff: n as f64,
            pool,
        };

        let mut beta_newton = vec![0.0; p];
        esteq::solve(
            &problem,
            &mut beta_newton,
            Solver::Newton,
            &SolveOptions::default(),
            &no_interrupt(),
        );

        let mut beta_fista = vec![0.0; p];
        fista::minimize(
            &mut beta_fista,
            &[0.0, 0.0],
            |b| problem.value(b).unwrap(),
            |b, g| problem.gradient(b, g),
            20000,
            1e-14,
            &no_interrupt(),
        );

        for j in 0..p {
            assert!(
                (beta_newton[j] - beta_fista[j]).abs() < 1e-6,
                "beta[{j}]: newton {} vs fista {}",
                beta_newton[j],
                beta_fista[j]
            );
        }
    }

    // The in-place rescale multiplies each block's psi columns, weight-derivative
    // columns, and Jacobian block by that block's scalar, leaving other blocks
    // and a unit scalar untouched.
    fn sample_result(n: usize, p: usize, n_groups: usize) -> EntropyResult {
        let total = n_groups * p;
        let psi: Vec<f64> = (0..n * total).map(|i| i as f64 + 1.0).collect();
        let dw: Vec<f64> = (0..n * total).map(|i| (i as f64 + 1.0) * 0.5).collect();
        let jac: Vec<f64> = (0..total * total).map(|i| i as f64 + 1.0).collect();
        EntropyResult {
            weights: vec![0.0; n],
            duals: vec![0.0; total],
            converged: true,
            interrupted: false,
            iterations: 0,
            grad_norm: 0.0,
            solver: "newton",
            psi: Some(psi),
            jac: Some(jac),
            dw_dbeta: Some(dw),
        }
    }

    #[test]
    fn scale_estimating_output_scales_each_block() {
        let (n, p, n_groups) = (4, 2, 2);
        let total = n_groups * p;
        let base = sample_result(n, p, n_groups);
        let mut scaled = sample_result(n, p, n_groups);
        let scales = [3.0, 5.0];
        scale_estimating_output(&mut scaled, p, &scales).unwrap();

        let (bp, sp) = (base.psi.unwrap(), scaled.psi.unwrap());
        let (bd, sd) = (base.dw_dbeta.unwrap(), scaled.dw_dbeta.unwrap());
        for (g, &scale) in scales.iter().enumerate().take(n_groups) {
            for col in (g * p)..(g * p + p) {
                for row in 0..n {
                    let k = col * n + row;
                    assert_eq!(sp[k], bp[k] * scale);
                    assert_eq!(sd[k], bd[k] * scale);
                }
            }
        }

        let (bj, sj) = (base.jac.unwrap(), scaled.jac.unwrap());
        for (g, &scale) in scales.iter().enumerate().take(n_groups) {
            for row in (g * p)..(g * p + p) {
                for col in (g * p)..(g * p + p) {
                    let k = col * total + row;
                    assert_eq!(sj[k], bj[k] * scale);
                }
            }
        }
    }

    #[test]
    fn scale_estimating_output_leaves_unit_blocks_untouched() {
        let (n, p, n_groups) = (3, 2, 2);
        let base = sample_result(n, p, n_groups);
        let mut scaled = sample_result(n, p, n_groups);
        // The first block is scaled; the second is left exactly as computed.
        scale_estimating_output(&mut scaled, p, &[2.0, 1.0]).unwrap();

        let (bp, sp) = (base.psi.unwrap(), scaled.psi.unwrap());
        for col in p..(2 * p) {
            for row in 0..n {
                let k = col * n + row;
                assert_eq!(sp[k].to_bits(), bp[k].to_bits());
            }
        }
    }

    #[test]
    fn scale_estimating_output_errors_on_length_mismatch() {
        // Two parameter blocks of size p, but only one scale supplied.
        let mut result = sample_result(4, 2, 2);
        let err = scale_estimating_output(&mut result, 2, &[3.0]).unwrap_err();
        assert!(err.contains("length 1"), "unexpected message: {err}");
        // The output is left untouched when the scale is rejected.
        let base = sample_result(4, 2, 2);
        assert_eq!(result.psi.unwrap(), base.psi.unwrap());
    }

    #[test]
    fn scale_estimating_output_errors_on_an_empty_parameter_set() {
        // A solve that produced no parameter blocks has no scale to apply. The
        // total parameter count is also the divisor the row count is recovered
        // from, so an empty set is rejected rather than divided by.
        let mut result = sample_result(4, 2, 0);
        let err = scale_estimating_output(&mut result, 2, &[]).unwrap_err();
        assert!(err.contains("length 0"), "unexpected message: {err}");
    }

    #[test]
    fn scale_estimating_output_ignores_scale_without_estimating_equations() {
        // The inexact problem carries no psi/jac/dw; a scale of any length is a
        // no-op rather than an error.
        let mut result = sample_result(3, 2, 1);
        result.psi = None;
        result.jac = None;
        result.dw_dbeta = None;
        scale_estimating_output(&mut result, 2, &[1.0, 2.0, 3.0]).unwrap();
    }

    // A two-group discrete design over two covariates, balanced to the pooled
    // mean. Entropy balancing carries no intercept column; both group hulls
    // contain the target, so the solve converges to exact balance.
    fn two_group_design() -> (Vec<f64>, Vec<i32>, Vec<f64>, usize, usize) {
        let n = 12;
        let p = 2;
        // Column-major n by p covariate matrix.
        let x1 = [
            0.3, 1.1, -0.4, 0.8, -1.0, 0.5, 1.3, -0.7, 0.2, -0.9, 0.6, -0.2,
        ];
        let x2 = [
            -0.5, 0.7, 1.2, -0.8, 0.4, -1.1, 0.9, 0.1, -0.6, 1.0, -0.3, 0.5,
        ];
        let mut covs = Vec::with_capacity(n * p);
        covs.extend_from_slice(&x1);
        covs.extend_from_slice(&x2);
        let group_idx = vec![0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1];
        // Pooled (unweighted) target of each column.
        let mut targets = vec![0.0; p];
        for j in 0..p {
            targets[j] = (0..n).map(|i| covs[j * n + i]).sum::<f64>() / n as f64;
        }
        (covs, group_idx, targets, n, p)
    }

    fn eval_inputs<'a>(
        covs: &'a [f64],
        targets: &'a [f64],
        base: &'a [f64],
        s: &'a [f64],
        n: usize,
        p: usize,
        tols: &'a [f64],
    ) -> EntropyInputs<'a> {
        EntropyInputs {
            covs,
            n,
            p,
            targets,
            tols,
            base,
            s,
            n_eff: 1.0,
            threads: 1,
            max_iter: 200,
            tol: 1e-12,
            solver: EntropySolver::Newton,
        }
    }

    // Re-evaluating psi at the solved duals reproduces the solve's own psi to
    // floating-point tolerance, and a central finite difference of the column
    // sums reproduces the analytic Jacobian.
    #[test]
    fn eval_psi_discrete_matches_solve_and_jacobian() {
        let (covs, group_idx, targets, n, p) = two_group_design();
        let base = vec![1.0; n];
        let s = vec![1.0; n];
        let tols = vec![0.0; p];
        let inputs = eval_inputs(&covs, &targets, &base, &s, n, p, &tols);
        let result = solve_discrete(&inputs, &group_idx, &no_interrupt());
        assert!(result.converged);

        let duals = &result.duals;
        let total_params = duals.len();
        let scales = vec![1.0; total_params / p];
        let stored = result.psi.as_ref().expect("exact problem has psi");
        let recomputed = eval_psi_discrete(&inputs, &group_idx, duals, &scales).unwrap();
        for (a, b) in stored.iter().zip(&recomputed) {
            assert!((a - b).abs() < 1e-10, "psi mismatch: {a} vs {b}");
        }

        let jac = result.jac.as_ref().expect("exact problem has jac");
        let eps = 1e-6;
        for col in 0..total_params {
            let mut up = duals.clone();
            let mut down = duals.clone();
            up[col] += eps;
            down[col] -= eps;
            let psi_up = eval_psi_discrete(&inputs, &group_idx, &up, &scales).unwrap();
            let psi_down = eval_psi_discrete(&inputs, &group_idx, &down, &scales).unwrap();
            for row in 0..total_params {
                let colsum_up: f64 = (0..n).map(|i| psi_up[row * n + i]).sum();
                let colsum_down: f64 = (0..n).map(|i| psi_down[row * n + i]).sum();
                let fd = (colsum_up - colsum_down) / (2.0 * eps);
                let analytic = jac[col * total_params + row];
                assert!(
                    (fd - analytic).abs() < 1e-4,
                    "jac[{row},{col}]: fd {fd} analytic {analytic}"
                );
            }
        }
    }

    // Non-unit sampling weights for the twelve-unit design. The tilt reads them
    // in its normalizing constant, so the weight map, not only the solved
    // duals, depends on them.
    fn sampling_weights() -> Vec<f64> {
        vec![0.7, 1.3, 0.5, 1.8, 1.1, 0.9, 1.4, 0.6, 1.2, 0.8, 1.5, 1.0]
    }

    // Re-evaluating the weights at the solved duals reproduces the solve's own
    // weight vector, and a central finite difference of that vector reproduces
    // the analytic weight derivative, under unit and non-unit sampling weights.
    #[test]
    fn eval_weights_discrete_matches_solve_and_weight_jacobian() {
        let (covs, group_idx, targets, n, p) = two_group_design();
        let base = vec![1.0; n];
        let tols = vec![0.0; p];
        for s in [vec![1.0; n], sampling_weights()] {
            let inputs = eval_inputs(&covs, &targets, &base, &s, n, p, &tols);
            let result = solve_discrete(&inputs, &group_idx, &no_interrupt());
            assert!(result.converged);

            let duals = &result.duals;
            let total_params = duals.len();
            let scales = vec![1.0; total_params / p];
            let dw = result.dw_dbeta.as_ref().expect("exact problem has dw");

            let recomputed = eval_weights_discrete(&inputs, &group_idx, duals, &scales).unwrap();
            assert_eq!(recomputed.len(), n);
            for (i, (a, b)) in result.weights.iter().zip(&recomputed).enumerate() {
                assert!(
                    (a - b).abs() < 1e-12,
                    "weight[{i}]: solve {a} versus eval {b}"
                );
            }

            let eps = 1e-6;
            for col in 0..total_params {
                let mut up = duals.clone();
                let mut down = duals.clone();
                up[col] += eps;
                down[col] -= eps;
                let w_up = eval_weights_discrete(&inputs, &group_idx, &up, &scales).unwrap();
                let w_down = eval_weights_discrete(&inputs, &group_idx, &down, &scales).unwrap();
                for i in 0..n {
                    let fd = (w_up[i] - w_down[i]) / (2.0 * eps);
                    let analytic = dw[col * n + i];
                    // The group-normalized weights are of order one over the
                    // group size, so the bound sits well below the derivative
                    // scale while staying orders above central-difference noise.
                    assert!(
                        (fd - analytic).abs() < 1e-6,
                        "dw[{i},{col}]: fd {fd} analytic {analytic}"
                    );
                }
            }
        }
    }

    // The per-group scale is the renormalization the reported weights carry, so
    // it multiplies the returned weight of every unit in its group and, being
    // free of the duals, its derivative too.
    #[test]
    fn eval_weights_discrete_applies_the_group_scale() {
        let (covs, group_idx, targets, n, p) = two_group_design();
        let base = vec![1.0; n];
        let s = vec![1.0; n];
        let tols = vec![0.0; p];
        let inputs = eval_inputs(&covs, &targets, &base, &s, n, p, &tols);
        let result = solve_discrete(&inputs, &group_idx, &no_interrupt());
        assert!(result.converged);

        let duals = &result.duals;
        let total_params = duals.len();
        let dw = result.dw_dbeta.as_ref().expect("exact problem has dw");
        let scales = vec![2.0, 3.0];

        let scaled = eval_weights_discrete(&inputs, &group_idx, duals, &scales).unwrap();
        for i in 0..n {
            let scale = scales[group_idx[i] as usize];
            let want = result.weights[i] * scale;
            assert!(
                (scaled[i] - want).abs() < 1e-12,
                "weight[{i}]: scaled {} expected {want}",
                scaled[i]
            );
        }

        let eps = 1e-6;
        for col in 0..total_params {
            let mut up = duals.clone();
            let mut down = duals.clone();
            up[col] += eps;
            down[col] -= eps;
            let w_up = eval_weights_discrete(&inputs, &group_idx, &up, &scales).unwrap();
            let w_down = eval_weights_discrete(&inputs, &group_idx, &down, &scales).unwrap();
            for i in 0..n {
                let fd = (w_up[i] - w_down[i]) / (2.0 * eps);
                let analytic = dw[col * n + i] * scales[group_idx[i] as usize];
                assert!(
                    (fd - analytic).abs() < 1e-6,
                    "dw[{i},{col}]: fd {fd} analytic {analytic}"
                );
            }
        }
    }

    // A focal estimand excludes its focal group from the solve, marking those
    // units with a negative group index. They carry no solved weight, and the
    // re-evaluation reproduces that zero alongside the solved group's weights.
    #[test]
    fn eval_weights_discrete_zeroes_the_unsolved_group() {
        let (covs, group_idx, targets, n, p) = two_group_design();
        // Shift the group labels down one so the first group becomes the focal
        // group, left out of the solve, and the second is solved on its own.
        let focal_idx: Vec<i32> = group_idx.iter().map(|&g| g - 1).collect();
        let base = vec![1.0; n];
        let s = vec![1.0; n];
        let tols = vec![0.0; p];
        let inputs = eval_inputs(&covs, &targets, &base, &s, n, p, &tols);
        let result = solve_discrete(&inputs, &focal_idx, &no_interrupt());
        assert!(result.converged);

        let scales = vec![1.0; result.duals.len() / p];
        let recomputed =
            eval_weights_discrete(&inputs, &focal_idx, &result.duals, &scales).unwrap();
        assert_eq!(recomputed.len(), n);
        for i in 0..n {
            if focal_idx[i] < 0 {
                assert_eq!(recomputed[i], 0.0, "unsolved unit {i} carries a weight");
            } else {
                assert!(recomputed[i] > 0.0, "solved unit {i} carries no weight");
            }
            assert!(
                (result.weights[i] - recomputed[i]).abs() < 1e-12,
                "weight[{i}]: solve {} versus eval {}",
                result.weights[i],
                recomputed[i]
            );
        }
    }

    // group_idx that implies more groups than coefs and scales cover is reported
    // rather than panicking on the block slice, so the R hook surfaces a
    // condition instead of a crash.
    #[test]
    fn eval_psi_discrete_rejects_inconsistent_group_count() {
        let (covs, group_idx, targets, n, p) = two_group_design();
        let base = vec![1.0; n];
        let s = vec![1.0; n];
        let tols = vec![0.0; p];
        let inputs = eval_inputs(&covs, &targets, &base, &s, n, p, &tols);
        // two_group_design implies two groups, so 2 * p coefficients and two
        // scales are required; a single group's worth must be rejected.
        let coefs = vec![0.0; p];
        let scales = vec![1.0];
        let err = eval_psi_discrete(&inputs, &group_idx, &coefs, &scales).unwrap_err();
        assert!(err.contains("group(s)"), "message was: {err}");
    }

    // The weight re-evaluation shares the argument check, so it reports the same
    // condition rather than slicing past the block.
    #[test]
    fn eval_weights_discrete_rejects_inconsistent_group_count() {
        let (covs, group_idx, targets, n, p) = two_group_design();
        let base = vec![1.0; n];
        let s = vec![1.0; n];
        let tols = vec![0.0; p];
        let inputs = eval_inputs(&covs, &targets, &base, &s, n, p, &tols);
        let coefs = vec![0.0; p];
        let scales = vec![1.0];
        let err = eval_weights_discrete(&inputs, &group_idx, &coefs, &scales).unwrap_err();
        assert!(err.contains("group(s)"), "message was: {err}");
    }

    /// Center and scale to the sample standard deviation, the standardization
    /// the R layer applies to a continuous exposure before it crosses the
    /// boundary.
    fn standardized(x: &[f64]) -> Vec<f64> {
        let n = x.len() as f64;
        let mean = x.iter().sum::<f64>() / n;
        let var = x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
        let sd = var.sqrt();
        x.iter().map(|v| (v - mean) / sd).collect()
    }

    // A continuous design in the shape the R layer builds: a standardized
    // exposure and two covariates as marginal columns, then the
    // exposure-covariate products. Each marginal column is held at its sample
    // mean and each product column at zero, the association the method removes.
    // The targets sit inside the hull of the constraint rows, so the single
    // group's tilt converges to exact balance.
    fn continuous_design() -> (Vec<f64>, Vec<f64>, Vec<i32>, usize, usize) {
        let n = 14;
        let p = 5;
        let x1 = [
            -1.0, -0.3, 0.3, -1.2, 0.2, 0.0, 0.1, 1.1, -1.2, 1.3, -0.7, -1.1, -0.7, 0.3,
        ];
        let x2 = [
            0.2, -0.3, -1.0, -0.6, 1.2, 0.2, -0.6, -0.9, -0.2, -1.7, -0.5, -0.7, 1.2, 1.0,
        ];
        let exposure = [
            -0.1, -1.1, 0.9, 0.9, 0.7, 0.7, -0.4, 0.7, 1.3, 0.0, -1.0, 0.8, 0.8, -0.3,
        ];
        let e = standardized(&exposure);

        let mut covs = Vec::with_capacity(n * p);
        covs.extend_from_slice(&e);
        covs.extend_from_slice(&x1);
        covs.extend_from_slice(&x2);
        covs.extend(x1.iter().zip(&e).map(|(x, ei)| x * ei));
        covs.extend(x2.iter().zip(&e).map(|(x, ei)| x * ei));

        let column_mean = |j: usize| covs[j * n..(j + 1) * n].iter().sum::<f64>() / n as f64;
        let targets = vec![column_mean(0), column_mean(1), column_mean(2), 0.0, 0.0];
        // The three marginal columns carry the reference distribution; the two
        // product columns are the relaxable association constraints.
        let dist_ind = vec![1, 1, 1, 0, 0];
        (covs, targets, dist_ind, n, p)
    }

    // The continuous fit normalizes its single group to the sampling-weight
    // total and reports it there, so the inputs carry that `n_eff`.
    fn continuous_inputs<'a>(
        covs: &'a [f64],
        targets: &'a [f64],
        base: &'a [f64],
        s: &'a [f64],
        n: usize,
        p: usize,
        tols: &'a [f64],
    ) -> EntropyInputs<'a> {
        let mut inputs = eval_inputs(covs, targets, base, s, n, p, tols);
        inputs.n_eff = s.iter().sum();
        inputs
    }

    // Non-unit sampling and base weights for the fourteen-unit continuous
    // design. Their product is the base measure the tilt anchors to, so both
    // move the solved duals and the weight map.
    fn continuous_sampling_weights() -> Vec<f64> {
        vec![
            0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,
        ]
    }

    fn continuous_base_weights() -> Vec<f64> {
        vec![
            1.3, 1.2, 1.2, 1.1, 1.1, 1.0, 1.0, 0.9, 0.9, 0.8, 0.8, 0.7, 0.7, 0.6,
        ]
    }

    // Re-evaluating psi at the solved duals reproduces the solve's own psi to
    // floating-point tolerance, its column sums vanish because the duals solve
    // the estimating equations, and a central finite difference of those column
    // sums reproduces the analytic Jacobian.
    #[test]
    fn eval_psi_continuous_matches_solve_and_jacobian() {
        let (covs, targets, dist_ind, n, p) = continuous_design();
        let tols = vec![0.0; p];
        for (s, base) in [
            (vec![1.0; n], vec![1.0; n]),
            (continuous_sampling_weights(), continuous_base_weights()),
        ] {
            let inputs = continuous_inputs(&covs, &targets, &base, &s, n, p, &tols);
            let result = solve_continuous(&inputs, &dist_ind, &no_interrupt());
            assert!(result.converged);
            assert_eq!(result.duals.len(), p);

            let duals = &result.duals;
            let stored = result.psi.as_ref().expect("exact problem has psi");
            let recomputed = eval_psi_continuous(&inputs, duals, 1.0).unwrap();
            assert_eq!(recomputed.len(), n * p);
            for (a, b) in stored.iter().zip(&recomputed) {
                assert!((a - b).abs() < 1e-10, "psi mismatch: {a} vs {b}");
            }

            // The solved duals satisfy the estimating equations, so each
            // column of psi sums to zero over the units.
            for j in 0..p {
                let colsum: f64 = (0..n).map(|i| recomputed[j * n + i]).sum();
                assert!(colsum.abs() < 1e-9, "psi column {j} sums to {colsum}");
            }

            let jac = result.jac.as_ref().expect("exact problem has jac");
            let eps = 1e-6;
            for col in 0..p {
                let mut up = duals.clone();
                let mut down = duals.clone();
                up[col] += eps;
                down[col] -= eps;
                let psi_up = eval_psi_continuous(&inputs, &up, 1.0).unwrap();
                let psi_down = eval_psi_continuous(&inputs, &down, 1.0).unwrap();
                for row in 0..p {
                    let colsum_up: f64 = (0..n).map(|i| psi_up[row * n + i]).sum();
                    let colsum_down: f64 = (0..n).map(|i| psi_down[row * n + i]).sum();
                    let fd = (colsum_up - colsum_down) / (2.0 * eps);
                    let analytic = jac[col * p + row];
                    // The Jacobian entries are of order the sampling-weight
                    // total, so the bound sits well below the derivative scale
                    // while staying orders above central-difference noise.
                    assert!(
                        (fd - analytic).abs() < 1e-6,
                        "jac[{row},{col}]: fd {fd} analytic {analytic}"
                    );
                }
            }
        }
    }

    // Re-evaluating the weights at the solved duals reproduces the solve's own
    // weight vector, and a central finite difference of that vector reproduces
    // the analytic weight derivative, under unit and non-unit measures.
    #[test]
    fn eval_weights_continuous_matches_solve_and_weight_jacobian() {
        let (covs, targets, dist_ind, n, p) = continuous_design();
        let tols = vec![0.0; p];
        for (s, base) in [
            (vec![1.0; n], vec![1.0; n]),
            (continuous_sampling_weights(), continuous_base_weights()),
        ] {
            let inputs = continuous_inputs(&covs, &targets, &base, &s, n, p, &tols);
            let result = solve_continuous(&inputs, &dist_ind, &no_interrupt());
            assert!(result.converged);

            let duals = &result.duals;
            let dw = result.dw_dbeta.as_ref().expect("exact problem has dw");
            let recomputed = eval_weights_continuous(&inputs, duals, 1.0).unwrap();
            assert_eq!(recomputed.len(), n);
            for (i, (a, b)) in result.weights.iter().zip(&recomputed).enumerate() {
                assert!(
                    (a - b).abs() < 1e-12,
                    "weight[{i}]: solve {a} versus eval {b}"
                );
            }

            let eps = 1e-6;
            for col in 0..p {
                let mut up = duals.clone();
                let mut down = duals.clone();
                up[col] += eps;
                down[col] -= eps;
                let w_up = eval_weights_continuous(&inputs, &up, 1.0).unwrap();
                let w_down = eval_weights_continuous(&inputs, &down, 1.0).unwrap();
                for i in 0..n {
                    let fd = (w_up[i] - w_down[i]) / (2.0 * eps);
                    let analytic = dw[col * n + i];
                    assert!(
                        (fd - analytic).abs() < 1e-6,
                        "dw[{i},{col}]: fd {fd} analytic {analytic}"
                    );
                }
            }
        }
    }

    // The reporting scale is the single group's renormalization, free of the
    // duals, so it multiplies every returned weight and every psi entry, the
    // same factor scale_estimating_output applies to the stored blocks.
    #[test]
    fn eval_continuous_applies_the_reporting_scale() {
        let (covs, targets, dist_ind, n, p) = continuous_design();
        let base = vec![1.0; n];
        let s = vec![1.0; n];
        let tols = vec![0.0; p];
        let inputs = continuous_inputs(&covs, &targets, &base, &s, n, p, &tols);
        let result = solve_continuous(&inputs, &dist_ind, &no_interrupt());
        assert!(result.converged);

        let duals = &result.duals;
        let scale = 2.5;
        let unit_weights = eval_weights_continuous(&inputs, duals, 1.0).unwrap();
        let scaled_weights = eval_weights_continuous(&inputs, duals, scale).unwrap();
        for i in 0..n {
            let want = unit_weights[i] * scale;
            assert!(
                (scaled_weights[i] - want).abs() < 1e-12,
                "weight[{i}]: scaled {} expected {want}",
                scaled_weights[i]
            );
        }

        let unit_psi = eval_psi_continuous(&inputs, duals, 1.0).unwrap();
        let scaled_psi = eval_psi_continuous(&inputs, duals, scale).unwrap();
        for (k, (a, b)) in unit_psi.iter().zip(&scaled_psi).enumerate() {
            assert!(
                (b - a * scale).abs() < 1e-12,
                "psi[{k}]: scaled {b} expected {}",
                a * scale
            );
        }
    }

    // A continuous solve carries one dual per constraint. A vector of any other
    // length is reported rather than left to read past the end of the tilt's
    // parameter block.
    #[test]
    fn eval_continuous_rejects_a_mismatched_dual_vector() {
        let (covs, targets, _dist_ind, n, p) = continuous_design();
        let base = vec![1.0; n];
        let s = vec![1.0; n];
        let tols = vec![0.0; p];
        let inputs = continuous_inputs(&covs, &targets, &base, &s, n, p, &tols);
        let coefs = vec![0.0; p - 1];

        let err = eval_psi_continuous(&inputs, &coefs, 1.0).unwrap_err();
        assert!(err.contains("length"), "message was: {err}");
        let err = eval_weights_continuous(&inputs, &coefs, 1.0).unwrap_err();
        assert!(err.contains("length"), "message was: {err}");
    }

    // A constraint matrix or a weight vector shorter than the dimensions it is
    // read at is reported rather than indexed past its end.
    #[test]
    fn eval_continuous_rejects_inputs_shorter_than_their_dimensions() {
        let (covs, targets, _dist_ind, n, p) = continuous_design();
        let base = vec![1.0; n];
        let s = vec![1.0; n];
        let tols = vec![0.0; p];
        let coefs = vec![0.0; p];

        let short_covs = continuous_inputs(&covs[..n * p - 1], &targets, &base, &s, n, p, &tols);
        let err = eval_psi_continuous(&short_covs, &coefs, 1.0).unwrap_err();
        assert!(err.contains("covs"), "message was: {err}");
        let err = eval_weights_continuous(&short_covs, &coefs, 1.0).unwrap_err();
        assert!(err.contains("covs"), "message was: {err}");

        let short_base = continuous_inputs(&covs, &targets, &base[..n - 1], &s, n, p, &tols);
        let err = eval_weights_continuous(&short_base, &coefs, 1.0).unwrap_err();
        assert!(err.contains("weight"), "message was: {err}");

        let short_s = continuous_inputs(&covs, &targets, &base, &s[..n - 1], n, p, &tols);
        let err = eval_psi_continuous(&short_s, &coefs, 1.0).unwrap_err();
        assert!(err.contains("weight"), "message was: {err}");
    }
}
