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
                let e = self.s[gi] * self.base[gi] * (-self.lin(gi, beta)).exp();
                acc.z += e;
                for j in 0..p {
                    let cj = self.covs[j * self.n + gi];
                    acc.m[j] += e * cj;
                    if want_smat {
                        for k in 0..p {
                            let ck = self.covs[k * self.n + gi];
                            acc.smat[j * p + k] += e * cj * ck;
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

/// Solve a discrete entropy problem whose units are partitioned into groups by
/// `group_idx` (values `0..G`), each balanced to the shared `targets`.
pub fn solve_discrete(
    inputs: &EntropyInputs<'_>,
    group_idx: &[i32],
    interrupt: &dyn Fn() -> bool,
) -> EntropyResult {
    let n_groups = group_idx.iter().copied().max().map_or(0, |g| g + 1).max(0) as usize;
    let mut groups: Vec<Vec<usize>> = vec![Vec::new(); n_groups];
    for (i, &g) in group_idx.iter().enumerate() {
        if g >= 0 {
            groups[g as usize].push(i);
        }
    }
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
    let all: Vec<usize> = (0..inputs.n).collect();
    let l1: Vec<f64> = inputs
        .tols
        .iter()
        .zip(dist_ind)
        .map(|(&t, &d)| if d != 0 { 0.0 } else { t })
        .collect();
    solve_groups(inputs, std::slice::from_ref(&all), &l1, interrupt)
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
                inputs,
                idx,
                g,
                &beta,
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

/// Fill the estimating-equation blocks for one group into the global arrays.
///
/// For unit `i` in group `g`, with parameter block starting at column `g * p`:
/// `psi_ij = s_i w_i (C_ij - target_j)`, `dw_ik = w_i (mbar_k - C_ik)`, and the
/// Jacobian block is `-sum_i s_i w_i (C_i - target)(C_i - target)'`.
#[allow(clippy::too_many_arguments)]
fn fill_estimating_output(
    inputs: &EntropyInputs<'_>,
    idx: &[usize],
    g: usize,
    _beta: &[f64],
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
    for (local, &gi) in idx.iter().enumerate() {
        let w = group_weights[local];
        let sw = inputs.s[gi] * w;
        for (j, &mbar_j) in mbar.iter().enumerate() {
            let dev = inputs.covs[j * n + gi] - inputs.targets[j];
            let col = offset + j;
            // psi and dw are n by P, column-major.
            psi[col * n + gi] = sw * dev;
            dw[col * n + gi] = w * (mbar_j - inputs.covs[j * n + gi]);
            for k in 0..p {
                let dev_k = inputs.covs[k * n + gi] - inputs.targets[k];
                let row = offset + j;
                let jcol = offset + k;
                jac[jcol * total_params + row] -= sw * dev * dev_k;
            }
        }
    }
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
}
