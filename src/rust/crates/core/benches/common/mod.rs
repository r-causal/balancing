//! Shared synthetic problem generator for the criterion benchmarks.
//!
//! The continuous generators build an exact entropy balancing problem with a
//! known solution: covariates are drawn from a fixed seed, a small target dual
//! vector `beta_star` is chosen, and the constraint targets are set to the
//! covariate means under the exponential tilt at `beta_star`. Solving from a
//! zero start therefore has to recover `beta_star`, so every solver does genuine
//! work and the achieved gradient can be checked against a common tolerance
//! before any timing counts. A correlation knob makes the covariance (and so the
//! Hessian) ill conditioned. A separate multi-group generator builds a discrete
//! problem whose groups are covariate dependent.
//!
//! This module lives under `benches/common/` so cargo does not treat it as a
//! benchmark target of its own.
//!
//! # Seed contract
//!
//! Regression comparisons are only meaningful when successive runs solve the
//! byte-identical problem, so the generated data is pinned in three ways and none
//! of them may change without re-baselining:
//!
//! - The pseudo-random source is [`SplitMix64`], reproduced inline so no external
//!   crate version can shift the stream. Its constants and the Box-Muller normal
//!   transform are fixed.
//! - The draw order is fixed: the shared factor is drawn before the columns, each
//!   column fills unit-major, and the group loadings and assignments follow in the
//!   order written here. Reordering any draw changes every downstream value.
//! - Every caller passes an explicit `seed`; a generator is never seeded from the
//!   clock or the environment. The seeds the benches use are literals in the
//!   benchmark files, so a given benchmark id always maps to one workload.
//!
//! A change to any of these produces a different problem at the same size, which
//! silently invalidates the criterion baselines. Bump the seeds or the baselines
//! deliberately rather than as a side effect of editing this module.

#![allow(dead_code)]

use balancing_core::methods::entropy::{EntropyInputs, EntropySolver};

/// A generated entropy problem and, for the continuous case, its known dual.
pub struct Problem {
    /// Column-major `n` by `p` constraint matrix.
    pub covs: Vec<f64>,
    /// Constraint targets.
    pub targets: Vec<f64>,
    /// Base weights, all one.
    pub base: Vec<f64>,
    /// Sampling weights, all one.
    pub s: Vec<f64>,
    /// Per-constraint tolerances, all zero (the exact problem).
    pub tols: Vec<f64>,
    /// Marginal-distribution constraint indicators, all zero: no column is a
    /// held-exact marginal, so the continuous entrypoint's L1 penalty logic is a
    /// no-op and the exact-problem hot path is measured.
    pub dist_ind: Vec<i32>,
    /// Group index per unit for the discrete entrypoint, empty for a continuous
    /// problem.
    pub group_idx: Vec<i32>,
    /// Number of units.
    pub n: usize,
    /// Number of constraints.
    pub p: usize,
    /// The dual vector the continuous solver must recover, empty for discrete.
    pub beta_star: Vec<f64>,
}

/// A small deterministic PRNG (SplitMix64) so benches need no external crate and
/// produce identical data on every run.
struct SplitMix64(u64);

impl SplitMix64 {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in the open interval (0, 1).
    fn unit(&mut self) -> f64 {
        // 53-bit mantissa, shifted off zero so the log in Box-Muller is finite.
        let bits = self.next_u64() >> 11;
        (bits as f64 + 0.5) * (1.0 / (1u64 << 53) as f64)
    }

    /// A standard normal draw by the Box-Muller transform.
    fn normal(&mut self) -> f64 {
        let u1 = self.unit();
        let u2 = self.unit();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

/// Build an exact continuous entropy problem of `n` units and `p` constraints.
///
/// `rho` in [0, 1) is the pairwise correlation induced across covariate columns
/// through a shared factor. `rho = 0` gives independent columns; a high `rho`
/// makes the columns near collinear, so the weighted covariance that forms the
/// Newton Hessian is ill conditioned.
pub fn make_problem_corr(n: usize, p: usize, seed: u64, rho: f64) -> Problem {
    let mut rng = SplitMix64(seed);
    let mut covs = vec![0.0; n * p];
    let a = rho.max(0.0).sqrt();
    let b = (1.0 - rho.max(0.0)).sqrt();
    // A shared factor per unit couples the columns; each column adds its own
    // independent part.
    let mut factor = vec![0.0; n];
    for f in factor.iter_mut() {
        *f = rng.normal();
    }
    for j in 0..p {
        for i in 0..n {
            covs[j * n + i] = a * factor[i] + b * rng.normal();
        }
    }

    // A small, structured dual so the tilt stays well conditioned even at large
    // p: |C . beta_star| has standard deviation near sqrt(p) * 0.03, moderate.
    let beta_star: Vec<f64> = (0..p).map(|j| 0.03 * (((j % 7) as f64) - 3.0)).collect();

    let (targets, _z) = tilted_means(&covs, n, p, &beta_star);

    Problem {
        covs,
        targets,
        base: vec![1.0; n],
        s: vec![1.0; n],
        tols: vec![0.0; p],
        dist_ind: vec![0; p],
        group_idx: Vec::new(),
        n,
        p,
        beta_star,
    }
}

/// Build a well-conditioned continuous problem (independent columns).
pub fn make_problem(n: usize, p: usize, seed: u64) -> Problem {
    make_problem_corr(n, p, seed, 0.0)
}

/// Build a discrete (multi-group) entropy problem with `n_groups` covariate
/// dependent groups over `n` units and `p` constraints.
///
/// Each unit is assigned to the group whose random loading best matches its
/// covariates, so the groups start imbalanced. Every group is reweighted to the
/// overall covariate means, the average-treatment-effect target.
pub fn make_discrete(n: usize, p: usize, n_groups: usize, seed: u64) -> Problem {
    let mut rng = SplitMix64(seed);
    let mut covs = vec![0.0; n * p];
    for value in covs.iter_mut() {
        *value = rng.normal();
    }

    // Random group loadings; a unit joins the group with the largest score. The
    // loadings are modest so the groups are imbalanced but not near-separable,
    // which keeps the per-group duals finite and the problem realistic.
    let loadings: Vec<f64> = (0..n_groups * p).map(|_| 0.2 * rng.normal()).collect();
    let mut group_idx = vec![0i32; n];
    for i in 0..n {
        let mut best_g = 0usize;
        let mut best_score = f64::NEG_INFINITY;
        for g in 0..n_groups {
            let mut score = 0.0;
            for j in 0..p {
                score += covs[j * n + i] * loadings[g * p + j];
            }
            // A small Gumbel-like perturbation keeps groups from degenerating.
            score += 0.3 * (-(-rng.unit().ln()).ln());
            if score > best_score {
                best_score = score;
                best_g = g;
            }
        }
        group_idx[i] = best_g as i32;
    }

    // Overall covariate means are the shared targets.
    let mut targets = vec![0.0; p];
    for (j, target) in targets.iter_mut().enumerate() {
        let mut acc = 0.0;
        for i in 0..n {
            acc += covs[j * n + i];
        }
        *target = acc / n as f64;
    }

    Problem {
        covs,
        targets,
        base: vec![1.0; n],
        s: vec![1.0; n],
        tols: vec![0.0; p],
        dist_ind: vec![0; p],
        group_idx,
        n,
        p,
        beta_star: Vec::new(),
    }
}

/// Weighted covariate means under the exponential tilt at `beta`, with the
/// normalizing constant.
fn tilted_means(covs: &[f64], n: usize, p: usize, beta: &[f64]) -> (Vec<f64>, f64) {
    let mut e = vec![0.0; n];
    let mut z = 0.0;
    for i in 0..n {
        let mut lin = 0.0;
        for (j, &b) in beta.iter().enumerate() {
            lin += covs[j * n + i] * b;
        }
        let ei = (-lin).exp();
        e[i] = ei;
        z += ei;
    }
    let mut means = vec![0.0; p];
    for (j, m) in means.iter_mut().enumerate() {
        let mut acc = 0.0;
        for i in 0..n {
            acc += e[i] * covs[j * n + i];
        }
        *m = acc / z;
    }
    (means, z)
}

/// Borrow a [`Problem`] as [`EntropyInputs`] for the given solver and threading.
pub fn inputs<'a>(
    prob: &'a Problem,
    threads: usize,
    solver: EntropySolver,
    max_iter: usize,
    tol: f64,
) -> EntropyInputs<'a> {
    EntropyInputs {
        covs: &prob.covs,
        n: prob.n,
        p: prob.p,
        targets: &prob.targets,
        tols: &prob.tols,
        base: &prob.base,
        s: &prob.s,
        n_eff: prob.n as f64,
        threads,
        max_iter,
        tol,
        solver,
    }
}
