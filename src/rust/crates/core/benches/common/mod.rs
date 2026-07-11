//! Shared synthetic problem generator for the criterion benchmarks.
//!
//! The generator builds an exact entropy balancing problem with a known
//! solution: covariates are drawn from a fixed seed, a small target dual vector
//! `beta_star` is chosen, and the constraint targets are set to the covariate
//! means under the exponential tilt at `beta_star`. Solving from a zero start
//! therefore has to recover `beta_star`, so every solver does genuine work and
//! the achieved gradient can be checked against a common tolerance before any
//! timing counts.
//!
//! This module lives under `benches/common/` so cargo does not treat it as a
//! benchmark target of its own.

#![allow(dead_code)]

use balancing_core::methods::entropy::{EntropyInputs, EntropySolver};

/// A generated exact entropy problem and its known dual solution.
pub struct Problem {
    /// Column-major `n` by `p` constraint matrix.
    pub covs: Vec<f64>,
    /// Constraint targets: the tilted covariate means at `beta_star`.
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
    /// Number of units.
    pub n: usize,
    /// Number of constraints.
    pub p: usize,
    /// The dual vector the solver must recover.
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

/// Build an exact entropy problem of `n` units and `p` constraints from `seed`.
pub fn make_problem(n: usize, p: usize, seed: u64) -> Problem {
    let mut rng = SplitMix64(seed);
    let mut covs = vec![0.0; n * p];
    // Column-major fill.
    for value in covs.iter_mut() {
        *value = rng.normal();
    }

    // A small, structured dual so the tilt stays well conditioned even at large
    // p: |C . beta_star| has standard deviation near sqrt(p) * 0.03, moderate.
    let beta_star: Vec<f64> = (0..p).map(|j| 0.03 * (((j % 7) as f64) - 3.0)).collect();

    // Tilted weights e_i = exp(-C_i . beta_star), then targets_j = weighted mean.
    let mut e = vec![0.0; n];
    let mut z = 0.0;
    for i in 0..n {
        let mut lin = 0.0;
        for (j, &b) in beta_star.iter().enumerate() {
            lin += covs[j * n + i] * b;
        }
        let ei = (-lin).exp();
        e[i] = ei;
        z += ei;
    }
    let mut targets = vec![0.0; p];
    for (j, target) in targets.iter_mut().enumerate() {
        let mut acc = 0.0;
        for i in 0..n {
            acc += e[i] * covs[j * n + i];
        }
        *target = acc / z;
    }

    Problem {
        covs,
        targets,
        base: vec![1.0; n],
        s: vec![1.0; n],
        tols: vec![0.0; p],
        dist_ind: vec![0; p],
        n,
        p,
        beta_star,
    }
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
