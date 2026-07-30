//! Criterion benchmark comparing the quadratic-program backends on the stable
//! balancing weights family.
//!
//! Stable balancing weights minimize a strictly convex diagonal quadratic form
//! over a simplex-type constraint set, so the quadratic term is positive
//! semidefinite by construction and both backends accept every spec. This is the
//! first shipped method whose real instances give a like-for-like OSQP versus
//! Clarabel comparison, unlike energy balancing, whose indefinite form Clarabel
//! rejects structurally.
//!
//! The spec assembly here mirrors `methods::sbw::solve_discrete` and
//! `solve_cont` exactly: the same group-normalized sampling weights, the same box
//! with a minimum-weight floor, the same group-sum rows, and the same moment or
//! correlation rows. Building the spec in the bench, rather than calling the
//! solver, is what lets the identical `QpSpec` cross into both backends through
//! the shared trait. The parity contract runs once outside the timing loop and
//! records, per instance, that both backends converged, that their objectives
//! agree to a relative 1e-6, and the worst constraint violation of each, so a
//! speed comparison is only ever drawn between solved, agreeing results.

use std::hint::black_box;
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};

use balancing_core::methods::qp_balance::{ConstraintBuilder, ZERO_SW, group_normalized};
use balancing_core::methods::sbw::SbwEstimand;
#[cfg(feature = "qp-clarabel")]
use balancing_core::qp::clarabel::Clarabel;
use balancing_core::qp::osqp::Osqp;
use balancing_core::qp::{Convexity, PMat, QpBackend, QpOptions, QpSpec};
#[cfg(feature = "qp-clarabel")]
use balancing_core::qp::{QpSolution, objective};

// -- Deterministic data generation ------------------------------------------

/// SplitMix64, so the bench needs no external RNG crate and produces identical
/// data on every run.
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

/// Standardize each column of a column-major `n` by `p` matrix in place to mean
/// zero and unit sample standard deviation, the scale the R layer hands the
/// solver for numeric covariates, so a moment tolerance reads as a standardized
/// mean difference.
fn standardize(covs: &mut [f64], n: usize, p: usize) {
    for j in 0..p {
        let col = &mut covs[j * n..(j + 1) * n];
        let mean = col.iter().sum::<f64>() / n as f64;
        let var = col.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0);
        let sd = if var > 0.0 { var.sqrt() } else { 1.0 };
        for v in col.iter_mut() {
            *v = (*v - mean) / sd;
        }
    }
}

/// The pooled sample mean of each standardized column, the average-treatment
/// target the group means are pulled toward.
fn pooled_targets(covs: &[f64], n: usize, p: usize) -> Vec<f64> {
    (0..p)
        .map(|j| covs[j * n..(j + 1) * n].iter().sum::<f64>() / n as f64)
        .collect()
}

/// The mean of each standardized column over the units at `level`, the focal
/// target a focal estimand pulls the other groups toward.
fn group_targets(covs: &[f64], levels: &[i32], n: usize, p: usize, level: i32) -> Vec<f64> {
    (0..p)
        .map(|j| {
            let mut acc = 0.0;
            let mut cnt = 0usize;
            for i in 0..n {
                if levels[i] == level {
                    acc += covs[j * n + i];
                    cnt += 1;
                }
            }
            if cnt > 0 { acc / cnt as f64 } else { 0.0 }
        })
        .collect()
}

/// A binary confounded instance: `p` standard-normal covariates and a binary
/// exposure whose propensity depends on a few of them, so the two groups start
/// imbalanced but overlap in covariate space and the moment balance is feasible.
/// Returns the standardized covariates and the zero-or-one level of each unit.
fn binary_data(n: usize, p: usize, seed: u64) -> (Vec<f64>, Vec<i32>) {
    let mut rng = SplitMix64(seed);
    let mut covs = vec![0.0; n * p];
    for v in covs.iter_mut() {
        *v = rng.normal();
    }
    let coef = [0.4, -0.3, 0.2, -0.2];
    let mut levels = vec![0i32; n];
    for i in 0..n {
        let mut lp = 0.0;
        for j in 0..p {
            lp += coef[j % coef.len()] * covs[j * n + i];
        }
        let pr = 1.0 / (1.0 + (-lp).exp());
        levels[i] = if rng.unit() < pr { 1 } else { 0 };
    }
    standardize(&mut covs, n, p);
    (covs, levels)
}

/// A categorical confounded instance in `n_levels` groups: the exposure level is
/// the argmax of a random linear score plus noise, so the groups differ in their
/// covariate distribution.
fn categorical_data(n: usize, p: usize, n_levels: usize, seed: u64) -> (Vec<f64>, Vec<i32>) {
    let mut rng = SplitMix64(seed);
    let mut covs = vec![0.0; n * p];
    for v in covs.iter_mut() {
        *v = rng.normal();
    }
    let loadings: Vec<f64> = (0..n_levels * p).map(|_| 0.5 * rng.normal()).collect();
    let mut levels = vec![0i32; n];
    for i in 0..n {
        let mut best_g = 0usize;
        let mut best = f64::NEG_INFINITY;
        for g in 0..n_levels {
            let mut score = 0.0;
            for j in 0..p {
                score += covs[j * n + i] * loadings[g * p + j];
            }
            score += 0.4 * (-(-rng.unit().ln()).ln());
            if score > best {
                best = score;
                best_g = g;
            }
        }
        levels[i] = best_g as i32;
    }
    standardize(&mut covs, n, p);
    (covs, levels)
}

/// A continuous confounded instance: `p` covariates and an exposure correlated
/// with them, so uniform weights leave a nonzero weighted exposure-covariate
/// correlation the constraint rows must pull inside the band. Covariates are
/// standardized; the exposure is returned raw, as the continuous solver
/// standardizes it internally.
fn continuous_data(n: usize, p: usize, seed: u64) -> (Vec<f64>, Vec<f64>) {
    let mut rng = SplitMix64(seed);
    let mut covs = vec![0.0; n * p];
    for v in covs.iter_mut() {
        *v = rng.normal();
    }
    let treat: Vec<f64> = (0..n)
        .map(|i| {
            let mut t = 0.6 * covs[i] + 0.3 * covs[n + i];
            t += rng.normal();
            t
        })
        .collect();
    standardize(&mut covs, n, p);
    (covs, treat)
}

// -- Spec assembly, mirroring methods::sbw ----------------------------------

/// The active variables and constrained group levels for a discrete estimand,
/// reproducing `methods::sbw::active_layout`: every present unit for the average
/// treatment effect, only the non-focal units for a focal estimand.
fn active_layout(
    n: usize,
    levels: &[i32],
    n_levels: usize,
    estimand: SbwEstimand,
) -> (Vec<usize>, Vec<usize>) {
    match estimand {
        SbwEstimand::Ate => (
            (0..n).filter(|&i| levels[i] >= 0).collect(),
            (0..n_levels).collect(),
        ),
        SbwEstimand::Focal { focal } => (
            (0..n)
                .filter(|&i| levels[i] >= 0 && levels[i] as usize != focal)
                .collect(),
            (0..n_levels).filter(|&t| t != focal).collect(),
        ),
    }
}

/// Assemble the discrete stable balancing spec exactly as `solve_discrete` does:
/// a doubled identity objective, a box with the minimum-weight floor and pinned
/// zero-sampling-weight units, one group-sum row per reweighted group, and one
/// moment row per group per covariate holding the group mean within the band.
#[allow(clippy::too_many_arguments)]
fn discrete_spec(
    covs: &[f64],
    levels: &[i32],
    n: usize,
    p: usize,
    n_levels: usize,
    estimand: SbwEstimand,
    s: &[f64],
    min_weight: f64,
    targets: &[f64],
    tols: &[f64],
) -> QpSpec {
    let s_norm = group_normalized(s, levels, n_levels);
    let mut n_t = vec![0usize; n_levels];
    for &g in levels {
        if g >= 0 {
            n_t[g as usize] += 1;
        }
    }
    let swnt: Vec<f64> = levels
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

    let (active, group_levels) = active_layout(n, levels, n_levels, estimand);
    let nvar = active.len();
    let p_mat = PMat::Diagonal(vec![2.0; nvar]);
    let q = vec![0.0; nvar];

    let mut builder = ConstraintBuilder::new(nvar);
    let pinned: Vec<bool> = active.iter().map(|&i| s[i].abs() < ZERO_SW).collect();
    builder.add_box(min_weight, &pinned);

    for &t in &group_levels {
        let coeffs: Vec<f64> = active
            .iter()
            .map(|&i| {
                if levels[i] as usize == t {
                    swnt[i]
                } else {
                    0.0
                }
            })
            .collect();
        builder.add_dense_row(&coeffs, 1.0, 1.0);
    }

    for &t in &group_levels {
        for c in 0..p {
            let coeffs: Vec<f64> = active
                .iter()
                .map(|&i| {
                    if levels[i] as usize == t {
                        covs[c * n + i] * swnt[i]
                    } else {
                        0.0
                    }
                })
                .collect();
            let band = tols[c].abs();
            builder.add_dense_row(&coeffs, targets[c] - band, targets[c] + band);
        }
    }

    let (m, indptr, indices, values, l, u) = builder.finish();
    QpSpec {
        n: nvar,
        m,
        p: p_mat,
        q,
        a_indptr: indptr,
        a_indices: indices,
        a_values: values,
        l,
        u,
        convexity: Convexity::Psd,
    }
}

/// Reliability-weighted variance, matching `methods::sbw::weighted_variance`, so
/// the standardized exposure the correlation rows use agrees with the solver.
fn weighted_variance(x: &[f64], w: &[f64]) -> f64 {
    let (mut sw, mut sw2, mut swx, mut swxx) = (0.0, 0.0, 0.0, 0.0);
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

/// Assemble the continuous stable balancing spec exactly as `solve_cont` does: a
/// doubled identity objective, a box, a single total-sum row fixing the weighted
/// total to `n`, and one bounded weighted-correlation row per covariate.
fn continuous_spec(
    covs: &[f64],
    treat: &[f64],
    n: usize,
    p: usize,
    s: &[f64],
    min_weight: f64,
    tols: &[f64],
) -> QpSpec {
    let nf = n as f64;
    let s_sum: f64 = s.iter().sum();
    let s_scaled: Vec<f64> = if s_sum > 0.0 {
        s.iter().map(|&si| si * nf / s_sum).collect()
    } else {
        vec![1.0; n]
    };
    let sw: f64 = s.iter().sum();
    let a_mean =
        s.iter().zip(treat).map(|(&wi, &ti)| wi * ti).sum::<f64>() / sw.max(f64::MIN_POSITIVE);
    let a_var = weighted_variance(treat, s);
    let a_sd = if a_var > 0.0 { a_var.sqrt() } else { 1.0 };
    let treat_std: Vec<f64> = treat.iter().map(|&t| (t - a_mean) / a_sd).collect();

    let p_mat = PMat::Diagonal(vec![2.0; n]);
    let q = vec![0.0; n];

    let mut builder = ConstraintBuilder::new(n);
    let pinned: Vec<bool> = s.iter().map(|&si| si.abs() < ZERO_SW).collect();
    builder.add_box(min_weight, &pinned);
    builder.add_dense_row(&s_scaled, nf, nf);

    let denom = (nf - 1.0).max(1.0);
    for c in 0..p {
        let coeffs: Vec<f64> = (0..n)
            .map(|i| covs[c * n + i] * treat_std[i] * s_scaled[i] / denom)
            .collect();
        let tol = tols[c].abs();
        builder.add_dense_row(&coeffs, -tol, tol);
    }

    let (m, indptr, indices, values, l, u) = builder.finish();
    QpSpec {
        n,
        m,
        p: p_mat,
        q,
        a_indptr: indptr,
        a_indices: indices,
        a_values: values,
        l,
        u,
        convexity: Convexity::Psd,
    }
}

// -- The promotion-rule grid -------------------------------------------------

/// One named quadratic-program instance both backends solve.
struct Instance {
    name: String,
    spec: QpSpec,
}

/// Build the stable balancing promotion dataset: binary average-treatment and
/// focal instances across the size and tolerance grid, a categorical
/// average-treatment instance, two continuous instances, and two stress
/// instances (a near-infeasible tight-tolerance binary problem and an ill-scaled
/// covariate problem).
fn build_grid() -> Vec<Instance> {
    let p = 4usize;
    let min_w = 1e-8;
    let mut out = Vec::new();

    let sizes = [500usize, 2000, 5000, 20000];
    let tols_grid = [0.01f64, 0.05, 0.1];

    for &n in &sizes {
        let (covs, levels) = binary_data(n, p, 0x5B_0000 + n as u64);
        for &tol in &tols_grid {
            let tols = vec![tol; p];
            // Average treatment effect: every group pulled to the pooled mean.
            let ate_targets = pooled_targets(&covs, n, p);
            let ate = discrete_spec(
                &covs,
                &levels,
                n,
                p,
                2,
                SbwEstimand::Ate,
                &vec![1.0; n],
                min_w,
                &ate_targets,
                &tols,
            );
            out.push(Instance {
                name: format!("binary_ate_n{n}_tol{tol}"),
                spec: ate,
            });

            // Focal (average treatment effect on the treated): level one held at
            // unit weight, level zero pulled to the level-one mean.
            let att_targets = group_targets(&covs, &levels, n, p, 1);
            let att = discrete_spec(
                &covs,
                &levels,
                n,
                p,
                2,
                SbwEstimand::Focal { focal: 1 },
                &vec![1.0; n],
                min_w,
                &att_targets,
                &tols,
            );
            out.push(Instance {
                name: format!("binary_att_n{n}_tol{tol}"),
                spec: att,
            });
        }
    }

    // Categorical average treatment effect, three levels, one representative size.
    {
        let n = 1000usize;
        let (covs, levels) = categorical_data(n, p, 3, 0x5B_CA71);
        for &tol in &tols_grid {
            let tols = vec![tol; p];
            let targets = pooled_targets(&covs, n, p);
            let spec = discrete_spec(
                &covs,
                &levels,
                n,
                p,
                3,
                SbwEstimand::Ate,
                &vec![1.0; n],
                min_w,
                &targets,
                &tols,
            );
            out.push(Instance {
                name: format!("categorical_ate_n{n}_tol{tol}"),
                spec,
            });
        }
    }

    // Continuous average treatment effect at two sizes, one representative
    // correlation tolerance.
    for &n in &[500usize, 2000] {
        let (covs, treat) = continuous_data(n, p, 0x5B_C047 + n as u64);
        let tols = vec![0.1; p];
        let spec = continuous_spec(&covs, &treat, n, p, &vec![1.0; n], min_w, &tols);
        out.push(Instance {
            name: format!("continuous_ate_n{n}_tol0.1"),
            spec,
        });
    }

    // Stress one: a near-infeasible tight tolerance. Strong confounding at a
    // 0.008 band leaves the feasible set a thin sliver, so the solve runs near
    // the feasibility boundary. Recorded honestly whatever status results.
    {
        let n = 2000usize;
        let mut rng = SplitMix64(0x5B_57E5);
        let mut covs = vec![0.0; n * p];
        for v in covs.iter_mut() {
            *v = rng.normal();
        }
        let mut levels = vec![0i32; n];
        for i in 0..n {
            // Strong single-covariate confounding so the groups separate sharply.
            let lp = 1.2 * covs[i] + 0.8 * covs[n + i];
            let pr = 1.0 / (1.0 + (-lp).exp());
            levels[i] = if rng.unit() < pr { 1 } else { 0 };
        }
        standardize(&mut covs, n, p);
        let targets = pooled_targets(&covs, n, p);
        let tols = vec![0.008; p];
        let spec = discrete_spec(
            &covs,
            &levels,
            n,
            p,
            2,
            SbwEstimand::Ate,
            &vec![1.0; n],
            min_w,
            &targets,
            &tols,
        );
        out.push(Instance {
            name: "stress_near_infeasible_n2000_tol0.008".to_string(),
            spec,
        });
    }

    // Stress two: ill-scaled covariates. The moment columns span several orders
    // of magnitude before any standardization, so the constraint matrix is badly
    // conditioned and both backends must cope with the scale spread.
    {
        let n = 2000usize;
        let mut rng = SplitMix64(0x5B_15CA);
        let mut covs = vec![0.0; n * p];
        // Column j is scaled by 10^(2j - 3): magnitudes from 1e-3 to 1e3.
        for j in 0..p {
            let scale = 10f64.powi(2 * j as i32 - 3);
            for i in 0..n {
                covs[j * n + i] = rng.normal() * scale;
            }
        }
        let mut levels = vec![0i32; n];
        for i in 0..n {
            // Recover each column's unit-scale draw so the propensity has real
            // signal regardless of the column magnitudes.
            let lp = 0.5 * (covs[i] / 1e-3) + 0.4 * (covs[n + i] / 0.1);
            let pr = 1.0 / (1.0 + (-lp).exp());
            levels[i] = if rng.unit() < pr { 1 } else { 0 };
        }
        // Deliberately do not standardize: the raw scale spread is the stressor.
        let targets = pooled_targets(&covs, n, p);
        // Bands proportional to each column's own spread so the rows stay active.
        let tols: Vec<f64> = (0..p)
            .map(|j| {
                let col = &covs[j * n..(j + 1) * n];
                let mean = col.iter().sum::<f64>() / n as f64;
                let sd =
                    (col.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0)).sqrt();
                0.05 * sd
            })
            .collect();
        let spec = discrete_spec(
            &covs,
            &levels,
            n,
            p,
            2,
            SbwEstimand::Ate,
            &vec![1.0; n],
            min_w,
            &targets,
            &tols,
        );
        out.push(Instance {
            name: "stress_ill_scaled_n2000".to_string(),
            spec,
        });
    }

    out
}

// -- Parity contract ---------------------------------------------------------

/// The worst constraint violation of `x` against the spec bounds: for every row,
/// how far `A x` falls below `l` or above `u`, maxed over rows and clamped at
/// zero. Read only by the parity contract, which compares the two backends.
// The indexed loops walk several parallel arrays (the CSC triplet, then the row
// bounds) where an index is the clearest form.
#[cfg(feature = "qp-clarabel")]
#[allow(clippy::needless_range_loop)]
fn max_violation(spec: &QpSpec, x: &[f64]) -> f64 {
    let mut ax = vec![0.0; spec.m];
    for j in 0..spec.n {
        let start = spec.a_indptr[j];
        let end = spec.a_indptr[j + 1];
        let xj = x[j];
        for k in start..end {
            ax[spec.a_indices[k]] += spec.a_values[k] * xj;
        }
    }
    let mut worst = 0.0f64;
    for r in 0..spec.m {
        let below = spec.l[r] - ax[r];
        let above = ax[r] - spec.u[r];
        worst = worst.max(below).max(above);
    }
    worst
}

/// Solve every grid instance through both backends outside the timing loop and
/// print a structured line per instance: statuses, objectives, their relative
/// gap, worst constraint violations, and iteration counts. The report reads
/// these lines to decide the parity precondition and the promotion inputs. No
/// assertion aborts the bench: a stress instance is allowed to fail, and the
/// failure is data the report records. The contract is a statement about the two
/// backends together, so it exists only where both are compiled in.
#[cfg(feature = "qp-clarabel")]
fn record_sbw_parity(instances: &[Instance]) {
    let opts = QpOptions::default();
    eprintln!("sbw_parity: begin ({} instances)", instances.len());
    for inst in instances {
        let os = Osqp.solve(&inst.spec, &opts, &|| false);
        let cl = Clarabel.solve(&inst.spec, &opts, &|| false);
        let (os_status, os_obj, os_iters, os_viol) = summarize(&inst.spec, &os);
        let (cl_status, cl_obj, cl_iters, cl_viol) = summarize(&inst.spec, &cl);
        let rel_gap = if os_obj.is_finite() && cl_obj.is_finite() {
            (os_obj - cl_obj).abs() / os_obj.abs().max(1.0)
        } else {
            f64::NAN
        };
        eprintln!(
            "sbw_parity: {name} osqp_status={os_status} clarabel_status={cl_status} \
             obj_osqp={os_obj:.10e} obj_clarabel={cl_obj:.10e} rel_gap={rel_gap:.3e} \
             osqp_maxviol={os_viol:.3e} clarabel_maxviol={cl_viol:.3e} \
             osqp_iters={os_iters} clarabel_iters={cl_iters}",
            name = inst.name,
        );
    }
    eprintln!("sbw_parity: end");
}

/// Status name, objective, iterations, and worst violation of a backend result.
#[cfg(feature = "qp-clarabel")]
fn summarize(
    spec: &QpSpec,
    res: &Result<QpSolution, balancing_core::qp::QpError>,
) -> (&'static str, f64, usize, f64) {
    match res {
        Ok(sol) => (
            sol.status.as_str(),
            objective(spec, &sol.x),
            sol.iterations,
            max_violation(spec, &sol.x),
        ),
        Err(_) => ("setup_error", f64::NAN, 0, f64::NAN),
    }
}

// -- Timing ------------------------------------------------------------------

fn bench_sbw(c: &mut Criterion) {
    let instances = build_grid();
    #[cfg(feature = "qp-clarabel")]
    record_sbw_parity(&instances);

    let opts = QpOptions::default();
    let mut group = c.benchmark_group("sbw_backends");
    group.sample_size(10);
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(3));

    for inst in &instances {
        let spec = &inst.spec;
        group.bench_with_input(BenchmarkId::new("osqp", &inst.name), spec, |b, spec| {
            b.iter(|| black_box(Osqp.solve(spec, &opts, &|| false).unwrap()));
        });
        #[cfg(feature = "qp-clarabel")]
        group.bench_with_input(BenchmarkId::new("clarabel", &inst.name), spec, |b, spec| {
            b.iter(|| black_box(Clarabel.solve(spec, &opts, &|| false).unwrap()));
        });
    }
    group.finish();
}

criterion_group!(benches, bench_sbw);
criterion_main!(benches);
