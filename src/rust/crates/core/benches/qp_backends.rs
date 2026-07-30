//! Criterion benchmarks comparing the quadratic-program backends.
//!
//! The shared backend trait lets the same spec run through OSQP and Clarabel so
//! their wall times are directly comparable. Two spec families are measured. The
//! positive-semidefinite simplex projection is solved by both backends, so it is
//! the like-for-like timing comparison on a problem Clarabel accepts. The energy
//! balancing family is driven through its public solve entry points; its
//! quadratic form is indefinite by construction, so it is an OSQP-only timing
//! path and Clarabel rejects it up front. The indefinite rejection, and the
//! numerical failure that justifies it, are recorded once outside the timing
//! loop rather than timed.

use std::hint::black_box;
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};

use balancing_core::dist::{Distance, pairwise};
use balancing_core::methods::energy::{
    EnergyContInputs, EnergyDiscreteInputs, EnergyEstimand, solve_cont, solve_discrete,
};
#[cfg(feature = "qp-clarabel")]
use balancing_core::qp::clarabel::Clarabel;
use balancing_core::qp::osqp::Osqp;
use balancing_core::qp::{Convexity, PMat, QpBackend, QpOptions, QpSpec};
#[cfg(feature = "qp-clarabel")]
use balancing_core::qp::{QpError, objective};

/// A positive-semidefinite simplex-projection spec in `n` variables: minimize
/// `0.5 sum x_i^2` on the unit simplex with a non-negativity box. Both backends
/// accept it, so it isolates their speed on a well-posed convex problem.
fn simplex_spec(n: usize) -> QpSpec {
    // Box rows 0..n (identity) plus a final sum row over every variable.
    let mut indptr = Vec::with_capacity(n + 1);
    let mut indices = Vec::new();
    let mut values = Vec::new();
    indptr.push(0);
    for j in 0..n {
        indices.push(j);
        values.push(1.0);
        indices.push(n);
        values.push(1.0);
        indptr.push(indices.len());
    }
    let mut l = vec![0.0; n];
    let mut u = vec![f64::INFINITY; n];
    l.push(1.0);
    u.push(1.0);
    QpSpec {
        n,
        m: n + 1,
        p: PMat::Diagonal(vec![1.0; n]),
        q: vec![0.0; n],
        a_indptr: indptr,
        a_indices: indices,
        a_values: values,
        l,
        u,
        convexity: Convexity::Psd,
    }
}

/// A deterministic confounded covariate matrix in `n` units and `p` columns,
/// column-major, with a covariate-dependent grouping into `n_levels`.
fn confounded(n: usize, p: usize, n_levels: usize) -> (Vec<f64>, Vec<i32>) {
    let mut covs = vec![0.0; n * p];
    let mut levels = vec![0i32; n];
    for i in 0..n {
        let t = i as f64;
        for j in 0..p {
            covs[j * n + i] = (t * (j as f64 + 1.0) * 0.017).sin() + 0.3 * (j as f64);
        }
        // A covariate-dependent score split into `n_levels` bands so the groups
        // differ in their covariate distribution.
        let score = covs[i] + 0.2 * covs[n + i];
        let band = ((score + 1.0) * 0.5 * n_levels as f64).floor();
        levels[i] = (band.max(0.0) as usize).min(n_levels - 1) as i32;
    }
    (covs, levels)
}

/// Borrow a confounded problem as discrete energy inputs under an estimand, with
/// optional first-moment constraints on the covariates themselves.
// The argument list mirrors the wide EnergyDiscreteInputs struct this helper fills.
#[allow(clippy::too_many_arguments)]
fn discrete_inputs<'a>(
    covs: &'a [f64],
    levels: &'a [i32],
    n: usize,
    p: usize,
    n_levels: usize,
    estimand: EnergyEstimand,
    s: &'a [f64],
    moment_covs: &'a [f64],
    targets: &'a [f64],
    tols: &'a [f64],
    n_moments: usize,
) -> EnergyDiscreteInputs<'a> {
    EnergyDiscreteInputs {
        covs,
        n,
        p,
        distance: Distance::ScaledEuclidean,
        levels,
        n_levels,
        estimand,
        s,
        min_weight: 1e-8,
        lambda: 1e-4,
        moment_covs,
        n_moments,
        targets,
        tols,
        threads: 1,
        qp: QpOptions::default(),
    }
}

/// Standardized first-moment targets and tolerance bands for the covariate
/// columns of a confounded problem, so a moment-constrained cell holds each
/// group's covariate means near the pooled sample mean.
fn moment_targets(covs: &[f64], n: usize, p: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut std = vec![0.0; n * p];
    let mut targets = vec![0.0; p];
    for j in 0..p {
        let col = &covs[j * n..(j + 1) * n];
        let mean = col.iter().sum::<f64>() / n as f64;
        let var = col.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0);
        let sd = if var > 0.0 { var.sqrt() } else { 1.0 };
        for i in 0..n {
            std[j * n + i] = (col[i] - mean) / sd;
        }
        targets[j] = std[j * n..(j + 1) * n].iter().sum::<f64>() / n as f64;
    }
    // A modest tolerance band, so the moment rows are active but feasible.
    let tols = vec![0.05; p];
    (std, targets, tols)
}

fn bench_backends(c: &mut Criterion) {
    record_backend_contract();

    let opts = QpOptions::default();
    let mut group = c.benchmark_group("qp_backends");
    group.measurement_time(Duration::from_secs(3));

    // Convex simplex projection: the like-for-like comparison on a spec both
    // backends accept.
    for &n in &[50usize, 200, 1000] {
        let spec = simplex_spec(n);
        group.bench_with_input(BenchmarkId::new("osqp_simplex", n), &n, |b, _| {
            b.iter(|| black_box(Osqp.solve(&spec, &opts, &|| false).unwrap()));
        });
        #[cfg(feature = "qp-clarabel")]
        group.bench_with_input(BenchmarkId::new("clarabel_simplex", n), &n, |b, _| {
            b.iter(|| black_box(Clarabel.solve(&spec, &opts, &|| false).unwrap()));
        });
    }
    group.finish();

    // Energy family, OSQP-only, driven through the public solve entry points.
    // Sizes stay in the range where the dense n by n solve is tractable for
    // repeated sampling; the larger end-to-end sizes are measured on the R side.
    let mut energy = c.benchmark_group("energy_osqp");
    energy.sample_size(10);
    energy.measurement_time(Duration::from_secs(4));
    let p = 4;

    for &n in &[250usize, 500, 1000] {
        let (covs, levels) = confounded(n, p, 2);
        let s = vec![1.0; n];
        let empty: Vec<f64> = Vec::new();

        let ate_impr = discrete_inputs(
            &covs,
            &levels,
            n,
            p,
            2,
            EnergyEstimand::Ate { improved: true },
            &s,
            &empty,
            &empty,
            &empty,
            0,
        );
        energy.bench_with_input(BenchmarkId::new("binary_ate_improved", n), &n, |b, _| {
            b.iter(|| black_box(solve_discrete(&ate_impr, &|| false).objective));
        });

        let ate_plain = discrete_inputs(
            &covs,
            &levels,
            n,
            p,
            2,
            EnergyEstimand::Ate { improved: false },
            &s,
            &empty,
            &empty,
            &empty,
            0,
        );
        energy.bench_with_input(BenchmarkId::new("binary_ate_plain", n), &n, |b, _| {
            b.iter(|| black_box(solve_discrete(&ate_plain, &|| false).objective));
        });

        let att = discrete_inputs(
            &covs,
            &levels,
            n,
            p,
            2,
            EnergyEstimand::Focal { focal: 1 },
            &s,
            &empty,
            &empty,
            &empty,
            0,
        );
        energy.bench_with_input(BenchmarkId::new("binary_att", n), &n, |b, _| {
            b.iter(|| black_box(solve_discrete(&att, &|| false).objective));
        });

        let (mcovs, mtargets, mtols) = moment_targets(&covs, n, p);
        let ate_moments = discrete_inputs(
            &covs,
            &levels,
            n,
            p,
            2,
            EnergyEstimand::Ate { improved: true },
            &s,
            &mcovs,
            &mtargets,
            &mtols,
            p,
        );
        energy.bench_with_input(BenchmarkId::new("binary_ate_moments", n), &n, |b, _| {
            b.iter(|| black_box(solve_discrete(&ate_moments, &|| false).objective));
        });
    }

    // Categorical average treatment effect at a single representative size.
    {
        let n = 1000usize;
        let (covs, levels) = confounded(n, p, 3);
        let s = vec![1.0; n];
        let empty: Vec<f64> = Vec::new();
        let cat = discrete_inputs(
            &covs,
            &levels,
            n,
            p,
            3,
            EnergyEstimand::Ate { improved: true },
            &s,
            &empty,
            &empty,
            &empty,
            0,
        );
        energy.bench_with_input(BenchmarkId::new("categorical_ate", n), &n, |b, _| {
            b.iter(|| black_box(solve_discrete(&cat, &|| false).objective));
        });
    }

    // Continuous average treatment effect: the distance-covariance form.
    {
        let n = 500usize;
        let (covs, _levels) = confounded(n, p, 2);
        let treat: Vec<f64> = (0..n)
            .map(|i| covs[i] + 0.3 * covs[n + i] + 0.1 * (i as f64 * 0.05).cos())
            .collect();
        let s = vec![1.0; n];
        let cont = EnergyContInputs {
            covs: &covs,
            n,
            p,
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
        energy.bench_with_input(BenchmarkId::new("continuous_ate", n), &n, |b, _| {
            b.iter(|| black_box(solve_cont(&cont, &|| false).objective));
        });
    }

    energy.finish();
}

/// An indefinite dense quadratic form of the energy kind: the negated pairwise
/// distance of a small confounded sample, doubled, on a simplex box. The energy
/// method tags exactly this shape `Indefinite`.
fn indefinite_energy_like_spec() -> QpSpec {
    let n = 40usize;
    let p = 3usize;
    let (covs, _levels) = confounded(n, p, 2);
    let dist = pairwise::euclidean(&covs, n, p, 1);
    // P = -2 * dist (already doubled). The negated distance matrix is indefinite.
    let pmat: Vec<f64> = dist.iter().map(|d| -2.0 * d).collect();

    // Box rows 0..n with a min-weight lower bound, then a sum row fixing the
    // total to n, so the indefinite objective is bounded below on the feasible
    // set exactly as the energy simplex constraints make it.
    let mut indptr = Vec::with_capacity(n + 1);
    let mut indices = Vec::new();
    let mut values = Vec::new();
    indptr.push(0);
    for j in 0..n {
        indices.push(j);
        values.push(1.0);
        indices.push(n);
        values.push(1.0);
        indptr.push(indices.len());
    }
    let mut l = vec![1e-8; n];
    let mut u = vec![f64::INFINITY; n];
    l.push(n as f64);
    u.push(n as f64);
    QpSpec {
        n,
        m: n + 1,
        p: PMat::Dense(pmat),
        q: vec![0.0; n],
        a_indptr: indptr,
        a_indices: indices,
        a_values: values,
        l,
        u,
        convexity: Convexity::Indefinite,
    }
}

/// Record the backend contract once, outside the timing loop: OSQP solves the
/// indefinite energy form, and where the Clarabel backend is compiled in, its
/// half of the contract is recorded alongside.
fn record_backend_contract() {
    let opts = QpOptions::default();
    let spec = indefinite_energy_like_spec();
    let mut unverified = spec.clone();
    unverified.convexity = Convexity::Unverified;

    #[cfg(feature = "qp-clarabel")]
    record_clarabel_contract(&opts, &spec, &unverified);

    let osqp_sol = Osqp
        .solve(&unverified, &opts, &|| false)
        .expect("osqp sets up the indefinite form");
    eprintln!(
        "qp_backends: osqp on the same form -> status {:?}, obj {:.6}",
        osqp_sol.status, osqp_sol.obj
    );
}

/// The Clarabel half of the backend contract: it rejects the indefinite energy
/// form up front on its convexity tag, the same form left unverified confirms
/// the rejection is not merely conservative because the interior-point solve
/// does not reach an optimal status, and the two backends agree on the convex
/// simplex projection their timings are compared on. That agreement is asserted
/// so the comparison is only ever drawn between solved, agreeing results.
#[cfg(feature = "qp-clarabel")]
fn record_clarabel_contract(opts: &QpOptions, spec: &QpSpec, unverified: &QpSpec) {
    match Clarabel.solve(spec, opts, &|| false) {
        Err(QpError::Indefinite) => {
            eprintln!("qp_backends: clarabel rejects the indefinite energy form (structured)");
        }
        other => panic!("expected an indefinite rejection from clarabel, got {other:?}"),
    }

    match Clarabel.solve(unverified, opts, &|| false) {
        Ok(sol) => eprintln!(
            "qp_backends: clarabel on the same form (unverified) -> status {:?}, obj {:.6}",
            sol.status, sol.obj
        ),
        Err(e) => {
            eprintln!("qp_backends: clarabel on the same form (unverified) -> setup error {e}")
        }
    }

    let psd = simplex_spec(200);
    let a = Osqp.solve(&psd, opts, &|| false).unwrap();
    let b = Clarabel.solve(&psd, opts, &|| false).unwrap();
    let oa = objective(&psd, &a.x);
    let ob = objective(&psd, &b.x);
    assert!(
        (oa - ob).abs() <= 1e-6 * oa.abs().max(1.0),
        "simplex objective parity: osqp {oa} vs clarabel {ob}"
    );
    eprintln!("qp_backends: simplex objective parity osqp {oa:.8} clarabel {ob:.8} (within 1e-6)");
}

criterion_group!(benches, bench_backends);
criterion_main!(benches);
