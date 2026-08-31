//! Criterion benchmarks for the covariate distance transforms.
//!
//! `dist::transform` is the standardizing pass that every energy and
//! covariate-function-of-distance fit runs before the pairwise matrix is
//! assembled. It is `O(n p)` in front of `O(n^2 p)` work, so it should never be
//! the visible cost, and the two benchmarks that touch it indirectly cannot say
//! so: `kernels` measures the kernel build the transform feeds, and
//! `qp_backends` reaches it through a full energy solve where it is a rounding
//! error. This file measures the pass on its own, which is what makes the price
//! of the second moment scan and of centering readable rather than inferred.
//!
//! Grid: ten covariates at n = 20000 and n = 50000, the two sizes where the
//! transform is large enough to time cleanly and which bracket the range an
//! energy fit can reach. Each size is measured on well-scaled covariates and on
//! the same covariates shifted to date-time magnitude. The offset copy is the
//! case the two-pass variance and the centering exist for, and it costs the same
//! arithmetic, so measuring both confirms the accuracy is not bought with a
//! branch that only the ordinary case avoids.

use std::time::Duration;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use balancing_core::dist::Distance;
use balancing_core::dist::transform::transform;

mod common;
use common::make_problem;

/// The magnitude `as.numeric()` gives a POSIXct, added to a well-scaled column
/// to reproduce a date-time covariate: a mean roughly half a million standard
/// deviations from zero.
const DATE_TIME_OFFSET: f64 = 1.7e9;

/// One covariate-matrix size to measure. The transform is linear in both, so `n`
/// dominates at the widths a balancing fit uses.
struct Point {
    n: usize,
    p: usize,
}

// Ordered fast to slow. Each benchmark's estimate is saved as it finishes, so a
// run stopped early still yields the smaller size.
const POINTS: &[Point] = &[Point { n: 20_000, p: 10 }, Point { n: 50_000, p: 10 }];

fn bench_transform(c: &mut Criterion) {
    let mut group = c.benchmark_group("dist_transform");
    // The pass is a few milliseconds at these sizes, so a modest sample gives a
    // stable median without a long run.
    group.sample_size(20);
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(3));

    for point in POINTS {
        let prob = make_problem(point.n, point.p, 0x7A05_5EED);
        let offset: Vec<f64> = prob.covs.iter().map(|&x| DATE_TIME_OFFSET + x).collect();
        let scalings: [(&str, &[f64]); 2] = [("scaled", &prob.covs), ("offset", &offset)];

        for (scaling, covs) in scalings {
            for (name, distance) in [
                ("scaled_euclidean", Distance::ScaledEuclidean),
                ("mahalanobis", Distance::Mahalanobis),
            ] {
                let id = format!("{name}_{scaling}_n{}_p{}", point.n, point.p);
                group.bench_with_input(BenchmarkId::from_parameter(&id), &covs, |b, covs| {
                    b.iter(|| {
                        let out =
                            transform(black_box(covs), point.n, point.p, distance, &prob.s, 1);
                        black_box(out.len())
                    });
                });
            }
        }
    }
    group.finish();
}

criterion_group!(benches, bench_transform);
criterion_main!(benches);
