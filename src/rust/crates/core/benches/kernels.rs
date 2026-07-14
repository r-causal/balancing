//! Criterion benchmarks for the entropy dual kernel.
//!
//! The per-unit fused gradient and Hessian accumulation is the hot path of the
//! entropy solve. That accumulation is a private implementation detail, so the
//! benchmark drives it through the public [`solve_continuous`] entry point,
//! whose per-iteration cost is dominated by that kernel. The generated problems
//! recover a known dual vector from a zero start, so the iteration count is
//! stable across runs and the per-pass kernel cost is recoverable from the
//! reported solve time and iteration count.
//!
//! Grid: the corners of the design grid (smallest and largest n by k) plus the
//! n = 50000, k = 200 reference point across thread counts. The full grid is too
//! long for a single run; these points bound the range and expose thread
//! scaling at the reference size.

use std::time::Duration;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use balancing_core::dist::kernels::{Kernel, KernelParams, MaternNu, build_kernel};
use balancing_core::methods::entropy::{EntropySolver, solve_continuous};

mod common;
use common::{inputs, make_problem};

/// Solver convergence tolerance on the gradient sup norm.
const GRAD_TOL: f64 = 1e-10;
/// Iteration ceiling; the exact Newton solve reaches the tolerance well below it.
const MAX_ITER: usize = 200;

/// One (n, k, threads) point to measure.
struct Point {
    n: usize,
    k: usize,
    threads: usize,
}

// Ordered fast to slow. Each benchmark's estimate is saved as it finishes, so a
// run stopped early still yields every completed point; the single-threaded
// reference, the most expensive point, is last.
const POINTS: &[Point] = &[
    Point {
        n: 1_000,
        k: 10,
        threads: 1,
    },
    Point {
        n: 1_000,
        k: 200,
        threads: 1,
    },
    Point {
        n: 10_000,
        k: 50,
        threads: 1,
    },
    Point {
        n: 10_000,
        k: 50,
        threads: 2,
    },
    Point {
        n: 10_000,
        k: 50,
        threads: 8,
    },
    Point {
        n: 50_000,
        k: 10,
        threads: 8,
    },
    Point {
        n: 50_000,
        k: 200,
        threads: 8,
    },
    Point {
        n: 50_000,
        k: 200,
        threads: 2,
    },
    Point {
        n: 50_000,
        k: 200,
        threads: 1,
    },
];

fn bench_entropy_solve(c: &mut Criterion) {
    let mut group = c.benchmark_group("entropy_solve");
    // These solves range from sub-millisecond to several seconds; a small sample
    // keeps the large points tractable while still giving a stable median.
    group.sample_size(10);
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(5));

    for point in POINTS {
        let prob = make_problem(point.n, point.k, 0x5EED_1234);
        let id = format!("n{}_k{}_t{}", point.n, point.k, point.threads);
        group.bench_with_input(BenchmarkId::from_parameter(&id), &prob, |b, prob| {
            b.iter(|| {
                let ins = inputs(
                    black_box(prob),
                    point.threads,
                    EntropySolver::Newton,
                    MAX_ITER,
                    GRAD_TOL,
                );
                let result = solve_continuous(&ins, &prob.dist_ind, &|| false);
                black_box(result.grad_norm)
            });
        });
    }
    group.finish();
}

/// One kernel-matrix size to measure. The kernel build is an `n` by `n`
/// assembly, so `n` dominates; `p` is held at a representative width.
struct KernelPoint {
    n: usize,
    p: usize,
}

// The kernel matrix is quadratic in n, so the two sizes bracket the range while
// keeping a single run modest; the larger point is last so an early stop still
// yields the smaller one.
const KERNEL_POINTS: &[KernelPoint] = &[
    KernelPoint { n: 1_000, p: 10 },
    KernelPoint { n: 10_000, p: 10 },
];

fn bench_kernel_matrix(c: &mut Criterion) {
    let mut group = c.benchmark_group("kernel_matrix");
    // The n = 10000 build is an order of magnitude larger than n = 1000, so a
    // small sample keeps the larger point tractable.
    group.sample_size(10);
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(5));

    for point in KERNEL_POINTS {
        let prob = make_problem(point.n, point.p, 0x9A17_5EED);
        let id = format!("gaussian_n{}_p{}", point.n, point.p);
        group.bench_with_input(BenchmarkId::from_parameter(&id), &prob, |b, prob| {
            b.iter(|| {
                let params = KernelParams {
                    kernel: Kernel::Gaussian,
                    bw_scale: 1.0,
                    matern_nu: MaternNu::ThreeHalves,
                    t_proj: &[],
                    n_draws: 0,
                };
                let k = build_kernel(
                    black_box(&prob.covs),
                    prob.n,
                    prob.p,
                    &prob.s,
                    &params,
                    &[],
                    1,
                );
                black_box(k.len())
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_entropy_solve, bench_kernel_matrix);
criterion_main!(benches);
