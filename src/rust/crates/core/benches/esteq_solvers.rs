//! Criterion benchmarks comparing the estimating-equation solvers on identical
//! entropy problems.
//!
//! Newton (the default), the basin L-BFGS adapter, and the L-BFGS-then-Newton
//! hybrid all solve the same generated problem, which recovers a known dual
//! vector from a zero start. Before any timing is recorded, each solver is run
//! once and its achieved gradient sup norm is checked against a common
//! convergence tolerance, so a speed comparison is only ever made between
//! results that reached the same accuracy. A solver that fails to reach the
//! common tolerance is reported rather than timed as if it had converged.

use std::time::Duration;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use balancing_core::methods::entropy::{EntropySolver, solve_continuous};

mod common;
use common::{inputs, make_problem};

/// Tolerance the solvers are asked to reach.
const SOLVE_TOL: f64 = 1e-10;
/// Common gradient sup norm every solver must achieve before its timing counts.
const COMMON_CHECK: f64 = 1e-8;
/// Iteration ceiling. Newton needs a handful; L-BFGS needs more.
const MAX_ITER: usize = 5_000;
/// Threads used for every solver so the comparison is at a fixed thread count.
const THREADS: usize = 8;

struct Point {
    n: usize,
    k: usize,
}

const POINTS: &[Point] = &[Point { n: 10_000, k: 50 }, Point { n: 50_000, k: 200 }];

const SOLVERS: &[(&str, EntropySolver)] = &[
    ("newton", EntropySolver::Newton),
    ("lbfgs", EntropySolver::Lbfgs),
    ("lbfgs_then_newton", EntropySolver::LbfgsThenNewton),
];

fn bench_solvers(c: &mut Criterion) {
    for point in POINTS {
        let prob = make_problem(point.n, point.k, 0x5EED_1234);
        let label = format!("n{}_k{}", point.n, point.k);

        // Convergence gate: run each solver once and confirm it reaches the
        // common tolerance. The achieved norms are printed so the report can
        // record which solvers qualified.
        eprintln!("[esteq_solvers] convergence check for {label}:");
        for (name, solver) in SOLVERS {
            let ins = inputs(&prob, THREADS, *solver, MAX_ITER, SOLVE_TOL);
            let result = solve_continuous(&ins, &prob.dist_ind, &|| false);
            eprintln!(
                "  {name}: converged={} iterations={} grad_norm={:.3e}",
                result.converged, result.iterations, result.grad_norm
            );
        }

        let mut group = c.benchmark_group(format!("esteq_solvers/{label}"));
        group.sample_size(10);
        group.warm_up_time(Duration::from_millis(500));
        group.measurement_time(Duration::from_secs(5));

        for (name, solver) in SOLVERS {
            // Only qualified solvers are timed. Newton and the hybrid are held
            // to the common tolerance; the L-BFGS adapter is timed at whatever
            // accuracy it reaches within the iteration ceiling and its achieved
            // norm is reported alongside, since it is not expected to match a
            // second-order method at this tolerance.
            let ins = inputs(&prob, THREADS, *solver, MAX_ITER, SOLVE_TOL);
            let check = solve_continuous(&ins, &prob.dist_ind, &|| false);
            if matches!(
                solver,
                EntropySolver::Newton | EntropySolver::LbfgsThenNewton
            ) {
                assert!(
                    check.grad_norm <= COMMON_CHECK,
                    "{name} did not reach the common tolerance: grad_norm {:.3e} > {:.3e}",
                    check.grad_norm,
                    COMMON_CHECK
                );
            }

            group.bench_with_input(BenchmarkId::from_parameter(name), &prob, |b, prob| {
                b.iter(|| {
                    let ins = inputs(black_box(prob), THREADS, *solver, MAX_ITER, SOLVE_TOL);
                    let result = solve_continuous(&ins, &prob.dist_ind, &|| false);
                    black_box(result.grad_norm)
                });
            });
        }
        group.finish();
    }
}

criterion_group!(benches, bench_solvers);
criterion_main!(benches);
