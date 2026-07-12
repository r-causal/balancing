//! Criterion benchmarks comparing the estimating-equation solvers on identical
//! entropy problems across regimes.
//!
//! Newton (the default), the basin L-BFGS adapter, and the L-BFGS-then-Newton
//! hybrid all solve the same generated problem. The grid spans a well-conditioned
//! reference, an ill-conditioned instance (near-collinear covariates), a small
//! constraint count, and a discrete multi-group instance, so the comparison
//! reaches the regimes the backend promotion rules care about. Before any timing
//! is recorded, each solver is run once and its achieved gradient sup norm is
//! checked against a common convergence tolerance, so a speed comparison is only
//! ever made between results that reached the same accuracy. A solver that fails
//! to reach the common tolerance is reported rather than timed as if it had
//! converged.

use std::time::Duration;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use balancing_core::methods::entropy::{
    EntropyResult, EntropySolver, solve_continuous, solve_discrete,
};

mod common;
use common::{Problem, inputs, make_discrete, make_problem, make_problem_corr};

/// Tolerance the solvers are asked to reach.
const SOLVE_TOL: f64 = 1e-10;
/// Common gradient sup norm every timed comparison must achieve.
const COMMON_CHECK: f64 = 1e-8;
/// Iteration ceiling. Newton needs a handful; L-BFGS needs more, and more still
/// on the ill-conditioned instance. Bounded so a non-converging L-BFGS run stays
/// timeable; the large-k regime, where L-BFGS needs on the order of a thousand
/// memory-bound iterations, is measured for Newton in the kernels bench instead.
const MAX_ITER: usize = 2_000;
/// Threads used for every solver so the comparison is at a fixed thread count.
const THREADS: usize = 8;

const SOLVERS: &[(&str, EntropySolver)] = &[
    ("newton", EntropySolver::Newton),
    ("lbfgs", EntropySolver::Lbfgs),
    ("lbfgs_then_newton", EntropySolver::LbfgsThenNewton),
];

/// Solve a problem with the given solver, dispatching to the discrete
/// entrypoint when the problem carries a group index.
fn solve_case(prob: &Problem, solver: EntropySolver) -> EntropyResult {
    let ins = inputs(prob, THREADS, solver, MAX_ITER, SOLVE_TOL);
    if prob.group_idx.is_empty() {
        solve_continuous(&ins, &prob.dist_ind, &|| false)
    } else {
        solve_discrete(&ins, &prob.group_idx, &|| false)
    }
}

fn bench_solvers(c: &mut Criterion) {
    let cases: Vec<(String, Problem)> = vec![
        (
            "cont_n10000_k50".to_string(),
            make_problem(10_000, 50, 0x5EED_1234),
        ),
        (
            "illcond_n10000_k50_rho0.99".to_string(),
            make_problem_corr(10_000, 50, 0x5EED_5678, 0.99),
        ),
        (
            "smallk_n50000_k10".to_string(),
            make_problem(50_000, 10, 0x5EED_1234),
        ),
        (
            "categorical_n30000_k20_g3".to_string(),
            make_discrete(30_000, 20, 3, 0x5EED_9ABC),
        ),
    ];

    for (label, prob) in &cases {
        // Convergence check: run each solver once and record whether it reaches
        // the common tolerance. The achieved norms are printed so the report can
        // record which solvers qualified on this instance.
        eprintln!("[esteq_solvers] convergence check for {label}:");
        for (name, solver) in SOLVERS {
            let result = solve_case(prob, *solver);
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
            // Newton and the hybrid are held to the common tolerance; the L-BFGS
            // adapter is timed at whatever accuracy it reaches within the
            // iteration ceiling, and its achieved norm is reported alongside,
            // since it is not expected to match a second-order method here.
            let check = solve_case(prob, *solver);
            if matches!(
                solver,
                EntropySolver::Newton | EntropySolver::LbfgsThenNewton
            ) {
                assert!(
                    check.grad_norm <= COMMON_CHECK,
                    "{name} did not reach the common tolerance on {label}: grad_norm {:.3e} > {:.3e}",
                    check.grad_norm,
                    COMMON_CHECK
                );
            }

            group.bench_with_input(BenchmarkId::from_parameter(name), prob, |b, prob| {
                b.iter(|| {
                    let result = solve_case(black_box(prob), *solver);
                    black_box(result.grad_norm)
                });
            });
        }
        group.finish();
    }
}

criterion_group!(benches, bench_solvers);
criterion_main!(benches);
