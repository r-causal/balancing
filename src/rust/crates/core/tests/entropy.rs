//! Fixture tests for entropy balancing against closed-form and invariant cases.

use balancing_core::methods::entropy::{EntropyInputs, EntropySolver, solve_discrete};

fn no_interrupt() -> impl Fn() -> bool {
    || false
}

// The argument list mirrors the wide EntropyInputs struct this helper fills.
#[allow(clippy::too_many_arguments)]
fn inputs<'a>(
    covs: &'a [f64],
    n: usize,
    p: usize,
    targets: &'a [f64],
    tols: &'a [f64],
    base: &'a [f64],
    s: &'a [f64],
    n_eff: f64,
    threads: usize,
) -> EntropyInputs<'a> {
    EntropyInputs {
        covs,
        n,
        p,
        targets,
        tols,
        base,
        s,
        n_eff,
        threads,
        max_iter: 200,
        tol: 1e-12,
        solver: EntropySolver::Newton,
    }
}

/// Build a standardized single-group problem with a nonzero dual, the shape the
/// R layer produces before crossing the boundary. A small LCG keeps the data
/// reproducible without a dependency.
fn standardized_problem(n: usize, p: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut covs = vec![0.0; n * p];
    let mut state: u64 = 0x2545F4914F6CDD1D;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) as f64) / ((1u64 << 53) as f64)
    };
    for j in 0..p {
        for i in 0..n {
            covs[j * n + i] = next();
        }
        let col = &mut covs[j * n..(j + 1) * n];
        let mean: f64 = col.iter().sum::<f64>() / n as f64;
        let var: f64 = col.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64;
        let sd = var.sqrt();
        for x in col.iter_mut() {
            *x = (*x - mean) / sd;
        }
    }
    let targets: Vec<f64> = (0..p).map(|j| 0.05 - 0.02 * j as f64).collect();
    let tols = vec![0.0; p];
    (covs, targets, tols)
}

/// The L-BFGS-then-Newton hybrid must reach a machine-precision solution: the
/// mandatory Newton polish drives the gradient far below what the L-BFGS warm
/// start alone reaches, so the estimating-equation output is evaluated at the
/// same accuracy Newton provides.
#[test]
fn hybrid_polish_reaches_machine_precision() {
    let n = 8_000;
    let p = 3;
    let (covs, targets, tols) = standardized_problem(n, p);
    let base = vec![1.0; n];
    let s = vec![1.0; n];
    let group_idx = vec![0_i32; n];

    let solve = |solver: EntropySolver| {
        let mut inp = inputs(&covs, n, p, &targets, &tols, &base, &s, n as f64, 1);
        inp.solver = solver;
        solve_discrete(&inp, &group_idx, &no_interrupt())
    };

    let lbfgs = solve(EntropySolver::Lbfgs);
    let hybrid = solve(EntropySolver::LbfgsThenNewton);

    assert_eq!(hybrid.solver, "lbfgs_then_newton");
    assert!(hybrid.converged);
    assert!(hybrid.psi.is_some());
    // The polish tightens the gradient well below the warm start's residual.
    assert!(
        hybrid.grad_norm <= 1e-10,
        "hybrid grad_norm {} not machine precision",
        hybrid.grad_norm
    );
    assert!(
        hybrid.grad_norm <= lbfgs.grad_norm,
        "polish did not improve on the warm start: hybrid {} vs lbfgs {}",
        hybrid.grad_norm,
        lbfgs.grad_norm
    );
}

/// An interrupt that becomes pending exactly at the hybrid handoff, after the
/// L-BFGS warm start's own interrupt check has passed, must still be surfaced:
/// the Newton polish never runs, so the estimate is only at the loose warm
/// tolerance and the solve must report `interrupted` and not `converged` rather
/// than inheriting the warm start's converged status.
#[test]
fn hybrid_handoff_interrupt_is_surfaced() {
    let n = 4_000;
    let p = 2;
    let (covs, targets, tols) = standardized_problem(n, p);
    let base = vec![1.0; n];
    let s = vec![1.0; n];
    let group_idx = vec![0_i32; n];

    // The basin warm start polls the interrupt once, at its end; the handoff
    // check polls it again. Returning `false` on the first probe and `true`
    // thereafter leaves the warm start not interrupted but makes the interrupt
    // pending precisely at the handoff.
    let calls = std::cell::Cell::new(0usize);
    let interrupt = || {
        let c = calls.get() + 1;
        calls.set(c);
        c >= 2
    };

    let mut inp = inputs(&covs, n, p, &targets, &tols, &base, &s, n as f64, 1);
    inp.solver = EntropySolver::LbfgsThenNewton;
    let result = solve_discrete(&inp, &group_idx, &interrupt);

    assert_eq!(result.solver, "lbfgs_then_newton");
    assert!(result.interrupted, "handoff interrupt was not surfaced");
    assert!(
        !result.converged,
        "solve claimed convergence at the loose warm-start point"
    );
}

/// The estimating-equation output, not only the weights, is bit-identical across
/// thread counts: the parallel psi and weight-derivative fill and the
/// deterministic Jacobian reduction all preserve the determinism contract.
#[test]
fn estimating_output_is_deterministic_across_thread_counts() {
    let n = 12_000;
    let p = 3;
    let (covs, targets, tols) = standardized_problem(n, p);
    let base = vec![1.0; n];
    let s = vec![1.0; n];
    let group_idx = vec![0_i32; n];

    let solve_with = |threads: usize| {
        let inp = inputs(&covs, n, p, &targets, &tols, &base, &s, n as f64, threads);
        solve_discrete(&inp, &group_idx, &no_interrupt())
    };

    let reference = solve_with(1);
    assert!(reference.converged);
    let ref_psi = reference.psi.as_ref().expect("exact problem returns psi");
    let ref_jac = reference.jac.as_ref().expect("exact problem returns jac");
    let ref_dw = reference
        .dw_dbeta
        .as_ref()
        .expect("exact problem returns dw_dbeta");

    for threads in [2, 4, 8] {
        let other = solve_with(threads);
        let psi = other.psi.as_ref().unwrap();
        let jac = other.jac.as_ref().unwrap();
        let dw = other.dw_dbeta.as_ref().unwrap();
        for i in 0..ref_psi.len() {
            assert_eq!(
                psi[i].to_bits(),
                ref_psi[i].to_bits(),
                "psi[{i}] differs at {threads} threads"
            );
            assert_eq!(
                dw[i].to_bits(),
                ref_dw[i].to_bits(),
                "dw_dbeta[{i}] differs at {threads} threads"
            );
        }
        for i in 0..ref_jac.len() {
            assert_eq!(
                jac[i].to_bits(),
                ref_jac[i].to_bits(),
                "jac[{i}] differs at {threads} threads"
            );
        }
    }
}

/// A pending interrupt stops the solve and is surfaced: the result reports
/// `interrupted` and is not marked converged, so the R layer can re-signal the
/// interrupt rather than warn about non-convergence.
#[test]
fn a_pending_interrupt_is_surfaced() {
    let covs = [1.0, 1.0, 1.0, 0.0, 0.0];
    let n = 5;
    let targets = [0.7];
    let tols = [0.0];
    let base = [1.0; 5];
    let s = [1.0; 5];
    let group_idx = [0, 0, 0, 0, 0];

    let inp = inputs(&covs, n, 1, &targets, &tols, &base, &s, 5.0, 1);
    let result = solve_discrete(&inp, &group_idx, &|| true);

    assert!(result.interrupted);
    assert!(!result.converged);
}

/// Entropy ATT with one binary covariate has a two-cell closed form: every
/// control with `x = 1` receives `n_eff * t / a` and every control with `x = 0`
/// receives `n_eff * (1 - t) / b`, where `t` is the treated share, `a` and `b`
/// are the control cell counts.
#[test]
fn att_binary_covariate_matches_closed_form() {
    // Three controls with x = 1, two with x = 0.
    let covs = [1.0, 1.0, 1.0, 0.0, 0.0];
    let n = 5;
    let a = 3.0;
    let b = 2.0;
    let t = 0.7;
    let targets = [t];
    let tols = [0.0];
    let base = [1.0; 5];
    let s = [1.0; 5];
    let n_eff = 5.0;
    let group_idx = [0, 0, 0, 0, 0];

    let inp = inputs(&covs, n, 1, &targets, &tols, &base, &s, n_eff, 1);
    let result = solve_discrete(&inp, &group_idx, &no_interrupt());

    assert!(result.converged);
    let w_x1 = n_eff * t / a;
    let w_x0 = n_eff * (1.0 - t) / b;
    for i in 0..3 {
        assert!(
            (result.weights[i] - w_x1).abs() < 1e-9,
            "weight[{i}] = {}",
            result.weights[i]
        );
    }
    for i in 3..5 {
        assert!(
            (result.weights[i] - w_x0).abs() < 1e-9,
            "weight[{i}] = {}",
            result.weights[i]
        );
    }

    // Weighted mean reproduces the target and total effective size is n_eff.
    let total: f64 = result.weights.iter().sum();
    assert!((total - n_eff).abs() < 1e-9);
    let weighted_mean: f64 = result
        .weights
        .iter()
        .zip(covs.iter())
        .map(|(w, x)| w * x)
        .sum::<f64>()
        / total;
    assert!((weighted_mean - t).abs() < 1e-9);

    // The exact problem returns estimating-equation output.
    assert!(result.psi.is_some());
    assert!(result.jac.is_some());
    assert_eq!(result.solver, "newton");
}

/// When the target already equals the base-weighted mean the dual solution is
/// zero and the weights collapse to the base weights (scaled so their sum is
/// `n_eff`). Here `n_eff` is the base-weight total, so weights equal base
/// weights exactly.
#[test]
fn weights_equal_base_weights_when_constraints_already_hold() {
    let covs = [0.2, 0.5, 0.8, 0.1, 0.9];
    let n = 5;
    let base = [1.0, 2.0, 1.0, 1.0, 2.0];
    let s = [1.0; 5];
    let base_total: f64 = base.iter().sum();
    let weighted_mean: f64 = covs
        .iter()
        .zip(base.iter())
        .map(|(x, q)| x * q)
        .sum::<f64>()
        / base_total;
    let targets = [weighted_mean];
    let tols = [0.0];
    let group_idx = [0, 0, 0, 0, 0];

    let inp = inputs(&covs, n, 1, &targets, &tols, &base, &s, base_total, 1);
    let result = solve_discrete(&inp, &group_idx, &no_interrupt());

    assert!(result.converged);
    for (i, (&weight, &b)) in result.weights.iter().zip(&base).enumerate() {
        assert!(
            (weight - b).abs() < 1e-9,
            "weight[{i}] = {weight} expected {b}"
        );
    }
    // The dual variable is zero.
    assert!(result.duals[0].abs() < 1e-8);
}

/// The same problem solved with different thread counts yields bit-identical
/// weights, the determinism guarantee of the reduction.
#[test]
fn weights_are_deterministic_across_thread_counts() {
    let n = 20_000;
    let p = 2;
    // Standardized covariates (weighted mean 0, sd 1 per column), which is what
    // the R layer produces before crossing the boundary. A small LCG keeps the
    // data reproducible without a dependency.
    let mut covs = vec![0.0; n * p];
    let mut state: u64 = 0x2545F4914F6CDD1D;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) as f64) / ((1u64 << 53) as f64)
    };
    for j in 0..p {
        for i in 0..n {
            covs[j * n + i] = next();
        }
        let col = &mut covs[j * n..(j + 1) * n];
        let mean: f64 = col.iter().sum::<f64>() / n as f64;
        let var: f64 = col.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64;
        let sd = var.sqrt();
        for x in col.iter_mut() {
            *x = (*x - mean) / sd;
        }
    }
    let base = vec![1.0; n];
    let s = vec![1.0; n];
    let group_idx = vec![0_i32; n];

    // Targets offset from the standardized mean of zero, so the dual is nonzero.
    let targets = [0.05, -0.03];
    let tols = [0.0, 0.0];

    let solve_with = |threads: usize| {
        let inp = inputs(&covs, n, p, &targets, &tols, &base, &s, n as f64, threads);
        solve_discrete(&inp, &group_idx, &no_interrupt())
    };

    let reference = solve_with(1);
    assert!(reference.converged);
    for threads in [2, 4, 8] {
        let other = solve_with(threads);
        for i in 0..n {
            assert_eq!(
                reference.weights[i].to_bits(),
                other.weights[i].to_bits(),
                "weight[{i}] differs at {threads} threads"
            );
        }
        assert_eq!(reference.duals[0].to_bits(), other.duals[0].to_bits());
        assert_eq!(reference.duals[1].to_bits(), other.duals[1].to_bits());
    }
}

/// A positive tolerance selects the inexact problem: FISTA runs, the weights
/// still satisfy the moment within the tolerance band, and no estimating
/// equations are returned.
#[test]
fn positive_tolerance_uses_the_inexact_problem() {
    let covs = [1.0, 1.0, 1.0, 0.0, 0.0];
    let n = 5;
    let targets = [0.7];
    let tols = [0.05];
    let base = [1.0; 5];
    let s = [1.0; 5];
    let group_idx = [0, 0, 0, 0, 0];

    let inp = inputs(&covs, n, 1, &targets, &tols, &base, &s, 5.0, 1);
    let result = solve_discrete(&inp, &group_idx, &no_interrupt());

    assert_eq!(result.solver, "fista");
    assert!(result.psi.is_none());
    assert!(result.jac.is_none());

    let total: f64 = result.weights.iter().sum();
    let weighted_mean: f64 = result
        .weights
        .iter()
        .zip(covs.iter())
        .map(|(w, x)| w * x)
        .sum::<f64>()
        / total;
    // The moment falls within the tolerance band around the target.
    assert!(
        (weighted_mean - 0.7).abs() <= 0.05 + 1e-6,
        "weighted mean {weighted_mean} outside tolerance"
    );
}
