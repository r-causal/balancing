//! Golden parity tests for entropy balancing.
//!
//! Fixtures live under `tests/golden/` as JSON files written by
//! `data-raw/make-golden.R` from a pinned reference implementation. Each fixture
//! captures the solver inputs at the matrix level (the constraint matrix,
//! targets, tolerances, and weights) so this harness compares solver parity in
//! isolation from covariate processing. When no fixtures are present the harness
//! passes and reports that they are generated in the R API step.
//!
//! Entropy fixture schema (all arrays are plain JSON numbers):
//!
//! ```json
//! {
//!   "kind": "discrete" | "continuous",
//!   "n": <int>, "p": <int>,
//!   "covs": [ column-major n*p ],
//!   "targets": [ p ], "tols": [ p ],
//!   "base": [ n ], "s": [ n ], "n_eff": <number>,
//!   "group_idx": [ n ],            // discrete only
//!   "expected_weights": [ n ],
//!   "rel_tol": <number>
//! }
//! ```
//!
//! The `ipt`, `cbps`, `cbps_multi`, and `cbps_cont` kinds carry the exposure and
//! method options instead of targets; each solver-specific handler below
//! documents its fields. The over-identified `cbps` fixture compares the GMM
//! criterion against `expected_obj` rather than the weights.

use std::path::{Path, PathBuf};

use balancing_core::dist::Distance;
use balancing_core::links::Link;
use balancing_core::methods::cbps::{
    CbpsContInputs, CbpsEstimand, CbpsInputs, CbpsMultiInputs, solve as solve_cbps,
    solve_cont as solve_cbps_cont, solve_multi as solve_cbps_multi,
};
use balancing_core::methods::energy::{
    EnergyContInputs, EnergyDiscreteInputs, EnergyEstimand, EnergyResult,
    solve_cont as solve_energy_cont, solve_discrete as solve_energy_discrete,
};
use balancing_core::methods::entropy::{
    EntropyInputs, EntropySolver, solve_continuous, solve_discrete,
};
use balancing_core::methods::ipt::{IptEstimand, IptInputs, solve as solve_ipt};
use balancing_core::methods::qp_balance::group_normalized;
use balancing_core::methods::sbw::{
    SbwContInputs, SbwDiscreteInputs, SbwEstimand, SbwNorm, SbwResult,
    solve_cont as solve_sbw_cont, solve_discrete as solve_sbw_discrete,
};
use balancing_core::qp::QpOptions;
use serde_json::Value;

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

fn nums(value: &Value, key: &str) -> Vec<f64> {
    value[key]
        .as_array()
        .unwrap_or_else(|| panic!("fixture field `{key}` is not an array"))
        .iter()
        .map(|v| v.as_f64().expect("fixture array holds numbers"))
        .collect()
}

fn scalar(value: &Value, key: &str) -> f64 {
    value[key]
        .as_f64()
        .unwrap_or_else(|| panic!("fixture field `{key}` is not a number"))
}

fn compare_weights(path: &Path, actual: &[f64], expected: &[f64], rel_tol: f64) {
    let scale = expected
        .iter()
        .fold(0.0_f64, |m, w| m.max(w.abs()))
        .max(1.0);
    for (i, &want) in expected.iter().enumerate() {
        let diff = (actual[i] - want).abs();
        assert!(
            diff <= rel_tol * scale,
            "{}: weight[{i}] = {} expected {want} (diff {diff})",
            path.display(),
            actual[i],
        );
    }
}

/// Solve an inverse probability tilting fixture and compare weights.
///
/// Schema adds `treat` (the zero-based level per unit), `link`, `estimand`, and
/// an optional `focal` (defaulting to the treated level `1`) to the shared
/// fields; `targets`, `tols`, `base`, and `n_eff` are absent.
fn check_ipt_fixture(path: &Path, f: &Value) {
    let n = scalar(f, "n") as usize;
    let p = scalar(f, "p") as usize;
    let covs = nums(f, "covs");
    let treat: Vec<i32> = nums(f, "treat").iter().map(|t| *t as i32).collect();
    let s = nums(f, "s");
    let expected = nums(f, "expected_weights");
    let rel_tol = scalar(f, "rel_tol");
    let link = match f["link"].as_str().expect("fixture link is a string") {
        "logit" => Link::Logit,
        "probit" => Link::Probit,
        "cloglog" => Link::Cloglog,
        other => panic!("unknown fixture link `{other}`"),
    };
    let n_levels = treat.iter().copied().max().map_or(0, |m| m + 1).max(0) as usize;
    let focal = f["focal"].as_f64().map(|x| x as usize).unwrap_or(1);
    let estimand = match f["estimand"]
        .as_str()
        .expect("fixture estimand is a string")
    {
        "ate" => IptEstimand::Ate,
        "att" | "atc" => IptEstimand::Focal(focal),
        other => panic!("unknown fixture estimand `{other}`"),
    };

    let inputs = IptInputs {
        covs: &covs,
        n,
        p,
        treat: &treat,
        n_levels,
        s: &s,
        link,
        estimand,
        threads: 1,
        max_iter: 500,
        tol: 1e-12,
    };
    let result = solve_ipt(&inputs, &|| false);
    compare_weights(path, &result.weights, &expected, rel_tol);
}

fn link_from(f: &Value) -> Link {
    match f["link"].as_str().unwrap_or("logit") {
        "logit" => Link::Logit,
        "probit" => Link::Probit,
        "cloglog" => Link::Cloglog,
        other => panic!("unknown fixture link `{other}`"),
    }
}

fn cbps_estimand(name: &str) -> CbpsEstimand {
    match name {
        "ate" => CbpsEstimand::Ate,
        "att" => CbpsEstimand::Att,
        "atc" => CbpsEstimand::Atc,
        "ato" => CbpsEstimand::Ato,
        other => panic!("unknown cbps estimand `{other}`"),
    }
}

/// Solve a binary covariate balancing propensity score fixture.
///
/// The just-identified fixtures compare weights; the over-identified fixtures
/// compare the generalized-method-of-moments criterion, which the design's
/// tolerance policy requires to sit at or below the reference plus `1e-8`.
/// Schema adds `treat` (zero/one), `estimand`, `link`, `over`, and `twostep` to
/// the shared fields, with `expected_obj` for the over-identified form.
fn check_cbps_fixture(path: &Path, f: &Value) {
    let n = scalar(f, "n") as usize;
    let p = scalar(f, "p") as usize;
    let covs = nums(f, "covs");
    let treat: Vec<i32> = nums(f, "treat").iter().map(|t| *t as i32).collect();
    let s = nums(f, "s");
    let over = f["over"].as_bool().unwrap_or(false);
    let twostep = f["twostep"].as_bool().unwrap_or(true);
    let estimand = cbps_estimand(f["estimand"].as_str().expect("cbps estimand is a string"));

    let inputs = CbpsInputs {
        covs_mod: &covs,
        covs_bal: &covs,
        n,
        p_mod: p,
        p_bal: p,
        treat: &treat,
        s: &s,
        link: link_from(f),
        estimand,
        over,
        twostep,
        threads: 1,
        max_iter: 500,
        tol: 1e-12,
    };
    let result = solve_cbps(&inputs, &|| false);
    if over {
        let reference = scalar(f, "expected_obj");
        let ours = result
            .gmm_obj
            .expect("over-identified reports the criterion");
        assert!(
            ours <= reference + 1e-8,
            "{}: gmm objective {ours} exceeds reference {reference} + 1e-8",
            path.display(),
        );
    } else {
        let expected = nums(f, "expected_weights");
        let rel_tol = scalar(f, "rel_tol");
        compare_weights(path, &result.weights, &expected, rel_tol);
    }
}

/// Solve a categorical covariate balancing propensity score fixture. Schema adds
/// `treat` (zero-based level), `estimand`, `link`, and an optional `focal`.
fn check_cbps_multi_fixture(path: &Path, f: &Value) {
    let n = scalar(f, "n") as usize;
    let p = scalar(f, "p") as usize;
    let covs = nums(f, "covs");
    let treat: Vec<i32> = nums(f, "treat").iter().map(|t| *t as i32).collect();
    let s = nums(f, "s");
    let expected = nums(f, "expected_weights");
    let rel_tol = scalar(f, "rel_tol");
    let n_levels = treat.iter().copied().max().map_or(0, |m| m + 1).max(0) as usize;
    let focal = f["focal"].as_f64().map(|x| x as usize).unwrap_or(0);
    let estimand = cbps_estimand(f["estimand"].as_str().expect("cbps estimand is a string"));

    let inputs = CbpsMultiInputs {
        covs: &covs,
        n,
        p,
        treat: &treat,
        n_levels,
        focal,
        s: &s,
        link: link_from(f),
        estimand,
        threads: 1,
        max_iter: 500,
        tol: 1e-12,
    };
    let result = solve_cbps_multi(&inputs, &|| false);
    compare_weights(path, &result.weights, &expected, rel_tol);
}

/// Solve a continuous covariate balancing propensity score fixture. Schema
/// carries `expo` (the exposure) in place of `treat`.
fn check_cbps_cont_fixture(path: &Path, f: &Value) {
    let n = scalar(f, "n") as usize;
    let p = scalar(f, "p") as usize;
    let covs = nums(f, "covs");
    let expo = nums(f, "expo");
    let s = nums(f, "s");
    let expected = nums(f, "expected_weights");
    let rel_tol = scalar(f, "rel_tol");

    let inputs = CbpsContInputs {
        covs: &covs,
        n,
        p,
        expo: &expo,
        s: &s,
        threads: 1,
        max_iter: 500,
        tol: 1e-12,
    };
    let result = solve_cbps_cont(&inputs, &|| false);
    compare_weights(path, &result.weights, &expected, rel_tol);
}

fn distance_from(f: &Value) -> Distance {
    match f["distance"].as_str().unwrap_or("scaled_euclidean") {
        "scaled_euclidean" => Distance::ScaledEuclidean,
        "mahalanobis" => Distance::Mahalanobis,
        "euclidean" => Distance::Euclidean,
        other => panic!("unknown fixture distance `{other}`"),
    }
}

/// Assert the energy solve succeeded before its objective is compared. A
/// one-sided objective comparison alone would accept a degenerate zero-weight
/// failure, whose objective of zero falls below every positive reference, so the
/// harness first requires a solved status (the backend reports this only when the
/// primal residual, including the group-sum rows, is within tolerance) and finite
/// non-negative weights.
fn assert_energy_solved(path: &Path, result: &EnergyResult) {
    assert!(
        result.converged,
        "{}: solve did not reach a solved status (status {})",
        path.display(),
        result.status,
    );
    assert!(
        result.weights.iter().all(|w| w.is_finite() && *w >= 0.0),
        "{}: solved weights are not all finite and non-negative",
        path.display(),
    );
    assert!(
        result.weights.iter().sum::<f64>() > 0.0,
        "{}: solved weights sum to zero",
        path.display(),
    );
}

/// Assert an achieved quadratic-program objective sits at or below the reference
/// plus a relative tolerance, the parity criterion for the energy family: both
/// implementations minimize the same energy loss, so ours must not exceed the
/// reference by more than the tolerance.
fn compare_objective(path: &Path, ours: f64, reference: f64, rel_tol: f64) {
    let scale = reference.abs().max(1.0);
    assert!(
        ours <= reference + rel_tol * scale,
        "{}: objective {ours} exceeds reference {reference} + {rel_tol} relative",
        path.display(),
    );
}

/// Solve a binary or multi-category energy balancing fixture and compare the
/// achieved objective. Schema adds `distance`, `estimand`, `improved`, `min_w`,
/// `lambda`, `treat` (zero/one for binary, zero-based level for multi), an
/// optional `focal`, and `expected_obj`; a multi fixture sets `kind` to
/// `energy_multi`.
fn check_energy_fixture(path: &Path, f: &Value, multi: bool) {
    let n = scalar(f, "n") as usize;
    let p = scalar(f, "p") as usize;
    let covs = nums(f, "covs");
    let s = nums(f, "s");
    let levels: Vec<i32> = nums(f, "treat").iter().map(|t| *t as i32).collect();
    let min_w = f["min_w"].as_f64().unwrap_or(1e-8);
    let lambda = f["lambda"].as_f64().unwrap_or(1e-4);
    let improved = f["improved"].as_bool().unwrap_or(true);
    let estimand_name = f["estimand"].as_str().unwrap_or("ate");

    let (n_levels, estimand) = if multi {
        let n_levels = levels.iter().copied().max().map_or(0, |m| m + 1).max(0) as usize;
        let focal = f["focal"].as_f64().map(|x| x as usize).unwrap_or(0);
        let estimand = match estimand_name {
            "ate" => EnergyEstimand::Ate { improved },
            "att" | "atc" => EnergyEstimand::Focal { focal },
            other => panic!("unknown energy estimand `{other}`"),
        };
        (n_levels, estimand)
    } else {
        let estimand = match estimand_name {
            "ate" => EnergyEstimand::Ate { improved },
            "att" => EnergyEstimand::Focal { focal: 1 },
            "atc" => EnergyEstimand::Focal { focal: 0 },
            other => panic!("unknown energy estimand `{other}`"),
        };
        (2usize, estimand)
    };

    let inputs = EnergyDiscreteInputs {
        covs: &covs,
        n,
        p,
        distance: distance_from(f),
        levels: &levels,
        n_levels,
        estimand,
        s: &s,
        min_weight: min_w,
        lambda,
        moment_covs: &[],
        n_moments: 0,
        targets: &[],
        tols: &[],
        threads: 1,
        qp: QpOptions::default(),
    };
    let result = solve_energy_discrete(&inputs, &|| false);
    assert_energy_solved(path, &result);
    let reference = scalar(f, "expected_obj");
    let rel_tol = f["rel_tol"].as_f64().unwrap_or(1e-6);
    compare_objective(path, result.objective, reference, rel_tol);
}

/// Solve a continuous-exposure energy balancing fixture and compare the achieved
/// objective. Schema carries `treat` as the exposure and `dimension_adj`.
fn check_energy_cont_fixture(path: &Path, f: &Value) {
    let n = scalar(f, "n") as usize;
    let p = scalar(f, "p") as usize;
    let covs = nums(f, "covs");
    let treat = nums(f, "treat");
    let s = nums(f, "s");
    let min_w = f["min_w"].as_f64().unwrap_or(1e-8);
    let lambda = f["lambda"].as_f64().unwrap_or(1e-4);
    let dimension_adj = f["dimension_adj"].as_bool().unwrap_or(true);

    let inputs = EnergyContInputs {
        covs: &covs,
        n,
        p,
        treat: &treat,
        distance: distance_from(f),
        s: &s,
        min_weight: min_w,
        lambda,
        dimension_adj,
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
    let result = solve_energy_cont(&inputs, &|| false);
    assert_energy_solved(path, &result);
    let reference = scalar(f, "expected_obj");
    let rel_tol = f["rel_tol"].as_f64().unwrap_or(1e-6);
    compare_objective(path, result.objective, reference, rel_tol);
}

/// Assert a stable balancing solve succeeded before its objective is compared, so
/// a degenerate failure whose objective is zero cannot pass the one-sided
/// comparison. The check requires a solved status and finite non-negative weights
/// that sum above zero.
fn assert_sbw_solved(path: &Path, result: &SbwResult) {
    assert!(
        result.converged,
        "{}: solve did not reach a solved status (status {})",
        path.display(),
        result.status,
    );
    assert!(
        result.weights.iter().all(|w| w.is_finite() && *w >= 0.0),
        "{}: solved weights are not all finite and non-negative",
        path.display(),
    );
    assert!(
        result.weights.iter().sum::<f64>() > 0.0,
        "{}: solved weights sum to zero",
        path.display(),
    );
}

/// Assert the returned discrete weights satisfy the fixture's constraint box: each
/// active group's group-normalized total is one and each column's group-normalized
/// weighted mean sits inside `target +/- tol`. The design's tolerance policy
/// (rust-architecture section 6) requires constraint satisfaction alongside the
/// one-sided objective bound, since a lower objective from a too-loose feasible set
/// is otherwise indistinguishable from solver quality.
fn assert_sbw_discrete_constraints(
    path: &Path,
    levels: &[i32],
    n_levels: usize,
    estimand: &SbwEstimand,
    s: &[f64],
    moment_covs: &[f64],
    targets: &[f64],
    tols: &[f64],
    weights: &[f64],
) {
    let n = levels.len();
    let q = targets.len();
    let s_norm = group_normalized(s, levels, n_levels);
    let mut n_t = vec![0usize; n_levels];
    for &g in levels {
        if g >= 0 {
            n_t[g as usize] += 1;
        }
    }
    let focal = match estimand {
        SbwEstimand::Focal { focal } => Some(*focal),
        SbwEstimand::Ate => None,
    };
    for t in 0..n_levels {
        if Some(t) == focal || n_t[t] == 0 {
            continue;
        }
        let swnt = |i: usize| s_norm[i] / n_t[t] as f64;
        let group_sum: f64 = (0..n)
            .filter(|&i| levels[i] == t as i32)
            .map(|i| swnt(i) * weights[i])
            .sum();
        assert!(
            (group_sum - 1.0).abs() < 1e-4,
            "{}: group {t} normalized total {group_sum} is not one",
            path.display(),
        );
        for c in 0..q {
            let mean: f64 = (0..n)
                .filter(|&i| levels[i] == t as i32)
                .map(|i| swnt(i) * weights[i] * moment_covs[c * n + i])
                .sum();
            assert!(
                mean >= targets[c] - tols[c] - 1e-4 && mean <= targets[c] + tols[c] + 1e-4,
                "{}: group {t} column {c} weighted mean {mean} outside {} +/- {}",
                path.display(),
                targets[c],
                tols[c],
            );
        }
    }
}

/// Reliability-weighted variance, matching the exposure standardization the
/// continuous solve uses so the checked correlation is the row the solver bounds.
fn sbw_weighted_variance(x: &[f64], w: &[f64]) -> f64 {
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

/// Assert the returned continuous weights satisfy the fixture's constraints: the
/// sampling-weighted total is n and each covariate's linearized weighted
/// correlation with the exposure sits inside its tolerance.
fn assert_sbw_cont_constraints(
    path: &Path,
    treat: &[f64],
    covs: &[f64],
    s: &[f64],
    tols: &[f64],
    weights: &[f64],
) {
    let n = treat.len();
    let q = tols.len();
    let s_sum: f64 = s.iter().sum();
    let s_scaled: Vec<f64> = if s_sum > 0.0 {
        s.iter().map(|&si| si * n as f64 / s_sum).collect()
    } else {
        vec![1.0; n]
    };
    let a_mean =
        s.iter().zip(treat).map(|(&wi, &ti)| wi * ti).sum::<f64>() / s_sum.max(f64::MIN_POSITIVE);
    let a_var = sbw_weighted_variance(treat, s);
    let a_sd = if a_var > 0.0 { a_var.sqrt() } else { 1.0 };

    let total: f64 = (0..n).map(|i| s_scaled[i] * weights[i]).sum();
    assert!(
        (total - n as f64).abs() < 1e-3 * n as f64,
        "{}: sampling-weighted total {total} is not n",
        path.display(),
    );
    for c in 0..q {
        let corr: f64 = (0..n)
            .map(|i| covs[c * n + i] * ((treat[i] - a_mean) / a_sd) * s_scaled[i] * weights[i])
            .sum::<f64>()
            / n as f64;
        assert!(
            corr.abs() <= tols[c] + 1e-4,
            "{}: column {c} correlation {corr} exceeds tolerance {}",
            path.display(),
            tols[c],
        );
    }
}

fn sbw_norm(f: &Value) -> SbwNorm {
    match f["norm"].as_str().unwrap_or("l2") {
        "l2" => SbwNorm::L2,
        "l1" => SbwNorm::L1,
        "linf" => SbwNorm::Linf,
        other => panic!("unknown sbw norm `{other}`"),
    }
}

/// Solve a binary or multi-category stable balancing fixture and compare the
/// achieved weight-dispersion objective. Schema adds `treat` (zero/one for binary,
/// zero-based level for multi), `estimand`, an optional `focal`, `norm`, `min_w`,
/// the standardized `moment_covs` with `targets` and `tols`, and `expected_obj`; a
/// multi fixture sets `kind` to `sbw_multi`.
fn check_sbw_fixture(path: &Path, f: &Value, multi: bool) {
    let n = scalar(f, "n") as usize;
    let p = scalar(f, "p") as usize;
    let s = nums(f, "s");
    let levels: Vec<i32> = nums(f, "treat").iter().map(|t| *t as i32).collect();
    let moment_covs = nums(f, "moment_covs");
    let targets = nums(f, "targets");
    let tols = nums(f, "tols");
    let min_w = f["min_w"].as_f64().unwrap_or(1e-8);
    let estimand_name = f["estimand"].as_str().unwrap_or("ate");

    let (n_levels, estimand) = if multi {
        let n_levels = levels.iter().copied().max().map_or(0, |m| m + 1).max(0) as usize;
        let focal = f["focal"].as_f64().map(|x| x as usize).unwrap_or(0);
        let estimand = match estimand_name {
            "ate" => SbwEstimand::Ate,
            "att" | "atc" => SbwEstimand::Focal { focal },
            other => panic!("unknown sbw estimand `{other}`"),
        };
        (n_levels, estimand)
    } else {
        let estimand = match estimand_name {
            "ate" => SbwEstimand::Ate,
            "att" => SbwEstimand::Focal { focal: 1 },
            "atc" => SbwEstimand::Focal { focal: 0 },
            other => panic!("unknown sbw estimand `{other}`"),
        };
        (2usize, estimand)
    };

    let inputs = SbwDiscreteInputs {
        n,
        levels: &levels,
        n_levels,
        estimand,
        s: &s,
        min_weight: min_w,
        norm: sbw_norm(f),
        moment_covs: &moment_covs,
        n_moments: p,
        targets: &targets,
        tols: &tols,
        qp: QpOptions::default(),
    };
    let result = solve_sbw_discrete(&inputs, &|| false).expect("supported norm");
    assert_sbw_solved(path, &result);
    assert_sbw_discrete_constraints(
        path,
        &levels,
        n_levels,
        &estimand,
        &s,
        &moment_covs,
        &targets,
        &tols,
        &result.weights,
    );
    let reference = scalar(f, "expected_obj");
    let rel_tol = f["rel_tol"].as_f64().unwrap_or(1e-6);
    compare_objective(path, result.objective, reference, rel_tol);
}

/// Solve a continuous-exposure stable balancing fixture and compare the achieved
/// dispersion objective. Schema carries `treat` as the exposure and the
/// standardized `covs` with per-column `tols`.
fn check_sbw_cont_fixture(path: &Path, f: &Value) {
    let n = scalar(f, "n") as usize;
    let p = scalar(f, "p") as usize;
    let treat = nums(f, "treat");
    let covs = nums(f, "covs");
    let tols = nums(f, "tols");
    let s = nums(f, "s");
    let min_w = f["min_w"].as_f64().unwrap_or(1e-8);

    let inputs = SbwContInputs {
        n,
        treat: &treat,
        covs: &covs,
        n_covs: p,
        s: &s,
        min_weight: min_w,
        norm: sbw_norm(f),
        tols: &tols,
        qp: QpOptions::default(),
    };
    let result = solve_sbw_cont(&inputs, &|| false).expect("supported norm");
    assert_sbw_solved(path, &result);
    assert_sbw_cont_constraints(path, &treat, &covs, &s, &tols, &result.weights);
    let reference = scalar(f, "expected_obj");
    let rel_tol = f["rel_tol"].as_f64().unwrap_or(1e-6);
    compare_objective(path, result.objective, reference, rel_tol);
}

fn check_fixture(path: &Path) {
    let text = std::fs::read_to_string(path).expect("read fixture");
    let f: Value = serde_json::from_str(&text).expect("parse fixture JSON");

    match f["kind"].as_str() {
        Some("ipt") => {
            check_ipt_fixture(path, &f);
            return;
        }
        Some("energy") => {
            check_energy_fixture(path, &f, false);
            return;
        }
        Some("energy_multi") => {
            check_energy_fixture(path, &f, true);
            return;
        }
        Some("energy_cont") => {
            check_energy_cont_fixture(path, &f);
            return;
        }
        Some("cbps") => {
            check_cbps_fixture(path, &f);
            return;
        }
        Some("cbps_multi") => {
            check_cbps_multi_fixture(path, &f);
            return;
        }
        Some("cbps_cont") => {
            check_cbps_cont_fixture(path, &f);
            return;
        }
        Some("sbw") => {
            check_sbw_fixture(path, &f, false);
            return;
        }
        Some("sbw_multi") => {
            check_sbw_fixture(path, &f, true);
            return;
        }
        Some("sbw_cont") => {
            check_sbw_cont_fixture(path, &f);
            return;
        }
        _ => {}
    }

    let n = scalar(&f, "n") as usize;
    let p = scalar(&f, "p") as usize;
    let covs = nums(&f, "covs");
    let targets = nums(&f, "targets");
    let tols = nums(&f, "tols");
    let base = nums(&f, "base");
    let s = nums(&f, "s");
    let n_eff = scalar(&f, "n_eff");
    let expected = nums(&f, "expected_weights");
    let rel_tol = scalar(&f, "rel_tol");

    let inputs = EntropyInputs {
        covs: &covs,
        n,
        p,
        targets: &targets,
        tols: &tols,
        base: &base,
        s: &s,
        n_eff,
        threads: 1,
        max_iter: 500,
        tol: 1e-12,
        solver: EntropySolver::Newton,
    };

    let result = match f["kind"].as_str().expect("fixture kind is a string") {
        "discrete" => {
            let group_idx: Vec<i32> = nums(&f, "group_idx").iter().map(|g| *g as i32).collect();
            solve_discrete(&inputs, &group_idx, &|| false)
        }
        "continuous" => {
            // Marginal-distribution columns are held exactly; the fixtures use
            // exact fits, so an all-zero `dist_ind` leaves every column relaxable
            // and the resulting L1 penalty is zero regardless.
            let dist_ind = f["dist_ind"]
                .as_array()
                .map(|a| a.iter().map(|v| v.as_f64().unwrap_or(0.0) as i32).collect())
                .unwrap_or_else(|| vec![0_i32; p]);
            solve_continuous(&inputs, &dist_ind, &|| false)
        }
        other => panic!("unknown fixture kind `{other}`"),
    };

    compare_weights(path, &result.weights, &expected, rel_tol);
}

#[test]
fn golden_fixtures_match() {
    let dir = golden_dir();
    let entries: Vec<PathBuf> = std::fs::read_dir(&dir)
        .map(|read| {
            read.filter_map(|e| e.ok().map(|e| e.path()))
                .filter(|p| p.extension().is_some_and(|ext| ext == "json"))
                .collect()
        })
        .unwrap_or_default();

    if entries.is_empty() {
        eprintln!(
            "no golden fixtures in {}; data-raw/make-golden.R generates them in the R API step",
            dir.display()
        );
        return;
    }

    for path in entries {
        check_fixture(&path);
    }
}
