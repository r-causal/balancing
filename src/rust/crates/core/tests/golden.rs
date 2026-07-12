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

use balancing_core::links::Link;
use balancing_core::methods::cbps::{
    CbpsContInputs, CbpsEstimand, CbpsInputs, CbpsMultiInputs, solve as solve_cbps,
    solve_cont as solve_cbps_cont, solve_multi as solve_cbps_multi,
};
use balancing_core::methods::entropy::{
    EntropyInputs, EntropySolver, solve_continuous, solve_discrete,
};
use balancing_core::methods::ipt::{IptEstimand, IptInputs, solve as solve_ipt};
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

fn check_fixture(path: &Path) {
    let text = std::fs::read_to_string(path).expect("read fixture");
    let f: Value = serde_json::from_str(&text).expect("parse fixture JSON");

    match f["kind"].as_str() {
        Some("ipt") => {
            check_ipt_fixture(path, &f);
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
