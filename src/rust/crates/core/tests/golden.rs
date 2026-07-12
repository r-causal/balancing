//! Golden parity tests for entropy balancing.
//!
//! Fixtures live under `tests/golden/` as JSON files written by
//! `data-raw/make-golden.R` from a pinned reference implementation. Each fixture
//! captures the solver inputs at the matrix level (the constraint matrix,
//! targets, tolerances, and weights) so this harness compares solver parity in
//! isolation from covariate processing. When no fixtures are present the harness
//! passes and reports that they are generated in the R API step.
//!
//! Fixture schema (all arrays are plain JSON numbers):
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

use std::path::{Path, PathBuf};

use balancing_core::links::Link;
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

fn check_fixture(path: &Path) {
    let text = std::fs::read_to_string(path).expect("read fixture");
    let f: Value = serde_json::from_str(&text).expect("parse fixture JSON");

    if matches!(f["kind"].as_str(), Some("ipt")) {
        check_ipt_fixture(path, &f);
        return;
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
