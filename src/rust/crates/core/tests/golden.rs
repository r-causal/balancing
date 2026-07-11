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

use balancing_core::methods::entropy::{
    EntropyInputs, EntropySolver, solve_continuous, solve_discrete,
};
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

fn check_fixture(path: &Path) {
    let text = std::fs::read_to_string(path).expect("read fixture");
    let f: Value = serde_json::from_str(&text).expect("parse fixture JSON");

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

    let scale = expected
        .iter()
        .fold(0.0_f64, |m, w| m.max(w.abs()))
        .max(1.0);
    for i in 0..n {
        let diff = (result.weights[i] - expected[i]).abs();
        assert!(
            diff <= rel_tol * scale,
            "{}: weight[{i}] = {} expected {} (diff {diff})",
            path.display(),
            result.weights[i],
            expected[i]
        );
    }
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
