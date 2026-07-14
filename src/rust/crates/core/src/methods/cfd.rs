//! Characteristic function distance balancing: weights that make the kernel mean
//! embeddings of the exposure groups agree.
//!
//! Each kernel embeds the covariate rows in a reproducing-kernel Hilbert space,
//! and the objective drives the reweighted mean embedding of every exposure group
//! toward a target embedding. For the average treatment effect every group is
//! pulled toward the full sample; the improved variant additionally pulls each
//! pair of groups toward each other. Both reduce to a quadratic form in the
//! weights built from the kernel matrix and the group-normalization outer product,
//! exactly as energy balancing does from the distance matrix, because the energy
//! kernel is the negative distance. Characteristic function distance balancing with
//! the energy kernel therefore reproduces energy balancing on the same inputs; the
//! other kernels share the assembly and differ only in the kernel matrix.
//!
//! The kernels other than energy are positive semidefinite by construction, so
//! their quadratic term is positive semidefinite and the interior-point backend is
//! eligible alongside the default ADMM backend. The energy kernel is only
//! conditionally positive semidefinite, so its quadratic term is indefinite and the
//! solve routes to the ADMM backend directly, matching energy balancing.
//!
//! The constraint set reuses the shared assembler: a box on each weight with a
//! minimum-weight floor and pinned zero-sampling-weight units, one group-sum row
//! per reweighted group, and optional moment rows holding a covariate's weighted
//! mean within a tolerance band.

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::dist::kernels::{KernelParams, build_kernel};
use crate::qp::osqp::Osqp;
use crate::qp::{
    Convexity, QpBackend, QpError, QpOptions, QpSolution, QpSpec, QpStatus, solve_psd,
};
use crate::threads::get_pool;

use super::qp_balance::{
    ConstraintBuilder, ZERO_SW, add_diagonal_penalty, doubled_dense, expand_and_floor,
    group_normalized,
};

/// The estimand a characteristic function distance solve targets.
#[derive(Debug, Clone, Copy)]
pub enum CfdEstimand {
    /// Average treatment effect: every group is reweighted toward the full
    /// sample, optionally adding the between-group term of the improved variant.
    Ate { improved: bool },
    /// A focal-group effect (treated or control): the non-focal groups are
    /// reweighted toward the focal group, whose units keep unit weight.
    Focal { focal: usize },
}

/// A reason a characteristic function distance solve cannot proceed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CfdError {
    /// An explicitly requested backend is not compiled into this build; the field
    /// names the missing backend.
    BackendUnavailable(&'static str),
}

impl std::fmt::Display for CfdError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CfdError::BackendUnavailable(name) => write!(
                f,
                "the `{name}` backend is not compiled into this build; rebuild with the `qp-clarabel` feature or choose another backend"
            ),
        }
    }
}

/// Outcome of a characteristic function distance solve.
#[derive(Debug, Clone)]
pub struct CfdResult {
    /// Balancing weights, one per unit, in the original unit order.
    pub weights: Vec<f64>,
    /// Constraint dual variables from the quadratic-program solve.
    pub duals: Vec<f64>,
    /// Whether the solve reached a solved status.
    pub converged: bool,
    /// Whether the solve stopped because a user interrupt was pending.
    pub interrupted: bool,
    /// Iterations performed by the backend.
    pub iterations: usize,
    /// Objective value at the solution.
    pub objective: f64,
    /// The backend that produced the result, `"osqp"` or `"clarabel"`, carried to
    /// R as the solver identity.
    pub backend: &'static str,
    /// Whether the automatic clarabel fallback engaged after an osqp
    /// primal-infeasibility certificate.
    pub fell_back: bool,
    /// The backend's terminal status name.
    pub status: &'static str,
    /// Primal residual from the backend.
    pub pri_res: f64,
    /// Dual residual from the backend.
    pub dua_res: f64,
}

/// Inputs for a binary or multi-category characteristic function distance solve.
pub struct CfdDiscreteInputs<'a> {
    /// Column-major `n` by `p` covariates for the kernel matrix.
    pub covs: &'a [f64],
    /// Number of units.
    pub n: usize,
    /// Number of covariate columns.
    pub p: usize,
    /// Kernel tuning.
    pub kernel: KernelParams<'a>,
    /// Zero-based exposure level of each unit; a negative level excludes the unit.
    pub levels: &'a [i32],
    /// Number of exposure levels.
    pub n_levels: usize,
    /// The estimand.
    pub estimand: CfdEstimand,
    /// Sampling weights.
    pub s: &'a [f64],
    /// Minimum allowable weight.
    pub min_weight: f64,
    /// Weight penalty coefficient.
    pub lambda: f64,
    /// Column-major `n` by `q` moment-constraint covariates, standardized on the
    /// R side; empty for no moment constraints.
    pub moment_covs: &'a [f64],
    /// Number of moment-constraint columns `q`.
    pub n_moments: usize,
    /// Target weighted mean for each moment column.
    pub targets: &'a [f64],
    /// Tolerance band half-width source for each moment column.
    pub tols: &'a [f64],
    /// Worker threads for the kernel assembly.
    pub threads: usize,
    /// Quadratic-program tuning.
    pub qp: QpOptions,
}

/// The interaction factor of the group-normalization outer product for two units
/// in levels `li` and `lj` under the estimand.
///
/// This mirrors energy balancing: for the average treatment effect the standard
/// variant pairs same-level units with weight one; the improved variant adds the
/// between-group term, which raises the same-level weight to the number of levels
/// and sets the cross-level weight to minus one. For a focal estimand only
/// same-level, non-focal units interact.
fn nn_factor(estimand: CfdEstimand, n_levels: usize, li: i32, lj: i32) -> f64 {
    match estimand {
        CfdEstimand::Ate { improved: false } => {
            if li == lj {
                1.0
            } else {
                0.0
            }
        }
        CfdEstimand::Ate { improved: true } => {
            if li == lj {
                n_levels as f64
            } else {
                -1.0
            }
        }
        CfdEstimand::Focal { .. } => {
            if li == lj {
                1.0
            } else {
                0.0
            }
        }
    }
}

/// Solve a binary or multi-category characteristic function distance problem.
pub fn solve_discrete(
    inputs: &CfdDiscreteInputs<'_>,
    interrupt: &dyn Fn() -> bool,
) -> Result<CfdResult, CfdError> {
    let n = inputs.n;
    let discarded: Vec<bool> = inputs.levels.iter().map(|&g| g < 0).collect();
    let kernel = build_kernel(
        inputs.covs,
        n,
        inputs.p,
        inputs.s,
        &inputs.kernel,
        &discarded,
        inputs.threads,
    );
    let s_norm = group_normalized(inputs.s, inputs.levels, inputs.n_levels);

    // Level sizes and each unit's own-column normalization value.
    let mut n_t = vec![0usize; inputs.n_levels];
    for &g in inputs.levels {
        if g >= 0 {
            n_t[g as usize] += 1;
        }
    }
    let swnt: Vec<f64> = inputs
        .levels
        .iter()
        .enumerate()
        .map(|(i, &g)| {
            if g >= 0 && n_t[g as usize] > 0 {
                s_norm[i] / n_t[g as usize] as f64
            } else {
                0.0
            }
        })
        .collect();

    // The active variables and the level rows that carry constraints. For the
    // average treatment effect every present unit is a variable and every level is
    // constrained; for a focal estimand only the non-focal units are variables and
    // only the non-focal levels are constrained.
    let (active, group_levels, focal) = match inputs.estimand {
        CfdEstimand::Ate { .. } => {
            let active: Vec<usize> = (0..n).filter(|&i| inputs.levels[i] >= 0).collect();
            let groups: Vec<usize> = (0..inputs.n_levels).collect();
            (active, groups, None)
        }
        CfdEstimand::Focal { focal } => {
            let active: Vec<usize> = (0..n)
                .filter(|&i| inputs.levels[i] >= 0 && inputs.levels[i] as usize != focal)
                .collect();
            let groups: Vec<usize> = (0..inputs.n_levels).filter(|&t| t != focal).collect();
            (active, groups, Some(focal))
        }
    };
    let nvar = active.len();

    // Quadratic term P = kernel * nn with the group-normalization outer product,
    // built over the active variables. Where energy balancing uses the negated
    // distance, kernel balancing uses the kernel directly; the two coincide for the
    // energy kernel, whose value is the negated distance.
    let mut pmat = vec![0.0; nvar * nvar];
    for a in 0..nvar {
        let ia = active[a];
        for b in 0..nvar {
            let ib = active[b];
            let factor = nn_factor(
                inputs.estimand,
                inputs.n_levels,
                inputs.levels[ia],
                inputs.levels[ib],
            );
            if factor == 0.0 {
                continue;
            }
            let k = kernel[ib * n + ia];
            pmat[b * nvar + a] = k * swnt[ia] * swnt[ib] * factor;
        }
    }
    // The penalty scale is each active unit's own group-normalization value.
    let scale: Vec<f64> = active.iter().map(|&i| swnt[i]).collect();
    add_diagonal_penalty(&mut pmat, nvar, inputs.lambda, &scale);

    // Linear term q_a = cross_a * scale_a, where cross_a is the negated
    // source-weighted kernel similarity of active column a to its target sample.
    // The source runs over the whole sample for the average treatment effect and
    // over the focal group for a focal estimand.
    let (src, mult): (Vec<f64>, f64) = match inputs.estimand {
        CfdEstimand::Ate { .. } => (s_norm.clone(), 2.0 / n as f64),
        CfdEstimand::Focal { .. } => {
            let f = focal.expect("focal estimand carries a focal level");
            let nf = n_t[f].max(1) as f64;
            let src: Vec<f64> = (0..n)
                .map(|i| {
                    if inputs.levels[i] >= 0 && inputs.levels[i] as usize == f {
                        s_norm[i] / nf
                    } else {
                        0.0
                    }
                })
                .collect();
            (src, 2.0)
        }
    };
    let cross = cross_similarity(&kernel, n, &active, &src, mult, inputs.threads);
    let q: Vec<f64> = (0..nvar).map(|a| cross[a] * scale[a]).collect();

    // Constraints: the box, the group sums, and the moment rows.
    let mut builder = ConstraintBuilder::new(nvar);
    let pinned: Vec<bool> = active
        .iter()
        .map(|&i| inputs.s[i].abs() < ZERO_SW)
        .collect();
    builder.add_box(inputs.min_weight, &pinned);

    for &t in &group_levels {
        let coeffs: Vec<f64> = active
            .iter()
            .map(|&i| {
                if inputs.levels[i] as usize == t {
                    swnt[i]
                } else {
                    0.0
                }
            })
            .collect();
        builder.add_dense_row(&coeffs, 1.0, 1.0);
    }

    // The tolerance is split across the constrained groups for the average
    // treatment effect, which bounds each group's mean around the shared target so
    // the between-group difference stays within the tolerance; a focal estimand
    // constrains one side and uses the full tolerance.
    let tol_half = match inputs.estimand {
        CfdEstimand::Ate { .. } => 0.5,
        CfdEstimand::Focal { .. } => 1.0,
    };
    for &t in &group_levels {
        for c in 0..inputs.n_moments {
            let coeffs: Vec<f64> = active
                .iter()
                .map(|&i| {
                    if inputs.levels[i] as usize == t {
                        inputs.moment_covs[c * n + i] * swnt[i]
                    } else {
                        0.0
                    }
                })
                .collect();
            let band = inputs.tols[c] * tol_half;
            builder.add_dense_row(&coeffs, inputs.targets[c] - band, inputs.targets[c] + band);
        }
    }

    let (m, indptr, indices, values, l, u) = builder.finish();
    let convexity = if inputs.kernel.kernel.is_psd() {
        Convexity::Psd
    } else {
        Convexity::Indefinite
    };
    let spec = QpSpec {
        n: nvar,
        m,
        p: doubled_dense(pmat),
        q,
        a_indptr: indptr,
        a_indices: indices,
        a_values: values,
        l,
        u,
        convexity,
    };

    let (solution, backend, fell_back) = route(&spec, &inputs.qp, interrupt)?;
    let weights = expand_and_floor(n, &active, &solution.x, inputs.min_weight);
    Ok(package(weights, &solution, backend, fell_back))
}

/// The negated source-weighted kernel similarity of each active column.
///
/// The value at active column `a` is `-mult` times the source-weighted sum of
/// kernel similarities from every unit to the unit at that column. For the energy
/// kernel, whose value is the negated distance, this reduces to energy balancing's
/// positively signed cross energy. Work is split over the active columns and is
/// deterministic by construction because each column is an independent sum in a
/// fixed index order.
fn cross_similarity(
    kernel: &[f64],
    n: usize,
    active: &[usize],
    src: &[f64],
    mult: f64,
    threads: usize,
) -> Vec<f64> {
    let pool = get_pool(threads);
    pool.install(|| {
        (0..active.len())
            .into_par_iter()
            .map(|a| {
                let j = active[a];
                let mut acc = 0.0;
                for i in 0..n {
                    acc += src[i] * kernel[j * n + i];
                }
                -mult * acc
            })
            .collect()
    })
}

/// Solve the spec through the backend policy. A positive-semidefinite kernel routes
/// through the shared osqp-primary policy with an eligible clarabel fallback; the
/// indefinite energy kernel routes to osqp directly, as energy balancing does,
/// since clarabel cannot accept it. An explicitly requested backend that is not
/// compiled in propagates as an error naming the missing feature; any other setup
/// failure falls back to a degenerate solution so the R layer raises the
/// convergence condition rather than propagating a bug as a panic.
fn route(
    spec: &QpSpec,
    opts: &QpOptions,
    interrupt: &dyn Fn() -> bool,
) -> Result<(QpSolution, &'static str, bool), CfdError> {
    if spec.convexity == Convexity::Indefinite {
        let solution = Osqp
            .solve(spec, opts, interrupt)
            .unwrap_or_else(|_| degenerate(spec));
        return Ok((solution, "osqp", false));
    }
    match solve_psd(spec, opts, interrupt) {
        Ok(routed) => Ok((routed.solution, routed.backend, routed.fell_back)),
        Err(QpError::BackendUnavailable(name)) => Err(CfdError::BackendUnavailable(name)),
        Err(_) => Ok((degenerate(spec), "osqp", false)),
    }
}

/// Pack floored weights and a solution into the result record.
fn package(
    weights: Vec<f64>,
    solution: &QpSolution,
    backend: &'static str,
    fell_back: bool,
) -> CfdResult {
    CfdResult {
        weights,
        duals: solution.duals.clone(),
        converged: solution.status.is_solved(),
        interrupted: solution.interrupted,
        iterations: solution.iterations,
        objective: solution.obj,
        backend,
        fell_back,
        status: solution.status.as_str(),
        pri_res: solution.pri_res,
        dua_res: solution.dua_res,
    }
}

/// A degenerate solution for a spec a backend refused to set up: zero weights
/// carrying the failure status so the R layer raises the convergence condition.
fn degenerate(spec: &QpSpec) -> QpSolution {
    QpSolution {
        x: vec![0.0; spec.n],
        duals: vec![0.0; spec.m],
        status: QpStatus::NonConvex,
        iterations: 0,
        obj: 0.0,
        pri_res: 0.0,
        dua_res: 0.0,
        interrupted: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::Distance;
    use crate::dist::kernels::Kernel;
    use crate::methods::energy::{EnergyDiscreteInputs, EnergyEstimand, solve_discrete as energy};

    fn kernel_params(kernel: Kernel) -> KernelParams<'static> {
        KernelParams {
            kernel,
            bw_scale: 1.0,
            smoothness: 1.5,
            t_proj: &[],
            n_draws: 0,
        }
    }

    fn discrete_inputs<'a>(
        covs: &'a [f64],
        levels: &'a [i32],
        n_levels: usize,
        s: &'a [f64],
        kernel: KernelParams<'a>,
        estimand: CfdEstimand,
    ) -> CfdDiscreteInputs<'a> {
        CfdDiscreteInputs {
            covs,
            n: levels.len(),
            p: covs.len() / levels.len(),
            kernel,
            levels,
            n_levels,
            estimand,
            s,
            min_weight: 1e-8,
            lambda: 1e-4,
            moment_covs: &[],
            n_moments: 0,
            targets: &[],
            tols: &[],
            threads: 1,
            qp: QpOptions::default(),
        }
    }

    // A confounded eight-by-two binary design shared by the property tests.
    fn confounded() -> (Vec<f64>, [i32; 8], Vec<f64>) {
        let covs = vec![
            1.5, 1.2, 0.9, 1.1, -1.0, -1.3, -0.7, -1.1, // covariate 1
            0.5, 0.2, -0.1, 0.3, -0.6, -0.2, 0.1, -0.4, // covariate 2
        ];
        let levels = [1, 1, 1, 1, 0, 0, 0, 0];
        let s = vec![1.0; 8];
        (covs, levels, s)
    }

    #[test]
    fn the_energy_kernel_reproduces_the_energy_method() {
        // With the energy kernel, whose value is the negated distance, the kernel
        // assembly is the energy assembly, so the weights match energy balancing on
        // the same inputs and constraints.
        let (covs, levels, s) = confounded();
        let cfd = solve_discrete(
            &discrete_inputs(
                &covs,
                &levels,
                2,
                &s,
                kernel_params(Kernel::Energy),
                CfdEstimand::Ate { improved: true },
            ),
            &|| false,
        )
        .unwrap();
        assert!(cfd.converged, "cfd status {}", cfd.status);

        let energy_inputs = EnergyDiscreteInputs {
            covs: &covs,
            n: 8,
            p: 2,
            distance: Distance::ScaledEuclidean,
            levels: &levels,
            n_levels: 2,
            estimand: EnergyEstimand::Ate { improved: true },
            s: &s,
            min_weight: 1e-8,
            lambda: 1e-4,
            moment_covs: &[],
            n_moments: 0,
            targets: &[],
            tols: &[],
            threads: 1,
            qp: QpOptions::default(),
        };
        let en = energy(&energy_inputs, &|| false);
        assert!(en.converged, "energy status {}", en.status);
        for (a, b) in cfd.weights.iter().zip(&en.weights) {
            assert!((a - b).abs() < 1e-6, "cfd {a} != energy {b}");
        }
    }

    #[test]
    fn a_gaussian_ate_normalizes_each_group_to_its_size() {
        let (covs, levels, s) = confounded();
        let result = solve_discrete(
            &discrete_inputs(
                &covs,
                &levels,
                2,
                &s,
                kernel_params(Kernel::Gaussian),
                CfdEstimand::Ate { improved: true },
            ),
            &|| false,
        )
        .unwrap();
        assert!(result.converged, "status {}", result.status);
        for level in [0, 1] {
            let sum: f64 = (0..8)
                .filter(|&i| levels[i] == level)
                .map(|i| result.weights[i])
                .sum();
            assert!((sum - 4.0).abs() < 1e-3, "level {level} sum {sum}");
        }
        assert!(result.weights.iter().all(|&w| w >= 1e-8 - 1e-12));
        assert_eq!(result.backend, "osqp");
    }

    #[test]
    fn an_att_pins_the_focal_group_at_unit_weight() {
        let (covs, levels, s) = confounded();
        let result = solve_discrete(
            &discrete_inputs(
                &covs,
                &levels,
                2,
                &s,
                kernel_params(Kernel::Gaussian),
                CfdEstimand::Focal { focal: 1 },
            ),
            &|| false,
        )
        .unwrap();
        assert!(result.converged, "status {}", result.status);
        for i in 0..8 {
            if levels[i] == 1 {
                assert!(
                    (result.weights[i] - 1.0).abs() < 1e-9,
                    "focal weight {}",
                    result.weights[i]
                );
            }
        }
        let sum_c: f64 = (0..8)
            .filter(|&i| levels[i] == 0)
            .map(|i| result.weights[i])
            .sum();
        assert!((sum_c - 4.0).abs() < 1e-3, "control sum {sum_c}");
    }

    #[test]
    fn the_improved_variant_differs_from_the_plain_variant() {
        let (covs, levels, s) = confounded();
        let improved = solve_discrete(
            &discrete_inputs(
                &covs,
                &levels,
                2,
                &s,
                kernel_params(Kernel::Gaussian),
                CfdEstimand::Ate { improved: true },
            ),
            &|| false,
        )
        .unwrap();
        let plain = solve_discrete(
            &discrete_inputs(
                &covs,
                &levels,
                2,
                &s,
                kernel_params(Kernel::Gaussian),
                CfdEstimand::Ate { improved: false },
            ),
            &|| false,
        )
        .unwrap();
        assert!(improved.converged && plain.converged);
        let differ = improved
            .weights
            .iter()
            .zip(&plain.weights)
            .any(|(a, b)| (a - b).abs() > 1e-4);
        assert!(differ, "improved and plain weights should differ");
    }

    #[test]
    fn distinct_kernels_produce_distinct_weights() {
        let (covs, levels, s) = confounded();
        let solve = |kernel: Kernel| {
            solve_discrete(
                &discrete_inputs(
                    &covs,
                    &levels,
                    2,
                    &s,
                    kernel_params(kernel),
                    CfdEstimand::Ate { improved: true },
                ),
                &|| false,
            )
            .unwrap()
            .weights
        };
        let gaussian = solve(Kernel::Gaussian);
        let laplace = solve(Kernel::Laplace);
        let energy = solve(Kernel::Energy);
        let differ = |a: &[f64], b: &[f64]| a.iter().zip(b).any(|(x, y)| (x - y).abs() > 1e-4);
        assert!(differ(&gaussian, &laplace), "gaussian equals laplace");
        assert!(differ(&gaussian, &energy), "gaussian equals energy");
        assert!(differ(&laplace, &energy), "laplace equals energy");
    }

    #[test]
    fn a_categorical_ate_normalizes_each_group_to_its_size() {
        let n = 9;
        let covs = vec![
            1.0, 1.2, 0.8, 0.1, -0.1, 0.2, -1.0, -1.2, -0.9, // covariate 1
            0.3, 0.1, 0.5, 1.1, 0.9, 1.3, -0.4, -0.2, -0.6, // covariate 2
        ];
        let levels = [0, 0, 0, 1, 1, 1, 2, 2, 2];
        let s = vec![1.0; n];
        let result = solve_discrete(
            &discrete_inputs(
                &covs,
                &levels,
                3,
                &s,
                kernel_params(Kernel::Gaussian),
                CfdEstimand::Ate { improved: true },
            ),
            &|| false,
        )
        .unwrap();
        assert!(result.converged, "status {}", result.status);
        for level in 0..3 {
            let sum: f64 = (0..n)
                .filter(|&i| levels[i] == level)
                .map(|i| result.weights[i])
                .sum();
            assert!((sum - 3.0).abs() < 1e-3, "level {level} sum {sum}");
        }
    }

    #[test]
    fn a_pending_interrupt_surfaces_on_the_result() {
        let (covs, levels, s) = confounded();
        let mut inputs = discrete_inputs(
            &covs,
            &levels,
            2,
            &s,
            kernel_params(Kernel::Gaussian),
            CfdEstimand::Ate { improved: true },
        );
        inputs.qp = QpOptions {
            chunk_iters: 1,
            ..QpOptions::default()
        };
        let result = solve_discrete(&inputs, &|| true).unwrap();
        assert!(result.interrupted);
        assert_eq!(result.status, "interrupted");
        assert!(!result.converged);
    }

    #[test]
    fn the_solve_is_deterministic() {
        let (covs, levels, s) = confounded();
        let solve = || {
            solve_discrete(
                &discrete_inputs(
                    &covs,
                    &levels,
                    2,
                    &s,
                    kernel_params(Kernel::Gaussian),
                    CfdEstimand::Ate { improved: true },
                ),
                &|| false,
            )
            .unwrap()
        };
        let a = solve();
        let b = solve();
        for (x, y) in a.weights.iter().zip(&b.weights) {
            assert_eq!(x.to_bits(), y.to_bits());
        }
    }

    #[test]
    #[cfg(feature = "qp-clarabel")]
    fn a_psd_kernel_solves_through_clarabel_when_requested() {
        use crate::qp::QpBackendChoice;
        let (covs, levels, s) = confounded();
        let mut inputs = discrete_inputs(
            &covs,
            &levels,
            2,
            &s,
            kernel_params(Kernel::Gaussian),
            CfdEstimand::Ate { improved: true },
        );
        inputs.qp.backend = QpBackendChoice::Clarabel;
        let result = solve_discrete(&inputs, &|| false).unwrap();
        assert!(result.converged, "status {}", result.status);
        assert_eq!(result.backend, "clarabel");
    }

    #[test]
    #[cfg(not(feature = "qp-clarabel"))]
    fn an_explicit_clarabel_choice_without_the_feature_errors() {
        use crate::qp::QpBackendChoice;
        let (covs, levels, s) = confounded();
        let mut inputs = discrete_inputs(
            &covs,
            &levels,
            2,
            &s,
            kernel_params(Kernel::Gaussian),
            CfdEstimand::Ate { improved: true },
        );
        inputs.qp.backend = QpBackendChoice::Clarabel;
        let err = solve_discrete(&inputs, &|| false).unwrap_err();
        assert_eq!(err, CfdError::BackendUnavailable("clarabel"));
    }
}
