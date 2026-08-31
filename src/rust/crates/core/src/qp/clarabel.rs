//! The Clarabel backend: an interior-point solver for convex problems.
//!
//! Interior-point methods require a positive-semidefinite quadratic form, so
//! this backend refuses a spec tagged `Indefinite` up front rather than
//! returning a meaningless iterate, and it refuses a spec with no decision
//! variables for the same reason. It is admissible for the diagonal and
//! kernel forms of the other quadratic-program methods; energy balancing always
//! routes to OSQP. Two-sided bounds are expressed in the conic form Clarabel
//! consumes: an equality row becomes a zero-cone row, and a finite lower or
//! upper bound becomes a nonnegative-cone row, with the lower bound negated so
//! `A x >= l` reads as `(-A) x <= -l`. That negation is undone when the conic
//! duals are gathered back into one multiplier per original constraint, so the
//! reported duals carry the same orientation OSQP reports.
//!
//! Clarabel takes a termination callback, so the interrupt is polled once per
//! interior-point iteration rather than at the chunk boundaries the OSQP backend
//! has to fall back on.

use std::ffi::{c_int, c_void};

use clarabel::algebra::CscMatrix;
use clarabel::solver::implementations::default::ffi::DefaultInfoFFI;
use clarabel::solver::{
    DefaultSettings, DefaultSolver, IPSolver, SolverStatus, SupportedConeT::NonnegativeConeT,
    SupportedConeT::ZeroConeT,
};

use super::{
    Convexity, QpBackend, QpError, QpOptions, QpSolution, QpSpec, QpStatus, objective,
    upper_triangular_csc,
};

/// The Clarabel backend.
pub struct Clarabel;

/// A single conic row assembled from an original two-sided constraint: the row
/// coefficients over the `n` variables, the right-hand side, which original
/// constraint it came from so the dual can be mapped back, and whether the
/// coefficients were negated to turn a lower bound into an upper one.
struct ConicRow {
    coeffs: Vec<(usize, f64)>,
    rhs: f64,
    source: usize,
    negated: bool,
}

impl QpBackend for Clarabel {
    fn name(&self) -> &'static str {
        "clarabel"
    }

    fn solve(
        &self,
        spec: &QpSpec,
        opts: &QpOptions,
        interrupt: &dyn Fn() -> bool,
    ) -> Result<QpSolution, QpError> {
        if spec.convexity == Convexity::Indefinite {
            return Err(QpError::Indefinite);
        }

        let n = spec.n;
        // A spec with no decision variables has nothing to solve for, and Clarabel's
        // own dimension check accepts it. With constraint rows still present it
        // returns an infeasibility certificate against rows no variable can satisfy;
        // with no rows either, the KKT system is empty and its factorization indexes
        // the first entry of nothing. The OSQP backend refuses the same spec,
        // and every method turns that refusal into a degenerate result carrying a
        // failure status, so refusing it here keeps the two backends interchangeable
        // and keeps the R layer raising one condition rather than two.
        if n == 0 {
            return Err(QpError::Setup(
                "the problem has no decision variables".to_string(),
            ));
        }

        let (p_indptr, p_indices, p_values) = upper_triangular_csc(&spec.p, n);
        let p_csc = CscMatrix::new(n, n, p_indptr, p_indices, p_values);

        // Row-major coefficients per constraint, gathered from the column-major
        // spec so each two-sided bound can be split into its conic rows.
        let rows_by_constraint = constraint_rows(spec);

        // Zero-cone (equality) rows first, then nonnegative-cone (inequality)
        // rows, as Clarabel expects cones grouped and ordered.
        let mut zero_rows: Vec<ConicRow> = Vec::new();
        let mut nonneg_rows: Vec<ConicRow> = Vec::new();
        for (c, coeffs) in rows_by_constraint.iter().enumerate() {
            let l = spec.l[c];
            let u = spec.u[c];
            if (l - u).abs() <= 0.0 {
                zero_rows.push(ConicRow {
                    coeffs: coeffs.clone(),
                    rhs: u,
                    source: c,
                    negated: false,
                });
                continue;
            }
            if u.is_finite() {
                nonneg_rows.push(ConicRow {
                    coeffs: coeffs.clone(),
                    rhs: u,
                    source: c,
                    negated: false,
                });
            }
            if l.is_finite() {
                let flipped: Vec<(usize, f64)> = coeffs.iter().map(|&(i, v)| (i, -v)).collect();
                nonneg_rows.push(ConicRow {
                    coeffs: flipped,
                    rhs: -l,
                    source: c,
                    negated: true,
                });
            }
        }

        let n_zero = zero_rows.len();
        let n_nonneg = nonneg_rows.len();
        let all_rows: Vec<ConicRow> = zero_rows.into_iter().chain(nonneg_rows).collect();
        let m = all_rows.len();

        let a_csc = rows_to_csc(&all_rows, m, n);
        let b: Vec<f64> = all_rows.iter().map(|r| r.rhs).collect();

        let mut cones = Vec::new();
        if n_zero > 0 {
            cones.push(ZeroConeT(n_zero));
        }
        if n_nonneg > 0 {
            cones.push(NonnegativeConeT(n_nonneg));
        }

        let settings = solver_settings(opts);

        // The callback reads the interrupt closure through this local, so it is
        // declared ahead of the solver and never moved while the solve runs.
        let mut poll: &dyn Fn() -> bool = interrupt;
        let mut solver = DefaultSolver::new(&p_csc, &spec.q, &a_csc, &b, &cones, settings)
            .map_err(|e| QpError::Setup(format!("{e:?}")))?;
        solver.set_termination_callback_c(poll_interrupt, (&raw mut poll).cast::<c_void>());
        solver.solve();

        let x = solver.solution.x.clone();
        // Map the conic duals back to one multiplier per original constraint. Every
        // conic dual is nonnegative, so the orientation has to be restored here: the
        // `l <= A x <= u` convention both backends report in signs a multiplier
        // positive when the upper side is active and negative when the lower side
        // is, and the lower side reached the cone through negated coefficients.
        // Combining the two sides as `z_upper - z_lower` is what makes the reported
        // dual satisfy `P x + q + A' y = 0`.
        let mut duals = vec![0.0; spec.m];
        for (row, z) in all_rows.iter().zip(solver.solution.z.iter()) {
            if row.negated {
                duals[row.source] -= *z;
            } else {
                duals[row.source] += *z;
            }
        }
        let obj = objective(spec, &x);
        // Clarabel reports a callback stop with its own status, which is the only
        // way the interrupted status arises here, so the flag reads off the mapping
        // rather than being tracked separately.
        let status = map_status(solver.solution.status);
        Ok(QpSolution {
            x,
            duals,
            status,
            iterations: solver.solution.iterations as usize,
            obj,
            pri_res: solver.info.res_primal,
            dua_res: solver.info.res_dual,
            interrupted: status == QpStatus::Interrupted,
        })
    }
}

/// Clarabel's per-iteration termination callback, which reports the caller's
/// interrupt closure so a pending interrupt stops the solve at the next iteration
/// rather than after it.
///
/// The closure reaches the callback through the C form rather than the Rust form
/// because the Rust form demands a `'static` callback, while the closure a backend
/// receives is borrowed for the length of the solve alone. `data` is the address
/// of the `&dyn Fn() -> bool` that [`Clarabel::solve`] holds in its own frame, and
/// the solver information record is not consulted: the decision to stop belongs to
/// the caller, not to the iterate.
extern "C" fn poll_interrupt(_info: *const DefaultInfoFFI<f64>, data: *mut c_void) -> c_int {
    // SAFETY: `solve` registers this callback with the address of a live
    // `&dyn Fn() -> bool` local that outlives the solver it registers on, Clarabel
    // stores that address and hands it back unchanged, and it calls the callback
    // synchronously from the solving thread, so no other reference to the local
    // exists while the solve runs.
    let interrupt = unsafe { &*(data as *const &dyn Fn() -> bool) };
    c_int::from(interrupt())
}

/// Translate the shared tuning into Clarabel settings.
///
/// A method's convergence tolerance is documented as setting both the absolute
/// and the relative solver tolerance, which is what OSQP receives through its
/// `eps_abs` and `eps_rel` pair. Clarabel splits the same idea across a duality
/// gap measured absolutely and relatively and a feasibility check measured only
/// absolutely, so the absolute request drives `tol_gap_abs` and `tol_feas` while
/// the relative request drives `tol_gap_rel`. Leaving the relative gap at its own
/// 1e-8 default would silently cap the request, because Clarabel's termination
/// test accepts either tolerance: a request for 1e-12 would terminate as soon as
/// the relative gap fell below 1e-8. Clarabel has no relative feasibility
/// tolerance to map, and the infeasibility-certificate tolerances stay at their
/// defaults, matching the OSQP backend leaving `eps_prim_inf` and `eps_dual_inf`
/// alone.
fn solver_settings(opts: &QpOptions) -> DefaultSettings<f64> {
    DefaultSettings::<f64> {
        verbose: false,
        max_iter: opts.max_iter.min(u32::MAX as usize) as u32,
        tol_gap_abs: opts.eps_abs,
        tol_gap_rel: opts.eps_rel,
        tol_feas: opts.eps_abs,
        ..Default::default()
    }
}

/// Gather the coefficients of each constraint (row of `A`) from the column-major
/// CSC spec.
fn constraint_rows(spec: &QpSpec) -> Vec<Vec<(usize, f64)>> {
    let mut rows = vec![Vec::new(); spec.m];
    for col in 0..spec.n {
        let start = spec.a_indptr[col];
        let end = spec.a_indptr[col + 1];
        for k in start..end {
            let row = spec.a_indices[k];
            rows[row].push((col, spec.a_values[k]));
        }
    }
    rows
}

/// Build a Clarabel CSC matrix from conic rows given as `(column, value)` lists.
fn rows_to_csc(rows: &[ConicRow], m: usize, n: usize) -> CscMatrix<f64> {
    // Collect column-major triplets, then compress.
    let mut cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for (r, row) in rows.iter().enumerate() {
        for &(col, val) in &row.coeffs {
            cols[col].push((r, val));
        }
    }
    let mut colptr = Vec::with_capacity(n + 1);
    let mut rowval = Vec::new();
    let mut nzval = Vec::new();
    colptr.push(0);
    for col in cols.iter_mut() {
        col.sort_by_key(|&(r, _)| r);
        for &(r, v) in col.iter() {
            rowval.push(r);
            nzval.push(v);
        }
        colptr.push(rowval.len());
    }
    CscMatrix::new(m, n, colptr, rowval, nzval)
}

/// Translate a Clarabel status into the backend-agnostic [`QpStatus`].
///
/// The match is exhaustive on purpose. The R layer chooses its advice from the
/// status, so a status folded into a neighbour sends the caller after the wrong
/// knob: a numerical breakdown or a stalled iterate reported as the iteration cap
/// invites raising `max_iterations`, which cannot help either. Listing every
/// variant also turns a new Clarabel status into a compile error rather than a
/// silent misclassification.
fn map_status(status: SolverStatus) -> QpStatus {
    match status {
        SolverStatus::Solved => QpStatus::Solved,
        SolverStatus::AlmostSolved => QpStatus::SolvedInaccurate,
        SolverStatus::PrimalInfeasible | SolverStatus::AlmostPrimalInfeasible => {
            QpStatus::PrimalInfeasible
        }
        SolverStatus::DualInfeasible | SolverStatus::AlmostDualInfeasible => {
            QpStatus::DualInfeasible
        }
        SolverStatus::MaxIterations => QpStatus::MaxIter,
        SolverStatus::MaxTime => QpStatus::TimeLimit,
        SolverStatus::NumericalError => QpStatus::NumericalError,
        SolverStatus::InsufficientProgress => QpStatus::InsufficientProgress,
        SolverStatus::Unsolved => QpStatus::NotSolved,
        SolverStatus::CallbackTerminated => QpStatus::Interrupted,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qp::PMat;
    use crate::qp::osqp::Osqp;

    fn simplex_spec(convexity: Convexity) -> QpSpec {
        QpSpec {
            n: 2,
            m: 3,
            p: PMat::Diagonal(vec![1.0, 1.0]),
            q: vec![0.0, 0.0],
            a_indptr: vec![0, 2, 4],
            a_indices: vec![0, 2, 1, 2],
            a_values: vec![1.0, 1.0, 1.0, 1.0],
            l: vec![0.0, 0.0, 1.0],
            u: vec![f64::INFINITY, f64::INFINITY, 1.0],
            convexity,
        }
    }

    #[test]
    fn rejects_an_indefinite_spec_up_front() {
        let spec = simplex_spec(Convexity::Indefinite);
        let err = Clarabel
            .solve(&spec, &QpOptions::default(), &|| false)
            .unwrap_err();
        assert!(matches!(err, QpError::Indefinite));
    }

    #[test]
    fn solves_a_psd_simplex_projection() {
        let spec = simplex_spec(Convexity::Psd);
        let sol = Clarabel
            .solve(&spec, &QpOptions::default(), &|| false)
            .expect("psd spec solves");
        assert!(sol.status.is_solved(), "status {:?}", sol.status);
        assert!((sol.x[0] - 0.5).abs() < 1e-6, "x0 = {}", sol.x[0]);
        assert!((sol.x[1] - 0.5).abs() < 1e-6, "x1 = {}", sol.x[1]);
        assert_eq!(sol.duals.len(), 3);
    }

    /// A spec whose optimum activates one lower bound and one upper bound, with a
    /// dual vector available in closed form.
    ///
    /// Minimize `0.5 (x0^2 + x1^2 + x2^2) - x0 - 2 x1 - 3 x2` subject to
    /// `x0 + x1 + x2 = 5`, `x0 >= 1.5`, and `x2 <= 2`. Dropping the two
    /// inequalities gives `x = (2/3, 5/3, 8/3)`, so both of them bind, and the
    /// optimum is `x = (1.5, 1.5, 2)`. Under the `l <= A x <= u` convention the
    /// stationarity condition `P x + q + A' y = 0` then fixes
    /// `y = (0.5, -1, 0.5)`: negative on the lower-active row, positive on the
    /// upper-active row.
    fn hand_kkt_spec() -> QpSpec {
        QpSpec {
            n: 3,
            m: 3,
            p: PMat::Diagonal(vec![1.0, 1.0, 1.0]),
            q: vec![-1.0, -2.0, -3.0],
            // Rows: the sum equality, the x0 floor, the x2 ceiling.
            a_indptr: vec![0, 2, 3, 5],
            a_indices: vec![0, 1, 0, 0, 2],
            a_values: vec![1.0, 1.0, 1.0, 1.0, 1.0],
            l: vec![5.0, 1.5, f64::NEG_INFINITY],
            u: vec![5.0, f64::INFINITY, 2.0],
            convexity: Convexity::Psd,
        }
    }

    /// A stable-balancing-shaped spec: one weight variable per unit under the
    /// diagonal dispersion objective, a group-sum equality row, a two-sided
    /// moment-tolerance row, and a weight floor per unit.
    ///
    /// Minimize `0.5 sum w^2` subject to `sum w = 4`, `z' w` within `0.05` of zero
    /// for `z = (-1.5, -0.5, 0.5, 1)`, and `w >= 0.85`. The moment row binds at its
    /// lower side and the first unit's floor binds, giving
    /// `w = (0.85, 0.925, 1.075, 1.15)` and, from stationarity,
    /// `y = (-1, -0.15, -0.075, 0, 0, 0)`.
    fn balancing_shaped_spec() -> QpSpec {
        QpSpec {
            n: 4,
            m: 6,
            p: PMat::Diagonal(vec![1.0, 1.0, 1.0, 1.0]),
            q: vec![0.0; 4],
            // Rows: the group sum, the moment tolerance, then one floor per unit.
            a_indptr: vec![0, 3, 6, 9, 12],
            a_indices: vec![0, 1, 2, 0, 1, 3, 0, 1, 4, 0, 1, 5],
            a_values: vec![1.0, -1.5, 1.0, 1.0, -0.5, 1.0, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0],
            l: vec![4.0, -0.05, 0.85, 0.85, 0.85, 0.85],
            u: vec![
                4.0,
                0.05,
                f64::INFINITY,
                f64::INFINITY,
                f64::INFINITY,
                f64::INFINITY,
            ],
            convexity: Convexity::Psd,
        }
    }

    /// The largest absolute entry of the stationarity residual `P x + q + A' y`,
    /// which vanishes at an optimum exactly when the duals carry the orientation
    /// the `l <= A x <= u` convention assigns them.
    fn stationarity_residual(spec: &QpSpec, x: &[f64], y: &[f64]) -> f64 {
        let n = spec.n;
        let mut grad: Vec<f64> = match &spec.p {
            PMat::Diagonal(d) => d.iter().zip(x).map(|(di, xi)| di * xi).collect(),
            PMat::Dense(mat) => (0..n)
                .map(|i| (0..n).map(|j| mat[j * n + i] * x[j]).sum::<f64>())
                .collect(),
        };
        for (g, qi) in grad.iter_mut().zip(spec.q.iter()) {
            *g += qi;
        }
        for (col, g) in grad.iter_mut().enumerate() {
            for k in spec.a_indptr[col]..spec.a_indptr[col + 1] {
                *g += y[spec.a_indices[k]] * spec.a_values[k];
            }
        }
        grad.iter().fold(0.0_f64, |worst, g| worst.max(g.abs()))
    }

    /// A spec with no decision variables and `rows` constraint rows. The
    /// no-variable shape is what a stable-balancing solve assembles when its active
    /// set is empty: the balance rows survive and demand a group total of one from
    /// nothing at all.
    fn empty_spec(rows: usize) -> QpSpec {
        QpSpec {
            n: 0,
            m: rows,
            p: PMat::Diagonal(vec![]),
            q: vec![],
            a_indptr: vec![0],
            a_indices: vec![],
            a_values: vec![],
            l: vec![1.0; rows],
            u: vec![1.0; rows],
            convexity: Convexity::Psd,
        }
    }

    #[test]
    fn a_spec_with_no_decision_variables_is_refused_by_both_backends() {
        // The OSQP backend refuses a spec with no variables, and the methods
        // turn that refusal into a degenerate result carrying a failure status.
        // Clarabel's own dimension check accepts it, so without a guard the same
        // spec would reach the R layer as an infeasibility certificate under one
        // backend and a solver failure under the other.
        let spec = empty_spec(1);
        let theirs = Osqp
            .solve(&spec, &QpOptions::default(), &|| false)
            .unwrap_err();
        assert!(matches!(theirs, QpError::Setup(_)), "osqp: {theirs:?}");
        let mine = Clarabel
            .solve(&spec, &QpOptions::default(), &|| false)
            .unwrap_err();
        assert!(matches!(mine, QpError::Setup(_)), "clarabel: {mine:?}");
    }

    #[test]
    fn a_spec_with_no_variables_and_no_rows_is_refused_before_the_factorization() {
        // With neither variables nor rows the KKT system Clarabel factors is empty,
        // and its factorization indexes the first entry of that empty system. The
        // refusal has to come before the solver is built.
        let spec = empty_spec(0);
        let err = Clarabel
            .solve(&spec, &QpOptions::default(), &|| false)
            .unwrap_err();
        assert!(matches!(err, QpError::Setup(_)), "clarabel: {err:?}");
    }

    #[test]
    fn a_pending_interrupt_stops_the_solve() {
        let spec = simplex_spec(Convexity::Psd);
        let sol = Clarabel
            .solve(&spec, &QpOptions::default(), &|| true)
            .expect("psd spec solves");
        assert!(sol.interrupted);
        assert_eq!(sol.status, QpStatus::Interrupted);
    }

    #[test]
    fn an_interrupt_that_never_fires_leaves_the_solve_alone() {
        // The interrupt is polled once per iteration, so a closure that stays false
        // must not disturb the trajectory or the reported flag.
        use std::sync::atomic::{AtomicUsize, Ordering};
        let spec = hand_kkt_spec();
        let polls = AtomicUsize::new(0);
        let sol = Clarabel
            .solve(&spec, &QpOptions::default(), &|| {
                polls.fetch_add(1, Ordering::Relaxed);
                false
            })
            .expect("psd spec solves");
        assert!(!sol.interrupted);
        assert!(sol.status.is_solved(), "status {:?}", sol.status);
        assert!(
            polls.load(Ordering::Relaxed) > 0,
            "the interrupt was never polled"
        );
        assert!((sol.x[0] - 1.5).abs() < 1e-6, "x0 = {}", sol.x[0]);
    }

    #[test]
    fn a_numerical_failure_is_not_reported_as_the_iteration_cap() {
        // Advice to raise the iteration cap cannot help a solve that broke down
        // numerically or stalled, so neither may arrive at the R layer wearing the
        // status that produces that advice.
        assert_eq!(
            map_status(SolverStatus::NumericalError),
            QpStatus::NumericalError
        );
        assert_eq!(
            map_status(SolverStatus::InsufficientProgress),
            QpStatus::InsufficientProgress
        );
        assert_eq!(map_status(SolverStatus::Unsolved), QpStatus::NotSolved);
    }

    #[test]
    fn map_status_translates_every_clarabel_status() {
        assert_eq!(map_status(SolverStatus::Solved), QpStatus::Solved);
        assert_eq!(
            map_status(SolverStatus::AlmostSolved),
            QpStatus::SolvedInaccurate
        );
        assert_eq!(
            map_status(SolverStatus::PrimalInfeasible),
            QpStatus::PrimalInfeasible
        );
        assert_eq!(
            map_status(SolverStatus::AlmostPrimalInfeasible),
            QpStatus::PrimalInfeasible
        );
        assert_eq!(
            map_status(SolverStatus::DualInfeasible),
            QpStatus::DualInfeasible
        );
        assert_eq!(
            map_status(SolverStatus::AlmostDualInfeasible),
            QpStatus::DualInfeasible
        );
        assert_eq!(map_status(SolverStatus::MaxIterations), QpStatus::MaxIter);
        assert_eq!(map_status(SolverStatus::MaxTime), QpStatus::TimeLimit);
        assert_eq!(
            map_status(SolverStatus::CallbackTerminated),
            QpStatus::Interrupted
        );
    }

    #[test]
    fn settings_map_both_convergence_tolerances() {
        let opts = QpOptions {
            eps_abs: 1e-7,
            eps_rel: 1e-5,
            max_iter: 137,
            ..QpOptions::default()
        };
        let settings = solver_settings(&opts);
        assert_eq!(settings.tol_gap_abs, 1e-7);
        assert_eq!(settings.tol_feas, 1e-7);
        assert_eq!(settings.tol_gap_rel, 1e-5);
        assert_eq!(settings.max_iter, 137);
        assert!(!settings.verbose);
    }

    #[test]
    fn a_relative_tolerance_below_the_clarabel_default_is_carried_through() {
        // The relative gap tolerance sits in an OR-condition termination test, so
        // leaving it at its own default would let a request far below that default
        // terminate at the default instead.
        let opts = QpOptions {
            eps_abs: 1e-12,
            eps_rel: 1e-12,
            ..QpOptions::default()
        };
        let settings = solver_settings(&opts);
        assert!(
            settings.tol_gap_rel <= 1e-12,
            "tol_gap_rel is {}",
            settings.tol_gap_rel
        );
    }

    #[test]
    fn a_lower_active_row_carries_a_negative_dual() {
        let spec = hand_kkt_spec();
        let sol = Clarabel
            .solve(&spec, &QpOptions::default(), &|| false)
            .expect("psd spec solves");
        assert!(sol.status.is_solved(), "status {:?}", sol.status);
        let expected = [0.5, -1.0, 0.5];
        for (i, want) in expected.iter().enumerate() {
            assert!(
                (sol.duals[i] - want).abs() < 1e-6,
                "dual {i} is {} but the convention fixes it at {want}",
                sol.duals[i]
            );
        }
        let residual = stationarity_residual(&spec, &sol.x, &sol.duals);
        assert!(residual < 1e-6, "stationarity residual is {residual}");
    }

    #[test]
    fn duals_agree_with_osqp_on_a_balancing_shaped_spec() {
        let spec = balancing_shaped_spec();
        let mine = Clarabel
            .solve(&spec, &QpOptions::default(), &|| false)
            .expect("psd spec solves");
        let theirs = Osqp
            .solve(&spec, &QpOptions::default(), &|| false)
            .expect("setup succeeds");
        assert!(mine.status.is_solved(), "clarabel status {:?}", mine.status);
        assert!(theirs.status.is_solved(), "osqp status {:?}", theirs.status);
        for i in 0..spec.m {
            assert!(
                (mine.duals[i] - theirs.duals[i]).abs() < 1e-4,
                "dual {i}: clarabel {} against osqp {}",
                mine.duals[i],
                theirs.duals[i]
            );
        }
        let residual = stationarity_residual(&spec, &mine.x, &mine.duals);
        assert!(residual < 1e-5, "stationarity residual is {residual}");
    }
}
