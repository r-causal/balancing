//! The Clarabel backend: an interior-point solver for convex problems.
//!
//! Interior-point methods require a positive-semidefinite quadratic form, so
//! this backend refuses a spec tagged `Indefinite` up front rather than
//! returning a meaningless iterate. It is admissible for the diagonal and
//! kernel forms of the other quadratic-program methods; energy balancing always
//! routes to OSQP. Two-sided bounds are expressed in the conic form Clarabel
//! consumes: an equality row becomes a zero-cone row, and a finite lower or
//! upper bound becomes a nonnegative-cone row, with the lower bound negated so
//! `A x >= l` reads as `(-A) x <= -l`.

use clarabel::algebra::CscMatrix;
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
/// coefficients over the `n` variables, the right-hand side, and which original
/// constraint it came from so the dual can be mapped back.
struct ConicRow {
    coeffs: Vec<(usize, f64)>,
    rhs: f64,
    source: usize,
}

impl QpBackend for Clarabel {
    fn name(&self) -> &'static str {
        "clarabel"
    }

    fn solve(
        &self,
        spec: &QpSpec,
        opts: &QpOptions,
        _interrupt: &dyn Fn() -> bool,
    ) -> Result<QpSolution, QpError> {
        if spec.convexity == Convexity::Indefinite {
            return Err(QpError::Indefinite);
        }

        let n = spec.n;
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
                });
                continue;
            }
            if u.is_finite() {
                zero_rows_push_upper(&mut nonneg_rows, coeffs, u, c);
            }
            if l.is_finite() {
                let negated: Vec<(usize, f64)> = coeffs.iter().map(|&(i, v)| (i, -v)).collect();
                nonneg_rows.push(ConicRow {
                    coeffs: negated,
                    rhs: -l,
                    source: c,
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

        let settings = DefaultSettings::<f64> {
            verbose: false,
            max_iter: opts.max_iter.min(u32::MAX as usize) as u32,
            tol_gap_abs: opts.eps_abs,
            tol_feas: opts.eps_abs,
            ..Default::default()
        };

        let mut solver = DefaultSolver::new(&p_csc, &spec.q, &a_csc, &b, &cones, settings)
            .map_err(|e| QpError::Setup(format!("{e:?}")))?;
        solver.solve();

        let x = solver.solution.x.clone();
        // Map the conic duals back to one multiplier per original constraint by
        // combining the upper and (negated) lower contributions.
        let mut duals = vec![0.0; spec.m];
        for (row, z) in all_rows.iter().zip(solver.solution.z.iter()) {
            duals[row.source] += *z;
        }
        let obj = objective(spec, &x);
        Ok(QpSolution {
            x,
            duals,
            status: map_status(solver.solution.status),
            iterations: solver.solution.iterations as usize,
            obj,
            pri_res: solver.info.res_primal,
            dua_res: solver.info.res_dual,
            interrupted: false,
        })
    }
}

/// Append the upper-bound inequality row `A x <= u` for constraint `c`.
fn zero_rows_push_upper(rows: &mut Vec<ConicRow>, coeffs: &[(usize, f64)], u: f64, c: usize) {
    rows.push(ConicRow {
        coeffs: coeffs.to_vec(),
        rhs: u,
        source: c,
    });
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
        SolverStatus::MaxTime => QpStatus::TimeLimit,
        _ => QpStatus::MaxIter,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qp::PMat;

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
}
