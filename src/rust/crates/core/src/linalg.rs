//! Dense linear algebra adapters over column-major slices.
//!
//! Inputs cross the R boundary as column-major `f64` slices; [`col_major`] wraps
//! one as a faer [`MatRef`] with no copy. The solver family also needs a small
//! symmetric solve for the Newton step, provided by [`solve_symmetric_ridge`].
//! Only the small parameter-by-parameter systems run through faer, and always
//! sequentially: the data-parallel work lives in the deterministic reductions.

use faer::linalg::solvers::Solve;
use faer::{Mat, MatRef, Side};

/// View a column-major slice as a faer matrix without copying.
///
/// The slice must hold exactly `nrows * ncols` elements in column-major order,
/// which is R's storage convention.
pub fn col_major(data: &[f64], nrows: usize, ncols: usize) -> MatRef<'_, f64> {
    MatRef::from_column_major_slice(data, nrows, ncols)
}

/// Solve `(H + ridge * I) x = rhs` for a symmetric `H`, writing `x` back into
/// `rhs`. `H` is a column-major `p * p` slice; only its lower triangle is read.
///
/// Returns `false` when the ridge-shifted matrix is not positive definite and
/// the LDLT factorization fails, leaving `rhs` untouched. Callers escalate the
/// ridge and retry.
pub fn solve_symmetric_ridge(h: &[f64], p: usize, ridge: f64, rhs: &mut [f64]) -> bool {
    debug_assert_eq!(h.len(), p * p);
    debug_assert_eq!(rhs.len(), p);
    let a = Mat::from_fn(p, p, |i, j| {
        let mut value = h[j * p + i];
        if i == j {
            value += ridge;
        }
        value
    });
    let factor = match a.ldlt(Side::Lower) {
        Ok(factor) => factor,
        Err(_) => return false,
    };
    let mut b = Mat::from_fn(p, 1, |i, _| rhs[i]);
    factor.solve_in_place(b.as_mut());
    let solution = b.col_as_slice(0);
    rhs.copy_from_slice(solution);
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn col_major_reads_column_order() {
        // Two columns [1,2,3] and [4,5,6].
        let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let m = col_major(&data, 3, 2);
        assert_eq!(*m.get(0, 0), 1.0);
        assert_eq!(*m.get(2, 0), 3.0);
        assert_eq!(*m.get(0, 1), 4.0);
        assert_eq!(*m.get(2, 1), 6.0);
    }

    #[test]
    fn solves_a_small_spd_system() {
        // H = [[2, 0], [0, 4]], rhs = [2, 8] -> x = [1, 2].
        let h = [2.0, 0.0, 0.0, 4.0];
        let mut rhs = [2.0, 8.0];
        assert!(solve_symmetric_ridge(&h, 2, 0.0, &mut rhs));
        assert!((rhs[0] - 1.0).abs() < 1e-12);
        assert!((rhs[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn ridge_shifts_the_diagonal() {
        // (H + 1*I) with H = [[1,0],[0,1]] gives 2 on the diagonal.
        let h = [1.0, 0.0, 0.0, 1.0];
        let mut rhs = [2.0, 2.0];
        assert!(solve_symmetric_ridge(&h, 2, 1.0, &mut rhs));
        assert!((rhs[0] - 1.0).abs() < 1e-12);
        assert!((rhs[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn reports_failure_on_zero_pivot() {
        // A zero leading pivot has no unpivoted LDLT, so the solve reports
        // failure and the caller escalates the ridge.
        let h = [0.0, 1.0, 1.0, 0.0];
        let mut rhs = [1.0, 1.0];
        assert!(!solve_symmetric_ridge(&h, 2, 0.0, &mut rhs));
    }

    #[test]
    fn a_positive_ridge_rescues_a_zero_pivot() {
        // Shifting the same matrix by a ridge restores a usable factorization.
        let h = [0.0, 1.0, 1.0, 0.0];
        let mut rhs = [1.0, 1.0];
        assert!(solve_symmetric_ridge(&h, 2, 4.0, &mut rhs));
    }
}
