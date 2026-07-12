//! Dense linear algebra adapters over column-major slices.
//!
//! Inputs cross the R boundary as column-major `f64` slices; [`col_major`] wraps
//! one as a faer [`MatRef`] with no copy. The solver family also needs a small
//! symmetric solve for the Newton step, provided by [`solve_symmetric_ridge`].
//! Only the small parameter-by-parameter systems run through faer, and always
//! sequentially: the data-parallel work lives in the deterministic reductions.

use faer::linalg::solvers::Solve;
use faer::{Mat, MatRef, Side};

/// Relative eigenvalue floor for the symmetric pseudo-inverse: directions whose
/// eigenvalue falls below this fraction of the largest are treated as null and
/// dropped, which is the standard generalized-inverse tolerance behavior.
const PINV_RCOND: f64 = 1e-12;

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

/// Moore-Penrose pseudo-inverse of a symmetric positive-semidefinite matrix.
///
/// `a` is a column-major `p * p` slice holding a symmetric matrix; only its
/// lower triangle is relied on through the self-adjoint eigendecomposition. The
/// result is the column-major `p * p` pseudo-inverse, formed as
/// `U diag(1 / lambda) U'` over the eigenvalues above a relative floor and zero
/// on the null directions. This is the two-step GMM weighting matrix inverse,
/// where the moment covariance can be singular when constraints are redundant.
///
/// A singular covariance is handled by the null-direction thresholding and still
/// returns `Some`. `None` signals a failed eigendecomposition, which is not
/// expected for a finite covariance matrix; the caller surfaces it rather than
/// degenerating the criterion to a silent zero weighting.
pub fn pseudo_inverse_symmetric(a: &[f64], p: usize) -> Option<Vec<f64>> {
    debug_assert_eq!(a.len(), p * p);
    if p == 0 {
        return Some(Vec::new());
    }
    let mat = Mat::from_fn(p, p, |i, j| a[j * p + i]);
    let eigen = mat.self_adjoint_eigen(Side::Lower).ok()?;
    let values = eigen.S();
    let vectors = eigen.U();

    let max_abs = (0..p).fold(0.0_f64, |m, i| m.max(values[i].abs()));
    let floor = PINV_RCOND * max_abs;

    // inv_lambda holds the reciprocal eigenvalues on the kept directions and zero
    // elsewhere, so the reconstruction below sums only the retained rank.
    let inv_lambda: Vec<f64> = (0..p)
        .map(|i| {
            let lambda = values[i];
            if lambda > floor { 1.0 / lambda } else { 0.0 }
        })
        .collect();

    let mut out = vec![0.0; p * p];
    for col in 0..p {
        for row in 0..p {
            let mut acc = 0.0;
            for (k, &il) in inv_lambda.iter().enumerate() {
                acc += *vectors.get(row, k) * il * *vectors.get(col, k);
            }
            out[col * p + row] = acc;
        }
    }
    Some(out)
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

    // Multiply two column-major p*p matrices for the pseudo-inverse checks.
    fn matmul(a: &[f64], b: &[f64], p: usize) -> Vec<f64> {
        let mut out = vec![0.0; p * p];
        for col in 0..p {
            for row in 0..p {
                let mut acc = 0.0;
                for k in 0..p {
                    acc += a[k * p + row] * b[col * p + k];
                }
                out[col * p + row] = acc;
            }
        }
        out
    }

    #[test]
    fn pseudo_inverse_is_the_inverse_for_a_full_rank_matrix() {
        // A well-conditioned symmetric positive-definite matrix: its pseudo-
        // inverse is the ordinary inverse, so the product is the identity.
        let a = [4.0, 1.0, 1.0, 3.0];
        let pinv = pseudo_inverse_symmetric(&a, 2).expect("decomposition succeeds");
        let prod = matmul(&a, &pinv, 2);
        let identity = [1.0, 0.0, 0.0, 1.0];
        for i in 0..4 {
            assert!(
                (prod[i] - identity[i]).abs() < 1e-12,
                "entry {i} = {}",
                prod[i]
            );
        }
    }

    #[test]
    fn pseudo_inverse_drops_the_null_direction() {
        // A rank-one matrix v v' with v = (1, 1): the pseudo-inverse is
        // v v' / (v'v)^2 = 0.25 * ones, and A A^+ A = A holds.
        let a = [1.0, 1.0, 1.0, 1.0];
        let pinv = pseudo_inverse_symmetric(&a, 2).expect("decomposition succeeds");
        for value in pinv {
            assert!((value - 0.25).abs() < 1e-12, "entry {value}");
        }
        // The defining relation A A^+ A = A confirms the generalized inverse.
        let a = [1.0, 1.0, 1.0, 1.0];
        let pinv = pseudo_inverse_symmetric(&a, 2).expect("decomposition succeeds");
        let recon = matmul(&matmul(&a, &pinv, 2), &a, 2);
        for i in 0..4 {
            assert!((recon[i] - a[i]).abs() < 1e-12);
        }
    }
}
