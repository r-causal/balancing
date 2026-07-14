//! Parallel, deterministic pairwise Euclidean distance matrices.
//!
//! The distance matrix is the dominant cost of energy balancing. Each entry is
//! an independent function of two covariate rows, so the assembly is parallel by
//! construction and bit-for-bit identical regardless of thread count: no
//! reduction crosses the boundary of a single output entry. Work is split over
//! output columns and each column is filled sequentially.

use rayon::iter::{IndexedParallelIterator, ParallelIterator};
use rayon::slice::ParallelSliceMut;

use crate::threads::get_pool;

/// Euclidean distance between rows `i` and `j` of a column-major `n` by `p`
/// matrix.
#[inline]
fn row_distance(x: &[f64], n: usize, p: usize, i: usize, j: usize) -> f64 {
    let mut acc = 0.0;
    for c in 0..p {
        let diff = x[c * n + i] - x[c * n + j];
        acc += diff * diff;
    }
    acc.sqrt()
}

/// Dense symmetric Euclidean distance matrix of the rows of `x`.
///
/// `x` is a column-major `n` by `p` matrix; the result is a column-major `n` by
/// `n` matrix with a zero diagonal. Columns are filled in parallel inside a pool
/// of `threads` workers. The exact diagonal zero is written directly rather than
/// computed, so it never carries a rounding artifact.
pub fn euclidean(x: &[f64], n: usize, p: usize, threads: usize) -> Vec<f64> {
    let mut out = vec![0.0; n * n];
    if n == 0 {
        return out;
    }
    let pool = get_pool(threads);
    pool.install(|| {
        out.par_chunks_mut(n).enumerate().for_each(|(j, col)| {
            for (i, slot) in col.iter_mut().enumerate() {
                *slot = if i == j {
                    0.0
                } else {
                    row_distance(x, n, p, i, j)
                };
            }
        });
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_dimensional_distances_are_absolute_differences() {
        // Points 0, 1, 4 on the line: distances are the absolute differences.
        let x = [0.0, 1.0, 4.0];
        let d = euclidean(&x, 3, 1, 1);
        // Column-major 3 by 3.
        assert_eq!(d[0], 0.0); // (0,0)
        assert_eq!(d[1], 1.0); // (1,0)
        assert_eq!(d[2], 4.0); // (2,0)
        assert_eq!(d[3], 1.0); // (0,1)
        assert_eq!(d[5], 3.0); // (2,1)
        assert_eq!(d[8], 0.0); // (2,2)
    }

    #[test]
    fn two_dimensional_distance_is_the_hypotenuse() {
        // Rows (0,0) and (3,4): distance 5.
        let x = [0.0, 3.0, 0.0, 4.0]; // column-major: col1 = [0,3], col2 = [0,4]
        let d = euclidean(&x, 2, 2, 1);
        assert!((d[1] - 5.0).abs() < 1e-12);
        assert!((d[2] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn the_matrix_is_symmetric_with_a_zero_diagonal() {
        let x = [0.2, 1.7, -3.1, 4.4, 0.9, -2.2];
        let n = 3;
        let d = euclidean(&x, n, 2, 1);
        for i in 0..n {
            assert_eq!(d[i * n + i], 0.0);
            for j in 0..n {
                assert!((d[j * n + i] - d[i * n + j]).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn assembly_is_bit_identical_across_thread_counts() {
        // A moderately sized problem so the parallel path splits work, checked
        // bit-for-bit against the single-threaded result.
        let n = 200;
        let p = 4;
        let x: Vec<f64> = (0..n * p)
            .map(|k| ((k as f64) * 0.3).sin() * 2.0 + ((k as f64) * 0.07).cos())
            .collect();
        let one = euclidean(&x, n, p, 1);
        for threads in [2, 4, 8] {
            let many = euclidean(&x, n, p, threads);
            for (a, b) in one.iter().zip(many.iter()) {
                assert_eq!(a.to_bits(), b.to_bits());
            }
        }
    }
}
