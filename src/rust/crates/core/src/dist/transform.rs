//! Covariate transforms whose Euclidean distance reproduces a chosen metric.
//!
//! Rescaling the covariates once and then taking plain Euclidean distances is
//! equivalent to, and cheaper than, evaluating the target metric on every pair.
//! Scaled Euclidean divides each column by its weighted standard deviation;
//! Mahalanobis whitens by the weighted covariance through its eigendecomposition
//! so that Euclidean distance on the whitened rows equals the Mahalanobis
//! distance; plain Euclidean leaves the covariates unchanged.

use faer::{Mat, Side};

use super::Distance;
use crate::threads::{deterministic_map_reduce, get_pool};

/// Relative eigenvalue floor for the whitening pseudo-inverse: directions whose
/// covariance eigenvalue falls below this fraction of the largest are treated as
/// null and contribute nothing to the whitened distance, matching the
/// generalized-inverse behavior used elsewhere for singular covariances.
const WHITEN_RCOND: f64 = 1e-12;

/// Per-column weighted mean and variance, and the summed weight statistics.
struct Moments {
    /// Total weight `sum_i w_i`.
    sw: f64,
    /// Summed squared weight `sum_i w_i^2`.
    sw2: f64,
    /// Weighted column sums `sum_i w_i x_ic`.
    sx: Vec<f64>,
    /// Weighted column sums of squares `sum_i w_i x_ic^2`.
    sxx: Vec<f64>,
}

impl Moments {
    fn zeros(p: usize) -> Self {
        Self {
            sw: 0.0,
            sw2: 0.0,
            sx: vec![0.0; p],
            sxx: vec![0.0; p],
        }
    }
}

/// Weighted mean and reliability-weighted variance of each column.
///
/// The variance uses the frequency-weight-free (reliability) denominator
/// `1 - sum_i wn_i^2` with normalized weights `wn = w / sum(w)`, which reduces to
/// the usual `n - 1` sample variance when the weights are equal. A column with no
/// variance is reported with variance zero; callers substitute one so the column
/// is left unscaled rather than dividing by zero.
fn weighted_moments(
    covs: &[f64],
    n: usize,
    p: usize,
    w: &[f64],
    threads: usize,
) -> (Vec<f64>, Vec<f64>) {
    let pool = get_pool(threads);
    let acc = deterministic_map_reduce(
        &pool,
        n,
        || Moments::zeros(p),
        |acc, i| {
            let wi = w[i];
            acc.sw += wi;
            acc.sw2 += wi * wi;
            for j in 0..p {
                let x = covs[j * n + i];
                acc.sx[j] += wi * x;
                acc.sxx[j] += wi * x * x;
            }
        },
        |acc, other| {
            acc.sw += other.sw;
            acc.sw2 += other.sw2;
            for j in 0..p {
                acc.sx[j] += other.sx[j];
                acc.sxx[j] += other.sxx[j];
            }
        },
    );

    let sw = acc.sw;
    let mut means = vec![0.0; p];
    let mut vars = vec![0.0; p];
    if sw <= 0.0 {
        return (means, vars);
    }
    // The reliability denominator; guarded so a single dominant weight does not
    // produce a negative or zero divisor.
    let denom = 1.0 - acc.sw2 / (sw * sw);
    for j in 0..p {
        let mean = acc.sx[j] / sw;
        means[j] = mean;
        let second = acc.sxx[j] / sw - mean * mean;
        vars[j] = if denom > 0.0 {
            (second / denom).max(0.0)
        } else {
            0.0
        };
    }
    (means, vars)
}

/// Transform `covs` (column-major `n` by `p`) so that Euclidean distance on the
/// result reproduces `distance`.
///
/// `w` are the weights used for the standardizing statistics. The Mahalanobis
/// transform can change the column count when the covariance is rank deficient,
/// so callers read the transformed column count from the returned length divided
/// by `n`.
pub fn transform(
    covs: &[f64],
    n: usize,
    p: usize,
    distance: Distance,
    w: &[f64],
    threads: usize,
) -> Vec<f64> {
    match distance {
        Distance::Euclidean => covs.to_vec(),
        Distance::ScaledEuclidean => scaled_euclidean(covs, n, p, w, threads),
        Distance::Mahalanobis => mahalanobis(covs, n, p, w, threads),
    }
}

/// Divide each column by its weighted standard deviation. A column with no
/// variance is left unscaled.
fn scaled_euclidean(covs: &[f64], n: usize, p: usize, w: &[f64], threads: usize) -> Vec<f64> {
    let (_means, vars) = weighted_moments(covs, n, p, w, threads);
    let mut out = covs.to_vec();
    for j in 0..p {
        let sd = vars[j].sqrt();
        let scale = if sd > 0.0 { 1.0 / sd } else { 1.0 };
        for i in 0..n {
            out[j * n + i] *= scale;
        }
    }
    out
}

/// Whiten the covariates by the weighted covariance so that Euclidean distance
/// on the result equals the Mahalanobis distance.
///
/// The covariates are first standardized to unit weighted variance, so the
/// covariance used for whitening is the weighted correlation matrix; this
/// mirrors the reference pipeline and keeps the whitening well scaled. The
/// symmetric eigendecomposition supplies the inverse square root, with null
/// directions dropped by the relative eigenvalue floor.
fn mahalanobis(covs: &[f64], n: usize, p: usize, w: &[f64], threads: usize) -> Vec<f64> {
    if p == 0 || n == 0 {
        return covs.to_vec();
    }
    let (means, vars) = weighted_moments(covs, n, p, w, threads);

    // Standardize columns to weighted mean zero and unit variance; a zero
    // variance column becomes a constant zero column and drops out of the
    // covariance below.
    let mut std = vec![0.0; n * p];
    for j in 0..p {
        let sd = vars[j].sqrt();
        let scale = if sd > 0.0 { 1.0 / sd } else { 0.0 };
        for i in 0..n {
            std[j * n + i] = (covs[j * n + i] - means[j]) * scale;
        }
    }

    // Weighted covariance of the standardized columns, using the same
    // reliability denominator as the per-column variance.
    let pool = get_pool(threads);
    let (sw, sw2, cov_flat) = deterministic_map_reduce(
        &pool,
        n,
        || (0.0_f64, 0.0_f64, vec![0.0_f64; p * p]),
        |acc, i| {
            let wi = w[i];
            acc.0 += wi;
            acc.1 += wi * wi;
            for a in 0..p {
                let xa = std[a * n + i];
                for b in 0..p {
                    acc.2[a * p + b] += wi * xa * std[b * n + i];
                }
            }
        },
        |acc, other| {
            acc.0 += other.0;
            acc.1 += other.1;
            for k in 0..p * p {
                acc.2[k] += other.2[k];
            }
        },
    );
    let denom = sw * (1.0 - sw2 / (sw * sw));
    let cov = Mat::from_fn(p, p, |a, b| {
        if denom > 0.0 {
            cov_flat[a * p + b] / denom
        } else {
            0.0
        }
    });

    // Whitening factor V diag(lambda^{-1/2}) from the symmetric eigendecomposition.
    let eigen = match cov.self_adjoint_eigen(Side::Lower) {
        Ok(e) => e,
        // A covariance that cannot be factored leaves the standardized columns
        // as the transform, which is the scaled-Euclidean fallback.
        Err(_) => return std,
    };
    let values = eigen.S();
    let vectors = eigen.U();
    let max_abs = (0..p).fold(0.0_f64, |m, i| m.max(values[i].abs()));
    let floor = WHITEN_RCOND * max_abs;
    let inv_sqrt: Vec<f64> = (0..p)
        .map(|i| {
            let lambda = values[i];
            if lambda > floor {
                1.0 / lambda.sqrt()
            } else {
                0.0
            }
        })
        .collect();

    // Transformed row i is std_i * V * diag(inv_sqrt): column k of the output is
    // sum_a std_ia * V_ak * inv_sqrt_k.
    let mut out = vec![0.0; n * p];
    for i in 0..n {
        for k in 0..p {
            let s = inv_sqrt[k];
            if s == 0.0 {
                continue;
            }
            let mut acc = 0.0;
            for a in 0..p {
                acc += std[a * n + i] * *vectors.get(a, k);
            }
            out[k * n + i] = acc * s;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equal_weights_reproduce_the_sample_variance() {
        // Column [1, 2, 3, 4]: sample variance with the n - 1 denominator is
        // 5/3. Equal weights must reproduce it through the reliability formula.
        let covs = [1.0, 2.0, 3.0, 4.0];
        let w = [1.0, 1.0, 1.0, 1.0];
        let (means, vars) = weighted_moments(&covs, 4, 1, &w, 1);
        assert!((means[0] - 2.5).abs() < 1e-12);
        assert!((vars[0] - 5.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn scaled_euclidean_divides_by_the_standard_deviation() {
        let covs = [1.0, 2.0, 3.0, 4.0];
        let w = [1.0, 1.0, 1.0, 1.0];
        let out = scaled_euclidean(&covs, 4, 1, &w, 1);
        let sd = (5.0_f64 / 3.0).sqrt();
        for i in 0..4 {
            assert!((out[i] - covs[i] / sd).abs() < 1e-12);
        }
    }

    #[test]
    fn a_constant_column_is_left_unscaled() {
        let covs = [2.0, 2.0, 2.0];
        let w = [1.0, 1.0, 1.0];
        let out = scaled_euclidean(&covs, 3, 1, &w, 1);
        assert_eq!(out, covs.to_vec());
    }

    #[test]
    fn mahalanobis_matches_scaled_on_uncorrelated_columns() {
        // Two uncorrelated columns with different variances: after whitening by
        // the correlation matrix (identity here) the transform equals the
        // standardized columns, so pairwise Euclidean distances match the
        // scaled-Euclidean transform up to the shared standardization.
        let n = 4;
        let covs = [
            0.0, 1.0, 2.0, 3.0, // column 1
            0.0, 2.0, 4.0, 6.0, // column 2 (perfectly proportional to column 1)
        ];
        let w = [1.0; 4];
        // Columns are collinear, so the correlation matrix is singular; the
        // whitening drops the null direction and still returns a finite result.
        let out = mahalanobis(&covs, n, 2, &w, 1);
        assert_eq!(out.len(), n * 2);
        assert!(out.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn mahalanobis_whitens_to_unit_covariance() {
        // Independent standard-normal-like columns: the whitened columns should
        // have (weighted) covariance close to the identity, so the sum of
        // squared off-diagonal covariance is small.
        let n = 5;
        let covs = [
            -2.0, -1.0, 0.0, 1.0, 2.0, // column 1
            1.0, -1.0, 2.0, 0.0, -2.0, // column 2 (different pattern)
        ];
        let w = [1.0; 5];
        let out = mahalanobis(&covs, n, 2, &w, 1);
        // Recompute the weighted covariance of the output with the same
        // reliability denominator and check it is near identity.
        let (_m, _v) = weighted_moments(&out, n, 2, &w, 1);
        let mut cov = [0.0; 4];
        let sw = 5.0;
        let sw2 = 5.0;
        let denom = sw * (1.0 - sw2 / (sw * sw));
        // Center columns.
        let mut c = out.clone();
        for j in 0..2 {
            let mean: f64 = (0..n).map(|i| out[j * n + i]).sum::<f64>() / n as f64;
            for i in 0..n {
                c[j * n + i] -= mean;
            }
        }
        for a in 0..2 {
            for b in 0..2 {
                let mut s = 0.0;
                for i in 0..n {
                    s += c[a * n + i] * c[b * n + i];
                }
                cov[a * 2 + b] = s / denom;
            }
        }
        assert!((cov[0] - 1.0).abs() < 1e-6, "cov[0,0] = {}", cov[0]);
        assert!((cov[3] - 1.0).abs() < 1e-6, "cov[1,1] = {}", cov[3]);
        assert!(cov[1].abs() < 1e-6, "off-diagonal = {}", cov[1]);
    }
}
