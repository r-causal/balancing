//! Covariate transforms whose Euclidean distance reproduces a chosen metric.
//!
//! Rescaling the covariates once and then taking plain Euclidean distances is
//! equivalent to, and cheaper than, evaluating the target metric on every pair.
//! Scaled Euclidean centers each column at its weighted mean and divides by its
//! weighted standard deviation; Mahalanobis centers as well and whitens by the
//! weighted covariance through its eigendecomposition so that Euclidean distance
//! on the whitened rows equals the Mahalanobis distance; plain Euclidean leaves
//! the covariates unchanged.
//!
//! Centering is invisible to the metric, since a common shift of every row
//! cancels in every pairwise difference, and it is what keeps the transformed
//! values at the spread's own scale rather than at the column's offset.

use faer::{Mat, Side};

use super::Distance;
use crate::stats::reliability_variance;
use crate::threads::{deterministic_map_reduce, get_pool};

/// Relative eigenvalue floor for the whitening pseudo-inverse: directions whose
/// covariance eigenvalue falls below this fraction of the largest are treated as
/// null and contribute nothing to the whitened distance, matching the
/// generalized-inverse behavior used elsewhere for singular covariances.
const WHITEN_RCOND: f64 = 1e-12;

/// Weight totals and per-column weighted sums, the first pass of the moments.
struct ColumnSums {
    /// Total weight `sum_i w_i`.
    sw: f64,
    /// Summed squared weight `sum_i w_i^2`.
    sw2: f64,
    /// Weighted column sums `sum_i w_i x_ic`.
    sx: Vec<f64>,
}

impl ColumnSums {
    fn zeros(p: usize) -> Self {
        Self {
            sw: 0.0,
            sw2: 0.0,
            sx: vec![0.0; p],
        }
    }
}

/// Weighted mean and reliability-weighted variance of each column.
///
/// The variance is accumulated in two scans of the covariates, the first
/// producing the weight totals and the column means and the second the weighted
/// squared deviations from those means. Forming the deviations explicitly is
/// what keeps the result accurate on a column whose mean is large against its
/// spread; the `stats` module records why that case is the ordinary one
/// rather than a corner. The alternative single-scan form, Welford's weighted online
/// update, would save the second read of the data but needs a division by the
/// running total inside the inner loop and a guarded pairwise merge to combine
/// the chunk accumulators, and the transform is `O(n p)` in front of the
/// `O(n^2 p)` pairwise work it feeds, so the second scan is not a cost worth
/// that complexity. Both scans stay column-major and both reduce through
/// [`deterministic_map_reduce`], so the result is still independent of the
/// thread count.
///
/// A column with no variance is reported with variance zero; callers substitute
/// one so the column is left unscaled rather than dividing by zero.
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
        || ColumnSums::zeros(p),
        |acc, i| {
            let wi = w[i];
            acc.sw += wi;
            acc.sw2 += wi * wi;
            for j in 0..p {
                acc.sx[j] += wi * covs[j * n + i];
            }
        },
        |acc, other| {
            acc.sw += other.sw;
            acc.sw2 += other.sw2;
            for j in 0..p {
                acc.sx[j] += other.sx[j];
            }
        },
    );

    let sw = acc.sw;
    if sw <= 0.0 {
        return (vec![0.0; p], vec![0.0; p]);
    }
    let means: Vec<f64> = acc.sx.iter().map(|&sxj| sxj / sw).collect();

    let ss = deterministic_map_reduce(
        &pool,
        n,
        || vec![0.0_f64; p],
        |acc, i| {
            let wi = w[i];
            for j in 0..p {
                let d = covs[j * n + i] - means[j];
                acc[j] += wi * d * d;
            }
        },
        |acc, other| {
            for j in 0..p {
                acc[j] += other[j];
            }
        },
    );

    let vars = ss
        .iter()
        .map(|&ssj| reliability_variance(ssj, sw, acc.sw2))
        .collect();
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

/// Subtract each column's weighted mean and divide by its weighted standard
/// deviation. A column with no variance is centered but left unscaled, so it
/// comes back as exactly zero.
///
/// Subtracting the mean changes no pairwise Euclidean distance: it shifts every
/// row by the same vector, which cancels in every difference. What it buys is
/// the range the output occupies. A date-time column, which `as.numeric()` puts
/// near 1.7e9, with an hour of spread divides to about 8e5 with deviations of
/// order one, and a unit in the last place there is about 1.2e-10, so the
/// distances built from it are quantized ten digits coarser than the doubles
/// carrying them. Centering leaves the output at the deviations' own scale and
/// costs one subtraction per entry, on a pass that already reads and writes
/// every entry. The Mahalanobis transform centers for the same reason.
fn scaled_euclidean(covs: &[f64], n: usize, p: usize, w: &[f64], threads: usize) -> Vec<f64> {
    let (means, vars) = weighted_moments(covs, n, p, w, threads);
    let mut out = vec![0.0; n * p];
    for j in 0..p {
        let sd = vars[j].sqrt();
        let scale = if sd > 0.0 { 1.0 / sd } else { 1.0 };
        for i in 0..n {
            out[j * n + i] = (covs[j * n + i] - means[j]) * scale;
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
        //
        // This arm is believed unreachable: the matrix is symmetric by
        // construction and the self-adjoint eigendecomposition of a finite
        // symmetric matrix converges, singular or not, which is why the rank
        // deficiency is handled by the eigenvalue floor below rather than here.
        // Non-finite input poisons the standardized columns as thoroughly as it
        // poisons the decomposition, so the fallback is no worse than the path it
        // replaces, and a distance definition is not worth a panic. Silently
        // changing the distance is the lesser cost of the two.
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
    use crate::dist::pairwise;
    use crate::threads::REDUCE_CHUNK;

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
    fn scaled_euclidean_centers_and_divides_by_the_standard_deviation() {
        let covs = [1.0, 2.0, 3.0, 4.0];
        let w = [1.0, 1.0, 1.0, 1.0];
        let out = scaled_euclidean(&covs, 4, 1, &w, 1);
        let sd = (5.0_f64 / 3.0).sqrt();
        for i in 0..4 {
            assert!((out[i] - (covs[i] - 2.5) / sd).abs() < 1e-12);
        }
    }

    #[test]
    fn a_constant_column_centers_to_zero_and_is_left_unscaled() {
        // Centering sends a column with no spread onto exactly zero, which is
        // also where the Mahalanobis standardization puts it. The scale stays at
        // one rather than dividing by a zero standard deviation.
        let covs = [2.0, 2.0, 2.0];
        let w = [1.0, 1.0, 1.0];
        let out = scaled_euclidean(&covs, 3, 1, &w, 1);
        assert_eq!(out, vec![0.0, 0.0, 0.0]);
    }

    #[test]
    fn mahalanobis_matches_scaled_on_uncorrelated_columns() {
        // Three deliberately orthogonal columns: each has mean zero and each pair
        // has zero cross product, so the weighted correlation matrix the whitening
        // reads is the identity and the whitening can only rotate the standardized
        // columns. A rotation leaves every pairwise Euclidean distance where it
        // was, so the transform's distances must equal the distances on the
        // per-column standardized data, which is what the Mahalanobis metric
        // reduces to when the covariance is diagonal.
        //
        // The standardizing scales are written out from the reliability variance
        // rather than taken from the transform: five equal weights give the
        // denominator 1 - 5/25, so a mean-zero column's variance is its mean
        // square over 0.8, and the mean squares here are 2, 2.8, and 2.
        let n = 5;
        let p = 3;
        let covs = [
            -2.0, -1.0, 0.0, 1.0, 2.0, // column 1
            2.0, -1.0, -2.0, -1.0, 2.0, // column 2
            -1.0, 2.0, 0.0, -2.0, 1.0, // column 3
        ];
        let w = [1.0; 5];
        let sds = [
            (2.0_f64 / 0.8).sqrt(),
            (2.8_f64 / 0.8).sqrt(),
            (2.0_f64 / 0.8).sqrt(),
        ];
        let mut standardized = vec![0.0; n * p];
        for j in 0..p {
            for i in 0..n {
                standardized[j * n + i] = covs[j * n + i] / sds[j];
            }
        }

        let distance = |x: &[f64], i: usize, k: usize| {
            (0..p)
                .map(|j| {
                    let d = x[j * n + i] - x[j * n + k];
                    d * d
                })
                .sum::<f64>()
                .sqrt()
        };

        // Standardizing has to be a real change on this fixture, or leaving the
        // covariates alone would satisfy the equality below.
        assert!(
            (distance(&covs, 0, 4) - distance(&standardized, 0, 4)).abs() > 0.5,
            "the fixture must standardize to something other than itself"
        );

        let out = mahalanobis(&covs, n, p, &w, 1);
        assert_eq!(out.len(), n * p);
        for i in 0..n {
            for k in (i + 1)..n {
                let got = distance(&out, i, k);
                let want = distance(&standardized, i, k);
                assert!(
                    (got - want).abs() < 1e-12,
                    "pair ({i}, {k}): whitened distance {got}, standardized distance {want}"
                );
            }
        }
    }

    #[test]
    fn mahalanobis_drops_the_null_direction_of_a_singular_covariance() {
        // Columns 1 and 2 are exactly collinear, so the correlation matrix the
        // whitening reads is singular and one eigenvalue is numerically zero.
        // The relative floor turns that direction's inverse square root into a
        // hard zero rather than clamping it at the floor, which is the
        // generalized-inverse reading: the null direction contributes nothing.
        // Whitened distance is then the Mahalanobis distance of the rank-
        // deficient data under the pseudo-inverse, and that equals the whitened
        // distance on any basis of the column space, here columns 1 and 3.
        //
        // The expectation is written out rather than taken from a second call.
        // Standardizing sends column 2 onto column 1 exactly, so the
        // correlation matrix is [[1, 1, 0], [1, 1, 0], [0, 0, 1]], whose
        // pseudo-inverse is [[1/4, 1/4, 0], [1/4, 1/4, 0], [0, 0, 1]]. A
        // difference of standardized rows is (da, da, dc), and the quadratic
        // form against that pseudo-inverse collapses to da^2 + dc^2: the plain
        // Euclidean distance on the two standardized basis columns. Their
        // scales come from the reliability variance, as in the test above: five
        // equal weights give the denominator 1 - 5/25, and both basis columns
        // are mean zero with mean square 2.
        let n = 5;
        let p = 3;
        let covs = [
            -2.0, -1.0, 0.0, 1.0, 2.0, // column 1
            -4.0, -2.0, 0.0, 2.0, 4.0, // column 2, twice column 1
            -1.0, 2.0, 0.0, -2.0, 1.0, // column 3, orthogonal to column 1
        ];
        let w = [1.0; 5];
        let sd = (2.0_f64 / 0.8).sqrt();
        let basis_distance = |i: usize, k: usize| {
            let da = (covs[i] - covs[k]) / sd;
            let dc = (covs[2 * n + i] - covs[2 * n + k]) / sd;
            (da * da + dc * dc).sqrt()
        };

        let out = mahalanobis(&covs, n, p, &w, 1);
        assert_eq!(out.len(), n * p);

        // Exactly one output column is identically zero: the floor drops the
        // null direction rather than clamping its eigenvalue at the floor and
        // keeping the direction. The eigenvalue here is a hard zero, so without
        // the floor its inverse square root is infinite and the direction comes
        // back as a column of NaN rather than as nothing at all.
        let zero_columns = (0..p)
            .filter(|k| (0..n).all(|i| out[k * n + i] == 0.0))
            .count();
        assert_eq!(zero_columns, 1, "the null direction must be dropped");

        let distance = |x: &[f64], i: usize, k: usize| {
            (0..p)
                .map(|j| {
                    let d = x[j * n + i] - x[j * n + k];
                    d * d
                })
                .sum::<f64>()
                .sqrt()
        };

        // Carrying the duplicated column as its own coordinate, which is what
        // the scaled-Euclidean fallback on the factorization failure does,
        // counts the collinear part twice inside the square root. The fixture
        // has to separate that reading from the pseudo-inverse one, or the
        // equality below would hold either way.
        let mut standardized = vec![0.0; n * p];
        for j in 0..p {
            let scale = if j == 1 { 2.0 * sd } else { sd };
            for i in 0..n {
                standardized[j * n + i] = covs[j * n + i] / scale;
            }
        }
        assert!(
            (distance(&standardized, 0, 4) - basis_distance(0, 4)).abs() > 0.5,
            "the fixture must separate the dropped direction from a kept one"
        );

        for i in 0..n {
            for k in (i + 1)..n {
                let got = distance(&out, i, k);
                let want = basis_distance(i, k);
                assert!(
                    (got - want).abs() < 1e-12,
                    "pair ({i}, {k}): whitened distance {got}, basis distance {want}"
                );
            }
        }
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

    /// A deterministic stream on the unit interval from a linear congruential
    /// generator. The offset fixtures below need many values with no structure,
    /// and generating them here keeps the numbers identical on every platform
    /// without reaching for a random-number dependency.
    fn lcg_unit(n: usize, seed: u64) -> Vec<f64> {
        let mut state = seed;
        (0..n)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                (state >> 11) as f64 / 9_007_199_254_740_992.0
            })
            .collect()
    }

    /// The offset the fixtures below sit at: the magnitude `as.numeric()` gives a
    /// POSIXct, with an hour of spread around it.
    const DATE_TIME_OFFSET: f64 = 1.7e9;

    /// `n` rows of two columns and their non-uniform weights. Column zero sits at
    /// date-time scale, `offset + spread * z` with the offset above and the
    /// spread one hour, which is the shape that makes a one-pass variance
    /// cancel: the mean is roughly half a million standard deviations from zero.
    /// Column one is an ordinary unit-scale column carried alongside so the same
    /// fixture covers the undemanding case.
    fn offset_fixture_n(n: usize) -> (Vec<f64>, Vec<f64>) {
        let spread = 3600.0;
        let z = lcg_unit(n, 20_250_828);
        let plain = lcg_unit(n, 12_345_677);
        let weight_draws = lcg_unit(n, 987_654_321);
        let mut covs = Vec::with_capacity(2 * n);
        covs.extend(
            z.iter()
                .map(|&zi| DATE_TIME_OFFSET + spread * (2.0 * zi - 1.0)),
        );
        covs.extend(plain.iter().map(|&pi| 2.0 * pi - 1.0));
        let w = weight_draws.iter().map(|&vi| 0.25 + 1.5 * vi).collect();
        (covs, w)
    }

    /// The fifty-row instance of that fixture, with its row count and offset, for
    /// the tests that center by the offset exactly.
    fn offset_fixture() -> (Vec<f64>, Vec<f64>, usize, f64) {
        let n = 50;
        let (covs, w) = offset_fixture_n(n);
        (covs, w, n, DATE_TIME_OFFSET)
    }

    /// Two-pass reliability-weighted variance: the weighted mean first, then the
    /// weighted squared deviations from it. The deviations are formed at the
    /// column's own scale rather than as a difference of two large sums, so this
    /// stays accurate at any offset and is the reference `weighted_moments` has to
    /// reproduce.
    fn two_pass_variance(x: &[f64], w: &[f64]) -> f64 {
        let sw: f64 = w.iter().sum();
        let sw2: f64 = w.iter().map(|wi| wi * wi).sum();
        let mean: f64 = x.iter().zip(w).map(|(&xi, &wi)| wi * xi).sum::<f64>() / sw;
        let ss: f64 = x
            .iter()
            .zip(w)
            .map(|(&xi, &wi)| {
                let d = xi - mean;
                wi * d * d
            })
            .sum();
        (ss / sw) / (1.0 - sw2 / (sw * sw))
    }

    #[test]
    fn the_standardizing_sd_survives_a_date_time_offset() {
        // Subtracting the offset from a column built as `offset + spread * z` is
        // exact in binary floating point, both values being within a factor of
        // two of each other, so the centered column holds precisely the
        // deviations the stored column has and its two-pass variance is the
        // variance the transform is being asked for.
        let (covs, w, n, offset) = offset_fixture();
        let centered: Vec<f64> = covs[..n].iter().map(|&xi| xi - offset).collect();
        let want = two_pass_variance(&centered, &w).sqrt();

        let (_means, vars) = weighted_moments(&covs, n, 2, &w, 1);
        let got = vars[0].sqrt();
        let rel = (got - want).abs() / want;
        assert!(
            rel < 1e-12,
            "standardizing sd {got} against the two-pass reference {want}, relative error {rel}"
        );
    }

    /// Weighted mean of a single column, the reference the centered transform's
    /// output is read against.
    fn two_pass_mean(x: &[f64], w: &[f64]) -> f64 {
        let sw: f64 = w.iter().sum();
        x.iter().zip(w).map(|(&xi, &wi)| wi * xi).sum::<f64>() / sw
    }

    #[test]
    fn a_scaled_column_has_zero_mean_and_unit_weighted_sd_at_a_date_time_offset() {
        // Scaled Euclidean subtracts each column's weighted mean and divides by
        // its weighted standard deviation, so every transformed column must read
        // back a weighted mean of zero and a weighted standard deviation of one.
        // The offset column is the demanding one; the plain column establishes
        // that the check itself is satisfiable.
        //
        // The standard deviation is held to 1e-12 because centering leaves the
        // output at the deviations' own scale. Without it the offset column would
        // come back around 8e5 with deviations of order one, so a unit in the
        // last place of the output would be about 1e-10 and any reading taken
        // from it would inherit that granularity, ten digits coarser than the
        // statistic the transform computed.
        //
        // The mean is held to 1e-9 instead, and the looser bound is a property of
        // the stored input rather than of the centering. A double near 1.7e9
        // resolves to about 2.4e-7, so the exact weighted mean of the stored
        // column cannot be named to better than half of that, which is 6e-11 of
        // the hour of spread the column is divided by. No arrangement of the
        // arithmetic reaches 1e-12 here; the measured residual, 9e-11, is already
        // inside one unit in the last place of the input read in standard
        // deviations. The plain column, which has no such floor, comes back at
        // 4e-17 and shows the check is not merely loose.
        let (covs, w, n, _offset) = offset_fixture();
        let out = scaled_euclidean(&covs, n, 2, &w, 1);
        for j in 0..2 {
            let column = &out[j * n..(j + 1) * n];
            let mean = two_pass_mean(column, &w);
            assert!(
                mean.abs() < 1e-9,
                "column {j} has weighted mean {mean} after scaling"
            );
            let sd = two_pass_variance(column, &w).sqrt();
            assert!(
                (sd - 1.0).abs() < 1e-12,
                "column {j} has weighted sd {sd} after scaling"
            );
        }
    }

    #[test]
    fn centering_leaves_the_scaled_euclidean_distances_where_they_were() {
        // Subtracting a per-column constant shifts every row by the same vector,
        // which cancels in every pairwise difference, so centering is a change of
        // representation and not of the metric. On well-scaled covariates, where
        // the uncentered form loses nothing to cancellation, the two must agree
        // to rounding: the reference divides by the same standard deviations
        // without subtracting the means.
        let n = 40;
        let p = 3;
        let draws = lcg_unit(n * p, 24_680_135);
        let covs: Vec<f64> = draws
            .iter()
            .enumerate()
            .map(|(k, &u)| 2.0 * u - 1.0 + 3.0 * ((k / n) as f64 + 1.0))
            .collect();
        let weight_draws = lcg_unit(n, 13_579_246);
        let w: Vec<f64> = weight_draws.iter().map(|&vi| 0.25 + 1.5 * vi).collect();

        let (means, vars) = weighted_moments(&covs, n, p, &w, 1);
        // Every column sits well away from zero, or centering would be a no-op
        // and the agreement below would hold for the wrong reason.
        for (j, mean) in means.iter().enumerate() {
            assert!(*mean > 1.0, "column {j} has weighted mean {mean}");
        }
        let mut uncentered = covs.clone();
        for (j, var) in vars.iter().enumerate() {
            let sd = var.sqrt();
            let scale = if sd > 0.0 { 1.0 / sd } else { 1.0 };
            for i in 0..n {
                uncentered[j * n + i] *= scale;
            }
        }

        let out = transform(&covs, n, p, Distance::ScaledEuclidean, &w, 1);
        let got = pairwise::euclidean(&out, n, p, 1);
        let want = pairwise::euclidean(&uncentered, n, p, 1);
        for (k, (&a, &b)) in got.iter().zip(&want).enumerate() {
            assert!(
                (a - b).abs() < 1e-12,
                "entry {k}: centered distance {a}, uncentered reference {b}"
            );
        }
    }

    #[test]
    fn the_moments_are_bit_identical_across_thread_counts() {
        // The chunked reduction fixes its summation tree independently of the
        // pool size, but a fixture below `REDUCE_CHUNK` never forms more than one
        // chunk and so never exercises the combine. These sizes do: one an exact
        // multiple of the chunk length and one that leaves a short final chunk.
        // The offset column is carried because a large mean is where a reordered
        // sum shows first.
        for n in [2 * REDUCE_CHUNK, 10_000] {
            let (covs, w) = offset_fixture_n(n);
            let (means_one, vars_one) = weighted_moments(&covs, n, 2, &w, 1);
            let (means_many, vars_many) = weighted_moments(&covs, n, 2, &w, 4);
            for j in 0..2 {
                assert_eq!(
                    means_one[j].to_bits(),
                    means_many[j].to_bits(),
                    "n {n}, column {j}: mean {} against {}",
                    means_one[j],
                    means_many[j]
                );
                assert_eq!(
                    vars_one[j].to_bits(),
                    vars_many[j].to_bits(),
                    "n {n}, column {j}: variance {} against {}",
                    vars_one[j],
                    vars_many[j]
                );
            }
        }
    }

    #[test]
    fn the_transforms_are_bit_identical_across_thread_counts() {
        // The same guarantee at the transform's own output. Mahalanobis reduces a
        // second time for the covariance, so it needs its own multi-chunk check
        // rather than inheriting the moments one.
        for n in [2 * REDUCE_CHUNK, 10_000] {
            let (covs, w) = offset_fixture_n(n);
            for distance in [Distance::ScaledEuclidean, Distance::Mahalanobis] {
                let one = transform(&covs, n, 2, distance, &w, 1);
                let many = transform(&covs, n, 2, distance, &w, 4);
                assert_eq!(one.len(), many.len());
                for (k, (&a, &b)) in one.iter().zip(&many).enumerate() {
                    assert_eq!(
                        a.to_bits(),
                        b.to_bits(),
                        "n {n}, {distance:?}, entry {k}: {a} against {b}"
                    );
                }
            }
        }
    }
}
