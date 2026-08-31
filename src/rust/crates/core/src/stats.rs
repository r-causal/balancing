//! Weighted moments shared by the standardizing statistics.
//!
//! Every method that standardizes a covariate or an exposure needs the same
//! reliability-weighted variance, and each one needs it computed in a way that
//! survives a large offset. The compact one-pass form, `sum_i w_i x_i^2 / sw`
//! minus the squared mean, subtracts two quantities that agree in their leading
//! digits whenever the mean is large against the spread, so the digits carrying
//! the variance are exactly the ones the subtraction destroys. A date-time
//! covariate is the ordinary case rather than an exotic one: `as.numeric()` on a
//! POSIXct is about 1.7e9, and an hour of spread around that leaves only a few
//! digits of the double intact, which cost the standardizing standard deviation
//! a relative error near 3e-4. Accumulating deviations from the weighted mean
//! instead forms every term at the spread's own scale, so the accuracy no longer
//! depends on where the column sits.
//!
//! Sharing the denominator as well as the deviations matters because the methods
//! compare against each other: a column standardized by the distance transforms
//! and the same column standardized inside a continuous solve have to agree, so
//! that a bounded weighted product of standardized columns reads as a bounded
//! correlation.

/// Finish a reliability-weighted variance from a weighted sum of squared
/// deviations and the weight totals.
///
/// The denominator is the frequency-weight-free (reliability) form
/// `1 - sum_i wn_i^2` with normalized weights `wn = w / sum(w)`, which reduces to
/// the usual `n - 1` sample variance when the weights are equal. A non-positive
/// total weight reports zero, as does a single dominant weight that drives the
/// denominator to zero or below, and a sum of squares that rounds negative is
/// clamped up.
pub(crate) fn reliability_variance(ss: f64, sw: f64, sw2: f64) -> f64 {
    if sw <= 0.0 {
        return 0.0;
    }
    let denom = 1.0 - sw2 / (sw * sw);
    if denom > 0.0 {
        ((ss / sw) / denom).max(0.0)
    } else {
        0.0
    }
}

/// Reliability-weighted variance of a single vector.
///
/// The weight totals and the weighted mean come from the first pass, the squared
/// deviations from that mean from the second. A vector of no positive total
/// weight has variance zero.
pub(crate) fn weighted_variance(x: &[f64], w: &[f64]) -> f64 {
    let mut sw = 0.0;
    let mut sw2 = 0.0;
    let mut swx = 0.0;
    for (&xi, &wi) in x.iter().zip(w) {
        sw += wi;
        sw2 += wi * wi;
        swx += wi * xi;
    }
    if sw <= 0.0 {
        return 0.0;
    }
    let mean = swx / sw;
    let mut ss = 0.0;
    for (&xi, &wi) in x.iter().zip(w) {
        let d = xi - mean;
        ss += wi * d * d;
    }
    reliability_variance(ss, sw, sw2)
}
