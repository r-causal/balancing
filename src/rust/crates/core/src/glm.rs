//! A small iteratively reweighted least squares fit for starting values.
//!
//! The tilting solver converges from any finite start because its Jacobian is
//! positive definite everywhere, but a propensity-model warm start shortens the
//! path and steadies the iterations when a link pushes propensities toward zero
//! or one. This module fits a binomial GLM with the requested link over a chosen
//! set of units and returns the coefficient vector; nothing here is exposed to R.

use crate::linalg::solve_symmetric_ridge;
use crate::links::Link;

/// Propensities are clamped away from the boundary so the IRLS working weights
/// stay finite when a saturated design drives a fitted value to zero or one.
const P_FLOOR: f64 = 1e-10;
/// Working weights are floored so a near-flat link region cannot make the normal
/// equations singular during the warm start.
const W_FLOOR: f64 = 1e-10;

/// Fit `p_i = link(x_i . beta)` to the binary `response` over the units in `idx`
/// by iteratively reweighted least squares, returning the coefficients.
///
/// `covs` is the full column-major `n` by `p` design; `idx` selects the fitting
/// units and `response` holds their zero/one outcomes in the same order. `s` is
/// the length-`n` sampling weight. The fit is a starting value, so it runs a
/// small fixed number of iterations and returns all-zero coefficients if it ever
/// leaves the finite range rather than propagating a non-finite start.
// The design matrix, its shape, the fitting subset, the response, the sampling
// weights, the link, and the iteration cap are each irreducible inputs to a GLM.
#[allow(clippy::too_many_arguments)]
pub fn irls(
    covs: &[f64],
    n: usize,
    p: usize,
    idx: &[usize],
    response: &[f64],
    s: &[f64],
    link: Link,
    max_iter: usize,
) -> Vec<f64> {
    debug_assert_eq!(response.len(), idx.len());
    let mut beta = vec![0.0; p];
    let mut xtwx = vec![0.0; p * p];
    let mut xtwz = vec![0.0; p];
    let mut row = vec![0.0; p];

    for _ in 0..max_iter {
        xtwx.iter_mut().for_each(|v| *v = 0.0);
        xtwz.iter_mut().for_each(|v| *v = 0.0);

        for (local, &i) in idx.iter().enumerate() {
            for (j, r) in row.iter_mut().enumerate() {
                *r = covs[j * n + i];
            }
            let eta: f64 = row.iter().zip(&beta).map(|(c, b)| c * b).sum();
            let mut prob = link.linkinv(eta);
            prob = prob.clamp(P_FLOOR, 1.0 - P_FLOOR);
            let mu_eta = link.mu_eta(eta);
            let var = prob * (1.0 - prob);
            // IRLS working weight and response for a binomial GLM.
            let w = (s[i] * mu_eta * mu_eta / var).max(W_FLOOR);
            let z = eta + (response[local] - prob) / mu_eta;
            for j in 0..p {
                let wxj = w * row[j];
                xtwz[j] += wxj * z;
                for k in 0..p {
                    xtwx[j * p + k] += wxj * row[k];
                }
            }
        }

        let mut next = xtwz.clone();
        // A ridge rescues a rank-deficient warm-start system; the tilting solver
        // that follows adds its own escalation, so a rough solve here suffices.
        if !solve_symmetric_ridge(&xtwx, p, 0.0, &mut next)
            && !solve_symmetric_ridge(&xtwx, p, 1e-8, &mut next)
        {
            return vec![0.0; p];
        }
        if next.iter().any(|v| !v.is_finite()) {
            return vec![0.0; p];
        }

        let delta: f64 = next
            .iter()
            .zip(&beta)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        beta = next;
        if delta < 1e-8 {
            break;
        }
    }

    if beta.iter().all(|v| v.is_finite()) {
        beta
    } else {
        vec![0.0; p]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // On a saturated single-covariate design the IRLS fit recovers the empirical
    // cell propensities, so the linear predictor reproduces them through the link.
    #[test]
    fn recovers_cell_propensities_on_a_saturated_design() {
        // Two cells via an intercept and a 0/1 covariate. Cell z = 0 has 2 of 4
        // positive; cell z = 1 has 1 of 4 positive.
        let n = 8;
        let z = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let y = [1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
        let mut covs = vec![1.0; n];
        covs.extend_from_slice(&z);
        let idx: Vec<usize> = (0..n).collect();
        let s = vec![1.0; n];
        let beta = irls(&covs, n, 2, &idx, &y, &s, Link::Logit, 50);

        let prob = |i: usize| Link::Logit.linkinv(beta[0] + beta[1] * z[i]);
        assert!(
            (prob(0) - 0.5).abs() < 1e-6,
            "cell 0 propensity {}",
            prob(0)
        );
        assert!(
            (prob(4) - 0.25).abs() < 1e-6,
            "cell 1 propensity {}",
            prob(4)
        );
    }

    #[test]
    fn returns_finite_coefficients_for_a_probit_fit() {
        let n = 6;
        let x = [0.0, 1.0, -1.0, 2.0, -2.0, 0.5];
        let y = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0];
        let mut covs = vec![1.0; n];
        covs.extend_from_slice(&x);
        let idx: Vec<usize> = (0..n).collect();
        let s = vec![1.0; n];
        let beta = irls(&covs, n, 2, &idx, &y, &s, Link::Probit, 50);
        assert!(beta.iter().all(|v| v.is_finite()));
    }
}
