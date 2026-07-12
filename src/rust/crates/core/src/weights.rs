//! Balancing weights from a propensity and their derivatives, per estimand.
//!
//! Inverse probability tilting reweights each treatment level to a target
//! population. The weight a level's units carry is a function of that level's
//! modeled propensity `p`, and the estimand fixes which function:
//!
//! - the average treatment effect tilts every level to the whole sample, and a
//!   level's units carry the inverse weight `1 / p`;
//! - a focal estimand (treatment effect on the treated or on the controls) tilts
//!   each non-focal level to the focal level, and a level's units carry
//!   `(1 - p) / p`, while the focal level's units carry weight one.
//!
//! The two forms differ by the constant one, so their derivative in the linear
//! predictor is identical, `dw/deta = -mu_eta / p^2`. That derivative feeds both
//! the tilting Jacobian and the weight-derivative block of the estimating
//! equations.

/// Which weight function a tilting block applies to its units.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WeightForm {
    /// `w = 1 / p`, tilting the level to the whole sample.
    Inverse,
    /// `w = (1 - p) / p`, tilting a non-focal level to the focal level.
    InverseComplement,
}

impl WeightForm {
    /// The weight at propensity `p`.
    pub fn weight(self, p: f64) -> f64 {
        match self {
            WeightForm::Inverse => 1.0 / p,
            WeightForm::InverseComplement => (1.0 - p) / p,
        }
    }
}

/// The weight derivative in the linear predictor, `dw/deta = -mu_eta / p^2`.
///
/// Both weight forms share this derivative because they differ only by an
/// additive constant. `mu_eta` is `dp/deta` from the link.
pub fn weight_deriv_eta(p: f64, mu_eta: f64) -> f64 {
    -mu_eta / (p * p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse_weight_is_reciprocal_propensity() {
        assert!((WeightForm::Inverse.weight(0.25) - 4.0).abs() < 1e-15);
        assert!((WeightForm::Inverse.weight(0.5) - 2.0).abs() < 1e-15);
    }

    #[test]
    fn complement_weight_is_the_odds_against() {
        // (1 - p) / p at p = 1/3 is 2.
        assert!((WeightForm::InverseComplement.weight(1.0 / 3.0) - 2.0).abs() < 1e-15);
    }

    #[test]
    fn both_forms_are_positive_across_the_unit_interval() {
        for k in 1..100 {
            let p = k as f64 / 100.0;
            assert!(WeightForm::Inverse.weight(p) > 0.0);
            assert!(WeightForm::InverseComplement.weight(p) > 0.0);
        }
    }

    #[test]
    fn derivative_matches_a_finite_difference() {
        // Logit at eta = 0.4: p and mu_eta by hand, checked against a central
        // difference of 1 / p through the link.
        let eta = 0.4;
        let link = crate::links::Link::Logit;
        let h = 1e-6;
        let w = |e: f64| 1.0 / link.linkinv(e);
        let fd = (w(eta + h) - w(eta - h)) / (2.0 * h);
        let analytic = weight_deriv_eta(link.linkinv(eta), link.mu_eta(eta));
        assert!((fd - analytic).abs() < 1e-6, "fd {fd} analytic {analytic}");
    }
}
