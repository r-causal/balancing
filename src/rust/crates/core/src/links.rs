//! Link functions for the propensity models behind inverse probability tilting.
//!
//! A link maps a linear predictor `eta` to a propensity `p` in `(0, 1)`. Each
//! link supplies three quantities the tilting solver and its warm-start GLM
//! need: `linkinv` (`p = G(eta)`), `mu_eta` (`dp/deta`), and `variance`
//! (`p (1 - p)`, the binomial variance the IRLS weights use).
//!
//! The probit link evaluates the standard normal distribution function through a
//! Cody rational approximation ported from the reference `pnorm` implementation,
//! which reaches full double precision including the tails. The unit tests check
//! it against tabulated constants at `1e-15`.
//!
//! The v1 exposures use logit, probit, and cloglog. The enum also names the
//! log-log, log, and complementary-log links so the set is fixed in one place;
//! their evaluators are wired when an exposure that needs them is added.

/// A propensity link.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Link {
    /// Logistic link, `p = 1 / (1 + exp(-eta))`.
    Logit,
    /// Probit link, `p = Phi(eta)`.
    Probit,
    /// Complementary log-log link, `p = 1 - exp(-exp(eta))`.
    Cloglog,
    /// Log-log link. Named for the fixed link set; not yet wired.
    #[allow(
        dead_code,
        reason = "part of the planned link set; wired when an exposure uses it"
    )]
    Loglog,
    /// Log link. Named for the fixed link set; not yet wired.
    #[allow(
        dead_code,
        reason = "part of the planned link set; wired when an exposure uses it"
    )]
    Log,
    /// Complementary log link. Named for the fixed link set; not yet wired.
    #[allow(
        dead_code,
        reason = "part of the planned link set; wired when an exposure uses it"
    )]
    Clog,
}

impl Link {
    /// Resolve a link by the name the R layer passes across the boundary.
    ///
    /// Only the wired links are accepted; an unknown or not-yet-wired name
    /// returns `None` so the boundary can raise a contract error.
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "logit" => Some(Link::Logit),
            "probit" => Some(Link::Probit),
            "cloglog" => Some(Link::Cloglog),
            _ => None,
        }
    }

    /// The propensity `p = G(eta)`.
    pub fn linkinv(self, eta: f64) -> f64 {
        match self {
            Link::Logit => {
                // Branch on the sign of eta so the exponential never overflows:
                // both forms are algebraically identical but each is evaluated
                // where its exponent is non-positive.
                if eta >= 0.0 {
                    1.0 / (1.0 + (-eta).exp())
                } else {
                    let e = eta.exp();
                    e / (1.0 + e)
                }
            }
            Link::Probit => pnorm(eta),
            Link::Cloglog => -(-eta.exp()).exp_m1(),
            Link::Loglog | Link::Log | Link::Clog => unimplemented!("link not yet wired"),
        }
    }

    /// The derivative `dp/deta = G'(eta)`.
    pub fn mu_eta(self, eta: f64) -> f64 {
        match self {
            Link::Logit => {
                // exp(-|eta|) / (1 + exp(-|eta|))^2, the symmetric form that
                // holds full precision where `p (1 - p)` would lose a bit.
                let e = (-eta.abs()).exp();
                let d = 1.0 + e;
                e / (d * d)
            }
            Link::Probit => dnorm(eta),
            Link::Cloglog => (eta - eta.exp()).exp(),
            Link::Loglog | Link::Log | Link::Clog => unimplemented!("link not yet wired"),
        }
    }

    /// The binomial variance `p (1 - p)` at `p = G(eta)`.
    pub fn variance(self, eta: f64) -> f64 {
        let p = self.linkinv(eta);
        p * (1.0 - p)
    }
}

/// `1 / sqrt(2 pi)`, the standard normal density's normalizing constant.
const M_1_SQRT_2PI: f64 = 0.398_942_280_401_432_677_939_946_059_934;
/// `sqrt(32)`, the boundary between the rational and tail regions of `pnorm`.
const M_SQRT_32: f64 = 5.656_854_249_492_380_195_206_754_896_838;

/// Standard normal density, `phi(x)`.
fn dnorm(x: f64) -> f64 {
    M_1_SQRT_2PI * (-0.5 * x * x).exp()
}

/// Standard normal distribution function `Phi(x)`, lower tail.
///
/// Cody's rational Chebyshev approximation over three regions, matching the
/// reference `pnorm` to full double precision. The tail region factors the
/// exponent as `exp(-xsq^2 / 2) exp(-del / 2)` with `xsq` truncated to a
/// multiple of `1/16`, which keeps the leading exponential exact and confines
/// rounding to the small correction `del`.
fn pnorm(x: f64) -> f64 {
    const A: [f64; 5] = [
        2.235_252_035_460_683_9e0,
        1.610_282_310_685_558_8e2,
        1.067_689_485_460_370_9e3,
        1.815_498_125_334_356e4,
        6.568_233_791_820_745e-2,
    ];
    const B: [f64; 4] = [
        4.720_258_190_468_824_2e1,
        9.760_985_517_377_767e2,
        1.026_093_220_861_897_8e4,
        4.550_778_933_502_672_9e4,
    ];
    const C: [f64; 9] = [
        3.989_415_120_881_346_7e-1,
        8.883_149_794_388_377,
        9.350_665_613_217_785e1,
        5.972_702_763_948_002e2,
        2.494_537_585_290_372_7e3,
        6.848_190_450_536_283e3,
        1.160_265_143_764_735e4,
        9.842_714_838_383_978e3,
        1.076_557_677_372_019_2e-8,
    ];
    const D: [f64; 8] = [
        2.226_668_804_432_811_6e1,
        2.353_879_017_826_25e2,
        1.519_377_599_407_554_8e3,
        6.485_558_298_266_761e3,
        1.861_557_164_088_509_8e4,
        3.490_095_272_114_598e4,
        3.891_200_328_609_327e4,
        1.968_542_967_685_999_1e4,
    ];
    const P: [f64; 6] = [
        2.158_985_340_579_57e-1,
        1.274_011_611_602_473_6e-1,
        2.223_527_787_064_980_7e-2,
        1.421_619_193_227_893_5e-3,
        2.911_287_495_116_879_2e-5,
        2.307_344_176_494_017_3e-2,
    ];
    const Q: [f64; 5] = [
        1.284_260_096_144_911,
        4.682_382_124_808_651e-1,
        6.598_813_786_892_856e-2,
        3.782_396_332_027_582_4e-3,
        7.297_515_550_839_662e-5,
    ];
    /// Half the machine epsilon, the threshold below which `Phi` is linear.
    const EPS: f64 = 1.110_223_024_625_156_5e-16;

    if x.is_nan() {
        return f64::NAN;
    }
    let y = x.abs();

    if y <= 0.674_489_75 {
        let mut xnum = 0.0;
        let mut xden = 0.0;
        if y > EPS {
            let xsq = x * x;
            xnum = A[4] * xsq;
            xden = xsq;
            for i in 0..3 {
                xnum = (xnum + A[i]) * xsq;
                xden = (xden + B[i]) * xsq;
            }
        }
        let temp = x * (xnum + A[3]) / (xden + B[3]);
        return 0.5 + temp;
    }

    let ccum = if y <= M_SQRT_32 {
        let mut xnum = C[8] * y;
        let mut xden = y;
        for i in 0..7 {
            xnum = (xnum + C[i]) * y;
            xden = (xden + D[i]) * y;
        }
        let temp = (xnum + C[7]) / (xden + D[7]);
        let xsq = (y * 16.0).trunc() / 16.0;
        let del = (y - xsq) * (y + xsq);
        (-xsq * xsq * 0.5).exp() * (-del * 0.5).exp() * temp
    } else {
        // Beyond this range one tail underflows to zero at double precision.
        let xsq = 1.0 / (x * x);
        let mut xnum = P[5] * xsq;
        let mut xden = xsq;
        for i in 0..4 {
            xnum = (xnum + P[i]) * xsq;
            xden = (xden + Q[i]) * xsq;
        }
        let mut temp = xsq * (xnum + P[4]) / (xden + Q[4]);
        temp = (M_1_SQRT_2PI - temp) / y;
        let xsq = (y * 16.0).trunc() / 16.0;
        let del = (y - xsq) * (y + xsq);
        (-xsq * xsq * 0.5).exp() * (-del * 0.5).exp() * temp
    };

    // `ccum` is the upper tail for the positive magnitude `y`; map it back to the
    // lower tail at the signed argument.
    if x > 0.0 { 1.0 - ccum } else { ccum }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Assert a relative match to a tabulated constant at the stated tolerance.
    fn close(actual: f64, expected: f64, rel: f64) {
        let scale = expected.abs().max(1e-300);
        assert!(
            (actual - expected).abs() <= rel * scale,
            "actual {actual:.17e} expected {expected:.17e} (rel {})",
            (actual - expected).abs() / scale
        );
    }

    // Probit propensities against pnorm to full double precision, including the
    // tails where a naive erfc loses digits.
    #[test]
    fn probit_matches_pnorm_at_tabulated_points() {
        let cases = [
            (0.0, 0.5),
            (0.5, 0.691_462_461_274_013),
            (1.0, 0.841_344_746_068_542_9),
            (-1.0, 0.158_655_253_931_457_05),
            (2.5, 0.993_790_334_674_223_8),
            (-2.5, 0.006_209_665_325_776_135),
            (3.0, 0.998_650_101_968_369_9),
            (-3.0, 0.001_349_898_031_630_094_6),
            (5.0, 0.999_999_713_348_428_1),
            (-5.0, 2.866_515_718_791_939e-7),
            (8.0, 0.999_999_999_999_999_3),
            (-8.0, 6.220_960_574_271_785e-16),
            (0.674_489_750_196_082, 0.750_000_000_000_000_1),
        ];
        for (eta, expected) in cases {
            close(Link::Probit.linkinv(eta), expected, 1e-15);
        }
    }

    #[test]
    fn probit_density_matches_dnorm() {
        let cases = [
            (0.0, 0.398_942_280_401_432_7),
            (1.0, 0.241_970_724_519_143_37),
            (-2.5, 0.017_528_300_493_568_54),
            (5.0, 1.486_719_514_734_297_7e-6),
            (8.0, 5.052_271_083_536_893e-15),
        ];
        for (eta, expected) in cases {
            close(Link::Probit.mu_eta(eta), expected, 1e-15);
        }
    }

    #[test]
    fn logit_matches_plogis() {
        let cases = [
            (0.0, 0.5),
            (1.0, 0.731_058_578_630_004_9),
            (-1.0, 0.268_941_421_369_995_1),
            (2.0, 0.880_797_077_977_882_3),
            (-2.0, 0.119_202_922_022_117_55),
            (4.0, 0.982_013_790_037_908_5),
        ];
        for (eta, expected) in cases {
            close(Link::Logit.linkinv(eta), expected, 1e-15);
        }
    }

    #[test]
    fn logit_density_matches_dlogis() {
        let cases = [
            (0.0, 0.25),
            (1.0, 0.196_611_933_241_481_88),
            (2.0, 0.104_993_585_403_506_5),
            (4.0, 0.017_662_706_213_291_114),
        ];
        for (eta, expected) in cases {
            close(Link::Logit.mu_eta(eta), expected, 1e-15);
        }
    }

    #[test]
    fn cloglog_matches_reference() {
        let inv = [
            (0.0, 0.632_120_558_828_557_7),
            (1.0, 0.934_011_964_154_687_5),
            (-1.0, 0.307_799_372_444_653_6),
            (0.5, 0.807_704_354_452_035),
            (-2.0, 0.126_576_981_506_883_35),
        ];
        for (eta, expected) in inv {
            close(Link::Cloglog.linkinv(eta), expected, 1e-15);
        }
        let deriv = [
            (0.0, 0.367_879_441_171_442_33),
            (1.0, 0.179_374_078_734_017_2),
            (-1.0, 0.254_646_380_043_582_5),
            (0.5, 0.317_041_921_077_942_16),
            (-2.0, 0.118_204_951_593_143_13),
        ];
        for (eta, expected) in deriv {
            close(Link::Cloglog.mu_eta(eta), expected, 1e-15);
        }
    }

    #[test]
    fn logit_is_stable_in_the_tails() {
        // No overflow at large magnitude, and the two symmetric forms agree.
        assert!(Link::Logit.linkinv(1000.0) == 1.0);
        assert!(Link::Logit.linkinv(-1000.0) == 0.0);
        close(
            Link::Logit.linkinv(30.0),
            1.0 - Link::Logit.linkinv(-30.0),
            1e-15,
        );
    }

    #[test]
    fn from_name_wires_only_the_v1_links() {
        assert_eq!(Link::from_name("logit"), Some(Link::Logit));
        assert_eq!(Link::from_name("probit"), Some(Link::Probit));
        assert_eq!(Link::from_name("cloglog"), Some(Link::Cloglog));
        assert_eq!(Link::from_name("loglog"), None);
        assert_eq!(Link::from_name("identity"), None);
    }
}
