//! Positive-semidefinite kernels for characteristic function distance balancing.
//!
//! Characteristic function distance balancing measures covariate balance through
//! a kernel embedding of the covariate rows. Each kernel is a function of the
//! pairwise Euclidean distances of the standardized covariates: the covariates
//! are first put on the scaled-Euclidean footing (each column divided by its
//! weighted standard deviation), the pairwise distance matrix is formed, and a
//! bandwidth-scaled transform is applied entrywise. Standardizing first makes the
//! kernel invariant to a uniform rescaling of the covariates, and the bandwidth,
//! the median of the pairwise distances times a scale factor, tracks the spread of
//! the data so a single kernel serves problems of different scales.
//!
//! Every kernel here yields a positive-semidefinite Gram matrix except the energy
//! kernel, the negative distance, which is only conditionally positive
//! semidefinite; that is why characteristic function distance balancing with the
//! energy kernel reproduces energy balancing, whose quadratic term is indefinite.
//! The Matern kernel is provided at the three half-integer smoothness values whose
//! closed forms avoid the Bessel-function evaluation of the general case. The t
//! kernel is a random-Fourier-feature approximation whose frequency vectors are
//! drawn on the R side under R's generator and passed in, so a fixed seed
//! reproduces the kernel exactly while all arithmetic stays in Rust.

use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, ParallelExtend, ParallelIterator,
};
use rayon::slice::ParallelSliceMut;

use super::{Distance, pairwise, transform};
use crate::threads::get_pool;

/// The kernel a characteristic function distance solve is built on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kernel {
    /// The negative pairwise distance. Only conditionally positive semidefinite,
    /// so the assembled quadratic term is indefinite and the solve reproduces
    /// energy balancing.
    Energy,
    /// The squared-exponential kernel `exp(-d^2 / (2 bw^2))`.
    Gaussian,
    /// The Laplacian kernel `exp(-d / bw)`.
    Laplace,
    /// The Matern kernel at half-integer smoothness (0.5, 1.5, or 2.5).
    Matern,
    /// A random-Fourier-feature approximation of the multivariate-t characteristic
    /// function kernel.
    T,
}

impl Kernel {
    /// Resolve the kernel named by the R layer, returning `None` for an
    /// unrecognized name.
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "energy" => Some(Kernel::Energy),
            "gaussian" => Some(Kernel::Gaussian),
            "laplace" => Some(Kernel::Laplace),
            "matern" => Some(Kernel::Matern),
            "t" => Some(Kernel::T),
            _ => None,
        }
    }

    /// Whether the kernel's Gram matrix is positive semidefinite. The energy
    /// kernel is only conditionally positive semidefinite, so its assembled
    /// quadratic term is indefinite; the others are positive semidefinite by
    /// construction, which makes the interior-point backend eligible.
    pub fn is_psd(self) -> bool {
        !matches!(self, Kernel::Energy)
    }
}

/// The Matern smoothness values with a closed form in this version.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MaternNu {
    Half,
    ThreeHalves,
    FiveHalves,
}

impl MaternNu {
    /// Resolve a smoothness value to a supported half-integer, returning `None`
    /// for an unsupported value.
    fn from_smoothness(smoothness: f64) -> Option<Self> {
        if (smoothness - 0.5).abs() < 1e-9 {
            Some(MaternNu::Half)
        } else if (smoothness - 1.5).abs() < 1e-9 {
            Some(MaternNu::ThreeHalves)
        } else if (smoothness - 2.5).abs() < 1e-9 {
            Some(MaternNu::FiveHalves)
        } else {
            None
        }
    }
}

/// Kernel tuning that crosses from the R layer.
pub struct KernelParams<'a> {
    /// Which kernel to build.
    pub kernel: Kernel,
    /// Bandwidth scale factor multiplying the median pairwise distance.
    pub bw_scale: f64,
    /// Matern smoothness; ignored unless `kernel` is `Matern`.
    pub smoothness: f64,
    /// Column-major `p` by `n_draws` frequency vectors for the t kernel, drawn on
    /// the R side; empty for the other kernels.
    pub t_proj: &'a [f64],
    /// Number of Monte Carlo draws (columns of `t_proj`) for the t kernel.
    pub n_draws: usize,
}

const SQRT3: f64 = 1.732_050_807_568_877_2;
const SQRT5: f64 = 2.236_067_977_499_79;

/// The squared-exponential kernel at scaled radius `r = d / bw`.
#[inline]
fn gaussian(r: f64) -> f64 {
    (-0.5 * r * r).exp()
}

/// The Laplacian kernel at scaled radius `r = d / bw`.
#[inline]
fn laplace(r: f64) -> f64 {
    (-r).exp()
}

/// The Matern kernel with smoothness one half: `exp(-r)`.
#[inline]
fn matern_half(r: f64) -> f64 {
    (-r).exp()
}

/// The Matern kernel with smoothness three halves: `(1 + sqrt(3) r) exp(-sqrt(3) r)`.
#[inline]
fn matern_three_halves(r: f64) -> f64 {
    let a = SQRT3 * r;
    (1.0 + a) * (-a).exp()
}

/// The Matern kernel with smoothness five halves:
/// `(1 + sqrt(5) r + 5 r^2 / 3) exp(-sqrt(5) r)`.
#[inline]
fn matern_five_halves(r: f64) -> f64 {
    let a = SQRT5 * r;
    (1.0 + a + a * a / 3.0) * (-a).exp()
}

/// The bandwidth: the median of the lower triangle of the distance matrix times
/// the scale factor.
///
/// Pairs where either unit is discarded do not enter the median, so a caliper's
/// discarded units cannot shift the bandwidth. The lower triangle is gathered in
/// parallel, and the median is found by a partial partition rather than a full
/// sort. A degenerate result (no admissible pairs, or a zero median) falls back to
/// one so the kernels never divide by zero.
fn median_bandwidth(
    dist: &[f64],
    n: usize,
    discarded: &[bool],
    bw_scale: f64,
    threads: usize,
) -> f64 {
    let disc = |i: usize| discarded.get(i).copied().unwrap_or(false);
    let pool = get_pool(threads);
    let mut lower: Vec<f64> = Vec::new();
    pool.install(|| {
        // The median does not depend on the gather order, so the parallel
        // extension is deterministic despite the unordered append.
        lower.par_extend((0..n).into_par_iter().flat_map_iter(|j| {
            ((j + 1)..n).filter_map(move |i| {
                if disc(i) || disc(j) {
                    None
                } else {
                    Some(dist[j * n + i])
                }
            })
        }));
    });

    let median = median_of(&mut lower);
    let bw = median * bw_scale;
    if bw > 0.0 { bw } else { 1.0 }
}

/// The median of a slice, found with a partial partition. The slice is reordered
/// in place. An empty slice has median one, the neutral bandwidth.
fn median_of(values: &mut [f64]) -> f64 {
    let m = values.len();
    if m == 0 {
        return 1.0;
    }
    if m % 2 == 1 {
        let (_, mid, _) = values.select_nth_unstable_by(m / 2, |a, b| a.total_cmp(b));
        *mid
    } else {
        let (lo_part, hi, _) = values.select_nth_unstable_by(m / 2, |a, b| a.total_cmp(b));
        let lo = lo_part.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        (lo + *hi) / 2.0
    }
}

/// Build the `n` by `n` kernel matrix for `covs` under `params`.
///
/// `covs` is a column-major `n` by `p` matrix and `s` the per-unit weights used
/// to standardize the covariates. `discarded` marks units excluded from the
/// bandwidth median; an empty slice discards none. The result is a column-major
/// `n` by `n` symmetric matrix.
pub fn build_kernel(
    covs: &[f64],
    n: usize,
    p: usize,
    s: &[f64],
    params: &KernelParams<'_>,
    discarded: &[bool],
    threads: usize,
) -> Vec<f64> {
    let z = transform::transform(covs, n, p, Distance::ScaledEuclidean, s, threads);
    let pt = if n == 0 { 0 } else { z.len() / n };

    if params.kernel == Kernel::T {
        return t_kernel(&z, n, pt, params, threads);
    }

    let dist = pairwise::euclidean(&z, n, pt, threads);
    if params.kernel == Kernel::Energy {
        return dist.iter().map(|&d| -d).collect();
    }

    let bw = median_bandwidth(&dist, n, discarded, params.bw_scale, threads);
    let matern_nu = MaternNu::from_smoothness(params.smoothness);
    let apply = |d: f64| -> f64 {
        let r = d / bw;
        match params.kernel {
            Kernel::Gaussian => gaussian(r),
            Kernel::Laplace => laplace(r),
            Kernel::Matern => match matern_nu.unwrap_or(MaternNu::ThreeHalves) {
                MaternNu::Half => matern_half(r),
                MaternNu::ThreeHalves => matern_three_halves(r),
                MaternNu::FiveHalves => matern_five_halves(r),
            },
            // The energy and t kernels return above.
            Kernel::Energy | Kernel::T => -d,
        }
    };

    let mut k = vec![0.0; n * n];
    let pool = get_pool(threads);
    pool.install(|| {
        k.par_chunks_mut(n.max(1)).enumerate().for_each(|(j, col)| {
            for (i, slot) in col.iter_mut().enumerate() {
                *slot = apply(dist[j * n + i]);
            }
        });
    });
    k
}

/// The random-Fourier-feature t kernel on the standardized covariates.
///
/// Each frequency vector (a column of `t_proj`) defines a cosine and sine feature
/// per unit; the kernel is the average inner product of those features across the
/// draws, which is a Gram matrix and therefore positive semidefinite. The feature
/// build parallelizes over draws and the assembly over columns, both in a fixed
/// per-entry order, so the result is deterministic.
fn t_kernel(z: &[f64], n: usize, p: usize, params: &KernelParams<'_>, threads: usize) -> Vec<f64> {
    let d = params.n_draws;
    if n == 0 || d == 0 {
        return vec![0.0; n * n];
    }
    let proj = params.t_proj;
    let pool = get_pool(threads);

    // Cosine and sine of the projection angle for every unit and draw, laid out
    // draw-major: feature `[draw * n + i]` is the angle of unit `i` under draw.
    let mut cos_f = vec![0.0; n * d];
    let mut sin_f = vec![0.0; n * d];
    pool.install(|| {
        cos_f
            .par_chunks_mut(n)
            .zip(sin_f.par_chunks_mut(n))
            .enumerate()
            .for_each(|(draw, (cos_col, sin_col))| {
                for i in 0..n {
                    let mut angle = 0.0;
                    for k in 0..p {
                        angle += proj[draw * p + k] * z[k * n + i];
                    }
                    cos_col[i] = angle.cos();
                    sin_col[i] = angle.sin();
                }
            });
    });

    let inv = 1.0 / d as f64;
    let mut kernel = vec![0.0; n * n];
    pool.install(|| {
        kernel.par_chunks_mut(n).enumerate().for_each(|(j, col)| {
            for (i, slot) in col.iter_mut().enumerate() {
                let mut acc = 0.0;
                for draw in 0..d {
                    acc += cos_f[draw * n + i] * cos_f[draw * n + j]
                        + sin_f[draw * n + i] * sin_f[draw * n + j];
                }
                *slot = acc * inv;
            }
        });
    });
    kernel
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::{Mat, Side};

    fn params(kernel: Kernel) -> KernelParams<'static> {
        KernelParams {
            kernel,
            bw_scale: 1.0,
            smoothness: 1.5,
            t_proj: &[],
            n_draws: 0,
        }
    }

    /// The smallest eigenvalue of a symmetric column-major matrix.
    fn min_eigenvalue(k: &[f64], n: usize) -> f64 {
        let mat = Mat::from_fn(n, n, |i, j| k[j * n + i]);
        let eigen = mat
            .self_adjoint_eigen(Side::Lower)
            .expect("eigendecomposition");
        let values = eigen.S();
        (0..n).fold(f64::INFINITY, |m, i| m.min(values[i]))
    }

    fn assert_symmetric_unit_diagonal(k: &[f64], n: usize, diag: f64) {
        for i in 0..n {
            assert!(
                (k[i * n + i] - diag).abs() < 1e-9,
                "diagonal {}",
                k[i * n + i]
            );
            for j in 0..n {
                assert!((k[j * n + i] - k[i * n + j]).abs() < 1e-12, "asymmetry");
            }
        }
    }

    #[test]
    fn scalar_kernels_match_their_closed_forms() {
        // Hand values at r = 1.
        assert!((gaussian(1.0) - (-0.5f64).exp()).abs() < 1e-15);
        assert!((laplace(1.0) - (-1.0f64).exp()).abs() < 1e-15);
        assert!((matern_half(1.0) - (-1.0f64).exp()).abs() < 1e-15);
        // (1 + sqrt(3)) exp(-sqrt(3)).
        let m32 = (1.0 + SQRT3) * (-SQRT3).exp();
        assert!((matern_three_halves(1.0) - m32).abs() < 1e-15);
        // (1 + sqrt(5) + 5/3) exp(-sqrt(5)).
        let m52 = (1.0 + SQRT5 + 5.0 / 3.0) * (-SQRT5).exp();
        assert!((matern_five_halves(1.0) - m52).abs() < 1e-15);
        // Every kernel is one at zero radius.
        for f in [
            gaussian as fn(f64) -> f64,
            laplace,
            matern_half,
            matern_three_halves,
            matern_five_halves,
        ] {
            assert!((f(0.0) - 1.0).abs() < 1e-15);
        }
    }

    #[test]
    fn the_median_bandwidth_is_the_lower_triangle_median() {
        // Points 0, 1, 2, 3 on a line: lower-triangle distances are
        // 1, 2, 3, 1, 2, 1, whose sorted median of six values is (1 + 2) / 2.
        let x = [0.0, 1.0, 2.0, 3.0];
        let dist = pairwise::euclidean(&x, 4, 1, 1);
        let bw = median_bandwidth(&dist, 4, &[], 1.0, 1);
        assert!((bw - 1.5).abs() < 1e-12, "bandwidth {bw}");
        // The scale factor multiplies the median.
        let scaled = median_bandwidth(&dist, 4, &[], 2.0, 1);
        assert!((scaled - 3.0).abs() < 1e-12, "scaled bandwidth {scaled}");
    }

    #[test]
    fn discarded_units_do_not_enter_the_bandwidth() {
        // Discarding the far point drops the large distances, lowering the median.
        let x = [0.0, 1.0, 2.0, 100.0];
        let dist = pairwise::euclidean(&x, 4, 1, 1);
        let full = median_bandwidth(&dist, 4, &[false, false, false, false], 1.0, 1);
        let dropped = median_bandwidth(&dist, 4, &[false, false, false, true], 1.0, 1);
        assert!(dropped < full, "dropped {dropped} !< full {full}");
        // Remaining pairs 0-1, 0-2, 1-2 have distances 1, 2, 1; median 1.
        assert!((dropped - 1.0).abs() < 1e-12, "dropped bandwidth {dropped}");
    }

    #[test]
    fn the_energy_kernel_is_the_negative_distance() {
        let covs = [0.0, 1.0, 2.0, 0.5, -0.5, 1.5];
        let n = 3;
        let s = vec![1.0; n];
        let k = build_kernel(&covs, n, 2, &s, &params(Kernel::Energy), &[], 1);
        let z = transform::transform(&covs, n, 2, Distance::ScaledEuclidean, &s, 1);
        let dist = pairwise::euclidean(&z, n, 2, 1);
        for (a, b) in k.iter().zip(&dist) {
            assert!((a + b).abs() < 1e-12, "energy kernel is not -distance");
        }
    }

    #[test]
    fn the_gaussian_kernel_is_symmetric_positive_semidefinite() {
        let n = 6;
        let covs: Vec<f64> = (0..n * 2).map(|k| ((k as f64) * 0.37).sin()).collect();
        let s = vec![1.0; n];
        let k = build_kernel(&covs, n, 2, &s, &params(Kernel::Gaussian), &[], 1);
        assert_symmetric_unit_diagonal(&k, n, 1.0);
        assert!(min_eigenvalue(&k, n) > -1e-9, "gaussian is not psd");
    }

    #[test]
    fn every_matern_smoothness_is_positive_semidefinite() {
        let n = 6;
        let covs: Vec<f64> = (0..n * 2).map(|k| ((k as f64) * 0.29).cos()).collect();
        let s = vec![1.0; n];
        for smoothness in [0.5, 1.5, 2.5] {
            let p = KernelParams {
                kernel: Kernel::Matern,
                bw_scale: 1.0,
                smoothness,
                t_proj: &[],
                n_draws: 0,
            };
            let k = build_kernel(&covs, n, 2, &s, &p, &[], 1);
            assert_symmetric_unit_diagonal(&k, n, 1.0);
            assert!(
                min_eigenvalue(&k, n) > -1e-9,
                "matern {smoothness} is not psd"
            );
        }
    }

    #[test]
    fn the_gaussian_kernel_uses_the_closed_form_with_the_median_bandwidth() {
        // One-dimensional standardized points let the entries be checked against
        // the closed form directly, confirming the bandwidth is wired in.
        let covs = [0.0, 1.0, 2.0, 3.0];
        let n = 4;
        let s = vec![1.0; n];
        let z = transform::transform(&covs, n, 1, Distance::ScaledEuclidean, &s, 1);
        let dist = pairwise::euclidean(&z, n, 1, 1);
        let bw = median_bandwidth(&dist, n, &[], 1.0, 1);
        let k = build_kernel(&covs, n, 1, &s, &params(Kernel::Gaussian), &[], 1);
        for j in 0..n {
            for i in 0..n {
                let expected = gaussian(dist[j * n + i] / bw);
                assert!((k[j * n + i] - expected).abs() < 1e-12, "gaussian entry");
            }
        }
    }

    #[test]
    fn the_t_kernel_is_symmetric_positive_semidefinite_and_reproducible() {
        // A fixed projection matrix stands in for the R-side draws, so the kernel
        // is deterministic and can be checked entrywise.
        let n = 5;
        let p = 2;
        let covs: Vec<f64> = (0..n * p)
            .map(|k| ((k as f64) * 0.41).sin() + 0.2)
            .collect();
        let s = vec![1.0; n];
        let d = 3;
        let proj = vec![0.3, -0.7, 1.1, 0.4, -0.2, 0.9]; // p by d column-major
        let tp = KernelParams {
            kernel: Kernel::T,
            bw_scale: 1.0,
            smoothness: 1.5,
            t_proj: &proj,
            n_draws: d,
        };
        let k = build_kernel(&covs, n, p, &s, &tp, &[], 1);
        assert_symmetric_unit_diagonal(&k, n, 1.0);
        assert!(min_eigenvalue(&k, n) > -1e-9, "t kernel is not psd");

        // A hand computation of one entry from the cosine-and-sine features.
        let z = transform::transform(&covs, n, p, Distance::ScaledEuclidean, &s, 1);
        let angle = |i: usize, draw: usize| -> f64 {
            (0..p).map(|k| proj[draw * p + k] * z[k * n + i]).sum()
        };
        let mut expected = 0.0;
        for draw in 0..d {
            let a0 = angle(0, draw);
            let a1 = angle(1, draw);
            expected += (a0 - a1).cos();
        }
        expected /= d as f64;
        assert!((k[1 * n + 0] - expected).abs() < 1e-12, "t kernel entry");

        // Rebuilding with the same projections gives a bit-identical matrix.
        let again = build_kernel(&covs, n, p, &s, &tp, &[], 1);
        for (a, b) in k.iter().zip(&again) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }

    #[test]
    fn the_kernel_is_invariant_to_a_uniform_covariate_rescaling() {
        // Standardizing first makes the gaussian kernel identical after scaling
        // every covariate by a constant.
        let n = 6;
        let covs: Vec<f64> = (0..n * 2).map(|k| ((k as f64) * 0.53).sin()).collect();
        let scaled: Vec<f64> = covs.iter().map(|v| v * 100.0).collect();
        let s = vec![1.0; n];
        let raw = build_kernel(&covs, n, 2, &s, &params(Kernel::Gaussian), &[], 1);
        let big = build_kernel(&scaled, n, 2, &s, &params(Kernel::Gaussian), &[], 1);
        for (a, b) in raw.iter().zip(&big) {
            assert!((a - b).abs() < 1e-10, "kernel is not scale invariant");
        }
    }

    #[test]
    fn the_kernel_build_is_deterministic_across_thread_counts() {
        let n = 40;
        let covs: Vec<f64> = (0..n * 3).map(|k| ((k as f64) * 0.19).sin()).collect();
        let s = vec![1.0; n];
        let one = build_kernel(&covs, n, 3, &s, &params(Kernel::Gaussian), &[], 1);
        for threads in [2, 4] {
            let many = build_kernel(&covs, n, 3, &s, &params(Kernel::Gaussian), &[], threads);
            for (a, b) in one.iter().zip(&many) {
                assert_eq!(a.to_bits(), b.to_bits());
            }
        }
    }
}
