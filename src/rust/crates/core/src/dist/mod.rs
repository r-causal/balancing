//! Distance matrices for the quadratic-program balancing family.
//!
//! Energy balancing measures covariate balance through the energy distance,
//! which is built from the pairwise distances of the covariate rows. Two steps
//! produce the matrix: [`transform`] rescales the covariates so that ordinary
//! Euclidean distance on the transformed rows reproduces the requested distance
//! (scaled Euclidean, Mahalanobis, or plain Euclidean), and [`pairwise`] forms
//! the dense symmetric distance matrix in parallel.

pub mod kernels;
pub mod pairwise;
pub mod transform;

/// The distance definition the energy objective is built on.
///
/// Each variant names a transform applied to the covariates before Euclidean
/// distances are taken. `ScaledEuclidean` centers each column at its weighted
/// mean and divides by its weighted standard deviation; `Mahalanobis` centers
/// and whitens by the weighted covariance; `Euclidean` uses the covariates
/// unchanged. Centering shifts every row alike and so leaves the distances
/// themselves untouched.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Distance {
    /// Euclidean distance on columns centered at their weighted mean and
    /// divided by their weighted standard deviation.
    ScaledEuclidean,
    /// Euclidean distance after centering and whitening by the weighted
    /// covariance.
    Mahalanobis,
    /// Euclidean distance on the covariates as supplied.
    Euclidean,
}

impl Distance {
    /// Resolve the distance named by the R layer, returning `None` for an
    /// unrecognized name.
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "scaled_euclidean" => Some(Distance::ScaledEuclidean),
            "mahalanobis" => Some(Distance::Mahalanobis),
            "euclidean" => Some(Distance::Euclidean),
            _ => None,
        }
    }
}

/// Build the `n` by `n` distance matrix for `covs` under `distance`.
///
/// `covs` is a column-major `n` by `p` matrix; `w` are the per-unit weights used
/// to compute the standardizing statistics (the sampling weights). The result is
/// a column-major `n` by `n` symmetric matrix with a zero diagonal. `threads`
/// sizes the pool for the pairwise assembly.
pub fn distance_matrix(
    covs: &[f64],
    n: usize,
    p: usize,
    distance: Distance,
    w: &[f64],
    threads: usize,
) -> Vec<f64> {
    let transformed = transform::transform(covs, n, p, distance, w, threads);
    let tp = transformed.len() / n.max(1);
    pairwise::euclidean(&transformed, n, tp, threads)
}
