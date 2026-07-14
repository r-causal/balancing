//! Pure computation core for the balancing package.
//!
//! This crate carries no R dependency. Unit tests, benchmarks, and the solver
//! comparison harness all run without an R toolchain, which keeps the numerical
//! work testable in isolation and leaves room for a second frontend against the
//! same core.

pub mod dist;
pub mod esteq;
pub mod glm;
pub mod linalg;
pub mod links;
pub mod methods;
pub mod qp;
pub mod threads;
pub mod weights;

/// Threads available for parallel work, and the constraint that set the count.
///
/// The count reflects the physical parallelism the process observes, lowered to
/// `OMP_THREAD_LIMIT` when that environment variable imposes a smaller cap. The
/// second element names the source: `"OMP_THREAD_LIMIT"` when that cap applied,
/// otherwise `"system"`. Automatic thread resolution for solver calls happens on
/// the R side and is passed down explicitly; this reports what the core sees so
/// the boundary can be exercised before any solver exists.
pub fn available_threads() -> (usize, &'static str) {
    let system = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let omp_limit = std::env::var("OMP_THREAD_LIMIT")
        .ok()
        .and_then(|value| value.parse::<usize>().ok());
    resolve_threads(system, omp_limit)
}

/// Resolve the effective thread count from the system parallelism and an
/// optional `OMP_THREAD_LIMIT`. Split out from environment access so the
/// resolution logic can be tested deterministically.
fn resolve_threads(system: usize, omp_limit: Option<usize>) -> (usize, &'static str) {
    let system = system.max(1);
    match omp_limit {
        Some(limit) if limit >= 1 && limit < system => (limit, "OMP_THREAD_LIMIT"),
        _ => (system, "system"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn available_threads_reports_at_least_one() {
        let (count, _source) = available_threads();
        assert!(count >= 1);
    }

    #[test]
    fn omp_limit_below_system_caps_the_count() {
        assert_eq!(resolve_threads(8, Some(2)), (2, "OMP_THREAD_LIMIT"));
    }

    #[test]
    fn omp_limit_at_or_above_system_is_ignored() {
        assert_eq!(resolve_threads(4, Some(4)), (4, "system"));
        assert_eq!(resolve_threads(4, Some(16)), (4, "system"));
    }

    #[test]
    fn missing_omp_limit_uses_system() {
        assert_eq!(resolve_threads(6, None), (6, "system"));
    }

    #[test]
    fn system_is_floored_at_one() {
        assert_eq!(resolve_threads(0, None), (1, "system"));
    }
}
