//! Thread pool management and deterministic reductions.
//!
//! Two guarantees drive this module. First, the core never installs a rayon
//! global pool: callers pass an explicit thread count resolved on the R side,
//! and every parallel region runs inside a cached [`rayon::ThreadPool`] of that
//! size. Second, every floating-point reduction is bit-for-bit identical across
//! thread counts on the same binary and inputs, which [`deterministic_map_reduce`]
//! provides by fixing the reduction order independent of how work is scheduled.

use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock};

use rayon::ThreadPool;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

/// Number of indices folded together before partial results are combined.
///
/// Determinism rests on this constant. Indices are partitioned into contiguous
/// chunks of this size; each chunk is folded sequentially in index order, and
/// the chunk partials are combined sequentially in chunk order. The partition
/// does not depend on the thread count, so the summation tree, and therefore the
/// exact floating-point result, is fixed regardless of how many threads run.
/// Changing this value perturbs results at the last bits and is a documented
/// behavior change.
pub const REDUCE_CHUNK: usize = 4096;

/// Cache of thread pools keyed by size, so repeated solver calls reuse workers.
fn pool_cache() -> &'static Mutex<HashMap<usize, Arc<ThreadPool>>> {
    static CACHE: OnceLock<Mutex<HashMap<usize, Arc<ThreadPool>>>> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(HashMap::new()))
}

/// Return a cached thread pool with `threads` workers, building it on first use.
///
/// A request for zero threads is treated as one. The pool count is clamped to at
/// least one so callers never have to special-case it.
pub fn get_pool(threads: usize) -> Arc<ThreadPool> {
    let threads = threads.max(1);
    let mut cache = pool_cache().lock().expect("thread pool cache poisoned");
    if let Some(pool) = cache.get(&threads) {
        return Arc::clone(pool);
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("failed to build thread pool");
    let pool = Arc::new(pool);
    cache.insert(threads, Arc::clone(&pool));
    pool
}

/// Fold `0..n` into a single accumulator with a fixed, thread-count-independent
/// reduction order.
///
/// `identity` builds a fresh zero accumulator, `fold` accumulates a single index
/// into a partial, and `combine` merges one partial into another. Chunks of
/// [`REDUCE_CHUNK`] indices are folded in parallel inside `pool`, but each chunk
/// is folded in ascending index order and the partials are merged in ascending
/// chunk order, so the result is identical for any pool size.
pub fn deterministic_map_reduce<T, Identity, Fold, Combine>(
    pool: &ThreadPool,
    n: usize,
    identity: Identity,
    fold: Fold,
    combine: Combine,
) -> T
where
    T: Send,
    Identity: Fn() -> T + Sync,
    Fold: Fn(&mut T, usize) + Sync,
    Combine: Fn(&mut T, T) + Sync,
{
    if n == 0 {
        return identity();
    }
    let n_chunks = n.div_ceil(REDUCE_CHUNK);
    // `into_par_iter().collect()` preserves chunk order, so `partials[c]` is the
    // fold of chunk `c` and the sequential merge below runs in chunk order.
    let partials: Vec<T> = pool.install(|| {
        (0..n_chunks)
            .into_par_iter()
            .map(|chunk| {
                let start = chunk * REDUCE_CHUNK;
                let end = (start + REDUCE_CHUNK).min(n);
                let mut acc = identity();
                for i in start..end {
                    fold(&mut acc, i);
                }
                acc
            })
            .collect()
    });
    let mut iter = partials.into_iter();
    let mut acc = iter.next().expect("at least one chunk exists when n > 0");
    for partial in iter {
        combine(&mut acc, partial);
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    // A reduction whose value is independent of association only up to rounding;
    // mixing large and small magnitudes makes a thread-count-dependent order
    // visible if one ever crept in.
    fn sample(i: usize) -> f64 {
        let x = i as f64;
        (x * 0.5).sin() * 1e6 + (x * 1e-3).cos() * 1e-9
    }

    fn sum_with(threads: usize, n: usize) -> f64 {
        let pool = get_pool(threads);
        deterministic_map_reduce(
            &pool,
            n,
            || 0.0_f64,
            |acc, i| *acc += sample(i),
            |acc, partial| *acc += partial,
        )
    }

    #[test]
    fn reduction_is_bit_identical_across_thread_counts() {
        let n = 4096 * 5 + 37;
        let one = sum_with(1, n).to_bits();
        for threads in [2, 3, 4, 8] {
            assert_eq!(sum_with(threads, n).to_bits(), one);
        }
    }

    #[test]
    fn reduction_matches_sequential_index_order() {
        let n = 4096 * 3 + 11;
        let mut expected = 0.0_f64;
        for chunk_start in (0..n).step_by(REDUCE_CHUNK) {
            let mut partial = 0.0_f64;
            for i in chunk_start..(chunk_start + REDUCE_CHUNK).min(n) {
                partial += sample(i);
            }
            expected += partial;
        }
        assert_eq!(sum_with(4, n).to_bits(), expected.to_bits());
    }

    #[test]
    fn empty_range_returns_identity() {
        assert_eq!(sum_with(4, 0), 0.0);
    }

    #[test]
    fn pools_are_cached_by_size() {
        let a = get_pool(2);
        let b = get_pool(2);
        assert!(Arc::ptr_eq(&a, &b));
    }
}
