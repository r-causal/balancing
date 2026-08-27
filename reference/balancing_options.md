# Package options

balancing reads a handful of global options, set with
[`options()`](https://rdrr.io/r/base/options.html), that tune how a fit
runs without changing any function's signature.

## Details

- `balancing.quiet`: when `TRUE`, suppresses the informational alerts
  [`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
  prints, such as the detected exposure type. Defaults to `FALSE`.

- `balancing.threads`: the number of worker threads the Rust core may
  use. When unset, the count is resolved automatically from the physical
  core count, capped by `OMP_THREAD_LIMIT` and `OMP_NUM_THREADS`, and
  forced to two under `R CMD check`.

- `balancing.entropy_solver`: the solver for the exact entropy problem,
  one of `"newton"` (the default), `"lbfgs"`, or `"lbfgs_then_newton"`.
  Newton is the only solver that drives the estimating equations to
  machine precision; the alternatives trade some precision for speed on
  large problems. When the option is unset and the Newton solve stops
  short of its tolerance or returns weights that are not finite, the fit
  retries once with `"lbfgs_then_newton"`, whose Newton polish restores
  that precision, and announces the retry. Setting the option pins the
  solver and disables the retry, so a fit pinned to `"newton"` reports
  the failed solve as it stands.

- `balancing.qp_backend`: the quadratic-program backend for the CFD and
  stable-balancing-weights methods, one of `"auto"` (the default),
  `"osqp"`, or `"clarabel"`. Under `"auto"`, the default solver runs
  first and the fit re-solves with the interior-point backend on a
  primal-infeasibility certificate; `"osqp"` disables that fallback and
  `"clarabel"` uses the interior-point backend directly. Energy
  balancing and the CFD energy kernel assemble an indefinite quadratic
  form, which the interior-point backend refuses, so they always solve
  through `"osqp"` and announce a `"clarabel"` pin as ignored.

Results are deterministic across thread counts: the same inputs produce
identical weights at any `balancing.threads` value.

## Examples

``` r
# Silence the exposure-type alert for a single fit.
withr::with_options(list(balancing.quiet = TRUE), {
  n <- 100
  x1 <- rnorm(n)
  df <- data.frame(exposure = rbinom(n, 1, plogis(x1)), x1 = x1)
  balance(df, exposure, x1, method = bw_entropy())
})
#> 
#> ── Entropy balancing ───────────────────────────────────────────────────────────
#> Exposure: "exposure" (binary)
#> Estimand: "ate"
#> Observations: 100
#> Solver: converged in 4 iterations
#> Constraints: 1 term (tolerance 0)
#> Largest imbalance: 0.0000 (standardized mean difference)
```
