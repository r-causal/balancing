# Entropy balancing

`bw_entropy()` specifies entropy balancing for
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md).
The weights minimize the Kullback-Leibler divergence from a set of base
weights subject to the covariate constraints, so among all reweightings
that achieve balance the solution stays as close as possible to the base
weights. Entropy balancing supports binary, categorical, and continuous
exposures.

## Usage

``` r
bw_entropy(
  ...,
  base_weights = NULL,
  distribution_moments = NULL,
  convergence_tolerance = 1e-10,
  max_iterations = NULL
)
```

## Arguments

- ...:

  Reserved for future extensions; must be empty. Tuning parameters must
  be passed by name.

- base_weights:

  A numeric vector of base weights, one per observation, or `NULL` for
  uniform base weights. The estimated weights minimize
  `sum(w * log(w / base_weights))`.

- distribution_moments:

  For continuous exposures, the number of exposure and covariate
  marginal moments held equal to the sample under the base measure, or
  `NULL` for the constraint moments. Raised automatically when smaller
  than the constraint moments. The base measure is the product of the
  sampling weights and `base_weights`, so without either the marginals
  are held equal to the unweighted sample.

- convergence_tolerance:

  The solver convergence tolerance. What it measures depends on which
  problem is solved. The exact problem, chosen when every tolerance in
  [`balance_terms()`](https://r-causal.github.io/balancing/reference/balance_terms.md)
  is zero, measures the gradient sup norm and responds to this value
  across its range. A positive tolerance selects the inexact problem,
  solved by FISTA against the relative change in the loss; that
  criterion is the weaker of the two, so the value is tightened to at
  most `1e-14` to hold the achieved balance inside the requested box,
  and anything above `1e-14` is inert there.

- max_iterations:

  The maximum solver iterations, or `NULL` for the core default.

## Value

An `bw_entropy` specification, a
[balance_method](https://r-causal.github.io/balancing/reference/balance_method.md).

## Details

For a binary exposure the average treatment effect reweights each
exposure group to the pooled covariate means, and the average treatment
effect on the treated reweights the control group to the treated
covariate means while the treated group is left unreweighted. Its
reported weights are its base weights carried to the group's
sampling-weighted total rather than the base weights themselves, so they
are proportional to the base weights and constant base weights come back
as ones whatever level they were set at. When every requested tolerance
is zero the constraints hold exactly and the weights solve smooth
estimating equations, which
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
records for the M-estimation variance in
[`ipw()`](https://r-causal.github.io/balancing/reference/ipw.balancing.md).
A positive `tolerance` in
[`balance_terms()`](https://r-causal.github.io/balancing/reference/balance_terms.md)
selects the inexact problem, which balances each constraint to within
the tolerance and does not produce estimating equations.

## References

Hainmueller, J. (2012). Entropy balancing for causal effects: A
multivariate reweighting method to produce balanced samples in
observational studies. *Political Analysis*, 20(1), 25-46.

## Examples

``` r
n <- 200
x1 <- rnorm(n)
x2 <- rnorm(n)
df <- data.frame(
  exposure = rbinom(n, 1, plogis(0.5 * x1 - 0.5 * x2)),
  x1 = x1,
  x2 = x2
)
fit <- balance(df, exposure, c(x1, x2), method = bw_entropy())
#> ℹ Treating `.exposure` as binary.
fit
#> 
#> ── Entropy balancing ───────────────────────────────────────────────────────────
#> Exposure: "exposure" (binary)
#> Estimand: "ate"
#> Observations: 200
#> Solver: converged in 4 iterations
#> Constraints: 2 terms (tolerance 0)
#> Largest imbalance: 0.0000 (standardized mean difference)
```
