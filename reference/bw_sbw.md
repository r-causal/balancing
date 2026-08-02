# Stable balancing weights

`bw_sbw()` specifies stable balancing weights for
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md).
Among all weightings that hold each reweighted exposure group's
covariate means inside a tolerance band, stable balancing weights select
the one of least dispersion, following Zubizarreta. The default `"l2"`
norm minimizes the sum of squared weights, so for a fixed per-group
total it minimizes the weight variance. Stable balancing weights support
binary, categorical, and continuous exposures.

## Usage

``` r
bw_sbw(
  ...,
  norm = c("l2", "l1", "linf"),
  min_weight = 1e-08,
  convergence_tolerance = NULL,
  max_iterations = NULL
)
```

## Arguments

- ...:

  Reserved for future extensions; must be empty. Tuning parameters must
  be passed by name.

- norm:

  The weight-dispersion norm to minimize, one of `"l2"` (the sum of
  squared weights, minimum variance), `"l1"` (the sum of absolute
  deviations from one), or `"linf"` (the largest absolute deviation from
  one).

- min_weight:

  The smallest permitted weight.

- convergence_tolerance:

  The quadratic-program solver tolerance, or `NULL` for the core
  default.

- max_iterations:

  The maximum solver iterations, or `NULL` for the core default.

## Value

An `bw_sbw` specification, a
[balance_method](https://r-causal.github.io/balancing/reference/balance_method.md).

## Details

The balance tolerance is the method's central tuning parameter. It is
set through `constraints = balance_terms(tolerance = ...)` in
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md),
not on the method spec, and it must be positive: with an exact (zero)
tolerance the problem reduces to exact moment balance, which abandons
the minimum-variance rationale and is prone to infeasibility, so a fit
without a positive tolerance is refused. Given a positive tolerance the
weights hold each group's weighted covariate means within the band while
minimizing the weight dispersion. For a discrete exposure no feasible
reweighting has smaller dispersion than the fit returns.

A single scalar tolerance applies to every covariate; a named vector
sets a tolerance per covariate. A covariate left at zero in a named
vector demands exact balance on that covariate while the others are
relaxed, the infeasibility-prone case, so a mixed specification is an
explicit choice rather than a convenience.

For a discrete exposure the constraints default to first-moment balance;
pass
[`balance_terms()`](https://r-causal.github.io/balancing/reference/balance_terms.md)
to balance higher moments, interactions, or quantiles, each inside its
tolerance band. For a focal estimand the focal group keeps unit weight
and the other groups are pulled to the focal group's covariate means.

For a continuous exposure the weighted exposure-covariate correlations
are held within the tolerance. The quadratic program bounds a linearized
correlation whose scales are fixed at the sample, so the fit tightens
that internal bound over a few passes until the reported weighted
correlation sits inside the requested band. The returned weights
therefore minimize dispersion over the tightened internal band rather
than over every weighting that meets the reported band, so the strict
minimum-dispersion guarantee is stated for discrete exposures only.

Stable balancing weights belong to the quadratic-program family, which
has no estimating equations, so a fit produces no estimating-equations
container.

The `norm` argument selects how the weight dispersion is measured,
always against the uniform baseline of one within each reweighted group.
`"l2"` minimizes the sum of squared weights, so for a fixed per-group
total it minimizes the weight variance. `"l1"` minimizes the sum of
absolute deviations from one, which tends to leave many weights
untouched and concentrate the reweighting on a few units. `"linf"`
minimizes the single largest absolute deviation from one, which spreads
the reweighting as evenly as the balance constraints allow. The `"l1"`
and `"linf"` problems are linear programs solved through the same
quadratic-program backends as `"l2"`; their solutions can be non-unique,
so a fit reports the achieved dispersion rather than promising a unique
weighting.

## References

Zubizarreta, J. R. (2015). Stable weights that balance covariates for
estimation with incomplete outcome data. *Journal of the American
Statistical Association*, 110(511), 910-922.

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
fit <- balance(
  df,
  exposure,
  c(x1, x2),
  method = bw_sbw(),
  constraints = balance_terms(tolerance = 0.05)
)
#> ℹ Treating `.exposure` as binary.
fit
#> 
#> ── Stable balancing weights ────────────────────────────────────────────────────
#> Exposure: "exposure" (binary)
#> Estimand: "ate"
#> Observations: 200
#> Solver: converged in 75 iterations
#> Constraints: 2 terms (tolerance 0.05)
#> Largest imbalance: 0.1000 (standardized mean difference)
```
