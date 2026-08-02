# Energy balancing

`bw_energy()` specifies energy balancing for
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md).
The weights minimize the energy distance between the reweighted exposure
groups and a target sample, subject to a simplex-type constraint set, so
the reweighting improves multivariate covariate balance without positing
a propensity model. Energy balancing supports binary, categorical, and
continuous exposures.

## Usage

``` r
bw_energy(
  ...,
  distance = c("scaled_euclidean", "mahalanobis", "euclidean"),
  improved = TRUE,
  weight_penalty = 1e-04,
  min_weight = 1e-08,
  distribution_moments = NULL,
  dimension_adjustment = TRUE,
  convergence_tolerance = NULL,
  max_iterations = NULL
)
```

## Arguments

- ...:

  Reserved for future extensions; must be empty. Tuning parameters must
  be passed by name.

- distance:

  The covariate distance definition the energy objective is built on,
  one of `"scaled_euclidean"` (each covariate divided by its standard
  deviation), `"mahalanobis"`, or `"euclidean"`.

- improved:

  Whether to add the between-group energy distance of the improved
  variant for the average treatment effect with a discrete exposure.

- weight_penalty:

  The L2 penalty on the weights, which stabilizes the quadratic program.

- min_weight:

  The smallest permitted weight. The reported weights average one within
  each exposure group, so a floor approaching one leaves almost no room
  above it: the weight spread shrinks in proportion to the headroom
  `1 - min_weight`, and the fit degenerates smoothly into uniform
  weights and reports the balance uniform weights achieve. Nothing warns
  at that boundary, because the problem stays feasible and the solution
  is a real one.
  [`bw_sbw()`](https://r-causal.github.io/balancing/reference/bw_sbw.md),
  whose tolerances are hard constraints rather than an objective,
  refuses the same floor as infeasible instead.

- distribution_moments:

  For a continuous exposure, the number of exposure and covariate
  marginal moments held equal to the sample under the base measure, or
  `NULL` for the constraint moments. Raised automatically when smaller
  than the constraint moments. Energy balancing carries no base weights,
  so the base measure is the sampling weights, and without them the
  marginals are held equal to the unweighted sample.

- dimension_adjustment:

  For a continuous exposure, whether to weight the covariate energy
  distance by the covariate dimensionality adjustment.

- convergence_tolerance:

  The quadratic-program solver tolerance, or `NULL` for the core
  default.

- max_iterations:

  The maximum solver iterations, or `NULL` for the core default.

## Value

An `bw_energy` specification, a
[balance_method](https://r-causal.github.io/balancing/reference/balance_method.md).

## Details

For a binary or categorical exposure the objective is the sum of each
group's energy distance to the target sample. The improved variant for
the average treatment effect adds the between-group energy distance,
which balances the groups against one another as well as against the
sample. A focal estimand reweights the non-focal groups toward the focal
group, whose units keep their base weight. The energy distance is built
from a pairwise covariate distance matrix; `distance` selects how that
matrix is formed.

For a continuous exposure the objective is the weighted distance
covariance between the exposure and the covariates, following Huling,
Greifer, and Chen, plus the marginal energy distances of the weighted
exposure and covariate distributions. `distribution_moments` sets how
many exposure and covariate marginal moments are held equal to the
sample under the base measure, and `dimension_adjustment` reweights the
covariate energy distance by the covariate dimensionality.

Energy balancing belongs to the quadratic-program family, which has no
estimating equations, so a fit produces no estimating-equations
container and the tolerance in
[`balance_terms()`](https://r-causal.github.io/balancing/reference/balance_terms.md)
relaxes any added moment constraints rather than selecting an inexact
solver. A tolerance supplied without moment constraints has nothing to
relax, so it is warned and ignored.

## References

Huling, J. D. and Mak, S. (2024). Energy balancing of covariate
distributions. *Journal of Causal Inference*, 12(1), 20220029.

Huling, J. D., Greifer, N., and Chen, G. (2024). Independence weights
for causal inference with continuous treatments. *Journal of the
American Statistical Association*, 119(546), 1657-1670.

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
fit <- balance(df, exposure, c(x1, x2), method = bw_energy())
#> ℹ Treating `.exposure` as binary.
fit
#> 
#> ── Energy balancing ────────────────────────────────────────────────────────────
#> Exposure: "exposure" (binary)
#> Estimand: "ate"
#> Observations: 200
#> Solver: converged in 75 iterations
#> Constraints: 2 terms (tolerance 0)
#> Largest imbalance: 0.0058 (standardized mean difference)
```
