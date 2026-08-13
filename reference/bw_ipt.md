# Inverse probability tilting

`bw_ipt()` specifies inverse probability tilting for
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md).
A propensity model is fit not by maximum likelihood but by a tilted
moment condition that forces each treatment group's weighted covariate
means to their estimand targets, so balance on the requested moments is
exact by construction. Inverse probability tilting supports binary and
categorical exposures.

## Usage

``` r
bw_ipt(
  ...,
  link = c("logit", "probit", "cloglog"),
  convergence_tolerance = 1e-10,
  max_iterations = NULL
)
```

## Arguments

- ...:

  Reserved for future extensions; must be empty. Tuning parameters must
  be passed by name.

- link:

  The propensity link, one of `"logit"`, `"probit"`, or `"cloglog"`.

- convergence_tolerance:

  The solver convergence tolerance on the tilting moment.

- max_iterations:

  The maximum solver iterations, or `NULL` for the core default.

## Value

An `bw_ipt` specification, a
[balance_method](https://r-causal.github.io/balancing/reference/balance_method.md).

## Details

For each treatment level the propensity `p_i = G(x_i' beta)` is
estimated so that the weighted covariate total matches the target
population total. The average treatment effect tilts every level to the
whole sample, weighting a unit by the inverse of its modeled propensity.
A focal estimand tilts each non-focal level to the focal level,
weighting a unit by `(1 - p_i) / p_i`, and leaves the focal units at
weight one. With mean balance and the logit link the tilt solves the
same treated-target problem as entropy balancing, so the two methods
produce the same average-treatment-effect-on-the-treated weights for a
binary exposure.

The weights solve smooth estimating equations regardless of the
requested tolerance, which
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
records for the M-estimation variance in
[`ipw()`](https://r-causal.github.io/balancing/reference/ipw.balancing.md).

## References

Graham, B. S., Pinto, C. C. de X., and Egel, D. (2012). Inverse
probability tilting for moment condition models with missing data. *The
Review of Economic Studies*, 79(3), 1053-1079.

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
fit <- balance(df, exposure, c(x1, x2), method = bw_ipt())
#> ℹ Treating `.exposure` as binary
fit
#> 
#> ── Inverse probability tilting ─────────────────────────────────────────────────
#> Exposure: "exposure" (binary)
#> Estimand: "ate"
#> Observations: 200
#> Solver: converged in 3 iterations
#> Constraints: 2 terms (tolerance 0)
#> Largest imbalance: 0.0000 (standardized mean difference)
```
