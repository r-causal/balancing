# balancing

balancing calculates optimization-based balancing weights for causal
inference. Rather than modeling the probability of exposure and then
hoping the resulting weights balance the covariates, these methods make
covariate balance the objective of a convex optimization problem and
solve for the weights that achieve it directly. The package covers six
methods (entropy balancing, inverse probability tilting, the covariate
balancing propensity score, energy balancing, characteristic function
distance balancing, and stable balancing weights) across binary,
categorical, and continuous exposures and a range of estimands. A Rust
core provides the numerical solvers.

## Installation

You can install the development version of balancing from
[r-causal.r-universe.dev](https://r-causal.r-universe.dev/) with:

``` r

install.packages(
  "balancing",
  repos = c("https://r-causal.r-universe.dev", getOption("repos"))
)
```

You can also install the development version of balancing from source
from [GitHub](https://github.com/r-causal/balancing) with:

``` r

# install.packages("pak")
pak::pak("r-causal/balancing")
```

Installing from source requires a Rust toolchain (`rustc` 1.88 or newer
and Cargo); we recommend installing one with
[rustup](https://rustup.rs/).

## Usage

You give
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
a data frame, name the exposure and covariates with tidyselect, choose a
method and an estimand, and it returns a fitted object carrying the
weights and a balance table.

``` r

library(balancing)

# Simulate data with two confounders and a binary exposure
set.seed(1)
n <- 500
x1 <- rnorm(n)
x2 <- rnorm(n)
z <- rbinom(n, 1, plogis(0.4 * x1 - 0.5 * x2))
y <- 1 + 0.9 * z + 0.6 * x1 - 0.4 * x2 + rnorm(n)
study <- data.frame(exposure = z, age = x1, score = x2, outcome = y)

# Fit entropy balancing weights for the ATT
fit <- balance(
  study,
  exposure,
  c(age, score),
  method = bw_entropy(),
  estimand = "att"
)

fit
#> 
#> ── Entropy balancing ───────────────────────────────────────────────────────────
#> Exposure: "exposure" (binary)
#> Estimand: "att" (focal level "1")
#> Observations: 500
#> Solver: converged in 3 iterations
#> Constraints: 2 terms (tolerance 0)
#> Largest imbalance: 0.0000 (standardized mean difference)
```

The printed summary reports the solver status and the largest imbalance
the weights leave behind on the constraint terms.

The weights are a `bw` vector, a sibling of
[`propensity::psw()`](https://r-causal.github.io/propensity/reference/psw.html).
Pass them to a weighted outcome model to estimate the effect.

``` r

study$w <- weights(fit)
outcome_mod <- lm(outcome ~ exposure, data = study, weights = w)
coef(outcome_mod)[["exposure"]]
#> [1] 1.131677
```

For the estimating-equation methods (entropy balancing, inverse
probability tilting, and the covariate balancing propensity score), pass
the fit and the weighted outcome model to
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
to get effect estimates with standard errors that account for having
estimated the weights. Binary, categorical, and continuous exposures are
supported; with a continuous exposure the reported effect is the
exposure coefficient of a weighted marginal structural model.

``` r

ipw(fit, outcome_mod)
#> Inverse Probability Weight Estimator
#> Estimand: ATT 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: balance(.data = study, .exposure = exposure, .covariates = c(age, 
#>     score), method = bw_entropy(), estimand = "att") 
#> 
#> Outcome Model:
#>   Call: lm(formula = outcome ~ exposure, data = study, weights = w) 
#> 
#> Marginal estimates:
#>      estimate std.err      z ci.lower ci.upper conf.level   p.value    
#> diff  1.13168 0.12345 9.1673  0.88972   1.3736       0.95 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## How balancing relates to the other r-causal packages

balancing is part of the [r-causal](https://github.com/r-causal) family
and shares its design language with three sibling packages.

- [propensity](https://r-causal.github.io/propensity/) calculates
  inverse probability weights from a fitted propensity score model and
  estimates effects with
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html).
  balancing produces a different kind of weight, one that targets
  balance directly rather than through a modeled score, but it plugs
  into the same effect-estimation workflow: the weights share the
  `causal_wts` vocabulary, and
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  dispatches on a balancing fit just as it does on a propensity model.
- [halfmoon](https://r-causal.github.io/halfmoon/) assesses covariate
  balance. balancing reports the balance achieved on its own constraint
  terms, but a full balance-assessment workflow, including variables the
  weights were not asked to balance, lives in halfmoon. Extract the
  weights with [`weights()`](https://rdrr.io/r/stats/weights.html) and
  pass them to
  [`halfmoon::check_balance()`](https://r-causal.github.io/halfmoon/reference/check_balance.html)
  and
  [`halfmoon::plot_balance()`](https://r-causal.github.io/halfmoon/reference/plot_balance.html).
- [positively](https://r-causal.github.io/positively/) diagnoses
  positivity violations and extrapolation across binary, categorical,
  and continuous exposures. balancing redistributes influence across
  units to achieve balance, and positivity problems are exactly where
  that redistribution becomes extrapolation, showing up as extreme
  weights and a collapsing effective sample size, so checking positivity
  with positively complements weighting with balancing.

## Learn more

- [`vignette("balancing")`](https://r-causal.github.io/balancing/articles/balancing.md)
  walks through a binary exposure workflow from weights to effects.
- [`vignette("choosing-a-method")`](https://r-causal.github.io/balancing/articles/choosing-a-method.md)
  compares the six methods and the constraint options.
- [`vignette("inference")`](https://r-causal.github.io/balancing/articles/inference.md)
  covers standard errors after balancing.
- [Causal Inference in R](https://www.r-causal.org/) is a book on causal
  inference methods in R.
