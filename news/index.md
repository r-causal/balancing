# Changelog

## balancing 0.0.0.9000

- Added
  [`balance()`](https://r-causal.github.io/balancing/reference/balance.md),
  which fits optimization-based balancing weights to a data frame with
  tidyselect covariate selection.

- Added six balancing methods, each a method specification passed to
  [`balance()`](https://r-causal.github.io/balancing/reference/balance.md):
  entropy balancing
  ([`bw_entropy()`](https://r-causal.github.io/balancing/reference/bw_entropy.md)),
  inverse probability tilting
  ([`bw_ipt()`](https://r-causal.github.io/balancing/reference/bw_ipt.md)),
  the covariate balancing propensity score
  ([`bw_cbps()`](https://r-causal.github.io/balancing/reference/bw_cbps.md)),
  energy balancing
  ([`bw_energy()`](https://r-causal.github.io/balancing/reference/bw_energy.md)),
  characteristic function distance balancing
  ([`bw_cfd()`](https://r-causal.github.io/balancing/reference/bw_cfd.md)),
  and stable balancing weights
  ([`bw_sbw()`](https://r-causal.github.io/balancing/reference/bw_sbw.md)).
  The methods support binary, categorical, and continuous exposures
  across the ATE, ATT, ATU, and (for the covariate balancing propensity
  score) ATO estimands.

- Added
  [`balance_terms()`](https://r-causal.github.io/balancing/reference/balance_terms.md)
  to specify the covariate functions a method balances, covering
  moments, pairwise interactions, quantile indicators, and per-covariate
  balance tolerances.

- [`bw_sbw()`](https://r-causal.github.io/balancing/reference/bw_sbw.md)
  supports three weight-dispersion norms: `"l2"` (the default, sum of
  squared weights), `"l1"` (sum of absolute deviations from uniform),
  and `"linf"` (the maximum deviation from uniform).

- Added a summary surface for fitted results:
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) report the method,
  estimand, solver status, constraint count, and largest imbalance;
  [`weights()`](https://rdrr.io/r/stats/weights.html) returns the `bw`
  weight vector; and the fit carries its balance table and solver dual
  variables for inspection. Broader balance assessment, including
  covariates the weights did not target, lives in the halfmoon package.

- Added an
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  method so a balancing fit drives the same effect-estimation workflow
  as a propensity score model. For the estimating-equation methods with
  a binary or categorical exposure, standard errors come from a stacked
  M-estimator that accounts for having estimated the weights, exposed
  through
  [`estimating_equations()`](https://r-causal.github.io/balancing/reference/estimating_equations.md).
  The outcome model may adjust for covariates alongside the exposure, in
  which case the marginal means are standardized over the estimand’s
  target population. A categorical exposure reports each non-reference
  level against the reference level, and the estimates table names each
  contrast in a `comparison` column.

- [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  also accepts a continuous exposure from an entropy balancing fit with
  exact balance, whose estimand is the average treatment effect. It
  reports the exposure coefficient of a weighted marginal structural
  model, which must carry exactly one term in the exposure, as a single
  row named for the outcome model’s link: `slope` for an identity link,
  `log(or)` for a logit, and `log(rr)` for a log link. The standard
  error comes from the same stacked M-estimator, the weight parameters
  above the outcome-model score, and is a large-sample one: at a few
  hundred observations it runs anticonservative.

- Added the `bw` weight vector class, a sibling of
  [`propensity::psw`](https://r-causal.github.io/propensity/reference/psw.html)
  under the shared `causal_wts` parent.

- The numerical solvers are implemented in a Rust core.

- Added a `NEWS.md` file to track changes to the package.
