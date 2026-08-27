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

- [`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
  resolves its exposure type through causalgenerics, which owns the
  detection heuristics, the announcement, and the refusal of a declared
  type across the r-causal packages. Four things change with the move.
  Missing values no longer count as a level of their own, so a
  two-valued exposure that also carries missingness reads as binary
  rather than continuous;
  [`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
  refuses a missing exposure before it reads a type, so no fit reaches
  the new reading. A factor or character exposure taking a single level
  reads as binary rather than categorical, and is refused for taking one
  level either way. The announcement of a detected type loses its
  trailing period. An `exposure_type` the data cannot carry is refused
  with the `causalgenerics_forced_exposure_type` condition rather than
  `balancing_exposure_type_error`, which stays the class of balancing’s
  own refusal of an exposure type the method cannot fit.

- The constraint expansion now tests for aliasing against the intercept
  every balancing method carries, rather than against the constraint
  columns alone. The level indicators of a factor sum to that constant,
  so one indicator per factor was redundant in the geometry the solver
  sees while surviving the old check; the entropy estimating equations
  were rank deficient by construction whenever a factor was balanced,
  and the flat direction that left behind produced order-dependent
  solver failures. The redundant indicator, which is the last level of
  each factor, is now dropped with an informational alert, and the
  balance table reports the surviving levels. Balance on the levels that
  remain implies balance on the omitted one, so no fit loses a
  constraint it previously met.

- A failed Newton solve of the exact entropy problem is now retried once
  with the L-BFGS-then-Newton hybrid. Newton starts from the base
  measure and is the only solver that reaches machine precision on the
  estimating equations, but a flat or badly scaled constraint set can
  leave that cold start short of its tolerance or send it out to weights
  that are not finite; the hybrid reaches a neighborhood with L-BFGS
  first and polishes it with Newton, so it clears such problems at the
  same precision. A successful retry announces itself and
  `@solver_status` records the solver the fit came from. Pinning the
  `balancing.entropy_solver` option disables the retry, and a fit that
  exhausts both solvers names them in its warning or error.

- When the variance engine cannot invert the stacked bread,
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  now reads the fit’s own estimating equations before it reports. A
  rank-deficient fit block is what makes the whole stack singular, so
  the refusal names the rank it found and points at the constraint
  columns to go and look at, in place of the reading that leaves the
  caller choosing between estimating functions that are not finite and a
  bread that is singular. A full-rank fit block reports as it did
  before, since the fit is then not what went wrong.

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
  contrast in a `contrast` column.

- [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  also accepts a continuous exposure from an entropy balancing fit with
  exact balance, whose estimand is the average treatment effect. It
  reports the dose response of a weighted marginal structural model. An
  exposure entering through one design column is the whole of that
  response, so its coefficient is the response’s slope and the table
  holds a single row named for the outcome model’s link: `slope` for an
  identity link, `log(or)` for a logit, and `log(rr)` for a log link.
  The standard error comes from the same stacked M-estimator, the weight
  parameters above the outcome-model score, and is a large-sample one:
  at a few hundred observations it runs anticonservative.

- A continuous marginal structural model may spread the exposure over
  several design columns, as in `y ~ exposure + I(exposure^2)` or
  `y ~ splines::ns(exposure, 3)`. What it has to keep is variable
  membership: every term reading the exposure reads the exposure alone,
  and a term reading a covariate alongside it is refused, since its
  coefficient is a change in the dose response per unit of that
  covariate rather than an effect a row could name. Such a model reports
  one row per exposure coefficient and gains a `contrast` column naming
  each row after the coefficient the fit names. The scale word steps
  back to `coef` at an identity link, since a curve has a different
  slope at every dose; a logit still reports `log(or)` and a log link
  `log(rr)`, because a coefficient of those models is a log ratio
  whatever column it multiplies.

- [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  takes a `.by` argument naming a modifier to report the effects within
  the levels of. A result carrying one reports the whole-sample effects,
  then those same effects within each subgroup, then each non-reference
  subgroup against the reference one, and the estimates table gains a
  `group` column naming the subgroup each row was estimated in. A
  subgroup row carries the collapsible measures alone, since the odds
  ratio over a sample is not an average of the odds ratios within its
  subgroups. Every added row is a parameter of the same stacked system,
  so the subgroups covary and a contrast of two of them is estimated
  rather than differenced after the fact. An outcome model with no term
  reading both the exposure and the modifier is reported with a warning
  rather than refused, and a continuous exposure, which reports
  coefficients rather than contrasts of standardized means, has no
  subgroup effect for the argument to name and refuses it.

- An exposure declared as a crossing of two treatments by
  causalgenerics’ `joint_exposure()` is reported in those two treatments
  rather than cell against cell: the counterfactual mean of each cell,
  the simple effects of each treatment within a fixed level of the
  other, and the interaction between the two treatments, reported once
  because it is symmetric in them. The cell-against-cell rows are
  replaced rather than supplemented, so their labels appear nowhere a
  row is named. Every row is a parameter of the same stacked system,
  which makes each interaction a double difference of cell means by
  construction. A crossing whose two treatments carry one name is
  refused, since every row is keyed by the treatment it contrasts and
  the level the other is held at, and a declared crossing is refused
  with `.by` or with an estimand other than the ATE.

- Missing covariate data is handled by imputing first and analyzing
  within each completed dataset. Balance, weight, and estimate once per
  imputation, then pool the results with
  [`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html),
  which balancing re-exports from causalgenerics. The complete-data
  degrees of freedom come from the outcome models automatically, since a
  balancing result reports none of its own.

- An
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  result reads through the accessors causalgenerics registers on the
  class: [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html), and
  [`weights()`](https://rdrr.io/r/stats/weights.html).
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the covariance
  of the reported effects, labeled the way the result labels its rows,
  so the covariances between effect measures that the estimates table
  cannot report are available to a caller combining them. The stored
  outcome model carries the outcome-model block of the same covariance,
  so [`vcov()`](https://rdrr.io/r/stats/vcov.html) on it accounts for
  having estimated the weights.

- [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  takes an `effects` argument recording the reading the result it builds
  presents: `"marginal"`, the default, for the population-averaged
  causal contrasts, or `"conditional"` for the outcome model’s
  coefficient surface. The stacked system is solved whichever value is
  named, so the argument settles which reading the result presents and
  nothing else. In the conditional reading,
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html), and
  [`confint()`](https://rdrr.io/r/stats/confint.html) report the outcome
  model’s coefficients against the outcome block of the stacked
  sandwich, so their standard errors account for having estimated the
  weights rather than treating the weights as fixed.
  [`as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
  and
  [`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html),
  which move a result between the two readings afterwards, are
  re-exported here, so a result is flipped without a second attachment.
  A printed result names its reading beside the estimand and again over
  the table it decides.

- A pooled result moves between the two readings as an unpooled one
  does.
  [`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html)
  pools both readings of a set of results, so
  [`as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
  and
  [`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
  flip the pooled result afterwards, and its
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html), and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) take an
  `effects` argument for a single call. The conditional reading is there
  to pool because
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  hands every outcome model over already wrapped with its block of the
  stacked covariance; the pooling itself comes from causalgenerics.

- The balancing fit an
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  result stores carries the leading block of the stacked covariance, the
  one belonging to the weight parameters alone, so
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) on that fit reports the
  covariance of the weights under their `theta_w` names. A fit that has
  not been through such a result carries no covariance and raises
  `balancing_vcov_error`, pointing at the
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  route when its weights solve estimating equations and at the bootstrap
  workflow of the inference vignette when they do not.

- Added the `bw` weight vector class, a sibling of
  [`propensity::psw`](https://r-causal.github.io/propensity/reference/psw.html)
  under the shared `causal_wts` parent.

- The numerical solvers are implemented in a Rust core.

- Added a `NEWS.md` file to track changes to the package.
