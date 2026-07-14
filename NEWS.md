# balancing 0.0.0.9000

* Added `balance()`, the entry point for fitting optimization-based balancing
  weights to a data frame with tidyselect covariate selection.

* Added six balancing methods, each a method specification passed to `balance()`:
  entropy balancing (`bw_entropy()`), inverse probability tilting (`bw_ipt()`),
  the covariate balancing propensity score (`bw_cbps()`), energy balancing
  (`bw_energy()`), characteristic function distance balancing (`bw_cfd()`), and
  stable balancing weights (`bw_sbw()`). The methods support binary,
  categorical, and continuous exposures across the ATE, ATT, ATU, and (for the
  covariate balancing propensity score) ATO estimands.

* Added `balance_terms()` to specify the covariate functions a method balances,
  covering moments, pairwise interactions, quantile indicators, and per-covariate
  balance tolerances.

* Added a diagnostics surface for fitted results: `print()`, `summary()`,
  `tidy()`, `glance()`, `ess()`, and `autoplot()` and `plot()` for love plots,
  weight distributions, and dual-variable charts.

* Added an `ipw()` method so a balancing fit drives the same effect-estimation
  workflow as a propensity score model. For the estimating-equation methods with
  a binary exposure, standard errors come from a stacked M-estimator that
  accounts for having estimated the weights, exposed through
  `estimating_equations()`.

* Added the `bw` weight vector class, a sibling of `propensity::psw` under the
  shared `causal_wts` parent.

* The numerical solvers are implemented in a Rust core.

* Added a `NEWS.md` file to track changes to the package.
