# Package index

## Estimate

Estimate balancing weights from a data frame.

- [`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
  : Estimate balancing weights

## Method specifications

The six balancing methods and the classes they share.

- [`bw_entropy()`](https://r-causal.github.io/balancing/reference/bw_entropy.md)
  : Entropy balancing
- [`bw_ipt()`](https://r-causal.github.io/balancing/reference/bw_ipt.md)
  : Inverse probability tilting
- [`bw_cbps()`](https://r-causal.github.io/balancing/reference/bw_cbps.md)
  : Covariate balancing propensity score
- [`bw_energy()`](https://r-causal.github.io/balancing/reference/bw_energy.md)
  : Energy balancing
- [`bw_cfd()`](https://r-causal.github.io/balancing/reference/bw_cfd.md)
  : Characteristic function distance balancing
- [`bw_sbw()`](https://r-causal.github.io/balancing/reference/bw_sbw.md)
  : Stable balancing weights
- [`balance_method()`](https://r-causal.github.io/balancing/reference/balance_method.md)
  [`estimating_equation_method()`](https://r-causal.github.io/balancing/reference/balance_method.md)
  [`quadratic_program_method()`](https://r-causal.github.io/balancing/reference/balance_method.md)
  : The balancing method specification classes
- [`supported_exposure_types()`](https://r-causal.github.io/balancing/reference/method_capabilities.md)
  [`supported_estimands()`](https://r-causal.github.io/balancing/reference/method_capabilities.md)
  [`supports_estimating_equations()`](https://r-causal.github.io/balancing/reference/method_capabilities.md)
  [`method_label()`](https://r-causal.github.io/balancing/reference/method_capabilities.md)
  [`fit_method()`](https://r-causal.github.io/balancing/reference/method_capabilities.md)
  : Method capability generics

## Balance constraints

Describe the covariate balance a method should achieve.

- [`balance_terms()`](https://r-causal.github.io/balancing/reference/balance_terms.md)
  : Balance constraints
- [`build_constraint_matrix()`](https://r-causal.github.io/balancing/reference/build_constraint_matrix.md)
  : Build the constraint matrix and recipe
- [`rebuild_constraint_matrix()`](https://r-causal.github.io/balancing/reference/rebuild_constraint_matrix.md)
  : Rebuild the constraint matrix from a recipe

## Inspect results

Summarize a fitted balancing result.

- [`balancing`](https://r-causal.github.io/balancing/reference/balancing.md)
  : A fitted balancing result
- [`weights,balancing::balancing-method`](https://r-causal.github.io/balancing/reference/weights-balancing-balancing-method.md)
  : Extract balancing weights
- [`new_bw()`](https://r-causal.github.io/balancing/reference/bw.md)
  [`bw()`](https://r-causal.github.io/balancing/reference/bw.md)
  [`as_bw()`](https://r-causal.github.io/balancing/reference/bw.md)
  [`is_bw()`](https://r-causal.github.io/balancing/reference/bw.md) :
  Balancing weight vectors
- [`balancing_options`](https://r-causal.github.io/balancing/reference/balancing_options.md)
  : Package options

## Effect estimation

Estimate causal effects and access the pieces inference needs.

- [`ipw.balancing`](https://r-causal.github.io/balancing/reference/ipw.balancing.md)
  : Inverse probability weighting for a balancing fit
- [`estimating_equations()`](https://r-causal.github.io/balancing/reference/estimating_equations.md)
  : Extract the estimating-equations container
- [`balancing_estimating_equations()`](https://r-causal.github.io/balancing/reference/balancing_estimating_equations.md)
  : The estimating-equations container

## Re-exports

- [`reexports`](https://r-causal.github.io/balancing/reference/reexports.md)
  [`ipw`](https://r-causal.github.io/balancing/reference/reexports.md)
  [`as_marginal`](https://r-causal.github.io/balancing/reference/reexports.md)
  [`as_conditional`](https://r-causal.github.io/balancing/reference/reexports.md)
  [`is_causal_wt`](https://r-causal.github.io/balancing/reference/reexports.md)
  [`estimand`](https://r-causal.github.io/balancing/reference/reexports.md)
  : Objects exported from other packages
