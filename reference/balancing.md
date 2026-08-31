# A fitted balancing result

`balancing` is the S7 object
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
returns. It carries the balancing weights, the resolved method and
estimand, the covariate expansion recipe, a balance table, and the
solver diagnostics. The raw data are not stored; the recipe carries what
is needed to rebuild the constraint matrix.

## Arguments

- weights:

  The balancing weights, a
  [bw](https://r-causal.github.io/balancing/reference/bw.md) vector.

- method:

  The fitted
  [balance_method](https://r-causal.github.io/balancing/reference/balance_method.md)
  specification.

- estimand:

  The resolved estimand string.

- exposure:

  The exposure column name.

- exposure_type:

  The resolved exposure type.

- exposure_levels:

  The exposure levels the fit weighted, as character strings in the
  order `levels(factor(exposure))` gives: a factor's declared order for
  a factor exposure, and the values sorted in their own type otherwise,
  so a numeric dose of 9 and 10 orders 9 first. Levels no observation
  takes are dropped and do not appear. The first element is the
  reference level every contrast in
  [`ipw()`](https://r-causal.github.io/balancing/reference/ipw.balancing.md)
  is measured against. A continuous exposure carries no levels, so the
  vector is empty.

- covariates:

  The covariate column names that contributed at least one retained
  constraint column, in selection order. The constraint expansion drops
  constant and aliased columns, and a covariate can be selected without
  contributing any column at all, so this records what the fit
  constrained rather than what was requested; the request itself stays
  in `call`. A fit whose balance comes from its objective rather than
  from constraints, such as energy or kernel balancing with no moment
  constraints, therefore records no covariates even though its objective
  reads every selected one.

- focal_level:

  The focal exposure level for `"att"` and `"atc"`, or `NULL`. This is
  the fitted object's property, set from the `.focal_level` argument of
  [`balance()`](https://r-causal.github.io/balancing/reference/balance.md).

- n:

  The number of observations.

- constraints:

  The resolved
  [balance_terms](https://r-causal.github.io/balancing/reference/balance_terms.md)
  specification, or `NULL`.

- recipe:

  The covariate expansion recipe, a list of per-column records.

- balance_table:

  The achieved balance, one row per constraint term.

- duals:

  Solver dual variables for diagnostics, or `NULL`.

- coefficients:

  The fitted coefficients or dual variables, or `NULL`.

- converged:

  Whether the solver met its convergence criterion.

- iterations:

  The solver iteration count. An energy fit that could not reach its
  tolerance re-solves at a reachable one, and when that re-solve
  converges this sums the original and the fallback solve, so it can
  exceed the requested `max_iterations`. When the re-solve does not
  converge the fit reports the original solve alone, so the count stays
  within the cap. A continuous energy or stable balancing fit with a
  positive tolerance refines the bound it hands the solver over several
  passes, each a solve of its own, and this sums every one of them. See
  [`bw_energy()`](https://r-causal.github.io/balancing/reference/bw_energy.md)
  and
  [`bw_sbw()`](https://r-causal.github.io/balancing/reference/bw_sbw.md)
  for the fuller account.

- objective:

  The solved objective value.

- solver_status:

  The solver that produced the result.

- estimating_equations:

  The
  [balancing_estimating_equations](https://r-causal.github.io/balancing/reference/balancing_estimating_equations.md)
  container, or `NULL`.

- vcov:

  The covariance of the fit's own weight parameters, a `p` by `p` matrix
  named for them, or `NULL`. A covariance for those parameters comes
  from the stacked system an
  [`ipw()`](https://r-causal.github.io/balancing/reference/ipw.balancing.md)
  result assembles, so
  [`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
  leaves this empty and the result fills it in on the copy of the fit it
  stores, where [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html)
  reads it back.

- sampling_weights:

  The sampling weights, or `NULL`.

- call:

  The originating call.

## Value

A `balancing` object.
