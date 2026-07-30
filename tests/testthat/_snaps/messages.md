# balancing_convergence_warning: the iteration cap is reached

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(max_iterations = 3L),
      estimand = "ate")
    Condition <balancing_convergence_warning>
      Warning in `balance()`:
      The solver did not reach its convergence tolerance.
      i Increase `max_iterations` or loosen `convergence_tolerance` in `balancing::bw_entropy()`.
    Output
      
      -- Entropy balancing -----------------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: did not converge in 2 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0000 (standardized mean difference)

# balancing_balance_warning: achieved balance exceeds the tolerance

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), estimand = "ate",
      constraints = balance_terms(tolerance = 0.1))
    Condition <balancing_balance_warning>
      Warning in `balance()`:
      The achieved balance exceeds the requested tolerance.
      x The largest imbalance is 0.1751.
      i Raise `tolerance` in `balance_terms()`, lower the moments, or drop interactions.
    Output
      
      -- Entropy balancing -----------------------------------------------------------
      Exposure: "exposure" (continuous)
      Estimand: "ate"
      Observations: 500
      Solver: converged in 30 iterations
      Constraints: 2 terms (tolerance 0.1)
      Largest imbalance: 0.1751 (correlation)

# balancing_ignored_argument_warning: two_step without over_identified

    Code
      balance(data, exposure, c(x1, x2), method = bw_cbps(two_step = FALSE,
        over_identified = FALSE), estimand = "ate")
    Condition <balancing_ignored_argument_warning>
      Warning in `method(fit_method, balancing::bw_cbps)`:
      `two_step` applies only to the over-identified fit and is ignored.
      i The two-step weighting matrix belongs to the over-identified criterion, which `bw_cbps()` fits for a binary exposure with `over_identified = TRUE`.
    Output
      
      -- Covariate balancing propensity score ----------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: converged in 4 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0000 (standardized mean difference)

# balancing_ignored_argument_warning: over_identified for a categorical exposure

    Code
      balance(data, exposure, c(x1, x2), method = bw_cbps(over_identified = TRUE),
      estimand = "ate")
    Condition <balancing_ignored_argument_warning>
      Warning in `method(fit_method, balancing::bw_cbps)`:
      `over_identified` applies only to a binary exposure and is ignored.
      i A categorical exposure has no over-identified criterion, so the fit balances its moment conditions exactly.
    Output
      
      -- Covariate balancing propensity score ----------------------------------------
      Exposure: "exposure" (categorical)
      Estimand: "ate"
      Observations: 500
      Solver: converged in 4 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0000 (standardized mean difference)

# balancing_ignored_argument_warning: over_identified for a continuous exposure

    Code
      balance(data, exposure, c(x1, x2), method = bw_cbps(over_identified = TRUE),
      estimand = "ate")
    Condition <balancing_ignored_argument_warning>
      Warning in `method(fit_method, balancing::bw_cbps)`:
      `over_identified` applies only to a binary exposure and is ignored.
      i A continuous exposure has no over-identified criterion, so the fit balances its moment conditions exactly.
    Output
      
      -- Covariate balancing propensity score ----------------------------------------
      Exposure: "exposure" (continuous)
      Estimand: "ate"
      Observations: 500
      Solver: converged in 6 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0000 (correlation)

# balancing_ignored_argument_warning: every argument a categorical fit ignores

    Code
      balance(data, exposure, c(x1, x2), method = bw_cbps(over_identified = TRUE,
        two_step = FALSE), estimand = "ate")
    Condition <balancing_ignored_argument_warning>
      Warning in `method(fit_method, balancing::bw_cbps)`:
      `over_identified` applies only to a binary exposure and is ignored.
      i A categorical exposure has no over-identified criterion, so the fit balances its moment conditions exactly.
      Warning in `method(fit_method, balancing::bw_cbps)`:
      `two_step` applies only to the over-identified fit and is ignored.
      i The two-step weighting matrix belongs to the over-identified criterion, which `bw_cbps()` fits for a binary exposure with `over_identified = TRUE`.
    Output
      
      -- Covariate balancing propensity score ----------------------------------------
      Exposure: "exposure" (categorical)
      Estimand: "ate"
      Observations: 500
      Solver: converged in 4 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0000 (standardized mean difference)

# balancing_ignored_argument_warning: a clarabel pin the energy kernel cannot honor

    Code
      balance(data, exposure, c(x1, x2), method = bw_cfd(kernel = "energy"),
      estimand = "ate")
    Condition <balancing_ignored_argument_warning>
      Warning in `method(fit_method, balancing::bw_cfd)`:
      The `balancing.qp_backend` option is "clarabel", which the "energy" kernel cannot use, and is ignored.
      x The energy kernel's quadratic form is indefinite, and the interior-point backend solves only positive-semidefinite forms.
      i The fit used "osqp" instead, which `@solver_status` records.
    Output
      
      -- Characteristic function distance balancing ----------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: converged in 75 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0118 (standardized mean difference)

# balancing_class_downgrade_warning: mismatched estimands

    Code
      vctrs::vec_c(x, y)
    Condition <balancing_class_downgrade_warning>
      Warning in `warn_bw_downgrade()`:
      Cannot combine <bw> weights with <bw>.
      i The result is a plain <numeric> vector without the estimand.
      i Set a common estimand on both vectors to keep the <bw> class.
    Output
      [1] 1 2 3 4

# alert: the detected exposure type is announced

    Code
      invisible(balance(data, exposure, c(x1, x2), method = bw_entropy(), estimand = "ate"))
    Message
      i Treating `.exposure` as binary.

# alert: the exposure is excluded from a covariate selection

    Code
      invisible(balance(data, exposure, everything(), method = bw_entropy(),
      estimand = "ate"))
    Message
      i Dropping the exposure "exposure" from `.covariates`.
      i Treating `.exposure` as binary.

# alert: aliased constraint columns are dropped

    Code
      invisible(balance(data, exposure, c(x1, x1_copy, x2), method = bw_entropy(),
      estimand = "ate", exposure_type = "binary"))
    Message
      i Dropping aliased constraint "x1_copy".

# alert: moments above one on a binary covariate are ignored

    Code
      invisible(balance(data, exposure, c(x2, flag), method = bw_entropy(), estimand = "ate",
      exposure_type = "binary", constraints = balance_terms(moments = 2L)))
    Message
      i Ignoring moments above one for the binary covariate "flag".

