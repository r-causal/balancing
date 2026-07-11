# balancing_convergence_warning: the iteration cap is reached

    Code
      expr
    Condition <balancing_convergence_warning>
      Warning in `balance()`:
      The solver did not reach its convergence tolerance.
      i Increase `max_iterations` or loosen `convergence_tolerance` in `balancing::entropy_balance()`.
    Output
      
      -- Entropy balancing -----------------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Effective sample size (0: 184.6 and 1: 247.6)
      Solver: did not converge in 2 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0000 (standardized mean difference)

# balancing_balance_warning: achieved balance exceeds the tolerance

    Code
      expr
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
      Effective sample size (overall: 371.0)
      Solver: converged in 30 iterations
      Constraints: 2 terms (tolerance 0.1)
      Largest imbalance: 0.1751 (correlation)

# balancing_class_downgrade_warning: mismatched estimands

    Code
      expr
    Condition <balancing_class_downgrade_warning>
      Warning in `warn_bw_downgrade()`:
      Cannot combine <bw> weights with <bw>.
      i The result is a plain <numeric> vector without the estimand.
      i Set a common estimand on both vectors to keep the <bw> class.
    Output
      [1] 1 2 3 4

# alert: the detected exposure type is announced

    Code
      invisible(balance(data, exposure, c(x1, x2), method = entropy_balance(),
      estimand = "ate"))
    Message
      i Treating `.exposure` as binary.

# alert: aliased constraint columns are dropped

    Code
      invisible(balance(data, exposure, c(x1, x1_copy, x2), method = entropy_balance(),
      estimand = "ate", exposure_type = "binary"))
    Message
      i Dropping aliased constraint "x1_copy".

# alert: moments above one on a binary covariate are ignored

    Code
      invisible(balance(data, exposure, c(x2, flag), method = entropy_balance(),
      estimand = "ate", exposure_type = "binary", constraints = balance_terms(
        moments = 2L)))
    Message
      i Ignoring moments above one for the binary covariate "flag".

