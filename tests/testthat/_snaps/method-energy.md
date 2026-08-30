# the ignored-tolerance warning records its class and message

    Code
      balance(data, exposure, c(x1, x2), method = bw_energy(), estimand = "ate",
      constraints = balance_terms(tolerance = 0.1))
    Condition <balancing_ignored_argument_warning>
      Warning in `method(fit_method, balancing::bw_energy)`:
      `tolerance` relaxes added constraints, but this fit has none to relax.
      i Drop `tolerance` from `balance_terms()`, or add constraints with `moments` or `interactions`, or with `quantiles` for a discrete exposure.
    Output
      -- Energy balancing ------------------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 150
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0437 (standardized mean difference)

# a continuous tolerance warns and is ignored

    Code
      balance(data, exposure, c(x1, x2), method = bw_energy(), estimand = "ate",
      constraints = balance_terms(tolerance = 0.1))
    Condition <balancing_ignored_argument_warning>
      Warning in `method(fit_method, balancing::bw_energy)`:
      `tolerance` relaxes added constraints, but this fit has none to relax.
      i Drop `tolerance` from `balance_terms()`, or add constraints with `moments` or `interactions`, or with `quantiles` for a discrete exposure.
    Output
      -- Energy balancing ------------------------------------------------------------
      Exposure: "exposure" (continuous)
      Estimand: "ate"
      Observations: 500
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.188 (correlation)

# an infeasible constraint set raises balancing_infeasible_error

    Code
      balance(data, exposure, c(x1, x2), method = bw_energy(), estimand = "ate",
      constraints = balance_terms(moments = 1L))
    Condition <balancing_infeasible_error>
      Error in `balance()`:
      ! The balancing problem is infeasible.
      x The solver reported that the constraints cannot be satisfied together.
      i Raise `tolerance` in `balance_terms()`, lower the moments, or drop interactions.

# an energy fit prints its summary block

    Code
      fit <- balance(data, exposure, c(x1, x2), method = bw_energy(), estimand = "ate")
      fit
    Output
      -- Energy balancing ------------------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0118 (standardized mean difference)

