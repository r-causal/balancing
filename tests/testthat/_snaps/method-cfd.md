# the ignored-tolerance warning records its class and printed tolerance

    Code
      balance(data, exposure, c(x1, x2), method = bw_cfd(), estimand = "ate",
      constraints = balance_terms(tolerance = 0.1))
    Condition <balancing_ignored_argument_warning>
      Warning in `method(fit_method, balancing::bw_cfd)`:
      `tolerance` relaxes added constraints, but this fit has none to relax.
      i Drop `tolerance` from `balance_terms()`, or add constraints with `moments` or `interactions`, or with `quantiles` for a discrete exposure.
    Output
      -- Characteristic function distance balancing ----------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 150
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0329 (standardized mean difference)

# a kernel balancing fit prints its summary block

    Code
      fit <- balance(data, exposure, c(x1, x2), method = bw_cfd(), estimand = "ate")
      fit
    Output
      -- Characteristic function distance balancing ----------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.00515 (standardized mean difference)

