# the required-tolerance message names the tuning parameter

    Code
      balance(data, exposure, c(x1, x2), method = bw_sbw(), estimand = "ate")
    Condition <balancing_constraints_error>
      Error in `method(fit_method, balancing::bw_sbw)`:
      ! Stable balancing weights require a positive balance tolerance.
      x No constraint carries a tolerance above zero.
      i Set `tolerance` in `balance_terms()` to a positive value, the central tuning parameter for `bw_sbw()`.

# an infeasible constraint set raises balancing_infeasible_error

    Code
      balance(data, exposure, c(x1, x2), method = bw_sbw(), estimand = "ate",
      constraints = balance_terms(tolerance = 0.01))
    Condition <balancing_infeasible_error>
      Error in `balance()`:
      ! The balancing problem is infeasible.
      x The solver reported that the constraints cannot be satisfied together.
      i Raise `tolerance` in `balance_terms()`, lower the moments, or drop interactions.

# the backend-fallback alert announces the switch

    Code
      alert_backend_fallback()
    Message <balancing_fallback_message>
      i The default solver certified this problem infeasible; re-solved with the "clarabel" backend.

# a stable balancing fit prints its summary block

    Code
      fit <- balance(data, exposure, c(x1, x2), method = bw_sbw(), estimand = "ate",
      constraints = balance_terms(tolerance = 0.05))
      fit
    Output
      
      -- Stable balancing weights ----------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: converged in 100 iterations
      Constraints: 2 terms (tolerance 0.05)
      Largest imbalance: 0.1000 (standardized mean difference)

