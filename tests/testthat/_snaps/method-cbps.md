# a cbps fit prints its summary block

    Code
      fit <- balance(data, exposure, c(x1, x2), method = cbps(), estimand = "ate")
      fit
    Output
      
      -- Covariate balancing propensity score ----------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Effective sample size (0: 180.1 and 1: 236.0)
      Solver: converged in 4 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0000 (standardized mean difference)

