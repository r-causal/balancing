# a kernel balancing fit prints its summary block

    Code
      fit <- balance(data, exposure, c(x1, x2), method = bw_cfd(), estimand = "ate")
      fit
    Output
      
      -- Characteristic function distance balancing ----------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: converged in 1625 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0051 (standardized mean difference)

