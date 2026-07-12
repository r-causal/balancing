# an ipt fit prints its summary block

    Code
      fit <- balance(data, exposure, c(x1, x2), method = ipt(), estimand = "ate")
      fit
    Output
      
      -- Inverse probability tilting -------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Effective sample size (0: 176.2 and 1: 240.4)
      Solver: converged in 4 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0000 (standardized mean difference)

