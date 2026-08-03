# an bw_ipt fit prints its summary block

    Code
      fit <- balance(data, exposure, c(x1, x2), method = bw_ipt(), estimand = "ate")
      fit
    Output
      -- Inverse probability tilting -------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0000 (standardized mean difference)

