# print() of a binary ate fit is stable

    Code
      print(fit)
    Output
      -- Entropy balancing -----------------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: <1e-7 (standardized mean difference)

# summary() of a binary ate fit is stable

    Code
      summary(fit)
    Output
      -- Entropy balancing -----------------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: <1e-7 (standardized mean difference)
      
      -- Weights --
      
      Range: 0.220 to 3.122
      Mean: 1.000
      Coefficient of variation: 0.398
      
      -- Balance --
      
        term   kind statistic group unweighted     weighted tolerance
      1   x1 moment       smd     1  0.552 <1e-7         0
      2   x2 moment       smd     1  0.416 <1e-7         0
        within_tolerance
      1             TRUE
      2             TRUE

# print() of a binary att fit renders the focal level

    Code
      print(fit)
    Output
      -- Entropy balancing -----------------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "att" (focal level "1")
      Observations: 500
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: <1e-7 (standardized mean difference)

# print() of a continuous ate fit is stable

    Code
      print(fit)
    Output
      -- Entropy balancing -----------------------------------------------------------
      Exposure: "exposure" (continuous)
      Estimand: "ate"
      Observations: 500
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: <1e-7 (correlation)

# print() of a categorical ate fit lists every level

    Code
      print(fit)
    Output
      -- Entropy balancing -----------------------------------------------------------
      Exposure: "exposure" (categorical)
      Estimand: "ate"
      Observations: 500
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: <1e-7 (standardized mean difference)

# print() of a fit whose factor lost a level is stable

    Code
      print(fit)
    Output
      -- Entropy balancing -----------------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: converged in <n> iterations
      Constraints: 4 terms (tolerance 0)
      Largest imbalance: <1e-7 (standardized mean difference)

# print() of an energy fit is stable

    Code
      print(fit)
    Output
      -- Energy balancing ------------------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0118 (standardized mean difference)

# summary() of an energy fit reports the weight floor count

    Code
      summary(fit)
    Output
      -- Energy balancing ------------------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0118 (standardized mean difference)
      
      -- Weights --
      
      Range: 0.027 to 7.351
      Mean: 1.000
      Coefficient of variation: 0.769
      Weights at the minimum-weight floor: <n> of 500
      
      -- Balance --
      
        term   kind statistic group unweighted   weighted tolerance within_tolerance
      1   x1 moment       smd     1  0.552 0.0118         0            FALSE
      2   x2 moment       smd     1  0.416 0.00511         0            FALSE

# summary() of a stable balancing fit reports the weight floor count

    Code
      summary(fit)
    Output
      -- Stable balancing weights ----------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0.05)
      Largest imbalance: 0.1 (standardized mean difference)
      
      -- Weights --
      
      Range: 0.000 to 1.957
      Mean: 1.000
      Coefficient of variation: 0.305
      Weights at the minimum-weight floor: <n> of 500
      
      -- Balance --
      
        term   kind statistic group unweighted weighted tolerance within_tolerance
      1   x1 moment       smd     1  0.552      0.1      0.05             TRUE
      2   x2 moment       smd     1  0.416      0.1      0.05             TRUE

# summary() of a cfd fit reports the weight floor count

    Code
      summary(fit)
    Output
      -- Characteristic function distance balancing ----------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Solver: converged in <n> iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.00515 (standardized mean difference)
      
      -- Weights --
      
      Range: 0.000 to 20.5
      Mean: 1.000
      Coefficient of variation: 1.732
      Weights at the minimum-weight floor: <n> of 500
      
      -- Balance --
      
        term   kind statistic group unweighted   weighted tolerance within_tolerance
      1   x1 moment       smd     1  0.552 0.00515         0            FALSE
      2   x2 moment       smd     1  0.416 0.00222         0            FALSE

