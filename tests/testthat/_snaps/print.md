# print() of a binary ate fit is stable

    Code
      print(fit)
    Output
      
      -- Entropy balancing -----------------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Effective sample size (0: 184.6 and 1: 247.6)
      Solver: converged in 6 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0000 (standardized mean difference)

# summary() of a binary ate fit is stable

    Code
      summary(fit)
    Output
      
      -- Entropy balancing -----------------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Effective sample size (0: 184.6 and 1: 247.6)
      Solver: converged in 6 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0000 (standardized mean difference)
      
      -- Weights --
      
      Range: 0.220 to 3.122
      Mean: 1.000
      Coefficient of variation: 0.398
      
      -- Balance --
      
        term   kind statistic group unweighted     weighted tolerance
      1   x1 moment       smd     1  0.5515682 6.438292e-11         0
      2   x2 moment       smd     1  0.4160432 1.014787e-10         0
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
      Effective sample size (0: 119.5 and 1: 277.0)
      Solver: converged in 5 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0000 (standardized mean difference)

# print() of a continuous ate fit is stable

    Code
      print(fit)
    Output
      
      -- Entropy balancing -----------------------------------------------------------
      Exposure: "exposure" (continuous)
      Estimand: "ate"
      Observations: 500
      Effective sample size (overall: 272.6)
      Solver: converged in 6 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0000 (correlation)

# print() of a categorical ate fit lists every level

    Code
      print(fit)
    Output
      
      -- Entropy balancing -----------------------------------------------------------
      Exposure: "exposure" (categorical)
      Estimand: "ate"
      Observations: 500
      Effective sample size (a: 139.0, b: 131.8, and c: 105.2)
      Solver: converged in 4 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0000 (standardized mean difference)

# print() of an energy fit is stable

    Code
      print(fit)
    Output
      
      -- Energy balancing ------------------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Effective sample size (0: 144.9 and 1: 169.9)
      Solver: converged in 75 iterations
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
      Effective sample size (0: 144.9 and 1: 169.9)
      Solver: converged in 75 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0118 (standardized mean difference)
      
      -- Weights --
      
      Range: 0.027 to 7.351
      Mean: 1.000
      Coefficient of variation: 0.769
      Weights at the minimum-weight floor: 0 of 500
      
      -- Balance --
      
        term   kind statistic group unweighted   weighted tolerance within_tolerance
      1   x1 moment       smd     1  0.5515682 0.01177578         0            FALSE
      2   x2 moment       smd     1  0.4160432 0.00511053         0            FALSE

# summary() of a stable balancing fit reports the weight floor count

    Code
      summary(fit)
    Output
      -- Stable balancing weights ----------------------------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Effective sample size (0: 198.3 and 1: 259.6)
      Solver: converged in 100 iterations
      Constraints: 2 terms (tolerance 0.05)
      Largest imbalance: 0.1000 (standardized mean difference)
      
      -- Weights --
      
      Range: 0.000 to 1.957
      Mean: 1.000
      Coefficient of variation: 0.305
      Weights at the minimum-weight floor: 1 of 500
      
      -- Balance --
      
        term   kind statistic group unweighted weighted tolerance within_tolerance
      1   x1 moment       smd     1  0.5515682      0.1      0.05             TRUE
      2   x2 moment       smd     1  0.4160432      0.1      0.05             TRUE

# summary() of a cfd fit reports the weight floor count

    Code
      summary(fit)
    Output
      -- Characteristic function distance balancing ----------------------------------
      Exposure: "exposure" (binary)
      Estimand: "ate"
      Observations: 500
      Effective sample size (0: 60.4 and 1: 65.3)
      Solver: converged in 1625 iterations
      Constraints: 2 terms (tolerance 0)
      Largest imbalance: 0.0051 (standardized mean difference)
      
      -- Weights --
      
      Range: 0.000 to 20.514
      Mean: 1.000
      Coefficient of variation: 1.732
      Weights at the minimum-weight floor: 13 of 500
      
      -- Balance --
      
        term   kind statistic group unweighted    weighted tolerance within_tolerance
      1   x1 moment       smd     1  0.5515682 0.005149743         0            FALSE
      2   x2 moment       smd     1  0.4160432 0.002219728         0            FALSE

