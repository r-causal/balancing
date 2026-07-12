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

