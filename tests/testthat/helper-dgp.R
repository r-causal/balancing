# Seeded, confounded data-generating processes for the property tests. Each
# returns a data frame whose covariates predict the exposure, so the unweighted
# imbalance between exposure groups is material and a balancing method has real
# work to do. `withr::local_seed()` keeps every draw reproducible without
# disturbing the session RNG.

# Binary exposure confounded by two continuous covariates and one factor.
sim_binary <- function(n = 500, seed = 2024) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  x3 <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  linear_predictor <- 0.8 *
    x1 -
    0.6 * x2 +
    0.4 * (x3 == "b") -
    0.5 * (x3 == "c")
  exposure <- stats::rbinom(n, 1L, stats::plogis(linear_predictor))
  data.frame(exposure = exposure, x1 = x1, x2 = x2, x3 = x3)
}

# Three-level categorical exposure confounded by two continuous covariates,
# drawn from a multinomial logit with reference level "a".
sim_categorical <- function(n = 500, seed = 2024) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  eta_b <- 0.7 * x1 - 0.3 * x2
  eta_c <- -0.5 * x1 + 0.6 * x2
  denom <- 1 + exp(eta_b) + exp(eta_c)
  prob_a <- 1 / denom
  prob_b <- exp(eta_b) / denom
  draw <- stats::runif(n)
  level <- ifelse(
    draw < prob_a,
    "a",
    ifelse(draw < prob_a + prob_b, "b", "c")
  )
  exposure <- factor(level, levels = c("a", "b", "c"))
  data.frame(exposure = exposure, x1 = x1, x2 = x2)
}

# Continuous exposure driven by one continuous covariate and one binary
# indicator. The indicator crosses the boundary as a raw zero/one column, so a
# fit over both covariates has to hold the indicator's marginal proportion as
# well as decorrelate it from the exposure.
sim_continuous_indicator <- function(n = 400, seed = 101) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  g <- stats::rbinom(n, 1L, 0.4)
  exposure <- 0.5 * x1 + 0.3 * g + stats::rnorm(n)
  data.frame(exposure = exposure, x1 = x1, g = g)
}

# Continuous exposure correlated with two continuous covariates, so the
# unweighted exposure-covariate correlations are materially nonzero.
sim_continuous <- function(n = 500, seed = 2024) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  exposure <- 0.9 * x1 - 0.6 * x2 + stats::rnorm(n)
  data.frame(exposure = exposure, x1 = x1, x2 = x2)
}
