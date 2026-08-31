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

# The two fixtures below are shaped for one route each rather than being general
# processes, and they live here because more than one file draws on them: a
# fixture reachable from only the file it was written in cannot be used by a
# spec that has to cross routes.

# A binary-exposure fixture whose effect differs across the levels of a
# two-level modifier. The modifier confounds the exposure as well as modifying
# its effect, so a fit that balances it has real work to do, and it rides along
# as a numeric indicator, `modifier_hi`, because that is the parameterization
# the weight parameters stay identified in.
#
# The modifier declares its levels in reverse alphabetical order on purpose. The
# reference subgroup every contrast of subgroups is measured against is the
# modifier's first level, which is `"lo"` here and would be `"hi"` for an
# implementation that sorted the levels itself.
ipw_by_fixture <- function(n = 400) {
  withr::with_seed(808, {
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    modifier <- factor(
      sample(c("lo", "hi"), n, replace = TRUE),
      levels = c("lo", "hi")
    )
    modifier_hi <- as.numeric(modifier == "hi")
    exposure <- stats::rbinom(
      n,
      1L,
      stats::plogis(0.7 * x1 - 0.5 * x2 + 0.6 * modifier_hi)
    )
    y <- stats::rbinom(
      n,
      1L,
      stats::plogis(
        -0.6 +
          0.2 * exposure +
          0.5 * x1 +
          0.3 * modifier_hi +
          1.4 * exposure * modifier_hi
      )
    )
    y_cont <- 1 +
      0.2 * exposure +
      0.5 * x1 -
      0.3 * x2 +
      1.2 * exposure * modifier_hi +
      stats::rnorm(n)
    data.frame(
      exposure = exposure,
      x1 = x1,
      x2 = x2,
      modifier = modifier,
      modifier_hi = modifier_hi,
      y = y,
      y_cont = y_cont
    )
  })
}

# Two binary treatments, the second depending on the first and both on a
# covariate, with an outcome carrying a real interaction between them. A binary
# and a gaussian outcome are drawn so each reported scale has something to read,
# and a modifier is drawn for the `.by` refusal.
ipw_joint_fixture <- function(n = 700) {
  withr::with_seed(4210, {
    x1 <- stats::rnorm(n)
    a <- stats::rbinom(n, 1L, stats::plogis(0.3 * x1))
    e <- stats::rbinom(n, 1L, stats::plogis(-0.2 + 0.5 * x1 - 0.4 * a))
    y <- stats::rbinom(
      n,
      1L,
      stats::plogis(-0.5 + 0.7 * a + 0.5 * e + 0.6 * x1 + 0.9 * a * e)
    )
    y_cont <- 1 +
      0.6 * a +
      0.4 * e +
      0.5 * x1 +
      0.8 * a * e +
      stats::rnorm(n)
    data <- data.frame(
      x1 = x1,
      y = y,
      y_cont = y_cont,
      a = factor(a, levels = c(0L, 1L)),
      e = factor(e, levels = c(0L, 1L)),
      modifier = factor(
        ifelse(x1 > 0, "hi", "lo"),
        levels = c("lo", "hi")
      )
    )
    # Assigned rather than built inside `data.frame()`, which would coerce the
    # crossing away before anything could read it.
    data$joint <- causalgenerics::joint_exposure(a = data$a, e = data$e)
    data
  })
}
