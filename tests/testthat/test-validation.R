# The sampling-weight validator guards type, length, finiteness, sign, and a
# degenerate total before a fit begins. Each branch raises a classed balancing
# error; the user-facing messages are snapshotted so a change in wording or class
# is caught. A valid vector passes through invisibly and unchanged.
#
# Finiteness is shared: `validate_finite()` is the single place the package
# decides what a missing value and an infinity mean for a numeric input, and the
# sampling-weight validator routes both checks through it.

test_that("validate_sampling_weights() returns valid weights invisibly", {
  weights <- c(0.5, 1, 1.5)
  expect_identical(validate_sampling_weights(weights, n = 3), weights)
})

test_that("validate_sampling_weights() rejects a non-numeric vector", {
  expect_balancing_error(
    validate_sampling_weights(c("a", "b"), n = 2)
  )
})

test_that("validate_sampling_weights() rejects a length mismatch", {
  expect_balancing_error(
    validate_sampling_weights(c(1, 2), n = 3)
  )
})

test_that("validate_sampling_weights() rejects missing values", {
  expect_balancing_error(
    validate_sampling_weights(c(1, NA, 1), n = 3)
  )
})

test_that("validate_sampling_weights() rejects negative values", {
  expect_balancing_error(
    validate_sampling_weights(c(1, -1, 1), n = 3)
  )
})

test_that("validate_sampling_weights() rejects infinite values of either sign", {
  expect_error(
    validate_sampling_weights(c(1, Inf, 1), n = 3),
    class = "balancing_range_error"
  )
  expect_error(
    validate_sampling_weights(c(1, -Inf, 1), n = 3),
    class = "balancing_range_error"
  )
  expect_balancing_error(
    validate_sampling_weights(c(1, Inf, 1), n = 3)
  )
})

test_that("validate_sampling_weights() rejects an all-zero vector", {
  expect_error(
    validate_sampling_weights(rep(0, 3), n = 3),
    class = "balancing_range_error"
  )
  expect_balancing_error(
    validate_sampling_weights(rep(0, 3), n = 3)
  )
})

# A unit given no sampling weight is pinned at zero rather than dropped, so an
# individual zero is a supported input and only a vector with no mass at all is
# degenerate.
test_that("validate_sampling_weights() accepts individual zero weights", {
  weights <- c(0, 1, 2)
  expect_identical(validate_sampling_weights(weights, n = 3), weights)
})

# ---- The base measure gate --------------------------------------------------

# The base measure is the product of the sampling weights and any base weights.
# Every group total a fit divides by is taken under it, so a group with no mass
# leaves that group's constraint targets and reported total undefined, and a
# measure with no mass at all leaves the pooled target undefined too.
test_that("validate_base_measure() passes a measure with mass in every group", {
  measure <- c(0, 1, 2, 1)
  groups <- list(a = c(1L, 2L), b = c(3L, 4L))
  expect_identical(validate_base_measure(measure, groups), measure)
})

test_that("validate_base_measure() names a group with no mass", {
  groups <- list(a = c(1L, 2L), b = c(3L, 4L))
  expect_error(
    validate_base_measure(c(0, 0, 2, 1), groups),
    class = "balancing_range_error"
  )
  expect_error(validate_base_measure(c(0, 0, 2, 1), groups), "\"a\"")
})

test_that("validate_base_measure() refuses a measure with no mass at all", {
  expect_error(
    validate_base_measure(rep(0, 4), NULL),
    class = "balancing_range_error"
  )
  expect_error(
    validate_base_measure(rep(0, 4), list(a = 1:2, b = 3:4)),
    class = "balancing_range_error"
  )
})

test_that("validate_base_measure() accepts a grouped fit with no groups", {
  measure <- c(0, 1, 2, 1)
  expect_identical(validate_base_measure(measure, NULL), measure)
})

# ---- The shared finiteness gate --------------------------------------------

test_that("validate_finite() returns a finite vector invisibly", {
  values <- c(0, 1, 2.5)
  expect_identical(validate_finite(values, "base_weights"), values)
})

test_that("validate_finite() rejects missing and NaN values", {
  expect_error(
    validate_finite(c(1, NA, 1), "base_weights"),
    class = "balancing_missing_error"
  )
  expect_error(
    validate_finite(c(1, NaN, 1), "base_weights"),
    class = "balancing_missing_error"
  )
})

test_that("validate_finite() rejects infinite values of either sign", {
  expect_error(
    validate_finite(c(1, Inf, 1), "base_weights"),
    class = "balancing_range_error"
  )
  expect_error(
    validate_finite(c(1, -Inf, 1), "base_weights"),
    class = "balancing_range_error"
  )
})

test_that("validate_finite() leaves a non-numeric vector alone", {
  # Character and factor columns reach the covariate gate too, and neither can
  # hold an infinity, so the check passes them through untouched.
  values <- c("a", "b")
  expect_identical(validate_finite(values, ".covariates"), values)
})

# ---- Solver boundary guards ------------------------------------------------

# The R-side validation above is the surface a user meets, and no fit can carry a
# non-finite value past it. These specs address the solver entry points directly,
# so the boundary refuses a non-finite input on its own rather than handing it to
# a solve whose arithmetic it would silently poison.

test_that("solve_entropy() refuses non-finite base and sampling weights", {
  covs <- c(-1, -1, 1, 1)
  group_idx <- c(0L, 0L, 1L, 1L)
  options <- list(threads = 1L)
  expect_error(
    solve_entropy(
      covs,
      group_idx,
      0,
      c(1, 1, 1, Inf),
      rep(1, 4),
      0,
      1,
      options
    ),
    "non-finite"
  )
  expect_error(
    solve_entropy(
      covs,
      group_idx,
      0,
      rep(1, 4),
      c(1, 1, 1, Inf),
      0,
      1,
      options
    ),
    "non-finite"
  )
})

test_that("solve_entropy_cont() refuses non-finite base and sampling weights", {
  covs <- c(-1, -0.5, 0.5, 1)
  options <- list(threads = 1L)
  expect_error(
    solve_entropy_cont(
      covs,
      0,
      0,
      1L,
      c(1, 1, 1, Inf),
      rep(1, 4),
      4,
      options
    ),
    "non-finite"
  )
  expect_error(
    solve_entropy_cont(
      covs,
      0,
      0,
      1L,
      rep(1, 4),
      c(1, 1, 1, Inf),
      4,
      options
    ),
    "non-finite"
  )
})

test_that("solve_cbps() refuses non-finite covariates and sampling weights", {
  covs <- c(-1, -0.5, 0.5, 1)
  treat <- c(0L, 0L, 1L, 1L)
  options <- list(threads = 1L)
  expect_error(
    solve_cbps(
      c(-1, -0.5, 0.5, Inf),
      covs,
      treat,
      rep(1, 4),
      "ate",
      "logit",
      FALSE,
      FALSE,
      options
    ),
    "`covs_mod` contains a non-finite value at position 4",
    fixed = TRUE
  )
  # The balance design is a second numeric block with a guard of its own, and the
  # just-identified fit passes the same matrix twice, which is what would hide a
  # missing guard on this one.
  expect_error(
    solve_cbps(
      covs,
      c(-1, -0.5, NaN, 1),
      treat,
      rep(1, 4),
      "ate",
      "logit",
      FALSE,
      FALSE,
      options
    ),
    "`covs_bal` contains a non-finite value at position 3",
    fixed = TRUE
  )
  expect_error(
    solve_cbps(
      covs,
      covs,
      treat,
      c(1, 1, 1, Inf),
      "ate",
      "logit",
      FALSE,
      FALSE,
      options
    ),
    "`s_weights` contains a non-finite value at position 4",
    fixed = TRUE
  )
})

test_that("solve_cbps_multi() refuses non-finite covariates and sampling weights", {
  covs <- c(-1, -0.5, 0.5, 1)
  treat_idx <- c(0L, 1L, 2L, 1L)
  options <- list(threads = 1L)
  expect_error(
    solve_cbps_multi(
      c(-1, -0.5, 0.5, Inf),
      treat_idx,
      0L,
      rep(1, 4),
      "ate",
      "logit",
      options
    ),
    "`covs` contains a non-finite value at position 4",
    fixed = TRUE
  )
  expect_error(
    solve_cbps_multi(
      covs,
      treat_idx,
      0L,
      c(1, 1, NaN, 1),
      "ate",
      "logit",
      options
    ),
    "`s_weights` contains a non-finite value at position 3",
    fixed = TRUE
  )
})

# The re-evaluation entry points take the same numeric blocks as the solves that
# produced them, at parameters a caller supplies, so they carry the same guards
# and are addressed the same way. Each spec plants one non-finite value in one
# guarded block and reads back which block the refusal names, since a guard
# pointed at the wrong argument would satisfy a bare check for the refusal
# alone.

test_that("the entropy re-evaluation entry points refuse non-finite blocks", {
  # One dual per group over a single constraint column: two groups of two units,
  # so `coefs` is one value per group and `covs` is the column itself.
  arguments <- list(
    coefs = c(0, 0),
    covs = c(-1, -1, 1, 1),
    group_idx = c(0L, 0L, 1L, 1L),
    targets = 0,
    base_weights = rep(1, 4),
    s_weights = rep(1, 4),
    n_eff = 4,
    esteq_scale = c(1, 1)
  )
  call_with <- function(entry, name, value) {
    modified <- arguments
    modified[[name]] <- value
    do.call(entry, modified)
  }

  for (entry in list(eval_psi_entropy, eval_weights_entropy)) {
    expect_error(
      call_with(entry, "covs", c(-1, -1, 1, Inf)),
      "`covs` contains a non-finite value at position 4",
      fixed = TRUE
    )
    expect_error(
      call_with(entry, "base_weights", c(1, NaN, 1, 1)),
      "`base_weights` contains a non-finite value at position 2",
      fixed = TRUE
    )
    expect_error(
      call_with(entry, "s_weights", c(1, 1, 1, Inf)),
      "`s_weights` contains a non-finite value at position 4",
      fixed = TRUE
    )
  }
})

test_that("the cbps re-evaluation entry points refuse non-finite blocks", {
  # The intercept column followed by one covariate, which is the design the
  # binary fit solves and re-evaluates at.
  arguments <- list(
    coefs = c(0, 0),
    covs = c(1, 1, 1, 1, -1, -0.5, 0.5, 1),
    treat = c(0L, 0L, 1L, 1L),
    s_weights = rep(1, 4),
    estimand = "ate",
    link = "logit"
  )
  call_with <- function(entry, name, value) {
    modified <- arguments
    modified[[name]] <- value
    do.call(entry, modified)
  }

  for (entry in list(eval_psi_cbps, eval_weights_cbps)) {
    expect_error(
      call_with(entry, "covs", c(1, 1, 1, 1, -1, -0.5, 0.5, Inf)),
      "`covs` contains a non-finite value at position 8",
      fixed = TRUE
    )
    expect_error(
      call_with(entry, "s_weights", c(1, 1, NaN, 1)),
      "`s_weights` contains a non-finite value at position 3",
      fixed = TRUE
    )
  }
})

test_that("solve_cbps_cont() refuses a non-finite exposure", {
  expect_error(
    solve_cbps_cont(
      c(-1, -0.5, 0.5, 1),
      c(-1, 0, 1, Inf),
      rep(1, 4),
      list(threads = 1L)
    ),
    "non-finite"
  )
})

# The quadratic-program entry points reached their solvers without a finiteness
# gate, so the boundary is made uniform with the estimating-equation family: every
# numeric block a solve reads is refused when it carries a non-finite value.

test_that("solve_energy() refuses non-finite covariates and sampling weights", {
  treat <- c(0L, 0L, 1L, 1L)
  empty <- numeric(0)
  options <- list(threads = 1L)
  expect_error(
    solve_energy(
      c(-1, -0.5, 0.5, Inf),
      treat,
      rep(1, 4),
      "scaled_euclidean",
      "ate",
      TRUE,
      empty,
      empty,
      empty,
      0,
      1e-4,
      options
    ),
    "non-finite"
  )
  expect_error(
    solve_energy(
      c(-1, -0.5, 0.5, 1),
      treat,
      c(1, 1, 1, NaN),
      "scaled_euclidean",
      "ate",
      TRUE,
      empty,
      empty,
      empty,
      0,
      1e-4,
      options
    ),
    "non-finite"
  )
})

test_that("solve_energy_cont() refuses a non-finite exposure", {
  empty <- numeric(0)
  bal_covs <- matrix(numeric(0), nrow = 4, ncol = 0)
  expect_error(
    solve_energy_cont(
      c(-1, -0.5, 0.5, 1),
      c(-1, 0, 1, Inf),
      rep(1, 4),
      "scaled_euclidean",
      TRUE,
      0,
      1e-4,
      matrix(c(-1, -0.5, 0.5, 1), ncol = 1),
      matrix(c(-1, 0, 1, 2), ncol = 1),
      bal_covs,
      empty,
      list(threads = 1L)
    ),
    "non-finite"
  )
})

test_that("solve_sbw() refuses non-finite sampling weights and moment columns", {
  treat <- c(0L, 0L, 1L, 1L)
  empty <- numeric(0)
  options <- list(threads = 1L)
  expect_error(
    solve_sbw(
      treat,
      c(1, 1, 1, Inf),
      "ate",
      "l2",
      empty,
      empty,
      empty,
      0,
      options
    ),
    "non-finite"
  )
  expect_error(
    solve_sbw(
      treat,
      rep(1, 4),
      "ate",
      "l2",
      c(-1, -0.5, 0.5, NaN),
      0,
      0.1,
      0,
      options
    ),
    "non-finite"
  )
})

test_that("solve_sbw_cont() refuses a non-finite exposure", {
  expect_error(
    solve_sbw_cont(
      c(-1, 0, 1, Inf),
      c(-1, -0.5, 0.5, 1),
      rep(1, 4),
      "l2",
      0.1,
      0,
      list(threads = 1L)
    ),
    "non-finite"
  )
})

test_that("solve_cfd() refuses non-finite covariates and projections", {
  treat <- c(0L, 0L, 1L, 1L)
  empty <- numeric(0)
  options <- list(threads = 1L, backend = "auto")
  expect_error(
    solve_cfd(
      c(-1, -0.5, 0.5, Inf),
      treat,
      rep(1, 4),
      "gaussian",
      1,
      1.5,
      empty,
      TRUE,
      "ate",
      empty,
      empty,
      empty,
      0,
      1e-4,
      options
    ),
    "non-finite"
  )
  # The t kernel's frequency projections are drawn on the R side, which is exactly
  # where an infinite degrees of freedom used to leave them all non-finite.
  expect_error(
    solve_cfd(
      c(-1, -0.5, 0.5, 1),
      treat,
      rep(1, 4),
      "t",
      1,
      1.5,
      c(NaN, NaN),
      TRUE,
      "ate",
      empty,
      empty,
      empty,
      0,
      1e-4,
      options
    ),
    "non-finite"
  )
})
