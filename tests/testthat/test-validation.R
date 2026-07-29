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
    "non-finite"
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
    "non-finite"
  )
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
