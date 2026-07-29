# Exposure-type detection resolves binary, categorical, and continuous
# exposures and honors an explicit type over the detection heuristics, refusing
# only a declaration the data cannot represent at all. These specs exercise the
# detection heuristics directly and the resolver's error and announcement
# branches.

# ---- detect_exposure_type -------------------------------------------------

test_that("a two-valued exposure is binary regardless of storage type", {
  expect_identical(detect_exposure_type(c(0L, 1L, 0L, 1L)), "binary")
  expect_identical(detect_exposure_type(c(TRUE, FALSE, TRUE)), "binary")
  expect_identical(
    detect_exposure_type(factor(c("a", "b", "a"))),
    "binary"
  )
})

test_that("a many-level factor or character exposure is categorical", {
  expect_identical(
    detect_exposure_type(factor(c("a", "b", "c", "a"))),
    "categorical"
  )
  expect_identical(
    detect_exposure_type(c("a", "b", "c", "a")),
    "categorical"
  )
})

test_that("a numeric exposure with few distinct values is categorical", {
  # Five distinct values in 100 observations is a 5 percent unique share, below
  # the 20 percent categorical threshold.
  withr::local_seed(1)
  exposure <- sample(1:5, 100, replace = TRUE)
  expect_identical(detect_exposure_type(exposure), "categorical")
})

test_that("a numeric exposure with many distinct values is continuous", {
  withr::local_seed(1)
  expect_identical(detect_exposure_type(stats::rnorm(100)), "continuous")
})

# ---- is_categorical guards ------------------------------------------------

test_that("is_categorical() is FALSE when every value is missing", {
  expect_false(is_categorical(c(NA_real_, NA_real_)))
})

# ---- resolve_exposure_type ------------------------------------------------

test_that("an explicit categorical type is honored on a factor exposure", {
  resolved <- resolve_exposure_type(
    "categorical",
    factor(c("a", "b", "c", "a")),
    bw_entropy()
  )
  expect_identical(resolved, "categorical")
})

test_that("an explicit continuous type is honored on a continuous exposure", {
  withr::local_seed(1)
  resolved <- resolve_exposure_type(
    "continuous",
    stats::rnorm(100),
    bw_entropy()
  )
  expect_identical(resolved, "continuous")
})

test_that("an explicit continuous type wins over the categorical heuristic", {
  # Ten distinct doses in 200 observations is a 5 percent unique share, which the
  # heuristic reads as categorical. An explicit type is the caller's declaration
  # of how the exposure is modeled, so it decides the fit.
  exposure <- withr::with_seed(
    1,
    sample(seq(10, 100, by = 10), 200, replace = TRUE)
  )
  expect_identical(detect_exposure_type(exposure), "categorical")
  expect_identical(
    resolve_exposure_type("continuous", exposure, bw_entropy()),
    "continuous"
  )
})

test_that("an explicit continuous type wins over two-level detection", {
  # A numeric exposure taking two values is still a dose the caller may model as
  # continuous, so the declaration stands.
  expect_identical(
    resolve_exposure_type("continuous", rep(c(0, 1), 50), bw_entropy()),
    "continuous"
  )
})

test_that("an explicit categorical type wins over continuous detection", {
  exposure <- withr::with_seed(1, stats::rnorm(20))
  expect_identical(
    resolve_exposure_type("categorical", exposure, bw_entropy()),
    "categorical"
  )
})

test_that("auto resolution announces the detected type", {
  withr::local_options(balancing.quiet = FALSE)
  expect_message(
    resolved <- resolve_exposure_type(
      "auto",
      c(0L, 1L, 0L, 1L),
      bw_entropy()
    ),
    "binary"
  )
  expect_identical(resolved, "binary")
})

test_that("a forced type the data contradict raises a classed error", {
  withr::local_seed(1)
  expect_error(
    resolve_exposure_type(
      "binary",
      stats::rnorm(100),
      bw_entropy()
    ),
    class = "balancing_exposure_type_error"
  )
})

test_that("a continuous type on a non-numeric exposure is a classed error", {
  # Structural impossibility rather than heuristic disagreement: a factor or
  # character exposure carries no dose to correlate the covariates against.
  expect_error(
    resolve_exposure_type(
      "continuous",
      factor(c("a", "b", "c", "a")),
      bw_entropy()
    ),
    class = "balancing_exposure_type_error"
  )
  expect_error(
    resolve_exposure_type(
      "continuous",
      c("a", "b", "c", "a"),
      bw_entropy()
    ),
    class = "balancing_exposure_type_error"
  )
})
