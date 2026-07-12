# Exposure-type detection resolves binary, categorical, and continuous
# exposures and honors an explicit type only when the data support it. These
# specs exercise the detection heuristics directly and the resolver's error and
# announcement branches.

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
    entropy_balance()
  )
  expect_identical(resolved, "categorical")
})

test_that("an explicit continuous type is honored on a continuous exposure", {
  withr::local_seed(1)
  resolved <- resolve_exposure_type(
    "continuous",
    stats::rnorm(100),
    entropy_balance()
  )
  expect_identical(resolved, "continuous")
})

test_that("auto resolution announces the detected type", {
  withr::local_options(balancing.quiet = FALSE)
  expect_message(
    resolved <- resolve_exposure_type(
      "auto",
      c(0L, 1L, 0L, 1L),
      entropy_balance()
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
      entropy_balance()
    ),
    class = "balancing_exposure_type_error"
  )
})
