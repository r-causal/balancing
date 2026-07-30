# Print and summary snapshots for a fitted entropy result. The snapshots record
# on the first successful run once the implementation exists.

test_that("print() of a binary ate fit is stable", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_snapshot(print(fit))
})

test_that("summary() of a binary ate fit is stable", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_snapshot(summary(fit))
})

test_that("print() of a binary att fit renders the focal level", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att"
  )
  expect_snapshot(print(fit))
})

test_that("print() of a continuous ate fit is stable", {
  data <- sim_continuous()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_snapshot(print(fit))
})

test_that("the balance table carries one row per constraint term", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  balance_table <- fit@balance_table
  expect_true(all(
    c(
      "term",
      "kind",
      "statistic",
      "group",
      "unweighted",
      "weighted",
      "tolerance",
      "within_tolerance"
    ) %in%
      names(balance_table)
  ))
})

# A categorical exposure reports one effective sample size per level, so its
# print block exercises the multi-group path.
test_that("print() of a categorical ate fit lists every level", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_snapshot(print(fit))
})

# The quadratic-program family reports a solver backend and a minimum-weight
# floor, so its print and summary blocks differ from the estimating-equation
# family. The weight summary names the count of weights resting on the floor.
test_that("print() of an energy fit is stable", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  expect_snapshot(print(fit))
})

test_that("summary() of an energy fit reports the weight floor count", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  expect_snapshot(summary(fit))
})

test_that("summary() of a stable balancing fit reports the weight floor count", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  expect_snapshot(summary(fit))
})

test_that("summary() of a cfd fit reports the weight floor count", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
    estimand = "ate"
  )
  expect_snapshot(summary(fit))
})

# ---- An imbalance that was never measured ----------------------------------

# A balance statistic that is not a number states no distance, so the display
# reports that the imbalance could not be assessed rather than offering NaN as the
# largest one. That is the ruling the balance warning already follows, applied to
# the print and summary blocks. No fit reaches this once every exposure level is
# required to carry base-measure mass, so the table is edited to reach it.
undefined_balance_fit <- function(value) {
  data <- sim_binary(n = 150)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  table <- fit@balance_table
  table$weighted[[1]] <- value
  fit@balance_table <- table
  fit
}

test_that("print() reports an unmeasurable imbalance rather than NaN", {
  fit <- undefined_balance_fit(NaN)
  printed <- utils::capture.output(expect_no_warning(print(fit)))
  expect_false(any(grepl("NaN", printed, fixed = TRUE)))
  expect_true(any(grepl("could not be assessed", printed, fixed = TRUE)))
})

test_that("print() reports a missing imbalance rather than NA", {
  fit <- undefined_balance_fit(NA_real_)
  printed <- utils::capture.output(expect_no_warning(print(fit)))
  expect_false(any(grepl("NA", printed, fixed = TRUE)))
  expect_true(any(grepl("could not be assessed", printed, fixed = TRUE)))
})

test_that("summary() reports an unmeasurable imbalance rather than NaN", {
  fit <- undefined_balance_fit(NaN)
  printed <- utils::capture.output(expect_no_warning(summary(fit)))
  expect_true(any(grepl("could not be assessed", printed, fixed = TRUE)))
  # The balance table itself still shows the term's value; only the headline
  # figure, which claims to measure how far the fit missed, is withheld.
  expect_false(any(grepl("Largest imbalance: NaN", printed, fixed = TRUE)))
})

test_that("summary() returns the fit invisibly", {
  data <- sim_binary(200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  returned <- withr::with_output_sink(tempfile(), summary(fit))
  expect_true(S7::S7_inherits(returned, balancing))
  expect_identical(returned@estimand, "ate")
})
