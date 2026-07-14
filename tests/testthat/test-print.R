# Print and summary snapshots for a fitted entropy result. The snapshots record
# on the first successful run once the implementation exists.

test_that("print() of a binary ate fit is stable", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bal_entropy(),
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
    method = bal_entropy(),
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
    method = bal_entropy(),
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
    method = bal_entropy(),
    estimand = "ate"
  )
  expect_snapshot(print(fit))
})

test_that("tidy() of a binary ate fit returns one row per constraint term", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bal_entropy(),
    estimand = "ate"
  )
  tidied <- generics::tidy(fit)
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
      names(tidied)
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
    method = bal_entropy(),
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
    method = bal_energy(),
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
    method = bal_energy(),
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
    method = bal_sbw(),
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
    method = bal_cfd(),
    estimand = "ate"
  )
  expect_snapshot(summary(fit))
})

test_that("summary() returns the fit invisibly", {
  data <- sim_binary(200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bal_entropy(),
    estimand = "ate"
  )
  returned <- withr::with_output_sink(tempfile(), summary(fit))
  expect_true(S7::S7_inherits(returned, balancing))
  expect_identical(returned@estimand, "ate")
})
