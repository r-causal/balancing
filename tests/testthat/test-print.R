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
