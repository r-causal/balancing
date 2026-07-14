# glance() reports a one-row fit summary. The imbalance column is named for the
# balance statistic the exposure type uses: a standardized mean difference for
# discrete exposures and a correlation for continuous ones.

test_that("glance() of a binary fit reports the max absolute SMD", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bal_entropy(),
    estimand = "ate"
  )
  glanced <- generics::glance(fit)
  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 1L)
  expect_true("max_absolute_smd" %in% names(glanced))
  expect_false("max_absolute_correlation" %in% names(glanced))
  expect_identical(glanced$method, "Entropy balancing")
  expect_identical(glanced$estimand, "ate")
  expect_identical(glanced$exposure_type, "binary")
  expect_gt(glanced$ess, 0)
})

test_that("glance() of a continuous fit reports the max absolute correlation", {
  data <- sim_continuous()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bal_entropy(),
    estimand = "ate"
  )
  glanced <- generics::glance(fit)
  expect_true("max_absolute_correlation" %in% names(glanced))
  expect_false("max_absolute_smd" %in% names(glanced))
  expect_identical(glanced$exposure_type, "continuous")
})
