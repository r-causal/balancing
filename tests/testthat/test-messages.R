# Snapshot the user-facing warnings and informational alerts the entropy slice
# raises. The classed fit-time warnings, the class-downgrade warning, and the
# three covariate-expansion alerts are recorded with their message text and
# condition class so a change in wording or class is caught.

# ---- Fit-time warnings ----------------------------------------------------

test_that("balancing_convergence_warning: the iteration cap is reached", {
  # Three Newton steps drive the balance essentially to zero but do not meet the
  # gradient tolerance, so the fit warns about convergence alone.
  data <- sim_binary()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = entropy_balance(max_iterations = 3L),
      estimand = "ate"
    )
  )
})

test_that("balancing_balance_warning: achieved balance exceeds the tolerance", {
  # A continuous tolerance without the second distribution moment leaves the
  # exposure variance free, so the weighted correlation exceeds the requested
  # bound and the fit warns.
  data <- sim_continuous()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = entropy_balance(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0.1)
    )
  )
})

# ---- Class-downgrade warning ----------------------------------------------

test_that("balancing_class_downgrade_warning: mismatched estimands", {
  x <- bw(c(1, 2), estimand = "ate")
  y <- bw(c(3, 4), estimand = "att")
  expect_balancing_warning(vctrs::vec_c(x, y))
})

# ---- Covariate-expansion alerts -------------------------------------------

test_that("alert: the detected exposure type is announced", {
  withr::local_options(balancing.quiet = FALSE)
  data <- sim_binary()
  expect_snapshot(
    invisible(balance(
      data,
      exposure,
      c(x1, x2),
      method = entropy_balance(),
      estimand = "ate"
    ))
  )
})

test_that("alert: aliased constraint columns are dropped", {
  withr::local_options(balancing.quiet = FALSE)
  data <- sim_binary()
  data$x1_copy <- data$x1
  expect_snapshot(
    invisible(balance(
      data,
      exposure,
      c(x1, x1_copy, x2),
      method = entropy_balance(),
      estimand = "ate",
      exposure_type = "binary"
    ))
  )
})

test_that("alert: moments above one on a binary covariate are ignored", {
  withr::local_options(balancing.quiet = FALSE)
  data <- sim_binary()
  data$flag <- as.integer(data$x1 > 0)
  expect_snapshot(
    invisible(balance(
      data,
      exposure,
      c(x2, flag),
      method = entropy_balance(),
      estimand = "ate",
      exposure_type = "binary",
      constraints = balance_terms(moments = 2L)
    ))
  )
})
