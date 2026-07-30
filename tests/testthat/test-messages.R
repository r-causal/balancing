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
      method = bw_entropy(max_iterations = 3L),
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
      method = bw_entropy(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0.1)
    )
  )
})

test_that("balancing_ignored_argument_warning: two_step without over_identified", {
  # The two-step weighting matrix belongs to the over-identified criterion, so
  # requesting it on a just-identified fit has no effect; the fit warns that the
  # argument is ignored and proceeds.
  data <- sim_binary()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(two_step = FALSE, over_identified = FALSE),
      estimand = "ate"
    )
  )
})

test_that("balancing_ignored_argument_warning: over_identified for a categorical exposure", {
  # The over-identified criterion is defined for a binary exposure alone, so a
  # categorical fit warns that the request is ignored and returns the exactly
  # balancing solution. The message names the exposure type it was raised for.
  data <- sim_categorical()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(over_identified = TRUE),
      estimand = "ate"
    )
  )
})

test_that("balancing_ignored_argument_warning: over_identified for a continuous exposure", {
  data <- sim_continuous()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(over_identified = TRUE),
      estimand = "ate"
    )
  )
})

test_that("balancing_ignored_argument_warning: every argument a categorical fit ignores", {
  # Once the over-identified request is ignored the fit is not over-identified,
  # so the two-step weighting matrix has no criterion to weight either. Each
  # ignored argument carries its own warning rather than the first standing in
  # for both.
  data <- sim_categorical()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(over_identified = TRUE, two_step = FALSE),
      estimand = "ate"
    )
  )
})

test_that("balancing_ignored_argument_warning: a clarabel pin the energy kernel cannot honor", {
  # The energy kernel's quadratic term is indefinite, which the interior-point
  # backend refuses, so a pinned clarabel request cannot be honored for it. The
  # fit names the request it dropped and the backend that ran instead.
  withr::local_options(balancing.qp_backend = "clarabel")
  data <- sim_binary()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cfd(kernel = "energy"),
      estimand = "ate"
    )
  )
})

test_that("balancing_ignored_argument_warning: a clarabel pin energy balancing cannot honor", {
  # Energy balancing assembles the same indefinite quadratic form as the energy
  # kernel, so the interior-point backend refuses it and a pinned clarabel
  # request cannot be honored. The fit names the request it dropped and the
  # backend that ran instead.
  withr::local_options(balancing.qp_backend = "clarabel")
  data <- sim_binary()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(),
      estimand = "ate"
    )
  )
})

test_that("balancing_ignored_argument_warning: focal_level with a pooled estimand", {
  # The average treatment effect reweights every exposure group rather than
  # holding one fixed, so it has no focal level to resolve and a supplied one is
  # never validated against the data. The fit names the estimand that ignores it.
  data <- sim_binary()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      estimand = "ate",
      focal_level = 1
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
      method = bw_entropy(),
      estimand = "ate"
    ))
  )
})

test_that("alert: the exposure is excluded from a covariate selection", {
  withr::local_options(balancing.quiet = FALSE)
  data <- sim_binary()
  expect_snapshot(
    invisible(balance(
      data,
      exposure,
      everything(),
      method = bw_entropy(),
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
      method = bw_entropy(),
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
      method = bw_entropy(),
      estimand = "ate",
      exposure_type = "binary",
      constraints = balance_terms(moments = 2L)
    ))
  )
})
