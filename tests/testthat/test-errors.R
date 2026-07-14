# Snapshot the message and condition class of every classed error reachable in
# the entropy slice. `expect_balancing_error()` records with cnd_class = TRUE so
# the subclass is captured alongside the text. Snapshots record on the first
# successful run once the implementation exists.

test_that("balancing_type_error: non-data-frame input", {
  expect_balancing_error(
    balance(list(a = 1), exposure, c(x1, x2), method = bal_entropy())
  )
})

test_that("balancing_selection_error: exposure selects two columns", {
  data <- sim_binary(n = 100)
  expect_balancing_error(
    balance(data, c(x1, x2), c(x1, x2), method = bal_entropy())
  )
})

test_that("balancing_missing_error: missing covariate values", {
  data <- sim_binary(n = 100)
  data$x1[1] <- NA
  expect_balancing_error(
    balance(data, exposure, c(x1, x2), method = bal_entropy())
  )
})

test_that("balancing_method_error: a bare-string method", {
  data <- sim_binary(n = 100)
  expect_balancing_error(
    balance(data, exposure, c(x1, x2), method = "entropy")
  )
})

test_that("balancing_estimand_error: an unsupported estimand", {
  data <- sim_binary(n = 100)
  expect_balancing_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bal_entropy(),
      estimand = "ato"
    )
  )
})

test_that("balancing_exposure_type_error: a method rejects an exposure type", {
  data <- sim_continuous(n = 150)
  expect_balancing_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bal_ipt(),
      estimand = "ate"
    )
  )
})

test_that("balancing_estimand_error: a categorical att without focal_level", {
  data <- sim_categorical(n = 150)
  expect_balancing_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bal_entropy(),
      estimand = "att"
    )
  )
})

test_that("balancing_exposure_type_error: a forced type contradicts the data", {
  data <- sim_continuous(n = 150)
  expect_balancing_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bal_entropy(),
      exposure_type = "binary"
    )
  )
})

test_that("balancing_constraints_error: quantiles with a continuous exposure", {
  data <- sim_continuous(n = 150)
  expect_balancing_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bal_entropy(),
      estimand = "ate",
      constraints = balance_terms(quantiles = 0.5)
    )
  )
})

test_that("balancing_constraints_error: an unnamed multi-element tolerance", {
  data <- sim_binary(n = 100)
  expect_balancing_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bal_entropy(),
      constraints = balance_terms(tolerance = c(0.1, 0.2))
    )
  )
})

test_that("balancing_constraints_error: a tolerance named for a non-covariate", {
  data <- sim_binary(n = 100)
  expect_balancing_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bal_entropy(),
      constraints = balance_terms(tolerance = c(nonesuch = 0.1))
    )
  )
})

test_that("balancing_range_error: base weights of the wrong length", {
  data <- sim_binary(n = 100)
  expect_balancing_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bal_entropy(base_weights = rep(1, 3))
    )
  )
})

test_that("balancing_range_error: an invalid entropy solver option", {
  data <- sim_binary(n = 100)
  withr::local_options(balancing.entropy_solver = "nope")
  expect_balancing_error(
    balance(data, exposure, c(x1, x2), method = bal_entropy())
  )
})

test_that("balancing_range_error: an invalid quadratic-program backend option", {
  data <- sim_binary(n = 100)
  withr::local_options(balancing.qp_backend = "nope")
  expect_balancing_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bal_sbw(),
      constraints = balance_terms(tolerance = 0.05)
    )
  )
})

test_that("balancing_ipw_unsupported_error: estimating_equations() when absent", {
  data <- sim_binary(n = 150)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bal_entropy(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  expect_balancing_error(estimating_equations(fit))
})

test_that("balancing_autoplot_duals_error: the duals view without dual variables", {
  data <- sim_binary(n = 150)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bal_entropy(),
    estimand = "ate"
  )
  # The estimating-equation family carries no dual variables, so the dual view
  # aborts before it touches ggplot2.
  expect_balancing_error(autoplot_duals(fit))
})
