# S7 property validators reject bad values at construction, and the abstract
# parent classes cannot be instantiated directly. Each test first exercises a
# valid construction so the spec depends on the real implementation, then
# asserts the rejection.

# ---- entropy_balance validators -------------------------------------------

test_that("entropy_balance() rejects a non-positive convergence tolerance", {
  expect_identical(entropy_balance()@convergence_tolerance, 1e-10)
  expect_error(entropy_balance(convergence_tolerance = -1e-10))
})

test_that("entropy_balance() rejects a negative iteration cap", {
  expect_null(entropy_balance()@max_iterations)
  expect_error(entropy_balance(max_iterations = -5L))
})

test_that("entropy_balance() rejects negative base weights", {
  expect_equal(entropy_balance(base_weights = c(1, 2))@base_weights, c(1, 2))
  expect_error(entropy_balance(base_weights = c(1, -1, 1)))
})

# ---- balance_terms validators ---------------------------------------------

test_that("balance_terms() rejects a negative tolerance", {
  expect_identical(balance_terms()@tolerance, 0)
  expect_error(balance_terms(tolerance = -0.5))
})

test_that("balance_terms() rejects negative moments", {
  expect_identical(balance_terms(moments = 2L)@moments, 2L)
  expect_error(balance_terms(moments = -2L))
})

test_that("balance_terms() rejects quantiles outside the open unit interval", {
  expect_identical(balance_terms(quantiles = 0.5)@quantiles, 0.5)
  expect_error(balance_terms(quantiles = c(0, 0.5)))
  expect_error(balance_terms(quantiles = c(0.5, 1)))
})

# ---- balance_terms: named and list variants -------------------------------

test_that("balance_terms() accepts a named tolerance vector", {
  terms <- balance_terms(tolerance = c(x1 = 0.1, x2 = 0.2))
  expect_identical(terms@tolerance, c(x1 = 0.1, x2 = 0.2))
})

test_that("balance_terms() accepts a quantile list and validates its elements", {
  terms <- balance_terms(quantiles = list(x1 = c(0.25, 0.75), x2 = 0.5))
  expect_type(terms@quantiles, "list")
  expect_identical(terms@quantiles$x1, c(0.25, 0.75))
  expect_error(balance_terms(quantiles = list(x1 = c(0.5, 1.5))))
})

# ---- Abstract parents ------------------------------------------------------

test_that("the abstract method parents are not constructible", {
  # A concrete subclass constructs; the abstract parents do not.
  expect_true(S7::S7_inherits(entropy_balance(), balance_method))
  expect_error(balance_method())
  expect_error(estimating_equation_method())
  expect_error(quadratic_program_method())
})
