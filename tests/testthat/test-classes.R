# S7 property validators reject bad values at construction, and the abstract
# parent classes cannot be instantiated directly. Each test first exercises a
# valid construction so the spec depends on the real implementation, then
# asserts the rejection.

# ---- bw_entropy validators -------------------------------------------

test_that("bw_entropy() rejects a non-positive convergence tolerance", {
  expect_identical(bw_entropy()@convergence_tolerance, 1e-10)
  expect_error(bw_entropy(convergence_tolerance = -1e-10))
})

test_that("bw_entropy() rejects a negative iteration cap", {
  expect_null(bw_entropy()@max_iterations)
  expect_error(bw_entropy(max_iterations = -5L))
})

test_that("bw_entropy() rejects negative base weights", {
  expect_equal(bw_entropy(base_weights = c(1, 2))@base_weights, c(1, 2))
  expect_error(bw_entropy(base_weights = c(1, -1, 1)))
})

# A missing base weight used to reach the sign comparison and stop the validator
# itself, and an infinity used to pass construction and fail mid-solve, so each
# now reports at construction with a message naming the property.
test_that("bw_entropy() rejects missing base weights", {
  expect_equal(bw_entropy(base_weights = c(1, 2))@base_weights, c(1, 2))
  expect_error(
    bw_entropy(base_weights = c(1, NA, 1)),
    "must not contain missing values"
  )
  expect_error(
    bw_entropy(base_weights = c(1, NaN, 1)),
    "must not contain missing values"
  )
})

test_that("bw_entropy() rejects infinite base weights", {
  expect_error(
    bw_entropy(base_weights = c(1, Inf, 1)),
    "must not contain infinite values"
  )
  expect_error(
    bw_entropy(base_weights = c(1, -Inf, 1)),
    "must not contain infinite values"
  )
})

# An all-zero base measure leaves every constraint target undefined, so it is
# refused at construction, which covers the discrete and the continuous fit paths
# alike. Individual zero base weights stay legal.
test_that("bw_entropy() rejects an all-zero base weight vector", {
  expect_equal(bw_entropy(base_weights = c(0, 1, 2))@base_weights, c(0, 1, 2))
  expect_error(bw_entropy(base_weights = rep(0, 3)), "must not be all zero")
})

# Entropy balancing carried no distribution-moments clause at all, so a missing
# value or a value below one passed construction and reached the continuous fit,
# which stopped on a base comparison. Energy balancing validates the same property
# and both now refuse the same set at construction.
test_that("bw_entropy() validates the distribution moments", {
  expect_identical(
    bw_entropy(distribution_moments = 2L)@distribution_moments,
    2L
  )
  expect_error(bw_entropy(distribution_moments = NA_integer_), "single")
  expect_error(bw_entropy(distribution_moments = c(1L, 2L)), "single")
  expect_error(bw_entropy(distribution_moments = 0L), "positive")
  expect_error(bw_entropy(distribution_moments = -1L), "positive")
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

# A missing value used to reach the range comparisons and stop the validator
# itself with a base error, so nothing named the property at fault. Each of the
# three range-checked properties reports the missingness instead.
test_that("balance_terms() rejects missing values in every range property", {
  expect_error(balance_terms(moments = NA), "missing values")
  expect_error(balance_terms(moments = c(x1 = 2L, x2 = NA_integer_)), "missing")
  expect_error(balance_terms(tolerance = NA_real_), "missing values")
  expect_error(balance_terms(tolerance = c(x1 = 0.1, x2 = NA)), "missing")
  expect_error(balance_terms(quantiles = NA_real_), "missing values")
  expect_error(balance_terms(quantiles = c(0.25, NA)), "missing values")
  expect_error(balance_terms(quantiles = list(x1 = c(0.5, NA))), "missing")
})

# The tolerance is the property callers tune most, and its three unusable inputs
# are turned away by three different mechanisms, none of them the package's own.
# `NULL` reaches the S7 property type check, which reports a class mismatch
# rather than a tolerance that has no value; `Inf` is neither missing nor
# negative, so it passes every check and constructs a specification whose box is
# unbounded; a missing value is caught, but by the validator, which raises the
# base S7 condition. Every other refusal a caller can trigger from a balance
# specification carries `balancing_constraints_error`, so these do too. Each is
# pinned by class rather than by wording, and separately, since a test stops at
# the first input it refuses the wrong way.
test_that("balance_terms() classes its refusal of an absent tolerance", {
  expect_error(
    balance_terms(tolerance = NULL),
    class = "balancing_constraints_error"
  )
})

test_that("balance_terms() classes its refusal of an infinite tolerance", {
  expect_error(
    balance_terms(tolerance = Inf),
    class = "balancing_constraints_error"
  )
  expect_error(
    balance_terms(tolerance = c(x1 = 0.1, x2 = Inf)),
    class = "balancing_constraints_error"
  )
})

test_that("balance_terms() classes its refusal of a missing tolerance", {
  expect_error(
    balance_terms(tolerance = NA_real_),
    class = "balancing_constraints_error"
  )
})

# ---- balance_method validators --------------------------------------------

# The optional solver tuning parameters are validated for being a single usable
# number, not only for their sign. A missing value or a vector of the wrong length
# used to steer the sign comparison and stop the validator with a base error, so
# the property was never named. Every method inherits the check from the abstract
# parent, so one specification covers all six.
test_that("a missing, non-finite, or multi-element tuning value is refused", {
  constructors <- list(
    bw_entropy = bw_entropy,
    bw_ipt = bw_ipt,
    bw_cbps = bw_cbps,
    bw_energy = bw_energy,
    bw_cfd = bw_cfd,
    bw_sbw = bw_sbw
  )
  for (name in names(constructors)) {
    constructor <- constructors[[name]]
    expect_identical(
      constructor(convergence_tolerance = 1e-8)@convergence_tolerance,
      1e-8
    )
    expect_error(constructor(convergence_tolerance = NA_real_), "finite")
    expect_error(constructor(convergence_tolerance = NaN), "finite")
    expect_error(constructor(convergence_tolerance = Inf), "finite")
    expect_error(constructor(convergence_tolerance = c(1e-8, 1e-9)), "single")
    expect_identical(constructor(max_iterations = 50L)@max_iterations, 50L)
    expect_error(constructor(max_iterations = NA_integer_), "single")
    expect_error(constructor(max_iterations = c(10L, 20L)), "single")
  }
})

# The five methods that cast their tuning arguments accept a bare numeric literal;
# entropy balancing used to be the exception, failing the S7 property type check
# on the integer iteration cap and the double tolerance alike. All six now take
# whichever numeric type the caller typed.
test_that("every method casts bare numeric tuning literals", {
  constructors <- list(
    bw_entropy = bw_entropy,
    bw_ipt = bw_ipt,
    bw_cbps = bw_cbps,
    bw_energy = bw_energy,
    bw_cfd = bw_cfd,
    bw_sbw = bw_sbw
  )
  for (name in names(constructors)) {
    constructor <- constructors[[name]]
    expect_identical(constructor(max_iterations = 200)@max_iterations, 200L)
    expect_identical(
      constructor(convergence_tolerance = 1L)@convergence_tolerance,
      1
    )
  }
})

# All six constructors document that their tuning parameters must be passed by
# name, and the position of the dots is what enforces it. With the dots declared
# last, a bare argument bound positionally to whichever tuning parameter came
# first, so `bw_ipt("probit")` set the link and `bw_energy("mahalanobis")` set the
# distance without either being named, while entropy balancing refused the same
# call. Dots-first makes all six refuse it through the classed method error.
test_that("no method constructor binds a positional argument", {
  constructors <- list(
    bw_entropy = bw_entropy,
    bw_ipt = bw_ipt,
    bw_cbps = bw_cbps,
    bw_energy = bw_energy,
    bw_cfd = bw_cfd,
    bw_sbw = bw_sbw
  )
  for (name in names(constructors)) {
    constructor <- constructors[[name]]
    expect_error(constructor(1e-4), class = "balancing_method_error")
    expect_error(constructor("probit"), class = "balancing_method_error")
    expect_error(constructor(TRUE), class = "balancing_method_error")
    # The named form of the same call still constructs.
    expect_identical(constructor(max_iterations = 50L)@max_iterations, 50L)
  }
})

# ---- supports_estimating_equations() context arguments ---------------------

# The generic takes its context through the dots, so every method has to declare
# the same context arguments for a caller to be able to ask any of them the same
# question. Entropy balancing and inverse probability tilting declared only
# `constraints`, so `exposure_type` fell through to their `check_dots_empty()` and
# the call errored for two of the six.
test_that("every method answers with either context argument", {
  methods <- list(
    bw_entropy(),
    bw_ipt(),
    bw_cbps(),
    bw_energy(),
    bw_cfd(),
    bw_sbw()
  )
  for (method in methods) {
    for (type in c("binary", "categorical", "continuous")) {
      answer <- supports_estimating_equations(method, exposure_type = type)
      expect_type(answer, "logical")
      expect_length(answer, 1L)
      expect_false(is.na(answer))
      combined <- supports_estimating_equations(
        method,
        exposure_type = type,
        constraints = balance_terms(tolerance = 0)
      )
      expect_type(combined, "logical")
      expect_length(combined, 1L)
    }
    expect_type(
      supports_estimating_equations(
        method,
        constraints = balance_terms(tolerance = 0)
      ),
      "logical"
    )
    # An argument no method consumes is still refused, so a misspelled context
    # name does not pass silently.
    expect_error(supports_estimating_equations(method, bogus = 1))
  }
})

# Entropy balancing's continuous fit solves smooth estimating equations exactly
# when every requested tolerance is zero, the same rule its discrete fit follows.
# No other method produces them for a continuous exposure at all.
test_that("only entropy balancing answers TRUE for a continuous exposure", {
  expect_true(supports_estimating_equations(
    bw_entropy(),
    exposure_type = "continuous"
  ))
  expect_true(supports_estimating_equations(
    bw_entropy(),
    exposure_type = "continuous",
    constraints = balance_terms(tolerance = 0)
  ))
  expect_false(supports_estimating_equations(
    bw_entropy(),
    exposure_type = "continuous",
    constraints = balance_terms(tolerance = 0.05)
  ))
  for (method in list(bw_ipt(), bw_cbps(), bw_energy(), bw_cfd(), bw_sbw())) {
    expect_false(supports_estimating_equations(
      method,
      exposure_type = "continuous"
    ))
  }
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
  expect_true(S7::S7_inherits(bw_entropy(), balance_method))
  expect_error(balance_method())
  expect_error(estimating_equation_method())
  expect_error(quadratic_program_method())
})
