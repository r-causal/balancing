# balance_terms() records the constraint set; build_constraint_matrix() turns a
# covariate selection plus a balance_terms() spec into the numeric constraint
# matrix and a serializable recipe following the expansion rules in the design;
# rebuild_constraint_matrix() reconstructs the matrix from the recipe alone.
#
# build_constraint_matrix() is an internal helper called after balance()
# resolves the covariate names, so these specs call it with a data frame, a
# character vector of covariate names, a balance_terms() spec, and the resolved
# exposure type, and expect a list with `matrix` and `recipe` components.

# ---- balance_terms() defaults and validation ------------------------------

test_that("balance_terms() carries its documented defaults", {
  terms <- balance_terms()
  expect_null(terms@moments)
  expect_false(terms@interactions)
  expect_null(terms@quantiles)
  expect_identical(terms@tolerance, 0)
})

test_that("balance_terms() stores supplied values", {
  terms <- balance_terms(
    moments = 2L,
    interactions = TRUE,
    quantiles = c(0.25, 0.75),
    tolerance = 0.1
  )
  expect_identical(terms@moments, 2L)
  expect_true(terms@interactions)
  expect_identical(terms@quantiles, c(0.25, 0.75))
  expect_identical(terms@tolerance, 0.1)
})

test_that("balance_terms() rejects negative moments", {
  expect_identical(balance_terms(moments = 1L)@moments, 1L)
  expect_error(balance_terms(moments = -1L))
})

test_that("balance_terms() rejects quantiles outside the unit interval", {
  expect_identical(balance_terms(quantiles = 0.75)@quantiles, 0.75)
  expect_error(balance_terms(quantiles = c(0.5, 1.5)))
})

test_that("balance_terms() rejects a negative tolerance", {
  expect_identical(balance_terms(tolerance = 0.1)@tolerance, 0.1)
  expect_error(balance_terms(tolerance = -0.1))
})

# ---- build_constraint_matrix(): base columns and factors ------------------

test_that("numeric covariates become one mean-balance column each", {
  # The two covariates are linearly independent of one another, so both survive
  # the aliasing check and each contributes exactly one column.
  data <- data.frame(x1 = c(-1, 0, 1, 2), x2 = c(2, 1, 0, 5))
  built <- build_constraint_matrix(
    data,
    c("x1", "x2"),
    balance_terms(),
    exposure_type = "binary"
  )

  expect_identical(ncol(built$matrix), 2L)
  kinds <- vapply(built$recipe, function(term) term$kind, character(1))
  expect_true(all(kinds == "moment"))
})

test_that("factor covariates expand to a full set of level indicators", {
  data <- data.frame(
    x1 = c(-1, 0, 1, 2, 0.5, -0.5),
    f = factor(c("a", "b", "c", "a", "b", "c"))
  )
  built <- build_constraint_matrix(
    data,
    c("x1", "f"),
    balance_terms(),
    exposure_type = "binary"
  )

  sources <- vapply(built$recipe, function(term) term$source, character(1))
  # One column for x1 plus one indicator per level of f (full set, not
  # reference-coded).
  expect_identical(sum(sources == "f"), 3L)
  expect_identical(ncol(built$matrix), 4L)
})

# ---- build_constraint_matrix(): powers and standardization ----------------

test_that("moments above one add centered raw powers", {
  data <- data.frame(x1 = c(-2, -1, 0, 1, 2, 3))
  built <- build_constraint_matrix(
    data,
    "x1",
    balance_terms(moments = 3L),
    exposure_type = "binary"
  )

  kinds <- vapply(built$recipe, function(term) term$kind, character(1))
  powers <- vapply(built$recipe, function(term) term$power, integer(1))
  expect_identical(ncol(built$matrix), 3L)
  expect_identical(sort(powers), c(1L, 2L, 3L))
  expect_identical(kinds[powers == 1L], "moment")
  expect_true(all(kinds[powers > 1L] == "power"))
})

test_that("the recipe records the standardization center and scale", {
  data <- data.frame(x1 = c(-2, -1, 0, 1, 2, 3))
  built <- build_constraint_matrix(
    data,
    "x1",
    balance_terms(),
    exposure_type = "binary"
  )

  record <- built$recipe[[1]]
  expect_true(is.numeric(record$center))
  expect_true(is.numeric(record$scale))
  # Columns cross the boundary standardized to unit scale.
  expect_equal(stats::sd(built$matrix[, 1]), 1, tolerance = 1e-8)
})

test_that("constraint columns standardize to the sampling-weighted scale", {
  # Numeric constraint columns are centered and scaled by the sampling-weighted
  # mean and standard deviation, matching the scale the reference implementations
  # measure standardized mean differences on. The reliability-weighted variance
  # denominator sum(w) - sum(w^2) / sum(w) reduces to n - 1 for equal weights.
  x <- c(-2, -1, 0, 1, 2, 8)
  data <- data.frame(x1 = x)
  weights <- c(5, 5, 5, 1, 1, 1)

  sw <- sum(weights)
  weighted_mean <- sum(weights * x) / sw
  denom <- sw - sum(weights^2) / sw
  weighted_sd <- sqrt(sum(weights * (x - weighted_mean)^2) / denom)

  built <- build_constraint_matrix(
    data,
    "x1",
    balance_terms(),
    exposure_type = "binary",
    sampling_weights = weights
  )
  record <- built$recipe[[1]]

  expect_equal(record$base_center, weighted_mean)
  expect_equal(record$scale, weighted_sd)
  # The weighted scale differs from the unweighted one when the weights are not
  # equal, so this is a genuine change of convention.
  expect_false(isTRUE(all.equal(record$scale, stats::sd(x))))
  # The standardized column has weighted mean zero and unit weighted variance.
  column <- built$matrix[, 1]
  expect_equal(sum(weights * column) / sw, 0, tolerance = 1e-12)
  expect_equal(sum(weights * column^2) / denom, 1, tolerance = 1e-12)

  # Without sampling weights the standardization is the unweighted sample scale.
  plain <- build_constraint_matrix(
    data,
    "x1",
    balance_terms(),
    exposure_type = "binary"
  )
  expect_equal(plain$recipe[[1]]$scale, stats::sd(x))
})

# ---- build_constraint_matrix(): interactions ------------------------------

test_that("interactions add pairwise products but skip within-factor pairs", {
  data <- data.frame(
    x1 = c(-1, 0, 1, 2, 0.5, -0.5),
    x2 = c(2, 1, 0, -1, 0.25, 0.75),
    f = factor(c("a", "b", "c", "a", "b", "c"))
  )
  built <- build_constraint_matrix(
    data,
    c("x1", "x2", "f"),
    balance_terms(interactions = TRUE),
    exposure_type = "binary"
  )

  kinds <- vapply(built$recipe, function(term) term$kind, character(1))
  interactions <- built$recipe[kinds == "interaction"]
  partners <- vapply(
    interactions,
    function(term) paste(sort(c(term$source, term$partner)), collapse = ":"),
    character(1)
  )
  # x1:x2, x1:f, x2:f cross-terms are present; no indicator of f is crossed
  # with another indicator of f.
  expect_true("x1:x2" %in% partners)
  expect_false(any(vapply(
    interactions,
    function(term) identical(term$source, "f") && identical(term$partner, "f"),
    logical(1)
  )))
})

# ---- build_constraint_matrix(): quantiles ---------------------------------

test_that("quantiles add indicator columns for a discrete exposure", {
  data <- data.frame(x1 = seq(-2, 2, length.out = 20))
  built <- build_constraint_matrix(
    data,
    "x1",
    balance_terms(quantiles = c(0.25, 0.5, 0.75)),
    exposure_type = "binary"
  )

  kinds <- vapply(built$recipe, function(term) term$kind, character(1))
  expect_identical(sum(kinds == "quantile"), 3L)
})

test_that("quantiles with a continuous exposure error", {
  data <- data.frame(x1 = seq(-2, 2, length.out = 20))
  expect_error(
    build_constraint_matrix(
      data,
      "x1",
      balance_terms(quantiles = 0.5),
      exposure_type = "continuous"
    ),
    class = "balancing_constraints_error"
  )
})

# ---- build_constraint_matrix(): aliased columns ---------------------------

test_that("aliased columns are dropped with an alert", {
  withr::local_options(balancing.quiet = FALSE)
  data <- data.frame(
    x1 = c(-1, 0, 1, 2, 0.5, -0.5)
  )
  data$x2 <- data$x1
  expect_message(
    built <- build_constraint_matrix(
      data,
      c("x1", "x2"),
      balance_terms(),
      exposure_type = "binary"
    )
  )
  expect_identical(ncol(built$matrix), 1L)
})

test_that("a shift-related covariate drops like an exact duplicate", {
  withr::local_options(balancing.quiet = FALSE)
  data <- data.frame(x1 = c(-1, 0, 1, 2, 0.5, -0.5))
  data$duplicate <- data$x1
  data$shifted <- data$x1 + 5

  expect_message(
    duplicated_pair <- build_constraint_matrix(
      data,
      c("x1", "duplicate"),
      balance_terms(),
      exposure_type = "binary"
    )
  )
  expect_message(
    shifted_pair <- build_constraint_matrix(
      data,
      c("x1", "shifted"),
      balance_terms(),
      exposure_type = "binary"
    )
  )

  # Centering removes the shift, so the two covariates contribute the same
  # constraint column and the shifted pair must drop exactly as the exact
  # duplicate does.
  expect_identical(ncol(duplicated_pair$matrix), 1L)
  expect_identical(ncol(shifted_pair$matrix), 1L)
  expect_equal(shifted_pair$matrix[, 1], duplicated_pair$matrix[, 1])
})

test_that("an affine-related covariate drops", {
  withr::local_options(balancing.quiet = FALSE)
  data <- data.frame(x1 = c(-1, 0, 1, 2, 0.5, -0.5))
  data$reversed <- 1 - data$x1

  expect_message(
    built <- build_constraint_matrix(
      data,
      c("x1", "reversed"),
      balance_terms(),
      exposure_type = "binary"
    )
  )
  expect_identical(ncol(built$matrix), 1L)
})

test_that("indicator, quantile, and factor columns keep their full set", {
  withr::local_seed(404)
  n <- 60
  data <- data.frame(
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n),
    b = rep(c(0, 1), length.out = n),
    f = factor(rep(c("a", "b", "c"), length.out = n))
  )
  built <- build_constraint_matrix(
    data,
    c("x1", "x2", "b", "f"),
    balance_terms(moments = 2L, quantiles = c(0.25, 0.75)),
    exposure_type = "binary"
  )
  terms <- vapply(built$recipe, function(record) record$term, character(1))

  # A zero/one indicator, a quantile indicator, and a full set of factor level
  # indicators are all linearly independent of one another, so nothing here is
  # aliased and every column survives.
  expect_identical(
    terms,
    c(
      "x1",
      "x1^2",
      "x2",
      "x2^2",
      "b",
      "f_a",
      "f_b",
      "f_c",
      "x1_q0.25",
      "x1_q0.75",
      "x2_q0.25",
      "x2_q0.75"
    )
  )
})

test_that("interactions with a factor keep their aliased-column drops", {
  withr::local_options(balancing.quiet = FALSE)
  withr::local_seed(505)
  n <- 60
  data <- data.frame(
    x1 = stats::rnorm(n),
    b = rep(c(0, 1), length.out = n),
    f = factor(rep(c("a", "b", "c"), length.out = n))
  )
  expect_message(
    built <- build_constraint_matrix(
      data,
      c("x1", "b", "f"),
      balance_terms(interactions = TRUE),
      exposure_type = "binary"
    )
  )
  terms <- vapply(built$recipe, function(record) record$term, character(1))

  # Each covariate crossed with the full set of factor indicators sums back to
  # the covariate itself, so one product per covariate is redundant.
  expect_identical(
    terms,
    c(
      "x1",
      "b",
      "f_a",
      "f_b",
      "f_c",
      "x1:b",
      "x1:f_a",
      "x1:f_b",
      "b:f_a",
      "b:f_b"
    )
  )
})

# ---- build_constraint_matrix(): constant columns --------------------------

test_that("a constant covariate is dropped with a constant-column alert", {
  withr::local_options(balancing.quiet = FALSE)
  data <- data.frame(x1 = c(-1, 0, 1, 2, 0.5, -0.5), fixed = 5)

  expect_message(
    built <- build_constraint_matrix(
      data,
      c("x1", "fixed"),
      balance_terms(),
      exposure_type = "continuous"
    ),
    "constant"
  )

  # A covariate with no spread carries no balance information, and the
  # statistics reported on it divide by its zero spread. It is reported as
  # constant rather than aliased: nothing else in the set stands in for it.
  terms <- vapply(built$recipe, function(record) record$term, character(1))
  expect_identical(terms, "x1")
  expect_identical(ncol(built$matrix), 1L)
})

test_that("a single-level factor is dropped with a constant-column alert", {
  withr::local_options(balancing.quiet = FALSE)
  data <- data.frame(
    x1 = c(-1, 0, 1, 2, 0.5, -0.5),
    f = factor(rep("a", 6))
  )

  # The single level's indicator is one everywhere, which is constant without
  # being zero, so the rank check leaves it in place.
  expect_message(
    built <- build_constraint_matrix(
      data,
      c("x1", "f"),
      balance_terms(),
      exposure_type = "continuous"
    ),
    "constant"
  )

  terms <- vapply(built$recipe, function(record) record$term, character(1))
  expect_identical(terms, "x1")
})

test_that("a constant interaction is dropped alongside its constant bases", {
  withr::local_options(balancing.quiet = FALSE)
  data <- data.frame(x1 = c(-1, 0, 1, 2, 0.5, -0.5), fixed = 5, held = 3)

  # The product of two constants is constant, so the interaction branch is
  # covered by the same detection; the products of the varying covariate with
  # each constant are proportional to it and remain aliased drops.
  expect_message(
    expect_message(
      built <- build_constraint_matrix(
        data,
        c("x1", "fixed", "held"),
        balance_terms(interactions = TRUE),
        exposure_type = "continuous"
      ),
      "constant"
    ),
    "aliased"
  )

  terms <- vapply(built$recipe, function(record) record$term, character(1))
  expect_identical(terms, "x1")
})

test_that("a wholly constant constraint set is a classed error", {
  data <- data.frame(fixed = rep(5, 6), held = rep(3, 6))
  expect_error(
    build_constraint_matrix(
      data,
      c("fixed", "held"),
      balance_terms(),
      exposure_type = "continuous"
    ),
    class = "balancing_constraints_error"
  )
})

# ---- rebuild_constraint_matrix(): round trip ------------------------------

test_that("rebuild_constraint_matrix() reproduces the built matrix", {
  data <- data.frame(
    x1 = c(-2, -1, 0, 1, 2, 3),
    f = factor(c("a", "b", "c", "a", "b", "c"))
  )
  built <- build_constraint_matrix(
    data,
    c("x1", "f"),
    balance_terms(moments = 2L),
    exposure_type = "binary"
  )
  rebuilt <- rebuild_constraint_matrix(built$recipe, data)

  expect_equal(rebuilt, built$matrix)
})

# ---- Per-covariate moments ------------------------------------------------

test_that("a scalar moments value applies to every covariate", {
  data <- withr::with_seed(31, data.frame(x1 = rnorm(20), x2 = rnorm(20)))
  built <- build_constraint_matrix(
    data,
    c("x1", "x2"),
    balance_terms(moments = 2L),
    exposure_type = "binary"
  )
  expect_identical(record_terms(built$recipe), c("x1", "x1^2", "x2", "x2^2"))
})

test_that("an unnamed multi-element moments vector is a classed error", {
  # An unnamed vector of length two recycled its first element to every
  # covariate and dropped the rest, which is the ambiguity to name rather than
  # resolve silently.
  data <- withr::with_seed(31, data.frame(x1 = rnorm(20), x2 = rnorm(20)))
  expect_error(
    build_constraint_matrix(
      data,
      c("x1", "x2"),
      balance_terms(moments = c(2L, 3L)),
      exposure_type = "binary"
    ),
    class = "balancing_constraints_error"
  )
})

test_that("moments named for a non-covariate is a classed error", {
  # A misspelled name previously left every covariate at first moments, so the
  # requested power never reached the fit.
  data <- withr::with_seed(31, data.frame(x1 = rnorm(20), x2 = rnorm(20)))
  expect_error(
    build_constraint_matrix(
      data,
      c("x1", "x2"),
      balance_terms(moments = c(x22 = 2L)),
      exposure_type = "binary"
    ),
    class = "balancing_constraints_error"
  )
})

# ---- Per-covariate tolerance ----------------------------------------------

test_that("a scalar tolerance applies to every constraint column", {
  data <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  built <- build_constraint_matrix(
    data,
    c("x1", "x2"),
    balance_terms(tolerance = 0.1),
    exposure_type = "binary"
  )
  tolerances <- vapply(built$recipe, function(r) r$tolerance, numeric(1))
  expect_equal(tolerances, c(0.1, 0.1))
})

test_that("a named tolerance sets per-covariate values and derived columns inherit", {
  data <- data.frame(x1 = rnorm(30), x2 = rnorm(30))
  built <- build_constraint_matrix(
    data,
    c("x1", "x2"),
    balance_terms(moments = c(x1 = 2L), tolerance = c(x1 = 0.2)),
    exposure_type = "binary"
  )
  terms <- vapply(built$recipe, function(r) r$term, character(1))
  tolerances <- stats::setNames(
    vapply(built$recipe, function(r) r$tolerance, numeric(1)),
    terms
  )
  # The x1 power column inherits x1's tolerance; the unnamed x2 defaults to exact.
  expect_equal(tolerances[["x1"]], 0.2)
  expect_equal(tolerances[["x1^2"]], 0.2)
  expect_equal(tolerances[["x2"]], 0)
})

test_that("an unnamed multi-element tolerance is a classed error", {
  data <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  expect_error(
    build_constraint_matrix(
      data,
      c("x1", "x2"),
      balance_terms(tolerance = c(0.1, 0.2)),
      exposure_type = "binary"
    ),
    class = "balancing_constraints_error"
  )
})

test_that("a tolerance named for a non-covariate is a classed error", {
  data <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  expect_error(
    build_constraint_matrix(
      data,
      c("x1", "x2"),
      balance_terms(tolerance = c(nonesuch = 0.1)),
      exposure_type = "binary"
    ),
    class = "balancing_constraints_error"
  )
})

# ---- List quantiles -------------------------------------------------------

test_that("a quantile list adds per-covariate quantile columns", {
  data <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  built <- build_constraint_matrix(
    data,
    c("x1", "x2"),
    balance_terms(quantiles = list(x1 = c(0.25, 0.75), x2 = 0.5)),
    exposure_type = "binary"
  )
  terms <- vapply(built$recipe, function(r) r$term, character(1))
  quantile_terms <- terms[grepl("_q", terms)]
  expect_setequal(quantile_terms, c("x1_q0.25", "x1_q0.75", "x2_q0.5"))
})

test_that("a covariate absent from the quantile list contributes no columns", {
  data <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  built <- build_constraint_matrix(
    data,
    c("x1", "x2"),
    balance_terms(quantiles = list(x1 = 0.5)),
    exposure_type = "binary"
  )
  terms <- vapply(built$recipe, function(r) r$term, character(1))
  expect_true("x1_q0.5" %in% terms)
  expect_false(any(grepl("x2_q", terms)))
})
