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

test_that("factor covariates expand to a level indicator set less the alias", {
  withr::local_options(balancing.quiet = FALSE)
  data <- data.frame(
    x1 = c(-1, 0, 1, 2, 0.5, -0.5),
    f = factor(c("a", "b", "c", "a", "b", "c"))
  )
  expect_message(
    built <- build_constraint_matrix(
      data,
      c("x1", "f"),
      balance_terms(),
      exposure_type = "binary"
    ),
    "aliased"
  )

  sources <- vapply(built$recipe, function(term) term$source, character(1))
  # The expansion is not reference-coded: every level gets an indicator. The
  # rank check against the constant the solvers carry then removes the one the
  # others determine, which is the last level rather than the first that
  # treatment contrasts drop.
  expect_identical(sum(sources == "f"), 2L)
  expect_identical(ncol(built$matrix), 3L)
  expect_identical(record_terms(built$recipe), c("x1", "f_a", "f_b"))
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

# Every solver carries an intercept: the entropy dual normalizes within each
# exposure group, IPT and CBPS bind an explicit column of ones, and the SBW,
# energy, and CFD programs carry a group-sum row. A constraint set that is
# affinely dependent on that constant is rank deficient in the geometry the
# solver sees even when the columns on their own are independent, so the rank
# check counts the constant alongside them. The shapes below are the ones where
# that matters: a full set of factor level indicators, a zero/one column with
# its complement, and a covariate that repeats a level indicator.

# Two continuous covariates, a two-level factor, and a binary exposure: the
# smallest shape whose level indicators sum to the constant.
two_level_factor_data <- function(n = 500) {
  withr::with_seed(
    1,
    data.frame(
      x1 = stats::rnorm(n),
      x2 = stats::rnorm(n),
      f = factor(sample(c("A", "B"), n, replace = TRUE)),
      exposure = stats::rbinom(n, 1L, 0.5)
    )
  )
}

# One continuous covariate and three factors of differing widths, so that the
# per-factor drop can be counted separately for each.
several_factor_data <- function(n = 300) {
  withr::with_seed(
    7,
    data.frame(
      x1 = stats::rnorm(n),
      f1 = factor(sample(c("A", "B"), n, replace = TRUE)),
      f2 = factor(sample(c("p", "q", "r"), n, replace = TRUE)),
      f3 = factor(sample(c("no", "yes"), n, replace = TRUE))
    )
  )
}

test_that("a two-level factor drops one level indicator", {
  withr::local_options(balancing.quiet = FALSE)
  data <- two_level_factor_data()

  expect_message(
    built <- build_constraint_matrix(
      data,
      c("x1", "x2", "f"),
      balance_terms(),
      exposure_type = "binary"
    ),
    "aliased"
  )

  # `f_A + f_B` is the constant column, so the second level constrains nothing
  # the first does not already constrain once the intercept is counted. The
  # constant comes first in the decomposition and LINPACK's dqrdc2 keeps column
  # order, moving only the deficient columns to the end, so the level dropped is
  # the last one.
  expect_identical(record_terms(built$recipe), c("x1", "x2", "f_A"))
  expect_equal(rebuild_constraint_matrix(built$recipe, data), built$matrix)
})

test_that("the built matrix is full rank against the intercept", {
  data <- two_level_factor_data()
  built <- build_constraint_matrix(
    data,
    c("x1", "x2", "f"),
    balance_terms(),
    exposure_type = "binary"
  )

  # The rank a solver sees is the rank of the constraint columns together with
  # the constant it carries, so that is the matrix the smallest singular value
  # has to be measured on. The bound sits far below the value this column set
  # reaches when it is full rank and far above the rounding-scale value a
  # deficient set leaves.
  expect_gt(min(svd(cbind(1, built$matrix))$d), 1e-6)
})

test_that("several factors each lose exactly one level indicator", {
  withr::local_options(balancing.quiet = FALSE)
  data <- several_factor_data()

  expect_message(
    built <- build_constraint_matrix(
      data,
      c("x1", "f1", "f2", "f3"),
      balance_terms(),
      exposure_type = "binary"
    ),
    "aliased"
  )

  # Each factor is affinely dependent on the constant on its own, so each loses
  # its last level and no more: the two-level factors keep one indicator each
  # and the three-level factor keeps two.
  expect_identical(
    record_terms(built$recipe),
    c("x1", "f1_A", "f2_p", "f2_q", "f3_no")
  )
  expect_gt(min(svd(cbind(1, built$matrix))$d), 1e-6)
})

test_that("the surviving level indicators do not depend on covariate order", {
  withr::local_options(balancing.quiet = FALSE)
  data <- several_factor_data()

  expect_message(
    built <- build_constraint_matrix(
      data,
      c("x1", "f1", "f2", "f3"),
      balance_terms(),
      exposure_type = "binary"
    ),
    "aliased"
  )
  expect_message(
    permuted <- build_constraint_matrix(
      data,
      c("f3", "f2", "x1", "f1"),
      balance_terms(),
      exposure_type = "binary"
    ),
    "aliased"
  )

  # Each factor's redundancy is with the constant alone, so which level survives
  # is settled inside the factor and no covariate borrows its drop from another.
  # The constrained set is therefore the same whatever order the covariates are
  # named in; only the column order changes.
  expect_setequal(record_terms(permuted$recipe), record_terms(built$recipe))
})

test_that("a zero/one covariate and its complement keep one column", {
  withr::local_options(balancing.quiet = FALSE)
  data <- data.frame(b = rep(c(0, 1), length.out = 20))
  data$b_complement <- 1 - data$b
  data$x1 <- withr::with_seed(3, stats::rnorm(20))

  expect_message(
    built <- build_constraint_matrix(
      data,
      c("x1", "b", "b_complement"),
      balance_terms(),
      exposure_type = "binary"
    ),
    "aliased"
  )

  # Indicator columns cross the boundary raw, so a zero/one column and its
  # complement are independent of one another and only the constant relates
  # them. Balancing either one balances the other.
  expect_identical(record_terms(built$recipe), c("x1", "b"))
})

test_that("a covariate repeating a factor level drops with the aliased level", {
  withr::local_options(balancing.quiet = FALSE)
  data <- data.frame(f = factor(rep(c("A", "B"), length.out = 20)))
  data$is_a <- as.numeric(data$f == "A")
  data$x1 <- withr::with_seed(4, stats::rnorm(20))

  expect_message(
    built <- build_constraint_matrix(
      data,
      c("x1", "f", "is_a"),
      balance_terms(),
      exposure_type = "binary"
    ),
    "aliased"
  )

  # Two redundancies are present at once: `f_B` repeats the constant less `f_A`,
  # and `is_a` repeats `f_A` outright. Both are later than the column that
  # stands in for them, so both are dropped and the single indicator that
  # carries the factor survives.
  expect_identical(record_terms(built$recipe), c("x1", "f_A"))
})

test_that("indicator and quantile columns keep their full set", {
  withr::local_options(balancing.quiet = FALSE)
  withr::local_seed(404)
  n <- 60
  data <- data.frame(
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n),
    b = rep(c(0, 1), length.out = n),
    f = factor(rep(c("a", "b", "c"), length.out = n))
  )
  expect_message(
    built <- build_constraint_matrix(
      data,
      c("x1", "x2", "b", "f"),
      balance_terms(moments = 2L, quantiles = c(0.25, 0.75)),
      exposure_type = "binary"
    ),
    "aliased"
  )
  terms <- vapply(built$recipe, function(record) record$term, character(1))

  # A zero/one indicator and a quantile indicator are each independent of the
  # constant the solver carries, so both keep every column they contribute. The
  # level indicators of a factor sum to that constant, so the factor loses one
  # level: the constant is first in the decomposition and LINPACK's dqrdc2 keeps
  # column order, moving only deficient columns to the end, so the last level of
  # each factor is the one dropped.
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

  # Three drops follow from the one rule, reading the columns left to right
  # after the constant. The factor's own indicators sum to the constant, so
  # `f_c` goes. Each covariate crossed with the full set of factor indicators
  # sums back to the covariate itself, and that covariate is already a column,
  # so the last product of each goes too: `x1:f_c` and `b:f_c`. What survives
  # spans the six cells the factor and the indicator cut the sample into,
  # through `f_a`, `f_b`, `b`, `b:f_a`, `b:f_b` and the constant, plus a slope
  # in `x1` for each of `1`, `f_a`, `f_b` and `b`.
  expect_identical(
    terms,
    c(
      "x1",
      "b",
      "f_a",
      "f_b",
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

# A name given twice used to take the first value and drop the rest, so the
# requested power silently never reached the fit.
test_that("a duplicated moments name is a classed error", {
  data <- withr::with_seed(31, data.frame(x1 = rnorm(20), x2 = rnorm(20)))
  expect_error(
    build_constraint_matrix(
      data,
      c("x1", "x2"),
      balance_terms(moments = c(x1 = 2L, x1 = 3L)),
      exposure_type = "binary"
    ),
    class = "balancing_constraints_error"
  )
})

# A partly named vector leaves its unnamed elements with no covariate to apply
# to. The empty name used to be reported as a covariate that does not exist,
# which named the wrong defect.
test_that("a partially named moments vector names the unnamed elements", {
  data <- withr::with_seed(31, data.frame(x1 = rnorm(20), x2 = rnorm(20)))
  expect_error(
    build_constraint_matrix(
      data,
      c("x1", "x2"),
      balance_terms(moments = c(x1 = 2L, 3L)),
      exposure_type = "binary"
    ),
    class = "balancing_constraints_error"
  )
  expect_error(
    build_constraint_matrix(
      data,
      c("x1", "x2"),
      balance_terms(moments = c(x1 = 2L, 3L)),
      exposure_type = "binary"
    ),
    "no name"
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

# The tolerance specification is resolved by the same rules as the moments, so a
# duplicated name and a partly named vector are refused there too.
test_that("a duplicated tolerance name is a classed error", {
  data <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  expect_error(
    build_constraint_matrix(
      data,
      c("x1", "x2"),
      balance_terms(tolerance = c(x1 = 0.1, x1 = 0.2)),
      exposure_type = "binary"
    ),
    class = "balancing_constraints_error"
  )
})

test_that("a partially named tolerance vector names the unnamed elements", {
  data <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  expect_error(
    build_constraint_matrix(
      data,
      c("x1", "x2"),
      balance_terms(tolerance = c(x1 = 0.1, 0.2)),
      exposure_type = "binary"
    ),
    "no name"
  )
})

# A name that is missing is a name the caller never supplied, so it belongs to
# the same refusal as an empty one. `nzchar()` reads a missing string as a name
# of some length, which walks a missing name past the unnamed-element check and
# leaves it to be reported as a covariate the data does not have. That names the
# wrong defect: nothing was misspelled, an element was left unnamed.
test_that("a missing tolerance name is reported as an element with no name", {
  data <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  tolerance <- c(0.1, 0.2)
  names(tolerance) <- c("x1", NA)
  expect_error(
    build_constraint_matrix(
      data,
      c("x1", "x2"),
      balance_terms(tolerance = tolerance),
      exposure_type = "binary"
    ),
    "no name",
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

# ---- An empty constraint set across the method families -------------------

# `moments = 0` asks for no moment columns, so a selection of continuous numeric
# covariates expands to a constraint matrix with no columns at all. The two
# families answer that differently by design. The estimating-equation family and
# stable balancing weights are defined by their constraints and have nothing left
# to solve, so they refuse the fit on the R side before the solver is called.
# Energy and characteristic function distance balancing are driven by their
# objective, which defines the solution on its own, so they fit as usual.
no_moment_constraints <- function() {
  balance_terms(moments = 0L)
}

test_that("no moments over continuous covariates leaves no constraint columns", {
  data <- sim_binary(n = 100)
  built <- build_constraint_matrix(
    data,
    c("x1", "x2"),
    no_moment_constraints(),
    exposure_type = "binary"
  )
  expect_identical(ncol(built$matrix), 0L)
  expect_length(built$recipe, 0L)
})

test_that("constraint-defined methods refuse an empty constraint set", {
  binary <- sim_binary(n = 150)
  categorical <- sim_categorical(n = 150)

  for (method in list(bw_entropy(), bw_ipt(), bw_cbps(), bw_sbw())) {
    expect_error(
      balance(
        binary,
        exposure,
        c(x1, x2),
        method = method,
        constraints = no_moment_constraints()
      ),
      class = "balancing_constraints_error"
    )
    expect_error(
      balance(
        categorical,
        exposure,
        c(x1, x2),
        method = method,
        constraints = no_moment_constraints()
      ),
      class = "balancing_constraints_error"
    )
  }
})

test_that("constraint-defined methods refuse an empty continuous constraint set", {
  continuous <- sim_continuous(n = 150)

  for (method in list(bw_entropy(), bw_cbps(), bw_sbw())) {
    expect_error(
      balance(
        continuous,
        exposure,
        c(x1, x2),
        method = method,
        constraints = no_moment_constraints()
      ),
      class = "balancing_constraints_error"
    )
  }
})

test_that("objective-driven methods fit with an empty constraint set", {
  binary <- sim_binary(n = 150)
  categorical <- sim_categorical(n = 150)

  for (method in list(bw_energy(), bw_cfd())) {
    fit <- balance(
      binary,
      exposure,
      c(x1, x2),
      method = method,
      constraints = no_moment_constraints()
    )
    expect_identical(nrow(fit@balance_table), 0L)
    expect_true(all(is.finite(as.numeric(stats::weights(fit)))))

    fit_categorical <- balance(
      categorical,
      exposure,
      c(x1, x2),
      method = method,
      constraints = no_moment_constraints()
    )
    expect_true(all(is.finite(as.numeric(stats::weights(fit_categorical)))))
  }

  fit_continuous <- balance(
    sim_continuous(n = 150),
    exposure,
    c(x1, x2),
    method = bw_energy(),
    constraints = no_moment_constraints()
  )
  expect_true(all(is.finite(as.numeric(stats::weights(fit_continuous)))))
})

test_that("a fit with no constraint terms prints and summarizes cleanly", {
  # An empty balance table has no largest imbalance and no tolerance to report,
  # which the maximum over an empty vector would render as -Inf behind a warning.
  data <- sim_binary(n = 150)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    constraints = no_moment_constraints()
  )
  printed <- utils::capture.output(expect_no_warning(print(fit)))
  expect_false(any(grepl("Inf", printed, fixed = TRUE)))
  expect_no_warning(utils::capture.output(summary(fit)))
})
