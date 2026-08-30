# compute_balance_table() assembles the balance table from the constraint
# matrix and the extracted weights. balance() has already built that matrix
# before it fits, so the table can be handed the built matrix rather than
# rebuilding it from the recipe. These specs pin that the two routes report the
# same table on every exposure type, estimand and weighting the fit supports,
# and pin the column statistics the table is assembled from against references
# written out in plain arithmetic, so a change of algorithm inside them has
# something independent to answer to.

# ---- Prebuilt constraint matrix --------------------------------------------

# Reassemble the arguments balance() passes to compute_balance_table(), so the
# equivalence specs run on the geometry a real fit reports on rather than on a
# plausible-looking substitute. The pieces mirror balance(): the constraint
# matrix is built with the sampling weights for the quadratic-program family and
# without them otherwise, the groups are the row indices of each exposure level,
# the weights are the fitted weights with the sampling weights composed on, and
# the reference measure is the product of the sampling weights and the method's
# base weights.
balance_table_pieces <- function(
  fit,
  data,
  covariates,
  constraint_target = "pooled"
) {
  constraint_sampling_weights <- if (
    S7::S7_inherits(fit@method, quadratic_program_method)
  ) {
    fit@sampling_weights
  } else {
    NULL
  }
  built <- build_constraint_matrix(
    data,
    covariates,
    fit@constraints,
    fit@exposure_type,
    sampling_weights = constraint_sampling_weights
  )

  exposure_vec <- data[[fit@exposure]]
  groups <- if (identical(fit@exposure_type, "continuous")) {
    NULL
  } else {
    exposure_key <- as.character(exposure_vec)
    stats::setNames(
      lapply(fit@exposure_levels, function(level) which(exposure_key == level)),
      fit@exposure_levels
    )
  }

  sampling_weights <- fit@sampling_weights %||% rep(1, fit@n)
  base_weights <- if ("base_weights" %in% S7::prop_names(fit@method)) {
    fit@method@base_weights %||% rep(1, fit@n)
  } else {
    rep(1, fit@n)
  }

  args <- list(
    recipe = built$recipe,
    data = data,
    exposure_vec = exposure_vec,
    exposure_type = fit@exposure_type,
    estimand = fit@estimand,
    focal_level = fit@focal_level,
    groups = groups,
    weights = as.numeric(stats::weights(fit)),
    tolerance = 0,
    reference = sampling_weights * base_weights,
    constraint_target = constraint_target,
    sampling_weights = fit@sampling_weights
  )

  list(args = args, matrix = built$matrix)
}

# The rebuild route is asserted against the fit's own table first, so the
# reassembled arguments are held to reproducing what balance() reported rather
# than only to agreeing with themselves.
expect_prebuilt_matrix_equivalence <- function(
  fit,
  data,
  covariates,
  constraint_target = "pooled"
) {
  pieces <- balance_table_pieces(fit, data, covariates, constraint_target)
  rebuilt <- do.call(compute_balance_table, pieces$args)
  expect_identical(rebuilt, fit@balance_table)

  prebuilt <- do.call(
    compute_balance_table,
    c(pieces$args, list(matrix = pieces$matrix))
  )
  expect_identical(prebuilt, rebuilt)
}

test_that("a prebuilt matrix reports the binary average-treatment-effect table", {
  data <- sim_binary(n = 200)
  fit <- balance(data, exposure, c(x1, x2, x3), method = bw_entropy())

  expect_prebuilt_matrix_equivalence(fit, data, c("x1", "x2", "x3"))
})

test_that("a prebuilt matrix reports the binary focal table", {
  # A focal estimand holds one arm fixed, so the table's target level and its
  # constraint geometry both differ from the average-treatment-effect case.
  data <- sim_binary(n = 200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2, x3),
    method = bw_entropy(),
    estimand = "att"
  )

  expect_prebuilt_matrix_equivalence(fit, data, c("x1", "x2", "x3"))
})

test_that("a prebuilt matrix reports the multi-group categorical table", {
  # Three exposure levels put two arms in the pairwise maximum, so the group
  # label the table reports comes from an argmax rather than a single contrast.
  data <- sim_categorical(n = 200)
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())

  expect_prebuilt_matrix_equivalence(fit, data, c("x1", "x2"))
})

test_that("a prebuilt matrix reports the continuous correlation table", {
  data <- sim_continuous(n = 200)
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())

  expect_prebuilt_matrix_equivalence(fit, data, c("x1", "x2"))
})

test_that("a prebuilt matrix reports the sampling-weighted table", {
  # Sampling weights move both the standardization the table reports on and the
  # reference measure the pooled constraint target is taken against.
  data <- sim_binary(n = 200)
  data$sw <- withr::with_seed(11, stats::runif(nrow(data), 0.5, 2))
  fit <- balance(
    data,
    exposure,
    c(x1, x2, x3),
    method = bw_entropy(),
    sampling_weights = sw
  )

  expect_prebuilt_matrix_equivalence(fit, data, c("x1", "x2", "x3"))
})

test_that("a prebuilt matrix reports the per-term tolerance table", {
  # Per-term tolerances reach the table through the recipe rather than the
  # `tolerance` argument, and the quadratic-program family builds its constraint
  # matrix on the sampling-weighted scale, so this case covers both.
  data <- sim_binary(n = 200)
  data$sw <- withr::with_seed(12, stats::runif(nrow(data), 0.5, 2))
  fit <- balance(
    data,
    exposure,
    c(x1, x2, x3),
    method = bw_sbw(),
    constraints = balance_terms(
      tolerance = c(x1 = 0.02, x2 = 0.05, x3 = 0.1)
    ),
    sampling_weights = sw
  )

  expect_prebuilt_matrix_equivalence(fit, data, c("x1", "x2", "x3"))
})

test_that("a prebuilt matrix reports the arm-to-arm constraint table", {
  # The covariate balancing propensity score equates the arms directly, so its
  # fit reports an arm-to-arm constraint target rather than the pooled one.
  data <- sim_binary(n = 200)
  fit <- balance(data, exposure, c(x1, x2, x3), method = bw_cbps())

  expect_prebuilt_matrix_equivalence(
    fit,
    data,
    c("x1", "x2", "x3"),
    constraint_target = "arms"
  )
})

test_that("a prebuilt matrix reports the base-weighted pooled target table", {
  # Entropy balancing is the only method carrying base weights, and they enter
  # the table as the reference measure the average-treatment-effect pooled
  # target is taken against. Non-uniform base weights move that target away from
  # the plain column mean, so this pins the branch a uniform measure hides.
  data <- sim_binary(n = 200)
  base_weights <- withr::with_seed(13, stats::runif(nrow(data), 0.5, 2))
  fit <- balance(
    data,
    exposure,
    c(x1, x2, x3),
    method = bw_entropy(base_weights = base_weights)
  )

  expect_prebuilt_matrix_equivalence(fit, data, c("x1", "x2", "x3"))
})

# ---- Column statistics -----------------------------------------------------

# A fixture with the features the column statistics have to survive: a column
# with no spread, which the zero-scale guard leaves unscaled; an indicator, whose
# scale is small next to the continuous columns; sampling weights that are not
# uniform, so the weighted and unweighted statistics differ; and a group
# membership that splits the rows unevenly, so a subset mean is not the whole
# column's mean.
#
# The column with no spread holds a power of two so that its weighted mean comes
# out exactly, in any summation order: scaling a weight by two is exact, so the
# weighted sum is exactly twice the sum of the weights and the centered column is
# exactly zero. A constant that is not a power of two leaves rounding noise in
# the centered column, the weighted scale picks that noise up instead of
# reporting zero, and the guard is never reached.
column_statistic_fixture <- function() {
  withr::with_seed(303, {
    n <- 40
    m <- cbind(
      x1 = stats::rnorm(n),
      x2 = stats::runif(n, -2, 3),
      constant = rep(2, n),
      indicator = as.numeric(stats::rbinom(n, 1L, 0.4))
    )
    list(
      m = m,
      sampling_weights = stats::runif(n, 0.3, 2.5),
      weights = stats::runif(n, 0.1, 4),
      idx = which(stats::rbinom(n, 1L, 0.35) == 1L)
    )
  })
}

test_that("standardize_columns() centers and scales on the unweighted scale", {
  fixture <- column_statistic_fixture()
  m <- fixture$m

  expected <- m
  for (j in seq_len(ncol(m))) {
    column <- m[, j]
    center <- sum(column) / length(column)
    variance <- sum((column - center)^2) / (length(column) - 1)
    scale <- sqrt(variance)
    if (scale == 0) {
      scale <- 1
    }
    expected[, j] <- (column - center) / scale
  }

  expect_equal(standardize_columns(m), expected, tolerance = 1e-12)
})

test_that("standardize_columns() centers and scales on the sampling-weighted scale", {
  fixture <- column_statistic_fixture()
  m <- fixture$m
  w <- fixture$sampling_weights

  expected <- m
  for (j in seq_len(ncol(m))) {
    column <- m[, j]
    center <- sum(w * column) / sum(w)
    denominator <- sum(w) - sum(w * w) / sum(w)
    variance <- sum(w * (column - center)^2) / denominator
    scale <- sqrt(variance)
    if (scale == 0) {
      scale <- 1
    }
    expected[, j] <- (column - center) / scale
  }

  expect_equal(
    standardize_columns(m, sampling_weights = w),
    expected,
    tolerance = 1e-12
  )
})

test_that("standardize_columns() reproduces the per-column weighted statistics", {
  # The vectorized weighted branch is bit for bit the per-column formulation it
  # replaced, so this pins identity rather than agreement to a tolerance. A
  # tolerance here would let a genuine change of accumulation order pass.
  fixture <- column_statistic_fixture()
  m <- fixture$m
  w <- fixture$sampling_weights

  centers <- apply(m, 2, weighted_center, w = w)
  scales <- apply(m, 2, weighted_scale, w = w)
  scales[scales == 0] <- 1
  expected <- sweep(sweep(m, 2, centers, "-"), 2, scales, "/")

  expect_identical(standardize_columns(m, sampling_weights = w), expected)
})

test_that("standardize_columns() agrees with stats::sd() to rounding", {
  # The unweighted scale carries the same corrected two-pass center `stats::sd()`
  # does, but `stats::sd()` carries the correction in long double, so the two
  # agree to a unit in the last place rather than exactly.
  fixture <- column_statistic_fixture()
  m <- fixture$m

  scales <- apply(m, 2, stats::sd)
  scales[scales == 0] <- 1
  expected <- sweep(sweep(m, 2, colMeans(m), "-"), 2, scales, "/")

  expect_equal(standardize_columns(m), expected, tolerance = 1e-14)

  # An offset of 1e14 on the non-constant columns is large enough against their
  # spread to separate the corrected two-pass center from the plain one: on this
  # input the plain form departs from `stats::sd()` by a relative 1.2e-4, well
  # outside the tolerance below, while the corrected form still tracks it.
  shifted <- m
  offset <- setdiff(colnames(m), "constant")
  shifted[, offset] <- shifted[, offset] + 1e14

  shifted_scales <- apply(shifted, 2, stats::sd)
  shifted_scales[shifted_scales == 0] <- 1
  shifted_expected <- sweep(
    sweep(shifted, 2, colMeans(shifted), "-"),
    2,
    shifted_scales,
    "/"
  )

  expect_equal(
    standardize_columns(shifted),
    shifted_expected,
    tolerance = 1e-14
  )
})

test_that("standardize_columns() leaves a column with no spread unscaled", {
  # The zero-scale guard divides by one rather than by zero, so the constant
  # column comes out as zeros instead of as missing values.
  fixture <- column_statistic_fixture()
  m <- fixture$m

  unweighted <- standardize_columns(m)
  weighted <- standardize_columns(
    m,
    sampling_weights = fixture$sampling_weights
  )

  expect_equal(unweighted[, "constant"], rep(0, nrow(m)), tolerance = 1e-12)
  expect_equal(weighted[, "constant"], rep(0, nrow(m)), tolerance = 1e-12)
})

# A constant column whose value is not exactly representable in the weighted
# center's arithmetic does not center to zero. `sum(w * x) / sum(w)` rounds
# twice, so the centered column holds a rounding residual rather than a zero,
# the weighted scale reports that residual as the column's spread, and dividing
# by it turns a column with nothing in it into order-one values. The column has
# no spread whatever the arithmetic says, so it standardizes to zero.
test_that("standardize_columns() flattens a constant column under sampling weights", {
  n <- 20L
  w <- withr::with_seed(1, stats::runif(n, 0.5, 2))
  m <- matrix(0.98, nrow = n, ncol = 1L, dimnames = list(NULL, "constant"))

  # The fixture is only worth having while the arithmetic still misses: the
  # weighted center of this column is a rounding step away from its value.
  expect_false(identical(sum(w * m[, 1L]) / sum(w), 0.98))

  expect_identical(
    standardize_columns(m, sampling_weights = w),
    matrix(0, nrow = n, ncol = 1L, dimnames = list(NULL, "constant"))
  )
  expect_identical(
    standardize_columns(m),
    matrix(0, nrow = n, ncol = 1L, dimnames = list(NULL, "constant"))
  )
})

# The other side of the same rule. A column that barely varies still varies, and
# a spread of 1e-6 around 0.98 is a real one rather than a rounding residual, so
# it is standardized rather than flattened.
test_that("standardize_columns() leaves a nearly constant column alone", {
  n <- 20L
  w <- withr::with_seed(1, stats::runif(n, 0.5, 2))
  m <- matrix(
    0.98 + 1e-6 * seq_len(n),
    nrow = n,
    ncol = 1L,
    dimnames = list(NULL, "nearly")
  )

  weighted <- standardize_columns(m, sampling_weights = w)
  unweighted <- standardize_columns(m)

  expect_equal(stats::sd(unweighted[, 1L]), 1, tolerance = 1e-8)
  expect_gt(diff(range(weighted[, 1L])), 1)
  expect_gt(diff(range(unweighted[, 1L])), 1)
})

# `column_is_constant()` reads the values rather than the computed scale, and it
# reads them with a comparison that a missing value answers with a missing
# value. A column carrying one would therefore make the reading itself missing,
# and the two places that consume it index with it: `centered[, constant] <- 0`
# and `scales[constant | scales == 0] <- 1` both refuse a missing subscript in an
# assignment, so a single missing value would stop the standardization with a
# base error rather than a classed one.
#
# A column with a missing value is read as not constant. That is the reading
# that leaves the rest of the standardization behaving exactly as it did: the
# column keeps its computed center and scale, and the missing value propagates
# through them into the standardized column, where it is visible for what it is.
# The constraint columns a fit builds carry no missing values, which the fit
# validates before a constraint matrix exists, so this governs only the direct
# callers.
test_that("column_is_constant() reads a column with a missing value as varying", {
  m <- cbind(
    constant = rep(0.98, 4),
    missing_constant = c(0.98, NA, 0.98, 0.98),
    all_missing = rep(NA_real_, 4),
    varying = c(1, 2, 3, 4)
  )

  expect_identical(
    column_is_constant(m),
    c(
      constant = TRUE,
      missing_constant = FALSE,
      all_missing = FALSE,
      varying = FALSE
    )
  )
})

test_that("column_is_constant() reads every column of a single row as constant", {
  m <- matrix(
    c(0.98, NA, 3),
    nrow = 1L,
    dimnames = list(NULL, c("value", "missing", "other"))
  )

  expect_identical(
    column_is_constant(m),
    c(value = TRUE, missing = FALSE, other = TRUE)
  )
})

test_that("standardize_columns() carries a missing value through the column", {
  n <- 20L
  w <- withr::with_seed(1, stats::runif(n, 0.5, 2))
  values <- rep(0.98, n)
  values[[3L]] <- NA_real_
  m <- matrix(values, nrow = n, ncol = 1L, dimnames = list(NULL, "constant"))

  weighted <- standardize_columns(m, sampling_weights = w)
  unweighted <- standardize_columns(m)

  expect_identical(dim(weighted), c(n, 1L))
  expect_identical(dim(unweighted), c(n, 1L))
  expect_true(is.na(weighted[3L, 1L]))
  expect_true(is.na(unweighted[3L, 1L]))
})

# The column arithmetic the weighted standardization and the solver's tolerance
# box both read their scales from. It takes an already-centered matrix because
# both callers have one in hand, and it reproduces the per-column
# `weighted_scale()` bit for bit rather than approximately.
test_that("centered_column_scales() reproduces the per-column weighted scale", {
  fixture <- column_statistic_fixture()
  m <- fixture$m
  w <- fixture$sampling_weights
  centered <- sweep(m, 2, apply(m, 2, weighted_center, w = w), "-")

  expect_identical(
    centered_column_scales(centered, w),
    apply(m, 2, weighted_scale, w = w)
  )
})

test_that("weighted_column_means() averages each column over the row subset", {
  fixture <- column_statistic_fixture()
  z <- standardize_columns(fixture$m, fixture$sampling_weights)
  w <- fixture$weights
  idx <- fixture$idx

  expected <- vapply(
    colnames(z),
    function(term) {
      column <- z[idx, term]
      sum(w[idx] * column) / sum(w[idx])
    },
    numeric(1)
  )

  expect_equal(
    weighted_column_means(z, idx, w),
    expected,
    tolerance = 1e-12
  )
})

test_that("weighted_column_means() averages over every row under unit weights", {
  # The table takes its unweighted arm means through the same helper with a unit
  # weight vector, so that path is pinned as well.
  fixture <- column_statistic_fixture()
  z <- standardize_columns(fixture$m)
  idx <- seq_len(nrow(z))
  w <- rep(1, nrow(z))

  expected <- vapply(
    colnames(z),
    function(term) sum(z[, term]) / nrow(z),
    numeric(1)
  )

  expect_equal(weighted_column_means(z, idx, w), expected, tolerance = 1e-12)
})
