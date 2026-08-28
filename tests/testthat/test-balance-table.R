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
