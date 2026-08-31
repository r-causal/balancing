# bw_energy() is the method spec; balance(..., method = bw_energy())
# fits it. Energy balancing (Huling and Mak) chooses weights that minimize the
# energy distance between the weighted exposure groups and the target sample; the
# improved variant adds the between-group energy term for the average treatment
# effect. For a continuous exposure the objective is the weighted distance
# covariance between the exposure and the covariates (Huling, Greifer, and Chen).
#
# Energy balancing belongs to the quadratic-program family, which carries a
# different set of promises from the estimating-equation family. Its central
# guarantee is that the achieved energy distance (or distance covariance)
# improves on the unweighted sample and that any added moment constraints are
# satisfied within tolerance; it does not promise exact moment balance when no
# moment constraints are requested. The family has no estimating equations, so
# estimating_equations() and ipw() raise the shared unsupported condition. These
# specs therefore verify the objective directly through independent energy-
# distance and distance-covariance oracles computed on the same distance
# definition the core uses, verify constraint satisfaction with the shared
# expect_balanced() oracle, and check the capability and diagnostic contracts.

# ---- Oracles --------------------------------------------------------------

# The scaled-euclidean distance matrix the default energy objective is built on:
# each covariate column is divided by its sample standard deviation before the
# pairwise Euclidean distances are formed. The comparison tests hold this matrix
# fixed across the weighted and unweighted evaluations, so the scaling constant
# never affects an inequality between two objective values.
scaled_distance_matrix <- function(covariates) {
  scales <- apply(covariates, 2, stats::sd)
  scales[scales == 0] <- 1
  scaled <- sweep(covariates, 2, scales, "/")
  as.matrix(stats::dist(scaled))
}

# The energy distance between a weighted subsample (rows `idx`, weights `w`
# normalized to sum to one) and the uniform full sample. The target self-energy
# term is constant across weightings and cancels in every comparison, but is
# kept so the value is a genuine energy distance.
energy_to_uniform <- function(dmat, idx, w) {
  n <- nrow(dmat)
  w <- w / sum(w)
  cross <- 2 * sum(w * rowSums(dmat[idx, , drop = FALSE])) / n
  within <- as.numeric(t(w) %*% dmat[idx, idx, drop = FALSE] %*% w)
  target <- sum(dmat) / n^2
  cross - within - target
}

# The between-group energy distance for the improved variant.
between_groups <- function(dmat, idx_a, idx_b, w_a, w_b) {
  w_a <- w_a / sum(w_a)
  w_b <- w_b / sum(w_b)
  cross <- 2 * as.numeric(t(w_a) %*% dmat[idx_a, idx_b, drop = FALSE] %*% w_b)
  within_a <- as.numeric(t(w_a) %*% dmat[idx_a, idx_a, drop = FALSE] %*% w_a)
  within_b <- as.numeric(t(w_b) %*% dmat[idx_b, idx_b, drop = FALSE] %*% w_b)
  cross - within_a - within_b
}

# The full improved energy objective for a binary average treatment effect: each
# group's energy distance to the target sample plus the between-group energy
# distance. The core minimizes this quantity subject to the group-sum and
# minimum-weight constraints, which uniform weights satisfy, so the solved
# weights attain a value at or below the unweighted value.
energy_ate_objective <- function(dmat, treated, w) {
  idx_t <- which(treated)
  idx_c <- which(!treated)
  energy_to_uniform(dmat, idx_t, w[idx_t]) +
    energy_to_uniform(dmat, idx_c, w[idx_c]) +
    between_groups(dmat, idx_t, idx_c, w[idx_t], w[idx_c])
}

# The weighted distance covariance between the exposure and the covariates, in
# the V-statistic form (Szekely, Rizzo, and Bakirov) with probability weights
# `p = w / sum(w)`. The continuous energy objective minimizes this quantity, so
# the balancing weights attain a smaller value than uniform weights.
weighted_distance_covariance <- function(exposure, covariates, w) {
  p <- w / sum(w)
  a <- as.matrix(stats::dist(matrix(exposure, ncol = 1)))
  b <- as.matrix(stats::dist(covariates))
  u <- as.numeric(a %*% p)
  v <- as.numeric(b %*% p)
  s1 <- as.numeric(t(p) %*% (a * b) %*% p)
  s2 <- as.numeric(t(p) %*% a %*% p) * as.numeric(t(p) %*% b %*% p)
  s3 <- sum(p * u * v)
  s1 + s2 - 2 * s3
}

# Weights normalized to mean one within each exposure group so comparisons test
# the weighting solution rather than a per-group reporting convention.
normalize_by_group <- function(w, g) {
  ave <- tapply(w, g, mean)
  as.numeric(w / ave[as.character(g)])
}

# ---- Constructor ----------------------------------------------------------

test_that("bw_energy() carries its documented defaults", {
  spec <- bw_energy()
  expect_true(S7::S7_inherits(spec, bw_energy))
  expect_true(S7::S7_inherits(spec, quadratic_program_method))
  expect_true(S7::S7_inherits(spec, balance_method))
  expect_identical(spec@distance, "scaled_euclidean")
  expect_true(spec@improved)
  expect_identical(spec@weight_penalty, 1e-4)
  expect_identical(spec@min_weight, 1e-8)
  expect_null(spec@distribution_moments)
  expect_true(spec@dimension_adjustment)
  expect_identical(spec@convergence_tolerance, 1e-6)
  expect_null(spec@max_iterations)
})

test_that("the loosened solver tolerance is energy's alone", {
  # Energy is the one quadratic program whose objective matrix is indefinite by
  # construction, and the negative curvature it carries grows as the sample
  # shrinks, so below a certain size the alternating-direction iteration stops
  # contracting: its residuals bottom out well above 1e-8 and then grow. The
  # method therefore names its own tolerance rather than taking the core default,
  # and the two positive-semidefinite quadratic programs still take it. The
  # contrast is pinned here so a later change to the shared default cannot move
  # energy with it, and so a change to energy's cannot leak into the others.
  expect_identical(bw_energy()@convergence_tolerance, 1e-6)
  expect_null(bw_sbw()@convergence_tolerance)
  expect_null(bw_cfd()@convergence_tolerance)
})

test_that("bw_energy() stores supplied tuning parameters", {
  spec <- bw_energy(
    distance = "mahalanobis",
    improved = FALSE,
    weight_penalty = 1e-3,
    min_weight = 1e-6,
    distribution_moments = 2L,
    dimension_adjustment = FALSE,
    convergence_tolerance = 1e-8,
    max_iterations = 500L
  )
  expect_identical(spec@distance, "mahalanobis")
  expect_false(spec@improved)
  expect_identical(spec@weight_penalty, 1e-3)
  expect_identical(spec@min_weight, 1e-6)
  expect_identical(spec@distribution_moments, 2L)
  expect_false(spec@dimension_adjustment)
  expect_identical(spec@convergence_tolerance, 1e-8)
  expect_identical(spec@max_iterations, 500L)
})

test_that("bw_energy() matches the distance argument", {
  expect_identical(bw_energy(distance = "euclidean")@distance, "euclidean")
  expect_error(bw_energy(distance = "manhattan"))
})

test_that("bw_energy() rejects unnamed and unknown extra arguments", {
  expect_true(S7::S7_inherits(bw_energy(), balance_method))
  expect_error(bw_energy(1e-4))
  expect_error(bw_energy(bogus = 1), class = "balancing_method_error")
})

# ---- Validators -----------------------------------------------------------

test_that("bw_energy() rejects a negative weight penalty", {
  expect_identical(bw_energy(weight_penalty = 1e-3)@weight_penalty, 1e-3)
  expect_error(bw_energy(weight_penalty = -1e-4))
})

test_that("bw_energy() rejects a negative minimum weight", {
  expect_identical(bw_energy(min_weight = 1e-6)@min_weight, 1e-6)
  expect_error(bw_energy(min_weight = -1e-8))
})

# The quadratic program pins each reweighted arm's mean weight at one, so a floor
# at one leaves the uniform weighting as the only feasible point and a floor above
# one leaves no feasible point at all. Both used to reach the solver and come back
# as an infeasibility blamed on the constraint set.
test_that("bw_energy() rejects a minimum weight at or above one", {
  expect_identical(bw_energy(min_weight = 0.5)@min_weight, 0.5)
  expect_error(bw_energy(min_weight = 1), "less than one")
  expect_error(bw_energy(min_weight = 2), "less than one")
  expect_error(bw_energy(min_weight = Inf), "less than one")
})

test_that("bw_energy() rejects a missing or infinite weight penalty", {
  expect_error(bw_energy(weight_penalty = NA_real_), "finite")
  expect_error(bw_energy(weight_penalty = Inf), "finite")
  expect_error(bw_energy(weight_penalty = c(1e-4, 1e-3)), "single")
})

# The distribution moments reach a comparison in the continuous fit path, so a
# missing value or a vector used to stop that fit with a base error rather than
# reporting at construction with the property named.
test_that("bw_energy() rejects missing or multi-element distribution moments", {
  expect_identical(
    bw_energy(distribution_moments = 2L)@distribution_moments,
    2L
  )
  expect_error(bw_energy(distribution_moments = NA_integer_), "single")
  expect_error(bw_energy(distribution_moments = c(1L, 2L)), "single")
})

test_that("bw_energy() rejects a non-positive convergence tolerance", {
  expect_identical(bw_energy()@convergence_tolerance, 1e-6)
  expect_error(bw_energy(convergence_tolerance = -1e-8))
})

test_that("bw_energy() rejects a negative iteration cap", {
  expect_null(bw_energy()@max_iterations)
  expect_error(bw_energy(max_iterations = -5L))
})

# ---- Capability methods ---------------------------------------------------

test_that("supported_exposure_types() lists every exposure type", {
  expect_setequal(
    supported_exposure_types(bw_energy()),
    c("binary", "categorical", "continuous")
  )
})

test_that("supported_estimands() depends on the exposure type", {
  binary <- supported_estimands(bw_energy(), "binary")
  expect_true(all(c("ate", "att") %in% binary))
  expect_true(any(c("atc", "atu") %in% binary))
  expect_false("ato" %in% binary)

  expect_setequal(
    supported_estimands(bw_energy(), "categorical"),
    c("ate", "att")
  )
  expect_setequal(
    supported_estimands(bw_energy(), "continuous"),
    "ate"
  )
})

test_that("supports_estimating_equations() is always FALSE for energy", {
  # The quadratic-program family has no estimating equations, whatever the
  # constraints or exposure type.
  expect_false(supports_estimating_equations(bw_energy()))
  expect_false(supports_estimating_equations(
    bw_energy(),
    constraints = balance_terms(moments = 1L)
  ))
  expect_false(supports_estimating_equations(
    bw_energy(),
    constraints = balance_terms(tolerance = 0.05)
  ))
  expect_false(supports_estimating_equations(
    bw_energy(),
    exposure_type = "binary"
  ))
  expect_false(supports_estimating_equations(
    bw_energy(),
    exposure_type = "continuous"
  ))
})

test_that("method_label() names the method", {
  expect_identical(method_label(bw_energy()), "Energy balancing")
})

# ---- Property tests: binary energy distance -------------------------------

test_that("energy balancing reduces the binary ate energy distance", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  expect_all(w, function(value) value >= 0)
  expect_all(w, function(value) value >= 1e-8)

  # Energy balancing drives balance through its objective rather than exact
  # moment constraints, so the achieved first-moment imbalance is verified
  # against the conventional good-balance ceiling rather than an exact oracle.
  expect_balanced(fit, data, tolerance = 0.1)

  dmat <- scaled_distance_matrix(as.matrix(data[c("x1", "x2")]))
  treated <- data$exposure == 1
  weighted <- energy_ate_objective(dmat, treated, w)
  unweighted <- energy_ate_objective(dmat, treated, rep(1, nrow(data)))
  expect_lt(weighted, unweighted)
})

test_that("a binary ate normalizes each group to its size", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  treated <- data$exposure == 1
  expect_equal(sum(w[treated]), sum(treated), tolerance = 1e-4)
  expect_equal(sum(w[!treated]), sum(!treated), tolerance = 1e-4)
})

test_that("a binary att targets the treated total in both groups", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "att"
  )
  w <- as.numeric(stats::weights(fit))
  treated <- data$exposure == 1
  n_treated <- sum(treated)
  expect_equal(sum(w[treated]), n_treated, tolerance = 1e-4)
  expect_equal(sum(w[!treated]), n_treated, tolerance = 1e-4)
  expect_all(w, function(value) value >= 0)
  expect_balanced(fit, data, tolerance = 0.1)
})

test_that("a binary atc fit produces non-negative weights", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "atc"
  )
  w <- as.numeric(stats::weights(fit))
  expect_all(w, function(value) value >= 0)
  expect_all(w, function(value) value >= 1e-8)
  expect_balanced(fit, data, tolerance = 0.1)
})

# A factor covariate crosses the covariate-distance boundary as one indicator
# column per level, since each level is a coordinate of that distance. The
# expansion is the only place the levels are read, so a regression in it would
# leave the objective running on the wrong coordinates while the fit still
# returned weights that normalize and balance the numeric covariates. The
# factor's own imbalance is therefore measured directly: the three level
# proportions differ between the exposure groups by up to 0.19 unweighted.

test_that("energy balancing balances a factor covariate", {
  data <- sim_binary()
  treated <- data$exposure == 1
  level_gaps <- function(w) {
    vapply(
      levels(data$x3),
      function(level) {
        indicator <- as.numeric(data$x3 == level)
        abs(
          stats::weighted.mean(indicator[treated], w[treated]) -
            stats::weighted.mean(indicator[!treated], w[!treated])
        )
      },
      numeric(1)
    )
  }
  unweighted <- level_gaps(rep(1, nrow(data)))
  expect_gt(max(unweighted), 0.15)

  for (estimand in c("ate", "att")) {
    fit <- balance(
      data,
      exposure,
      c(x1, x2, x3),
      method = bw_energy(),
      estimand = estimand
    )
    w <- as.numeric(stats::weights(fit))
    control_target <- if (identical(estimand, "ate")) {
      sum(!treated)
    } else {
      sum(treated)
    }
    expect_true(fit@converged)
    expect_equal(sum(w[treated]), sum(treated), tolerance = 1e-4)
    expect_equal(sum(w[!treated]), control_target, tolerance = 1e-4)
    expect_all(w, function(value) value >= 1e-8)

    # Every level's gap closes by at least a factor of four and lands inside a
    # ceiling no unweighted level clears.
    weighted <- level_gaps(w)
    expect_all(weighted, function(value) value < unweighted / 4)
    expect_lt(max(weighted), 0.01)
    expect_balanced(fit, data, tolerance = 0.1)
  }
})

# ---- ESS ------------------------------------------------------------------

# Kish effective sample size, computed inline from the reported weights since
# balance assessment moved to halfmoon. Bounds alone are no test of it: every
# strictly positive weight vector sits between zero and its own group size by
# Cauchy-Schwarz, so the specs below pin where the figure lands, what moves it,
# and the one case where it is exact.
kish_ess <- function(w) sum(w)^2 / sum(w^2)

# The groups come from the data, which is where the level a row belongs to is
# recorded.
exposure_groups <- function(data) {
  split(seq_len(nrow(data)), as.character(data$exposure))
}

test_that("a binary ate spends part of each group on balance", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  w <- as.numeric(weights(fit))
  for (idx in exposure_groups(data)) {
    group_ess <- kish_ess(w[idx])
    # Driving the energy distance down on this confounded exposure costs
    # precision, so the figure is strictly below the group size by a wide margin
    # rather than merely bounded by it: uniform weights, which balance nothing,
    # would sit at the size exactly. It stays well clear of the floor a handful
    # of dominating weights would leave, which is the other way a fit can fail
    # while still reporting positive weights.
    expect_lt(group_ess, 0.85 * length(idx))
    expect_gt(group_ess, 0.35 * length(idx))
  }
})

test_that("the per-group effective sample size rises with the weight penalty", {
  # The weight penalty is the L2 term that pulls the solution toward the base
  # weights, so raising it buys precision back from the energy objective. Tying
  # the figure to the knob that moves it is what distinguishes these weights
  # from any other positive vector: an arbitrary one has no reason to be ordered
  # this way in both groups at once. The penalties are spaced by decades, which
  # is what it takes for each step to move the figure by more than a rounding.
  data <- sim_binary()
  ess_at <- function(penalty) {
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(weight_penalty = penalty),
      estimand = "ate"
    )
    w <- as.numeric(weights(fit))
    vapply(exposure_groups(data), function(idx) kish_ess(w[idx]), numeric(1))
  }
  penalties <- c(1e-2, 1e-1, 1)
  curve <- lapply(penalties, ess_at)

  for (step in seq_len(length(curve) - 1L)) {
    expect_all(curve[[step + 1L]], function(value) value > curve[[step]])
  }
})

test_that("a focal group keeps its base weights, so its ESS is its size", {
  # The focal group of a treated estimand is not a variable of the quadratic
  # program at all, and renormalizing it to its own total leaves every one of
  # its weights equal to one. Its effective sample size is therefore its size
  # exactly, an equality the reweighted group cannot meet: that group carries
  # the whole tilt and comes back at a fraction of its own size.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "att"
  )
  w <- as.numeric(weights(fit))
  groups <- exposure_groups(data)
  focal <- groups[[fit@focal_level]]
  reweighted <- groups[[setdiff(names(groups), fit@focal_level)]]

  expect_equal(w[focal], rep(1, length(focal)), tolerance = 1e-10)
  expect_equal(kish_ess(w[focal]), length(focal), tolerance = 1e-8)
  expect_lt(kish_ess(w[reweighted]), 0.6 * length(reweighted))
  expect_gt(kish_ess(w[reweighted]), 0.2 * length(reweighted))
})

# ---- Categorical ----------------------------------------------------------

test_that("categorical ate energy balancing produces valid weights", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  expect_all(w, function(value) value >= 0)
  for (level in levels(data$exposure)) {
    idx <- data$exposure == level
    expect_equal(sum(w[idx]), sum(idx), tolerance = 1e-4)
  }
  expect_balanced(fit, data, tolerance = 0.1)
})

test_that("categorical att energy balancing produces valid weights", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "att",
    .focal_level = "b"
  )
  w <- as.numeric(stats::weights(fit))
  expect_all(w, function(value) value >= 0)
  expect_balanced(fit, data, tolerance = 0.1)
})

# ---- Moment constraints ---------------------------------------------------

test_that("moment constraints are satisfied within tolerance", {
  # Without moment constraints the quadratic-program family does not promise
  # exact moment balance; with them the constraint rows hold to numerical
  # precision, which the shared oracle verifies.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate",
    constraints = balance_terms(moments = 1L)
  )
  expect_balanced(fit, data)
})

test_that("a positive tolerance relaxes the moment constraints", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate",
    constraints = balance_terms(moments = 1L, tolerance = 0.1)
  )
  expect_balanced(fit, data, tolerance = 0.1)
})

test_that("a tolerance without moment constraints warns and is ignored", {
  # The tolerance relaxes added moment constraints; with none present it has
  # nothing to act on, so energy balancing warns and proceeds. The reference
  # fit without the tolerance produces the same weights.
  data <- sim_binary()
  expect_warning(
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0.1)
    ),
    class = "balancing_warning"
  )
  reference <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  expect_equal(
    as.numeric(stats::weights(fit)),
    as.numeric(stats::weights(reference)),
    tolerance = 1e-6
  )
})

test_that("quantile constraints balance a discrete exposure", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate",
    constraints = balance_terms(moments = 1L, quantiles = 0.5)
  )
  expect_balanced(fit, data)
})

# ---- Distance options -----------------------------------------------------

test_that("the distance definitions produce different weights", {
  data <- sim_binary()
  fit_of <- function(distance) {
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(distance = distance),
      estimand = "ate"
    )
  }
  w_scaled <- as.numeric(stats::weights(fit_of("scaled_euclidean")))
  w_maha <- as.numeric(stats::weights(fit_of("mahalanobis")))
  w_eucl <- as.numeric(stats::weights(fit_of("euclidean")))
  expect_false(isTRUE(all.equal(w_scaled, w_maha)))
  expect_false(isTRUE(all.equal(w_scaled, w_eucl)))
  expect_false(isTRUE(all.equal(w_maha, w_eucl)))
})

test_that("the improved variant differs from the plain variant for a binary ate", {
  data <- sim_binary()
  fit_improved <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(improved = TRUE),
    estimand = "ate"
  )
  fit_plain <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(improved = FALSE),
    estimand = "ate"
  )
  w_improved <- normalize_by_group(
    as.numeric(stats::weights(fit_improved)),
    data$exposure
  )
  w_plain <- normalize_by_group(
    as.numeric(stats::weights(fit_plain)),
    data$exposure
  )
  expect_false(isTRUE(all.equal(w_improved, w_plain)))
})

# ---- Diagnostics ----------------------------------------------------------

test_that("the fit reports dual variables and the quadratic-program backend", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  expect_s3_class(fit@duals, "data.frame")
  expect_gt(nrow(fit@duals), 0)
  # The default backend is osqp; the solver status reflects it.
  expect_identical(fit@solver_status, "osqp")
})

# ---- Backend routing ------------------------------------------------------

# Energy balancing assembles an indefinite quadratic form, so it always solves
# through the ADMM backend whatever the option asks for. What it owes the caller
# is the account its quadratic-program siblings give: an unknown value is an
# error rather than a silent default, a pinned interior-point backend is
# announced as dropped, and the recorded backend names what actually ran.

test_that("an unknown quadratic-program backend option is refused", {
  data <- sim_binary()
  withr::local_options(balancing.qp_backend = "bogus")
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(),
      estimand = "ate"
    ),
    class = "balancing_range_error"
  )
})

test_that("an unknown backend option is refused for a continuous fit", {
  data <- sim_continuous()
  withr::local_options(balancing.qp_backend = "bogus")
  expect_error(
    balance(data, exposure, c(x1, x2), method = bw_energy()),
    class = "balancing_range_error"
  )
})

test_that("a clarabel pin energy balancing cannot honor is announced", {
  withr::local_options(balancing.qp_backend = "clarabel")
  data <- sim_binary()
  expect_warning(
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(),
      estimand = "ate"
    ),
    class = "balancing_ignored_argument_warning"
  )
  # The recorded backend names what actually ran, not what was asked for.
  expect_identical(fit@solver_status, "osqp")
  expect_true(fit@converged)
})

test_that("a continuous fit announces a clarabel pin as well", {
  withr::local_options(balancing.qp_backend = "clarabel")
  data <- sim_continuous()
  expect_warning(
    fit <- balance(data, exposure, c(x1, x2), method = bw_energy()),
    class = "balancing_ignored_argument_warning"
  )
  expect_identical(fit@solver_status, "osqp")
})

test_that("energy balancing is silent under the automatic and osqp backends", {
  data <- sim_binary()
  for (backend in c("auto", "osqp")) {
    fit <- withr::with_options(
      list(balancing.qp_backend = backend),
      expect_no_warning(balance(
        data,
        exposure,
        c(x1, x2),
        method = bw_energy(),
        estimand = "ate"
      ))
    )
    expect_identical(fit@solver_status, "osqp")
  }
})

test_that("an energy fit has no estimating equations", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  expect_null(fit@estimating_equations)
  expect_error(
    estimating_equations(fit),
    class = "balancing_ipw_unsupported_error"
  )
})

test_that("ipw() rejects an energy fit", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  outcome <- stats::lm(x1 ~ exposure, data = data)
  expect_error(
    ipw(fit, outcome),
    class = "balancing_ipw_unsupported_error"
  )
})

# ---- Continuous exposure --------------------------------------------------

test_that("continuous energy balancing reduces the distance covariance", {
  data <- sim_continuous()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  expect_all(w, function(value) value >= 0)

  covariates <- as.matrix(data[c("x1", "x2")])
  weighted <- weighted_distance_covariance(data$exposure, covariates, w)
  unweighted <- weighted_distance_covariance(
    data$exposure,
    covariates,
    rep(1, nrow(data))
  )
  expect_lt(weighted, unweighted)
})

test_that("distribution_moments changes the continuous weights and holds the exposure variance", {
  # distribution_moments adds marginal-moment constraints to the continuous
  # objective, so raising it is not inert and the constrained moment is held at
  # the unweighted sample value.
  data <- sim_continuous()
  fit_default <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  fit_moments <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(distribution_moments = 2L),
    estimand = "ate"
  )
  w_default <- as.numeric(stats::weights(fit_default))
  w_moments <- as.numeric(stats::weights(fit_moments))
  expect_false(isTRUE(all.equal(w_default, w_moments)))

  sample_variance <- sum((data$exposure - mean(data$exposure))^2) /
    length(data$exposure)
  weighted_variance <- function(w) {
    center <- stats::weighted.mean(data$exposure, w)
    sum(w * (data$exposure - center)^2) / sum(w)
  }
  expect_equal(weighted_variance(w_moments), sample_variance, tolerance = 1e-3)
})

test_that("a continuous fit preserves an indicator covariate's marginal", {
  # An indicator covariate crosses the boundary as a raw zero/one column, so the
  # marginal it is held to is the sample proportion. A distribution row that
  # drives that column's weighted mean to zero instead empties the indicated
  # stratum, which balances nothing.
  data <- sim_continuous_indicator()
  fit <- balance(
    data,
    exposure,
    c(x1, g),
    method = bw_energy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  expect_equal(
    stats::weighted.mean(data$g, w),
    mean(data$g),
    tolerance = 1e-6
  )
  # The stratum keeps its share of the total weight rather than being annihilated.
  expect_equal(sum(w[data$g == 1]), sum(data$g == 1), tolerance = 1e-4)
  expect_all(w, function(value) value >= 0)
})

test_that("a continuous fit holds the base-measure marginals under sampling weights", {
  # Every distribution row of one fit targets one population. Under informative
  # sampling that population is the sampling-weighted sample, so the exposure
  # marginal moves to the sampling-weighted exposure mean alongside the
  # covariate marginals rather than staying at the unweighted mean, which would
  # leave the fit targeting two populations at once.
  withr::local_seed(42)
  n <- 300
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  exposure <- 0.8 * x1 - 0.5 * x2 + stats::rnorm(n)
  data <- data.frame(
    exposure = exposure,
    x1 = x1,
    x2 = x2,
    sw = ifelse(x1 + exposure > 0, 3, 0.5)
  )
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate",
    sampling_weights = sw
  )
  # The accessor already composes the sampling weights, so the achieved means
  # read off it directly.
  w <- as.numeric(stats::weights(fit))
  expect_equal(
    stats::weighted.mean(data$exposure, w),
    stats::weighted.mean(data$exposure, data$sw),
    tolerance = 1e-6
  )
  expect_equal(
    stats::weighted.mean(data$x1, w),
    stats::weighted.mean(data$x1, data$sw),
    tolerance = 1e-6
  )
  expect_equal(
    stats::weighted.mean(data$x2, w),
    stats::weighted.mean(data$x2, data$sw),
    tolerance = 1e-6
  )
})

test_that("a continuous fit holds an indicator's base-measure proportion", {
  # The two defects meet on an indicator covariate under sampling weights: the
  # proportion the fit preserves is the sampling-weighted one, and the exposure
  # mean is preserved on that same measure.
  data <- sim_continuous_indicator()
  withr::local_seed(31)
  data$sw <- stats::runif(nrow(data), 0.3, 3)
  fit <- balance(
    data,
    exposure,
    c(x1, g),
    method = bw_energy(),
    estimand = "ate",
    sampling_weights = sw
  )
  w <- as.numeric(stats::weights(fit))
  expect_equal(
    stats::weighted.mean(data$g, w),
    stats::weighted.mean(data$g, data$sw),
    tolerance = 1e-6
  )
  expect_equal(
    stats::weighted.mean(data$exposure, w),
    stats::weighted.mean(data$exposure, data$sw),
    tolerance = 1e-6
  )
})

test_that("the higher distribution moments are held on the base measure", {
  # The rows above the first moment take the same measure as the first. With the
  # exposure mean and its second moment both held at their sampling-weighted
  # values, the weighted exposure variance is the sampling-weighted variance.
  data <- sim_continuous()
  withr::local_seed(17)
  data$sw <- stats::runif(nrow(data), 0.3, 3)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(distribution_moments = 2L),
    estimand = "ate",
    sampling_weights = sw
  )
  w <- as.numeric(stats::weights(fit))
  weighted_variance <- function(weights) {
    center <- stats::weighted.mean(data$exposure, weights)
    sum(weights * (data$exposure - center)^2) / sum(weights)
  }
  expect_equal(
    stats::weighted.mean(data$exposure, w),
    stats::weighted.mean(data$exposure, data$sw),
    tolerance = 1e-6
  )
  expect_equal(
    weighted_variance(w),
    weighted_variance(data$sw),
    tolerance = 1e-4
  )
})

test_that("a third distribution moment is held on the base measure", {
  # The distribution rows go up in order together, so a third-moment fit holds
  # the first three central moments of the exposure and of every numeric
  # covariate at their sampling-weighted values. The third moment is what
  # separates this from the second-moment case, which matches the mean and the
  # variance and leaves the skewness where the objective puts it.
  data <- sim_continuous()
  withr::local_seed(17)
  data$sw <- stats::runif(nrow(data), 0.3, 3)
  fit_of <- function(moments) {
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(distribution_moments = moments),
      estimand = "ate",
      sampling_weights = sw
    )
  }
  central_moment <- function(values, weights, order) {
    center <- stats::weighted.mean(values, weights)
    sum(weights * (values - center)^order) / sum(weights)
  }

  w <- as.numeric(stats::weights(fit_of(3L)))
  for (column in c("exposure", "x1", "x2")) {
    values <- data[[column]]
    expect_equal(
      stats::weighted.mean(values, w),
      stats::weighted.mean(values, data$sw),
      tolerance = 1e-8
    )
    for (order in 2:3) {
      expect_equal(
        central_moment(values, w, order),
        central_moment(values, data$sw, order),
        tolerance = 1e-8
      )
    }
  }

  # The third row does real work: without it the exposure's third moment is not
  # held anywhere near the target the row pins it to.
  w_two <- as.numeric(stats::weights(fit_of(2L)))
  target <- central_moment(data$exposure, data$sw, 3)
  expect_gt(
    abs(central_moment(data$exposure, w_two, 3) - target),
    1e-4
  )
})

test_that("dimension_adjustment toggles the continuous solution", {
  data <- sim_continuous()
  fit_on <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(dimension_adjustment = TRUE),
    estimand = "ate"
  )
  fit_off <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(dimension_adjustment = FALSE),
    estimand = "ate"
  )
  expect_false(isTRUE(all.equal(
    as.numeric(stats::weights(fit_on)),
    as.numeric(stats::weights(fit_off))
  )))
})

# ---- Continuous correlation constraints -----------------------------------

# The weighted exposure-covariate Pearson correlation the balance table reports,
# computed here so the specs below judge the fit on the statistic a reader sees
# rather than on the solver's own row.
weighted_correlation <- function(exposure, column, weights) {
  abs(stats::cov.wt(cbind(exposure, column), wt = weights, cor = TRUE)$cor[
    1,
    2
  ])
}

test_that("moments requests exposure-covariate correlation constraints", {
  # For a continuous exposure `balance_terms(moments = k)` asks for the weighted
  # correlation of the exposure with each covariate power up to k to be held at
  # zero, the meaning `moments` carries for a discrete exposure and the meaning
  # WeightIt gives its own `moments` argument. The energy objective alone leaves
  # a residual correlation of a tenth or more at these sample sizes, so a fit
  # that meets this really did add the rows.
  data <- sim_continuous(n = 350)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate",
    constraints = balance_terms(moments = 1L)
  )
  w <- as.numeric(stats::weights(fit))
  for (column in c("x1", "x2")) {
    expect_lt(weighted_correlation(data$exposure, data[[column]], w), 1e-6)
  }

  default <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  w_default <- as.numeric(stats::weights(default))
  expect_gt(weighted_correlation(data$exposure, data$x1, w_default), 0.05)
})

test_that("a second moment constrains the correlation with the covariate squares", {
  # Each power of a covariate is its own constraint column, so the second moment
  # adds the correlation of the exposure with the squares alongside the first.
  data <- sim_continuous(n = 350)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate",
    constraints = balance_terms(moments = 2L)
  )
  w <- as.numeric(stats::weights(fit))
  for (column in c("x1", "x2")) {
    values <- data[[column]]
    expect_lt(weighted_correlation(data$exposure, values, w), 1e-6)
    expect_lt(
      weighted_correlation(data$exposure, (values - mean(values))^2, w),
      1e-6
    )
  }
})

test_that("interactions constrain the correlation with the product column", {
  # `interactions = TRUE` adds the product of two covariates as a constraint
  # column, which on the continuous path is one more correlation row.
  data <- sim_continuous(n = 350)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate",
    constraints = balance_terms(interactions = TRUE)
  )
  w <- as.numeric(stats::weights(fit))
  expect_lt(
    weighted_correlation(data$exposure, data$x1 * data$x2, w),
    1e-6
  )
})

test_that("the balance table reports the correlation rows the fit constrained", {
  # The rows a continuous fit constrains are the rows the table reports, so a
  # second-moment request shows four correlation rows, each met within its
  # tolerance.
  data <- sim_continuous(n = 350)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate",
    constraints = balance_terms(moments = 2L)
  )
  table <- as.data.frame(fit@balance_table)
  expect_identical(table$term, c("x1", "x1^2", "x2", "x2^2"))
  expect_true(all(table$statistic == "correlation"))
  expect_lt(max(table$weighted), 1e-6)
  expect_true(all(table$within_tolerance))
})

# A continuous fit's achieved weighted exposure-covariate correlations, the
# statistic the correlation rows are judged on and the one the balance table
# reports.
achieved_correlations <- function(fit, data, covariates) {
  w <- as.numeric(stats::weights(fit))
  vapply(
    covariates,
    function(covariate) {
      weighted_correlation(data$exposure, data[[covariate]], w)
    },
    numeric(1)
  )
}

test_that("a continuous tolerance is honored as a band rather than held exactly", {
  # The correlation rows are relaxable. The quadratic program bounds a
  # linearized correlation whose exposure and covariate scales are fixed at the
  # sample, so a single solve at the requested band overshoots it; the fit
  # tightens the bound it hands the program over a few passes until the reported
  # correlation sits inside the band. What that has to show is both halves: no
  # correlation above the band, and at least one well inside it rather than at
  # the zero exact rows would produce, so the band is used rather than ignored.
  data <- sim_continuous(n = 350)
  fit <- expect_no_warning(balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate",
    constraints = balance_terms(moments = 1L, tolerance = 0.05)
  ))
  achieved <- achieved_correlations(fit, data, c("x1", "x2"))

  expect_all(achieved, function(r) r <= 0.05 + balance_margin(0.05))
  expect_gt(max(achieved), 0.5 * 0.05)
})

test_that("a continuous tolerance of zero still holds the rows exactly", {
  # Exact balance is the tolerance the refinement has nothing to tighten, so it
  # reaches the same solution it always did, in a single solve.
  data <- sim_continuous(n = 350)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate",
    constraints = balance_terms(moments = 1L)
  )
  achieved <- achieved_correlations(fit, data, c("x1", "x2"))

  expect_all(achieved, function(r) r < 1e-6)
})

test_that("the refinement takes several passes and sums their iterations", {
  # The pass count is what separates an honored band from a single overshooting
  # solve, so it is counted at the solver rather than inferred from the weights,
  # and the reported iterations have to account for every pass rather than for
  # the last one alone. Exact balance is the control: it has nothing to tighten
  # and takes one solve.
  data <- sim_continuous(n = 350)
  solves <- 0L
  per_solve <- integer()
  original <- solve_energy_cont
  testthat::local_mocked_bindings(
    solve_energy_cont = function(...) {
      result <- original(...)
      solves <<- solves + 1L
      per_solve <<- c(per_solve, as.integer(result$iterations))
      result
    }
  )

  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate",
    constraints = balance_terms(moments = 1L, tolerance = 0.05)
  )
  expect_gt(solves, 1L)
  expect_identical(fit@iterations, sum(per_solve))

  solves <- 0L
  balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate",
    constraints = balance_terms(moments = 1L)
  )
  expect_identical(solves, 1L)
})

test_that("a tightened pass that stops at the iteration cap keeps the last converged iterate", {
  # A tightened bound is harder than the one before it, so a pass can spend the
  # iteration cap on a band an earlier pass met comfortably. What the fit reports
  # then is the last iterate that converged rather than the unsettled one the
  # tightened pass left behind: the weights are the earlier pass's, the fit
  # reports itself converged, and the correlations that iterate actually achieved
  # sit above the requested band, so the ordinary balance warning judges them
  # instead of a convergence warning claiming the solve failed. The reported
  # iterations still account for every pass, the failed one included, because
  # each pass cost a whole solve.
  data <- sim_continuous(n = 350)
  solves <- 0L
  per_solve <- integer()
  first_weights <- NULL
  original <- solve_energy_cont
  testthat::local_mocked_bindings(
    solve_energy_cont = function(...) {
      result <- original(...)
      solves <<- solves + 1L
      per_solve <<- c(per_solve, as.integer(result$iterations))
      if (solves == 1L) {
        first_weights <<- as.numeric(result$weights)
      } else {
        result$converged <- FALSE
        result$status <- "max_iter"
      }
      result
    }
  )

  # Every warning is collected rather than one being matched, so a convergence
  # warning raised alongside the balance warning fails here instead of passing
  # under an expectation that looked only for the one it wanted.
  seen <- character()
  fit <- withCallingHandlers(
    balance(
      data,
      exposure,
      c(x1, x2),
      # Pinned at the tolerance the problem is taken to reach, so the mocked
      # failure cannot draw in the re-solve fallback and the two solves counted
      # here are the two refinement passes.
      method = bw_energy(convergence_tolerance = 1e-6),
      estimand = "ate",
      constraints = balance_terms(moments = 1L, tolerance = 0.05)
    ),
    warning = function(cnd) {
      seen <<- c(seen, class(cnd)[[1]])
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(seen, "balancing_balance_warning")
  expect_identical(solves, 2L)
  expect_true(fit@converged)
  expect_identical(fit@iterations, sum(per_solve))

  # The reported weights renormalize the solver's, so the restored iterate shows
  # as proportionality to the first pass's raw weights rather than as equality.
  w <- as.numeric(stats::weights(fit))
  expect_equal(
    w / sum(w),
    first_weights / sum(first_weights),
    tolerance = 1e-10
  )

  table <- as.data.frame(fit@balance_table)
  expect_all(table$weighted, function(r) r > 0.05 + balance_margin(0.05))
  expect_column_all(table, "within_tolerance", function(value) !value)
})

test_that("a tightened pass certified infeasible raises rather than restoring", {
  # The other half of the same guard. An infeasibility certificate is a claim
  # about the constraint set, not about the iteration cap, so restoring the
  # looser iterate would answer a band the solver said cannot be met with
  # weights that do not meet it. The certificate surfaces as the infeasible
  # condition instead.
  data <- sim_continuous(n = 350)
  solves <- 0L
  original <- solve_energy_cont
  testthat::local_mocked_bindings(
    solve_energy_cont = function(...) {
      result <- original(...)
      solves <<- solves + 1L
      if (solves > 1L) {
        result$converged <- FALSE
        result$status <- "primal_infeasible"
      }
      result
    }
  )

  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(),
      estimand = "ate",
      constraints = balance_terms(moments = 1L, tolerance = 0.05)
    ),
    class = "balancing_infeasible_error"
  )
  expect_identical(solves, 2L)
})

test_that("moments no longer sets the marginal distribution moments", {
  # `distribution_moments` is the only route to the marginal rows. A
  # second-moment constraint request therefore leaves the weighted exposure
  # variance where the objective puts it, while `distribution_moments = 2` pins
  # it at the sample value.
  data <- sim_continuous(n = 350)
  sample_variance <- sum((data$exposure - mean(data$exposure))^2) / nrow(data)
  weighted_variance <- function(w) {
    center <- stats::weighted.mean(data$exposure, w)
    sum(w * (data$exposure - center)^2) / sum(w)
  }
  fit_terms <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate",
    constraints = balance_terms(moments = 2L)
  )
  fit_distribution <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(distribution_moments = 2L),
    estimand = "ate"
  )
  w_terms <- as.numeric(stats::weights(fit_terms))
  w_distribution <- as.numeric(stats::weights(fit_distribution))
  expect_gt(abs(weighted_variance(w_terms) - sample_variance), 1e-3)
  expect_equal(
    weighted_variance(w_distribution),
    sample_variance,
    tolerance = 1e-3
  )
})

test_that("a covariate left out of the constraint set keeps its marginal rows", {
  # The marginal rows belong to the covariates, not to the constraint set, so
  # `moments` must not reach them from either side. Excluding x1 from the
  # constraint set drops x1's correlation row and nothing else: its weighted mean
  # and variance stay at the sample values `distribution_moments` pins, exactly
  # as x2's and the exposure's do. Reading the marginal columns off the
  # constraint recipe dropped x1's marginal rows along with its correlation row,
  # and its weighted variance then floated to wherever the objective put it.
  data <- sim_continuous(n = 350)
  central_moment <- function(values, weights, order) {
    center <- stats::weighted.mean(values, weights)
    sum(weights * (values - center)^order) / sum(weights)
  }
  uniform <- rep(1, nrow(data))

  # The constraint sets differ in what they ask of the correlation rows and
  # agree in what they leave to the marginals, so the marginals must come out
  # the same under all three. The middle one is the case that failed: x1 has no
  # constraint record to read a marginal off. The last one is its mirror, where
  # x1's record reaches past the distribution moments.
  constraint_sets <- list(
    balance_terms(moments = 1L),
    balance_terms(moments = c(x1 = 0L, x2 = 1L)),
    balance_terms(moments = 3L)
  )

  for (constraints in constraint_sets) {
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(distribution_moments = 2L),
      estimand = "ate",
      constraints = constraints
    )
    w <- as.numeric(stats::weights(fit))
    for (column in c("exposure", "x1", "x2")) {
      values <- data[[column]]
      expect_equal(
        stats::weighted.mean(values, w),
        mean(values),
        tolerance = 1e-8
      )
      expect_equal(
        central_moment(values, w, 2),
        central_moment(values, uniform, 2),
        tolerance = 1e-8
      )
    }
  }
})

test_that("the default continuous fit keeps the objective-driven solution", {
  # The correlation rows are added only when the constraint set asks for them,
  # so the default fit solves the same program it solved before: the weights
  # average one, the correlation improves on the unweighted sample, and the
  # residual the energy objective leaves behind is still there.
  data <- sim_continuous(n = 350)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  expect_equal(mean(w), 1, tolerance = 1e-8)
  expect_gt(weighted_correlation(data$exposure, data$x1, w), 0.05)
  expect_lt(
    weighted_correlation(data$exposure, data$x1, w),
    abs(stats::cor(data$exposure, data$x1))
  )
})

# ---- Unsupported estimands ------------------------------------------------

# The overlap estimand is legal only for the covariate balancing propensity
# score, so energy balancing rejects it for every exposure type with the shared
# unsupported-estimand condition, whose message the entropy slice already
# snapshots.

test_that("the ato estimand raises balancing_estimand_error for a binary exposure", {
  data <- sim_binary()
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(),
      estimand = "ato"
    ),
    class = "balancing_estimand_error"
  )
})

test_that("the ato estimand raises balancing_estimand_error for a categorical exposure", {
  data <- sim_categorical()
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(),
      estimand = "ato"
    ),
    class = "balancing_estimand_error"
  )
})

# ---- Newly reachable classed condition ------------------------------------

test_that("the ignored-tolerance warning records its class and message", {
  data <- sim_binary(n = 150)
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0.1)
    )
  )
})

test_that("a continuous tolerance warns and is ignored", {
  # A continuous fit holds its distribution moments exactly and adds no relaxable
  # moment constraints, so a positive tolerance has nothing to act on and warns.
  data <- sim_continuous()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0.1)
    )
  )
})

test_that("the balance table reports the tolerance a continuous fit enforced", {
  # The correlation rows are held inside the band the specification asked for,
  # so that band is what the table reports and what its verdict is judged
  # against.
  data <- sim_continuous(n = 350)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate",
    constraints = balance_terms(moments = 1L, tolerance = 0.05)
  )
  table <- as.data.frame(fit@balance_table)
  expect_column_all(table, "tolerance", function(value) value == 0.05)
  expect_column_all(table, "within_tolerance", function(value) value)
})

test_that("a fit that added no constraint rows reports the tolerance it enforced", {
  # A tolerance with no constraint rows to relax reaches no row of the program,
  # so the table must not report it as the fit's tolerance: the fit enforced
  # nothing, which is a tolerance of zero. Both exposure types answer the same
  # way, the rows being absent for the same reason in each.
  for (data in list(sim_binary(n = 200), sim_continuous(n = 200))) {
    expect_warning(
      fit <- balance(
        data,
        exposure,
        c(x1, x2),
        method = bw_energy(),
        estimand = "ate",
        constraints = balance_terms(tolerance = 0.1)
      ),
      class = "balancing_ignored_argument_warning"
    )
    expect_column_all(
      as.data.frame(fit@balance_table),
      "tolerance",
      function(value) value == 0
    )
  }
})

# ---- Infeasible constraint set --------------------------------------------

test_that("an infeasible constraint set raises balancing_infeasible_error", {
  # A covariate that perfectly separates the exposure groups cannot meet exact
  # moment balance: the treated group's weighted mean of the separating column is
  # fixed at one while the pooled target is a proportion below one, so the solver
  # reports primal infeasibility, which maps to the infeasible condition rather
  # than the generic convergence warning.
  n <- 60L
  data <- data.frame(
    exposure = rep(0:1, each = n / 2),
    x1 = rep(0:1, each = n / 2),
    x2 = rep_len(c(-1, 0, 1), n)
  )
  expect_balancing_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(),
      estimand = "ate",
      constraints = balance_terms(moments = 1L)
    )
  )
})

# ---- Small-sample convergence ---------------------------------------------

# A frame with the shape the discrete energy objective destabilizes on: a
# three-level character covariate, a coarse numeric near 8e4 whose spread is small
# beside its mean, a numeric near 80, and a binary exposure with roughly a third
# of the sample treated, at n = 354. The energy quadratic form is the negated
# distance matrix, conditionally positive semidefinite, so it is indefinite and
# its negative curvature scales as 1/n. At this size the alternating-direction
# iteration is no longer a contraction: its residuals reach a floor above 1e-8,
# and a tolerance below that floor keeps the run going past the optimum until the
# iterate it carries is renormalized back to uniform weights, which balance
# nothing. The seed belongs to the fixture because whether a draw destabilizes
# before it clears its tolerance depends on the draw.
make_energy_frame <- function(seed, n = 354) {
  withr::with_seed(seed, {
    season <- sample(
      c("peak", "regular", "value"),
      n,
      TRUE,
      prob = c(0.22, 0.55, 0.23)
    )
    shift <- c(peak = 1.4, regular = 0, value = -1.1)[season]
    close <- sample(
      c(59400, 64800, 72000, 75600, 79200, 82800, 86400, 90000),
      n,
      TRUE,
      prob = c(0.02, 0.03, 0.16, 0.16, 0.2, 0.2, 0.15, 0.08)
    )
    close <- pmin(pmax(close + 3600 * round(shift), 59400), 90000)
    temp <- 82 + 4 * shift + stats::rnorm(n, 0, 8)
    lp <- -0.9 +
      0.55 * (season == "peak") -
      0.35 * (season == "value") +
      0.9 * scale(close)[, 1] -
      0.5 * scale(temp)[, 1]
    data.frame(
      z = stats::rbinom(n, 1, stats::plogis(lp)),
      season = season,
      close = close,
      temp = temp
    )
  })
}

test_that("the default tolerance fits a small indefinite energy problem", {
  # The default has to be a tolerance the objective can actually reach on an
  # ordinary sample of this size, so the fit converges without warning and moves
  # every covariate a long way toward balance. Under a tolerance below the
  # solver's residual floor the same fit spends its whole iteration cap instead.
  data <- make_energy_frame(9)
  fit <- expect_no_warning(
    balance(data, z, c(season, close, temp), method = bw_energy()),
    class = "balancing_convergence_warning"
  )
  expect_true(fit@converged)
  table <- as.data.frame(fit@balance_table)
  expect_gt(max(abs(table$unweighted)), 0.5)
  expect_lt(max(abs(table$weighted)), 0.1)
})

test_that("a tolerance below the residual floor spends the iteration cap", {
  # The counterpart of the spec above. The tolerance the default used to carry on
  # this frame was 1e-8, which sits so close to the solver's residual floor here
  # that scaling one covariate column by 1 + 5e-16 flips the verdict; a platform
  # whose compiler contracts a multiply-add differently would disagree with this
  # machine. The spec therefore asks for 1e-14, far below any residual floor the
  # iteration reaches, so every platform spends the cap for the same reason. The
  # cap is set well above the count a reachable tolerance converges in on this
  # frame, so what the run fails on is the tolerance rather than the budget; the
  # default cap of 200000 reaches the same verdict and costs two orders of
  # magnitude more time.
  data <- make_energy_frame(9)
  expect_warning(
    fit <- balance(
      data,
      z,
      c(season, close, temp),
      method = bw_energy(convergence_tolerance = 1e-14, max_iterations = 1000L)
    ),
    class = "balancing_convergence_warning"
  )
  expect_false(fit@converged)
})

test_that("a fit that cannot reach its tolerance returns usable weights", {
  # A run whose iterate walked away from the optimum must not hand that iterate
  # back: the stable balancing weights continuous path keeps the last iterate that
  # met its criterion rather than the failed one, and energy owes the same. The
  # promises are the ones a caller can check on the returned object. The weights
  # are finite, sit at or above the documented floor, and carry each group at its
  # estimand target total. Their effective sample size is strictly inside the
  # range a real fit occupies, below the group size that uniform weights sit at
  # exactly and above the collapse a few dominating weights would leave. And the
  # balance the fit reports improves on the unweighted sample rather than
  # reproducing it. The tolerance is the 1e-14 of the spec above, and for the
  # same reason: at the 1e-8 the default used to carry, whether this frame
  # reaches its tolerance turns on the last bits of the draw.
  data <- make_energy_frame(9)
  expect_warning(
    fit <- balance(
      data,
      z,
      c(season, close, temp),
      method = bw_energy(convergence_tolerance = 1e-14, max_iterations = 1000L)
    ),
    class = "balancing_convergence_warning"
  )

  w <- as.numeric(weights(fit))
  expect_all(w, is.finite)
  expect_all(w, function(value) value >= bw_energy()@min_weight)
  groups <- split(seq_len(nrow(data)), as.character(data$z))
  for (idx in groups) {
    expect_equal(sum(w[idx]), length(idx))
    group_ess <- kish_ess(w[idx])
    expect_lt(group_ess, 0.95 * length(idx))
    expect_gt(group_ess, 0.25 * length(idx))
  }

  table <- as.data.frame(fit@balance_table)
  expect_lt(max(abs(table$weighted)), 0.5 * max(abs(table$unweighted)))
})

test_that("a continuous fit that spends its cap reports the re-solve", {
  # The re-solve at a reachable tolerance is not a discrete-path device: the
  # continuous solve routes through the same fallback, and a run that ends at
  # the iteration cap with a tolerance below the reachable 1e-6 is retried
  # there. The fit still calls itself unconverged, because the tolerance asked
  # for was not met, and its reported iterations are the two solves added
  # together, so they exceed the cap the caller set.
  #
  # A small cap is what reaches this. The benchmark sweep that measured the
  # continuous path against WeightIt never reached it at realistic settings:
  # not one of eighty cells, up to a thousand observations at tolerances from
  # 1e-5 to 1e-8, spent its cap. Driving it with the cap is legitimate all the
  # same, since the fallback keys on the terminal status and the tolerance
  # rather than on how the cap was reached.
  data <- sim_continuous(n = 200)
  cap <- 100L
  expect_warning(
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(
        convergence_tolerance = 1e-14,
        max_iterations = cap
      ),
      estimand = "ate"
    ),
    class = "balancing_convergence_warning"
  )
  expect_false(fit@converged)
  expect_gt(fit@iterations, cap)

  # The iterate the fit reports is the re-solve's, which is a real answer: the
  # weights are finite, sit at or above the documented floor, carry the sample
  # at its target total, and improve on the unweighted correlation.
  w <- as.numeric(stats::weights(fit))
  expect_all(w, is.finite)
  expect_all(w, function(value) value >= bw_energy()@min_weight)
  expect_equal(mean(w), 1, tolerance = 1e-6)
  table <- as.data.frame(fit@balance_table)
  expect_lt(max(abs(table$weighted)), 0.5 * max(abs(table$unweighted)))
})

test_that("a continuous cap too small for the re-solve keeps the first iterate", {
  # The retry is given the same cap, so a cap below what the reachable
  # tolerance needs leaves it unconverged too, and the fit reports the original
  # solve rather than a second failed one. The iteration count is then the cap
  # itself rather than the sum.
  data <- sim_continuous(n = 200)
  cap <- 25L
  expect_warning(
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(
        convergence_tolerance = 1e-14,
        max_iterations = cap
      ),
      estimand = "ate"
    ),
    class = "balancing_convergence_warning"
  )
  expect_false(fit@converged)
  expect_identical(fit@iterations, cap)
})

# ---- Live consistency against WeightIt ------------------------------------

test_that("energy weights meet the objective tolerance against WeightIt for a binary ate", {
  skip_on_cran()
  skip_if_not_installed("WeightIt")

  # The quadratic-program family compares objectives rather than weights: both
  # implementations minimize the energy distance, so ours must attain a value at
  # or below WeightIt's on the same distance matrix. The objective is evaluated
  # in R for both weight vectors from the shared oracle.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  reference <- WeightIt::weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "energy",
    estimand = "ATE"
  )

  dmat <- scaled_distance_matrix(as.matrix(data[c("x1", "x2")]))
  treated <- data$exposure == 1
  ours <- energy_ate_objective(
    dmat,
    treated,
    as.numeric(stats::weights(fit))
  )
  theirs <- energy_ate_objective(dmat, treated, reference$weights)
  expect_lte(ours, theirs + 1e-6)
})

test_that("energy weights meet the objective tolerance against WeightIt for a binary att", {
  skip_on_cran()
  skip_if_not_installed("WeightIt")

  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "att"
  )
  reference <- WeightIt::weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "energy",
    estimand = "ATT",
    focal = "1"
  )

  # The focal fit minimizes the energy distance between the reweighted controls
  # and the treated sample, not the improved average-treatment-effect functional,
  # so both weight vectors are scored on that focal objective with the treated
  # group held at unit weight.
  dmat <- scaled_distance_matrix(as.matrix(data[c("x1", "x2")]))
  idx_t <- which(data$exposure == 1)
  idx_c <- which(data$exposure == 0)
  focal_objective <- function(w) {
    between_groups(dmat, idx_c, idx_t, w[idx_c], rep(1, length(idx_t)))
  }
  ours <- focal_objective(as.numeric(stats::weights(fit)))
  theirs <- focal_objective(reference$weights)
  expect_lte(ours, theirs + 1e-6)
})

# ---- Print snapshot -------------------------------------------------------

test_that("an energy fit prints its summary block", {
  # Records on the first successful run once the fit path exists.
  data <- sim_binary()
  expect_balancing_snapshot({
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(),
      estimand = "ate"
    )
    fit
  })
})
