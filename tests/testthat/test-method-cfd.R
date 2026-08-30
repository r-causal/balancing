# bw_cfd() is the method spec; balance(..., method = bw_cfd()) fits it.
# Characteristic function distance balancing (kernel balancing) chooses weights
# that minimize a kernel measure of the distance between the reweighted exposure
# groups and the target sample. Every kernel but energy is positive semidefinite
# by construction, so the objective is a convex quadratic program with a
# simplex-type constraint set; the energy kernel is conditionally positive
# semidefinite alone, so its quadratic form is indefinite and the solve routes to
# the ADMM backend. The method spec carries only tuning
# parameters. The improved variant adds the mean-embedding target adjustment of
# the CFD literature. The energy kernel is the negative pairwise distance, so
# bw_cfd(kernel = "energy") reproduces bw_energy() on the same data
# and constraints; that equivalence is the parity anchor for this slice, since
# there is no external reference implementation for the other kernels.
#
# Characteristic function distance balancing supports binary and categorical
# exposures only; a continuous exposure is unsupported. Like the rest of the
# quadratic-program family it has no estimating equations, so
# estimating_equations() and ipw() raise the shared unsupported condition, and
# its central guarantee is objective-level rather than exact moment balance:
# without moment constraints the kernel objective drives balance, and with them
# the constraint rows hold within tolerance. These specs therefore verify the
# energy-kernel equivalence directly, verify constraint satisfaction with the
# shared expect_balanced() oracle, exercise every kernel for convergent, valid,
# distinct weights, and check the capability and diagnostic contracts.

# ---- Oracles --------------------------------------------------------------

# Weights normalized to mean one within each exposure group so comparisons test
# the weighting solution rather than a per-group reporting convention.
normalize_by_group <- function(w, g) {
  ave <- tapply(w, g, mean)
  as.numeric(w / ave[as.character(g)])
}

# TRUE when two weight vectors are not numerically equal, used to assert that a
# tuning change actually moves the solution.
weights_differ <- function(a, b) {
  !isTRUE(all.equal(a, b))
}

# ---- Constructor ----------------------------------------------------------

test_that("bw_cfd() carries its documented defaults", {
  spec <- bw_cfd()
  expect_true(S7::S7_inherits(spec, bw_cfd))
  expect_true(S7::S7_inherits(spec, quadratic_program_method))
  expect_true(S7::S7_inherits(spec, balance_method))
  expect_identical(spec@kernel, "gaussian")
  expect_identical(spec@smoothness, 1.5)
  expect_equal(spec@degrees_of_freedom, 5)
  expect_equal(spec@simulation_draws, 5000)
  expect_true(spec@improved)
  expect_identical(spec@weight_penalty, 1e-4)
  expect_identical(spec@min_weight, 1e-8)
  expect_null(spec@convergence_tolerance)
  expect_null(spec@max_iterations)
})

test_that("bw_cfd() stores supplied tuning parameters", {
  spec <- bw_cfd(
    kernel = "matern",
    smoothness = 2.5,
    degrees_of_freedom = 8,
    simulation_draws = 1000,
    improved = FALSE,
    weight_penalty = 1e-3,
    min_weight = 1e-6,
    convergence_tolerance = 1e-8,
    max_iterations = 500L
  )
  expect_identical(spec@kernel, "matern")
  expect_identical(spec@smoothness, 2.5)
  expect_equal(spec@degrees_of_freedom, 8)
  expect_equal(spec@simulation_draws, 1000)
  expect_false(spec@improved)
  expect_identical(spec@weight_penalty, 1e-3)
  expect_identical(spec@min_weight, 1e-6)
  expect_identical(spec@convergence_tolerance, 1e-8)
  expect_identical(spec@max_iterations, 500L)
})

test_that("bw_cfd() matches the kernel argument", {
  expect_identical(bw_cfd(kernel = "laplace")@kernel, "laplace")
  expect_identical(bw_cfd(kernel = "t")@kernel, "t")
  expect_identical(bw_cfd(kernel = "energy")@kernel, "energy")
  expect_error(bw_cfd(kernel = "cauchy"))
})

test_that("bw_cfd() rejects unnamed and unknown extra arguments", {
  expect_true(S7::S7_inherits(bw_cfd(), balance_method))
  expect_error(bw_cfd(1e-4))
  expect_error(bw_cfd(bogus = 1), class = "balancing_method_error")
})

# ---- Validators -----------------------------------------------------------

test_that("bw_cfd() accepts only the supported Matern smoothness values", {
  expect_identical(bw_cfd(smoothness = 0.5)@smoothness, 0.5)
  expect_identical(bw_cfd(smoothness = 1.5)@smoothness, 1.5)
  expect_identical(bw_cfd(smoothness = 2.5)@smoothness, 2.5)
  expect_error(bw_cfd(smoothness = 1))
  expect_error(bw_cfd(smoothness = 3.5))
})

test_that("bw_cfd() rejects degrees of freedom at or below two", {
  expect_equal(bw_cfd(degrees_of_freedom = 3)@degrees_of_freedom, 3)
  expect_error(bw_cfd(degrees_of_freedom = 2))
  expect_error(bw_cfd(degrees_of_freedom = 1))
})

# An infinity used to pass the sign and missingness checks, since `is.na(Inf)` is
# FALSE and every infinity is greater than two. The t kernel's radial draw then
# took a chi-square on infinite degrees of freedom, which is a missing value, so
# every Monte Carlo frequency projection came out non-finite and the fit failed
# blaming collinearity. All three non-finite values report at construction.
test_that("bw_cfd() rejects a non-finite degrees of freedom", {
  expect_equal(bw_cfd(degrees_of_freedom = 3)@degrees_of_freedom, 3)
  expect_error(bw_cfd(degrees_of_freedom = Inf), "finite")
  expect_error(bw_cfd(degrees_of_freedom = NA_real_), "finite")
  expect_error(bw_cfd(degrees_of_freedom = NaN), "finite")
})

test_that("bw_cfd() rejects a multi-element degrees of freedom", {
  expect_error(bw_cfd(degrees_of_freedom = c(3, 5)), "single")
})

# The t-kernel projections are what a non-finite degrees of freedom poisons, so
# the rejection is also pinned on the draw itself: every projection a legal
# specification produces is finite.
test_that("cfd_projection() draws finite projections", {
  spec <- bw_cfd(kernel = "t", degrees_of_freedom = 5, simulation_draws = 50)
  projection <- withr::with_seed(414, cfd_projection(spec, 3))
  expect_length(projection, 3 * 50)
  expect_true(all(is.finite(projection)))
})

# The validator and the kernel's own boundary have to admit exactly the same set.
# The validator compared with `all.equal()`, whose relative tolerance is wider
# than the kernel's absolute 1e-9, so a value in between passed construction and
# stopped the fit at the boundary with an unclassed error. The constructor now
# snaps a request inside the validator's neighborhood onto the order it names, so
# whatever construction accepts the kernel accepts.
test_that("bw_cfd() snaps a near-canonical smoothness onto the exact order", {
  expect_identical(bw_cfd(smoothness = 1.5 + 1e-8)@smoothness, 1.5)
  expect_identical(bw_cfd(smoothness = 0.5 - 1e-9)@smoothness, 0.5)
  expect_identical(bw_cfd(smoothness = 2.5 + 1e-9)@smoothness, 2.5)
  expect_identical(bw_cfd(smoothness = 3 / 2)@smoothness, 1.5)
})

test_that("bw_cfd() refuses a smoothness outside every order's neighborhood", {
  expect_error(bw_cfd(smoothness = 1.5 + 1e-6), "0.5, 1.5, or 2.5")
  expect_error(bw_cfd(smoothness = 2), "0.5, 1.5, or 2.5")
  expect_error(bw_cfd(smoothness = NA_real_), "0.5, 1.5, or 2.5")
  expect_error(bw_cfd(smoothness = c(0.5, 1.5)), "0.5, 1.5, or 2.5")
})

test_that("a smoothness the constructor accepts always reaches the kernel", {
  data <- sim_binary(n = 150)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(kernel = "matern", smoothness = 1.5 + 1e-8),
    estimand = "ate"
  )
  expect_identical(fit@method@smoothness, 1.5)
  expect_all(as.numeric(stats::weights(fit)), is.finite)
})

test_that("bw_cfd() rejects a negative weight penalty", {
  expect_identical(bw_cfd(weight_penalty = 1e-3)@weight_penalty, 1e-3)
  expect_error(bw_cfd(weight_penalty = -1e-4))
})

test_that("bw_cfd() rejects a negative minimum weight", {
  expect_identical(bw_cfd(min_weight = 1e-6)@min_weight, 1e-6)
  expect_error(bw_cfd(min_weight = -1e-8))
})

# The quadratic program pins each reweighted arm's mean weight at one, so a floor
# at one leaves the uniform weighting as the only feasible point and a floor above
# one leaves no feasible point at all. Both used to reach the solver and come back
# as an infeasibility blamed on the constraint set.
test_that("bw_cfd() rejects a minimum weight at or above one", {
  expect_identical(bw_cfd(min_weight = 0.5)@min_weight, 0.5)
  expect_error(bw_cfd(min_weight = 1), "less than one")
  expect_error(bw_cfd(min_weight = 2), "less than one")
  expect_error(bw_cfd(min_weight = Inf), "less than one")
})

test_that("bw_cfd() rejects a non-positive number of simulation draws", {
  expect_equal(bw_cfd(simulation_draws = 2000)@simulation_draws, 2000)
  expect_error(bw_cfd(simulation_draws = 0))
  expect_error(bw_cfd(simulation_draws = -100))
})

test_that("bw_cfd() rejects a non-positive convergence tolerance", {
  expect_null(bw_cfd()@convergence_tolerance)
  expect_error(bw_cfd(convergence_tolerance = -1e-8))
})

test_that("bw_cfd() rejects a negative iteration cap", {
  expect_null(bw_cfd()@max_iterations)
  expect_error(bw_cfd(max_iterations = -5L))
})

# ---- Capability methods ---------------------------------------------------

test_that("supported_exposure_types() lists binary and categorical only", {
  expect_setequal(
    supported_exposure_types(bw_cfd()),
    c("binary", "categorical")
  )
  expect_false("continuous" %in% supported_exposure_types(bw_cfd()))
})

test_that("supported_estimands() depends on the exposure type", {
  binary <- supported_estimands(bw_cfd(), "binary")
  expect_true(all(c("ate", "att") %in% binary))
  expect_true(any(c("atc", "atu") %in% binary))
  expect_false("ato" %in% binary)

  expect_setequal(
    supported_estimands(bw_cfd(), "categorical"),
    c("ate", "att")
  )
})

test_that("supports_estimating_equations() is always FALSE for cfd", {
  # The quadratic-program family has no estimating equations, whatever the
  # constraints or exposure type.
  expect_false(supports_estimating_equations(bw_cfd()))
  expect_false(supports_estimating_equations(
    bw_cfd(),
    constraints = balance_terms(moments = 1L)
  ))
  expect_false(supports_estimating_equations(
    bw_cfd(),
    constraints = balance_terms(tolerance = 0.05)
  ))
  expect_false(supports_estimating_equations(
    bw_cfd(),
    exposure_type = "binary"
  ))
})

test_that("method_label() names the method", {
  expect_identical(
    method_label(bw_cfd()),
    "Characteristic function distance balancing"
  )
})

# ---- Energy-kernel equivalence: the parity anchor -------------------------

test_that("the energy kernel reproduces energy balancing for a binary ate", {
  # bw_cfd(kernel = "energy") is the negative pairwise distance, so with the
  # matched defaults (improved variant, weight penalty, minimum weight) it solves
  # the same quadratic program as bw_energy() on the same data and
  # constraints and returns the same weights. This is the parity anchor for the
  # slice: the other kernels share the assembly and differ only in the kernel
  # matrix.
  data <- sim_binary()
  fit_cfd <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(kernel = "energy"),
    estimand = "ate"
  )
  fit_energy <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  expect_equal(
    as.numeric(stats::weights(fit_cfd)),
    as.numeric(stats::weights(fit_energy)),
    tolerance = 1e-5
  )
})

test_that("the energy kernel reproduces energy balancing for a binary att", {
  data <- sim_binary()
  fit_cfd <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(kernel = "energy"),
    estimand = "att"
  )
  fit_energy <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "att"
  )
  expect_equal(
    as.numeric(stats::weights(fit_cfd)),
    as.numeric(stats::weights(fit_energy)),
    tolerance = 1e-5
  )
})

# ---- Property tests: binary --------------------------------------------------

test_that("a binary ate normalizes each group to its size", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  expect_all(w, function(value) value >= 0)
  expect_all(w, function(value) value >= 1e-8)
  treated <- data$exposure == 1
  expect_equal(sum(w[treated]), sum(treated), tolerance = 1e-4)
  expect_equal(sum(w[!treated]), sum(!treated), tolerance = 1e-4)
  expect_true(fit@converged)
  # Characteristic-function-distance balancing drives balance through its
  # objective rather than exact moment constraints, so the achieved first-moment
  # imbalance is verified against the conventional good-balance ceiling.
  expect_balanced(fit, data, tolerance = 0.1)
})

test_that("a binary att targets the treated total in both groups", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
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

test_that("a binary atc fit produces non-negative floored weights", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
    estimand = "atc"
  )
  w <- as.numeric(stats::weights(fit))
  expect_all(w, function(value) value >= 0)
  expect_all(w, function(value) value >= 1e-8)
  expect_balanced(fit, data, tolerance = 0.1)
})

# The kernel is built from the same covariate matrix energy balancing forms, so
# a factor covariate reaches it as one indicator column per level. The expansion
# is the only place the levels are read, and a fit that lost them would still
# normalize and still balance the numeric covariates, so the factor's own
# imbalance is measured directly: the three level proportions differ between the
# exposure groups by up to 0.19 unweighted.

test_that("kernel balancing balances a factor covariate", {
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
      method = bw_cfd(),
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

# Kish effective sample size, computed inline since balance assessment moved to
# halfmoon. Bounds alone are no test of it: every strictly positive weight
# vector sits between zero and its own group size by Cauchy-Schwarz, so the
# specs below pin where the figure lands, what moves it, and the one case where
# it is exact.
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
    method = bw_cfd(),
    estimand = "ate"
  )
  w <- as.numeric(weights(fit))
  for (idx in exposure_groups(data)) {
    group_ess <- kish_ess(w[idx])
    # Balancing this confounded exposure costs precision, so the figure is
    # strictly below the group size by a wide margin rather than merely bounded
    # by it: uniform weights, which balance nothing, would sit at the size
    # exactly. It stays well clear of the floor a handful of dominating weights
    # would leave, which is the other way a fit can fail while still reporting
    # positive weights.
    expect_lt(group_ess, 0.9 * length(idx))
    expect_gt(group_ess, 0.1 * length(idx))
  }
})

test_that("the per-group effective sample size rises with the weight penalty", {
  # The weight penalty is the L2 term that pulls the solution toward the base
  # weights, so raising it buys precision back from balance. Tying the figure to
  # the knob that moves it is what distinguishes these weights from any other
  # positive vector: an arbitrary one has no reason to be ordered this way in
  # both groups at once.
  data <- sim_binary()
  ess_at <- function(penalty) {
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cfd(weight_penalty = penalty),
      estimand = "ate"
    )
    w <- as.numeric(weights(fit))
    vapply(exposure_groups(data), function(idx) kish_ess(w[idx]), numeric(1))
  }
  penalties <- c(1e-4, 1e-2, 1e-1)
  curve <- lapply(penalties, ess_at)

  for (step in seq_len(length(curve) - 1L)) {
    expect_all(curve[[step + 1L]], function(value) value > curve[[step]])
  }
})

test_that("a focal group keeps its base weights, so its ESS is its size", {
  # The focal group of a treated estimand is held at its base weights and
  # renormalized to its own total, which leaves every one of its weights equal
  # to one. Its effective sample size is therefore its size exactly, an equality
  # the reweighted group cannot meet: that group carries the whole tilt and
  # comes back at a fraction of its own size.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
    estimand = "att"
  )
  w <- as.numeric(weights(fit))
  groups <- exposure_groups(data)
  focal <- groups[[fit@focal_level]]
  reweighted <- groups[[setdiff(names(groups), fit@focal_level)]]

  expect_equal(kish_ess(w[focal]), length(focal), tolerance = 1e-8)
  expect_lt(kish_ess(w[reweighted]), 0.5 * length(reweighted))
  expect_gt(kish_ess(w[reweighted]), 0.05 * length(reweighted))
})

# ---- Categorical ----------------------------------------------------------

test_that("categorical ate kernel balancing produces valid weights", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
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

test_that("categorical att kernel balancing produces valid weights", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
    estimand = "att",
    .focal_level = "b"
  )
  w <- as.numeric(stats::weights(fit))
  expect_all(w, function(value) value >= 0)
  expect_balanced(fit, data, tolerance = 0.1)
})

# ---- Kernels --------------------------------------------------------------

test_that("every kernel converges to valid, floored weights", {
  data <- sim_binary()
  fit_of <- function(method) {
    balance(data, exposure, c(x1, x2), method = method, estimand = "ate")
  }
  methods <- list(
    gaussian = bw_cfd(kernel = "gaussian"),
    matern05 = bw_cfd(kernel = "matern", smoothness = 0.5),
    matern15 = bw_cfd(kernel = "matern", smoothness = 1.5),
    matern25 = bw_cfd(kernel = "matern", smoothness = 2.5),
    laplace = bw_cfd(kernel = "laplace")
  )
  for (method in methods) {
    fit <- fit_of(method)
    w <- as.numeric(stats::weights(fit))
    expect_true(fit@converged)
    expect_all(w, function(value) value >= 0)
    expect_all(w, function(value) value >= 1e-8)
  }
})

test_that("distinct kernels produce distinct weight sets", {
  data <- sim_binary()
  weights_for <- function(method) {
    normalize_by_group(
      as.numeric(stats::weights(
        balance(data, exposure, c(x1, x2), method = method, estimand = "ate")
      )),
      data$exposure
    )
  }
  w_gaussian <- weights_for(bw_cfd(kernel = "gaussian"))
  w_laplace <- weights_for(bw_cfd(kernel = "laplace"))
  w_energy <- weights_for(bw_cfd(kernel = "energy"))
  expect_true(weights_differ(w_gaussian, w_laplace))
  expect_true(weights_differ(w_gaussian, w_energy))
  expect_true(weights_differ(w_laplace, w_energy))
})

test_that("the Matern smoothness values produce distinct weight sets", {
  data <- sim_binary()
  weights_for <- function(smoothness) {
    normalize_by_group(
      as.numeric(stats::weights(
        balance(
          data,
          exposure,
          c(x1, x2),
          method = bw_cfd(kernel = "matern", smoothness = smoothness),
          estimand = "ate"
        )
      )),
      data$exposure
    )
  }
  w_05 <- weights_for(0.5)
  w_15 <- weights_for(1.5)
  w_25 <- weights_for(2.5)
  expect_true(weights_differ(w_05, w_15))
  expect_true(weights_differ(w_15, w_25))
  expect_true(weights_differ(w_05, w_25))
})

test_that("a distance-based kernel is invariant to a uniform covariate rescaling", {
  # The bandwidth is the median of the pairwise distances times a scale factor,
  # so multiplying the whole covariate matrix by a constant scales every distance
  # and the bandwidth together and leaves the kernel matrix, and therefore the
  # weights, unchanged.
  data <- sim_binary()
  scaled <- data
  scaled$x1 <- data$x1 * 100
  scaled$x2 <- data$x2 * 100
  # Both solves are tightened well past the 1e-8 solver default so that each one
  # settles inside its own tolerance ball rather than wherever its iteration path
  # happened to stop. At the default the two agree only to about 4e-5, the same
  # order as the assertion below, so the comparison was reading solver
  # reproducibility rather than the kernel identity it is about. At 2e-11 they
  # agree to about 5e-7. Tighter is possible but not by much: this problem stops
  # converging below roughly 1e-11, and a fit that gives up warns.
  fit <- function(df) {
    balance(
      df,
      exposure,
      c(x1, x2),
      method = bw_cfd(kernel = "gaussian", convergence_tolerance = 2e-11),
      estimand = "ate"
    )
  }
  w_raw <- normalize_by_group(
    as.numeric(stats::weights(fit(data))),
    data$exposure
  )
  w_scaled <- normalize_by_group(
    as.numeric(stats::weights(fit(scaled))),
    scaled$exposure
  )
  expect_equal(w_raw, w_scaled, tolerance = 1e-4)
})

# ---- t kernel reproducibility ---------------------------------------------

test_that("t-kernel weights are reproducible under a fixed seed", {
  # The t-kernel Monte Carlo projections are drawn on the R side under R's RNG,
  # so a fixed seed reproduces the projection matrix and therefore the weights.
  data <- sim_binary()
  fit_once <- function() {
    withr::with_seed(
      99,
      as.numeric(stats::weights(
        balance(
          data,
          exposure,
          c(x1, x2),
          method = bw_cfd(kernel = "t"),
          estimand = "ate"
        )
      ))
    )
  }
  expect_equal(fit_once(), fit_once())
})

test_that("t-kernel weights respond to the number of simulation draws", {
  data <- sim_binary()
  weights_with <- function(draws) {
    withr::with_seed(
      99,
      as.numeric(stats::weights(
        balance(
          data,
          exposure,
          c(x1, x2),
          method = bw_cfd(kernel = "t", simulation_draws = draws),
          estimand = "ate"
        )
      ))
    )
  }
  expect_true(weights_differ(
    normalize_by_group(weights_with(500), data$exposure),
    normalize_by_group(weights_with(5000), data$exposure)
  ))
})

# ---- Improved variant -----------------------------------------------------

test_that("the improved variant differs from the plain variant and both converge", {
  data <- sim_binary()
  fit_of <- function(improved) {
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cfd(improved = improved),
      estimand = "ate"
    )
  }
  fit_improved <- fit_of(TRUE)
  fit_plain <- fit_of(FALSE)
  expect_true(fit_improved@converged)
  expect_true(fit_plain@converged)
  w_improved <- normalize_by_group(
    as.numeric(stats::weights(fit_improved)),
    data$exposure
  )
  w_plain <- normalize_by_group(
    as.numeric(stats::weights(fit_plain)),
    data$exposure
  )
  expect_true(weights_differ(w_improved, w_plain))
})

# ---- Sampling weights -----------------------------------------------------

test_that("a binary ate balances under non-uniform sampling weights", {
  data <- sim_binary()
  sw <- withr::with_seed(11, stats::runif(nrow(data), 0.3, 3))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
    estimand = "ate",
    sampling_weights = sw
  )
  # The reported weights compose the sampling weights, so each group's
  # sampling-weighted total returns to the group's sampling-weight sum.
  w <- as.numeric(stats::weights(fit))
  treated <- data$exposure == 1
  expect_equal(sum(w[treated]), sum(sw[treated]), tolerance = 1e-3)
  expect_equal(sum(w[!treated]), sum(sw[!treated]), tolerance = 1e-3)
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
    method = bw_cfd(),
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
    method = bw_cfd(),
    estimand = "ate",
    constraints = balance_terms(moments = 1L, tolerance = 0.1)
  )
  expect_balanced(fit, data, tolerance = 0.1)
})

test_that("a tolerance without moment constraints warns and is ignored", {
  # The tolerance relaxes added moment constraints; with none present it has
  # nothing to act on, so kernel balancing warns and proceeds. The reference fit
  # without the tolerance produces the same weights.
  data <- sim_binary()
  expect_warning(
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cfd(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0.1)
    ),
    class = "balancing_warning"
  )
  reference <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
    estimand = "ate"
  )
  expect_equal(
    as.numeric(stats::weights(fit)),
    as.numeric(stats::weights(reference)),
    tolerance = 1e-6
  )
})

# ---- Diagnostics ----------------------------------------------------------

test_that("the fit reports dual variables and the quadratic-program backend", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
    estimand = "ate"
  )
  expect_s3_class(fit@duals, "data.frame")
  expect_gt(nrow(fit@duals), 0)
  # The default backend is osqp; the solver status reflects it.
  expect_identical(fit@solver_status, "osqp")
})

test_that("a kernel balancing fit has no estimating equations", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
    estimand = "ate"
  )
  expect_null(fit@estimating_equations)
  expect_error(
    estimating_equations(fit),
    class = "balancing_ipw_unsupported_error"
  )
})

test_that("ipw() rejects a kernel balancing fit", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
    estimand = "ate"
  )
  outcome <- stats::lm(x1 ~ exposure, data = data)
  expect_error(
    ipw(fit, outcome),
    class = "balancing_ipw_unsupported_error"
  )
})

# ---- Backend routing ------------------------------------------------------

test_that("an explicit clarabel backend solves and records itself", {
  # The default kernel, like every kernel but energy, is positive semidefinite by
  # construction, so the interior-point backend is eligible and honors the option.
  withr::local_options(balancing.qp_backend = "clarabel")
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
    estimand = "ate"
  )
  expect_identical(fit@solver_status, "clarabel")
  expect_true(fit@converged)
})

# The energy kernel is only conditionally positive semidefinite, so the quadratic
# term it assembles is indefinite and the interior-point backend refuses it up
# front. The routing to the ADMM backend is therefore correct rather than
# optional, and what the fit owed the caller was to say so: a pinned clarabel
# request used to be dropped in silence while the recorded backend read osqp.
test_that("a clarabel pin the energy kernel cannot honor is announced", {
  withr::local_options(balancing.qp_backend = "clarabel")
  data <- sim_binary()
  expect_warning(
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cfd(kernel = "energy"),
      estimand = "ate"
    ),
    class = "balancing_ignored_argument_warning"
  )
  # The recorded backend names what actually ran, not what was asked for.
  expect_identical(fit@solver_status, "osqp")
  expect_true(fit@converged)
})

test_that("the energy kernel is silent under the automatic and osqp backends", {
  data <- sim_binary()
  for (backend in c("auto", "osqp")) {
    fit <- withr::with_options(
      list(balancing.qp_backend = backend),
      expect_no_warning(balance(
        data,
        exposure,
        c(x1, x2),
        method = bw_cfd(kernel = "energy"),
        estimand = "ate"
      ))
    )
    expect_identical(fit@solver_status, "osqp")
  }
})

test_that("a clarabel pin a positive-semidefinite kernel honors is silent", {
  data <- sim_binary()
  for (kernel in c("gaussian", "laplace", "matern")) {
    fit <- withr::with_options(
      list(balancing.qp_backend = "clarabel"),
      expect_no_warning(balance(
        data,
        exposure,
        c(x1, x2),
        method = bw_cfd(kernel = kernel),
        estimand = "ate"
      ))
    )
    expect_identical(fit@solver_status, "clarabel")
  }
})

# ---- Unsupported exposure type --------------------------------------------

test_that("a continuous exposure raises balancing_exposure_type_error", {
  # Characteristic function distance balancing supports discrete exposures only;
  # a continuous exposure is refused with the shared exposure-type condition.
  data <- sim_continuous()
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cfd(),
      estimand = "ate"
    ),
    class = "balancing_exposure_type_error"
  )
})

# ---- Unsupported estimands ------------------------------------------------

# The overlap estimand is legal only for the covariate balancing propensity
# score, so kernel balancing rejects it for every exposure type with the shared
# unsupported-estimand condition.

test_that("the ato estimand raises balancing_estimand_error for a binary exposure", {
  data <- sim_binary()
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cfd(),
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
      method = bw_cfd(),
      estimand = "ato"
    ),
    class = "balancing_estimand_error"
  )
})

# ---- Print snapshot -------------------------------------------------------

test_that("a kernel balancing fit prints its summary block", {
  # Records on the first successful run once the fit path exists.
  data <- sim_binary()
  expect_balancing_snapshot({
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cfd(),
      estimand = "ate"
    )
    fit
  })
})

# ---- Solver tolerance box -------------------------------------------------

# The moment-constraint band a kernel balancing fit hands the solver is built by
# `solver_box()` (R/method-entropy.R), which the CFD assembly calls directly at
# R/method-cfd.R with the prepared matrix, the requested tolerances, and the
# sampling weights. A tolerance is written on the standardized scale, so the box
# converts it by the column's standard deviation, and a column holding one value
# repeated has no spread to convert against: its box is the tolerance itself.
#
# It does not arrive that way on its own. The weighted center divides a sum of
# products by a sum of weights and need not return the repeated value exactly,
# so the centered column carries a rounding residual and the computed scale
# reports that residual as the column's spread. Left alone, the constant 0.98
# column below shrinks its own band by roughly fifteen orders of magnitude and
# the fit is constrained against rounding. This pins the guard on the call the
# CFD path makes rather than only on the entropy one.
cfd_solver_box_fixture <- function() {
  withr::with_seed(707, {
    n <- 300L
    z <- cbind(
      stats::rnorm(n),
      stats::runif(n, -2, 3),
      rep(0.98, n),
      as.numeric(stats::rbinom(n, 1L, 0.4))
    )
    list(z = z, sampling_weights = stats::runif(n, 0.3, 2.5))
  })
}

test_that("the kernel balancing tolerance box leaves a constant column raw", {
  fixture <- cfd_solver_box_fixture()
  z <- fixture$z
  w <- fixture$sampling_weights
  tolerances <- seq_len(ncol(z)) / 100
  constant <- 3L

  expect_all(z[, constant], function(value) value == 0.98)
  expect_identical(
    solver_box(z, tolerances, w)[[constant]],
    tolerances[[constant]]
  )
  expect_identical(
    solver_box(z, tolerances)[[constant]],
    tolerances[[constant]]
  )
})
