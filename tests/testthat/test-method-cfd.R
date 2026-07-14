# bw_cfd() is the method spec; balance(..., method = bw_cfd()) fits it.
# Characteristic function distance balancing (kernel balancing) chooses weights
# that minimize a kernel measure of the distance between the reweighted exposure
# groups and the target sample. Each kernel is positive semidefinite by
# construction, so the objective is a convex quadratic program with a
# simplex-type constraint set, and the method spec carries only tuning
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

test_that("bw_cfd() rejects a negative weight penalty", {
  expect_identical(bw_cfd(weight_penalty = 1e-3)@weight_penalty, 1e-3)
  expect_error(bw_cfd(weight_penalty = -1e-4))
})

test_that("bw_cfd() rejects a negative minimum weight", {
  expect_identical(bw_cfd(min_weight = 1e-6)@min_weight, 1e-6)
  expect_error(bw_cfd(min_weight = -1e-8))
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
  expect_true(all(w >= 0))
  expect_true(all(w >= 1e-8))
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
  expect_true(all(w >= 0))
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
  expect_true(all(w >= 0))
  expect_true(all(w >= 1e-8))
  expect_balanced(fit, data, tolerance = 0.1)
})

# ---- ESS ------------------------------------------------------------------

test_that("the effective sample size is bounded by n within each group", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
    estimand = "ate"
  )
  ess_tbl <- ess(fit)
  expect_true(all(ess_tbl$ess <= ess_tbl$n + 1e-8))
  expect_true(all(ess_tbl$ess > 0))
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
  expect_true(all(w >= 0))
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
    focal_level = "b"
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
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
    expect_true(all(w >= 0))
    expect_true(all(w >= 1e-8))
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
  fit <- function(df) {
    balance(
      df,
      exposure,
      c(x1, x2),
      method = bw_cfd(kernel = "gaussian"),
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
    propensity::ipw(fit, outcome),
    class = "balancing_ipw_unsupported_error"
  )
})

# ---- Backend routing ------------------------------------------------------

test_that("an explicit clarabel backend solves and records itself", {
  # The kernels are positive semidefinite by construction, so the interior-point
  # backend is eligible and honors the option.
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
  expect_snapshot({
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
