# bw_entropy() is the method spec; balance(..., method = bw_entropy())
# fits it. These specs cover the constructor, the capability methods, and the
# statistical promises: achieved balance through expect_balanced(), non-negative
# weights, estimand-correct group sums, ESS bounded by n, the fixed points where
# the weights reduce to the base weights, and the inexact (tolerance) problem.
#
# Group-sum conventions follow entropy balancing's usual normalization: for the
# ate every exposure group's weights sum to its group size (mean weight one);
# for the att the treated group keeps its base weights and the control group is
# reweighted to the treated total.

# ---- Constructor ----------------------------------------------------------

test_that("bw_entropy() carries its documented defaults", {
  spec <- bw_entropy()
  expect_true(S7::S7_inherits(spec, bw_entropy))
  expect_true(S7::S7_inherits(spec, balance_method))
  expect_null(spec@base_weights)
  expect_null(spec@distribution_moments)
  expect_identical(spec@convergence_tolerance, 1e-10)
  expect_null(spec@max_iterations)
})

test_that("bw_entropy() stores supplied tuning parameters", {
  spec <- bw_entropy(
    base_weights = rep(1, 5),
    convergence_tolerance = 1e-8,
    max_iterations = 200L
  )
  expect_equal(spec@base_weights, rep(1, 5))
  expect_identical(spec@convergence_tolerance, 1e-8)
  expect_identical(spec@max_iterations, 200L)
})

test_that("bw_entropy() rejects unnamed extra arguments", {
  expect_true(S7::S7_inherits(bw_entropy(), balance_method))
  expect_error(bw_entropy(1e-8), class = "balancing_method_error")
  expect_error(bw_entropy(bogus = 1), class = "balancing_method_error")
})

# ---- Capability methods ---------------------------------------------------

test_that("supported_exposure_types() lists every exposure type", {
  expect_setequal(
    supported_exposure_types(bw_entropy()),
    c("binary", "categorical", "continuous")
  )
})

test_that("supported_estimands() depends on the exposure type", {
  binary <- supported_estimands(bw_entropy(), "binary")
  expect_true(all(c("ate", "att") %in% binary))
  expect_true(any(c("atc", "atu") %in% binary))
  expect_false("ato" %in% binary)

  expect_setequal(
    supported_estimands(bw_entropy(), "categorical"),
    c("ate", "att")
  )
  expect_setequal(
    supported_estimands(bw_entropy(), "continuous"),
    "ate"
  )
})

test_that("supports_estimating_equations() follows the tolerance rule", {
  expect_true(supports_estimating_equations(bw_entropy()))
  expect_true(supports_estimating_equations(
    bw_entropy(),
    constraints = balance_terms(tolerance = 0)
  ))
  expect_false(supports_estimating_equations(
    bw_entropy(),
    constraints = balance_terms(tolerance = 0.05)
  ))
})

test_that("method_label() names the method", {
  expect_identical(method_label(bw_entropy()), "Entropy balancing")
})

# ---- Statistical promises: binary -----------------------------------------

test_that("entropy balancing balances a binary ate", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

test_that("entropy balancing balances a binary att", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

test_that("entropy balancing balances a binary atc", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "atc"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

test_that("entropy balancing balances a factor covariate for a binary ate", {
  # A factor covariate crosses the boundary as one raw indicator column per
  # level, whose balance target is the pooled level proportion rather than the
  # zero a standardized numeric column carries. Every exposure group's weighted
  # level proportions must therefore land on the pooled proportions.
  #
  # The level indicators of a factor sum to a constant column, so the entropy
  # dual is exactly flat along that direction and its Hessian is singular. The
  # damped Newton solve reaches the dual's numerical optimum without reaching a
  # 1e-10 gradient, so the gradient tolerance is set to the scale this constraint
  # set can reach; the achieved balance below is exact well inside the assertion.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, x3),
    method = bw_entropy(convergence_tolerance = 1e-8),
    estimand = "ate"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))

  w <- as.numeric(stats::weights(fit))
  for (level in levels(data$x3)) {
    indicator <- as.numeric(data$x3 == level)
    for (group in unique(data$exposure)) {
      idx <- data$exposure == group
      expect_equal(
        stats::weighted.mean(indicator[idx], w[idx]),
        mean(indicator),
        tolerance = 1e-6
      )
    }
  }
})

test_that("entropy balancing balances a factor covariate for a binary att", {
  # For the treated target the control group's weighted level proportions match
  # the treated group's own proportions.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, x3),
    method = bw_entropy(),
    estimand = "att"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))

  w <- as.numeric(stats::weights(fit))
  treated <- data$exposure == 1
  for (level in levels(data$x3)) {
    indicator <- as.numeric(data$x3 == level)
    expect_equal(
      stats::weighted.mean(indicator[!treated], w[!treated]),
      mean(indicator[treated]),
      tolerance = 1e-6
    )
  }
})

test_that("a binary ate normalizes each group to its size", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  treated <- data$exposure == 1
  expect_equal(sum(w[treated]), sum(treated), tolerance = 1e-4)
  expect_equal(sum(w[!treated]), sum(!treated), tolerance = 1e-4)
})

test_that("a binary att keeps treated base weights and matches the control sum", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att"
  )
  w <- as.numeric(stats::weights(fit))
  treated <- data$exposure == 1
  n_treated <- sum(treated)
  expect_equal(w[treated], rep(1, n_treated), tolerance = 1e-6)
  expect_equal(sum(w[!treated]), n_treated, tolerance = 1e-4)
})

# ---- Statistical promises: categorical ------------------------------------

test_that("entropy balancing balances a categorical ate", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

test_that("entropy balancing balances a categorical att", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att",
    focal_level = "b"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

# ---- Statistical promises: continuous -------------------------------------

test_that("entropy balancing balances a continuous ate", {
  data <- sim_continuous()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

test_that("a continuous fit preserves an indicator covariate's marginal", {
  # An indicator covariate crosses the boundary as a raw zero/one column, so the
  # marginal it is held to is the sample proportion. A marginal target of zero
  # would instead drive every unit in the indicated stratum to no weight at all,
  # which balances nothing and empties the stratum.
  data <- sim_continuous_indicator()
  fit <- balance(
    data,
    exposure,
    c(x1, g),
    method = bw_entropy(),
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
  expect_balanced(fit, data)
  expect_true(all(w >= 0))
})

test_that("a continuous fit holds the base-measure marginals under sampling weights", {
  # Under informative sampling the reference is the sampling-weighted sample, so
  # the marginals the fit preserves are the sampling-weighted ones: the
  # indicator's proportion and the exposure mean alike.
  data <- sim_continuous_indicator()
  withr::local_seed(31)
  data$sw <- stats::runif(nrow(data), 0.3, 3)
  fit <- balance(
    data,
    exposure,
    c(x1, g),
    method = bw_entropy(),
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
  expect_balanced(fit, data)
})

test_that("a continuous fit returns the base weights when balance already holds", {
  # The exposure distribution is identical in the two covariate strata and the
  # base weights depend on the stratum alone, so the base measure already
  # decorrelates the exposure from the covariate. The marginal targets are taken
  # under that same base measure, which leaves the base weights themselves as
  # the solution. Targets read off the sample instead would pull the marginals
  # away from the base measure and tilt a solution that needs no tilt.
  grid <- seq(-2, 2, length.out = 50)
  data <- data.frame(
    exposure = rep(grid, 2),
    x1 = rep(c(-1, 1), each = length(grid))
  )
  base <- rep(c(1, 3), each = length(grid))
  fit <- balance(
    data,
    exposure,
    x1,
    method = bw_entropy(base_weights = base),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  expect_equal(w, base * nrow(data) / sum(base), tolerance = 1e-6)
})

# ---- ESS ------------------------------------------------------------------

test_that("the effective sample size is bounded by n within each group", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  # Kish effective sample size computed inline within each exposure group, since
  # balance assessment moved to halfmoon; each group's figure stays positive and
  # bounded by that group's size.
  w <- as.numeric(weights(fit))
  groups <- attr(fit@weights, "groups")
  for (idx in groups) {
    group_ess <- sum(w[idx])^2 / sum(w[idx]^2)
    expect_gt(group_ess, 0)
    expect_lte(group_ess, length(idx) + 1e-8)
  }
})

# ---- Fixed points ---------------------------------------------------------

test_that("weights reduce to the base weights when balance already holds", {
  # Treated and control share an identical covariate distribution, so no
  # reweighting is required and the entropy solution is the base weights.
  data <- data.frame(
    exposure = c(0L, 0L, 1L, 1L),
    x1 = c(-1, 1, -1, 1)
  )
  fit <- balance(
    data,
    exposure,
    x1,
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  expect_equal(w, rep(1, 4), tolerance = 1e-6)
})

test_that("base weights anchor the solution when balance already holds", {
  data <- data.frame(
    exposure = c(0L, 0L, 1L, 1L),
    x1 = c(-1, 1, -1, 1)
  )
  base <- c(2, 1, 2, 1)
  fit <- balance(
    data,
    exposure,
    x1,
    method = bw_entropy(base_weights = base),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  # With balance already satisfied, the KL solution is proportional to the
  # base weights within each exposure group.
  expect_equal(w[1] / w[2], base[1] / base[2], tolerance = 1e-6)
  expect_equal(w[3] / w[4], base[3] / base[4], tolerance = 1e-6)
})

# ---- Inexact (tolerance) problem ------------------------------------------

test_that("a positive tolerance balances within tolerance", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("the inexact tolerance box binds on the weighted scale under informative sampling weights", {
  # The inexact entropy path carries a tolerance box like the quadratic-program
  # family, so its box must bind on the sampling-weighted scale the balance table
  # reports on. Sampling weights correlated with x1 shrink its weighted spread
  # below its unweighted one; a box measured on the unweighted scale would bind at
  # the wrong width and warn against its own converged, in-box solution.
  data <- sim_binary()
  sw <- 0.3 + 2 * (abs(data$x1) < 0.5)
  fit <- expect_no_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0.01),
      sampling_weights = sw
    ),
    class = "balancing_balance_warning"
  )
  expect_true(all(fit@balance_table$within_tolerance))
  expect_balanced(fit, data, tolerance = 0.01)
})

test_that("a positive tolerance disables the estimating equations", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  expect_null(fit@estimating_equations)
})

test_that("the exact problem produces estimating equations", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att"
  )
  expect_false(is.null(fit@estimating_equations))
})

# ---- Estimating equations from the core -----------------------------------

test_that("the estimating equations satisfy the moment-condition identity", {
  # The container is built from the matrices the Rust core returns, so the
  # per-unit estimating functions must sum to zero at the solution: the weighted
  # constraint mean equals the target in every solved group. This verifies the
  # boundary rather than trusting it.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  ee <- estimating_equations(fit)
  expect_lt(max(abs(colSums(ee@psi))), 1e-6)
  # The Jacobian is the negative weighted covariance of the constraints, which
  # is symmetric.
  expect_equal(ee@jacobian, t(ee@jacobian), tolerance = 1e-8)
  # The container has one parameter block per solved group.
  n_groups <- length(unique(data$exposure))
  expect_identical(ncol(ee@psi), 2L * n_groups)
  expect_identical(dim(ee@weight_jacobian), dim(ee@psi))
})

test_that("an att fit reports estimating equations only for the control group", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att"
  )
  ee <- estimating_equations(fit)
  # Only the reweighted control group carries parameters (two covariates).
  expect_identical(ncol(ee@psi), 2L)
  expect_lt(max(abs(colSums(ee@psi))), 1e-6)
})

# ---- Solver option routing ------------------------------------------------

test_that("the entropy solver option routes the solver and holds the solution", {
  # The exact-problem solver is the promotion knob: it must route to the core
  # solver named in the option and reach the same solution as the Newton
  # default, which the quasi-Newton alternatives do to well within tolerance.
  data <- sim_binary()
  fit_newton <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_identical(fit_newton@solver_status, "newton")
  w_newton <- as.numeric(stats::weights(fit_newton))

  for (solver in c("lbfgs", "lbfgs_then_newton")) {
    fit <- withr::with_options(
      list(balancing.entropy_solver = solver),
      balance(
        data,
        exposure,
        c(x1, x2),
        method = bw_entropy(),
        estimand = "ate"
      )
    )
    expect_identical(fit@solver_status, solver)
    expect_equal(
      as.numeric(stats::weights(fit)),
      w_newton,
      tolerance = 1e-6
    )
  }
})

# ---- distribution_moments (continuous) ------------------------------------

test_that("distribution_moments holds the exposure variance", {
  data <- sim_continuous()
  fit_default <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  fit_moments <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(distribution_moments = 2L),
    estimand = "ate"
  )

  w_default <- as.numeric(stats::weights(fit_default))
  w_moments <- as.numeric(stats::weights(fit_moments))
  # The parameter is no longer inert: raising the distribution moments changes
  # the weights.
  expect_false(isTRUE(all.equal(w_default, w_moments)))

  sample_variance <- sum((data$exposure - mean(data$exposure))^2) /
    length(data$exposure)
  weighted_variance <- function(w) {
    center <- stats::weighted.mean(data$exposure, w)
    sum(w * (data$exposure - center)^2) / sum(w)
  }
  # The second-moment fit holds the weighted exposure variance at the sample
  # value; the first-moment default lets it drift.
  expect_equal(weighted_variance(w_moments), sample_variance, tolerance = 1e-4)
  expect_gt(
    abs(weighted_variance(w_default) - sample_variance),
    0.1 * sample_variance
  )
})

test_that("distribution_moments is raised to the constraint moments with an alert", {
  withr::local_options(balancing.quiet = FALSE)
  data <- sim_continuous()
  # With alerts on, the fit also announces the detected exposure type; capture
  # that outer alert so it does not leak into the test console.
  expect_message(
    expect_message(
      balance(
        data,
        exposure,
        c(x1, x2),
        method = bw_entropy(distribution_moments = 1L),
        estimand = "ate",
        constraints = balance_terms(moments = 2L)
      ),
      regexp = "distribution_moments"
    ),
    regexp = "continuous"
  )
})

# ---- Continuous with tolerance --------------------------------------------

test_that("a continuous tolerance relaxes correlations but holds the marginals", {
  data <- sim_continuous()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(distribution_moments = 2L),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.1)
  )
  # The inexact continuous path runs FISTA and produces no estimating equations.
  expect_identical(fit@solver_status, "fista")
  expect_null(fit@estimating_equations)
  # With the marginal variances held, the tolerance binds on the correlation
  # scale the design specifies, up to the small structural gap between the
  # constrained weighted association and the reported Pearson correlation.
  expect_lte(max(fit@balance_table$weighted), 0.105)
  expect_true(all(fit@balance_table$within_tolerance))

  # The exposure marginals stay at the unweighted sample: the mean is preserved.
  w <- as.numeric(stats::weights(fit))
  expect_equal(
    stats::weighted.mean(data$exposure, w),
    mean(data$exposure),
    tolerance = 1e-6
  )
})

# ---- Property tests under sampling and base weights -----------------------

test_that("entropy balancing balances a binary ate under sampling weights", {
  data <- sim_binary()
  withr::local_seed(11)
  data$sw <- stats::runif(nrow(data), 0.3, 3)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate",
    sampling_weights = sw
  )
  # Balance holds against the sampling-weighted pooled reference, which
  # expect_balanced() derives from the fit's own sampling weights.
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

test_that("entropy balancing balances a binary att under sampling weights", {
  data <- sim_binary()
  withr::local_seed(12)
  data$sw <- stats::runif(nrow(data), 0.3, 3)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att",
    sampling_weights = sw
  )
  expect_balanced(fit, data)
})

test_that("entropy balancing balances a binary ate under base weights", {
  data <- sim_binary()
  withr::local_seed(13)
  base <- stats::runif(nrow(data), 0.5, 2)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(base_weights = base),
    estimand = "ate"
  )
  # The base measure moves the pooled target; expect_balanced() reads the base
  # weights from the fitted method.
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

test_that("entropy balancing balances a binary atu under base weights", {
  data <- sim_binary()
  withr::local_seed(14)
  base <- stats::runif(nrow(data), 0.5, 2)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(base_weights = base),
    estimand = "atc"
  )
  expect_balanced(fit, data)
})

# ---- Live consistency against WeightIt ------------------------------------

test_that("entropy weights match WeightIt for a binary ate", {
  skip_on_cran()
  skip_if_not_installed("WeightIt")

  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  reference <- WeightIt::weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "ebal",
    estimand = "ATE"
  )

  ours <- as.numeric(stats::weights(fit))
  theirs <- reference$weights
  # Compare on a common normalization (mean one within each exposure group) so
  # only the weighting solution is under test, not the reporting convention.
  normalize <- function(w, g) {
    ave <- tapply(w, g, mean)
    # Return a bare numeric vector: subsetting the tapply() array by name would
    # otherwise leak dim and dimnames into the normalized weights.
    as.numeric(w / ave[as.character(g)])
  }
  ours_n <- normalize(ours, data$exposure)
  theirs_n <- normalize(theirs, data$exposure)
  expect_equal(ours_n, unname(theirs_n), tolerance = 1e-6)
})

test_that("entropy weights match WeightIt for a binary att", {
  skip_on_cran()
  skip_if_not_installed("WeightIt")

  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att"
  )
  reference <- WeightIt::weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "ebal",
    estimand = "ATT",
    focal = "1"
  )

  ours <- as.numeric(stats::weights(fit))
  theirs <- reference$weights
  normalize <- function(w, g) {
    ave <- tapply(w, g, mean)
    # Return a bare numeric vector: subsetting the tapply() array by name would
    # otherwise leak dim and dimnames into the normalized weights.
    as.numeric(w / ave[as.character(g)])
  }
  ours_n <- normalize(ours, data$exposure)
  theirs_n <- normalize(theirs, data$exposure)
  expect_equal(ours_n, unname(theirs_n), tolerance = 1e-6)
})
