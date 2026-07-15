# bw_ipt() is the method spec; balance(..., method = bw_ipt()) fits it. Inverse
# probability tilting fits a propensity model whose tilted score equations force
# the weighted covariate means to their estimand targets, so the achieved
# balance is exact on the requested moments. These specs cover the constructor,
# the property validators, the capability methods, and the statistical promises:
# achieved balance through expect_balanced(), non-negative weights,
# estimand-correct group sums, ESS bounded by n, the stored link coefficients,
# and the estimating-equations container. The link set is logit, probit, and
# cloglog. IPT supports binary and categorical exposures only.
#
# The identity that entropy balancing and just-identified inverse probability
# tilting produce the same average-treatment-effect-on-the-treated weights for a
# binary exposure is tested directly against entropy balancing, which is already
# implemented.

# Normalize weights to mean one within each exposure group so comparisons test
# the weighting solution rather than a reporting convention.
normalize_by_group <- function(w, g) {
  ave <- tapply(w, g, mean)
  as.numeric(w / ave[as.character(g)])
}

# ---- Constructor ----------------------------------------------------------

test_that("bw_ipt() carries its documented defaults", {
  spec <- bw_ipt()
  expect_true(S7::S7_inherits(spec, bw_ipt))
  expect_true(S7::S7_inherits(spec, estimating_equation_method))
  expect_true(S7::S7_inherits(spec, balance_method))
  expect_identical(spec@link, "logit")
  expect_identical(spec@convergence_tolerance, 1e-10)
  expect_null(spec@max_iterations)
})

test_that("bw_ipt() stores supplied tuning parameters", {
  spec <- bw_ipt(
    link = "probit",
    convergence_tolerance = 1e-8,
    max_iterations = 200L
  )
  expect_identical(spec@link, "probit")
  expect_identical(spec@convergence_tolerance, 1e-8)
  expect_identical(spec@max_iterations, 200L)
})

test_that("bw_ipt() matches the link argument", {
  expect_identical(bw_ipt(link = "cloglog")@link, "cloglog")
  expect_error(bw_ipt(link = "identity"))
})

test_that("bw_ipt() rejects unnamed extra arguments", {
  expect_true(S7::S7_inherits(bw_ipt(), balance_method))
  expect_error(bw_ipt(bogus = 1), class = "balancing_method_error")
})

# ---- Validators -----------------------------------------------------------

test_that("bw_ipt() rejects a non-positive convergence tolerance", {
  expect_identical(bw_ipt()@convergence_tolerance, 1e-10)
  expect_error(bw_ipt(convergence_tolerance = -1e-10))
})

test_that("bw_ipt() rejects a negative iteration cap", {
  expect_null(bw_ipt()@max_iterations)
  expect_error(bw_ipt(max_iterations = -5L))
})

# ---- Capability methods ---------------------------------------------------

test_that("supported_exposure_types() excludes continuous", {
  expect_setequal(
    supported_exposure_types(bw_ipt()),
    c("binary", "categorical")
  )
  expect_false("continuous" %in% supported_exposure_types(bw_ipt()))
})

test_that("supported_estimands() depends on the exposure type", {
  binary <- supported_estimands(bw_ipt(), "binary")
  expect_true(all(c("ate", "att") %in% binary))
  expect_true(any(c("atc", "atu") %in% binary))
  expect_false("ato" %in% binary)

  expect_setequal(
    supported_estimands(bw_ipt(), "categorical"),
    c("ate", "att")
  )
})

test_that("supports_estimating_equations() is always TRUE for bw_ipt", {
  expect_true(supports_estimating_equations(bw_ipt()))
  expect_true(supports_estimating_equations(
    bw_ipt(),
    constraints = balance_terms(tolerance = 0)
  ))
  # Unlike entropy balancing, inverse probability tilting keeps its estimating
  # equations regardless of the requested tolerance.
  expect_true(supports_estimating_equations(
    bw_ipt(),
    constraints = balance_terms(tolerance = 0.05)
  ))
})

test_that("method_label() names the method", {
  expect_identical(method_label(bw_ipt()), "Inverse probability tilting")
})

# ---- Statistical promises: binary -----------------------------------------

test_that("bw_ipt balances a binary ate", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

test_that("bw_ipt balances a binary att", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "att"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

test_that("bw_ipt balances a binary atc", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "atc"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

# ---- Statistical promises: categorical ------------------------------------

test_that("bw_ipt balances a categorical ate", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

test_that("bw_ipt balances a categorical att", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "att",
    focal_level = "b"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

test_that("a capped iteration count warns about convergence", {
  # Three Newton steps drive the balance essentially to zero but do not meet the
  # gradient tolerance, so the fit returns a usable iterate and warns about
  # convergence alone rather than erroring.
  data <- sim_binary()
  expect_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_ipt(max_iterations = 3L),
      estimand = "ate"
    ),
    class = "balancing_convergence_warning"
  )
})

# ---- Group sums -----------------------------------------------------------

test_that("a binary ate normalizes each group to its size", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
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
    method = bw_ipt(),
    estimand = "att"
  )
  w <- as.numeric(stats::weights(fit))
  treated <- data$exposure == 1
  n_treated <- sum(treated)
  expect_equal(w[treated], rep(1, n_treated), tolerance = 1e-6)
  expect_equal(sum(w[!treated]), n_treated, tolerance = 1e-4)
})

# ---- Property tests under sampling weights --------------------------------

test_that("bw_ipt balances a binary ate under sampling weights", {
  data <- sim_binary()
  withr::local_seed(21)
  data$sw <- stats::runif(nrow(data), 0.3, 3)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate",
    sampling_weights = sw
  )
  # Balance holds against the sampling-weighted pooled reference, which
  # expect_balanced() derives from the fit's own sampling weights.
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))

  # weights() folds the sampling weights in, so each group's total matches its
  # sampling-weighted size, the mean-one-per-group convention under weights.
  w <- as.numeric(stats::weights(fit))
  treated <- data$exposure == 1
  expect_equal(sum(w[treated]), sum(data$sw[treated]), tolerance = 1e-4)
  expect_equal(sum(w[!treated]), sum(data$sw[!treated]), tolerance = 1e-4)
})

test_that("bw_ipt balances a binary att under sampling weights", {
  data <- sim_binary()
  withr::local_seed(22)
  data$sw <- stats::runif(nrow(data), 0.3, 3)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "att",
    sampling_weights = sw
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))

  # The focal (treated) units keep base weight one, so their reported weight is
  # the sampling weight alone; the control group matches the focal total.
  w <- as.numeric(stats::weights(fit))
  treated <- data$exposure == 1
  expect_equal(w[treated], data$sw[treated], tolerance = 1e-6)
  expect_equal(sum(w[!treated]), sum(data$sw[treated]), tolerance = 1e-4)
})

test_that("bw_ipt balances a categorical att under sampling weights", {
  data <- sim_categorical()
  withr::local_seed(23)
  data$sw <- stats::runif(nrow(data), 0.3, 3)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "att",
    focal_level = "b",
    sampling_weights = sw
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))

  # Each non-focal level is tilted to the focal total, and the focal level keeps
  # base weight one.
  w <- as.numeric(stats::weights(fit))
  focal <- data$exposure == "b"
  focal_total <- sum(data$sw[focal])
  expect_equal(w[focal], data$sw[focal], tolerance = 1e-6)
  for (level in c("a", "c")) {
    idx <- data$exposure == level
    expect_equal(sum(w[idx]), focal_total, tolerance = 1e-4)
  }
})

# ---- ESS ------------------------------------------------------------------

test_that("the effective sample size is bounded by n within each group", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
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

# ---- Stored link coefficients ---------------------------------------------

test_that("the fit stores the link coefficients", {
  # The result class stores the tilt's link coefficients; propensity scores
  # inform the weights but are not stored on the object.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  expect_false(is.null(fit@coefficients))
  expect_type(fit@coefficients, "double")
  expect_gt(length(fit@coefficients), 0)
  expect_true(all(is.finite(fit@coefficients)))
})

# ---- Estimating equations from the core -----------------------------------

test_that("the estimating equations are populated with consistent dimensions", {
  # Inverse probability tilting always supplies estimating equations. The
  # container carries the matrices the Rust core returns, so the per-unit
  # estimating functions must sum to zero at the solution.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  ee <- estimating_equations(fit)
  n <- nrow(data)
  p <- ncol(ee@psi)
  expect_equal(nrow(ee@psi), n)
  expect_equal(nrow(ee@jacobian), p)
  expect_equal(ncol(ee@jacobian), p)
  expect_equal(dim(ee@weight_jacobian), dim(ee@psi))
  expect_lt(max(abs(colSums(ee@psi))), 1e-6)
})

# ---- Link functions -------------------------------------------------------

test_that("each link function fits and balances a binary ate", {
  data <- sim_binary()
  for (link in c("logit", "probit", "cloglog")) {
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_ipt(link = link),
      estimand = "ate"
    )
    expect_balanced(fit, data)
    expect_true(all(stats::weights(fit) >= 0))
  }
})

test_that("the link changes the weights", {
  data <- sim_binary()
  fit_logit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(link = "logit"),
    estimand = "ate"
  )
  fit_probit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(link = "probit"),
    estimand = "ate"
  )
  w_logit <- as.numeric(stats::weights(fit_logit))
  w_probit <- as.numeric(stats::weights(fit_probit))
  expect_false(isTRUE(all.equal(w_logit, w_probit)))
})

# ---- Cross-method identity ------------------------------------------------

test_that("bw_ipt att weights equal entropy balancing att weights (binary)", {
  # With mean balance and the logit link, inverse probability tilting and
  # entropy balancing solve the same treated-target problem, so their
  # average-treatment-effect-on-the-treated weights agree.
  data <- sim_binary()
  fit_ipt <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "att"
  )
  fit_ebal <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att"
  )
  w_ipt <- normalize_by_group(
    as.numeric(stats::weights(fit_ipt)),
    data$exposure
  )
  w_ebal <- normalize_by_group(
    as.numeric(stats::weights(fit_ebal)),
    data$exposure
  )
  expect_equal(w_ipt, w_ebal, tolerance = 1e-6)
})

# ---- Unsupported exposure type and estimand -------------------------------

test_that("a continuous exposure raises balancing_exposure_type_error", {
  data <- sim_continuous()
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_ipt(),
      estimand = "ate"
    ),
    class = "balancing_exposure_type_error"
  )
})

test_that("the ato estimand raises balancing_estimand_error", {
  data <- sim_binary()
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_ipt(),
      estimand = "ato"
    ),
    class = "balancing_estimand_error"
  )
})

# ---- Live consistency against WeightIt ------------------------------------

test_that("bw_ipt weights match WeightIt for a binary ate", {
  skip_on_cran()
  skip_if_not_installed("WeightIt")
  # WeightIt's ipt method solves its moment equations through rootSolve, a
  # transitive backend balancing does not depend on directly.
  skip_if_not_installed("rootSolve")

  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  reference <- WeightIt::weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "ipt",
    estimand = "ATE"
  )

  ours <- normalize_by_group(
    as.numeric(stats::weights(fit)),
    data$exposure
  )
  theirs <- normalize_by_group(reference$weights, data$exposure)
  expect_equal(ours, theirs, tolerance = 1e-6)
})

test_that("bw_ipt weights match WeightIt for a binary att", {
  skip_on_cran()
  skip_if_not_installed("WeightIt")
  # WeightIt's ipt method solves its moment equations through rootSolve, a
  # transitive backend balancing does not depend on directly.
  skip_if_not_installed("rootSolve")

  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "att"
  )
  reference <- WeightIt::weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "ipt",
    estimand = "ATT",
    focal = "1"
  )

  ours <- normalize_by_group(
    as.numeric(stats::weights(fit)),
    data$exposure
  )
  theirs <- normalize_by_group(reference$weights, data$exposure)
  expect_equal(ours, theirs, tolerance = 1e-6)
})

# ---- Print snapshot -------------------------------------------------------

test_that("an bw_ipt fit prints its summary block", {
  # Records on the first successful run once the fit path exists.
  data <- sim_binary()
  expect_snapshot({
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_ipt(),
      estimand = "ate"
    )
    fit
  })
})
