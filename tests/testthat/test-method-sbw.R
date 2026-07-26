# bw_sbw() is the method spec; balance(..., method = bw_sbw()) fits it. Stable
# balancing weights (Zubizarreta) minimize the dispersion of the weights subject
# to approximate covariate balance: the weighted covariate means are held within
# a tolerance band rather than balanced exactly, and among all weightings that
# meet that band the method selects the one of least variance. The default "l2"
# norm minimizes the sum of squared weights, so for fixed per-group totals it
# minimizes the weight variance.
#
# Stable balancing weights belong to the quadratic-program family, which carries
# a different set of promises from the estimating-equation family. The balance
# tolerance is not an optional relaxation here but the method's central tuning
# parameter: with an exact (zero) tolerance the problem reduces to exact moment
# balance, which defeats the minimum-variance rationale and is infeasible-prone,
# so the design requires a positive tolerance and the fit refuses a fit without
# one. Given a positive tolerance the central guarantee is that the achieved
# standardized mean differences sit inside the band and that no feasible
# reweighting has smaller weight dispersion. The family has no estimating
# equations, so estimating_equations() and ipw() raise the shared unsupported
# condition.
#
# These specs therefore verify the minimum-variance-under-constraints
# characterization directly (achieved balance within tolerance through the shared
# expect_balanced() oracle, and monotone weight dispersion as the tolerance
# tightens, since a smaller feasible set cannot lower the minimized dispersion),
# check the estimand-correct group sums and the minimum-weight floor, and confirm
# the capability and diagnostic contracts. The live comparison against optweight
# is objective-level: both implementations solve the same strictly convex program,
# so ours must attain weight dispersion at or below theirs on the same tolerance
# box, the contract the design specifies for the quadratic-program family.

# ---- Oracles --------------------------------------------------------------

# The weight dispersion the default "l2" objective minimizes, up to an additive
# constant. Within each exposure group the weights are normalized to mean one and
# the summed squared deviation from one is accumulated. Because each group is
# reported at a fixed total, this equals the minimized sum of squared weights
# shifted by a constant, so it is monotone in the true objective and can order two
# fits on the same constraint set.
weight_dispersion <- function(w, g) {
  parts <- split(w, g)
  sum(vapply(
    parts,
    function(v) {
      v <- v / mean(v)
      sum((v - 1)^2)
    },
    numeric(1)
  ))
}

# Weights normalized to mean one within each exposure group so comparisons test
# the weighting solution rather than a per-group reporting convention.
normalize_by_group <- function(w, g) {
  ave <- tapply(w, g, mean)
  as.numeric(w / ave[as.character(g)])
}

# The weight dispersion the "l1" objective minimizes: the summed absolute
# departure from one within each group, normalized to mean one so it measures the
# weighting solution rather than a per-group total. Monotone in the true objective
# and comparable across two fits on the same constraint set.
l1_dispersion <- function(w, g) {
  parts <- split(w, g)
  sum(vapply(
    parts,
    function(v) {
      v <- v / mean(v)
      sum(abs(v - 1))
    },
    numeric(1)
  ))
}

# The weight dispersion the "linf" objective minimizes: the single largest
# absolute departure from one across every reweighted unit, normalized to mean one
# within each group, matching the solver's single shared deviation variable.
linf_dispersion <- function(w, g) {
  parts <- split(w, g)
  max(vapply(
    parts,
    function(v) {
      v <- v / mean(v)
      max(abs(v - 1))
    },
    numeric(1)
  ))
}

# The achieved absolute standardized mean difference of each constraint column,
# arm to pooled, for the average treatment effect. The constraint matrix is
# rebuilt from the recipe and standardized to the unweighted sample, the pooled
# target is the uniform-measure column mean, and each column reports the largest
# arm-to-target deviation. This resolves per-covariate tolerances, which the
# scalar expect_balanced() oracle cannot.
per_column_smd <- function(fit, .data) {
  constraint_matrix <- rebuild_constraint_matrix(fit@recipe, .data)
  centers <- colMeans(constraint_matrix)
  scales <- apply(constraint_matrix, 2, stats::sd)
  scales[scales == 0] <- 1
  z <- sweep(sweep(constraint_matrix, 2, centers, "-"), 2, scales, "/")

  w <- as.numeric(stats::weights(fit))
  exposure <- as.character(.data[[fit@exposure]])
  arms <- split(seq_along(exposure), exposure)
  target <- colMeans(z)

  achieved <- vapply(
    seq_len(ncol(z)),
    function(j) {
      max(vapply(
        arms,
        function(idx) abs(stats::weighted.mean(z[idx, j], w[idx]) - target[j]),
        numeric(1)
      ))
    },
    numeric(1)
  )
  stats::setNames(
    achieved,
    vapply(fit@recipe, function(r) r$term, character(1))
  )
}

# ---- Constructor ----------------------------------------------------------

test_that("bw_sbw() carries its documented defaults", {
  spec <- bw_sbw()
  expect_true(S7::S7_inherits(spec, bw_sbw))
  expect_true(S7::S7_inherits(spec, quadratic_program_method))
  expect_true(S7::S7_inherits(spec, balance_method))
  expect_identical(spec@norm, "l2")
  expect_identical(spec@min_weight, 1e-8)
  expect_null(spec@convergence_tolerance)
  expect_null(spec@max_iterations)
})

test_that("bw_sbw() stores supplied tuning parameters", {
  spec <- bw_sbw(
    norm = "l1",
    min_weight = 1e-6,
    convergence_tolerance = 1e-8,
    max_iterations = 500L
  )
  expect_identical(spec@norm, "l1")
  expect_identical(spec@min_weight, 1e-6)
  expect_identical(spec@convergence_tolerance, 1e-8)
  expect_identical(spec@max_iterations, 500L)
})

test_that("bw_sbw() matches the norm argument", {
  expect_identical(bw_sbw(norm = "l2")@norm, "l2")
  expect_identical(bw_sbw(norm = "l1")@norm, "l1")
  expect_identical(bw_sbw(norm = "linf")@norm, "linf")
  expect_error(bw_sbw(norm = "huber"))
})

test_that("bw_sbw() rejects unnamed and unknown extra arguments", {
  expect_true(S7::S7_inherits(bw_sbw(), balance_method))
  expect_error(bw_sbw(1e-8))
  expect_error(bw_sbw(bogus = 1), class = "balancing_method_error")
})

# ---- Validators -----------------------------------------------------------

test_that("bw_sbw() rejects a negative minimum weight", {
  expect_identical(bw_sbw(min_weight = 1e-6)@min_weight, 1e-6)
  expect_error(bw_sbw(min_weight = -1e-8))
})

test_that("bw_sbw() rejects a non-positive convergence tolerance", {
  expect_null(bw_sbw()@convergence_tolerance)
  expect_error(bw_sbw(convergence_tolerance = -1e-8))
})

test_that("bw_sbw() rejects a negative iteration cap", {
  expect_null(bw_sbw()@max_iterations)
  expect_error(bw_sbw(max_iterations = -5L))
})

# ---- Capability methods ---------------------------------------------------

test_that("supported_exposure_types() lists every exposure type", {
  expect_setequal(
    supported_exposure_types(bw_sbw()),
    c("binary", "categorical", "continuous")
  )
})

test_that("supported_estimands() depends on the exposure type", {
  binary <- supported_estimands(bw_sbw(), "binary")
  expect_true(all(c("ate", "att") %in% binary))
  expect_true(any(c("atc", "atu") %in% binary))
  expect_false("ato" %in% binary)

  expect_setequal(
    supported_estimands(bw_sbw(), "categorical"),
    c("ate", "att")
  )
  expect_setequal(
    supported_estimands(bw_sbw(), "continuous"),
    "ate"
  )
})

test_that("supports_estimating_equations() is always FALSE for bw_sbw", {
  # The quadratic-program family has no estimating equations, whatever the
  # constraints or exposure type.
  expect_false(supports_estimating_equations(bw_sbw()))
  expect_false(supports_estimating_equations(
    bw_sbw(),
    constraints = balance_terms(tolerance = 0.05)
  ))
  expect_false(supports_estimating_equations(
    bw_sbw(),
    exposure_type = "binary"
  ))
  expect_false(supports_estimating_equations(
    bw_sbw(),
    exposure_type = "continuous"
  ))
})

test_that("method_label() names the method", {
  expect_identical(method_label(bw_sbw()), "Stable balancing weights")
})

# ---- Property tests: minimum variance under constraints -------------------

test_that("a binary ate fit meets the tolerance and floors the weights", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  expect_true(all(w >= 1e-8))
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("a binary ate normalizes each group to its size", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
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
    method = bw_sbw(),
    estimand = "att",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  treated <- data$exposure == 1
  n_treated <- sum(treated)
  expect_equal(sum(w[treated]), n_treated, tolerance = 1e-4)
  expect_equal(sum(w[!treated]), n_treated, tolerance = 1e-4)
  expect_true(all(w >= 0))
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("a binary atc fit produces non-negative floored weights", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "atc",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  expect_true(all(w >= 1e-8))
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("tightening the tolerance cannot lower the weight dispersion", {
  # Stable balancing weights minimize the weight dispersion subject to the
  # tolerance band. A smaller tolerance is a smaller feasible set, so its
  # minimized dispersion is at least that of a larger tolerance. This is the
  # operational form of the minimum-variance-under-constraints characterization,
  # checked without an external reference.
  data <- sim_binary()
  fit_of <- function(tolerance) {
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_sbw(),
      estimand = "ate",
      constraints = balance_terms(tolerance = tolerance)
    )
  }
  tight <- fit_of(0.01)
  loose <- fit_of(0.1)
  expect_balanced(tight, data, tolerance = 0.01)
  expect_balanced(loose, data, tolerance = 0.1)

  disp_tight <- weight_dispersion(
    as.numeric(stats::weights(tight)),
    data$exposure
  )
  disp_loose <- weight_dispersion(
    as.numeric(stats::weights(loose)),
    data$exposure
  )
  expect_gte(disp_tight, disp_loose - 1e-6)
})

test_that("the minimum-weight floor holds on the reported scale", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(min_weight = 1e-3),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 1e-3))
})

# ---- Sampling weights -----------------------------------------------------

test_that("a binary ate balances under non-uniform sampling weights", {
  data <- sim_binary()
  sw <- withr::with_seed(11, stats::runif(nrow(data), 0.3, 3))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05),
    sampling_weights = sw
  )
  # The reported weights compose the sampling weights, so each group's
  # sampling-weighted total returns to the group's sampling-weight sum and the
  # balance oracle holds on that scale.
  expect_balanced(fit, data, tolerance = 0.05)
  w <- as.numeric(stats::weights(fit))
  treated <- data$exposure == 1
  expect_equal(sum(w[treated]), sum(sw[treated]), tolerance = 1e-3)
  expect_equal(sum(w[!treated]), sum(sw[!treated]), tolerance = 1e-3)
})

test_that("the tolerance box binds on the weighted scale under informative sampling weights", {
  # Sampling weights correlated with a covariate make its weighted standard
  # deviation differ from its unweighted one, so a box measured on the wrong scale
  # binds at the wrong width. The weights below shrink the weighted spread of x1
  # to well below its unweighted spread; the fit must still bind its box, and
  # report balance, on the weighted scale, so no balance warning fires on a
  # converged in-box solution.
  data <- sim_binary()
  sw <- 0.3 + 2 * (abs(data$x1) < 0.5)
  fit <- expect_no_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_sbw(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0.01),
      sampling_weights = sw
    ),
    class = "balancing_balance_warning"
  )
  expect_true(fit@converged)
  expect_true(all(fit@balance_table$within_tolerance))
  # The achieved balance binds at the tolerance on the weighted scale, not the
  # looser unweighted scale a mis-scaled box would have allowed.
  expect_balanced(fit, data, tolerance = 0.01)
})

test_that("a continuous ate meets the correlation tolerance under sampling weights", {
  data <- sim_continuous()
  sw <- withr::with_seed(12, stats::runif(nrow(data), 0.3, 3))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05),
    sampling_weights = sw
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  expect_equal(sum(w), sum(sw), tolerance = 1e-3)
  expect_balanced(fit, data, tolerance = 0.05)
})

# ---- ESS ------------------------------------------------------------------

test_that("the effective sample size is bounded by n within each group", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
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

# ---- Categorical ----------------------------------------------------------

test_that("categorical ate stable balancing produces valid weights", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  for (level in levels(data$exposure)) {
    idx <- data$exposure == level
    expect_equal(sum(w[idx]), sum(idx), tolerance = 1e-4)
  }
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("categorical att stable balancing produces valid weights", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "att",
    focal_level = "b",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  expect_balanced(fit, data, tolerance = 0.05)
})

# ---- Continuous -----------------------------------------------------------

test_that("continuous ate stable balancing meets the correlation tolerance", {
  data <- sim_continuous()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  expect_true(all(w >= 1e-8))
  expect_balanced(fit, data, tolerance = 0.05)
})

# ---- Tolerance semantics --------------------------------------------------

test_that("a fit without a positive tolerance raises balancing_constraints_error", {
  # The balance tolerance is the method's central tuning parameter. An exact
  # (zero) tolerance reduces stable balancing weights to exact moment balance,
  # which abandons the minimum-variance rationale and is infeasible-prone, so the
  # fit refuses to proceed and names the knob to set.
  data <- sim_binary()
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_sbw(),
      estimand = "ate"
    ),
    class = "balancing_constraints_error"
  )
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_sbw(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0)
    ),
    class = "balancing_constraints_error"
  )
})

test_that("the required-tolerance message names the tuning parameter", {
  data <- sim_binary()
  expect_balancing_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_sbw(),
      estimand = "ate"
    )
  )
})

# ---- Absolute-deviation norms ---------------------------------------------

test_that("an l1 binary ate meets the tolerance and normalizes each group", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(norm = "l1"),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  expect_true(all(w >= 1e-8))
  treated <- data$exposure == 1
  expect_equal(sum(w[treated]), sum(treated), tolerance = 1e-3)
  expect_equal(sum(w[!treated]), sum(!treated), tolerance = 1e-3)
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("an l1 binary att targets the treated total and meets the tolerance", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(norm = "l1"),
    estimand = "att",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  treated <- data$exposure == 1
  n_treated <- sum(treated)
  expect_equal(sum(w[treated]), n_treated, tolerance = 1e-3)
  expect_equal(sum(w[!treated]), n_treated, tolerance = 1e-3)
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("a linf binary ate meets the tolerance and normalizes each group", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(norm = "linf"),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  expect_true(all(w >= 1e-8))
  treated <- data$exposure == 1
  expect_equal(sum(w[treated]), sum(treated), tolerance = 1e-3)
  expect_equal(sum(w[!treated]), sum(!treated), tolerance = 1e-3)
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("a linf binary att meets the tolerance", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(norm = "linf"),
    estimand = "att",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("an l1 categorical ate produces valid balanced weights", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(norm = "l1"),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  for (level in levels(data$exposure)) {
    idx <- data$exposure == level
    expect_equal(sum(w[idx]), sum(idx), tolerance = 1e-3)
  }
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("a linf categorical ate produces valid balanced weights", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(norm = "linf"),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("an l1 continuous ate meets the correlation tolerance", {
  data <- sim_continuous()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(norm = "l1"),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  expect_true(all(w >= 1e-8))
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("a linf continuous ate meets the correlation tolerance", {
  data <- sim_continuous()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(norm = "linf"),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("an l1 fit balances under non-uniform sampling weights", {
  data <- sim_binary()
  sw <- withr::with_seed(11, stats::runif(nrow(data), 0.3, 3))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(norm = "l1"),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05),
    sampling_weights = sw
  )
  w <- as.numeric(stats::weights(fit))
  treated <- data$exposure == 1
  expect_equal(sum(w[treated]), sum(sw[treated]), tolerance = 1e-3)
  expect_equal(sum(w[!treated]), sum(sw[!treated]), tolerance = 1e-3)
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("a linf fit balances under non-uniform sampling weights", {
  # The supremum norm measures deviation on the group-normalized scale, so the
  # sampling weights enter only through the constraints. Unlike the reference
  # implementation the fit is therefore well defined with sampling weights, and it
  # returns balanced, group-normalized weights.
  data <- sim_binary()
  sw <- withr::with_seed(11, stats::runif(nrow(data), 0.3, 3))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(norm = "linf"),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05),
    sampling_weights = sw
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))
  treated <- data$exposure == 1
  expect_equal(sum(w[treated]), sum(sw[treated]), tolerance = 1e-3)
  expect_equal(sum(w[!treated]), sum(sw[!treated]), tolerance = 1e-3)
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("tightening the tolerance cannot lower the l1 dispersion", {
  # A smaller tolerance is a smaller feasible set, so the minimized summed
  # absolute deviation cannot fall, the operational form of the
  # minimum-dispersion characterization for the l1 norm.
  data <- sim_binary()
  fit_of <- function(tolerance) {
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_sbw(norm = "l1"),
      estimand = "ate",
      constraints = balance_terms(tolerance = tolerance)
    )
  }
  tight <- fit_of(0.01)
  loose <- fit_of(0.1)
  expect_balanced(tight, data, tolerance = 0.01)
  expect_balanced(loose, data, tolerance = 0.1)
  disp_tight <- l1_dispersion(as.numeric(stats::weights(tight)), data$exposure)
  disp_loose <- l1_dispersion(as.numeric(stats::weights(loose)), data$exposure)
  expect_gte(disp_tight, disp_loose - 1e-3)
})

# ---- Live consistency against optweight for the deviation norms ------------

test_that("l1 stable balancing meets the objective tolerance against optweight", {
  skip_on_cran()
  skip_if_not_installed("optweight")

  # Both implementations minimize the same summed absolute deviation on the same
  # tolerance box, so ours attains an l1 dispersion at or below optweight's while
  # both satisfy the band. The comparison is at the objective level because the
  # linear program's solution can be non-unique.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(norm = "l1"),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  # The reference solver emits its own convergence chatter on this linear program;
  # its returned weights are still feasible, which is all the comparison needs.
  reference <- suppressWarnings(optweight::optweight(
    exposure ~ x1 + x2,
    data = data,
    tols = 0.05,
    estimand = "ATE",
    norm = "l1"
  ))

  expect_balanced(fit, data, tolerance = 0.05)
  ours <- l1_dispersion(as.numeric(stats::weights(fit)), data$exposure)
  theirs <- l1_dispersion(reference$weights, data$exposure)
  expect_lte(ours, theirs + 0.01 * theirs + 1e-6)
})

test_that("linf stable balancing meets the objective tolerance against optweight", {
  skip_on_cran()
  skip_if_not_installed("optweight")

  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(norm = "linf"),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  # The reference solver emits its own convergence chatter on this linear program;
  # its returned weights are still feasible, which is all the comparison needs.
  reference <- suppressWarnings(optweight::optweight(
    exposure ~ x1 + x2,
    data = data,
    tols = 0.05,
    estimand = "ATE",
    norm = "linf"
  ))

  expect_balanced(fit, data, tolerance = 0.05)
  ours <- linf_dispersion(as.numeric(stats::weights(fit)), data$exposure)
  theirs <- linf_dispersion(reference$weights, data$exposure)
  expect_lte(ours, theirs + 0.01 * theirs + 1e-6)
})

test_that("a per-covariate tolerance binds each covariate to its own band", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = c(x1 = 0.01, x2 = 0.2))
  )
  achieved <- per_column_smd(fit, data)
  expect_lte(achieved[["x1"]], 0.01 + 1e-6)
  expect_lte(achieved[["x2"]], 0.2 + 1e-6)
})

test_that("derived columns inherit the source covariate tolerance", {
  # A second moment adds a power column whose tolerance is inherited from its
  # source covariate, so the whole column set is balanced within the same band
  # and the recipe records the inherited value on the derived term.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(moments = 2L, tolerance = 0.05)
  )
  power_records <- Filter(function(r) identical(r$kind, "power"), fit@recipe)
  expect_gt(length(power_records), 0)
  expect_true(all(vapply(
    power_records,
    function(r) isTRUE(all.equal(r$tolerance, 0.05)),
    logical(1)
  )))
  expect_balanced(fit, data, tolerance = 0.05)
})

# ---- Constant covariates --------------------------------------------------

test_that("the weighted correlations read a constant column as zero", {
  # The refinement loop compares each achieved correlation against its target.
  # A column with no weighted spread has an undefined correlation, and reporting
  # it as zero keeps the comparison that decides which tolerances still bind
  # from resolving to a missing value.
  exposure <- c(-1, 0, 1, 2, 0.5, -0.5)
  z <- cbind(varying = exposure, constant = rep(0, 6))
  achieved <- sbw_weighted_correlations(exposure, z, rep(1, 6))

  expect_equal(achieved, c(1, 0))
})

test_that("a constant covariate leaves a continuous stable-balancing fit intact", {
  data <- sim_continuous(n = 200)
  data$fixed <- 5
  fit <- balance(
    data,
    exposure,
    c(x1, x2, fixed),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )

  expect_false("fixed" %in% fit@balance_table$term)
  expect_true(all(is.finite(as.numeric(stats::weights(fit)))))
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("a single-level factor leaves a continuous stable-balancing fit intact", {
  data <- sim_continuous(n = 200)
  data$f <- factor(rep("a", nrow(data)))
  fit <- balance(
    data,
    exposure,
    c(x1, x2, f),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )

  expect_false("f_a" %in% fit@balance_table$term)
  expect_true(all(is.finite(as.numeric(stats::weights(fit)))))
  expect_balanced(fit, data, tolerance = 0.05)
})

# ---- Infeasible constraint set --------------------------------------------

test_that("an infeasible constraint set raises balancing_infeasible_error", {
  # A covariate that perfectly separates the exposure groups cannot meet a tight
  # balance band: the treated group's weighted mean of the separating column is
  # fixed at one while the control group's is fixed at zero, so no reweighting
  # brings the standardized mean difference inside a small tolerance and the
  # solver reports primal infeasibility, which maps to the infeasible condition
  # rather than the generic convergence warning. The tolerance is positive so the
  # fit passes the required-tolerance gate and reaches the solver.
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
      method = bw_sbw(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0.01)
    )
  )
})

test_that("the clarabel fallback does not mask genuine infeasibility", {
  # A perfectly-separating covariate is infeasible for both backends. Under the
  # default routing the primary solver certifies infeasibility, the interior-point
  # backend is consulted, and it certifies infeasibility too, so the user still
  # receives one clear infeasible error rather than a silent rescue.
  n <- 60L
  data <- data.frame(
    exposure = rep(0:1, each = n / 2),
    x1 = rep(0:1, each = n / 2),
    x2 = rep_len(c(-1, 0, 1), n)
  )
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_sbw(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0.01)
    ),
    class = "balancing_infeasible_error"
  )
})

# ---- Backend routing and fallback -----------------------------------------

# An ill-scaled binary average-treatment-effect instance: four covariate columns
# spanning 1e-3 to 1e3, left unstandardized so the raw constraint matrix is badly
# conditioned. The default solver falsely certifies this feasible instance
# infeasible; the interior-point backend solves it. This is the small form of the
# misfire the n=20000 instance shows, exercised through the internal solver so the
# routing and the boundary fields can be tested without the large, slow shape.
ill_scaled_instance <- function(n = 2000L, seed = 101) {
  withr::local_seed(seed)
  p <- 4L
  covs <- matrix(0, n, p)
  for (j in seq_len(p)) {
    covs[, j] <- stats::rnorm(n) * 10^(2 * (j - 1) - 3)
  }
  linear_predictor <- 0.5 * (covs[, 1] / 1e-3) + 0.4 * (covs[, 2] / 0.1)
  treat <- as.integer(stats::runif(n) < stats::plogis(linear_predictor))
  list(
    covs = covs,
    treat = treat,
    targets = colMeans(covs),
    tols = 0.05 * apply(covs, 2, stats::sd),
    n = n
  )
}

solve_sbw_internal <- function(instance, backend) {
  solve_sbw <- getFromNamespace("solve_sbw", "balancing")
  solve_sbw(
    instance$treat,
    rep(1, instance$n),
    "ate",
    "l2",
    instance$covs,
    instance$targets,
    instance$tols,
    1e-8,
    list(threads = 1L, backend = backend)
  )
}

test_that("the default routing falls back to clarabel on an osqp infeasibility certificate", {
  instance <- ill_scaled_instance()
  auto <- solve_sbw_internal(instance, "auto")
  reference <- solve_sbw_internal(instance, "clarabel")

  expect_true(auto$converged)
  expect_true(auto$fell_back)
  expect_identical(auto$solver_status, "clarabel")
  weights <- as.numeric(auto$weights)
  expect_true(all(is.finite(weights)))
  expect_true(all(weights >= 0))

  # The rescue matches a direct clarabel solve: the same strictly convex program,
  # solved to the same objective, so the fallback adds no accuracy cost.
  expect_equal(auto$objective, reference$objective, tolerance = 1e-6)

  # The returned weights satisfy the group sums and every moment band, so the
  # rescue is a genuine feasible solution rather than a passed status.
  for (level in c(0L, 1L)) {
    idx <- instance$treat == level
    group_mean <- mean(weights[idx])
    expect_equal(group_mean, 1, tolerance = 1e-4)
    for (j in seq_len(ncol(instance$covs))) {
      achieved <- stats::weighted.mean(instance$covs[idx, j], weights[idx])
      expect_lte(
        abs(achieved - instance$targets[j]),
        instance$tols[j] + 1e-8
      )
    }
  }
})

test_that("an explicit osqp backend does not fall back", {
  # Pinning the primary backend honors the user's choice: the false infeasibility
  # certificate is returned as it stands, with no interior-point rescue, so a user
  # who wants the primary solver's answer gets exactly that.
  instance <- ill_scaled_instance()
  osqp <- solve_sbw_internal(instance, "osqp")
  expect_false(osqp$converged)
  expect_false(osqp$fell_back)
  expect_identical(osqp$status, "primal_infeasible")
  expect_identical(osqp$solver_status, "osqp")
})

test_that("an explicit osqp backend surfaces the false infeasibility through balance", {
  # With the primary backend pinned the fit surfaces the certificate as the
  # infeasible error, the pre-fallback behavior, so the option is a genuine escape
  # hatch. The tiny instance below is genuinely infeasible, standing in for the
  # large shape the primary solver misfires on, since both reach the same path.
  withr::local_options(balancing.qp_backend = "osqp")
  n <- 60L
  data <- data.frame(
    exposure = rep(0:1, each = n / 2),
    x1 = rep(0:1, each = n / 2),
    x2 = rep_len(c(-1, 0, 1), n)
  )
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_sbw(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0.01)
    ),
    class = "balancing_infeasible_error"
  )
})

test_that("an explicit clarabel backend solves and records itself", {
  withr::local_options(balancing.qp_backend = "clarabel")
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  expect_identical(fit@solver_status, "clarabel")
  expect_true(fit@converged)
  expect_balanced(fit, data, tolerance = 0.05)
})

test_that("the backend-fallback alert announces the switch", {
  withr::local_options(balancing.quiet = FALSE)
  expect_snapshot(alert_backend_fallback(), cnd_class = TRUE)
})

# ---- Diagnostics ----------------------------------------------------------

test_that("the fit reports dual variables and the quadratic-program backend", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  expect_s3_class(fit@duals, "data.frame")
  expect_gt(nrow(fit@duals), 0)
  # The default backend is osqp; the solver status reflects it.
  expect_identical(fit@solver_status, "osqp")
})

test_that("a stable balancing fit has no estimating equations", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  expect_null(fit@estimating_equations)
  expect_error(
    estimating_equations(fit),
    class = "balancing_ipw_unsupported_error"
  )
})

test_that("ipw() rejects a stable balancing fit", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  outcome <- stats::lm(x1 ~ exposure, data = data)
  expect_error(
    ipw(fit, outcome),
    class = "balancing_ipw_unsupported_error"
  )
})

# ---- Unsupported estimands ------------------------------------------------

# The overlap estimand is legal only for the covariate balancing propensity
# score, so stable balancing weights reject it for every exposure type with the
# shared unsupported-estimand condition, whose message the entropy slice already
# snapshots.

test_that("the ato estimand raises balancing_estimand_error for a binary exposure", {
  data <- sim_binary()
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_sbw(),
      estimand = "ato",
      constraints = balance_terms(tolerance = 0.05)
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
      method = bw_sbw(),
      estimand = "ato",
      constraints = balance_terms(tolerance = 0.05)
    ),
    class = "balancing_estimand_error"
  )
})

# ---- Live consistency against optweight -----------------------------------

test_that("stable balancing weights meet the objective tolerance against optweight for a binary ate", {
  skip_on_cran()
  skip_if_not_installed("optweight")

  # The quadratic-program family compares objectives rather than weights: both
  # implementations solve the same strictly convex minimum-dispersion program on
  # the same tolerance box, so ours must attain weight dispersion at or below
  # optweight's while both satisfy the band. The dispersion is evaluated in R for
  # both weight vectors from the shared oracle.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  reference <- optweight::optweight(
    exposure ~ x1 + x2,
    data = data,
    tols = 0.05,
    estimand = "ATE"
  )

  expect_balanced(fit, data, tolerance = 0.05)
  ours <- weight_dispersion(as.numeric(stats::weights(fit)), data$exposure)
  theirs <- weight_dispersion(reference$weights, data$exposure)
  expect_lte(ours, theirs + 1e-8)
})

test_that("stable balancing weights meet the objective tolerance against optweight for a binary att", {
  skip_on_cran()
  skip_if_not_installed("optweight")

  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "att",
    constraints = balance_terms(tolerance = 0.05)
  )
  reference <- optweight::optweight(
    exposure ~ x1 + x2,
    data = data,
    tols = 0.05,
    estimand = "ATT",
    focal = 1
  )

  expect_balanced(fit, data, tolerance = 0.05)
  ours <- weight_dispersion(as.numeric(stats::weights(fit)), data$exposure)
  theirs <- weight_dispersion(reference$weights, data$exposure)
  expect_lte(ours, theirs + 1e-8)
})

# ---- Print snapshot -------------------------------------------------------

test_that("a stable balancing fit prints its summary block", {
  # Records on the first successful run once the fit path exists.
  data <- sim_binary()
  expect_snapshot({
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_sbw(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0.05)
    )
    fit
  })
})
