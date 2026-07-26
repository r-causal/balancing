# bw_cbps() is the method spec; balance(..., method = bw_cbps()) fits it. The
# covariate balancing propensity score fits a propensity model whose parameters
# satisfy covariate balancing moment conditions. In the just-identified form the
# number of moment conditions equals the number of parameters, so the balancing
# conditions hold exactly and the achieved balance matches the requested
# moments. In the over-identified form the model score equations are stacked onto
# the balancing conditions and a GMM criterion is minimized, so balance is
# approximate and the method reports the criterion value rather than estimating
# equations.
#
# These specs cover the constructor, the property validators, the capability
# methods (including the overlap estimand that only the covariate balancing
# propensity score supports), and the statistical promises: achieved balance
# through expect_balanced(), non-negative weights, estimand-correct group sums,
# ESS bounded by n, the overlap-weight form of the ato estimand, the
# estimating-equations container for the just-identified form, and the design's
# tier-one promise that entropy balancing, inverse probability tilting, and
# just-identified covariate balancing all produce the same
# average-treatment-effect-on-the-treated weights for a binary exposure.

# Normalize weights to mean one within each exposure group so comparisons test
# the weighting solution rather than a reporting convention.
normalize_by_group <- function(w, g) {
  ave <- tapply(w, g, mean)
  as.numeric(w / ave[as.character(g)])
}

# Reconstruct the fitted propensity from the stored logit coefficients. The
# recipe rebuilds the exact standardized constraint matrix the core solved
# against, and the intercept leads the coefficient vector, matching the model
# matrix the fit prepends. This performs only generic matrix algebra with the
# stored pieces; it never re-derives the method's moment conditions.
fitted_propensity <- function(fit, .data) {
  z <- rebuild_constraint_matrix(fit@recipe, .data)
  coefs <- fit@coefficients
  eta <- coefs[[1]] + as.numeric(z %*% coefs[-1])
  stats::plogis(eta)
}

# Just-identified covariate balancing for the average treatment effect equates
# the two exposure arms' weighted covariate means to each other rather than
# balancing each arm to the pooled sample mean, so it is checked arm to arm on
# the standardized scale. This is the exact balance the moment conditions deliver
# and the weighting the reference implementation produces; the shared
# expect_balanced() oracle instead measures each arm against the pooled mean,
# which the arm-to-arm moment conditions do not target.
expect_arms_balanced <- function(fit, .data) {
  constraint_matrix <- rebuild_constraint_matrix(fit@recipe, .data)
  centers <- colMeans(constraint_matrix)
  scales <- apply(constraint_matrix, 2, stats::sd)
  scales[scales == 0] <- 1
  z <- sweep(sweep(constraint_matrix, 2, centers, "-"), 2, scales, "/")

  w <- as.numeric(stats::weights(fit))
  exposure <- as.character(.data[[fit@exposure]])
  arms <- split(seq_along(exposure), exposure)
  arm_means <- lapply(arms, function(idx) {
    apply(z[idx, , drop = FALSE], 2, stats::weighted.mean, w = w[idx])
  })
  testthat::expect_lte(max(abs(arm_means[[2]] - arm_means[[1]])), 1e-6)
}

# The over-identified generalized-method-of-moments criterion the core minimizes,
# evaluated at a coefficient vector. The moment stack is the score residual on the
# model covariates followed by the average-treatment-effect balancing factor times
# the balance covariates, weighted by the two-step pseudo-inverse of the moment
# covariance anchored at the maximum-likelihood fit. This mirrors the core so the
# design's objective-level tolerance can be checked live against WeightIt, whose
# fit does not expose its criterion. The design matrix is the intercept-prefixed
# standardized constraint matrix the core solved against.
cbps_gmm_objective_at <- function(design, treat, beta) {
  n <- nrow(design)
  p <- ncol(design)
  m_total <- 2 * p
  clamp <- function(prob) pmin(pmax(prob, 1e-8), 1 - 1e-8)
  moments <- function(coef) {
    prob <- clamp(stats::plogis(as.vector(design %*% coef)))
    factor <- (prob - treat) / (prob * (1 - prob))
    g <- matrix(0, n, m_total)
    for (j in seq_len(p)) {
      g[, j] <- (treat - prob) * design[, j]
      g[, p + j] <- factor * design[, j]
    }
    g
  }
  anchor <- stats::glm.fit(
    design,
    treat,
    family = stats::binomial()
  )$coefficients
  g_anchor <- moments(anchor)
  covariance <- crossprod(g_anchor) / n
  eig <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  floor <- 1e-12 * max(abs(eig$values))
  inv <- ifelse(eig$values > floor, 1 / eig$values, 0)
  weighting <- eig$vectors %*% diag(inv, m_total, m_total) %*% t(eig$vectors)
  m <- colSums(moments(beta)) / n
  as.numeric(t(m) %*% weighting %*% m)
}

# ---- Constructor ----------------------------------------------------------

test_that("bw_cbps() carries its documented defaults", {
  spec <- bw_cbps()
  expect_true(S7::S7_inherits(spec, bw_cbps))
  expect_true(S7::S7_inherits(spec, estimating_equation_method))
  expect_true(S7::S7_inherits(spec, balance_method))
  expect_false(spec@over_identified)
  expect_true(spec@two_step)
  expect_identical(spec@link, "logit")
  expect_identical(spec@convergence_tolerance, 1e-10)
  expect_null(spec@max_iterations)
})

test_that("bw_cbps() stores supplied tuning parameters", {
  spec <- bw_cbps(
    over_identified = TRUE,
    two_step = FALSE,
    link = "probit",
    convergence_tolerance = 1e-8,
    max_iterations = 200L
  )
  expect_true(spec@over_identified)
  expect_false(spec@two_step)
  expect_identical(spec@link, "probit")
  expect_identical(spec@convergence_tolerance, 1e-8)
  expect_identical(spec@max_iterations, 200L)
})

test_that("bw_cbps() matches the link argument", {
  expect_identical(bw_cbps(link = "cloglog")@link, "cloglog")
  expect_error(bw_cbps(link = "identity"))
})

test_that("bw_cbps() rejects unnamed extra arguments", {
  expect_true(S7::S7_inherits(bw_cbps(), balance_method))
  expect_error(bw_cbps(bogus = 1), class = "balancing_method_error")
})

# ---- Validators -----------------------------------------------------------

test_that("bw_cbps() rejects a non-positive convergence tolerance", {
  expect_identical(bw_cbps()@convergence_tolerance, 1e-10)
  expect_error(bw_cbps(convergence_tolerance = -1e-10))
})

test_that("bw_cbps() rejects a negative iteration cap", {
  expect_null(bw_cbps()@max_iterations)
  expect_error(bw_cbps(max_iterations = -5L))
})

test_that("bw_cbps() rejects non-logical flags", {
  expect_false(bw_cbps(over_identified = FALSE)@over_identified)
  expect_error(bw_cbps(over_identified = "yes"))
  expect_error(bw_cbps(two_step = 1))
})

# ---- Capability methods ---------------------------------------------------

test_that("supported_exposure_types() lists every exposure type", {
  expect_setequal(
    supported_exposure_types(bw_cbps()),
    c("binary", "categorical", "continuous")
  )
})

test_that("supported_estimands() depends on the exposure type", {
  binary <- supported_estimands(bw_cbps(), "binary")
  expect_true(all(c("ate", "att", "ato") %in% binary))
  expect_true(any(c("atc", "atu") %in% binary))

  # The overlap estimand is legal only for a binary exposure.
  expect_setequal(
    supported_estimands(bw_cbps(), "categorical"),
    c("ate", "att")
  )
  expect_false("ato" %in% supported_estimands(bw_cbps(), "categorical"))

  expect_setequal(
    supported_estimands(bw_cbps(), "continuous"),
    "ate"
  )
  expect_false("ato" %in% supported_estimands(bw_cbps(), "continuous"))
})

test_that("supports_estimating_equations() follows the design's rules", {
  # The just-identified discrete form supplies estimating equations.
  expect_true(supports_estimating_equations(bw_cbps()))
  expect_true(supports_estimating_equations(
    bw_cbps(),
    exposure_type = "binary"
  ))
  expect_true(supports_estimating_equations(
    bw_cbps(),
    exposure_type = "categorical"
  ))

  # The over-identified form minimizes a GMM criterion and has none. Only a
  # binary exposure fits that criterion, so with no exposure type supplied the
  # answer covers the binary reading.
  expect_false(supports_estimating_equations(bw_cbps(over_identified = TRUE)))
  expect_false(supports_estimating_equations(
    bw_cbps(over_identified = TRUE),
    exposure_type = "binary"
  ))

  # A continuous exposure has none either.
  expect_false(supports_estimating_equations(
    bw_cbps(),
    exposure_type = "continuous"
  ))
})

test_that("method_label() names the method", {
  expect_identical(
    method_label(bw_cbps()),
    "Covariate balancing propensity score"
  )
})

# ---- Statistical promises: binary just-identified -------------------------

test_that("bw_cbps balances a binary ate", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ate"
  )
  expect_arms_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

test_that("a binary ate fit reports arm-to-arm balance without a balance warning", {
  # Just-identified bw_cbps for the average treatment effect equates the two arms'
  # weighted means, which the balance table reports on the arm-to-arm
  # standardized-mean-difference convention. The achieved imbalance is therefore
  # zero and no balance warning fires, unlike an arm-to-pooled report.
  data <- sim_binary()
  fit <- expect_no_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(),
      estimand = "ate"
    ),
    class = "balancing_balance_warning"
  )
  expect_lt(max(fit@balance_table$weighted), 1e-6)
})

test_that("bw_cbps balances a binary att", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "att"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

test_that("bw_cbps balances a binary atc", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "atc"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

# ---- Statistical promises: categorical ------------------------------------

test_that("bw_cbps balances a categorical ate", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ate"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

test_that("bw_cbps balances a categorical att", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "att",
    focal_level = "b"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

# ---- Statistical promises: continuous -------------------------------------

test_that("bw_cbps balances a continuous ate on the correlation scale", {
  # A continuous exposure balances the weighted exposure-covariate covariance;
  # expect_balanced() reads that on the correlation scale.
  data <- sim_continuous()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ate"
  )
  expect_balanced(fit, data)
  expect_true(all(stats::weights(fit) >= 0))
})

# ---- Group sums -----------------------------------------------------------

test_that("a binary ate normalizes each group to its size", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
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
    method = bw_cbps(),
    estimand = "att"
  )
  w <- as.numeric(stats::weights(fit))
  treated <- data$exposure == 1
  n_treated <- sum(treated)
  expect_equal(w[treated], rep(1, n_treated), tolerance = 1e-6)
  expect_equal(sum(w[!treated]), n_treated, tolerance = 1e-4)
})

# ---- ESS ------------------------------------------------------------------

test_that("the effective sample size is bounded by n within each group", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
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

# ---- Overlap (ato) weights ------------------------------------------------

test_that("bw_cbps ato weights take the overlap form", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ato"
  )
  w <- as.numeric(stats::weights(fit))
  expect_true(all(w >= 0))

  # Overlap weighting tilts by the overlap factor h = p(1 - p): a treated unit
  # is weighted by (1 - p) and a control unit by p, evaluated at the fitted
  # propensity. Within each group the weight is therefore proportional to the
  # corresponding factor, which a per-group reporting rescale cannot disturb.
  p <- fitted_propensity(fit, data)
  treated <- data$exposure == 1
  ratio_treated <- w[treated] / (1 - p[treated])
  ratio_control <- w[!treated] / p[!treated]
  expect_lt(stats::sd(ratio_treated) / mean(ratio_treated), 1e-4)
  expect_lt(stats::sd(ratio_control) / mean(ratio_control), 1e-4)
})

test_that("bw_cbps ato weights differ from the ate weights", {
  data <- sim_binary()
  fit_ato <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ato"
  )
  fit_ate <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ate"
  )
  w_ato <- normalize_by_group(
    as.numeric(stats::weights(fit_ato)),
    data$exposure
  )
  w_ate <- normalize_by_group(
    as.numeric(stats::weights(fit_ate)),
    data$exposure
  )
  expect_false(isTRUE(all.equal(w_ato, w_ate)))
})

# ---- Stored link coefficients ---------------------------------------------

test_that("the fit stores the link coefficients", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ate"
  )
  expect_false(is.null(fit@coefficients))
  expect_type(fit@coefficients, "double")
  expect_gt(length(fit@coefficients), 0)
  expect_true(all(is.finite(fit@coefficients)))
})

# ---- Estimating equations from the core -----------------------------------

test_that("a just-identified binary fit populates consistent estimating equations", {
  # The just-identified moment conditions match the parameter count, so the
  # per-unit estimating functions sum to zero column by column at the solution.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
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

test_that("a just-identified overlap fit populates consistent estimating equations", {
  # The overlap estimand's just-identified conditions are smooth, so it supplies
  # estimating equations like the other just-identified discrete estimands; only
  # the over-identified form and a continuous exposure do not.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ato"
  )
  expect_true(supports_estimating_equations(
    bw_cbps(),
    exposure_type = "binary"
  ))
  ee <- estimating_equations(fit)
  n <- nrow(data)
  p <- ncol(ee@psi)
  expect_equal(nrow(ee@psi), n)
  expect_equal(nrow(ee@jacobian), p)
  expect_equal(ncol(ee@jacobian), p)
  expect_equal(dim(ee@weight_jacobian), dim(ee@psi))
  expect_lt(max(abs(colSums(ee@psi))), 1e-6)
})

test_that("a just-identified categorical fit populates consistent estimating equations", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
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
      method = bw_cbps(link = link),
      estimand = "ate"
    )
    expect_arms_balanced(fit, data)
    expect_true(all(stats::weights(fit) >= 0))
  }
})

# ---- Cross-method identity ------------------------------------------------

test_that("just-identified bw_cbps att weights equal entropy and bw_ipt att weights (binary)", {
  # The design's tier-one promise: entropy balancing, inverse probability
  # tilting, and just-identified covariate balancing solve the same
  # treated-target moment conditions with the logit link, so their
  # average-treatment-effect-on-the-treated weights agree.
  data <- sim_binary()
  fit_cbps <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "att"
  )
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
  w_cbps <- normalize_by_group(
    as.numeric(stats::weights(fit_cbps)),
    data$exposure
  )
  w_ipt <- normalize_by_group(
    as.numeric(stats::weights(fit_ipt)),
    data$exposure
  )
  w_ebal <- normalize_by_group(
    as.numeric(stats::weights(fit_ebal)),
    data$exposure
  )
  expect_equal(w_cbps, w_ipt, tolerance = 1e-6)
  expect_equal(w_cbps, w_ebal, tolerance = 1e-6)
})

# ---- two_step interaction warning -----------------------------------------

test_that("two_step is warned and ignored without over_identified", {
  # The two-step weighting matrix belongs to the over-identified GMM criterion.
  # Setting it while the fit is just-identified has no effect, so the fit warns
  # and proceeds.
  data <- sim_binary()
  expect_warning(
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(two_step = FALSE, over_identified = FALSE),
      estimand = "ate"
    ),
    class = "balancing_warning"
  )
  # The setting is ignored: the just-identified solution does not depend on it.
  reference <- suppressWarnings(balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(two_step = TRUE, over_identified = FALSE),
    estimand = "ate"
  ))
  expect_equal(
    as.numeric(stats::weights(fit)),
    as.numeric(stats::weights(reference)),
    tolerance = 1e-6
  )
})

# ---- over_identified outside a binary exposure ----------------------------

test_that("over_identified is warned and ignored for a categorical exposure", {
  # The generalized-method-of-moments criterion stacks the propensity model's
  # score equations onto the balancing conditions, which the core minimizes for a
  # binary exposure alone. A categorical fit cannot honor the request, so it
  # announces the setting as ignored on the convention two_step follows and
  # returns the just-identified solution.
  data <- sim_categorical()
  expect_warning(
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(over_identified = TRUE),
      estimand = "ate"
    ),
    class = "balancing_ignored_argument_warning"
  )
  reference <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ate"
  )
  expect_equal(
    as.numeric(stats::weights(fit)),
    as.numeric(stats::weights(reference)),
    tolerance = 1e-12
  )
})

test_that("over_identified is warned and ignored for a continuous exposure", {
  # The continuous form balances the exposure-covariate covariance through an
  # exponential tilt and has no criterion to over-identify, so the request is
  # ignored the same way.
  data <- sim_continuous()
  expect_warning(
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(over_identified = TRUE),
      estimand = "ate"
    ),
    class = "balancing_ignored_argument_warning"
  )
  reference <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ate"
  )
  expect_equal(
    as.numeric(stats::weights(fit)),
    as.numeric(stats::weights(reference)),
    tolerance = 1e-12
  )
})

test_that("a binary fit honors over_identified without the ignored warning", {
  # The binary path is the one that fits the criterion, so nothing is ignored
  # there and the solution genuinely departs from the just-identified one.
  data <- sim_binary()
  fit <- expect_no_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(over_identified = TRUE),
      estimand = "ate"
    ),
    class = "balancing_ignored_argument_warning"
  )
  reference <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ate"
  )
  expect_false(isTRUE(all.equal(
    as.numeric(stats::weights(fit)),
    as.numeric(stats::weights(reference))
  )))
})

test_that("the categorical over-identified capability matches its container", {
  # A caller reads supports_estimating_equations() to decide whether
  # estimating_equations() will answer, so the two must agree. The categorical
  # path ignores the over-identified request, which leaves a just-identified fit
  # whose container is real.
  data <- sim_categorical()
  fit <- suppressWarnings(balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(over_identified = TRUE),
    estimand = "ate"
  ))
  expect_true(supports_estimating_equations(
    bw_cbps(over_identified = TRUE),
    exposure_type = "categorical"
  ))
  ee <- estimating_equations(fit)
  expect_true(S7::S7_inherits(ee, balancing_estimating_equations))
  expect_equal(nrow(ee@psi), nrow(data))
  expect_lt(max(abs(colSums(ee@psi))), 1e-6)
  expect_false(is.null(ee@psi_fn))
  expect_false(is.null(ee@weights_fn))
})

test_that("the continuous over-identified capability matches its container", {
  data <- sim_continuous()
  fit <- suppressWarnings(balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(over_identified = TRUE),
    estimand = "ate"
  ))
  expect_false(supports_estimating_equations(
    bw_cbps(over_identified = TRUE),
    exposure_type = "continuous"
  ))
  expect_null(fit@estimating_equations)
  expect_error(
    estimating_equations(fit),
    class = "balancing_ipw_unsupported_error"
  )
})

# ---- Over-identified GMM --------------------------------------------------

test_that("an over-identified fit succeeds and records its criterion", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(over_identified = TRUE),
    estimand = "ate"
  )
  expect_true(all(stats::weights(fit) >= 0))
  # The GMM criterion is recorded on the objective slot.
  expect_type(fit@objective, "double")
  expect_true(is.finite(fit@objective))
  expect_gte(fit@objective, 0)
})

test_that("an over-identified fit has no estimating equations", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(over_identified = TRUE),
    estimand = "ate"
  )
  expect_null(fit@estimating_equations)
  expect_error(
    estimating_equations(fit),
    class = "balancing_ipw_unsupported_error"
  )
})

test_that("the two-step weighting changes the over-identified solution", {
  data <- sim_binary()
  fit_twostep <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(over_identified = TRUE, two_step = TRUE),
    estimand = "ate"
  )
  fit_full <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(over_identified = TRUE, two_step = FALSE),
    estimand = "ate"
  )
  w_twostep <- as.numeric(stats::weights(fit_twostep))
  w_full <- as.numeric(stats::weights(fit_full))
  expect_false(isTRUE(all.equal(w_twostep, w_full)))
})

# ---- Continuous exposure: no estimating equations -------------------------

test_that("a continuous fit has no estimating equations", {
  data <- sim_continuous()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ate"
  )
  expect_null(fit@estimating_equations)
  expect_error(
    estimating_equations(fit),
    class = "balancing_ipw_unsupported_error"
  )
})

# ---- Unsupported estimands ------------------------------------------------

# The overlap estimand is legal only for a binary exposure, so it is unsupported
# for the other two exposure types. Both cases raise the shared unsupported
# estimand condition, whose message the entropy slice already snapshots, so the
# classed assertion is the substance here.

test_that("the ato estimand raises balancing_estimand_error for a categorical exposure", {
  data <- sim_categorical()
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(),
      estimand = "ato"
    ),
    class = "balancing_estimand_error"
  )
})

test_that("the ato estimand raises balancing_estimand_error for a continuous exposure", {
  data <- sim_continuous()
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(),
      estimand = "ato"
    ),
    class = "balancing_estimand_error"
  )
})

# ---- Live consistency against WeightIt ------------------------------------

test_that("just-identified bw_cbps weights match WeightIt for a binary ate", {
  skip_on_cran()
  skip_if_not_installed("WeightIt")

  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ate"
  )
  # over = FALSE selects WeightIt's just-identified (exactly balancing) CBPS.
  reference <- WeightIt::weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "cbps",
    estimand = "ATE",
    over = FALSE
  )

  ours <- normalize_by_group(
    as.numeric(stats::weights(fit)),
    data$exposure
  )
  theirs <- normalize_by_group(reference$weights, data$exposure)
  expect_equal(ours, theirs, tolerance = 1e-6)
})

test_that("just-identified bw_cbps weights match WeightIt for a binary att", {
  skip_on_cran()
  skip_if_not_installed("WeightIt")

  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "att"
  )
  reference <- WeightIt::weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "cbps",
    estimand = "ATT",
    focal = "1",
    over = FALSE
  )

  ours <- normalize_by_group(
    as.numeric(stats::weights(fit)),
    data$exposure
  )
  theirs <- normalize_by_group(reference$weights, data$exposure)
  expect_equal(ours, theirs, tolerance = 1e-6)
})

test_that("over-identified bw_cbps meets the design objective tolerance against WeightIt", {
  skip_on_cran()
  skip_if_not_installed("WeightIt")

  # WeightIt does not expose the generalized-method-of-moments criterion value,
  # and its over-identified fit weights the moments by a different matrix than the
  # two-step pseudo-inverse this core uses, so the two implementations minimize
  # different criteria and their fitted propensity scores need not agree. The
  # design's tolerance policy for the over-identified form is therefore
  # objective-level: our criterion must sit at or below the criterion at
  # WeightIt's fit. That criterion is reconstructed here from WeightIt's fitted
  # propensity using the core's moment and weighting definitions, the live analog
  # of the Rust golden check.
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(over_identified = TRUE, two_step = TRUE),
    estimand = "ate"
  )
  # over = TRUE selects WeightIt's over-identified CBPS; twostep matches ours.
  reference <- WeightIt::weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "cbps",
    estimand = "ATE",
    over = TRUE,
    twostep = TRUE
  )

  design <- cbind(1, rebuild_constraint_matrix(fit@recipe, data))
  levels <- sort(unique(as.character(data$exposure)))
  treat <- as.integer(as.character(data$exposure) == levels[[2]])
  # The coefficients that reproduce WeightIt's fitted propensity in the
  # standardized design, so the criterion is evaluated at WeightIt's fit.
  eta_reference <- stats::qlogis(pmin(
    pmax(as.numeric(reference$ps), 1e-8),
    1 - 1e-8
  ))
  beta_reference <- solve(crossprod(design), crossprod(design, eta_reference))
  reference_objective <- cbps_gmm_objective_at(design, treat, beta_reference)

  expect_lte(fit@objective, reference_objective + 1e-8)
})

# ---- Print snapshot -------------------------------------------------------

test_that("a bw_cbps fit prints its summary block", {
  # Records on the first successful run once the fit path exists.
  data <- sim_binary()
  expect_snapshot({
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(),
      estimand = "ate"
    )
    fit
  })
})
