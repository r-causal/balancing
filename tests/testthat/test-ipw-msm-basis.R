# Specs for the dose-response marginal structural models `ipw()` reports for a
# continuous exposure: a weighted outcome model whose exposure enters through a
# transformation or a basis rather than through one bare term.
#
# The boundary is variable membership. A term is an exposure term when every
# variable it reads is the exposure, so `exposure + I(exposure^2)`,
# `exposure + sin(exposure)`, `poly(exposure, 2)`, `splines::ns(exposure, 3)`
# and `splines::bs(exposure, 3)` are all models of the exposure alone, however
# many design columns they expand to. A term reading a covariate as well
# contributes a coefficient that is a change in the effect per unit of that
# covariate, so no row could name the value it is held at, and such a model
# keeps the refusal it has always had.
#
# What an admitted model reports is one row per exposure-reading coefficient,
# named in `contrast` by the coefficient the fit names, and read straight off
# the weighted fit: nothing is standardized here, so each row is exactly a
# coefficient of the outcome model and the sandwich is what the stack adds.
#
# The vocabulary follows from that. An exposure entering through one column is
# the whole of the dose response, so its coefficient is that response's slope
# everywhere and the row keeps the word `slope`, naming nothing further. An
# exposure entering through several columns has no such row, since a curve has a
# different slope at every dose, so the scale word steps back to `coef` at an
# identity link. A logit reports `log(or)` and a log link `log(rr)` either way,
# because a coefficient of those models is a log ratio whatever column it
# multiplies.
#
# None of this needs a frame. The outcome design is read off the fit rather than
# rebuilt from data, and a continuous exposure standardizes nothing, so a basis
# model is reported from the fitted objects alone even though its model frame
# records the basis and carries no exposure column at all. That is pinned here
# rather than left to hold by accident.

# ---- Fixtures --------------------------------------------------------------

# A continuous exposure with a curved response, so a second basis coefficient is
# a real number rather than noise around zero, and a binary outcome carrying the
# same curve for the logit variants.
#
# The exposure and the balanced covariates come from `sim_continuous_indicator()`
# rather than from `sim_continuous()`, which the single-term continuous specs
# use. Its weaker confounding leaves the largest entropy weight a few hundred
# times the smallest rather than tens of millions of times it, and that is what
# keeps a three-column spline in the outcome model a well-conditioned fit whose
# coefficients and standard errors are worth comparing anything against. The
# outcomes are drawn here under their own seed so the fixture is fixed without
# changing the shared helper.
msm_basis_fixture <- function(n = 400) {
  data <- sim_continuous_indicator(n)
  withr::with_seed(6120, {
    data$y_cont <- 1 +
      0.5 * data$exposure +
      0.25 * data$exposure^2 +
      0.4 * data$x1 +
      stats::rnorm(n)
    data$y <- stats::rbinom(
      n,
      1L,
      stats::plogis(
        -0.3 + 0.6 * data$exposure - 0.25 * data$exposure^2 + 0.3 * data$x1
      )
    )
  })
  data
}

# The fit every model in this file is weighted by. Entropy balancing at exact
# balance is the only route a continuous exposure has to `ipw()`, and the
# average treatment effect is the only estimand it targets.
msm_basis_fit <- function(data) {
  balance(
    data,
    exposure,
    c(x1, g),
    method = bw_entropy(),
    estimand = "ate"
  )
}

# A weighted marginal structural model of whatever shape a test wants. A `NULL`
# family fits a linear model, which is the identity-link case; every other
# family goes through `glm()`, wrapped because balancing weights are not counts
# and a binomial fit says so at every call. The weights ride along as a column
# so the model frame resolves them, as they do everywhere else in this suite.
fit_basis_msm <- function(formula, data, wts, family = NULL) {
  data[[".wts"]] <- wts
  if (is.null(family)) {
    return(stats::lm(formula, data = data, weights = .wts))
  }
  suppressWarnings(
    stats::glm(formula, data = data, family = family, weights = .wts)
  )
}

# One row of a coefficient surface, read by the coefficient it is named after.
# The lookup says how many rows it matched before it returns a value: a key
# matching nothing gives an empty vector, and a comparison between two of those
# holds, so a reader that did not report the count could report agreement
# between rows that are not there.
basis_msm_value <- function(estimates, contrast, column = "estimate") {
  value <- estimates[[column]][estimates$contrast == contrast]
  testthat::expect_identical(length(value), 1L, label = contrast)
  value
}

# The column contract of a multi-row coefficient surface. The coefficient each
# row reports is named in `contrast`, immediately after the measure it
# qualifies, which is where every other surface in this package puts that
# column. Nothing about a dose-response curve is evaluated within a subgroup, so
# there is no `group` column.
expect_basis_msm_columns <- function(estimates) {
  testthat::expect_identical(
    names(estimates),
    c(
      "effect",
      "contrast",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
  testthat::expect_null(estimates[["group"]])
  invisible(estimates)
}

# Everything a coefficient surface has to keep in lockstep with its estimates
# table, written once because the claim is the same for every basis.
#
# Such a result records the conditional reading and supports no other, so its
# accessors key by the outcome model's own coefficient names in one order, and
# the covariance they report is the outcome block of the stacked sandwich rather
# than the one the weighted fit computed for itself while treating its weights
# as fixed. The coefficients come out of one fit of one model, so they covary; a
# surface assembled row by row would report zeros off the diagonal.
#
# The stored estimates frame is checked against that same block. It is not a
# reading the result presents, but `pool_ipw()` reads it rather than the
# accessors, so its standard errors reach a pooled result whether or not
# anything else reports them, and nothing else pins them to a covariance.
expect_basis_msm_accessors <- function(result, effect, outcome_mod) {
  estimates <- result$estimates
  labels <- paste(estimates$effect, estimates$contrast)
  coefficients <- names(stats::coef(outcome_mod))

  testthat::expect_identical(estimates$effect, rep(effect, nrow(estimates)))
  testthat::expect_identical(anyDuplicated(labels), 0L)
  expect_finite_column(estimates, "std.err")
  expect_column_all(estimates, "std.err", function(x) x > 0)
  expect_column_all(estimates, "ci.lower", function(x) x < estimates$ci.upper)

  testthat::expect_identical(result$effects, "conditional")
  testthat::expect_identical(result$readings, "conditional")
  testthat::expect_identical(stats::coef(result), stats::coef(outcome_mod))
  testthat::expect_identical(
    dimnames(stats::vcov(result)),
    list(coefficients, coefficients)
  )
  testthat::expect_identical(rownames(stats::confint(result)), coefficients)

  # Each row of the stored frame reports the standard error of the coefficient
  # it names, which is the square root of that coefficient's diagonal entry in
  # the same block, read at the coefficient rather than by position. A frame
  # built from the wrong rows of the sandwich, or from the covariance the
  # weighted fit computed for itself, differs here by more than rounding.
  testthat::expect_equal(
    unname(sqrt(diag(stats::vcov(result)))[estimates$contrast]),
    estimates$std.err,
    tolerance = 1e-12
  )

  covariance <- stats::vcov(result)
  testthat::expect_equal(covariance, t(covariance), tolerance = 1e-12)
  off_diagonal <- covariance[upper.tri(covariance)]
  testthat::expect_true(all(is.finite(off_diagonal)))
  testthat::expect_gt(max(abs(off_diagonal)), 1e-8)

  # The stacked parameter vector names the stored entries by the estimates
  # table's own labels, which is what lets the stored estimates and standard
  # errors be read back out of the variance system the result carries.
  theta <- result$fit$theta
  testthat::expect_true(all(labels %in% names(theta)))
  testthat::expect_equal(
    unname(theta[labels]),
    estimates$estimate,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    unname(sqrt(diag(result$fit$vcov))[labels]),
    estimates$std.err,
    tolerance = 1e-12
  )

  printed <- paste(capture.output(print(result)), collapse = "\n")
  for (coefficient in coefficients) {
    testthat::expect_match(printed, coefficient, fixed = TRUE)
  }

  invisible(result)
}

# The standard errors a weighted marginal structural model carries when the
# weights are treated as fixed: the empirical sandwich of the model's own score
# with no weight-parameter block above it. It is what `naive_msm_se()` in
# test-ipw.R reports for a single exposure column, read at every column of the
# design instead, which is what a basis needs, since a basis has no one column
# to name. The bread is a central finite difference of the score's column sums,
# so the same function serves a linear model and a logit one.
naive_msm_ses <- function(outcome_mod, wts) {
  design <- stats::model.matrix(outcome_mod)
  n <- nrow(design)
  family <- stats::family(outcome_mod)
  response <- as.numeric(
    stats::model.response(stats::model.frame(outcome_mod))
  )
  beta <- stats::coef(outcome_mod)

  score <- function(coefficients) {
    eta <- as.numeric(design %*% coefficients)
    mu <- family$linkinv(eta)
    (wts * (response - mu) * family$mu.eta(eta) / family$variance(mu)) * design
  }

  meat <- crossprod(score(beta)) / n
  bread <- vapply(
    seq_along(beta),
    function(j) {
      step <- 1e-6 * max(1, abs(beta[[j]]))
      up <- beta
      down <- beta
      up[[j]] <- up[[j]] + step
      down[[j]] <- down[[j]] - step
      (colSums(score(up)) - colSums(score(down))) / (2 * step * n)
    },
    numeric(length(beta))
  )

  bread_inv <- solve(bread)
  covariance <- bread_inv %*% meat %*% t(bread_inv) / n
  stats::setNames(sqrt(diag(covariance)), colnames(design))
}

# The five shapes a dose response is written in, which every claim about the
# surface is made across. Each is a model of the exposure alone: two write the
# curve out term by term and three hand it to a basis constructor.
msm_basis_designs <- list(
  quadratic = list(
    formula = y_cont ~ exposure + I(exposure^2),
    contrast = c("exposure", "I(exposure^2)")
  ),
  sine = list(
    formula = y_cont ~ exposure + sin(exposure),
    contrast = c("exposure", "sin(exposure)")
  ),
  polynomial = list(
    formula = y_cont ~ poly(exposure, 2),
    contrast = c("poly(exposure, 2)1", "poly(exposure, 2)2")
  ),
  natural_spline = list(
    formula = y_cont ~ splines::ns(exposure, 3),
    contrast = c(
      "splines::ns(exposure, 3)1",
      "splines::ns(exposure, 3)2",
      "splines::ns(exposure, 3)3"
    )
  ),
  b_spline = list(
    formula = y_cont ~ splines::bs(exposure, 3),
    contrast = c(
      "splines::bs(exposure, 3)1",
      "splines::bs(exposure, 3)2",
      "splines::bs(exposure, 3)3"
    )
  )
)

# ---- The single-coefficient surface stays where it is ----------------------

# The whole of today's continuous surface, pinned by value rather than by shape,
# because relaxing what the outcome model may look like is not permission to
# change what the models it already accepted report. A model whose exposure
# enters through one bare term keeps the eight-column table, keeps the word
# `slope`, names no coefficient, and returns the same numbers it returns today.

test_that("a bare exposure term reports the slope surface unchanged", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ exposure,
    data,
    as.numeric(stats::weights(fit))
  )

  result <- ipw(fit, outcome_mod)
  estimates <- result$estimates

  expect_named(
    estimates,
    c(
      "effect",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
  expect_identical(nrow(estimates), 1L)
  expect_identical(estimates$effect, "slope")
  expect_null(estimates[["contrast"]])
  expect_null(estimates[["group"]])

  expect_equal(estimates$estimate, 0.5624256945180591, tolerance = 1e-8)
  expect_equal(estimates$std.err, 0.06319511226385538, tolerance = 1e-8)
  expect_equal(estimates$z, 8.899829027437935, tolerance = 1e-8)
  expect_equal(estimates$ci.lower, 0.4385655504819371, tolerance = 1e-8)
  expect_equal(estimates$ci.upper, 0.6862858385541811, tolerance = 1e-8)
  expect_identical(estimates$conf.level, 0.95)
  expect_equal(estimates$p.value, 0)

  # The one reported entry is still named for the link in the stacked parameter
  # vector, and the accessors still label the single row that way.
  expect_identical(names(stats::coef(result)), "slope")
  expect_identical(
    unname(result$fit$theta[["slope"]]),
    estimates$estimate
  )
})

# A basis of a covariate is a covariate term, however many columns it expands
# to, so it contributes no row and leaves the exposure's own coefficient the
# whole of the surface. This holds today and has to keep holding: the relaxation
# is about which terms read the exposure, not about how wide the design is.

test_that("a covariate basis leaves the single slope row alone", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ exposure + splines::ns(x1, 3),
    data,
    as.numeric(stats::weights(fit))
  )

  estimates <- ipw(fit, outcome_mod)$estimates

  expect_identical(nrow(estimates), 1L)
  expect_identical(estimates$effect, "slope")
  expect_null(estimates[["contrast"]])
  expect_length(names(estimates), 8L)
  expect_equal(
    estimates$estimate,
    stats::coef(outcome_mod)[["exposure"]],
    tolerance = 1e-10
  )
})

# ---- Terms reading a covariate keep the refusal ----------------------------

# A term reading the exposure and a covariate together contributes a coefficient
# that is a change in the dose response per unit of that covariate. There is no
# one effect for a row to report and no covariate value a label could name, so
# such a model is refused, whether the mixing is written as a crossing, as an
# explicit interaction, as arithmetic inside `I()`, or as a basis expanded
# against a covariate. The refusal is the one the continuous path already makes.

test_that("a term reading the exposure and a covariate is refused", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  w <- as.numeric(stats::weights(fit))

  refused <- list(
    crossing = fit_basis_msm(y_cont ~ exposure * x1, data, w),
    interaction = fit_basis_msm(y_cont ~ exposure + exposure:x1, data, w),
    product = fit_basis_msm(y_cont ~ I(exposure * x1), data, w),
    basis = fit_basis_msm(y_cont ~ poly(exposure, 2):x1, data, w),
    logit = fit_basis_msm(y ~ exposure * x1, data, w, stats::binomial())
  )

  for (outcome_mod in refused) {
    expect_error(
      ipw(fit, outcome_mod),
      class = "balancing_ipw_input_error"
    )
  }
})

# ---- One row per exposure coefficient --------------------------------------

test_that("a quadratic dose response reports one row per exposure coefficient", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ exposure + I(exposure^2),
    data,
    as.numeric(stats::weights(fit))
  )

  result <- ipw(fit, outcome_mod)
  estimates <- result$estimates

  expect_basis_msm_columns(estimates)
  expect_identical(nrow(estimates), 2L)
  expect_identical(estimates$effect, c("coef", "coef"))

  # Named by the coefficient the fit names, which for a curve written out term
  # by term is the term's own label.
  expect_identical(estimates$contrast, c("exposure", "I(exposure^2)"))
  expect_equal(
    basis_msm_value(estimates, "exposure"),
    stats::coef(outcome_mod)[["exposure"]],
    tolerance = 1e-10
  )
  expect_equal(
    basis_msm_value(estimates, "I(exposure^2)"),
    stats::coef(outcome_mod)[["I(exposure^2)"]],
    tolerance = 1e-10
  )
  expect_basis_msm_accessors(result, "coef", outcome_mod)
})

test_that("a transformed dose term reports one row per exposure coefficient", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ exposure + sin(exposure),
    data,
    as.numeric(stats::weights(fit))
  )

  result <- ipw(fit, outcome_mod)
  estimates <- result$estimates

  # The boundary is which variables a term reads, not which functions it is
  # built from, so a transformation with no polynomial or spline structure is
  # admitted on exactly the same terms as one that has it.
  expect_basis_msm_columns(estimates)
  expect_identical(nrow(estimates), 2L)
  expect_identical(estimates$contrast, c("exposure", "sin(exposure)"))
  expect_equal(
    basis_msm_value(estimates, "sin(exposure)"),
    stats::coef(outcome_mod)[["sin(exposure)"]],
    tolerance = 1e-10
  )
  expect_basis_msm_accessors(result, "coef", outcome_mod)
})

test_that("an orthogonal polynomial basis reports one row per basis coefficient", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ poly(exposure, 2),
    data,
    as.numeric(stats::weights(fit))
  )

  result <- ipw(fit, outcome_mod)
  estimates <- result$estimates

  # One term produced both columns, so the term names neither of them on its
  # own and the row is named after the coefficient instead.
  expect_basis_msm_columns(estimates)
  expect_identical(nrow(estimates), 2L)
  expect_identical(
    estimates$contrast,
    c("poly(exposure, 2)1", "poly(exposure, 2)2")
  )
  expect_equal(
    basis_msm_value(estimates, "poly(exposure, 2)1"),
    stats::coef(outcome_mod)[["poly(exposure, 2)1"]],
    tolerance = 1e-10
  )
  expect_equal(
    basis_msm_value(estimates, "poly(exposure, 2)2"),
    stats::coef(outcome_mod)[["poly(exposure, 2)2"]],
    tolerance = 1e-10
  )
  expect_basis_msm_accessors(result, "coef", outcome_mod)
})

test_that("a natural spline basis reports one row per basis coefficient", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)

  # No bare exposure term anywhere in this model. The boundary is variable
  # membership rather than the presence of a linear term, so a model that is
  # nothing but a spline in the exposure is as reportable as one that is not.
  outcome_mod <- fit_basis_msm(
    y_cont ~ splines::ns(exposure, 3),
    data,
    as.numeric(stats::weights(fit))
  )

  result <- ipw(fit, outcome_mod)
  estimates <- result$estimates

  expect_basis_msm_columns(estimates)
  expect_identical(nrow(estimates), 3L)
  expect_identical(
    estimates$contrast,
    c(
      "splines::ns(exposure, 3)1",
      "splines::ns(exposure, 3)2",
      "splines::ns(exposure, 3)3"
    )
  )
  expect_equal(
    estimates$estimate,
    unname(stats::coef(outcome_mod)[estimates$contrast]),
    tolerance = 1e-10
  )
  expect_basis_msm_accessors(result, "coef", outcome_mod)
})

test_that("a B-spline basis reports one row per basis coefficient", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ splines::bs(exposure, 3),
    data,
    as.numeric(stats::weights(fit))
  )

  result <- ipw(fit, outcome_mod)
  estimates <- result$estimates

  expect_basis_msm_columns(estimates)
  expect_identical(nrow(estimates), 3L)
  expect_identical(
    estimates$contrast,
    c(
      "splines::bs(exposure, 3)1",
      "splines::bs(exposure, 3)2",
      "splines::bs(exposure, 3)3"
    )
  )
  expect_equal(
    estimates$estimate,
    unname(stats::coef(outcome_mod)[estimates$contrast]),
    tolerance = 1e-10
  )
  expect_basis_msm_accessors(result, "coef", outcome_mod)
})

# The intercept is not a causal coefficient and a covariate's columns are not
# either, so the count of rows follows the exposure basis rather than the width
# of the design. Two terms expanding to several columns each, only one of them
# reading the exposure, is what separates the two.

test_that("an exposure basis beside a covariate basis reports only its own rows", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ poly(exposure, 2) + splines::ns(x1, 3),
    data,
    as.numeric(stats::weights(fit))
  )

  result <- ipw(fit, outcome_mod)
  estimates <- result$estimates

  expect_basis_msm_columns(estimates)
  expect_identical(nrow(estimates), 2L)
  expect_identical(
    estimates$contrast,
    c("poly(exposure, 2)1", "poly(exposure, 2)2")
  )
  expect_identical(length(stats::coef(outcome_mod)), 6L)
  expect_equal(
    estimates$estimate,
    unname(stats::coef(outcome_mod)[estimates$contrast]),
    tolerance = 1e-10
  )
  expect_basis_msm_accessors(result, "coef", outcome_mod)
})

# ---- The vocabulary --------------------------------------------------------

# A transformation that expands to one column leaves the exposure entering the
# model through one column, so its coefficient is the slope of the dose response
# on that column everywhere and the row keeps the word `slope` and names no
# coefficient. `coef` is for the surfaces where no row is a slope. Written twice
# because the two arrive differently: one is arithmetic on the exposure and the
# other is a function of it.

test_that("a single transformed exposure column keeps the slope row", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  w <- as.numeric(stats::weights(fit))

  models <- list(
    scaled = list(
      model = fit_basis_msm(y_cont ~ I(exposure / 2), data, w),
      coefficient = "I(exposure/2)"
    ),
    transformed = list(
      model = fit_basis_msm(y_cont ~ sin(exposure), data, w),
      coefficient = "sin(exposure)"
    )
  )

  for (spec in models) {
    estimates <- ipw(fit, spec$model)$estimates

    expect_identical(nrow(estimates), 1L)
    expect_identical(estimates$effect, "slope")
    expect_null(estimates[["contrast"]])
    expect_length(names(estimates), 8L)
    expect_equal(
      estimates$estimate,
      unname(stats::coef(spec$model)[[spec$coefficient]]),
      tolerance = 1e-10
    )
  }
})

# The scale word follows the outcome link, which is what the single-term
# continuous surface already does. A coefficient of a logit model is a log odds
# ratio per unit of whatever column it multiplies, so the word is honest of a
# basis coefficient as it stands and does not step back the way `slope` does.

test_that("a logit marginal structural model reports its coefficients as log odds ratios", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  w <- as.numeric(stats::weights(fit))

  models <- list(
    quadratic = list(
      model = fit_basis_msm(
        y ~ exposure + I(exposure^2),
        data,
        w,
        stats::binomial()
      ),
      contrast = c("exposure", "I(exposure^2)")
    ),
    polynomial = list(
      model = fit_basis_msm(y ~ poly(exposure, 2), data, w, stats::binomial()),
      contrast = c("poly(exposure, 2)1", "poly(exposure, 2)2")
    ),
    natural_spline = list(
      model = fit_basis_msm(
        y ~ splines::ns(exposure, 3),
        data,
        w,
        stats::binomial()
      ),
      contrast = c(
        "splines::ns(exposure, 3)1",
        "splines::ns(exposure, 3)2",
        "splines::ns(exposure, 3)3"
      )
    )
  )

  for (spec in models) {
    result <- ipw(fit, spec$model)
    estimates <- result$estimates

    expect_basis_msm_columns(estimates)
    expect_identical(nrow(estimates), length(spec$contrast))
    expect_identical(estimates$contrast, spec$contrast)
    expect_identical(
      estimates$effect,
      rep("log(or)", length(spec$contrast))
    )
    expect_equal(
      estimates$estimate,
      unname(stats::coef(spec$model)[spec$contrast]),
      tolerance = 1e-10
    )
    expect_basis_msm_accessors(result, "log(or)", spec$model)
  }

  # A logit model whose exposure enters through one bare term keeps `log(or)`
  # and keeps the eight-column table, so the scale word is unchanged by the
  # relaxation and only the naming of the rows is new.
  bare <- ipw(fit, fit_basis_msm(y ~ exposure, data, w, stats::binomial()))
  expect_identical(nrow(bare$estimates), 1L)
  expect_identical(bare$estimates$effect, "log(or)")
  expect_null(bare$estimates[["contrast"]])
})

# ---- The point estimates ---------------------------------------------------

# An identity-link marginal structural model is the estimator: nothing is
# standardized, so each row is exactly a coefficient of the weighted fit and the
# comparison is a closed form rather than an agreement to within a
# g-computation. Read row by row, so a reordering of the surface fails here
# rather than passing on a vector that happens to line up.

test_that("the reported coefficients are the weighted fit's own", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  w <- as.numeric(stats::weights(fit))

  for (spec in msm_basis_designs) {
    outcome_mod <- fit_basis_msm(spec$formula, data, w)
    estimates <- ipw(fit, outcome_mod)$estimates
    beta <- stats::coef(outcome_mod)

    expect_identical(nrow(estimates), length(spec$contrast))
    for (contrast in spec$contrast) {
      expect_equal(
        basis_msm_value(estimates, contrast),
        unname(beta[[contrast]]),
        tolerance = 1e-10
      )
    }

    # The intercept is not a causal coefficient and is not reported, so the
    # surface is one row shorter than the coefficient vector.
    expect_identical(nrow(estimates), length(beta) - 1L)
  }
})

# The same curve written two ways is the same fit: `poly(exposure, 2, raw =
# TRUE)` spans the exposure and its square, exactly as writing the two terms out
# does, so the two models have identical coefficients and identical covariances
# and differ only in how their columns are named. What that pins is that the
# surface reads the fit rather than the formula, and it reaches the whole
# reported covariance rather than its diagonal.

test_that("two parameterizations of one dose response report the same numbers", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  w <- as.numeric(stats::weights(fit))

  written <- ipw(fit, fit_basis_msm(y_cont ~ exposure + I(exposure^2), data, w))
  basis <- ipw(
    fit,
    fit_basis_msm(y_cont ~ poly(exposure, 2, raw = TRUE), data, w)
  )

  expect_identical(
    basis$estimates$contrast,
    c("poly(exposure, 2, raw = TRUE)1", "poly(exposure, 2, raw = TRUE)2")
  )
  expect_equal(
    basis$estimates$estimate,
    written$estimates$estimate,
    tolerance = 1e-10
  )
  expect_equal(
    basis$estimates$std.err,
    written$estimates$std.err,
    tolerance = 1e-10
  )
  expect_equal(
    unname(stats::vcov(basis)),
    unname(stats::vcov(written)),
    tolerance = 1e-10
  )
})

# ---- The standard errors ---------------------------------------------------

# Every reported standard error is finite and positive, and none of them is the
# one a weighted regression reports when it treats its weights as a design
# quantity. Accounting for having estimated the weights moves each coefficient
# of each shape, which is what says the stack reaches the whole coefficient
# block rather than the one column the single-term surface reported.

test_that("the basis standard errors account for having estimated the weights", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  w <- as.numeric(stats::weights(fit))

  for (spec in msm_basis_designs) {
    outcome_mod <- fit_basis_msm(spec$formula, data, w)
    naive <- naive_msm_ses(outcome_mod, w)[spec$contrast]
    estimates <- ipw(fit, outcome_mod)$estimates

    reported <- vapply(
      spec$contrast,
      function(contrast) basis_msm_value(estimates, contrast, "std.err"),
      numeric(1)
    )

    expect_true(all(is.finite(reported)))
    expect_true(all(reported > 0))
    expect_true(all(abs(reported / naive - 1) > 1e-6))
  }
})

# Everything above reads the reported standard errors against the surface's own
# arithmetic: they are the diagonal of the covariance the accessors carry, they
# are not the ones a weighted regression reports when it treats its weights as
# fixed, and two parameterizations of one curve report the same pair of them.
# None of that says they are the right size. A nonparametric bootstrap does, and
# it says it for a basis the way the single-term specs say it for a slope: the
# whole procedure is repeated on resampled data, weights and curve together, and
# the spread of the coefficients it produces is what the reported standard errors
# are held against.
#
# The tolerance is set by two measured quantities rather than a guess. The
# resampling standard deviation carries its own noise, of order five percent at
# this replicate count, and the stacked standard error is anticonservative at a
# few hundred observations, which `ipw()` documents as a ratio of 0.836 at 300
# observations for a slope. Over the seeds 2024, 11 and 777 the reported standard
# error came to 0.92, 0.99 and 1.02 of the bootstrap's for the linear
# coefficient and to 0.87, 0.88 and 0.91 of it for the curvature one. Twenty
# percent is what those two together leave.
#
# A tolerance that wide separates a standard error of the wrong magnitude rather
# than one of the wrong calibration, which is the whole of what an external
# anchor is asked for here: whether the spread the surface reports is the spread
# repeating the procedure produces. The weights-fixed sandwich the test above
# rules out is not ruled out by this one, since on this fixture it sits inside
# the same twenty percent, which is why that comparison is made separately and
# exactly. Each coefficient is compared on its own, so neither of them passes on
# the other's agreement.

test_that("the quadratic basis standard errors track a nonparametric bootstrap", {
  skip_on_cran()
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ exposure + I(exposure^2),
    data,
    as.numeric(stats::weights(fit))
  )
  estimates <- ipw(fit, outcome_mod)$estimates

  n <- nrow(data)
  boot <- withr::with_seed(2024, {
    vapply(
      seq_len(200),
      function(b) {
        idx <- sample.int(n, n, replace = TRUE)
        resampled <- data[idx, , drop = FALSE]
        # A resampled data set can legitimately fail to converge; that replicate
        # drops out through the error handler, and its convergence warning is
        # suppressed so it does not leak into the suite output.
        tryCatch(
          suppressWarnings({
            boot_fit <- msm_basis_fit(resampled)
            boot_mod <- fit_basis_msm(
              y_cont ~ exposure + I(exposure^2),
              resampled,
              as.numeric(stats::weights(boot_fit))
            )
            stats::coef(boot_mod)[c("exposure", "I(exposure^2)")]
          }),
          error = function(e) c(NA_real_, NA_real_)
        )
      },
      numeric(2)
    )
  })
  boot_se <- apply(boot, 1L, stats::sd, na.rm = TRUE)

  expect_identical(estimates$contrast, c("exposure", "I(exposure^2)"))
  for (index in seq_len(nrow(estimates))) {
    expect_equal(
      estimates$std.err[[index]],
      unname(boot_se[[index]]),
      tolerance = 0.2,
      label = estimates$contrast[[index]]
    )
  }
})

# ---- The call and the readings ---------------------------------------------

# A model frame records the term rather than the variables inside it, so the
# frame of a `poly(exposure, 2)` fit carries no exposure column and a route that
# read the exposure off it would have to ask the caller for the data. This one
# never reads it: the outcome design comes off the fit and a continuous exposure
# standardizes nothing, so the fitted objects are the whole of what a basis
# model is reported from. Supplying a frame anyway changes nothing, so neither
# call is the privileged one.

test_that("a basis marginal structural model is reported without a frame", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  w <- as.numeric(stats::weights(fit))

  for (spec in msm_basis_designs) {
    outcome_mod <- fit_basis_msm(spec$formula, data, w)

    result <- ipw(fit, outcome_mod)
    expect_s3_class(result, "ipw")
    expect_identical(result$estimates$contrast, spec$contrast)
    expect_identical(
      ipw(fit, outcome_mod, .data = data)$estimates,
      result$estimates
    )
  }
})

# The conditional reading of any continuous fit is the outcome model's own
# coefficient surface, which for a basis is the whole vector, intercept
# included, rather than the rows the stored estimates table names. The
# covariance it reports is the outcome block of the stacked one, so it accounts
# for the weights having been estimated as well.

test_that("the conditional reading of a basis fit is the whole coefficient vector", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ splines::ns(exposure, 3),
    data,
    as.numeric(stats::weights(fit))
  )

  result <- ipw(fit, outcome_mod)

  # The reading such a result records is the one it presents, so asking for it
  # again is the result that went in rather than a rebuild of it.
  expect_identical(causalgenerics::as_conditional(result), result)

  expect_identical(stats::coef(result), stats::coef(outcome_mod))
  covariance <- stats::vcov(result)
  expect_identical(
    dimnames(covariance),
    list(names(stats::coef(outcome_mod)), names(stats::coef(outcome_mod)))
  )
  expect_true(all(is.finite(covariance)))
  expect_true(all(diag(covariance) > 0))

  # The whole block rather than its diagonal: what the reading reports is the
  # outcome block of the stacked covariance, off-diagonal entries included, so
  # a route that rebuilt it from the model alone would agree nowhere. The block
  # is addressed here by the names the stack carries, which are not the model's:
  # the intercept keeps the `beta_` prefix every column that reads no exposure
  # keeps, and each exposure-reading column is renamed to the label its
  # estimates row is read under. The dimnames differ by construction, since the
  # reading presents the model's own names, so the values are what is compared.
  labels <- paste(result$estimates$effect, result$estimates$contrast)
  keys <- c("beta_(Intercept)", labels)
  expect_true(all(keys %in% names(result$fit$theta)))
  # `keys` spans the whole vector only because this outcome model reads no
  # covariate: every column of it is the intercept or an exposure basis term.
  expect_identical(unname(covariance), unname(result$fit$vcov[keys, keys]))

  # The stored table is the exposure-reading coefficients alone, so it is one
  # row shorter than the vector the reading presents: the intercept is a
  # coefficient of the model and no row of the surface the stack estimated.
  expect_identical(
    nrow(result$estimates),
    length(stats::coef(outcome_mod)) - 1L
  )
})

# An exposure entering the outcome model through several columns has no
# coefficient that is a causal effect: a curve has a different slope at every
# dose, so no row of the table answers the question the marginal reading is
# asked. Such a result therefore declares the conditional reading as the only
# one it supports, and every door into the marginal one is shut, at the
# constructor and at each accessor alike.
#
# The declaration is a default rather than something the caller asked for, so it
# is announced once at construction. The announcement says which reading was
# recorded, why there is no other, where the marginalization this package does
# not compute belongs, and how to stop being told.

test_that("a basis fit announces the reading it records", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ poly(exposure, 2),
    data,
    as.numeric(stats::weights(fit))
  )

  withr::local_options(balancing.quiet = FALSE)
  expect_snapshot(invisible(ipw(fit, outcome_mod)))
})

test_that("a basis fit records the conditional reading and supports no other", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  w <- as.numeric(stats::weights(fit))

  for (spec in msm_basis_designs) {
    outcome_mod <- fit_basis_msm(spec$formula, data, w)
    result <- ipw(fit, outcome_mod)

    # The refusals come from causalgenerics, which owns the reading contract, so
    # what this package settles is the pair of fields it declares. The class is
    # the shared one: nothing here is a subclass carrying refusals of its own.
    expect_identical(class(result), "ipw")
    expect_identical(result$effects, "conditional")
    expect_identical(result$readings, "conditional")
  }
})

test_that("naming the conditional reading builds the same result silently", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ splines::bs(exposure, 3),
    data,
    as.numeric(stats::weights(fit))
  )

  # A caller who named the reading has been told, so the announcement is a
  # default being explained rather than a fact being reported, and it stops.
  named <- withr::with_options(
    list(balancing.quiet = FALSE),
    expect_no_message(ipw(fit, outcome_mod, effects = "conditional"))
  )

  # Naming the reading the default records builds the result the default builds,
  # whole: the same stacked system, the same stored table, and the same wrapped
  # component models.
  expect_identical(named, ipw(fit, outcome_mod))
})

test_that("a basis fit refuses the marginal reading at construction", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ poly(exposure, 2),
    data,
    as.numeric(stats::weights(fit))
  )

  # Asking for the marginal reading of a model that has none is a question
  # rather than a preference, so it is answered rather than quietly given the
  # other reading. It is answered the same way whether or not the announcement
  # would have been printed, which is why the messages are left on: the recorded
  # entry is the error alone, and an announcement reaching a caller who is about
  # to be refused would grow a Message block into it.
  withr::local_options(balancing.quiet = FALSE)
  expect_error(
    ipw(fit, outcome_mod, effects = "marginal"),
    class = "balancing_ipw_input_error"
  )
  expect_balancing_error(ipw(fit, outcome_mod, effects = "marginal"))
})

# A caller who wrote no reading into their own call has asked for nothing, and a
# wrapper forwarding its own `effects` default hands this method the whole
# vector rather than one element of it. That vector is the default, so a result
# whose reading has to change treats it as the default being overridden rather
# than as a marginal request being refused: the announcement is the right answer
# to a caller who never named a reading, and a refusal would be a wrapper's
# argument list refusing calls its author never wrote. Only a single reading is
# a request.

test_that("a forwarded default reading counts as no request", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_basis_msm(y_cont ~ poly(exposure, 2), data, w)

  wrapper <- function(..., effects = c("marginal", "conditional")) {
    ipw(..., effects = effects)
  }

  withr::local_options(balancing.quiet = FALSE)

  # The announcement is compared against the one the bare default prints rather
  # than snapshotted again: what this pins is that the two calls are the same
  # call, and the wording is pinned where it is recorded.
  forwarded <- testthat::evaluate_promise(wrapper(fit, outcome_mod))
  direct <- testthat::evaluate_promise(ipw(fit, outcome_mod))

  expect_gt(length(direct$messages), 0L)
  expect_identical(forwarded$messages, direct$messages)
  expect_identical(forwarded$result$effects, "conditional")
  expect_identical(forwarded$result$readings, "conditional")
  expect_identical(forwarded$result$estimates, direct$result$estimates)

  # A wrapper forwards a reading its own caller did name, so asking for the
  # marginal one through a wrapper meets the refusal that asking for it
  # directly meets.
  expect_error(
    wrapper(fit, outcome_mod, effects = "marginal"),
    class = "balancing_ipw_input_error"
  )

  # A single-column dose has the marginal reading, so the same forwarded vector
  # resolves through the match to the default it always resolved to, with
  # nothing announced.
  single <- fit_basis_msm(y_cont ~ exposure, data, w)
  result <- expect_no_message(wrapper(fit, single))
  expect_identical(result$effects, "marginal")
  expect_identical(result$readings, c("marginal", "conditional"))
})

test_that("the marginal reading is refused wherever it is asked for", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ splines::ns(exposure, 3),
    data,
    as.numeric(stats::weights(fit))
  )

  result <- ipw(fit, outcome_mod)

  # The refusal belongs to causalgenerics, which reads the declared readings off
  # the result, so the classes are asserted rather than the wording: the text is
  # that package's to pin. The specific class comes first so a caller can catch
  # the missing reading by name, and the general one is what everything that
  # refuses any reading shares.
  specific <- "causalgenerics_unsupported_reading_marginal"
  general <- "causalgenerics_unsupported_reading"

  for (cls in c(specific, general)) {
    expect_error(causalgenerics::as_marginal(result), class = cls)
    expect_error(stats::coef(result, effects = "marginal"), class = cls)
    expect_error(stats::vcov(result, effects = "marginal"), class = cls)
    expect_error(stats::confint(result, effects = "marginal"), class = cls)
    expect_error(
      as.data.frame(result, effects = "marginal"),
      class = cls
    )
  }

  # Asking for the reading the result already records is answered rather than
  # refused, which is what makes the refusal about the reading rather than about
  # the argument.
  expect_identical(
    stats::coef(result, effects = "conditional"),
    stats::coef(result)
  )
  expect_identical(
    as.data.frame(result, effects = "conditional"),
    as.data.frame(result)
  )
})

# A single-column dose is unchanged by any of this. Its coefficient is the slope
# of the dose response everywhere, so the marginal reading is a reading it has,
# both readings are declared, the default is the marginal one it always
# reported, and nothing is announced. Written across the two shapes that reach
# it, a model of the exposure alone and one carrying a covariate as well, so the
# width of the design is not what the boundary reads.

test_that("a single-column dose fit keeps both readings and its marginal default", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  w <- as.numeric(stats::weights(fit))

  formulas <- list(y_cont ~ exposure, y_cont ~ exposure + x1)

  for (formula in formulas) {
    outcome_mod <- fit_basis_msm(formula, data, w)

    result <- withr::with_options(
      list(balancing.quiet = FALSE),
      expect_no_message(ipw(fit, outcome_mod))
    )

    expect_identical(class(result), "ipw")
    expect_identical(result$effects, "marginal")
    expect_identical(result$readings, c("marginal", "conditional"))
    expect_identical(result$estimates$effect, "slope")
    expect_null(result$estimates[["contrast"]])

    # Both readings answer, and the round trip between them is the result that
    # went in rather than a rebuild of it.
    conditional <- causalgenerics::as_conditional(result)
    expect_identical(conditional$effects, "conditional")
    expect_identical(causalgenerics::as_marginal(conditional), result)
    expect_identical(names(stats::coef(result)), "slope")
    expect_identical(
      names(stats::coef(result, effects = "conditional")),
      names(stats::coef(outcome_mod))
    )
  }
})

# The modifier is refused for a continuous exposure before the outcome model is
# looked at, so a basis fit meets that refusal rather than the reading one and
# is told nothing about readings on the way.

test_that(".by on a basis fit keeps the continuous-exposure refusal", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ poly(exposure, 2),
    data,
    as.numeric(stats::weights(fit))
  )

  withr::local_options(balancing.quiet = FALSE)
  expect_no_message(expect_error(
    ipw(fit, outcome_mod, .by = g),
    class = "balancing_ipw_unsupported_error"
  ))
})

# The variance system a basis fit reports is the same stacked one every
# continuous fit reports, widened by the basis. The weight parameters lead it
# and keep their own names, and the reported rows follow under the labels the
# estimates table carries.

test_that("a basis fit reports the stacked variance system it was read from", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  outcome_mod <- fit_basis_msm(
    y_cont ~ poly(exposure, 2),
    data,
    as.numeric(stats::weights(fit))
  )

  result <- ipw(fit, outcome_mod)
  parameters <- length(estimating_equations(fit)@parameters)

  expect_identical(result$se_method, "mestimation")
  expect_named(result$fit, c("theta", "vcov"))
  expect_identical(
    length(result$fit$theta),
    parameters + length(stats::coef(outcome_mod))
  )
  expect_identical(
    names(result$fit$theta)[seq_len(parameters)],
    paste0("theta_w", seq_len(parameters))
  )
  expect_identical(
    dimnames(result$fit$vcov),
    list(names(result$fit$theta), names(result$fit$theta))
  )
})

# The surface is the one description of what the exposure contributed to the
# outcome design, and the sandwich has to work it out before it can name the
# stacked parameters it returns. Handing it back is what keeps the caller from
# deriving the same description a second time, so the two cannot disagree about
# which columns carry the dose response or what their rows are called. Both
# shapes are pinned, since a lone exposure column and a basis describe the
# surface differently and only one of them is exercised by the naming above.

test_that("the msm sandwich returns the coefficient surface it named from", {
  data <- msm_basis_fixture()
  fit <- msm_basis_fit(data)
  w <- as.numeric(stats::weights(fit))
  container <- estimating_equations(fit)

  formulas <- list(
    bare = y_cont ~ exposure,
    basis = y_cont ~ poly(exposure, 2)
  )

  for (formula in formulas) {
    outcome_mod <- fit_basis_msm(formula, data, w)
    result <- ipw_deli_msm_sandwich(
      container = container,
      outcome_mod = outcome_mod,
      exposure_name = "exposure",
      sampling_weights = fit@sampling_weights
    )

    expect_true("surface" %in% names(result))
    expect_identical(
      result$surface,
      msm_coefficient_identity(outcome_mod, "exposure")
    )
  }
})
