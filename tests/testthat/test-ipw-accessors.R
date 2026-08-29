# Specs for the covariance a balancing `ipw()` result hands to the shared
# accessors. causalgenerics owns every S3 method on class `ipw`, and those
# methods read only the contract fields of the result plus one attribute:
# `ipw_vcov` on the `estimates` data frame, the covariance of the reported
# effects, labeled the way every surface of the result labels its rows. The
# fitting package attaches it. These specs pin what balancing attaches, on real
# `balance()` plus `ipw()` fits, and pin that balancing registers no method on
# the class itself.
#
# The stacked variance system already carries the numbers. What the attribute
# adds is the block of it that describes the reported effects, under the effect
# labels rather than under the stacked theta names, which are seeds of the
# solver rather than names a caller ever sees. The component outcome model is
# carried the same way: the corrected covariance travels with the model, so
# `vcov(result$outcome_mod)` accounts for having estimated the weights while a
# bare refit of the same weighted model does not.
#
# The same division decides the presentation mode specs at the foot of this
# file. What balancing owes the two readings of a result is the `effects`
# argument that records one at construction and the corrected block the
# conditional reading reports. The reading itself, the accessors that report it,
# and the generics that move a result between the two are causalgenerics'.

# ---- Fixtures --------------------------------------------------------------

# A binary-exposure fixture: the shared data-generating process plus a binary
# and a continuous outcome, drawn under their own seed so the arithmetic is
# fixed without changing the shared helper.
accessor_binary_fixture <- function(n = 200) {
  data <- sim_binary(n)
  withr::with_seed(505, {
    data$y <- stats::rbinom(
      n,
      1L,
      stats::plogis(-0.3 + 0.5 * data$exposure + 0.4 * data$x1)
    )
    data$y_cont <- 1 +
      0.6 * data$exposure +
      0.5 * data$x1 -
      0.3 * data$x2 +
      stats::rnorm(n)
  })
  data
}

# A three-level categorical fixture, whose exposure is the shared process's and
# whose binary outcome depends on both non-reference levels, so every
# contrast carries a distinct effect.
accessor_categorical_fixture <- function(n = 200) {
  data <- sim_categorical(n)
  withr::with_seed(606, {
    linear_predictor <- -0.4 +
      0.7 * (data$exposure == "b") +
      1.1 * (data$exposure == "c") +
      0.5 * data$x1
    data$y <- stats::rbinom(n, 1L, stats::plogis(linear_predictor))
  })
  data
}

# A continuous-exposure fixture with a gaussian outcome, which is the marginal
# structural model whose exposure coefficient the result reports as `slope`.
accessor_continuous_fixture <- function(n = 300) {
  data <- sim_continuous(n)
  withr::with_seed(707, {
    data$y_cont <- 1 +
      0.4 * data$exposure +
      0.5 * data$x1 -
      0.3 * data$x2 +
      stats::rnorm(n)
  })
  data
}

# A weighted outcome model of the form `ipw()` expects. The weights ride along
# as a data column so the model frame resolves them, and they are the fit's own
# `bw` vector rather than a bare numeric: `weights()` on the result reads that
# frame back, and the weight class is part of what it returns.
fit_accessor_outcome <- function(formula, data, weights, family) {
  data[[".wts"]] <- weights
  suppressWarnings(
    stats::glm(formula, data = data, family = family, weights = .wts)
  )
}

# ---- Shared expectations ---------------------------------------------------

# The effect covariance contract. `labels` are the display labels the result's
# own rows carry, spelled out by each caller so the label rule is pinned
# against literal strings rather than recomputed here; `keys` name the same
# entries in the stacked parameter vector, which the attribute is a relabeled
# block of.
expect_effect_covariance <- function(result, labels, keys) {
  covariance <- attr(result$estimates, "ipw_vcov", exact = TRUE)

  # Everything below indexes the matrix, so an absent attribute would raise an
  # error in place of each assertion. Reporting it once on its own is what
  # keeps the reason readable.
  if (!is.matrix(covariance) || !is.numeric(covariance)) {
    testthat::fail(paste0(
      "The estimates carry no numeric `ipw_vcov` matrix. They carry an object ",
      "of class ",
      paste(class(covariance), collapse = "/"),
      "."
    ))
    return(invisible(NULL))
  }

  size <- length(labels)
  testthat::expect_identical(dim(covariance), c(size, size))
  testthat::expect_identical(dimnames(covariance), list(labels, labels))
  testthat::expect_equal(covariance, t(covariance))

  # The reported standard errors are the diagonal of this matrix. The
  # off-diagonals are what the estimates table alone cannot say: the effect
  # measures are transformations of one another's marginal means, so they
  # covary, and a caller combining two of them needs that covariance.
  testthat::expect_equal(
    diag(covariance),
    stats::setNames(result$estimates$std.err^2, labels)
  )
  testthat::expect_equal(
    unname(covariance),
    unname(result$fit$vcov[keys, keys, drop = FALSE])
  )

  # The generic reads the attribute rather than recomputing anything.
  testthat::expect_identical(stats::vcov(result), covariance)
  invisible(NULL)
}

# What the accessors return for a balancing result. These read the contract
# fields alone, so they describe the result the method already builds as much
# as the covariance it gains.
expect_accessor_contract <- function(result, labels, n, weights) {
  estimates <- result$estimates

  testthat::expect_identical(
    stats::coef(result),
    stats::setNames(estimates$estimate, labels)
  )

  # Asked for the level the result was fitted at, the accessor returns the
  # stored bounds rather than recomputing them, so the comparison is exact.
  bounds <- stats::confint(result)
  testthat::expect_identical(
    dimnames(bounds),
    list(labels, c("2.5 %", "97.5 %"))
  )
  testthat::expect_identical(unname(bounds[, 1L]), estimates$ci.lower)
  testthat::expect_identical(unname(bounds[, 2L]), estimates$ci.upper)

  testthat::expect_identical(stats::nobs(result), n)
  testthat::expect_identical(stats::weights(result), weights)

  # `fit` here is the stacked parameter vector and its covariance, a bare list,
  # which reaches `df.residual.default` and yields nothing. The contract maps
  # that to a missing value rather than to an error, so a balancing result
  # reports no residual degrees of freedom and says so.
  testthat::expect_identical(stats::df.residual(result), NA_integer_)
  invisible(NULL)
}

# The component outcome model carries the corrected covariance with it, which
# is the block of the stacked sandwich belonging to the outcome-model
# coefficients, under the names that model calls them by. The block is
# addressed by position because the stacked names are not the model's: they
# carry a `beta_` prefix, and a continuous exposure renames the entry its
# reported effect is read from.
expect_outcome_model_wrap <- function(result, outcome_mod, weight_parameters) {
  wrapped <- result$outcome_mod
  testthat::expect_s3_class(wrapped, "ipw_model")

  coefficients <- names(stats::coef(outcome_mod))
  index <- weight_parameters + seq_along(coefficients)
  block <- result$fit$vcov[index, index, drop = FALSE]
  corrected <- stats::vcov(wrapped)

  testthat::expect_identical(
    dimnames(corrected),
    list(coefficients, coefficients)
  )
  testthat::expect_equal(unname(corrected), unname(block))

  # A bare refit of the same weighted model treats the weights as fixed, which
  # is the variance the correction exists to replace.
  testthat::expect_false(isTRUE(all.equal(
    unname(corrected),
    unname(stats::vcov(outcome_mod))
  )))

  # Wrapping prepends a class and attaches an attribute, so everything else
  # about the model is inherited.
  testthat::expect_equal(
    stats::predict(wrapped, type = "response"),
    stats::predict(outcome_mod, type = "response")
  )
  invisible(NULL)
}

# ---- Binary exposure -------------------------------------------------------

test_that("a binary ate result carries the covariance of its means and effects", {
  data <- accessor_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- stats::weights(fit)
  outcome_mod <- fit_accessor_outcome(
    y ~ exposure,
    data,
    w,
    stats::quasibinomial()
  )

  result <- ipw(fit, outcome_mod)

  labels <- c(
    "mean 0",
    "mean 1",
    "rd 1 vs 0",
    "log(rr) 1 vs 0",
    "log(or) 1 vs 0"
  )
  keys <- c("mu0", "mu1", "rd", "log(rr)", "log(or)")
  expect_effect_covariance(result, labels = labels, keys = keys)
  expect_accessor_contract(
    result,
    labels = labels,
    n = nrow(data),
    weights = w
  )
  expect_outcome_model_wrap(
    result,
    outcome_mod,
    weight_parameters = length(estimating_equations(fit)@parameters)
  )
})

test_that("a binary att result carries the covariance of its means and effects", {
  data <- accessor_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att",
    .focal_level = "1"
  )
  w <- stats::weights(fit)
  outcome_mod <- fit_accessor_outcome(
    y ~ exposure + x1,
    data,
    w,
    stats::quasibinomial()
  )

  result <- ipw(fit, outcome_mod)

  labels <- c(
    "mean 0",
    "mean 1",
    "rd 1 vs 0",
    "log(rr) 1 vs 0",
    "log(or) 1 vs 0"
  )
  keys <- c("mu0", "mu1", "rd", "log(rr)", "log(or)")
  expect_effect_covariance(result, labels = labels, keys = keys)
  expect_accessor_contract(
    result,
    labels = labels,
    n = nrow(data),
    weights = w
  )
  expect_outcome_model_wrap(
    result,
    outcome_mod,
    weight_parameters = length(estimating_equations(fit)@parameters)
  )
})

# A gaussian outcome reports one contrast measure rather than three, so its
# covariance is the block naming the two level means and the difference of them.
test_that("a gaussian-outcome result carries the covariance of its difference", {
  data <- accessor_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- stats::weights(fit)
  outcome_mod <- fit_accessor_outcome(
    y_cont ~ exposure,
    data,
    w,
    stats::gaussian()
  )

  result <- ipw(fit, outcome_mod)

  labels <- c("mean 0", "mean 1", "diff 1 vs 0")
  expect_effect_covariance(
    result,
    labels = labels,
    keys = c("mu0", "mu1", "diff")
  )
  expect_accessor_contract(
    result,
    labels = labels,
    n = nrow(data),
    weights = w
  )
  expect_outcome_model_wrap(
    result,
    outcome_mod,
    weight_parameters = length(estimating_equations(fit)@parameters)
  )
})

# ---- Categorical exposure --------------------------------------------------

# A categorical exposure reports one block of measures per non-reference level,
# so the effect measure alone no longer identifies a row and the labels join
# the measure to the contrast. The stacked names suffix the measure with the
# level instead, which is what the attribute relabels.
test_that("a categorical ate result labels its covariance by contrast", {
  data <- accessor_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- stats::weights(fit)
  outcome_mod <- fit_accessor_outcome(
    y ~ exposure,
    data,
    w,
    stats::quasibinomial()
  )

  result <- ipw(fit, outcome_mod)

  labels <- c(
    "mean a",
    "mean b",
    "mean c",
    "rd b vs a",
    "log(rr) b vs a",
    "log(or) b vs a",
    "rd c vs a",
    "log(rr) c vs a",
    "log(or) c vs a"
  )
  keys <- c(
    "mu_a",
    "mu_b",
    "mu_c",
    "rd_b",
    "log(rr)_b",
    "log(or)_b",
    "rd_c",
    "log(rr)_c",
    "log(or)_c"
  )
  expect_effect_covariance(result, labels = labels, keys = keys)
  expect_accessor_contract(
    result,
    labels = labels,
    n = nrow(data),
    weights = w
  )
  expect_outcome_model_wrap(
    result,
    outcome_mod,
    weight_parameters = length(estimating_equations(fit)@parameters)
  )
})

# ---- Continuous exposure ---------------------------------------------------

# A continuous exposure reports one effect, the exposure coefficient of the
# marginal structural model, so its covariance is a one by one matrix named for
# the model's link. The stack for that fit is the weight parameters and the
# model's coefficients alone, so the outcome-model block is the whole of what
# follows the weight parameters.
test_that("a continuous result carries a one by one effect covariance", {
  data <- accessor_continuous_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- stats::weights(fit)
  outcome_mod <- fit_accessor_outcome(
    y_cont ~ exposure,
    data,
    w,
    stats::gaussian()
  )

  result <- ipw(fit, outcome_mod)

  expect_effect_covariance(result, labels = "slope", keys = "slope")
  expect_accessor_contract(
    result,
    labels = "slope",
    n = nrow(data),
    weights = w
  )
  expect_outcome_model_wrap(
    result,
    outcome_mod,
    weight_parameters = length(estimating_equations(fit)@parameters)
  )
})

# ---- The vcov() generic ----------------------------------------------------

# The accessor reads the attribute and nothing else, so a result that carries
# none raises rather than substituting a diagonal built from the standard
# errors. That the generic succeeds here is the whole of what makes the
# attribute reachable.
test_that("vcov() returns the effect covariance and agrees with coef()", {
  data <- accessor_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- stats::weights(fit)
  outcome_mod <- fit_accessor_outcome(
    y ~ exposure,
    data,
    w,
    stats::quasibinomial()
  )

  result <- ipw(fit, outcome_mod)
  covariance <- stats::vcov(result)

  expect_identical(covariance, attr(result$estimates, "ipw_vcov", exact = TRUE))
  expect_identical(names(stats::coef(result)), colnames(covariance))
  expect_identical(rownames(covariance), colnames(covariance))
})

# ---- The nobs() generic ----------------------------------------------------

# The accessor delegates to the stored outcome model, and a weighted glm counts
# only the rows it was fitted on that carry a nonzero weight. A unit given no
# sampling weight is pinned at zero rather than dropped, so the two numbers a
# caller might read as the sample size come apart: the weight vector is still
# the length of the data the fit saw, while the outcome model counted one row
# fewer. The single zero is what separates those two counts. It is not what
# separates this fit's weights from the binary ate case above: sampling weights
# enter the entropy solve, so every unit's weight moves with them.
test_that("nobs() counts the outcome model's nonzero-weight rows", {
  data <- accessor_binary_fixture()
  data$sw <- rep(1, nrow(data))
  data$sw[7L] <- 0

  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate",
    sampling_weights = sw
  )
  w <- stats::weights(fit)
  expect_identical(as.numeric(w)[7L], 0)

  outcome_mod <- fit_accessor_outcome(
    y ~ exposure,
    data,
    w,
    stats::quasibinomial()
  )
  result <- ipw(fit, outcome_mod)

  expect_identical(stats::nobs(result), nrow(data) - 1L)
  expect_length(stats::weights(result), nrow(data))

  # The count is the outcome model's own, read through the delegation rather
  # than recomputed from the weights.
  expect_identical(stats::nobs(result), as.integer(stats::nobs(outcome_mod)))
})

# ---- The fit's own covariance ----------------------------------------------

# The stacked system a result is built from covers every parameter in it, and
# its leading block belongs to the weight parameters alone. That block is the
# fit's covariance, the one a caller reads to report uncertainty in the weights
# rather than in the effects, so the result stores it on the copy of the fit it
# carries and `vcov()` on that copy reads it back. The block keeps the stacked
# parameter names, `theta_w1` through `theta_wp`, because those name the fit's
# own parameters and neither route renames them: a continuous stack relabels one
# entry, and that entry is the marginal structural model's exposure coefficient,
# which follows the block rather than sitting in it.

test_that("a binary result carries the weight covariance on the fit it stores", {
  data <- accessor_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- stats::weights(fit)
  outcome_mod <- fit_accessor_outcome(
    y ~ exposure,
    data,
    w,
    stats::quasibinomial()
  )

  result <- ipw(fit, outcome_mod)

  # The block covers the parameters the fit solved for, so the container is what
  # says how large it is. The literal is pinned beside it because a block
  # compared only against a size read off the same fit would move with any
  # change to the expansion and report nothing.
  p <- length(estimating_equations(fit)@parameters)
  expect_identical(p, 4L)

  covariance <- stats::vcov(result$wt_mod)
  expect_identical(
    covariance,
    result$fit$vcov[seq_len(p), seq_len(p), drop = FALSE]
  )
  expect_identical(
    dimnames(covariance),
    list(paste0("theta_w", seq_len(p)), paste0("theta_w", seq_len(p)))
  )
})

test_that("a continuous result carries the weight covariance on its fit", {
  data <- accessor_continuous_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- stats::weights(fit)
  outcome_mod <- fit_accessor_outcome(
    y_cont ~ exposure,
    data,
    w,
    stats::gaussian()
  )

  result <- ipw(fit, outcome_mod)

  p <- length(estimating_equations(fit)@parameters)
  expect_identical(p, 5L)

  covariance <- stats::vcov(result$wt_mod)
  expect_identical(
    covariance,
    result$fit$vcov[seq_len(p), seq_len(p), drop = FALSE]
  )
  expect_identical(
    dimnames(covariance),
    list(paste0("theta_w", seq_len(p)), paste0("theta_w", seq_len(p)))
  )

  # The entry this route renames is the effect it reports, which is a
  # coefficient of the marginal structural model. It sits past the weight
  # parameters, so the block carries the same names on either route.
  expect_identical(names(result$fit$theta)[[p + 2L]], "slope")
})

# Nothing in balancing computes a covariance for a fit on its own. The block
# arrives only as a by-product of the stacked assembly an `ipw()` result
# performs, so a fit that has not been through one has no block to report and
# refuses, rather than returning an empty matrix.
#
# One class covers both ways of arriving at that refusal, so both are pinned
# here: a fit whose weights solve estimating equations could reach a covariance
# through an `ipw()` result, and a fit from a method that solves none has no
# such route at all, since `ipw()` turns it away too. What separates them is the
# guidance rather than the class, and test-errors.R records the two wordings.
test_that("a fit on its own reports no covariance", {
  data <- accessor_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )

  expect_error(stats::vcov(fit), class = "balancing_vcov_error")

  no_equations <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  expect_null(no_equations@estimating_equations)
  expect_error(stats::vcov(no_equations), class = "balancing_vcov_error")
})

# S7 objects have value semantics, so filling the covariance in on the fit a
# result stores cannot reach the object the caller passed. That is the whole of
# what this pins: the caller's fit still refuses after a result has been built
# from it, and the stored copy differs from it in the covariance and in nothing
# else, which is what tells one property being filled in apart from the fit
# being rebuilt.
test_that("building a result leaves the fit that went into it alone", {
  data <- accessor_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- stats::weights(fit)
  outcome_mod <- fit_accessor_outcome(
    y ~ exposure,
    data,
    w,
    stats::quasibinomial()
  )

  result <- ipw(fit, outcome_mod)

  expect_true(is.matrix(stats::vcov(result$wt_mod)))
  expect_error(stats::vcov(fit), class = "balancing_vcov_error")

  shared <- setdiff(S7::prop_names(fit), "vcov")
  expect_identical(S7::props(result$wt_mod)[shared], S7::props(fit)[shared])
})

# ---- The presentation mode -------------------------------------------------

# A result reports its effects in one of two readings, recorded in the `effects`
# field the result class contracts to hold. The marginal reading is the causal
# contrast estimates every route reported before the field existed; the
# conditional reading presents the outcome model's coefficient surface. Both
# surfaces exist on every result balancing builds, since the stacked system is
# solved whichever reading is asked for, so the field says which one is
# presented rather than which one was computed, and `ipw()` takes an `effects`
# argument saying which one the result it builds records.

# The binary route's fit and its weighted outcome model, which the specs below
# need together rather than one at a time. The fit is the entropy one every
# binary spec in this file uses.
accessor_binary_models <- function() {
  data <- accessor_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  list(
    fit = fit,
    outcome_mod = fit_accessor_outcome(
      y ~ exposure,
      data,
      stats::weights(fit),
      stats::quasibinomial()
    )
  )
}

# The result an `effects` argument builds, against the same result built without
# the argument and moved to that reading afterwards. The field records which
# reading is presented rather than which one was computed, so naming a reading
# at construction settles the field and nothing else: the same estimates, the
# same covariance attached to them, and the same component models, wrapper
# included. Asserting only the field and the surface it selects would leave a
# construction that computed just the named reading, skipping the marginal
# estimates or the covariance because they are not the ones on show, passing
# every test here.
#
# Whole-object identity is available because the variance system a balancing
# result stores is the stacked parameter vector and its covariance, both plain
# numeric, and the fits and models going in are the same objects. Nothing in it
# carries the closures of the call that produced it, which is what would make
# two solves of one system agree in every number and be identical in none of
# them. The narrower assertions come first so that a mismatch names the field or
# the shape rather than the whole result.
expect_ipw_built_as <- function(result, expected) {
  testthat::expect_identical(result$effects, expected$effects)
  testthat::expect_identical(result$readings, expected$readings)
  testthat::expect_identical(names(result), names(expected))
  testthat::expect_identical(result, expected)
  invisible(result)
}

test_that("every route defaults to the marginal reading and round-trips", {
  binary <- accessor_binary_models()

  categorical_data <- accessor_categorical_fixture()
  categorical_fit <- balance(
    categorical_data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  categorical_mod <- fit_accessor_outcome(
    y ~ exposure,
    categorical_data,
    stats::weights(categorical_fit),
    stats::quasibinomial()
  )

  continuous_data <- accessor_continuous_fixture()
  continuous_fit <- balance(
    continuous_data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  continuous_mod <- fit_accessor_outcome(
    y_cont ~ exposure,
    continuous_data,
    stats::weights(continuous_fit),
    stats::gaussian()
  )

  results <- list(
    binary = ipw(binary$fit, binary$outcome_mod),
    categorical = ipw(categorical_fit, categorical_mod),
    continuous = ipw(continuous_fit, continuous_mod)
  )

  # Marginal is what every route reported before the mode was a field, so a call
  # that names no mode still reports it. Each of these results is one the whole
  # surface exists on, which the declared readings say and the round trip below
  # relies on.
  expect_identical(
    vapply(results, function(res) res$effects, character(1)),
    stats::setNames(rep("marginal", length(results)), names(results))
  )
  for (res in results) {
    expect_identical(res$readings, c("marginal", "conditional"))
  }

  # Both readings exist on every result, so moving to the other one records the
  # move and reads nothing else, and moving back is the result that went in
  # rather than a rebuild of it.
  flipped <- lapply(results, causalgenerics::as_conditional)
  expect_identical(
    vapply(flipped, function(res) res$effects, character(1)),
    stats::setNames(rep("conditional", length(results)), names(results))
  )
  expect_identical(lapply(flipped, causalgenerics::as_marginal), results)

  # Asking for the reading a result already records says what asking once said.
  expect_identical(lapply(flipped, causalgenerics::as_conditional), flipped)
  expect_identical(lapply(results, causalgenerics::as_marginal), results)

  # A continuous fit whose outcome model reads the exposure through several
  # columns is the one result the round trip does not hold for, because the
  # marginal reading is not a reading it has. It declares the conditional one
  # alone and refuses the flip rather than answering it, which is what says the
  # loop above ran on the results that support both rather than on every result
  # the package builds.
  basis_mod <- fit_accessor_outcome(
    y_cont ~ poly(exposure, 2),
    continuous_data,
    stats::weights(continuous_fit),
    stats::gaussian()
  )
  basis <- ipw(continuous_fit, basis_mod)

  expect_identical(basis$effects, "conditional")
  expect_identical(basis$readings, "conditional")
  expect_error(
    causalgenerics::as_marginal(basis),
    class = "causalgenerics_unsupported_reading_marginal"
  )
})

test_that("a binary result records the reading it was built in", {
  models <- accessor_binary_models()

  base <- ipw(models$fit, models$outcome_mod)
  result <- ipw(models$fit, models$outcome_mod, effects = "conditional")

  # Building in the conditional reading is building the result and recording the
  # reading, so it is the default build moved to that reading and nothing else.
  # The stacked system is solved either way, so the corrected block still reaches
  # the outcome model and the marginal estimates are still there to present.
  expect_ipw_built_as(result, causalgenerics::as_conditional(base))
  expect_s3_class(result$outcome_mod, "ipw_model")
  expect_identical(stats::vcov(result), stats::vcov(result$outcome_mod))

  # The reading is a field of the result rather than something the accessors are
  # told at each call, so a result built in the conditional one reports it with
  # nothing named at the call site.
  expect_identical(stats::coef(result), stats::coef(models$outcome_mod))

  # Naming the default is the other half of the argument, and it has to build the
  # result that leaving the argument out builds.
  named <- ipw(models$fit, models$outcome_mod, effects = "marginal")
  expect_ipw_built_as(named, base)
  expect_identical(
    names(stats::coef(named)),
    c("mean 0", "mean 1", "rd 1 vs 0", "log(rr) 1 vs 0", "log(or) 1 vs 0")
  )
})

test_that("a categorical result records the reading it was built in", {
  data <- accessor_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  outcome_mod <- fit_accessor_outcome(
    y ~ exposure,
    data,
    stats::weights(fit),
    stats::quasibinomial()
  )

  base <- ipw(fit, outcome_mod)
  result <- ipw(fit, outcome_mod, effects = "conditional")

  # The conditional reading is of the outcome model, which has one coefficient
  # per non-reference exposure level rather than the six rows the marginal
  # reading of this fit reports.
  expect_ipw_built_as(result, causalgenerics::as_conditional(base))
  expect_identical(stats::coef(result), stats::coef(outcome_mod))
  expect_identical(
    names(stats::coef(result)),
    c("(Intercept)", "exposureb", "exposurec")
  )
})

test_that("a continuous result records the reading it was built in", {
  data <- accessor_continuous_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  outcome_mod <- fit_accessor_outcome(
    y_cont ~ exposure,
    data,
    stats::weights(fit),
    stats::gaussian()
  )

  base <- ipw(fit, outcome_mod)
  result <- ipw(fit, outcome_mod, effects = "conditional")

  # A continuous exposure takes the other branch of the method, whose stack is
  # the weight parameters and the marginal structural model alone, so the
  # argument has to reach the construction from there as well.
  expect_ipw_built_as(result, causalgenerics::as_conditional(base))
  expect_s3_class(result$outcome_mod, "ipw_model")
  expect_identical(stats::coef(result), stats::coef(outcome_mod))
  expect_identical(names(stats::coef(result)), c("(Intercept)", "exposure"))
})

test_that("an invalid effects value errors", {
  models <- accessor_binary_models()

  err <- expect_error(
    ipw(models$fit, models$outcome_mod, effects = "banana"),
    class = "rlang_error"
  )

  # The message names both readings alongside the value, which is what tells a
  # rejection of the value apart from a rejection of the argument itself. This
  # method's dots absorb a name it does not recognize, so a value that reached
  # them would be reported, if at all, as a dots problem naming neither reading.
  expect_match(conditionMessage(err), "marginal", fixed = TRUE)
  expect_match(conditionMessage(err), "conditional", fixed = TRUE)
  expect_match(conditionMessage(err), "banana", fixed = TRUE)

  # The value is settled before the outcome model is looked at, so a call that
  # is wrong in both places reports the value rather than the model. Every
  # refusal this method makes carries the same parent class, so the reading
  # names are what separate them.
  early <- expect_error(
    ipw(models$fit, list(), effects = "banana"),
    class = "rlang_error"
  )
  expect_match(conditionMessage(early), "conditional", fixed = TRUE)
})

test_that("the effects argument reports the other reading for one call", {
  models <- accessor_binary_models()
  result <- ipw(models$fit, models$outcome_mod)

  expect_identical(result$effects, "marginal")
  expect_identical(
    stats::coef(result, effects = "conditional"),
    stats::coef(models$outcome_mod)
  )
  expect_identical(
    stats::vcov(result, effects = "conditional"),
    stats::vcov(result$outcome_mod)
  )

  # Naming a reading at the call site answers in it and leaves the result where
  # it was, so the next call with nothing named answers in the stored one.
  expect_identical(result$effects, "marginal")
  expect_identical(
    stats::coef(result),
    stats::setNames(
      result$estimates$estimate,
      c("mean 0", "mean 1", "rd 1 vs 0", "log(rr) 1 vs 0", "log(or) 1 vs 0")
    )
  )

  # A value that names neither reading is refused at the call site too.
  expect_error(
    stats::coef(result, effects = "banana"),
    class = "causalgenerics_invalid_argument_effects"
  )
})

test_that("a conditional result reports the corrected outcome block", {
  models <- accessor_binary_models()
  result <- ipw(models$fit, models$outcome_mod)
  conditional <- causalgenerics::as_conditional(result)

  expect_identical(stats::coef(conditional), stats::coef(models$outcome_mod))

  # The covariance the conditional reading reports is the outcome block of the
  # stacked sandwich, which balancing attaches to the model itself, rather than
  # anything derived from it here.
  covariance <- stats::vcov(conditional)
  expect_identical(covariance, stats::vcov(result$outcome_mod))
  expect_identical(
    dimnames(covariance),
    list(
      names(stats::coef(models$outcome_mod)),
      names(stats::coef(models$outcome_mod))
    )
  )

  # It is not the covariance the outcome model computed for itself, which treats
  # the estimated weights as fixed and reports an uncertainty the coefficients do
  # not have.
  expect_false(isTRUE(all.equal(
    covariance,
    stats::vcov(models$outcome_mod),
    check.attributes = FALSE
  )))

  # The limits are built from that block at every level, since the ones the
  # result stores belong to the effects the other reading reports.
  interval <- stats::confint(conditional)
  expect_identical(
    rownames(interval),
    names(stats::coef(models$outcome_mod))
  )
  expect_identical(colnames(interval), c("2.5 %", "97.5 %"))
})

test_that("a conditional result prints the outcome model's coefficients", {
  # testthat 3e pins the output width but not the number of significant digits,
  # and `printCoefmat()` wraps its table past 80 columns under a larger `digits`,
  # which splits the rows this test reads by position.
  withr::local_options(digits = 7)

  models <- accessor_binary_models()
  result <- ipw(models$fit, models$outcome_mod)
  out <- capture.output(print(causalgenerics::as_conditional(result)))

  # The reading is named twice, because the two readings are different tables of
  # different numbers: once beside the estimand and once over the table.
  expect_true(any(grepl(
    "Effects: conditional (outcome model)",
    out,
    fixed = TRUE
  )))
  expect_false(any(grepl("Marginal estimates:", out, fixed = TRUE)))

  header <- which(startsWith(out, "Conditional estimates (outcome model):"))
  expect_length(header, 1L)

  # The corrected block is there, so the coefficients are reported beside the
  # standard errors it implies rather than on their own.
  expect_true(grepl("Std. Error", out[[header + 1L]], fixed = TRUE))
  expect_false(any(grepl(
    "Standard errors are not reported",
    out,
    fixed = TRUE
  )))

  # The rows are the outcome model's coefficients rather than the effect
  # measures the marginal reading of the same result tabulates.
  labels <- names(stats::coef(models$outcome_mod))
  rows <- out[seq(header + 2L, header + 1L + length(labels))]
  expect_identical(sub(" .*$", "", rows), labels)
})

# ---- Where the accessors come from -----------------------------------------

# The readers these specs use to name the package a method comes from live in
# helper-s3-registration.R.

test_that("the ipw accessors are the ones causalgenerics registers", {
  specs <- c(
    "coef.ipw",
    "vcov.ipw",
    "confint.ipw",
    "nobs.ipw",
    "df.residual.ipw",
    "weights.ipw",
    "vcov.ipw_model"
  )
  sources <- vapply(
    specs,
    method_source,
    character(1),
    where = asNamespace("stats")
  )
  expect_identical(
    sources,
    stats::setNames(rep("causalgenerics", length(specs)), specs)
  )

  # The generics that move a result between its readings are causalgenerics'
  # own rather than stats', so their methods are recorded in that package's
  # table instead.
  modes <- c("as_marginal.ipw", "as_conditional.ipw")
  mode_sources <- vapply(
    modes,
    method_source,
    character(1),
    where = asNamespace("causalgenerics")
  )
  expect_identical(
    mode_sources,
    stats::setNames(rep("causalgenerics", length(modes)), modes)
  )
})

test_that("balancing registers no S3 method on the ipw result classes", {
  # `pool_ipw()` returns an `ipw_pooled` object, so the pooled result carries the
  # same borrowed accessor surface as a single fit and is swept alongside it.
  specs <- c(
    "coef.ipw",
    "vcov.ipw",
    "confint.ipw",
    "nobs.ipw",
    "df.residual.ipw",
    "weights.ipw",
    "print.ipw",
    "as.data.frame.ipw",
    "vcov.ipw_model",
    "as_marginal.ipw",
    "as_conditional.ipw",
    "estimand.ipw",
    "coef.ipw_pooled",
    "vcov.ipw_pooled",
    "confint.ipw_pooled",
    "nobs.ipw_pooled",
    "df.residual.ipw_pooled",
    "weights.ipw_pooled",
    "print.ipw_pooled",
    "as.data.frame.ipw_pooled",
    "as_marginal.ipw_pooled",
    "as_conditional.ipw_pooled",
    "estimand.ipw_pooled"
  )
  defined <- vapply(specs, defined_in_balancing, logical(1))
  expect_identical(
    defined,
    stats::setNames(rep(FALSE, length(specs)), specs)
  )

  # A method under a generic this list does not name would still be recorded in
  # the table belonging to that generic's package, so every table the ipw
  # methods live in is read whole rather than only at the expected names.
  pattern <- "\\.ipw(_model|_pooled)?$"
  sources <- c(
    registered_sources(baseenv(), pattern),
    registered_sources(asNamespace("stats"), pattern),
    registered_sources(asNamespace("causalgenerics"), pattern)
  )
  expect_false("balancing" %in% sources)
})

# The methods balancing does register belong to the fit class rather than to the
# result, and their generics are `stats`' own, so S7 records them at load time in
# that package's table under the class's fully qualified name. Reading the table
# is what pins the registration rather than the dispatch alone: an S7 method for
# an external generic that never reaches the table is never dispatched to, and
# the table is where that shows.
test_that("the fit's own accessors are balancing's methods", {
  specs <- c(
    "getCall.balancing::balancing",
    "weights.balancing::balancing",
    "vcov.balancing::balancing"
  )
  sources <- vapply(
    specs,
    method_source,
    character(1),
    where = asNamespace("stats")
  )
  expect_identical(
    sources,
    stats::setNames(rep("balancing", length(specs)), specs)
  )
})
