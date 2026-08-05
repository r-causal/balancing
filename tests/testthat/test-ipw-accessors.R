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
# comparison carries a distinct effect.
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

test_that("a binary ate result carries the covariance of its three effects", {
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

  labels <- c("rd", "log(rr)", "log(or)")
  expect_effect_covariance(result, labels = labels, keys = labels)
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

test_that("a binary att result carries the covariance of its three effects", {
  data <- accessor_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att",
    focal_level = "1"
  )
  w <- stats::weights(fit)
  outcome_mod <- fit_accessor_outcome(
    y ~ exposure + x1,
    data,
    w,
    stats::quasibinomial()
  )

  result <- ipw(fit, outcome_mod)

  labels <- c("rd", "log(rr)", "log(or)")
  expect_effect_covariance(result, labels = labels, keys = labels)
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

# A gaussian outcome reports one measure rather than three, so its covariance
# is the one by one block naming it.
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

  expect_effect_covariance(result, labels = "diff", keys = "diff")
  expect_accessor_contract(
    result,
    labels = "diff",
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
# the measure to the comparison. The stacked names suffix the measure with the
# level instead, which is what the attribute relabels.
test_that("a categorical ate result labels its covariance by comparison", {
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
    "rd b vs a",
    "log(rr) b vs a",
    "log(or) b vs a",
    "rd c vs a",
    "log(rr) c vs a",
    "log(or) c vs a"
  )
  keys <- c(
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

# ---- Where the accessors come from -----------------------------------------

# R records a registered S3 method in the `.__S3MethodsTable__.` of the
# environment where its generic is defined, so reading that table names the
# package a method actually comes from. `getS3method()` is not a substitute: it
# returns `NULL` for a generic that is not visible on the search path, which
# under `R CMD check` would make an absence assertion pass without testing
# anything. test-dependencies.R reads the same tables for the printers; the two
# readers are kept separate because a test file's definitions are local to it.
method_source <- function(name, where) {
  table <- get(".__S3MethodsTable__.", envir = where)
  if (!exists(name, envir = table, inherits = FALSE)) {
    return(NA_character_)
  }
  environmentName(environment(get(name, envir = table, inherits = FALSE)))
}

# `UseMethod()` searches the environment its generic was called from as well as
# the registration table, so a method left behind in balancing's namespace
# would shadow an inherited one even with no registration behind it.
defined_in_balancing <- function(name) {
  exists(name, envir = asNamespace("balancing"), inherits = FALSE)
}

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
})

test_that("balancing registers no S3 method on class ipw", {
  specs <- c(
    "coef.ipw",
    "vcov.ipw",
    "confint.ipw",
    "nobs.ipw",
    "df.residual.ipw",
    "weights.ipw",
    "print.ipw",
    "as.data.frame.ipw",
    "vcov.ipw_model"
  )
  defined <- vapply(specs, defined_in_balancing, logical(1))
  expect_identical(
    defined,
    stats::setNames(rep(FALSE, length(specs)), specs)
  )

  # A method under a generic this list does not name would still be recorded in
  # the table belonging to that generic's package, so both tables the ipw
  # methods live in are read whole rather than only at the expected names.
  registered_sources <- function(where) {
    table <- get(".__S3MethodsTable__.", envir = where)
    names <- grep(
      "\\.ipw(_model)?$",
      ls(table, all.names = TRUE),
      value = TRUE
    )
    vapply(
      names,
      function(name) environmentName(environment(get(name, envir = table))),
      character(1)
    )
  }
  sources <- c(
    registered_sources(baseenv()),
    registered_sources(asNamespace("stats"))
  )
  expect_false("balancing" %in% sources)
})
