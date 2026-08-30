# Specs for the counterfactual mean at each exposure level, which a binary or
# categorical `ipw()` result reports as a row of its own.
#
# The stacked system has always estimated those means: the contrasts it reports
# are transformations of them and each mean is a parameter of the stack, so the
# sandwich already covers them. What was missing was reporting them, which left
# a reader holding the difference between two risks and no way to say what
# either risk was.
#
# The rows lead the table, one per exposure level in the fit's own level order
# with the reference level first, under the effect measure `"mean"` and keyed in
# the `contrast` column by the level they belong to. The contrast rows that
# follow name the pair they compare in that same column, so a binary result is
# read by the rule a categorical one already was rather than leaving its single
# comparison unnamed. A `.by` request repeats the means within each stratum,
# after the whole-sample rows and ahead of the stratum contrasts.
#
# What these pin is that contract: the row set and its order, the labels, the
# arithmetic behind the estimates, and that every surface built from a result
# carries the new rows. The estimates themselves are checked against the stacked
# parameters they are read from and against a g-computation plug-in built from
# the fitted outcome model, rather than against recorded values, so nothing here
# has to be re-recorded when a fixture changes.

# ---- Fixtures --------------------------------------------------------------

# A binary-exposure fixture carrying the exposure twice: as the numeric 0/1
# column the rest of the suite uses, and as a factor whose levels are declared
# in reverse alphabetical order. The reference level is the fit's first, so an
# implementation that sorted the levels itself would name `"no"` the reference
# for the factor arm and every label would come back the other way round.
#
# The modifier enters the outcome model, so a `.by` request has the
# exposure-by-modifier term whose absence it otherwise warns about, and it
# confounds the exposure as well, so the fit has real work to do.
level_means_binary_fixture <- function(n = 600) {
  withr::with_seed(4321, {
    x1 <- stats::rnorm(n)
    x2 <- stats::rbinom(n, 1L, 0.4)
    exposure <- stats::rbinom(n, 1L, stats::plogis(0.4 * x1 - 0.5 * x2))
    y <- stats::rbinom(
      n,
      1L,
      stats::plogis(-0.3 + 0.9 * exposure + 0.6 * x1 - 0.4 * x2)
    )
    y_cont <- 12 +
      0.7 * exposure +
      0.5 * x1 -
      0.3 * x2 +
      stats::rnorm(n)
    data.frame(
      exposure = exposure,
      arm = factor(
        ifelse(exposure == 1L, "yes", "no"),
        levels = c("yes", "no")
      ),
      x1 = x1,
      x2 = x2,
      modifier = factor(ifelse(x2 == 1L, "hi", "lo"), levels = c("lo", "hi")),
      y = y,
      y_cont = y_cont
    )
  })
}

# A three-level categorical fixture whose levels all carry distinct marginal
# means, built on the shared categorical process with an outcome drawn here
# under its own seed.
level_means_categorical_fixture <- function(n = 600) {
  data <- sim_categorical(n)
  withr::with_seed(4322, {
    linear_predictor <- -0.3 +
      0.4 * (data$exposure == "b") +
      0.8 * (data$exposure == "c") +
      0.5 * data$x1
    data$y <- stats::rbinom(n, 1L, stats::plogis(linear_predictor))
  })
  data
}

# A weighted outcome model of the shape `ipw()` expects, with the weights riding
# along as a column so the model frame resolves them.
fit_level_means_outcome <- function(formula, data, wts, family) {
  data[[".wts"]] <- wts
  suppressWarnings(
    stats::glm(formula, data = data, family = family, weights = .wts)
  )
}

# ---- Oracles ---------------------------------------------------------------

# The counterfactual mean at one exposure level: the outcome model predicted
# with the exposure column set to that level for every unit, averaged over the
# target population. `tilt` is that population's weight, one per unit, which is
# uniform for a pooled estimand and the focal group's indicator for a focal one,
# and it takes a stratum's indicator on top for a subgroup's mean. The oracle
# reads the fitted model and the data rather than anything the result carries.
level_mean_plugin <- function(
  outcome_mod,
  data,
  exposure_name,
  value,
  tilt = NULL
) {
  counterfactual <- data
  counterfactual[[exposure_name]] <- value
  stats::weighted.mean(
    stats::predict(outcome_mod, newdata = counterfactual, type = "response"),
    tilt %||% rep(1, nrow(data))
  )
}

# Every level of an exposure, in the order the counterfactual designs fix them,
# which is the factor's own level order and the sorted values of a numeric
# column.
level_mean_values <- function(values) {
  if (is.factor(values)) {
    return(lapply(levels(values), function(level) {
      factor(level, levels = levels(values))
    }))
  }
  as.list(sort(unique(values)))
}

# The rows of a stored estimates table on either side of the boundary this file
# is about.
level_mean_rows <- function(estimates) {
  estimates[estimates$effect == "mean", , drop = FALSE]
}

level_contrast_rows <- function(estimates) {
  estimates[estimates$effect != "mean", , drop = FALSE]
}

# ---- The row set a binary result reports -----------------------------------

test_that("a binary result leads with one mean row per exposure level", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure,
    data,
    w,
    stats::binomial()
  )

  estimates <- ipw(fit, outcome_mod)$estimates

  # The means come first, in level order with the reference level leading, and
  # each is keyed by the level it belongs to rather than by a pair of levels.
  expect_identical(
    estimates$effect,
    c("mean", "mean", "rd", "log(rr)", "log(or)")
  )
  expect_identical(estimates$contrast, c("0", "1", rep("1 vs 0", 3L)))
  expect_identical(names(estimates)[1:2], c("effect", "contrast"))
})

test_that("a binary factor exposure names its means for the fit's own levels", {
  data <- level_means_binary_fixture()
  fit <- balance(data, arm, c(x1, x2), method = bw_entropy(), estimand = "ate")
  expect_identical(fit@exposure_levels, c("yes", "no"))

  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(y ~ arm, data, w, stats::binomial())

  estimates <- ipw(fit, outcome_mod)$estimates

  # The factor declares its levels out of alphabetical order, so the reference
  # level is `"yes"` and every row names the level the fit gave it rather than
  # the one a sort would have.
  expect_identical(estimates$contrast, c("yes", "no", rep("no vs yes", 3L)))
  expect_equal(
    level_mean_rows(estimates)$estimate,
    vapply(
      level_mean_values(data$arm),
      function(value) level_mean_plugin(outcome_mod, data, "arm", value),
      numeric(1)
    ),
    tolerance = 1e-8
  )
})

test_that("a binary gaussian outcome reports its means beside the difference", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y_cont ~ exposure,
    data,
    w,
    stats::gaussian()
  )

  estimates <- ipw(fit, outcome_mod)$estimates

  # A continuous outcome reports one contrast rather than three, and the means
  # are on the outcome's own scale, so nothing confines them to the unit
  # interval the risks of a binomial model sit in.
  expect_identical(estimates$effect, c("mean", "mean", "diff"))
  expect_identical(estimates$contrast, c("0", "1", "1 vs 0"))
  expect_true(all(level_mean_rows(estimates)$estimate > 1))
})

# ---- The arithmetic behind the estimates -----------------------------------

test_that("the binary mean rows are the mean block of the stacked system", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod)

  # The stack carries the marginal means as the parameters `mu0` and `mu1`. The
  # reported rows are those parameters rather than a second computation of the
  # same quantity, and their standard errors are the diagonal of the same
  # sandwich the contrast rows are read off.
  means <- level_mean_rows(result$estimates)
  expect_equal(
    means$estimate,
    unname(result$fit$theta[c("mu0", "mu1")]),
    tolerance = 1e-12
  )
  expect_equal(
    means$std.err,
    unname(sqrt(diag(result$fit$vcov))[c("mu0", "mu1")]),
    tolerance = 1e-12
  )
})

test_that("the binary mean rows match a g-computation plug-in", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))

  # An adjusted model predicts a value per unit, so the standardization is part
  # of the estimand rather than a formality a marginal model would hide.
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure + x1,
    data,
    w,
    stats::binomial()
  )

  estimates <- ipw(fit, outcome_mod)$estimates

  expect_equal(
    level_mean_rows(estimates)$estimate,
    vapply(
      level_mean_values(data$exposure),
      function(value) level_mean_plugin(outcome_mod, data, "exposure", value),
      numeric(1)
    ),
    tolerance = 1e-8
  )
})

test_that("a focal estimand standardizes its mean rows over the focal group", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att",
    .focal_level = "1"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure + x1,
    data,
    w,
    stats::binomial()
  )

  estimates <- ipw(fit, outcome_mod)$estimates
  values <- level_mean_values(data$exposure)
  tilt <- as.numeric(data$exposure == 1L)

  focal <- vapply(
    values,
    function(value) {
      level_mean_plugin(outcome_mod, data, "exposure", value, tilt = tilt)
    },
    numeric(1)
  )
  pooled <- vapply(
    values,
    function(value) level_mean_plugin(outcome_mod, data, "exposure", value),
    numeric(1)
  )

  expect_equal(
    level_mean_rows(estimates)$estimate,
    focal,
    tolerance = 1e-8
  )

  # The two standardizations disagree on this fixture, so the assertion above
  # is a check on the population averaged over rather than on the predictions
  # alone.
  expect_false(isTRUE(all.equal(focal, pooled)))
})

test_that("each binary contrast is the transform of the means above it", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure,
    data,
    w,
    stats::binomial()
  )

  estimates <- ipw(fit, outcome_mod)$estimates
  mu <- stats::setNames(
    level_mean_rows(estimates)$estimate,
    level_mean_rows(estimates)$contrast
  )
  effects <- stats::setNames(
    level_contrast_rows(estimates)$estimate,
    level_contrast_rows(estimates)$effect
  )

  # A counterfactual risk is a probability, so both means lie strictly inside
  # the unit interval and every transform below is defined.
  expect_true(all(mu > 0 & mu < 1))

  expect_equal(effects[["rd"]], mu[["1"]] - mu[["0"]], tolerance = 1e-8)
  expect_equal(
    effects[["log(rr)"]],
    log(mu[["1"]]) - log(mu[["0"]]),
    tolerance = 1e-8
  )
  expect_equal(
    effects[["log(or)"]],
    stats::qlogis(mu[["1"]]) - stats::qlogis(mu[["0"]]),
    tolerance = 1e-8
  )
})

test_that("every binary mean row carries usable inference", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure,
    data,
    w,
    stats::binomial()
  )

  means <- level_mean_rows(ipw(fit, outcome_mod)$estimates)

  expect_finite_column(means, "std.err")
  expect_column_all(means, "std.err", function(x) x > 0)
  expect_column_all(means, "ci.lower", function(x) x < means$estimate)
  expect_column_all(means, "estimate", function(x) x < means$ci.upper)
  expect_equal(means$z, means$estimate / means$std.err, tolerance = 1e-12)

  # The interval is the normal approximation the rest of the table is built on.
  expect_equal(
    (means$ci.upper - means$estimate) / means$std.err,
    rep(stats::qnorm(0.975), 2L),
    tolerance = 1e-8
  )
})

# ---- The covariance the estimates carry ------------------------------------

test_that("the stored covariance covers the binary mean rows", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))

  # The model adjusts so that the two arms share a covariate coefficient. A
  # marginal model at a pooled estimand fits each arm from its own units and
  # nothing else, which leaves the pair of means orthogonal to working
  # precision, so the coupling asserted below would be a property of the
  # adjustment rather than of the block.
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure + x1,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod)
  covariance <- attr(result$estimates, "ipw_vcov", exact = TRUE)
  labels <- c(
    "mean 0",
    "mean 1",
    "rd 1 vs 0",
    "log(rr) 1 vs 0",
    "log(or) 1 vs 0"
  )

  expect_identical(dim(covariance), c(5L, 5L))
  expect_identical(dimnames(covariance), list(labels, labels))
  expect_equal(
    sqrt(diag(covariance)),
    stats::setNames(result$estimates$std.err, labels),
    tolerance = 1e-12
  )

  # Both means read the same adjusted model over the same units, so they covary
  # rather than being independent estimates. A block assembled per row would
  # report an exact zero there.
  expect_true(is.finite(covariance[["mean 0", "mean 1"]]))
  expect_gt(abs(covariance[["mean 0", "mean 1"]]), 1e-8)

  # A risk difference is the exposed arm's mean less the reference arm's, so it
  # moves with the one and against the other.
  expect_lt(covariance[["mean 0", "rd 1 vs 0"]], 0)
  expect_gt(covariance[["mean 1", "rd 1 vs 0"]], 0)
})

# ---- The surfaces built from the result ------------------------------------

test_that("the accessors report the binary mean rows", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod)
  labels <- c(
    "mean 0",
    "mean 1",
    "rd 1 vs 0",
    "log(rr) 1 vs 0",
    "log(or) 1 vs 0"
  )

  expect_identical(names(stats::coef(result)), labels)
  expect_identical(
    unname(stats::coef(result)),
    result$estimates$estimate
  )
  expect_identical(dimnames(stats::vcov(result)), list(labels, labels))

  bounds <- stats::confint(result)
  expect_identical(rownames(bounds), labels)
  expect_identical(unname(bounds[, 1L]), result$estimates$ci.lower)
  expect_identical(unname(bounds[, 2L]), result$estimates$ci.upper)
})

test_that("a coerced binary result reports the mean rows under term", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod)
  coerced <- as.data.frame(result, conf.int = TRUE)

  expect_identical(names(coerced)[1:2], c("term", "contrast"))
  expect_identical(coerced$term, result$estimates$effect)
  expect_identical(coerced$contrast, result$estimates$contrast)
  expect_identical(coerced$estimate, result$estimates$estimate)
  expect_identical(coerced$conf.low, result$estimates$ci.lower)
})

test_that("exponentiating a coerced binary result leaves the mean rows alone", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod)
  plain <- as.data.frame(result, conf.int = TRUE)
  exponentiated <- as.data.frame(
    result,
    conf.int = TRUE,
    exponentiate = TRUE
  )

  # A counterfactual risk is not a ratio, so moving the ratio rows to their
  # natural scale must leave it where it is and must not rename it to something
  # it is not.
  is_mean <- plain$term == "mean"
  expect_identical(exponentiated$term[is_mean], rep("mean", 2L))
  expect_identical(exponentiated$contrast[is_mean], plain$contrast[is_mean])
  expect_identical(exponentiated$estimate[is_mean], plain$estimate[is_mean])
  expect_identical(exponentiated$conf.low[is_mean], plain$conf.low[is_mean])
  expect_identical(exponentiated$conf.high[is_mean], plain$conf.high[is_mean])

  # The ratio rows still move, so the agreement above is with rows that were
  # left alone rather than with a table nothing happened to.
  ratio <- plain$term %in% c("log(rr)", "log(or)")
  expect_identical(exponentiated$term[ratio], c("rr", "or"))
  expect_equal(
    exponentiated$estimate[ratio],
    exp(plain$estimate[ratio]),
    tolerance = 1e-12
  )
})

# ---- A categorical exposure ------------------------------------------------

test_that("a categorical result reports one mean row per level in level order", {
  data <- level_means_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod)
  estimates <- result$estimates

  expect_identical(
    estimates$effect,
    c(rep("mean", 3L), rep(c("rd", "log(rr)", "log(or)"), times = 2L))
  )
  expect_identical(
    estimates$contrast,
    c("a", "b", "c", rep(c("b vs a", "c vs a"), each = 3L))
  )

  # The rows are the categorical mean block, one parameter per level, in the
  # order the fit codes its levels.
  expect_equal(
    level_mean_rows(estimates)$estimate,
    unname(result$fit$theta[c("mu_a", "mu_b", "mu_c")]),
    tolerance = 1e-12
  )
})

test_that("the categorical mean rows match a g-computation plug-in", {
  data <- level_means_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure + x1,
    data,
    w,
    stats::binomial()
  )

  estimates <- ipw(fit, outcome_mod)$estimates

  expect_equal(
    level_mean_rows(estimates)$estimate,
    vapply(
      level_mean_values(data$exposure),
      function(value) level_mean_plugin(outcome_mod, data, "exposure", value),
      numeric(1)
    ),
    tolerance = 1e-8
  )
})

test_that("a categorical focal estimand standardizes its mean rows over the focal group", {
  data <- level_means_categorical_fixture()
  # The focal level is neither the reference level nor the last one, so a mean
  # block standardized over the wrong group could not pass by coincidence of
  # position.
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "att",
    .focal_level = "b"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure + x1,
    data,
    w,
    stats::binomial()
  )

  estimates <- ipw(fit, outcome_mod)$estimates
  values <- level_mean_values(data$exposure)
  tilt <- as.numeric(data$exposure == "b")

  focal <- vapply(
    values,
    function(value) {
      level_mean_plugin(outcome_mod, data, "exposure", value, tilt = tilt)
    },
    numeric(1)
  )
  pooled <- vapply(
    values,
    function(value) level_mean_plugin(outcome_mod, data, "exposure", value),
    numeric(1)
  )

  expect_equal(
    level_mean_rows(estimates)$estimate,
    focal,
    tolerance = 1e-8
  )

  # The two standardizations disagree on this fixture, so the assertion above
  # is a check on the population averaged over rather than on the predictions
  # alone.
  expect_false(isTRUE(all.equal(focal, pooled)))
})

test_that("each categorical contrast is the transform of the two means it names", {
  data <- level_means_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure,
    data,
    w,
    stats::binomial()
  )

  estimates <- ipw(fit, outcome_mod)$estimates
  mu <- stats::setNames(
    level_mean_rows(estimates)$estimate,
    level_mean_rows(estimates)$contrast
  )
  contrasts <- level_contrast_rows(estimates)

  for (level in c("b", "c")) {
    label <- paste(level, "vs a")
    rows <- contrasts[contrasts$contrast == label, , drop = FALSE]
    values <- stats::setNames(rows$estimate, rows$effect)
    expect_equal(
      values[["rd"]],
      mu[[level]] - mu[["a"]],
      tolerance = 1e-8,
      label = label
    )
    expect_equal(
      values[["log(rr)"]],
      log(mu[[level]]) - log(mu[["a"]]),
      tolerance = 1e-8,
      label = label
    )
    expect_equal(
      values[["log(or)"]],
      stats::qlogis(mu[[level]]) - stats::qlogis(mu[["a"]]),
      tolerance = 1e-8,
      label = label
    )
  }
})

test_that("a categorical covariance keys its mean rows by level", {
  data <- level_means_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod)
  covariance <- attr(result$estimates, "ipw_vcov", exact = TRUE)

  expect_identical(dim(covariance), c(9L, 9L))
  expect_identical(
    rownames(covariance)[1:3],
    c("mean a", "mean b", "mean c")
  )
  expect_identical(anyDuplicated(rownames(covariance)), 0L)
  expect_equal(
    unname(sqrt(diag(covariance))),
    result$estimates$std.err,
    tolerance = 1e-12
  )

  # Both contrasts are measured against the reference level, so both covary
  # with the reference level's mean.
  expect_gt(abs(covariance[["mean a", "rd b vs a"]]), 1e-8)
  expect_gt(abs(covariance[["mean a", "rd c vs a"]]), 1e-8)
})

# ---- A `.by` request -------------------------------------------------------

test_that("a .by fit reports the means overall and within every stratum", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  estimates <- ipw(fit, outcome_mod, .by = modifier)$estimates
  means <- level_mean_rows(estimates)

  # Two means over the whole sample and two within each of the two strata. A
  # contrast of strata compares two effects and has no mean of its own.
  expect_identical(nrow(means), 6L)
  expect_identical(
    means$group,
    c(
      rep("overall", 2L),
      rep("modifier = lo", 2L),
      rep("modifier = hi", 2L)
    )
  )
  expect_identical(means$contrast, rep(c("0", "1"), times = 3L))
  expect_false(any(means$group == "modifier = hi vs modifier = lo"))

  # The whole-sample means lead the table and the stratum means follow the
  # whole-sample contrasts, so a block of means is never split by the contrasts
  # built from it.
  expect_identical(
    which(estimates$effect == "mean"),
    c(1L, 2L, 6L, 7L, 8L, 9L)
  )
})

test_that("a .by stratum mean is the g-computation over that stratum", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure * modifier + x1,
    data,
    w,
    stats::binomial()
  )

  estimates <- ipw(fit, outcome_mod, .by = modifier)$estimates
  values <- level_mean_values(data$exposure)

  for (level in c("lo", "hi")) {
    group <- paste0("modifier = ", level)
    reported <- estimates$estimate[
      estimates$effect == "mean" & estimates$group == group
    ]
    stratum <- as.numeric(data$modifier == level)
    expect_equal(
      reported,
      vapply(
        values,
        function(value) {
          level_mean_plugin(
            outcome_mod,
            data,
            "exposure",
            value,
            tilt = stratum
          )
        },
        numeric(1)
      ),
      tolerance = 1e-8,
      label = group
    )
  }
})

test_that("a .by fit gives every mean row a usable standard error", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure * modifier + x1,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod, .by = modifier)
  means <- level_mean_rows(result$estimates)

  expect_finite_column(means, "std.err")
  expect_column_all(means, "std.err", function(x) x > 0)

  # The stratum means are parameters of the same stacked system, so their
  # standard errors are the diagonal of the sandwich the rest of the table is
  # read off.
  expect_equal(
    unname(sqrt(diag(stats::vcov(result)))),
    result$estimates$std.err,
    tolerance = 1e-12
  )

  # A stratum mean is estimated from that stratum's units alone, so it carries
  # more uncertainty than the whole-sample mean of the same arm. The stratum
  # rows run stratum-major over the same pair of levels, so recycling the
  # whole-sample pair over them pairs each mean with its own arm.
  overall <- means$std.err[means$group == "overall"]
  strata <- means$std.err[means$group != "overall"]
  expect_length(overall, 2L)
  expect_length(strata, 4L)
  expect_true(all(strata > rep(overall, times = 2L)))
})

test_that("a .by fit labels every mean row by level and stratum together", {
  data <- level_means_binary_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod, .by = modifier)
  labels <- paste(
    result$estimates$effect,
    result$estimates$contrast,
    result$estimates$group
  )

  expect_identical(anyDuplicated(labels), 0L)
  expect_identical(names(stats::coef(result)), labels)
  expect_true(
    all(c("mean 0 modifier = lo", "mean 1 modifier = hi") %in% labels)
  )
})

# ---- The exposure types the rows do not reach ------------------------------

# A continuous exposure has no levels to fix, so there is no counterfactual mean
# to report and the surface is the marginal structural model's own coefficients.
# The boundary is stated here because everything above is written from the level
# set, and an implementation reading a level set off a continuous column would
# report means for values no unit repeats.

test_that("a continuous-exposure result reports no mean row", {
  data <- sim_continuous(400)
  data$y_cont <- withr::with_seed(4323, {
    1 +
      0.4 * data$exposure +
      0.5 * data$x1 -
      0.3 * data$x2 +
      stats::rnorm(nrow(data))
  })
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_level_means_outcome(
    y_cont ~ exposure,
    data,
    w,
    stats::gaussian()
  )

  estimates <- ipw(fit, outcome_mod)$estimates

  expect_identical(estimates$effect, "slope")
  expect_false(any(estimates$effect == "mean"))
})
