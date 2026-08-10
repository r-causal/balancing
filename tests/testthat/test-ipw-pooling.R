# Specs for pooling balancing `ipw()` results across multiply imputed datasets.
# causalgenerics owns the pooling itself, including Rubin's rules and the
# Barnard-Rubin degrees of freedom; what balancing owes is the entrypoint under
# an unqualified call and results the pooling can actually read. These specs
# pin the workflow end to end rather than the arithmetic, which is that
# package's to prove, except where the arithmetic depends on something
# balancing's results supply.
#
# The workflow is a within-imputation one: each completed dataset is balanced,
# weighted, and estimated on its own, and the per-imputation results are pooled
# afterward. `balance()` cannot run inside `with(imp, ...)`, since it takes its
# data first and the completed frame has no name inside that call, so the
# imputations are walked with `lapply()` over `mice::complete(imp, "all")`.

# A binary-exposure fixture carrying missing values in one covariate. The
# missingness depends on `x2`, which is fully observed, so the values are
# missing at random given the data rather than completely at random, which is
# the condition multiple imputation addresses. The size is small on purpose:
# three imputations of 250 rows is enough to exercise every field the pooled
# object reports without making the suite wait.
ipw_pooling_fixture <- function(n = 250) {
  withr::with_seed(707, {
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    exposure <- stats::rbinom(n, 1L, stats::plogis(0.6 * x1 - 0.4 * x2))
    y <- stats::rbinom(
      n,
      1L,
      stats::plogis(-0.3 + 0.5 * exposure + 0.4 * x1)
    )
    data <- data.frame(exposure = exposure, x1 = x1, x2 = x2, y = y)
    missing <- stats::rbinom(n, 1L, stats::plogis(-1.2 + 0.5 * x2)) == 1L
    data$x1[missing] <- NA_real_
    data
  })
}

# The per-imputation analysis, written the way the documentation shows it: one
# completed dataset in, one `ipw` result out. `quasibinomial()` solves the same
# estimating equation as `binomial()` and does not warn that the weights are
# not counts.
fit_pooling_ipw <- function(data) {
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  data$.wts <- stats::weights(fit)
  ipw(
    fit,
    stats::glm(
      y ~ exposure,
      data = data,
      family = stats::quasibinomial(),
      weights = .wts
    )
  )
}

test_that("pool_ipw() pools balancing results across imputations", {
  skip_if_not_installed("mice")
  data <- ipw_pooling_fixture()
  imp <- withr::with_seed(4321, mice::mice(data, m = 3, print = FALSE))
  fits <- lapply(mice::complete(imp, "all"), fit_pooling_ipw)

  # Unqualified, since the entrypoint being reachable after
  # `library(balancing)` alone is half of what this pins.
  pooled <- pool_ipw(fits)

  expect_s3_class(pooled, "ipw_pooled")
  expect_identical(pooled$m, 3L)
  expect_identical(nrow(pooled$estimates), 3L)
  expect_identical(pooled$estimates$effect, c("rd", "log(rr)", "log(or)"))
  expect_true(all(is.finite(pooled$estimates$estimate)))
  expect_true(all(is.finite(pooled$estimates$std.err)))

  # A finite pooled degrees of freedom is what says the Barnard-Rubin
  # adjustment ran rather than falling back to a normal reference.
  expect_true(all(is.finite(pooled$estimates$df)))
})

# The complete-data degrees of freedom are looked for on the results first and
# on their outcome models second. A balancing result reports none of its own,
# since the stacked system is not a fit with residual degrees of freedom, so
# this workflow always reaches the second place. Pinning the value it lands on
# rather than only its finiteness is what says the fallback ran: `NA` would mean
# the chain broke, and `Inf` would mean it gave up and warned instead.
test_that("the pooled degrees of freedom fall back to the outcome models", {
  skip_if_not_installed("mice")
  data <- ipw_pooling_fixture()
  imp <- withr::with_seed(4321, mice::mice(data, m = 3, print = FALSE))
  fits <- lapply(mice::complete(imp, "all"), fit_pooling_ipw)

  pooled <- pool_ipw(fits)

  result_df <- vapply(
    fits,
    function(fit) as.numeric(stats::df.residual(fit)),
    numeric(1)
  )
  outcome_df <- vapply(
    fits,
    function(fit) as.numeric(stats::df.residual(fit$outcome_mod)),
    numeric(1)
  )

  expect_true(all(is.na(result_df)))
  expect_false(is.na(pooled$dfcom))
  expect_false(is.infinite(pooled$dfcom))
  expect_equal(pooled$dfcom, min(outcome_df))
})

# Rubin's point estimate is the mean of the per-imputation ones. Recomputing it
# from the stored tables rather than from anything the pooled object carries is
# what makes this a check of the pooling arithmetic instead of a restatement
# of it.
test_that("the pooled estimate is the mean of the per-imputation estimates", {
  skip_if_not_installed("mice")
  data <- ipw_pooling_fixture()
  imp <- withr::with_seed(4321, mice::mice(data, m = 3, print = FALSE))
  fits <- lapply(mice::complete(imp, "all"), fit_pooling_ipw)

  pooled <- pool_ipw(fits)

  # The recompute reads each result's estimates by position, so the row
  # correspondence it assumes is pinned rather than trusted: every fit has to
  # report the same measures in the same order as the pooled table, which is
  # also the agreement the pooling itself requires.
  for (fit in fits) {
    expect_identical(fit$estimates$effect, pooled$estimates$effect)
  }

  for (row in seq_len(nrow(pooled$estimates))) {
    per_imputation <- vapply(
      fits,
      function(fit) fit$estimates$estimate[[row]],
      numeric(1)
    )
    expect_equal(pooled$estimates$estimate[[row]], mean(per_imputation))
  }
})

# One call pools both readings of a set of results and stores the one it was not
# asked for beside the one it was, so a pooled result moves between them
# afterwards. What makes the second reading reachable at all is balancing's:
# `ipw()` hands every outcome model over already wrapped, and the conditional
# reading is that wrapped model's coefficient surface under the corrected
# covariance. A set of results carrying no such surface pools on the marginal
# reading alone and refuses the flip, so pinning that a balancing pool is never
# that kind belongs here rather than upstream, where the wrapping is not in
# view.
test_that("a pooled result flips to the outcome models' coefficients", {
  skip_if_not_installed("mice")
  data <- ipw_pooling_fixture()
  imp <- withr::with_seed(4321, mice::mice(data, m = 3, print = FALSE))
  fits <- lapply(mice::complete(imp, "all"), fit_pooling_ipw)

  # Unqualified, for the reason `pool_ipw()` is: a caller who has attached
  # balancing alone flips a result that balancing's `ipw()` built.
  pooled <- pool_ipw(fits)
  flipped <- as_conditional(pooled)

  # Every imputation is analyzed with the same outcome formula, so one fit's
  # coefficient names are the set's, which the loop pins rather than assumes.
  # The literal is pinned beside them because a comparison of two names read the
  # same way would be satisfied by two empty vectors.
  coefficients <- names(stats::coef(fits[[1L]]$outcome_mod))
  expect_identical(coefficients, c("(Intercept)", "exposure"))
  for (fit in fits) {
    expect_identical(names(stats::coef(fit$outcome_mod)), coefficients)
  }

  expect_identical(flipped$effects, "conditional")
  expect_identical(flipped$estimates$effect, coefficients)

  # The flip exchanges two readings the pooling already computed rather than
  # computing either, so the result taken out to the other one and back is the
  # result that went in, the pooling diagnostics and the stored reading
  # included.
  expect_identical(as_marginal(flipped), pooled)
})

# Which reading a call pools actively says what the result presents rather than
# what was estimated, so the two routes to one pooled result have to arrive at
# the same object. The per-imputation results are the same objects on both
# routes, and everything the pooled result reports about the analyses rather
# than about a reading of them is read off those. The complete-data degrees of
# freedom are among them, which is why this leaves `dfcom` to be resolved on
# each route instead of naming it to hold it fixed.
test_that("pooling either reading of one set gives the same pooled result", {
  skip_if_not_installed("mice")
  data <- ipw_pooling_fixture()
  imp <- withr::with_seed(4321, mice::mice(data, m = 3, print = FALSE))
  fits <- lapply(mice::complete(imp, "all"), fit_pooling_ipw)

  pooled <- pool_ipw(fits)
  from_conditional <- pool_ipw(lapply(fits, as_conditional))

  # The narrower assertions come first so that a mismatch names the field or the
  # shape rather than the whole pooled result.
  expect_identical(from_conditional$effects, "conditional")
  expect_identical(names(from_conditional), names(pooled))
  expect_identical(from_conditional$dfcom, pooled$dfcom)
  expect_identical(from_conditional$m, pooled$m)

  expect_identical(as_marginal(from_conditional), pooled)
})

# The pooled accessors take a reading for one call, as the accessors on an
# unpooled result do. What is pinned here is shape alone: which rows each
# reading reports, under what names, and that naming one leaves the pooled
# result presenting the reading it stored. The numbers in either table are
# causalgenerics' to prove, and this file leaves them there.
test_that("the pooled accessors report either reading for one call", {
  skip_if_not_installed("mice")
  data <- ipw_pooling_fixture()
  imp <- withr::with_seed(4321, mice::mice(data, m = 3, print = FALSE))
  fits <- lapply(mice::complete(imp, "all"), fit_pooling_ipw)

  pooled <- pool_ipw(fits)
  coefficients <- names(stats::coef(fits[[1L]]$outcome_mod))

  expect_identical(
    names(stats::coef(pooled, effects = "conditional")),
    coefficients
  )

  # The tidied frame names its rows in a `term` column, which is the column both
  # readings of a pooled result use and the one an unpooled result uses too.
  conditional <- as.data.frame(pooled, effects = "conditional")
  expect_identical(conditional$term, coefficients)

  # Naming a reading answers in it and leaves the result where it was, so a
  # following call with nothing named answers in the stored one.
  expect_identical(pooled$effects, "marginal")
  expect_identical(names(stats::coef(pooled)), c("rd", "log(rr)", "log(or)"))
  expect_identical(as.data.frame(pooled)$term, c("rd", "log(rr)", "log(or)"))
})

# Only the estimating-equation methods carry the container `ipw()`
# differentiates, so a bootstrap-only method refuses before there is anything to
# pool. The refusal belongs to `ipw()` rather than to the pooling: a method
# without estimating equations never produces a result to put in the list, so
# there is no pooled surface for a spec to reach.
test_that("a bootstrap-only method refuses at the per-imputation step", {
  data <- ipw_pooling_fixture()
  complete <- data[!is.na(data$x1), , drop = FALSE]
  fit <- balance(
    complete,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  complete$.wts <- stats::weights(fit)
  outcome_mod <- stats::glm(
    y ~ exposure,
    data = complete,
    family = stats::quasibinomial(),
    weights = .wts
  )

  expect_false(supports_estimating_equations(bw_energy()))
  expect_error(
    ipw(fit, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  )
})

# ---- Pooling effects reported by subgroup ----------------------------------
#
# A `.by` result reports the effects over the whole sample, then within each
# subgroup of a modifier, then for each non-reference subgroup against the
# reference one, and names every row by the subgroup it belongs to. The pooling
# keys rows on that name along with the effect measure, so a subgroup's rows
# combine with the same subgroup's rows in the other analyses rather than with
# whichever rows happen to sit at the same position.
#
# The keying is causalgenerics', which reads the identity columns of the frames
# it is handed. What balancing owes is frames that carry the subgroup column at
# all, and carry the same subgroups in every analysis, so these run the whole
# path on real fits rather than restating the arithmetic.

# The same missingness as the ungrouped fixture, plus a modifier that is
# complete in every row. The modifier has to be observed rather than imputed:
# an imputed one would put a unit in different subgroups in different
# analyses, so the subgroups themselves would differ and there would be nothing
# to align.
ipw_pooling_by_fixture <- function(n = 250) {
  withr::with_seed(717, {
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    modifier <- factor(
      sample(c("lo", "hi"), n, replace = TRUE),
      levels = c("lo", "hi")
    )
    modifier_hi <- as.numeric(modifier == "hi")
    exposure <- stats::rbinom(
      n,
      1L,
      stats::plogis(0.6 * x1 - 0.4 * x2 + 0.5 * modifier_hi)
    )
    y <- stats::rbinom(
      n,
      1L,
      stats::plogis(
        -0.3 + 0.2 * exposure + 0.4 * x1 + 1.2 * exposure * modifier_hi
      )
    )
    # The indicator the fit balances on is derived in the analysis rather than
    # carried here. A column that is a deterministic function of another is
    # collinear with it, which the imputation reports as a logged event.
    data <- data.frame(
      exposure = exposure,
      x1 = x1,
      x2 = x2,
      modifier = modifier,
      y = y
    )
    missing <- stats::rbinom(n, 1L, stats::plogis(-1.2 + 0.5 * x2)) == 1L
    data$x1[missing] <- NA_real_
    data
  })
}

# The per-imputation analysis, the ungrouped one plus the modifier: balanced on
# it, interacted with the exposure in the outcome model, and named to `.by`.
fit_pooling_by_ipw <- function(data) {
  data$modifier_hi <- as.numeric(data$modifier == "hi")
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  data$.wts <- stats::weights(fit)
  ipw(
    fit,
    stats::glm(
      y ~ exposure * modifier,
      data = data,
      family = stats::quasibinomial(),
      weights = .wts
    ),
    .by = modifier
  )
}

test_that("pool_ipw() keys grouped balancing results by effect and subgroup", {
  skip_if_not_installed("mice")
  data <- ipw_pooling_by_fixture()
  # The modifier is observed, so it enters the imputation model as a predictor
  # and is never imputed itself.
  imp <- withr::with_seed(4321, mice::mice(data, m = 3, print = FALSE))
  fits <- lapply(mice::complete(imp, "all"), fit_pooling_by_ipw)

  pooled <- pool_ipw(fits)

  expect_s3_class(pooled, "ipw_pooled")
  expect_identical(pooled$m, 3L)
  expect_identical(names(pooled$estimates)[1:2], c("effect", "group"))
  expect_identical(nrow(pooled$estimates), 9L)
  expect_identical(
    unique(pooled$estimates$group),
    c(
      "overall",
      "modifier = lo",
      "modifier = hi",
      "modifier = hi vs modifier = lo"
    )
  )

  # Every analysis reports the same subgroups in the same order, which is both
  # what the pooling requires and what makes the pooled frame's key the key of
  # the frames it pooled.
  for (fit in fits) {
    expect_identical(fit$estimates$effect, pooled$estimates$effect)
    expect_identical(fit$estimates$group, pooled$estimates$group)
  }

  expect_true(all(is.finite(pooled$estimates$estimate)))
  expect_true(all(is.finite(pooled$estimates$std.err)))
  expect_true(all(is.finite(pooled$estimates$df)))

  # The pooled accessors label their rows by measure and subgroup together, the
  # way each analysis labels its own.
  labels <- paste(pooled$estimates$effect, pooled$estimates$group)
  expect_identical(names(stats::coef(pooled)), labels)
  expect_identical(as.data.frame(pooled)$group, pooled$estimates$group)
})

# Rubin's point estimate is the mean of the per-imputation ones, and it is a
# mean within a subgroup rather than across subgroups. Recomputing one cell of
# each block from the stored tables, after pinning that the row means the same
# subgroup in every analysis, is what says the alignment held.
test_that("pool_ipw() applies Rubin's rules within each subgroup", {
  skip_if_not_installed("mice")
  data <- ipw_pooling_by_fixture()
  imp <- withr::with_seed(4321, mice::mice(data, m = 3, print = FALSE))
  fits <- lapply(mice::complete(imp, "all"), fit_pooling_by_ipw)

  pooled <- pool_ipw(fits)

  cells <- list(
    c("rd", "modifier = hi"),
    c("log(rr)", "modifier = hi vs modifier = lo")
  )
  labels <- paste(pooled$estimates$effect, pooled$estimates$group)
  for (cell in cells) {
    label <- paste(cell[[1L]], cell[[2L]])
    row <- match(label, labels)

    # The cell has to be in the pooled frame before anything can be read at its
    # position, and its absence is the whole of what a missing subgroup column
    # would show, so the loop reports that and moves on rather than indexing
    # past the end of every frame it holds.
    expect_false(is.na(row), label = paste0("pooled row for ", label))
    if (is.na(row)) {
      next
    }

    per_imputation <- vapply(
      fits,
      function(fit) {
        expect_identical(fit$estimates$group[[row]], cell[[2L]])
        fit$estimates$estimate[[row]]
      },
      numeric(1)
    )
    expect_equal(pooled$estimates$estimate[[row]], mean(per_imputation))
  }
})

# Analyses whose subgroups differ have no common set of rows to average, and
# the pooling refuses rather than lining up whatever sits at each position. The
# refusal is causalgenerics', keyed on the row labels, and the labels carry the
# subgroup because balancing writes the subgroup column: a grouped result whose
# rows were keyed by the effect measure alone would look poolable against any
# other, so reaching that refusal from balancing results is the thing pinned.
#
# No imputation is involved, since the disagreement is between two analyses
# rather than between two completed datasets, and the fixture's complete rows
# are enough to build both.
test_that("pool_ipw() refuses grouped results whose subgroups disagree", {
  data <- ipw_pooling_by_fixture()
  complete <- data[!is.na(data$x1), , drop = FALSE]
  # The same analysis on the same rows, with one subgroup split off the top of
  # `x2`, so the two results agree about everything except which subgroups they
  # report.
  widened <- complete
  widened$modifier <- factor(
    ifelse(complete$x2 > 1, "top", as.character(complete$modifier)),
    levels = c("lo", "hi", "top")
  )

  narrow_result <- fit_pooling_by_ipw(complete)
  wide_result <- fit_pooling_by_ipw(widened)

  expect_identical(
    unique(narrow_result$estimates$group),
    c(
      "overall",
      "modifier = lo",
      "modifier = hi",
      "modifier = hi vs modifier = lo"
    )
  )
  expect_true("modifier = top" %in% wide_result$estimates$group)

  expect_error(
    pool_ipw(list(narrow_result, wide_result)),
    class = "causalgenerics_pool_mismatch_labels"
  )
})
