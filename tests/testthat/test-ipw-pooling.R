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
  imp <- mice::mice(data, m = 3, print = FALSE, seed = 4321)
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
  imp <- mice::mice(data, m = 3, print = FALSE, seed = 4321)
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
  imp <- mice::mice(data, m = 3, print = FALSE, seed = 4321)
  fits <- lapply(mice::complete(imp, "all"), fit_pooling_ipw)

  pooled <- pool_ipw(fits)

  for (row in seq_len(nrow(pooled$estimates))) {
    per_imputation <- vapply(
      fits,
      function(fit) fit$estimates$estimate[[row]],
      numeric(1)
    )
    expect_equal(pooled$estimates$estimate[[row]], mean(per_imputation))
  }
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
