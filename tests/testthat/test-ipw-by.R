# Specs for effect modification on balancing's `ipw()` method: the `.by`
# argument, which reports the effects again within the levels of a modifier and
# contrasts each level against the reference one. The rows a request adds are
# parameters of the same stacked system the whole-sample rows come from, so they
# carry standard errors that account for having estimated the weights and
# covariances that couple one subgroup to another.
#
# What these pin is the reported surface and the arithmetic behind it: which
# rows a grouped result carries and under what names, that the whole-sample rows
# are the rows the ungrouped fit already reported, that each subgroup's means are
# the g-computation means over that subgroup alone, that the covariance is one
# joint block rather than per-subgroup fits stitched together, and the
# configurations the argument refuses.

# ---- Fixtures --------------------------------------------------------------

# A binary-exposure fixture whose effect differs across the levels of a
# two-level modifier. The modifier confounds the exposure as well as modifying
# its effect, so a fit that balances it has real work to do, and it rides along
# as a numeric indicator, `modifier_hi`, because that is the parameterization
# the weight parameters stay identified in.
#
# The modifier declares its levels in reverse alphabetical order on purpose. The
# reference subgroup every contrast of subgroups is measured against is the
# modifier's first level, which is `"lo"` here and would be `"hi"` for an
# implementation that sorted the levels itself.
ipw_by_fixture <- function(n = 400) {
  withr::with_seed(808, {
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
      stats::plogis(0.7 * x1 - 0.5 * x2 + 0.6 * modifier_hi)
    )
    y <- stats::rbinom(
      n,
      1L,
      stats::plogis(
        -0.6 +
          0.2 * exposure +
          0.5 * x1 +
          0.3 * modifier_hi +
          1.4 * exposure * modifier_hi
      )
    )
    y_cont <- 1 +
      0.2 * exposure +
      0.5 * x1 -
      0.3 * x2 +
      1.2 * exposure * modifier_hi +
      stats::rnorm(n)
    data.frame(
      exposure = exposure,
      x1 = x1,
      x2 = x2,
      modifier = modifier,
      modifier_hi = modifier_hi,
      y = y,
      y_cont = y_cont
    )
  })
}

# A three-level categorical exposure crossed with the same two-level modifier.
# The exposure comes from the shared `sim_categorical()` process; the modifier
# and the outcome are drawn here under their own seed, with the interaction
# concentrated on the `"c"` level so the two subgroups disagree about one
# contrast and agree about the other.
ipw_by_categorical_fixture <- function(n = 450) {
  data <- sim_categorical(n)
  withr::with_seed(909, {
    data$modifier <- factor(
      sample(c("lo", "hi"), n, replace = TRUE),
      levels = c("lo", "hi")
    )
    data$modifier_hi <- as.numeric(data$modifier == "hi")
    is_b <- as.numeric(data$exposure == "b")
    is_c <- as.numeric(data$exposure == "c")
    linear_predictor <- -0.5 +
      0.3 * is_b +
      0.2 * is_c +
      0.5 * data$x1 +
      0.2 * data$modifier_hi +
      1.1 * is_b * data$modifier_hi +
      0.9 * is_c * data$modifier_hi
    data$y <- stats::rbinom(n, 1L, stats::plogis(linear_predictor))
  })
  data
}

# A weighted outcome model of the shape `ipw()` expects, with the weights riding
# along as a column so the model frame resolves them.
fit_by_outcome <- function(formula, data, wts, family) {
  data[[".wts"]] <- wts
  suppressWarnings(
    stats::glm(formula, data = data, family = family, weights = .wts)
  )
}

# The marginal means a subgroup's rows report: the outcome model predicted with
# the exposure fixed to each level in turn, averaged over the units of that
# subgroup alone. `stratum` is the subgroup's indicator and `tilt` the
# standardization weight the whole-sample means already use, which is uniform
# for a pooled estimand and the focal group's indicator for a focal one. The
# oracle reads the fitted model and the data rather than anything the result
# carries, so it anchors the subgroup blocks independently of the stack that
# reports them.
by_stratum_means <- function(
  outcome_mod,
  data,
  stratum,
  exposure_name = "exposure",
  tilt = NULL
) {
  values <- data[[exposure_name]]
  levels <- if (is.factor(values)) levels(values) else sort(unique(values))
  weight <- (tilt %||% rep(1, nrow(data))) * as.numeric(stratum)
  vapply(
    levels,
    function(level) {
      counterfactual <- data
      counterfactual[[exposure_name]] <- if (is.factor(values)) {
        factor(level, levels = levels)
      } else {
        level
      }
      stats::weighted.mean(
        stats::predict(
          outcome_mod,
          newdata = counterfactual,
          type = "response"
        ),
        weight
      )
    },
    numeric(1)
  )
}

# One row's estimate, read out of a coerced result by the columns that name it.
# The lookup says how many rows it matched before it returns: a key matching
# nothing gives an empty vector, and a comparison between two of those holds,
# so a reader that did not report the count could report agreement between two
# rows that are not there.
by_estimate <- function(estimates, effect, group, contrast = NULL) {
  rows <- estimates$term == effect & estimates$group == group
  if (!is.null(contrast)) {
    rows <- rows & estimates$contrast == contrast
  }
  value <- estimates$estimate[rows]
  testthat::expect_identical(
    length(value),
    1L,
    label = paste(c(effect, contrast, group), collapse = " ")
  )
  value
}

# ---- Row identity ----------------------------------------------------------

# The column is added by a request rather than carried by every result. An
# ungrouped fit reports one set of effects over the whole sample, so a column
# naming the subgroup each row belongs to would repeat a single value down the
# table and read as a subgroup that was named.

test_that("an ungrouped binary fit names no subgroups", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  estimates <- ipw(fit, outcome_mod)$estimates

  expect_false("group" %in% names(estimates))
  expect_identical(
    estimates$effect,
    c("mean", "mean", "rd", "log(rr)", "log(or)")
  )
})

test_that("a .by fit reports the whole sample, each stratum, then their contrast", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  # The stored table, for the reason every other column contract in this suite
  # reads it: the coerced frame is causalgenerics' contract, and this is the one
  # balancing fills in.
  estimates <- ipw(fit, outcome_mod, .by = modifier)$estimates

  # The subgroup column sits after the contrast, which is where the shared
  # contract places it: a contrast names the level a row belongs to or the pair
  # it compares, and the subgroup qualifies the whole of that.
  expect_named(
    estimates,
    c(
      "effect",
      "contrast",
      "group",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
  expect_identical(
    estimates$effect,
    c(
      "mean",
      "mean",
      "rd",
      "log(rr)",
      "log(or)",
      rep("mean", 4L),
      rep(c("rd", "log(rr)"), times = 3)
    )
  )
  expect_identical(
    estimates$contrast,
    c(
      "0",
      "1",
      rep("1 vs 0", 3L),
      rep(c("0", "1"), times = 2),
      rep("1 vs 0", 6L)
    )
  )
  expect_identical(
    estimates$group,
    c(
      rep("overall", 5),
      rep("modifier = lo", 2),
      rep("modifier = hi", 2),
      rep("modifier = lo", 2),
      rep("modifier = hi", 2),
      rep("modifier = hi vs modifier = lo", 2)
    )
  )
  expect_identical(nrow(estimates), 15L)
  expect_true(all(is.finite(estimates$estimate)))

  # A block of means is never split by the contrasts built from it: the
  # whole-sample pair leads the table and the stratum pairs sit together after
  # the whole-sample contrasts and ahead of the stratum ones.
  expect_identical(
    which(estimates$effect == "mean"),
    c(1L, 2L, 6L, 7L, 8L, 9L)
  )
})

# An odds ratio is noncollapsible: the odds ratio over a whole sample is not an
# average of the odds ratios within its subgroups, and the difference of two of
# them is not the difference in effect it reads as. The whole-sample rows keep
# it, since nothing there averages anything, and the subgroup rows report the
# collapsible measures alone.

test_that("a .by fit reports no odds ratio outside its whole-sample rows", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  estimates <- ipw(fit, outcome_mod, .by = modifier)$estimates
  odds <- estimates$effect == "log(or)"

  expect_identical(sum(odds), 1L)
  expect_identical(estimates$group[odds], "overall")
  expect_identical(
    sort(unique(estimates$effect[estimates$group != "overall"])),
    c("log(rr)", "mean", "rd")
  )
})

# The blocks a request appends come after every block the ungrouped fit already
# solves, and no earlier equation reads a parameter of theirs, so the leading
# block of the sandwich is the one that fit produces. That is what lets a
# grouped result report the whole-sample rows it grew out of rather than
# recomputing them at a different variance.

test_that("a .by fit leaves the whole-sample rows it already reported alone", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  ungrouped <- ipw(fit, outcome_mod)$estimates
  grouped <- ipw(fit, outcome_mod, .by = modifier)$estimates
  overall <- grouped[grouped$group == "overall", , drop = FALSE]

  expect_identical(overall$effect, ungrouped$effect)
  expect_equal(overall$estimate, ungrouped$estimate, tolerance = 1e-10)
  expect_equal(overall$std.err, ungrouped$std.err, tolerance = 1e-10)
})

# `.by = NULL` is the default written out, so it has to be the absence of a
# request rather than a request naming nothing. The whole stored frame is
# compared, the covariance it carries included, since an implementation that
# added an empty subgroup block would still agree on the numbers.

test_that(".by = NULL reports the frame an ungrouped fit reports", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  baseline <- ipw(fit, outcome_mod)
  explicit <- ipw(fit, outcome_mod, .by = NULL)

  expect_identical(explicit$estimates, baseline$estimates)
  expect_identical(names(explicit$fit$theta), names(baseline$fit$theta))
})

# ---- The subgroup arithmetic -----------------------------------------------

test_that("a .by fit standardizes each stratum's means over that stratum", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  # An adjusted model, so the subgroup means are a real g-computation average
  # rather than the pair of weighted cell means a saturated model collapses to.
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier + x1,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod, .by = modifier)

  for (level in levels(data$modifier)) {
    means <- by_stratum_means(outcome_mod, data, data$modifier == level)
    expect_equal(
      unname(result$fit$theta[[paste0("mu0_modifier = ", level)]]),
      unname(means[[1L]]),
      tolerance = 1e-8
    )
    expect_equal(
      unname(result$fit$theta[[paste0("mu1_modifier = ", level)]]),
      unname(means[[2L]]),
      tolerance = 1e-8
    )
  }

  estimates <- as.data.frame(result)
  for (level in levels(data$modifier)) {
    means <- by_stratum_means(outcome_mod, data, data$modifier == level)
    group <- paste0("modifier = ", level)
    expect_equal(
      by_estimate(estimates, "rd", group),
      means[[2L]] - means[[1L]],
      tolerance = 1e-8
    )
    expect_equal(
      by_estimate(estimates, "log(rr)", group),
      log(means[[2L]]) - log(means[[1L]]),
      tolerance = 1e-8
    )
  }
})

# A model saturated in the exposure and the modifier predicts one value per
# cell, so its subgroup g-computation means reduce to the weighted outcome means
# of the four cells. That is the reading a caller can compute by hand from the
# weights alone, and it anchors the subgroup block against an arithmetic that
# borrows nothing from the outcome model at all.

test_that("saturated stratum means are that stratum's weighted group means", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod, .by = modifier)

  for (level in levels(data$modifier)) {
    for (exposed in c(0L, 1L)) {
      cell <- data$modifier == level & data$exposure == exposed
      expect_equal(
        unname(result$fit$theta[[
          paste0("mu", exposed, "_modifier = ", level)
        ]]),
        stats::weighted.mean(data$y[cell], w[cell]),
        tolerance = 1e-8
      )
    }
  }
})

test_that("a .by fit contrasts each stratum against the reference stratum", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier + x1,
    data,
    w,
    stats::binomial()
  )

  estimates <- as.data.frame(ipw(fit, outcome_mod, .by = modifier))

  # The reference stratum is the modifier's first level, which this fixture
  # declares out of alphabetical order, so a sorted implementation would report
  # the contrast the other way around and carry the opposite sign.
  for (effect in c("rd", "log(rr)")) {
    expect_equal(
      by_estimate(estimates, effect, "modifier = hi vs modifier = lo"),
      by_estimate(estimates, effect, "modifier = hi") -
        by_estimate(estimates, effect, "modifier = lo"),
      tolerance = 1e-10
    )
  }

  # The modification the fixture was drawn with is large, so the contrast is
  # not a rounding artifact of two subgroup effects that agree.
  expect_gt(
    abs(by_estimate(estimates, "rd", "modifier = hi vs modifier = lo")),
    0.1
  )
})

test_that("a .by fit on a continuous outcome reports one difference per stratum", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y_cont ~ exposure * modifier,
    data,
    w,
    stats::gaussian()
  )

  estimates <- ipw(fit, outcome_mod, .by = modifier)$estimates

  expect_identical(
    estimates$effect,
    c("mean", "mean", "diff", rep("mean", 4L), rep("diff", 3L))
  )
  expect_identical(
    estimates$group,
    c(
      rep("overall", 3L),
      rep("modifier = lo", 2L),
      rep("modifier = hi", 2L),
      "modifier = lo",
      "modifier = hi",
      "modifier = hi vs modifier = lo"
    )
  )

  # A saturated linear model puts the difference in differences in one of its
  # own coefficients, so the contrast of subgroups is a number the fitted model
  # already reports and the row can be checked against it exactly.
  interaction <- stats::coef(outcome_mod)[["exposure:modifierhi"]]
  contrast <- estimates$estimate[
    estimates$group == "modifier = hi vs modifier = lo"
  ]
  expect_equal(contrast, interaction, tolerance = 1e-8)
})

# A focal estimand standardizes over the focal group rather than over everyone,
# and a subgroup's means standardize over the focal units inside that subgroup.
# The tilt therefore has to reach the subgroup rows as well as the whole-sample
# ones: an implementation that averaged over every unit of the subgroup would
# report the pooled effect within it while every other assertion still held.

test_that("a .by att fit standardizes each stratum over its treated units", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "att"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier + x1,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod, .by = modifier)

  for (level in levels(data$modifier)) {
    means <- by_stratum_means(
      outcome_mod,
      data,
      data$modifier == level,
      tilt = as.numeric(data$exposure == 1L)
    )
    expect_equal(
      unname(result$fit$theta[[paste0("mu0_modifier = ", level)]]),
      unname(means[[1L]]),
      tolerance = 1e-8
    )
    expect_equal(
      unname(result$fit$theta[[paste0("mu1_modifier = ", level)]]),
      unname(means[[2L]]),
      tolerance = 1e-8
    )
  }
})

# A factor carries its own level order and the reference stratum is the first of
# them. A character column carries none, so it is read as a factor on the way in
# and the reference stratum is whichever value sorts first. That is a rule about
# the column's type rather than about the data, and the two readings disagree
# whenever the values appear in an order other than their sorted one, so the
# fixture below makes them disagree.

test_that(".by measures a character modifier against its first sorted value", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))

  # The modifier relabelled so that the value appearing first sorts last. Which
  # of the two labels lands on the first row is read off the fixture rather than
  # assumed, so the disagreement is built rather than hoped for.
  first_seen <- as.character(data$modifier)[[1L]]
  labelled <- data
  labelled$letters <- ifelse(data$modifier == first_seen, "zebra", "alpha")
  expect_identical(labelled$letters[[1L]], "zebra")

  outcome_mod <- fit_by_outcome(
    y ~ exposure * letters,
    labelled,
    w,
    stats::binomial()
  )

  result <- expect_no_warning(ipw(fit, outcome_mod, .by = letters))
  estimates <- result$estimates

  expect_identical(
    unique(estimates$group),
    c(
      "overall",
      "letters = alpha",
      "letters = zebra",
      "letters = zebra vs letters = alpha"
    )
  )

  # The contrast of strata is the non-reference stratum minus the reference one,
  # so a reference read off first appearance rather than off sorted order would
  # carry the opposite sign.
  coerced <- as.data.frame(result)
  expect_equal(
    by_estimate(coerced, "rd", "letters = zebra vs letters = alpha"),
    by_estimate(coerced, "rd", "letters = zebra") -
      by_estimate(coerced, "rd", "letters = alpha"),
    tolerance = 1e-10
  )
})

# ---- The stacked system ----------------------------------------------------

test_that("a .by fit appends a mean and a contrast block for every stratum", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod, .by = modifier)

  parameters <- length(estimating_equations(fit)@parameters)
  coefficients <- length(stats::coef(outcome_mod))
  strata <- nlevels(data$modifier)
  # The widths come from the fit rather than from literals so that the
  # accounting reads as the rule it is. The weight parameters and the outcome
  # model's coefficients lead, as they do without a request; then the two
  # whole-sample means and the three measures contrasting them; then, per
  # stratum, a mean at each exposure level and the two collapsible measures;
  # then those measures once more for each non-reference stratum against the
  # reference one.
  expected <- parameters +
    coefficients +
    2L +
    3L +
    strata * 2L +
    strata * 2L +
    (strata - 1L) * 2L
  expect_identical(length(result$fit$theta), expected)

  expect_identical(
    utils::tail(names(result$fit$theta), 15L),
    c(
      "mu0",
      "mu1",
      "rd",
      "log(rr)",
      "log(or)",
      "mu0_modifier = lo",
      "mu1_modifier = lo",
      "mu0_modifier = hi",
      "mu1_modifier = hi",
      "rd_modifier = lo",
      "log(rr)_modifier = lo",
      "rd_modifier = hi",
      "log(rr)_modifier = hi",
      "rd_modifier = hi vs modifier = lo",
      "log(rr)_modifier = hi vs modifier = lo"
    )
  )
  expect_identical(
    dimnames(result$fit$vcov),
    list(names(result$fit$theta), names(result$fit$theta))
  )

  # Every reported row is read off the stack at its own name, so the estimates
  # table and the variance system cannot drift apart. A mean row is keyed by the
  # marginal-mean parameter of the level it belongs to, a contrast row by the
  # contrast itself, and either is suffixed with its subgroup outside the
  # whole-sample block.
  estimates <- as.data.frame(result)
  whole_sample <- ifelse(
    estimates$term == "mean",
    paste0("mu", sub(" vs .*$", "", estimates$contrast)),
    estimates$term
  )
  keys <- ifelse(
    estimates$group == "overall",
    whole_sample,
    paste0(whole_sample, "_", estimates$group)
  )
  expect_equal(
    estimates$estimate,
    unname(result$fit$theta[keys]),
    tolerance = 1e-12
  )
  expect_equal(
    estimates$std.error,
    unname(sqrt(diag(result$fit$vcov))[keys]),
    tolerance = 1e-12
  )
})

test_that("a .by fit reports a usable standard error for every row", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier + x1,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod, .by = modifier)
  estimates <- result$estimates

  expect_identical(nrow(estimates), 15L)
  expect_true(all(is.finite(estimates$std.err)))
  expect_true(all(estimates$std.err > 0))
  expect_true(all(estimates$ci.lower < estimates$estimate))
  expect_true(all(estimates$ci.upper > estimates$estimate))
  expect_equal(
    unname(sqrt(diag(stats::vcov(result)))),
    estimates$std.err,
    tolerance = 1e-12
  )
})

# The subgroups share the weight parameters and the outcome model's
# coefficients, so their effects covary. Per-subgroup fits assembled into a
# block-diagonal matrix would report exact zeros in the cross-subgroup entries
# while agreeing with themselves everywhere on the diagonal, so the coupling is
# what separates one joint system from that assembly.

test_that("a .by fit's covariance couples the subgroups it reports", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier + x1,
    data,
    w,
    stats::binomial()
  )

  covariance <- stats::vcov(ipw(fit, outcome_mod, .by = modifier))

  couples <- list(
    c("rd 1 vs 0 modifier = lo", "rd 1 vs 0 modifier = hi"),
    c("log(rr) 1 vs 0 modifier = lo", "log(rr) 1 vs 0 modifier = hi"),
    c("rd 1 vs 0 overall", "rd 1 vs 0 modifier = hi"),
    c(
      "rd 1 vs 0 modifier = hi",
      "rd 1 vs 0 modifier = hi vs modifier = lo"
    )
  )
  for (pair in couples) {
    entry <- covariance[pair[[1L]], pair[[2L]]]
    expect_true(is.finite(entry), label = paste(pair, collapse = " with "))
    expect_gt(abs(entry), 1e-8, label = paste(pair, collapse = " with "))
  }
  expect_equal(covariance, t(covariance), tolerance = 1e-12)

  # The contrast of two subgroups is a parameter of the same system, so its
  # variance is the joint variance of the difference. A stitched assembly would
  # have to read the sum of the two variances instead, which the coupling above
  # makes a different number.
  #
  # The identity is exact algebra, but every entry it combines comes out of the
  # finite-difference bread inversion, so the comparison cancels independently
  # rounded numbers whose last bits follow whichever BLAS and compiler the
  # platform provides. Release Linux and Windows builds have fallen outside
  # 1e-10 relative where this machine reads 2e-16, so every assertion of this
  # shape in the suite reads at 1e-6. That still sits three orders below the
  # smallest gap to the stitched sum any of them refuses, which is 1.5e-3.
  variance_lo <- covariance[
    "rd 1 vs 0 modifier = lo",
    "rd 1 vs 0 modifier = lo"
  ]
  variance_hi <- covariance[
    "rd 1 vs 0 modifier = hi",
    "rd 1 vs 0 modifier = hi"
  ]
  coupling <- covariance[
    "rd 1 vs 0 modifier = lo",
    "rd 1 vs 0 modifier = hi"
  ]
  contrast <- covariance[
    "rd 1 vs 0 modifier = hi vs modifier = lo",
    "rd 1 vs 0 modifier = hi vs modifier = lo"
  ]
  expect_equal(
    contrast,
    variance_lo + variance_hi - 2 * coupling,
    tolerance = 1e-6
  )
  expect_false(isTRUE(all.equal(contrast, variance_lo + variance_hi)))
})

# A focal estimand standardizes each stratum over its focal units alone, and the
# tilt that says so enters the mean rows of the stacked system rather than being
# applied to them afterwards. So it reaches the covariance as well as the
# estimates: an implementation that tilted the seed values and left the rows
# pooled would report the focal means correctly and a covariance belonging to a
# different pair of estimators, which the estimate assertions above could not
# see. Reading the same coupling and the same variance identity off an att fit
# is what pins the tilt into the sandwich.

test_that("a .by att fit couples its subgroups through the focal tilt", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "att"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier + x1,
    data,
    w,
    stats::binomial()
  )

  covariance <- stats::vcov(ipw(fit, outcome_mod, .by = modifier))

  couples <- list(
    c("rd 1 vs 0 modifier = lo", "rd 1 vs 0 modifier = hi"),
    c("rd 1 vs 0 overall", "rd 1 vs 0 modifier = hi"),
    c(
      "rd 1 vs 0 modifier = hi",
      "rd 1 vs 0 modifier = hi vs modifier = lo"
    )
  )
  for (pair in couples) {
    entry <- covariance[pair[[1L]], pair[[2L]]]
    expect_true(is.finite(entry), label = paste(pair, collapse = " with "))
    expect_gt(abs(entry), 1e-8, label = paste(pair, collapse = " with "))
  }

  variance_lo <- covariance[
    "rd 1 vs 0 modifier = lo",
    "rd 1 vs 0 modifier = lo"
  ]
  variance_hi <- covariance[
    "rd 1 vs 0 modifier = hi",
    "rd 1 vs 0 modifier = hi"
  ]
  coupling <- covariance[
    "rd 1 vs 0 modifier = lo",
    "rd 1 vs 0 modifier = hi"
  ]
  contrast <- covariance[
    "rd 1 vs 0 modifier = hi vs modifier = lo",
    "rd 1 vs 0 modifier = hi vs modifier = lo"
  ]
  # The identity is exact in the system and approximate in the sandwich, since
  # the row carrying it is differenced rather than written down, and a focal
  # estimand standardizes each stratum over a quarter of the sample, which
  # leaves the mean block's diagonal smaller and the difference less accurate.
  # It reads 1e-10 relative here against the 1e-16 a pooled estimand reaches.
  # Both are asserted at the tolerance the whole class carries, which is sized
  # for the spread across platforms rather than for either measurement.
  expect_equal(
    contrast,
    variance_lo + variance_hi - 2 * coupling,
    tolerance = 1e-6
  )
  expect_false(isTRUE(all.equal(contrast, variance_lo + variance_hi)))

  # The att covariance is its own, not the ate one relabelled, so the tilt is
  # doing work here rather than cancelling.
  pooled_fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  pooled_w <- as.numeric(stats::weights(pooled_fit))
  pooled_mod <- fit_by_outcome(
    y ~ exposure * modifier + x1,
    data,
    pooled_w,
    stats::binomial()
  )
  pooled <- stats::vcov(ipw(pooled_fit, pooled_mod, .by = modifier))
  expect_false(isTRUE(all.equal(
    covariance["rd 1 vs 0 modifier = lo", "rd 1 vs 0 modifier = hi"],
    pooled["rd 1 vs 0 modifier = lo", "rd 1 vs 0 modifier = hi"]
  )))
})

# The six blocks the stack is assembled from are only all present under a
# request: an ungrouped fit leaves the last two absent. The assembly of those
# blocks into the stacked matrix claims agreement with `rbind()` to the bit, and
# this is the route that states the claim over a full stack. What the
# expectation compares is the assembled matrix itself, at every evaluation the
# finite difference asks the closure for.

test_that("a .by fit stacks its psi blocks as rbind would", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  expect_stacked_psi_matches_rbind(ipw(fit, outcome_mod, .by = modifier))
})

# ---- Labels ----------------------------------------------------------------

# The measure repeats across subgroups and names no row on its own, so every
# label a grouped result carries joins the measure to the subgroup. balancing
# builds those labels for the covariance it attaches to the estimates table, and
# causalgenerics builds them again for the accessors and the printed table. Both
# have to grow the subgroup segment together: a label built from the measure
# alone leaves two measures standing for nine rows, which the covariance carries
# as duplicated dimnames rather than as an error.

test_that("a .by fit labels its coefficients, covariance, and printed rows alike", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod, .by = modifier)
  labels <- c(
    "mean 0 overall",
    "mean 1 overall",
    "rd 1 vs 0 overall",
    "log(rr) 1 vs 0 overall",
    "log(or) 1 vs 0 overall",
    "mean 0 modifier = lo",
    "mean 1 modifier = lo",
    "mean 0 modifier = hi",
    "mean 1 modifier = hi",
    "rd 1 vs 0 modifier = lo",
    "log(rr) 1 vs 0 modifier = lo",
    "rd 1 vs 0 modifier = hi",
    "log(rr) 1 vs 0 modifier = hi",
    "rd 1 vs 0 modifier = hi vs modifier = lo",
    "log(rr) 1 vs 0 modifier = hi vs modifier = lo"
  )

  expect_identical(anyDuplicated(labels), 0L)
  expect_identical(
    dimnames(attr(result$estimates, "ipw_vcov", exact = TRUE)),
    list(labels, labels)
  )
  expect_identical(names(stats::coef(result)), labels)
  expect_identical(dimnames(stats::vcov(result)), list(labels, labels))
  expect_identical(rownames(stats::confint(result)), labels)

  output <- utils::capture.output(print(result))
  for (label in labels) {
    expect_true(
      any(grepl(label, output, fixed = TRUE)),
      label = paste0("printed row ", label)
    )
  }
})

test_that("a coerced .by result heads its subgroup column after the term", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  coerced <- as.data.frame(ipw(fit, outcome_mod, .by = modifier))

  expect_identical(names(coerced)[1:3], c("term", "contrast", "group"))
  expect_identical(
    coerced$group,
    c(
      rep("overall", 5),
      rep("modifier = lo", 2),
      rep("modifier = hi", 2),
      rep("modifier = lo", 2),
      rep("modifier = hi", 2),
      rep("modifier = hi vs modifier = lo", 2)
    )
  )
})

# ---- The conditional reading -----------------------------------------------

# A request reports the marginal effects within subgroups. The conditional
# reading is the outcome model's own coefficient surface under the corrected
# covariance, and a model fitted once has one such surface however the marginal
# side is broken up, so `.by` must leave it exactly where an ungrouped result
# leaves it. The reading is read through the accessors rather than through the
# coerced frame, which reports the stored marginal table whichever reading a
# result records.

test_that("a .by result's conditional reading is the outcome model's surface", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod, .by = modifier)
  conditional <- causalgenerics::as_conditional(result)
  coefficients <- names(stats::coef(outcome_mod))

  expect_identical(conditional$effects, "conditional")
  expect_identical(stats::coef(conditional), stats::coef(outcome_mod))
  expect_identical(
    dimnames(stats::vcov(conditional)),
    list(coefficients, coefficients)
  )
  expect_identical(rownames(stats::confint(conditional)), coefficients)

  # No subgroup reaches the coefficient names, so nothing there carries the
  # `"var = value"` segment the marginal rows are keyed by.
  expect_false(any(grepl("modifier = ", coefficients, fixed = TRUE)))

  # The corrected covariance is the block the result already carries on the
  # stored outcome model, which is what says the reading reports the same
  # numbers a grouped result computed rather than recomputing anything.
  expect_identical(stats::vcov(conditional), stats::vcov(result$outcome_mod))

  # The flip exchanges two readings the result already holds, so a grouped
  # result taken out to the other reading and back is the result that went in.
  expect_identical(causalgenerics::as_marginal(conditional), result)
})

# ---- Categorical exposures -------------------------------------------------

test_that("a .by categorical fit crosses its contrasts with its subgroups", {
  data <- ipw_by_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  estimates <- ipw(fit, outcome_mod, .by = modifier)$estimates

  # The subgroup column follows the contrast column, since the contrast names
  # one side of the comparison and the subgroup qualifies the whole of it.
  expect_named(
    estimates,
    c(
      "effect",
      "contrast",
      "group",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )

  # Blocks run subgroup-major; inside a block the ordering is the one an
  # ungrouped categorical fit uses, contrast-major with the measures of one
  # contrast together.
  contrasts <- c("b vs a", "c vs a")
  overall <- c("rd", "log(rr)", "log(or)")
  stratum <- c("rd", "log(rr)")
  groups <- c(
    "modifier = lo",
    "modifier = hi",
    "modifier = hi vs modifier = lo"
  )

  levels <- c("a", "b", "c")
  strata <- c("modifier = lo", "modifier = hi")

  expect_identical(
    estimates$effect,
    c(
      rep("mean", length(levels)),
      rep(overall, times = length(contrasts)),
      rep("mean", length(levels) * length(strata)),
      rep(rep(stratum, times = length(contrasts)), times = length(groups))
    )
  )
  expect_identical(
    estimates$contrast,
    c(
      levels,
      rep(contrasts, each = length(overall)),
      rep(levels, times = length(strata)),
      rep(rep(contrasts, each = length(stratum)), times = length(groups))
    )
  )
  expect_identical(
    estimates$group,
    c(
      rep("overall", length(levels) + length(overall) * length(contrasts)),
      rep(strata, each = length(levels)),
      rep(groups, each = length(stratum) * length(contrasts))
    )
  )
  expect_identical(nrow(estimates), 27L)
})

test_that("a .by categorical fit keeps its odds ratios among the whole-sample rows", {
  data <- ipw_by_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  estimates <- ipw(fit, outcome_mod, .by = modifier)$estimates
  odds <- estimates$effect == "log(or)"

  expect_identical(sum(odds), 2L)
  expect_identical(estimates$group[odds], rep("overall", 2L))
  expect_identical(estimates$contrast[odds], c("b vs a", "c vs a"))
})

test_that("a .by categorical fit reports each contrast within each subgroup", {
  data <- ipw_by_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier + x1,
    data,
    w,
    stats::binomial()
  )

  estimates <- as.data.frame(ipw(fit, outcome_mod, .by = modifier))

  for (level in levels(data$modifier)) {
    means <- by_stratum_means(outcome_mod, data, data$modifier == level)
    group <- paste0("modifier = ", level)
    for (compared in c("b", "c")) {
      contrast <- paste0(compared, " vs a")
      expect_equal(
        by_estimate(estimates, "rd", group, contrast),
        means[[compared]] - means[["a"]],
        tolerance = 1e-8
      )
      expect_equal(
        by_estimate(estimates, "log(rr)", group, contrast),
        log(means[[compared]]) - log(means[["a"]]),
        tolerance = 1e-8
      )
    }
  }

  # Each subgroup contrast is measured against the reference subgroup within
  # the exposure contrast it belongs to, so the crossing carries both
  # references at once.
  for (compared in c("b", "c")) {
    contrast <- paste0(compared, " vs a")
    expect_equal(
      by_estimate(
        estimates,
        "rd",
        "modifier = hi vs modifier = lo",
        contrast
      ),
      by_estimate(estimates, "rd", "modifier = hi", contrast) -
        by_estimate(estimates, "rd", "modifier = lo", contrast),
      tolerance = 1e-10
    )
  }
})

test_that("a .by categorical fit labels its rows by measure, contrast, and subgroup", {
  data <- ipw_by_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  result <- ipw(fit, outcome_mod, .by = modifier)
  estimates <- result$estimates
  labels <- paste(estimates$effect, estimates$contrast, estimates$group)

  expect_identical(length(labels), 27L)
  expect_identical(anyDuplicated(labels), 0L)
  expect_identical(
    dimnames(attr(estimates, "ipw_vcov", exact = TRUE)),
    list(labels, labels)
  )
  expect_identical(names(stats::coef(result)), labels)
  expect_identical(rownames(stats::confint(result)), labels)
})

# ---- Subgroups no unit belongs to ------------------------------------------

# A factor may declare a level nothing carries, and such a level names an empty
# subgroup: no means to standardize, no contrast to report, and a row of missing
# values if one were reported anyway. It is dropped rather than refused, which
# is propensity's behavior and the reading that lets a modifier subset without
# being recoded first. The declaration reaches the argument only through
# `.data`, since a model frame drops a factor's unused levels on the way in.

test_that(".by drops a modifier level no unit carries", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )
  supplied <- data
  supplied$modifier <- factor(
    as.character(data$modifier),
    levels = c("lo", "hi", "none")
  )
  expect_identical(sum(supplied$modifier == "none"), 0L)

  # The outcome model still carries the interaction, since the supplied column
  # is the modifier the model was fitted on under one more declared level, so
  # nothing here is warned about.
  result <- expect_no_warning(
    ipw(fit, outcome_mod, .data = supplied, .by = modifier)
  )
  estimates <- result$estimates

  expect_identical(nrow(estimates), 15L)
  expect_identical(
    unique(estimates$group),
    c(
      "overall",
      "modifier = lo",
      "modifier = hi",
      "modifier = hi vs modifier = lo"
    )
  )
  expect_false(any(grepl("none", estimates$group, fixed = TRUE)))
})

# ---- The rows a supplied frame holds ---------------------------------------

# `.by` makes `.data` the channel the modifier arrives on, so what that frame
# holds decides which units each subgroup is built from. Half of the stacked
# system reads it and half reads the fit's own order: the counterfactual
# predictions, the stratum indicators, and a focal estimand's standardization
# come from `.data`, while the weight equations, the outcome-model score, and
# every cross term of the meat come from the fit. A frame carrying the fit's
# rows in another order pairs each unit's prediction with another unit's weight,
# which a weighted mean is blind to and a sandwich is not: the effects come back
# unchanged and the standard errors are wrong. Nothing downstream can notice, so
# the frame is checked against the outcome model's own, which the weight
# preflight has already pinned to the fit's order.

test_that(".data must hold the fit's rows in the fit's order", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )
  # Sorted by the modifier, which is the reordering a caller is most likely to
  # arrive at and the one that moves the standard errors furthest, since it
  # groups the units each stratum indicator selects.
  sorted <- data[order(data$modifier), , drop = FALSE]

  expect_error(
    ipw(fit, outcome_mod, .data = sorted, .by = modifier),
    class = "balancing_ipw_input_error"
  )

  # The mechanism is not `.by`'s. An ungrouped result reads the same frame for
  # its counterfactual designs, so the refusal covers that call as well.
  expect_error(
    ipw(fit, outcome_mod, .data = sorted),
    class = "balancing_ipw_input_error"
  )
})

# Reordering is one way a frame stops describing the fit's rows and replacing a
# value is the other. A modifier edited in one row moves that unit between
# strata while every count and every level stays as it was, so the check has to
# compare values rather than shapes.

test_that(".data must hold the values the models were fitted on", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )
  edited <- data
  first_lo <- which(data$modifier == "lo")[[1L]]
  edited$modifier[[first_lo]] <- "hi"

  expect_identical(nrow(edited), nrow(data))
  expect_identical(levels(edited$modifier), levels(data$modifier))

  expect_error(
    ipw(fit, outcome_mod, .data = edited, .by = modifier),
    class = "balancing_ipw_input_error"
  )
})

# The check compares the columns the two frames name in common, so a frame
# carrying columns the outcome model never saw is still the frame it was fitted
# on. That is the shape every `.by` workflow supplies, since the modifier
# arrives beside the model's own variables.

test_that(".data carrying columns the model never saw is still aligned", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )
  supplied <- data
  supplied$unused <- withr::with_seed(31, stats::rnorm(nrow(data)))

  supplied_result <- ipw(fit, outcome_mod, .data = supplied, .by = modifier)
  frame_result <- ipw(fit, outcome_mod, .by = modifier)

  expect_identical(nrow(supplied_result$estimates), 15L)
  expect_identical(
    supplied_result$estimates$group,
    frame_result$estimates$group
  )
  expect_equal(
    supplied_result$estimates$estimate,
    frame_result$estimates$estimate,
    tolerance = 1e-12
  )
  expect_equal(
    supplied_result$estimates$std.err,
    frame_result$estimates$std.err,
    tolerance = 1e-12
  )
})

# ---- Refusals --------------------------------------------------------------

# A continuous exposure reports a coefficient of its marginal structural model
# rather than a contrast of standardized means, so there is no subgroup effect
# for the argument to name. Fitting the model on each subgroup would report a
# coefficient apiece and no covariance between them, which is the workflow the
# refusal has to point at rather than approximate.

test_that(".by refuses a continuous exposure", {
  data <- ipw_by_fixture()
  data$dose <- withr::with_seed(
    12,
    0.5 * data$x1 - 0.3 * data$x2 + stats::rnorm(nrow(data))
  )
  fit <- balance(
    data,
    dose,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  data$.wts <- as.numeric(stats::weights(fit))
  outcome_mod <- stats::lm(y_cont ~ dose, data = data, weights = .wts)

  # No warning first: the configuration is refused rather than fitted with a
  # diagnostic about the outcome model's terms.
  expect_no_warning(expect_error(
    ipw(fit, outcome_mod, .by = modifier),
    class = "balancing_ipw_unsupported_error"
  ))
})

# A missing value names no subgroup, so the units carrying one belong to none of
# the strata the effects would be reported within. Dropping them silently would
# report every subgroup over a different sample than the whole-sample rows use.

test_that(".by refuses a modifier with missing values", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )
  # The gap is carried on a column supplied through `.data`, since a modifier
  # the outcome model was fitted on could not have missing values and keep the
  # fit's row count.
  supplied <- data
  supplied$patchy <- supplied$modifier
  supplied$patchy[c(3L, 17L, 42L)] <- NA

  expect_no_warning(expect_error(
    ipw(fit, outcome_mod, .data = supplied, .by = patchy),
    class = "balancing_ipw_input_error"
  ))
})

# An effect within a subgroup contrasts the exposure levels inside it, so a
# subgroup holding one level alone has no contrast to report there. Refitting
# either model does not help, which is what makes this a refusal rather than a
# diagnostic: the data hold no comparison in that subgroup.

test_that(".by refuses a subgroup missing an exposure level", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )
  supplied <- data
  supplied$thin <- factor(
    ifelse(data$exposure == 1L & data$x2 > 0.8, "narrow", "wide"),
    levels = c("wide", "narrow")
  )

  # The subgroup is not empty and it is not tiny; it simply holds treated units
  # alone, which is the failure a count of its rows would not show.
  expect_gt(sum(supplied$thin == "narrow"), 10L)

  expect_no_warning(expect_error(
    ipw(fit, outcome_mod, .data = supplied, .by = thin),
    class = "balancing_ipw_input_error"
  ))
})

# The effects are reported within the levels of one variable, so a selection
# resolving to any other number of columns names no set of subgroups. Both ends
# are refused: nothing selected is not the same request as `.by = NULL`, which
# is the argument's default and asks for no subgroups at all, and two columns
# describe a crossing the caller has to build for themselves.

test_that(".by refuses a selection that is not exactly one column", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier + x1,
    data,
    w,
    stats::binomial()
  )

  expect_error(
    ipw(fit, outcome_mod, .by = tidyselect::starts_with("no_such_column")),
    class = "balancing_ipw_input_error"
  )
  expect_error(
    ipw(fit, outcome_mod, .by = c(modifier, x1)),
    class = "balancing_ipw_input_error"
  )
})

# The subgroups are the levels of the modifier, so the modifier has to name a
# fixed set of them. A numeric column names one subgroup per distinct value, and
# a logical one names subgroups whose labels would read as the conditions that
# produced them rather than as levels a caller declared, so both are refused
# with the recoding named rather than guessed at.

test_that(".by refuses a modifier that is not a factor or a character vector", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier + x1,
    data,
    w,
    stats::binomial()
  )
  supplied <- data
  supplied$flag <- data$x2 > 0

  expect_error(
    ipw(fit, outcome_mod, .by = x1),
    class = "balancing_ipw_input_error"
  )
  expect_error(
    ipw(fit, outcome_mod, .data = supplied, .by = flag),
    class = "balancing_ipw_input_error"
  )

  # The indicator the fit balances on holds the same subgroup structure as the
  # modifier, written as a number, so the refusal is about the column's type
  # rather than about what it holds and the recoding has to be asked for rather
  # than guessed at.
  expect_error(
    ipw(fit, outcome_mod, .data = supplied, .by = modifier_hi),
    class = "balancing_ipw_input_error"
  )
})

# The subgroup effects are g-computation on the outcome model as it was
# specified. A model with no term reading both the exposure and the modifier
# forces one and the same effect on every subgroup, so its subgroup rows agree
# by construction rather than by evidence. That is a modeling choice a caller
# may have made deliberately, so it is a warning and the result is still built.

test_that(".by warns when the outcome model has no exposure-by-modifier term", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure + modifier,
    data,
    w,
    stats::binomial()
  )

  expect_warning(
    result <- ipw(fit, outcome_mod, .by = modifier),
    class = "balancing_ipw_by_interaction_warning"
  )

  expect_s3_class(result, "ipw")
  expect_true("group" %in% names(result$estimates))
  expect_identical(nrow(result$estimates), 15L)
})

# The diagnostic reads the terms of the outcome model, and a model may carry the
# effect modification through a column derived from the modifier rather than
# through the modifier itself. `y ~ exposure * modifier_hi` reported by
# `.by = modifier` is that case: no term names `modifier`, so the diagnostic
# reports, while the subgroup effects it reports differ by as much as the
# fixture was drawn with. What the message says has to stay true here, which is
# why it reports which terms were read rather than announcing that the effect is
# the same in every subgroup.

test_that(".by reports the terms it read for a modifier carried by an indicator", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier_hi,
    data,
    w,
    stats::binomial()
  )

  expect_warning(
    result <- ipw(fit, outcome_mod, .data = data, .by = modifier),
    class = "balancing_ipw_by_interaction_warning"
  )

  estimates <- as.data.frame(result)
  expect_equal(
    by_estimate(estimates, "rd", "modifier = hi vs modifier = lo"),
    by_estimate(estimates, "rd", "modifier = hi") -
      by_estimate(estimates, "rd", "modifier = lo"),
    tolerance = 1e-10
  )
  expect_gt(
    abs(by_estimate(estimates, "rd", "modifier = hi vs modifier = lo")),
    0.1
  )
})

# ---- Recorded messages -----------------------------------------------------

# The wording and the condition class of everything `.by` refuses, and of the
# one thing it diagnoses. The warning is recorded off an invisible call rather
# than off the call itself: what the snapshot is for is the message, and letting
# the result print would carry nine rows of estimates into it alongside.

test_that("balancing_ipw_input_error: .by selects no column", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  expect_balancing_error(
    ipw(fit, outcome_mod, .by = tidyselect::starts_with("no_such_column"))
  )
})

test_that("balancing_ipw_input_error: .by selects two columns", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  # The model adjusts for `x1` so that the selection below reaches two columns
  # of the model frame. Named against a frame holding one of them, the selection
  # is tidyselect's own out-of-bounds failure rather than this refusal.
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier + x1,
    data,
    w,
    stats::binomial()
  )

  expect_balancing_error(ipw(fit, outcome_mod, .by = c(modifier, x1)))
})

test_that("balancing_ipw_input_error: .by names a modifier with missing values", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )
  supplied <- data
  supplied$patchy <- supplied$modifier
  supplied$patchy[c(3L, 17L, 42L)] <- NA

  expect_balancing_error(
    ipw(fit, outcome_mod, .data = supplied, .by = patchy)
  )
})

test_that("balancing_ipw_input_error: .by names a numeric modifier", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  expect_balancing_error(ipw(fit, outcome_mod, .data = data, .by = modifier_hi))
})

test_that("balancing_ipw_input_error: a stratum holds one exposure level", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )
  supplied <- data
  supplied$thin <- factor(
    ifelse(data$exposure == 1L & data$x2 > 0.8, "narrow", "wide"),
    levels = c("wide", "narrow")
  )

  expect_balancing_error(ipw(fit, outcome_mod, .data = supplied, .by = thin))
})

test_that("balancing_ipw_input_error: .data holds the fit's rows out of order", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )
  sorted <- data[order(data$modifier), , drop = FALSE]

  expect_balancing_error(ipw(fit, outcome_mod, .data = sorted, .by = modifier))
})

test_that("balancing_ipw_unsupported_error: .by on a continuous exposure", {
  data <- ipw_by_fixture()
  data$dose <- withr::with_seed(
    12,
    0.5 * data$x1 - 0.3 * data$x2 + stats::rnorm(nrow(data))
  )
  fit <- balance(
    data,
    dose,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  data$.wts <- as.numeric(stats::weights(fit))
  outcome_mod <- stats::lm(y_cont ~ dose, data = data, weights = .wts)

  expect_balancing_error(ipw(fit, outcome_mod, .by = modifier))
})

test_that("balancing_ipw_by_interaction_warning: no term reads both columns", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure + modifier,
    data,
    w,
    stats::binomial()
  )

  expect_balancing_warning(
    invisible(ipw(fit, outcome_mod, .by = modifier))
  )
})

# ---- The analytic contrast block ------------------------------------------

# A request appends two more deterministic blocks to the stack: each stratum's
# contrasts, written from that stratum's means, and each non-reference stratum's
# contrasts against the reference stratum's, written from the stratum contrast
# parameters. Both are constant across units, so both are candidates for an
# analytic bread row alongside the whole-sample contrasts, and a grouped fit is
# where the saving is largest.
#
# The reference system differences every one of those rows, which is what the
# package does today, and the reported system has to stay identical to it to the
# bit. The evaluation count beside it is red until they leave the differenced
# system.

test_that("a .by fit reports the fully differenced stacked system", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  frame <- stats::model.frame(outcome_mod)
  by <- ipw_resolve_by(
    rlang::quo(modifier),
    frame = frame,
    exposure = frame[["exposure"]],
    exposure_levels = fit@exposure_levels,
    exposure_name = "exposure",
    outcome_mod = outcome_mod
  )
  reference <- ipw_reference_stack(
    container = estimating_equations(fit),
    outcome_mod = outcome_mod,
    frame = frame,
    exposure_name = "exposure",
    levels = fit@exposure_levels,
    by = by
  )

  expect_ipw_matches_reference_stack(
    ipw(fit, outcome_mod, .by = modifier),
    reference,
    keys = c(
      "mu0",
      "mu1",
      "rd",
      "log(rr)",
      "log(or)",
      "mu0_modifier = lo",
      "mu1_modifier = lo",
      "mu0_modifier = hi",
      "mu1_modifier = hi",
      "rd_modifier = lo",
      "log(rr)_modifier = lo",
      "rd_modifier = hi",
      "log(rr)_modifier = hi",
      "rd_modifier = hi vs modifier = lo",
      "log(rr)_modifier = hi vs modifier = lo"
    )
  )
})

test_that("a .by fit differences no contrast row", {
  data <- ipw_by_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  frame <- stats::model.frame(outcome_mod)
  by <- ipw_resolve_by(
    rlang::quo(modifier),
    frame = frame,
    exposure = frame[["exposure"]],
    exposure_levels = fit@exposure_levels,
    exposure_name = "exposure",
    outcome_mod = outcome_mod
  )
  reference <- ipw_reference_stack(
    container = estimating_equations(fit),
    outcome_mod = outcome_mod,
    frame = frame,
    exposure_name = "exposure",
    levels = fit@exposure_levels,
    by = by
  )

  expect_stacked_evaluations(
    ipw(fit, outcome_mod, .by = modifier),
    2L * (reference$width - reference$deterministic) + 1L
  )
})

# The two cases above hold the ate surface of an entropy fit. What they cannot
# see is whether the analytic block still agrees once the rows around it change
# shape: a focal estimand standardizes every mean over the treated units alone,
# and non-uniform sampling weights enter the tilt and the reported weight scale
# both. Neither reaches the deterministic rows directly, since those rows read
# mean and contrast parameters and nothing else, and that is exactly why the
# case is worth pinning. An implementation that let the tilt leak into the
# rows it fills in analytically would still agree with the reference on the ate
# surface and disagree here.

test_that("a focal .by fit under sampling weights matches the differenced system", {
  data <- ipw_by_fixture()
  sampling <- withr::with_seed(2718, stats::runif(nrow(data), 0.4, 2.6))
  fit <- balance(
    data,
    exposure,
    c(x1, x2, modifier_hi),
    method = bw_ipt(),
    estimand = "att",
    sampling_weights = sampling
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_by_outcome(
    y ~ exposure * modifier,
    data,
    w,
    stats::binomial()
  )

  frame <- stats::model.frame(outcome_mod)
  by <- ipw_resolve_by(
    rlang::quo(modifier),
    frame = frame,
    exposure = frame[["exposure"]],
    exposure_levels = fit@exposure_levels,
    exposure_name = "exposure",
    outcome_mod = outcome_mod
  )
  reference <- ipw_reference_stack(
    container = estimating_equations(fit),
    outcome_mod = outcome_mod,
    frame = frame,
    exposure_name = "exposure",
    levels = fit@exposure_levels,
    by = by,
    sampling_weights = fit@sampling_weights,
    focal_level = fit@focal_level
  )

  expect_stacked_evaluations(
    expect_ipw_matches_reference_stack(
      ipw(fit, outcome_mod, .by = modifier),
      reference,
      keys = c(
        "mu0",
        "mu1",
        "rd",
        "log(rr)",
        "log(or)",
        "mu0_modifier = lo",
        "mu1_modifier = lo",
        "mu0_modifier = hi",
        "mu1_modifier = hi",
        "rd_modifier = lo",
        "log(rr)_modifier = lo",
        "rd_modifier = hi",
        "log(rr)_modifier = hi",
        "rd_modifier = hi vs modifier = lo",
        "log(rr)_modifier = hi vs modifier = lo"
      )
    ),
    2L * (reference$width - reference$deterministic) + 1L
  )
})
