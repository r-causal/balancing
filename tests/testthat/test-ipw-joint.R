# Specs for reporting a joint exposure through balancing's `ipw()` method: two
# treatments intervened on at once, crossed into one categorical exposure by
# `causalgenerics::joint_exposure()`.
#
# The crossing is a factor subclass, so the categorical path already fits it and
# already reports what it reports for any four-level exposure: each cell against
# the reference cell, under labels such as `"a = 1, e = 0 vs a = 0, e = 0"`.
# Those rows are arithmetically right and they answer a question nobody asked. A
# declared crossing is reported in the two treatments instead:
#
#   * the counterfactual mean of each cell as a row of its own, under the effect
#     label `"mean"`, keyed by the cell;
#   * the simple effects, each treatment's effect within a fixed level of the
#     other, keyed by the contrast naming the treatment and the group naming the
#     level it is held at, including the comparisons vs-reference reporting
#     cannot express at all;
#   * the interaction, the difference between two of the first treatment's
#     simple effects, on each collapsible scale.
#
# The odds ratio is off every contrast row here for the reason it is off a
# stratum row: it is noncollapsible, so neither a simple effect reported beside
# one nor a difference of two of them says what it appears to.
#
# What these pin is that grammar, the arithmetic behind it, the route by which
# the declaration reaches `ipw()`, and that the flat cell-against-cell labels
# become unreachable rather than staying available by accident.

# ---- Fixtures --------------------------------------------------------------

# The two treatments of `ipw_joint_fixture()` in helper-dgp.R crossed under one
# name, assembled from the parts
# rather than declared through `causalgenerics::joint_exposure()`, which now
# refuses two components sharing a name.
#
# That refusal covers the public constructor alone. The internal one the class's
# own subsetting and coercion methods call is unvalidated by design, so a
# crossing carrying one name twice stays constructible by any package holding
# the pieces, and balancing reads its rows back by a key written from the
# component names. Building the object here is what keeps balancing's own
# refusal under test against a crossing that never passed the upstream check.
#
# The layout is the one `joint_exposure()` writes: an integer core indexing the
# cells, the cells as `levels`, each component's levels as `components` named
# for the treatment it belongs to, and a class vector placing `"factor"` ahead
# of `"vctrs_vctr"` so the formula machinery treats the result natively. The
# cells vary the first component fastest, which is what puts the reference cell
# first. What the crossing carries is asserted where it is used, so the assembly
# cannot quietly build something else.
ipw_joint_collision_fixture <- function() {
  data <- ipw_joint_fixture()
  first <- levels(data$a)
  second <- levels(data$e)
  cells <- paste0(
    "a = ",
    rep(first, times = length(second)),
    ", a = ",
    rep(second, each = length(first))
  )
  codes <- as.integer(data$a) + (as.integer(data$e) - 1L) * length(first)
  data$joint <- structure(
    codes,
    levels = cells,
    components = stats::setNames(list(first, second), c("a", "a")),
    class = c("joint_exposure", "factor", "vctrs_vctr", "integer")
  )
  data
}

# The cells, in the order the crossing declares them: the first component varies
# fastest, so the reference cell crosses the two reference levels and comes
# first.
joint_cells <- c(
  "a = 0, e = 0",
  "a = 1, e = 0",
  "a = 0, e = 1",
  "a = 1, e = 1"
)

# The labels the categorical path reports when nothing is declared, which is
# what a declared crossing must never produce.
joint_flat_contrasts <- paste(joint_cells[-1], "vs", joint_cells[[1]])

# The weights every joint fixture is fitted with, solved on a plain factor
# carrying exactly the same cells under exactly the same labels.
#
# Two things follow from fitting the flattened copy. The weights are provably
# the plain categorical path's weights, so any difference between the two
# reported surfaces comes from the declaration and from nothing else, which is
# what makes the agreement comparisons below a test of the declaration. And the
# specs for the reported surface do not rest on the separate question of whether
# `balance()` itself takes a declared column, which is pinned on its own below.
fit_joint_weights <- function(data, estimand = "ate", focal_level = NULL) {
  flattened <- data
  flattened$joint <- factor(
    as.character(data$joint),
    levels = levels(data$joint)
  )
  balance(
    flattened,
    joint,
    c(x1),
    method = bw_ipt(),
    estimand = estimand,
    focal_level = focal_level
  )
}

# A weighted outcome model of the shape `ipw()` expects, with the weights riding
# along as a column so the model frame resolves them.
fit_joint_outcome <- function(formula, data, wts, family) {
  data[[".wts"]] <- wts
  suppressWarnings(
    stats::glm(formula, data = data, family = family, weights = .wts)
  )
}

# Every call that reaches a declared crossing in this file goes through this
# gate, whether it is expected to return a result or to refuse.
#
# The counterfactual designs are built by fixing the exposure to one cell at a
# time, and a joint exposure that loses cells gives up its declaration and says
# so, so designs built by writing into the declared column warn once per cell.
# The gate turns that into a failure rather than into noise the suite would
# carry, and it is the assertion that says the designs are built from a plain
# factor over the same cells.
#
# The handler collects rather than stopping at the first, because
# `expect_no_warning()` muffles only the one it reports and lets the rest
# through to the reporter, which is how four warnings become one failure and
# three the suite counts. Collecting them reports the whole set at once and
# leaves the value there to read either way.
expect_joint_quiet <- function(expr) {
  observed <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(cnd) {
      observed <<- c(observed, conditionMessage(cnd))
      rlang::cnd_muffle(cnd)
    }
  )
  testthat::expect_identical(observed, character())
  value
}

joint_ipw <- function(...) {
  expect_joint_quiet(ipw(...))
}

# The four counterfactual means: predict the outcome model with the exposure
# fixed to each cell in turn and average over the target population. `tilt` is
# that population's weight, uniform for a pooled estimand and the focal cell's
# indicator for a focal one.
#
# The counterfactual column is written as a plain factor rather than assigned
# into the declared one, for the reason the gate above exists: setting every row
# to one cell leaves the other three unpopulated, and the declaration would be
# given up before the prediction was taken.
joint_cell_means <- function(
  outcome_mod,
  data,
  tilt = NULL,
  exposure_name = "joint"
) {
  weight <- tilt %||% rep(1, nrow(data))
  vapply(
    joint_cells,
    function(cell) {
      counterfactual <- data
      counterfactual[[exposure_name]] <- factor(cell, levels = joint_cells)
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

# One row's estimate, read out of a coerced result by the three columns that
# name it. The lookup says how many rows it matched before it returns: a key
# matching nothing gives an empty vector, and a comparison between two of those
# holds, so a reader that did not report the count could report agreement
# between rows that are not there.
joint_estimate <- function(estimates, effect, contrast, group) {
  rows <- estimates$term == effect &
    estimates$contrast == contrast &
    estimates$group == group
  value <- estimates$estimate[rows]
  testthat::expect_identical(
    length(value),
    1L,
    label = paste(effect, contrast, group)
  )
  value
}

# ---- The declaration -------------------------------------------------------

# The fixture rests on the crossing, so what it declares is pinned against
# causalgenerics before anything is read through it.

test_that("the fixture declares the crossing the joint surface is written in", {
  data <- ipw_joint_fixture()

  expect_true(causalgenerics::is_joint_exposure(data$joint))
  expect_identical(levels(data$joint), joint_cells)
  expect_identical(
    causalgenerics::joint_components(data$joint),
    list(a = c("0", "1"), e = c("0", "1"))
  )
  expect_identical(
    causalgenerics::joint_reference(data$joint),
    joint_cells[[1]]
  )

  # Every cell is populated, which is what a crossing needs to be identified and
  # what the four mean rows each stand for.
  expect_true(all(table(data$joint) > 0L))
})

# A declared column is a factor over the cells, so weighting it is weighting
# them: the crossing changes which effects are reported and not which population
# is balanced. `balance()` therefore has to take one and return the fit the
# flattened column returns. Today it neither weights it nor turns it away: the
# exposure's finiteness check reaches a vctrs method the class does not carry,
# and the caller gets an error about `is.infinite()` that names nothing they
# did.

test_that("balance() weights a declared crossing as it weights its cells", {
  data <- ipw_joint_fixture()

  fit <- expect_no_error(
    balance(data, joint, c(x1), method = bw_ipt(), estimand = "ate")
  )
  reference <- fit_joint_weights(data)

  expect_identical(fit@exposure_type, "categorical")
  expect_identical(fit@exposure_levels, joint_cells)
  expect_equal(
    as.numeric(stats::weights(fit)),
    as.numeric(stats::weights(reference)),
    tolerance = 1e-10
  )
})

# ---- The reported surface --------------------------------------------------

test_that("a declared crossing reports cell means, simple effects, and their interaction", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )

  # The stored table, for the reason every other column contract in this suite
  # reads it: the coerced frame is causalgenerics' contract, and this is the one
  # balancing fills in.
  estimates <- joint_ipw(fit, outcome_mod)$estimates

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

  measures <- c("rd", "log(rr)")
  expect_identical(
    estimates$effect,
    c(rep("mean", 4L), rep(measures, times = 5L))
  )
  expect_identical(
    estimates$contrast,
    c(
      joint_cells,
      rep("a: 1 vs 0", 4L),
      rep("e: 1 vs 0", 4L),
      rep("a: 1 vs 0", 2L)
    )
  )
  expect_identical(
    estimates$group,
    c(
      rep("overall", 4L),
      rep(c("e = 0", "e = 1"), each = 2L),
      rep(c("a = 0", "a = 1"), each = 2L),
      rep("e = 1 vs e = 0", 2L)
    )
  )
  expect_identical(nrow(estimates), 14L)
  expect_true(all(is.finite(estimates$estimate)))
})

# The whole point of the declaration is that the cells stop being the vocabulary
# the effects are reported in. A surface still carrying the cell-against-cell
# strings anywhere a caller reads is the surface the declaration was supposed to
# replace, so every place a row can be named is checked rather than the stored
# contrast column alone.

test_that("a declared crossing never names a cell against the reference cell", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )

  result <- joint_ipw(fit, outcome_mod)
  written <- c(
    result$estimates$contrast,
    result$estimates$group,
    names(stats::coef(result)),
    rownames(stats::vcov(result)),
    colnames(stats::vcov(result)),
    rownames(stats::confint(result)),
    names(result$fit$theta)
  )

  for (flat in joint_flat_contrasts) {
    expect_false(
      any(grepl(flat, written, fixed = TRUE)),
      label = paste0("the flat contrast ", flat)
    )
  }

  # The printed table is where a caller meets the labels first, so it is read
  # as well rather than trusted to follow from the accessors.
  output <- utils::capture.output(print(result))
  for (flat in joint_flat_contrasts) {
    expect_false(
      any(grepl(flat, output, fixed = TRUE)),
      label = paste0("the printed flat contrast ", flat)
    )
  }
})

test_that("a declared crossing reports no odds ratio", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )

  estimates <- joint_ipw(fit, outcome_mod)$estimates

  expect_false("log(or)" %in% estimates$effect)
  expect_identical(
    sort(unique(estimates$effect)),
    c("log(rr)", "mean", "rd")
  )
})

test_that("a gaussian outcome reports one difference per simple effect", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y_cont ~ joint + x1,
    data,
    w,
    stats::gaussian()
  )

  estimates <- joint_ipw(fit, outcome_mod)$estimates

  expect_identical(estimates$effect, c(rep("mean", 4L), rep("diff", 5L)))
  expect_identical(
    estimates$contrast,
    c(
      joint_cells,
      "a: 1 vs 0",
      "a: 1 vs 0",
      "e: 1 vs 0",
      "e: 1 vs 0",
      "a: 1 vs 0"
    )
  )
  expect_identical(
    estimates$group,
    c(
      rep("overall", 4L),
      "e = 0",
      "e = 1",
      "a = 0",
      "a = 1",
      "e = 1 vs e = 0"
    )
  )
  expect_identical(nrow(estimates), 9L)
})

# ---- The arithmetic --------------------------------------------------------

test_that("the cell mean rows are the weighted g-computation means", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  # An adjusted model, so the cell means are a real g-computation average rather
  # than the weighted cell outcomes a saturated model collapses to. Under
  # `y ~ joint` alone the counterfactual predictions are constant within a cell,
  # and an implementation averaging the wrong population would match anyway.
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )
  mu <- joint_cell_means(outcome_mod, data)

  estimates <- as.data.frame(joint_ipw(fit, outcome_mod))

  for (index in seq_along(joint_cells)) {
    expect_equal(
      joint_estimate(estimates, "mean", joint_cells[[index]], "overall"),
      unname(mu[[index]]),
      tolerance = 1e-8
    )
  }
})

test_that("each simple effect contrasts two cell means", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )
  mu <- joint_cell_means(outcome_mod, data)

  estimates <- as.data.frame(joint_ipw(fit, outcome_mod))

  # `mu` is indexed in cell order, so the first is the reference cell, the
  # second moves the first treatment, the third moves the second, and the fourth
  # moves both. A simple effect fixes one treatment and contrasts the other.
  simple <- list(
    list(contrast = "a: 1 vs 0", group = "e = 0", hi = 2L, lo = 1L),
    list(contrast = "a: 1 vs 0", group = "e = 1", hi = 4L, lo = 3L),
    list(contrast = "e: 1 vs 0", group = "a = 0", hi = 3L, lo = 1L),
    list(contrast = "e: 1 vs 0", group = "a = 1", hi = 4L, lo = 2L)
  )

  for (effect in simple) {
    expect_equal(
      joint_estimate(estimates, "rd", effect$contrast, effect$group),
      unname(mu[[effect$hi]] - mu[[effect$lo]]),
      tolerance = 1e-8
    )
    expect_equal(
      joint_estimate(estimates, "log(rr)", effect$contrast, effect$group),
      unname(log(mu[[effect$hi]]) - log(mu[[effect$lo]])),
      tolerance = 1e-8
    )
  }
})

# The interaction is one number reported once, under the first treatment's
# framing. Under the second treatment's framing it is the same number, which is
# an identity of the four means rather than a second row: reporting both would
# put one quantity in the table twice under two names.

test_that("the interaction is the difference of two simple effects either way", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )
  mu <- joint_cell_means(outcome_mod, data)

  result <- joint_ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)
  reported <- joint_estimate(
    estimates,
    "rd",
    "a: 1 vs 0",
    "e = 1 vs e = 0"
  )

  first_framing <- (mu[[4]] - mu[[3]]) - (mu[[2]] - mu[[1]])
  second_framing <- (mu[[4]] - mu[[2]]) - (mu[[3]] - mu[[1]])
  expect_equal(unname(first_framing), unname(second_framing), tolerance = 1e-12)
  expect_equal(reported, unname(first_framing), tolerance = 1e-8)

  # It is also the difference of the two reported simple-effect rows, which is
  # what says the row is a parameter of the same system rather than a number
  # computed beside it.
  expect_equal(
    reported,
    joint_estimate(estimates, "rd", "a: 1 vs 0", "e = 1") -
      joint_estimate(estimates, "rd", "a: 1 vs 0", "e = 0"),
    tolerance = 1e-10
  )

  # The fixture carries a real interaction, so the row is not a rounding
  # artifact of two simple effects that agree.
  expect_gt(abs(reported), 0.1)

  # Once, under one framing. The second treatment's framing names no row.
  expect_identical(sum(estimates$group == "e = 1 vs e = 0"), 2L)
  expect_false(any(grepl("a = 1 vs a = 0", estimates$group, fixed = TRUE)))
})

# The declaration changes which rows are reported and not what any of them is
# worth. Two of the vs-reference contrasts are simple effects under another
# name: the cell that moves the first treatment alone against the reference is
# that treatment's effect at the second's reference level, and likewise the
# other way. Those rows have to agree on the estimate and on the standard error,
# since both surfaces read the same stacked system.

test_that("the joint surface agrees with the categorical path on the rows both report", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  flattened <- data
  flattened$joint <- factor(
    as.character(data$joint),
    levels = levels(data$joint)
  )
  declared_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )
  flat_mod <- fit_joint_outcome(
    y ~ joint + x1,
    flattened,
    w,
    stats::binomial()
  )

  declared <- as.data.frame(joint_ipw(fit, declared_mod))
  plain <- as.data.frame(ipw(fit, flat_mod))

  # The undeclared surface is the one the categorical path has always reported,
  # pinned here so the comparison is against a known table rather than against
  # whatever the other call happened to produce.
  expect_identical(
    plain$contrast,
    c(joint_cells, rep(joint_flat_contrasts, each = 3L))
  )
  expect_false("group" %in% names(plain))

  shared <- list(
    list(
      flat = joint_flat_contrasts[[1]],
      contrast = "a: 1 vs 0",
      group = "e = 0"
    ),
    list(
      flat = joint_flat_contrasts[[2]],
      contrast = "e: 1 vs 0",
      group = "a = 0"
    )
  )
  for (row in shared) {
    for (effect in c("rd", "log(rr)")) {
      index <- which(plain$term == effect & plain$contrast == row$flat)
      expect_identical(length(index), 1L)
      expect_equal(
        joint_estimate(declared, effect, row$contrast, row$group),
        plain$estimate[index],
        tolerance = 1e-10
      )
      expect_equal(
        declared$std.error[
          declared$term == effect &
            declared$contrast == row$contrast &
            declared$group == row$group
        ],
        plain$std.error[index],
        tolerance = 1e-10
      )
    }
  }
})

test_that("a declared crossing reports a usable standard error for every row", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )

  result <- joint_ipw(fit, outcome_mod)
  estimates <- result$estimates

  expect_identical(nrow(estimates), 14L)
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

# Every row of the surface is a parameter of one stacked system, so the rows
# covary. Simple effects computed on their own subsets and assembled into a
# block-diagonal matrix would report exact zeros where the coupling belongs,
# and the interaction's standard error would come out as the sum of two
# variances rather than as the variance of their difference.

test_that("a declared crossing's covariance couples the rows it reports", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )

  covariance <- stats::vcov(joint_ipw(fit, outcome_mod))
  at_zero <- "rd a: 1 vs 0 e = 0"
  at_one <- "rd a: 1 vs 0 e = 1"
  interaction <- "rd a: 1 vs 0 e = 1 vs e = 0"

  expect_true(all(c(at_zero, at_one, interaction) %in% rownames(covariance)))
  expect_equal(covariance, t(covariance), tolerance = 1e-12)

  coupling <- covariance[at_zero, at_one]
  expect_true(is.finite(coupling))
  expect_gt(abs(coupling), 1e-8)

  # The entries combined here are finite-difference results, so the identity is
  # exact algebra read off independently rounded numbers, and the last bits of
  # the cancellation follow the platform's BLAS and compiler. The tolerance is
  # the one every assertion of this shape carries: loose enough for that
  # spread, three orders below the gap to the stitched sum the next assertion
  # refuses.
  expect_equal(
    covariance[interaction, interaction],
    covariance[at_zero, at_zero] +
      covariance[at_one, at_one] -
      2 * coupling,
    tolerance = 1e-6
  )
  expect_false(isTRUE(all.equal(
    covariance[interaction, interaction],
    covariance[at_zero, at_zero] + covariance[at_one, at_one]
  )))
})

# A declared crossing takes the contrast block over rather than sitting beside
# it, so the assembly of the stacked matrix meets a block of a shape no other
# exposure produces. What the expectation compares is the assembled matrix
# against the `rbind()` of the same blocks, at every evaluation the finite
# difference asks the closure for.

test_that("a declared crossing stacks its psi blocks as rbind would", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )

  expect_stacked_psi_matches_rbind(joint_ipw(fit, outcome_mod))
})

# ---- Labels ----------------------------------------------------------------

# The measure repeats across the simple effects and the cells, so a row is named
# by all three of its identity columns. balancing builds those labels for the
# covariance it attaches to the estimates table and causalgenerics builds them
# again for the accessors and the printed table, and the two have to agree or a
# caller reading a covariance out by the name `coef()` gave gets a different
# entry from the one `print()` showed.

test_that("a declared crossing labels its rows by measure, contrast, and cell", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )

  result <- joint_ipw(fit, outcome_mod)
  labels <- c(
    paste("mean", joint_cells, "overall"),
    "rd a: 1 vs 0 e = 0",
    "log(rr) a: 1 vs 0 e = 0",
    "rd a: 1 vs 0 e = 1",
    "log(rr) a: 1 vs 0 e = 1",
    "rd e: 1 vs 0 a = 0",
    "log(rr) e: 1 vs 0 a = 0",
    "rd e: 1 vs 0 a = 1",
    "log(rr) e: 1 vs 0 a = 1",
    "rd a: 1 vs 0 e = 1 vs e = 0",
    "log(rr) a: 1 vs 0 e = 1 vs e = 0"
  )

  expect_identical(length(labels), 14L)
  expect_identical(anyDuplicated(labels), 0L)
  expect_identical(
    dimnames(attr(result$estimates, "ipw_vcov", exact = TRUE)),
    list(labels, labels)
  )
  expect_identical(names(stats::coef(result)), labels)
  expect_identical(dimnames(stats::vcov(result)), list(labels, labels))
  expect_identical(rownames(stats::confint(result)), labels)
})

# ---- How the declaration reaches ipw() -------------------------------------

# The fit records the levels it weighted as plain strings and keeps no memory of
# the crossing, so the declaration cannot come from there. It comes from the
# exposure column of the frame the discrete path already resolves: the outcome
# model's own frame, or `.data` when the caller supplied one. Both arms are
# pinned, since the frame is resolved once and either arm can fill it.

test_that("ipw() reads the declaration from the outcome model's own frame", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )

  # The fit was solved on the flattened copy, so nothing it stores says the
  # exposure was a crossing. That is what makes the frame the only source.
  expect_identical(fit@exposure_levels, joint_cells)
  expect_false(causalgenerics::is_joint_exposure(fit@exposure_levels))

  frame <- stats::model.frame(outcome_mod)
  expect_true(causalgenerics::is_joint_exposure(frame$joint))

  estimates <- joint_ipw(fit, outcome_mod)$estimates
  expect_identical(estimates$effect[1:4], rep("mean", 4L))
  expect_identical(estimates$contrast[1:4], joint_cells)
})

test_that("ipw() reads the declaration from .data when one is supplied", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  flattened <- data
  flattened$joint <- factor(
    as.character(data$joint),
    levels = levels(data$joint)
  )
  # The outcome model is fitted on the flattened column, so its own frame
  # declares nothing and the supplied frame is the only place the crossing is
  # written.
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    flattened,
    w,
    stats::binomial()
  )
  expect_false(
    causalgenerics::is_joint_exposure(stats::model.frame(outcome_mod)$joint)
  )

  estimates <- joint_ipw(fit, outcome_mod, .data = data)$estimates

  expect_identical(nrow(estimates), 14L)
  expect_identical(estimates$contrast[1:4], joint_cells)
  expect_identical(estimates$group[[5]], "e = 0")
})

# The declaration is what changes the surface, so its absence has to leave the
# categorical path exactly where it was: one mean per cell, then each cell
# against the reference cell over all three measures including the odds ratio,
# and no subgroup column. Both surfaces lead with the same cell means; what the
# declaration changes is the contrasts written from them, so the labels are
# where the two part company.

test_that("an undeclared crossing still reports each cell against the reference", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  flattened <- data
  flattened$joint <- factor(
    as.character(data$joint),
    levels = levels(data$joint)
  )
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    flattened,
    w,
    stats::binomial()
  )

  estimates <- ipw(fit, outcome_mod)$estimates

  expect_identical(nrow(estimates), 13L)
  expect_identical(
    estimates$effect,
    c(
      rep("mean", length(joint_cells)),
      rep(c("rd", "log(rr)", "log(or)"), times = 3L)
    )
  )
  expect_identical(
    estimates$contrast,
    c(joint_cells, rep(joint_flat_contrasts, each = 3L))
  )
  expect_false("group" %in% names(estimates))
})

# ---- Refusals --------------------------------------------------------------

# Every row of the surface is keyed by the treatment it contrasts and the level
# the other treatment is held at, and both of those are written from the
# component names. Two treatments under one name therefore write one key over
# two different effects, and since the rows are read back by name, the reported
# table would carry the first of each colliding pair twice in place of the two
# effects the crossing holds. The declaration is turned away instead.

test_that("a declared crossing refuses two treatments under one name", {
  data <- ipw_joint_collision_fixture()

  # The crossing itself is well formed, which is why the refusal is balancing's
  # to carry: the cells are distinct and every one of them is populated, and it
  # is only reporting in the two treatments that the shared name is a defect.
  expect_true(causalgenerics::is_joint_exposure(data$joint))
  expect_identical(
    names(causalgenerics::joint_components(data$joint)),
    c("a", "a")
  )
  expect_identical(anyDuplicated(levels(data$joint)), 0L)
  expect_true(all(table(data$joint) > 0L))

  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )

  expect_joint_quiet(expect_error(
    ipw(fit, outcome_mod),
    class = "balancing_ipw_input_error"
  ))
})

# A joint exposure already reports an interaction between two treatments.
# Reporting it again within the levels of a modifier is a three-way question,
# and the surface answers a two-way one, so the combination is turned away
# rather than answered in some reading of it the labels would not distinguish.

test_that("a declared crossing refuses .by", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint * modifier + x1,
    data,
    w,
    stats::binomial()
  )

  expect_joint_quiet(expect_error(
    ipw(fit, outcome_mod, .by = modifier),
    class = "balancing_ipw_unsupported_error"
  ))
})

# Every cell mean on the joint surface standardizes to one population. A tilted
# estimand standardizes each of them to the focal cell's population, over which
# the simple effects and the interaction are not defined, so the surface is
# reported for the pooled estimand alone.

test_that("a declared crossing is reported for the ate estimand alone", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(
    data,
    estimand = "att",
    focal_level = "a = 1, e = 1"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )

  expect_joint_quiet(expect_error(
    ipw(fit, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  ))
})

# ---- Recorded messages -----------------------------------------------------

# The wording and the condition class of the three configurations a declared
# crossing turns away. Each of them refuses before any counterfactual design is
# built, so no snapshot here carries a declaration-loss warning alongside the
# refusal, and a change that moved any of the checks after the designs would show
# up here as well as in the gates above.

test_that("balancing_ipw_input_error: two treatments under one name", {
  data <- ipw_joint_collision_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )

  expect_balancing_error(ipw(fit, outcome_mod))
})

test_that("balancing_ipw_unsupported_error: .by on a declared crossing", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint * modifier + x1,
    data,
    w,
    stats::binomial()
  )

  expect_balancing_error(ipw(fit, outcome_mod, .by = modifier))
})

test_that("balancing_ipw_unsupported_error: a focal estimand on a crossing", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(
    data,
    estimand = "att",
    focal_level = "a = 1, e = 1"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(
    y ~ joint + x1,
    data,
    w,
    stats::binomial()
  )

  expect_balancing_error(ipw(fit, outcome_mod))
})

# ---- Components that cannot be crossed -------------------------------------

# A continuous component cannot be declared in the first place. Crossing a
# treatment with a continuous variable puts one unit in almost every cell and
# none in the rest, and a crossing with an empty cell is not identified, so
# causalgenerics refuses it at construction. Nothing of the kind can reach
# `ipw()`, which is why balancing carries no refusal of its own for it. The
# control is here rather than upstream because it is the reason for that
# absence.

test_that("a continuous component cannot be declared at all", {
  data <- ipw_joint_fixture()

  expect_error(
    causalgenerics::joint_exposure(a = data$a, dose = data$x1),
    class = "causalgenerics_invalid_joint_exposure"
  )

  # The whole crossing has to be populated, so a component with more levels than
  # the data support is refused on the same terms rather than on its type.
  coarse <- cut(data$x1, breaks = c(-Inf, 0, Inf), labels = c("lo", "hi"))
  expect_true(causalgenerics::is_joint_exposure(
    causalgenerics::joint_exposure(a = data$a, dose = coarse)
  ))
})

# ---- The analytic contrast block ------------------------------------------

# A declared crossing replaces the vs-reference contrast block with the simple
# effects and the interaction, and those rows are deterministic on the same
# terms: a simple effect is a contrast of two mean parameters and an interaction
# row is the difference of two simple-effect parameters, so both are constant
# across units and both have bread rows that are known without differencing
# anything. This is the widest contrast block any surface reports, so it is the
# one an analytic bread block saves the most on.
#
# The reference system differences every one of those rows, which is what the
# package does today, and the reported system has to stay identical to it to the
# bit. The evaluation count beside it is red until they leave the differenced
# system.

test_that("a declared crossing reports the fully differenced stacked system", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(y ~ joint + x1, data, w, stats::binomial())

  frame <- stats::model.frame(outcome_mod)
  joint <- ipw_joint_plan(
    frame[["joint"]],
    fit@exposure_levels,
    is_gaussian_outcome(outcome_mod)
  )
  frame[["joint"]] <- ipw_joint_bare(frame[["joint"]])
  reference <- ipw_reference_stack(
    container = estimating_equations(fit),
    outcome_mod = outcome_mod,
    frame = frame,
    exposure_name = "joint",
    levels = fit@exposure_levels,
    categorical = TRUE,
    joint = joint
  )

  expect_ipw_matches_reference_stack(
    expect_joint_quiet(ipw(fit, outcome_mod)),
    reference,
    keys = c(
      paste0("mu_", joint_cells),
      "rd_a: 1 vs 0 e = 0",
      "log(rr)_a: 1 vs 0 e = 0",
      "rd_a: 1 vs 0 e = 1",
      "log(rr)_a: 1 vs 0 e = 1",
      "rd_e: 1 vs 0 a = 0",
      "log(rr)_e: 1 vs 0 a = 0",
      "rd_e: 1 vs 0 a = 1",
      "log(rr)_e: 1 vs 0 a = 1",
      "rd_a: 1 vs 0 e = 1 vs e = 0",
      "log(rr)_a: 1 vs 0 e = 1 vs e = 0"
    )
  )
})

test_that("a declared crossing differences no contrast row", {
  data <- ipw_joint_fixture()
  fit <- fit_joint_weights(data)
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_joint_outcome(y ~ joint + x1, data, w, stats::binomial())

  frame <- stats::model.frame(outcome_mod)
  joint <- ipw_joint_plan(
    frame[["joint"]],
    fit@exposure_levels,
    is_gaussian_outcome(outcome_mod)
  )
  frame[["joint"]] <- ipw_joint_bare(frame[["joint"]])
  reference <- ipw_reference_stack(
    container = estimating_equations(fit),
    outcome_mod = outcome_mod,
    frame = frame,
    exposure_name = "joint",
    levels = fit@exposure_levels,
    categorical = TRUE,
    joint = joint
  )

  expect_stacked_evaluations(
    expect_joint_quiet(ipw(fit, outcome_mod)),
    2L * (reference$width - reference$deterministic) + 1L
  )
})
