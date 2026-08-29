# Print and summary snapshots for a fitted entropy result. The snapshots record
# on the first successful run once the implementation exists.

test_that("print() of a binary ate fit is stable", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_balancing_snapshot(print(fit))
})

test_that("summary() of a binary ate fit is stable", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_balancing_snapshot(summary(fit))
})

test_that("print() of a binary att fit renders the focal level", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att"
  )
  expect_balancing_snapshot(print(fit))
})

test_that("print() of a continuous ate fit is stable", {
  data <- sim_continuous()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_balancing_snapshot(print(fit))
})

test_that("the balance table carries one row per constraint term", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  balance_table <- fit@balance_table
  expect_true(all(
    c(
      "term",
      "kind",
      "statistic",
      "group",
      "unweighted",
      "weighted",
      "tolerance",
      "within_tolerance"
    ) %in%
      names(balance_table)
  ))
})

# A categorical exposure reports one effective sample size per level, so its
# print block exercises the multi-group path.
test_that("print() of a categorical ate fit lists every level", {
  data <- sim_categorical()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_balancing_snapshot(print(fit))
})

# A factor covariate contributes one constraint per level, and the full set is
# affine with the intercept, so the aliasing check drops the last level. The
# print block is the record that a fit still describes itself correctly once a
# constraint has gone: the term count is the surviving four, not the five the
# formula named. The block does not carry the balance table, so the surviving
# terms are asserted directly below it.
test_that("print() of a fit whose factor lost a level is stable", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2, x3),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_balancing_snapshot(print(fit))
  expect_identical(fit@balance_table$term, c("x1", "x2", "x3_a", "x3_b"))
})

# The quadratic-program family reports a solver backend and a minimum-weight
# floor, so its print and summary blocks differ from the estimating-equation
# family. The weight summary names the count of weights resting on the floor.
test_that("print() of an energy fit is stable", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  expect_balancing_snapshot(print(fit))
})

test_that("summary() of an energy fit reports the weight floor count", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  expect_balancing_snapshot(summary(fit))
})

test_that("summary() of a stable balancing fit reports the weight floor count", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  expect_balancing_snapshot(summary(fit))
})

test_that("summary() of a cfd fit reports the weight floor count", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
    estimand = "ate"
  )
  expect_balancing_snapshot(summary(fit))
  # The snapshot scrubs the count, which moves with the solver's floating-point
  # path, so the claim that this fixture puts weights on the floor at all is made
  # here instead. Recomputed with the rule the summary method applies, on the
  # same reported-weight scale. The count itself is not asserted: that is the
  # platform-volatile number.
  floor <- fit@method@min_weight
  at_floor <- sum(as.numeric(fit@weights) <= floor * (1 + 1e-6) + 1e-12)
  expect_gt(at_floor, 0)
})

# ---- The width of the largest-imbalance figure ------------------------------

# The headline figure and the balance-exceeded warning report the same quantity
# and render it the same way, to three significant digits. A fixed four decimal
# places spends its width on leading zeros, so it loses the digits that
# distinguish one small imbalance from another and collapses everything below
# half a ten-thousandth to the same "0.0000". The two fixtures below are the
# cases where the formats disagree: one imbalance a few thousandths wide, and
# one the fit drove to zero. Both are read off the fit rather than written in,
# because the low-order digits move with the platform's floating-point path.
#
# The `<1e-8` placeholder in helper-snapshot.R applies to recorded snapshots
# only. These tests capture the printed block directly, so they see the value
# the print method rendered.
test_that("print() renders the largest imbalance to three significant digits", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(),
    estimand = "ate"
  )
  largest <- max(abs(fit@balance_table$weighted))
  significant <- formatC(largest, format = "g", digits = 3)
  fixed <- formatC(largest, format = "f", digits = 4)
  # The fixture only says anything if the two formats disagree on it.
  expect_false(identical(significant, fixed))

  printed <- utils::capture.output(print(fit))
  line <- grep("Largest imbalance", printed, value = TRUE, fixed = TRUE)
  expect_length(line, 1L)
  # The trailing space and parenthesis pin the whole rendered figure. Without
  # them the fixed rendering is a prefix of the significant one, so a bare
  # substring test would accept either.
  expect_match(
    line,
    paste0("Largest imbalance: ", significant, " ("),
    fixed = TRUE
  )
  expect_no_match(
    line,
    paste0("Largest imbalance: ", fixed, " ("),
    fixed = TRUE
  )
})

test_that("print() keeps a largest imbalance below the fourth decimal legible", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  largest <- max(abs(fit@balance_table$weighted))
  expect_lt(largest, 1e-8)

  printed <- utils::capture.output(print(fit))
  line <- grep("Largest imbalance", printed, value = TRUE, fixed = TRUE)
  expect_length(line, 1L)
  expect_match(
    line,
    paste0(
      "Largest imbalance: ",
      formatC(largest, format = "g", digits = 3),
      " ("
    ),
    fixed = TRUE
  )
  expect_no_match(line, "Largest imbalance: 0.0000", fixed = TRUE)
})

# ---- An imbalance that was never measured ----------------------------------

# A balance statistic that is not a number states no distance, so the display
# reports that the imbalance could not be assessed rather than offering NaN as the
# largest one. That is the ruling the balance warning already follows, applied to
# the print and summary blocks. No fit reaches this once every exposure level is
# required to carry base-measure mass, so the table is edited to reach it.
undefined_balance_fit <- function(value) {
  data <- sim_binary(n = 150)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  table <- fit@balance_table
  table$weighted[[1]] <- value
  fit@balance_table <- table
  fit
}

test_that("print() reports an unmeasurable imbalance rather than NaN", {
  fit <- undefined_balance_fit(NaN)
  printed <- utils::capture.output(expect_no_warning(print(fit)))
  expect_false(any(grepl("NaN", printed, fixed = TRUE)))
  expect_true(any(grepl("could not be assessed", printed, fixed = TRUE)))
})

test_that("print() reports a missing imbalance rather than NA", {
  fit <- undefined_balance_fit(NA_real_)
  printed <- utils::capture.output(expect_no_warning(print(fit)))
  expect_false(any(grepl("NA", printed, fixed = TRUE)))
  expect_true(any(grepl("could not be assessed", printed, fixed = TRUE)))
})

test_that("summary() reports an unmeasurable imbalance rather than NaN", {
  fit <- undefined_balance_fit(NaN)
  printed <- utils::capture.output(expect_no_warning(summary(fit)))
  expect_true(any(grepl("could not be assessed", printed, fixed = TRUE)))
  # The balance table itself still shows the term's value; only the headline
  # figure, which claims to measure how far the fit missed, is withheld.
  expect_false(any(grepl("Largest imbalance: NaN", printed, fixed = TRUE)))
})

test_that("summary() returns the fit invisibly", {
  data <- sim_binary(200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  returned <- withr::with_output_sink(tempfile(), summary(fit))
  expect_true(S7::S7_inherits(returned, balancing))
  expect_identical(returned@estimand, "ate")
})
