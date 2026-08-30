# Specs for the assembly of the stacked estimating function that `ipw()`'s
# variance is read from. The stack is built one block at a time, and the blocks
# themselves are pinned elsewhere: what these cover is only the step that puts
# them together into the S-by-n matrix the sandwich differentiates.
#
# The step is `stack_psi_blocks()`, which fills a preallocated matrix rather
# than growing one through `rbind()`. Its whole contract is that the two agree
# to the bit, values and dimnames alike, so both halves of that contract are
# stated here: the helper against `rbind()` on blocks written out by hand,
# including the shapes a stack actually presents it with, and the helper against
# `rbind()` on the blocks a real fit builds.
#
# Beside it sits `sum_psi_blocks()`, which the bread reaches instead, since the
# bread reads only row sums and has no use for the matrix. Its contract is the
# assembly's read through `rowSums()`, and it is stated the same two ways.

# ---- Fixtures --------------------------------------------------------------

# A binary-exposure fixture with an outcome that depends on the exposure and on
# one covariate, so an adjusted outcome model has something to adjust for and
# the mean rows standardize over predictions that vary by unit. The exposure and
# the covariates come from the shared `sim_binary()` process; the outcome is
# drawn here under its own seed.
ipw_deli_fixture <- function(n = 200) {
  data <- sim_binary(n)
  withr::with_seed(717, {
    data$y <- stats::rbinom(
      n,
      1L,
      stats::plogis(-0.3 + 0.5 * data$exposure + 0.4 * data$x1)
    )
  })
  data
}

# The same shape at three exposure levels, which is what puts more than one
# contrast row in the stack.
ipw_deli_categorical_fixture <- function(n = 200) {
  data <- sim_categorical(n)
  withr::with_seed(818, {
    linear_predictor <- -0.4 +
      0.3 * (data$exposure == "b") +
      0.6 * (data$exposure == "c") +
      0.5 * data$x1
    data$y <- stats::rbinom(n, 1L, stats::plogis(linear_predictor))
  })
  data
}

# A weighted outcome model of the shape `ipw()` expects, with the weights riding
# along as a column so the model frame resolves them.
fit_deli_outcome <- function(formula, data, wts, family) {
  data[[".wts"]] <- wts
  suppressWarnings(
    stats::glm(formula, data = data, family = family, weights = .wts)
  )
}

# ---- The assembly itself ---------------------------------------------------

# The blocks a stack presents are not all matrices of the same shape. A route
# that reports no contrasts of its own passes nothing for that block, and a
# request that adds a block the fit turned out to have no rows for passes an
# empty one, so both have to stack to what `rbind()` stacks them to: dropped
# entirely, contributing no row and no row name.
#
# Every block a stack carries holds doubles, so the blocks written out here hold
# doubles too. Nothing states what the assembly does with an integer block,
# because nothing builds one.

test_that("stack_psi_blocks() stacks present, absent, and empty blocks as rbind does", {
  n <- 4L
  weight_block <- matrix(
    seq_len(2L * n) / 10,
    nrow = 2L,
    dimnames = list(c("theta_w1", "theta_w2"), NULL)
  )
  score_block <- matrix(seq_len(3L * n) * 2.5, nrow = 3L)
  empty_block <- matrix(numeric(), nrow = 0L, ncol = n)
  contrast_block <- matrix(
    -1,
    nrow = 1L,
    ncol = n,
    dimnames = list("rd", NULL)
  )

  blocks <- list(
    weight_block,
    score_block,
    NULL,
    empty_block,
    contrast_block
  )

  expect_identical(
    stack_psi_blocks(blocks, n),
    rbind(weight_block, score_block, empty_block, contrast_block)
  )
})

# Row names are the half of the contract a preallocated matrix is most likely to
# get wrong, and it can be got wrong in either direction: naming rows `rbind()`
# left unnamed, or leaving unnamed the rows it fills in with an empty string
# because some other block carried names. Both are stated here rather than left
# to the fits, whose blocks carry no row names at all.

test_that("stack_psi_blocks() names no rows when no block does", {
  n <- 3L
  first <- matrix(seq_len(2L * n) + 0.5, nrow = 2L)
  second <- matrix(seq_len(n) - 0.5, nrow = 1L)

  stacked <- stack_psi_blocks(list(first, second), n)

  expect_identical(stacked, rbind(first, second))
  expect_null(rownames(stacked))
})

test_that("stack_psi_blocks() pads unnamed rows the way rbind does", {
  n <- 2L
  named <- matrix(
    c(1, 2, 3, 4),
    nrow = 2L,
    dimnames = list(c("mu0", "mu1"), NULL)
  )
  unnamed <- matrix(c(5, 6), nrow = 1L)

  stacked <- stack_psi_blocks(list(named, unnamed), n)

  expect_identical(stacked, rbind(named, unnamed))
  expect_identical(rownames(stacked), c("mu0", "mu1", ""))
})

# The buffer is allocated at the width the caller declares, so a block that is
# not that wide is a mistake in the block rather than something to fit in. R
# would recycle it into the rows instead, and the stack would carry a psi matrix
# whose values belong to no unit, so the width is checked rather than trusted.
#
# The refusal is worded from what the entry is and reported at the assembly the
# entries were handed to. A matrix entry is the only kind with columns to count;
# a bare vector is one row, so what is wrong with it is how many values it
# holds. Neither is refused from the `vapply()` closure the count is taken in,
# which is a frame no caller wrote and which names no block.

test_that("balancing_internal_error: a matrix block of the wrong width", {
  n <- 4L
  wide_enough <- matrix(seq_len(n) + 0.5, nrow = 1L)
  too_narrow <- matrix(c(1, 2), nrow = 1L)

  expect_balancing_error(stack_psi_blocks(list(wide_enough, too_narrow), n))
})

test_that("balancing_internal_error: a bare-vector block of the wrong length", {
  n <- 4L
  wide_enough <- matrix(seq_len(n) + 0.5, nrow = 1L)

  expect_balancing_error(stack_psi_blocks(list(wide_enough, c(1, 2)), n))
})

test_that("balancing_internal_error: a bare-vector block the reduction refuses", {
  n <- 4L
  wide_enough <- matrix(seq_len(n) + 0.5, nrow = 1L)

  expect_balancing_error(sum_psi_blocks(list(wide_enough, c(1, 2)), n))
})

# The two readings agree only while every block holds doubles. `rowSums()` on
# the assembled stack returns a double whatever the block's storage was, while
# the reduction takes a per-observation row with `sum()`, which returns an
# integer for an integer row and can overflow it to `NA`. No route builds such a
# row today, so the storage is refused where the width is rather than left to
# make the two readings differ.

test_that("balancing_internal_error: a block that does not hold doubles", {
  n <- 4L
  wide_enough <- matrix(seq_len(n) + 0.5, nrow = 1L)

  expect_balancing_error(sum_psi_blocks(list(wide_enough, seq_len(n)), n))
})

test_that("the assembly refuses a block that does not hold doubles", {
  n <- 4L
  wide_enough <- matrix(seq_len(n) + 0.5, nrow = 1L)

  expect_error(
    stack_psi_blocks(list(wide_enough, matrix(seq_len(n), nrow = 1L)), n),
    class = "balancing_internal_error"
  )
})

# The blocks a stack carries are not all matrices, and two kinds of them used to
# be. A route's mean rows were stacked into a block of their own before that
# block was copied into the destination, and its contrast rows were expanded
# into a block holding one value repeated across every column. Both are written
# a row at a time now, and the deterministic rows are written as the single
# value they repeat, so what the assembly returns has to be what stacking those
# blocks returned, to the bit. All four of the kinds a `.by` request puts in the
# stack are present here: the whole-sample mean and contrast rows, and the
# stratum rows of each.

test_that("stack_psi_blocks() writes rows where blocks were stacked before", {
  n <- 5L
  weight_block <- matrix(seq_len(2L * n) / 10, nrow = 2L)
  score_block <- matrix(seq_len(3L * n) * 2.5, nrow = 3L)
  mean_rows <- list(seq_len(n) + 0.25, seq_len(n) - 0.75)
  contrast_values <- c(rd = -0.5, "log(rr)" = 0.25, "log(or)" = 1.5)
  by_mean_rows <- list(seq_len(n) * 0.5, seq_len(n) * -0.5)
  by_contrast_values <- c("rd_g = lo" = 0.1, "rd_g = hi" = -0.2)

  stacked <- stack_psi_blocks(
    c(
      list(weight_block, score_block),
      mean_rows,
      constant_psi_rows(contrast_values),
      by_mean_rows,
      constant_psi_rows(by_contrast_values)
    ),
    n
  )

  expect_identical(
    stacked,
    rbind(
      weight_block,
      score_block,
      do.call(rbind, mean_rows),
      matrix(contrast_values, nrow = 3L, ncol = n),
      do.call(rbind, by_mean_rows),
      matrix(by_contrast_values, nrow = 2L, ncol = n)
    )
  )
})

# ---- The reduction the bread reads -----------------------------------------

# The bread never looks at the stack itself, only at its row sums, so the sums
# are taken block by block where the blocks are built and the S-by-n
# destination is never allocated on that path. That is worth having only while
# the two readings agree to the bit, which is what these state: the reduction of
# a list of blocks has to be the row sums of the matrix the same list assembles
# to, over the shapes a stack presents and over the absent and empty entries a
# route that carries no block of some kind passes.

test_that("sum_psi_blocks() reduces to the row sums of the assembled stack", {
  n <- 5L
  weight_block <- matrix(seq_len(2L * n) / 10, nrow = 2L)
  score_block <- matrix(seq_len(3L * n) * 2.5, nrow = 3L)
  mean_rows <- list(seq_len(n) + 0.25, seq_len(n) - 0.75)
  contrast_values <- c(rd = -0.5, "log(rr)" = 0.25, "log(or)" = 1.5)
  by_mean_rows <- list(seq_len(n) * 0.5, seq_len(n) * -0.5)
  by_contrast_values <- c("rd_g = lo" = 0.1, "rd_g = hi" = -0.2)

  blocks <- c(
    list(weight_block, score_block),
    mean_rows,
    constant_psi_rows(contrast_values),
    by_mean_rows,
    constant_psi_rows(by_contrast_values)
  )

  expect_identical(
    sum_psi_blocks(blocks, n),
    unname(rowSums(stack_psi_blocks(blocks, n)))
  )
})

test_that("sum_psi_blocks() drops absent and empty blocks as the assembly does", {
  n <- 4L
  weight_block <- matrix(seq_len(2L * n) / 10, nrow = 2L)
  empty_block <- matrix(numeric(0), nrow = 0L, ncol = n)
  blocks <- list(weight_block, NULL, empty_block, seq_len(n) + 0.5, 1.25)

  expect_identical(
    sum_psi_blocks(blocks, n),
    unname(rowSums(stack_psi_blocks(blocks, n)))
  )
})

# ---- The assembly a real fit performs --------------------------------------

# Four routes build four different sets of blocks: a binary exposure two mean
# rows and one block of contrasts, a categorical one a mean row per level and a
# block of contrasts per non-reference level, a `.by` request two further blocks
# after those, and a declared crossing the same means under contrasts written in
# the two treatments. The whole-sample routes leave the two subgroup blocks
# absent, so between them the four cover the blocks that are always there and
# the ones that are only sometimes.
#
# The two whole-sample routes are stated here. The other two are stated beside
# the fixtures they need, in test-ipw-by.R and test-ipw-joint.R.

test_that("a binary fit stacks its psi blocks as rbind would", {
  data <- ipw_deli_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_deli_outcome(
    y ~ exposure + x1,
    data,
    w,
    stats::binomial()
  )

  expect_stacked_psi_matches_rbind(ipw(fit, outcome_mod))
})

test_that("a categorical fit stacks its psi blocks as rbind would", {
  data <- ipw_deli_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_deli_outcome(
    y ~ exposure + x1,
    data,
    w,
    stats::binomial()
  )

  expect_stacked_psi_matches_rbind(ipw(fit, outcome_mod))
})
