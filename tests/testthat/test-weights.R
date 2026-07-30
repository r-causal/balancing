# The bw class is a sibling of propensity::psw under the shared causal_wts
# parent. These specs cover construction, the vctrs prototype labels, the
# coercion lattice with its classed downgrade warning, and arithmetic
# preservation.

# ---- Construction ---------------------------------------------------------

test_that("new_bw() builds a bw vector with the expected class stack", {
  w <- new_bw(c(0.5, 1, 1.5), estimand = "ate")
  expect_s3_class(
    w,
    c("bw", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )
  expect_equal(vctrs::vec_data(w), c(0.5, 1, 1.5))
  expect_identical(estimand(w), "ate")
})

test_that("bw() coerces and validates the input", {
  w <- bw(c(1L, 2L, 3L), estimand = "att")
  expect_s3_class(w, "bw")
  expect_type(vctrs::vec_data(w), "double")
  expect_identical(estimand(w), "att")
})

test_that("as_bw() coerces a plain numeric vector", {
  w <- as_bw(c(1, 2), estimand = "atu")
  expect_true(is_bw(w))
  expect_identical(estimand(w), "atu")
})

test_that("is_bw() identifies bw vectors only", {
  expect_true(is_bw(bw(c(1, 2), estimand = "ate")))
  expect_false(is_bw(c(1, 2)))
})

test_that("is_causal_wt() recognizes bw vectors", {
  expect_true(is_causal_wt(bw(c(1, 2), estimand = "ate")))
})

# ---- Prototype labels -----------------------------------------------------

test_that("the prototype abbreviation includes the estimand", {
  w <- bw(c(1, 2), estimand = "att")
  expect_identical(vctrs::vec_ptype_abbr(w), "bw{att}")
})

test_that("the full prototype includes the estimand", {
  w <- bw(c(1, 2), estimand = "att")
  expect_match(vctrs::vec_ptype_full(w), "att")
})

# ---- Coercion lattice -----------------------------------------------------

test_that("combining bw vectors with matching estimands preserves the class", {
  x <- bw(c(1, 2), estimand = "ate")
  y <- bw(c(3, 4), estimand = "ate")
  combined <- vctrs::vec_c(x, y)
  expect_true(is_bw(combined))
  expect_identical(estimand(combined), "ate")
})

test_that("combining bw vectors with mismatched estimands warns and downgrades", {
  x <- bw(c(1, 2), estimand = "ate")
  y <- bw(c(3, 4), estimand = "att")
  expect_warning(
    combined <- vctrs::vec_c(x, y),
    class = "balancing_class_downgrade_warning"
  )
  expect_type(combined, "double")
  expect_false(is_bw(combined))
})

test_that("combining a bw with a psw warns and downgrades", {
  skip_if_not_installed("propensity")
  x <- bw(c(1, 2), estimand = "ate")
  y <- propensity::psw(c(3, 4), estimand = "ate")
  expect_warning(
    combined <- vctrs::vec_c(x, y),
    class = "balancing_class_downgrade_warning"
  )
  expect_type(combined, "double")
})

test_that("combining a bw with character warns and downgrades to character", {
  x <- bw(c(1, 2), estimand = "ate")
  expect_warning(
    combined <- vctrs::vec_c(x, "label"),
    class = "balancing_class_downgrade_warning"
  )
  expect_type(combined, "character")
  expect_false(is_bw(combined))
})

test_that("bw casts to character", {
  x <- bw(c(1, 2, 3), estimand = "ate")
  as_character <- vctrs::vec_cast(x, character())
  expect_type(as_character, "character")
  expect_identical(as_character, c("1", "2", "3"))
})

test_that("the character common type with a bw is character in both orders", {
  x <- bw(c(1, 2), estimand = "ate")
  suppressWarnings({
    expect_identical(vctrs::vec_ptype2(x, character()), character())
    expect_identical(vctrs::vec_ptype2(character(), x), character())
  })
})

test_that("bw casts to and from double", {
  x <- bw(c(1, 2, 3), estimand = "ate")
  as_double <- vctrs::vec_cast(x, double())
  expect_type(as_double, "double")
  expect_equal(as_double, c(1, 2, 3))

  back <- vctrs::vec_cast(c(1, 2, 3), x)
  expect_true(is_bw(back))
  expect_identical(estimand(back), "ate")
})

# ---- Arithmetic -----------------------------------------------------------

test_that("arithmetic with a scalar preserves the bw class and estimand", {
  w <- bw(c(1, 2, 3), estimand = "att")
  scaled <- w * 2
  expect_true(is_bw(scaled))
  expect_identical(estimand(scaled), "att")
  expect_equal(vctrs::vec_data(scaled), c(2, 4, 6))
})

test_that("normalizing a bw vector preserves the class", {
  w <- bw(c(1, 2, 3), estimand = "ate")
  normalized <- w / sum(w)
  expect_true(is_bw(normalized))
})

# ---- Prototype labels without an estimand ---------------------------------

test_that("the prototype labels report an unknown estimand when it is NULL", {
  w <- bw(c(1, 2))
  expect_identical(vctrs::vec_ptype_abbr(w), "bw")
  expect_identical(vctrs::vec_ptype_full(w), "bw{estimand = unknown}")
})

# ---- Coercion lattice with double and integer -----------------------------

test_that("the common type of a bw with double or integer is double", {
  w <- bw(c(1, 2), estimand = "ate")
  expect_identical(vctrs::vec_ptype2(w, double()), double())
  expect_identical(vctrs::vec_ptype2(double(), w), double())
  expect_identical(vctrs::vec_ptype2(w, integer()), double())
  expect_identical(vctrs::vec_ptype2(integer(), w), double())
})

test_that("bw casts to and from integer", {
  w <- bw(c(1, 2, 3), estimand = "ate")
  as_integer <- vctrs::vec_cast(w, integer())
  expect_identical(as_integer, c(1L, 2L, 3L))

  from_integer <- vctrs::vec_cast(c(1L, 2L, 3L), w)
  expect_true(is_bw(from_integer))
  expect_identical(estimand(from_integer), "ate")
})

test_that("combining an integer with a bw yields a plain double", {
  w <- bw(c(1, 2), estimand = "ate")
  combined <- vctrs::vec_c(1L, w)
  expect_type(combined, "double")
  expect_false(is_bw(combined))
  expect_equal(combined, c(1, 1, 2))
})

test_that("combining a psw before a bw warns and downgrades", {
  skip_if_not_installed("propensity")
  x <- propensity::psw(c(1, 2), estimand = "ate")
  y <- bw(c(3, 4), estimand = "ate")
  expect_warning(
    combined <- vctrs::vec_c(x, y),
    class = "balancing_class_downgrade_warning"
  )
  expect_type(combined, "double")
  expect_false(is_bw(combined))
})

# ---- Arithmetic branches --------------------------------------------------

test_that("adding two bw vectors with a matching estimand keeps the estimand", {
  x <- bw(c(1, 2, 3), estimand = "ate")
  y <- bw(c(0.5, 1, 1.5), estimand = "ate")
  sum <- x + y
  expect_true(is_bw(sum))
  expect_identical(estimand(sum), "ate")
  expect_equal(vctrs::vec_data(sum), c(1.5, 3, 4.5))
})

test_that("adding two bw vectors with different estimands records both", {
  x <- bw(c(1, 2), estimand = "ate")
  y <- bw(c(3, 4), estimand = "att")
  combined <- x + y
  expect_true(is_bw(combined))
  expect_identical(estimand(combined), "ate, att")
  expect_equal(vctrs::vec_data(combined), c(4, 6))
})

test_that("arithmetic with a bare numeric preserves the class in either order", {
  w <- bw(c(1, 2, 3), estimand = "att")
  left <- w * 2
  right <- 2 * w
  expect_true(is_bw(left))
  expect_true(is_bw(right))
  expect_equal(vctrs::vec_data(left), c(2, 4, 6))
  expect_equal(vctrs::vec_data(right), c(2, 4, 6))
})

test_that("arithmetic with an integer preserves the class", {
  w <- bw(c(1, 2, 3), estimand = "ate")
  shifted <- w + 1L
  expect_true(is_bw(shifted))
  expect_equal(vctrs::vec_data(shifted), c(2, 3, 4))
})

test_that("unary minus and plus preserve the bw class", {
  w <- bw(c(1, -2, 3), estimand = "ate")
  negated <- -w
  expect_true(is_bw(negated))
  expect_equal(vctrs::vec_data(negated), c(-1, 2, -3))
  expect_identical(+w, w)
})

test_that("arithmetic with an unsupported type errors", {
  w <- bw(c(1, 2), estimand = "ate")
  expect_error(w + "a")
})

# ---- Math, Summary, min, max ----------------------------------------------

test_that("cumulative math preserves the class while reductions strip it", {
  w <- bw(c(1, 2, 3), estimand = "ate")
  running <- cumsum(w)
  expect_true(is_bw(running))
  expect_equal(vctrs::vec_data(running), c(1, 3, 6))

  expect_false(is_bw(sqrt(w)))
  expect_false(is_bw(sum(w)))
  expect_equal(sum(w), 6)
})

test_that("Summary, min, and max operate on the underlying data", {
  w <- bw(c(3, 1, 2), estimand = "ate")
  expect_equal(min(w), 1)
  expect_equal(max(w), 3)
  expect_equal(range(w), c(1, 3))
  expect_false(is_bw(min(w)))
})

test_that("median() and quantile() operate on the underlying data", {
  w <- bw(c(3, 1, 2, 4), estimand = "ate")
  expect_equal(median(w), median(c(3, 1, 2, 4)))
  expect_false(is_bw(median(w)))
  expect_equal(
    quantile(w, probs = c(0.25, 0.75)),
    quantile(c(3, 1, 2, 4), probs = c(0.25, 0.75))
  )
  expect_false(is_bw(quantile(w)))
})

# ---- Subsetting -----------------------------------------------------------

test_that("subsetting a bw preserves the class and metadata", {
  w <- bw(c(1, 2, 3, 4), estimand = "ate")
  first_two <- w[1:2]
  expect_true(is_bw(first_two))
  expect_identical(estimand(first_two), "ate")
  expect_equal(vctrs::vec_data(first_two), c(1, 2))
})

test_that("subsetting with no index returns the whole vector", {
  w <- bw(c(1, 2, 3), estimand = "ate")
  expect_true(is_bw(w[]))
  expect_equal(vctrs::vec_data(w[]), c(1, 2, 3))
})

test_that("matrix subsetting drops to the underlying data", {
  w <- bw(c(10, 20, 30), estimand = "ate")
  index <- matrix(c(1L, 3L), ncol = 1)
  picked <- w[index]
  expect_false(is_bw(picked))
  expect_equal(picked, c(10, 30))
})

# ---- Restoration and the groups attribute ---------------------------------

# `balance()` records which rows belong to each exposure level on the weight
# vector it builds, as a `groups` attribute, and `ipw()` reads the level order
# back off `fit@weights` to name the reference level every contrast is measured
# against. Those entries are row positions into that one vector, so they belong
# to the fit's own weights and to nothing derived from them. Restoration drops
# them in every case, and the tests below pin the rule at each operation that
# routes through `vec_restore()`, so a regression names the route it came in on.

test_that("vec_restore drops the groups attribute", {
  w <- new_bw(
    c(1, 2, 3, 4),
    estimand = "ate",
    groups = list(`0` = 1:2, `1` = 3:4)
  )
  restored <- vctrs::vec_restore(c(2, 4, 6, 8), w)
  expect_true(is_bw(restored))
  expect_identical(estimand(restored), "ate")
  expect_null(attr(restored, "groups"))
})

test_that("arithmetic and cumulative math drop the groups attribute", {
  w <- new_bw(
    c(1, 2, 3, 4),
    estimand = "ate",
    groups = list(`0` = 1:2, `1` = 3:4)
  )
  expect_identical(estimand(w * 2), "ate")
  expect_null(attr(w * 2, "groups"))
  expect_null(attr(-w, "groups"))
  expect_null(attr(cumsum(w), "groups"))
})

# Reordering is the case the rule is unconditional for. These come back at the
# size they started at, so a restoration that kept the attribute whenever the
# size held would re-attach positions the rows no longer match: `rev(w)` would
# claim level `0` owns what are now the last two weights. `vec_restore()` is
# handed the restored data and the object it came from, never the index that
# reordered them, so it cannot recognize this and rescale the positions.
test_that("reordering drops the groups attribute at the original size", {
  w <- new_bw(
    c(4, 3, 2, 1),
    estimand = "ate",
    groups = list(`0` = 1:2, `1` = 3:4)
  )

  reversed <- rev(w)
  expect_true(is_bw(reversed))
  expect_identical(vctrs::vec_size(reversed), 4L)
  expect_null(attr(reversed, "groups"))

  expect_null(attr(sort(w), "groups"))

  # Repeating positions holds the size while duplicating rows.
  expect_null(attr(w[c(1, 1, 2, 2)], "groups"))
})

# Re-attaching row positions to a shorter vector would describe rows that are no
# longer there, which is the same defect the reordering case has, arrived at
# from the other direction.
test_that("slicing drops the groups attribute rather than keeping stale rows", {
  w <- new_bw(
    c(1, 2, 3, 4),
    estimand = "ate",
    groups = list(`0` = 1:2, `1` = 3:4)
  )
  sliced <- w[1:2]
  expect_true(is_bw(sliced))
  expect_identical(estimand(sliced), "ate")
  expect_null(attr(sliced, "groups"))
})

# Combining restores onto the common prototype rather than onto either input,
# and `vec_ptype2.bw.bw()` builds that prototype from the estimand alone, so the
# `to` a combination restores through carries no `groups` to begin with. Under a
# uniform drop the two agree, and the test holds the agreement in place: a
# prototype that started carrying `groups` would be caught here.
test_that("combining or repeating a bw drops the groups attribute", {
  w <- new_bw(
    c(1, 2, 3, 4),
    estimand = "ate",
    groups = list(`0` = 1:2, `1` = 3:4)
  )

  combined <- vctrs::vec_c(w, w)
  expect_true(is_bw(combined))
  expect_identical(estimand(combined), "ate")
  expect_null(attr(combined, "groups"))

  single <- vctrs::vec_c(w)
  expect_identical(vctrs::vec_size(single), 4L)
  expect_null(attr(single, "groups"))

  expect_null(attr(rep(w, 2), "groups"))
})

# A fit carries no row positions on its weight vector at all: what the exposure
# levels are is recorded as a property of the fit, and which rows take each level
# is in the data. Both branches of `weights()` are pinned here because the
# symmetry between them is the contract: the composing branch would drop any extra
# attribute by way of `vec_restore.bw()`, and the branch that returns the
# balancing weights alone hands back the fit's own vector, so the accessor has one
# return shape whether or not the fit was given sampling weights.
test_that("weights() carries the estimand and no fit metadata", {
  data <- sim_binary(200)
  data$sampling <- rep(c(0.8, 1.2), length.out = nrow(data))
  sampled <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate",
    sampling_weights = sampling
  )
  plain <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )

  # The exposure levels are on the fit, not on the weights.
  expect_identical(sampled@exposure_levels, c("0", "1"))
  expect_identical(plain@exposure_levels, c("0", "1"))
  expect_null(attr(sampled@weights, "groups"))
  expect_null(attr(plain@weights, "groups"))

  # The class and the estimand are the whole of what a weight vector carries, on
  # the fit and out of the accessor alike.
  expect_setequal(names(attributes(sampled@weights)), c("class", "estimand"))
  expect_setequal(
    names(attributes(stats::weights(sampled))),
    c("class", "estimand")
  )
  expect_setequal(
    names(attributes(stats::weights(
      sampled,
      include_sampling_weights = FALSE
    ))),
    c("class", "estimand")
  )
  expect_null(attr(stats::weights(sampled), "groups"))
  expect_null(
    attr(stats::weights(sampled, include_sampling_weights = FALSE), "groups")
  )
  expect_null(attr(stats::weights(plain), "groups"))

  # What the accessor does keep.
  expect_true(is_bw(stats::weights(sampled)))
  expect_identical(estimand(stats::weights(sampled)), "ate")
  expect_identical(vctrs::vec_size(stats::weights(plain)), nrow(data))
})
