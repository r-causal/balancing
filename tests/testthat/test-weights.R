# The bw class is a sibling of propensity::psw under the shared causal_wts
# parent. These specs cover construction, the vctrs prototype labels, the
# coercion lattice with its classed downgrade warning, arithmetic preservation,
# and the ess() accessor shape.

# ---- Construction ---------------------------------------------------------

test_that("new_bw() builds a bw vector with the expected class stack", {
  w <- new_bw(c(0.5, 1, 1.5), estimand = "ate")
  expect_s3_class(
    w,
    c("bw", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )
  expect_equal(vctrs::vec_data(w), c(0.5, 1, 1.5))
  expect_identical(propensity::estimand(w), "ate")
})

test_that("bw() coerces and validates the input", {
  w <- bw(c(1L, 2L, 3L), estimand = "att")
  expect_s3_class(w, "bw")
  expect_type(vctrs::vec_data(w), "double")
  expect_identical(propensity::estimand(w), "att")
})

test_that("as_bw() coerces a plain numeric vector", {
  w <- as_bw(c(1, 2), estimand = "atu")
  expect_true(is_bw(w))
  expect_identical(propensity::estimand(w), "atu")
})

test_that("is_bw() identifies bw vectors only", {
  expect_true(is_bw(bw(c(1, 2), estimand = "ate")))
  expect_false(is_bw(c(1, 2)))
})

test_that("propensity::is_causal_wt() recognizes bw vectors", {
  expect_true(propensity::is_causal_wt(bw(c(1, 2), estimand = "ate")))
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
  expect_identical(propensity::estimand(combined), "ate")
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
  expect_identical(propensity::estimand(back), "ate")
})

# ---- Arithmetic -----------------------------------------------------------

test_that("arithmetic with a scalar preserves the bw class and estimand", {
  w <- bw(c(1, 2, 3), estimand = "att")
  scaled <- w * 2
  expect_true(is_bw(scaled))
  expect_identical(propensity::estimand(scaled), "att")
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
  expect_identical(propensity::estimand(from_integer), "ate")
})

test_that("combining an integer with a bw yields a plain double", {
  w <- bw(c(1, 2), estimand = "ate")
  combined <- vctrs::vec_c(1L, w)
  expect_type(combined, "double")
  expect_false(is_bw(combined))
  expect_equal(combined, c(1, 1, 2))
})

test_that("combining a psw before a bw warns and downgrades", {
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
  expect_identical(propensity::estimand(sum), "ate")
  expect_equal(vctrs::vec_data(sum), c(1.5, 3, 4.5))
})

test_that("adding two bw vectors with different estimands records both", {
  x <- bw(c(1, 2), estimand = "ate")
  y <- bw(c(3, 4), estimand = "att")
  combined <- x + y
  expect_true(is_bw(combined))
  expect_identical(propensity::estimand(combined), "ate, att")
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

# ---- Subsetting -----------------------------------------------------------

test_that("subsetting a bw preserves the class and metadata", {
  w <- bw(c(1, 2, 3, 4), estimand = "ate")
  first_two <- w[1:2]
  expect_true(is_bw(first_two))
  expect_identical(propensity::estimand(first_two), "ate")
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

# ---- ess() ----------------------------------------------------------------

test_that("ess() returns a group, n, ess tibble", {
  data <- sim_binary(n = 200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = entropy_balance(),
    estimand = "ate"
  )
  ess_tbl <- ess(fit)
  expect_true(all(c("group", "n", "ess") %in% names(ess_tbl)))
  expect_true(all(ess_tbl$ess <= ess_tbl$n + 1e-8))
})
