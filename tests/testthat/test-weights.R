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
