# The sampling-weight validator guards type, length, missingness, and sign
# before a fit begins. Each branch raises a classed balancing error; the
# user-facing messages are snapshotted so a change in wording or class is
# caught. A valid vector passes through invisibly and unchanged.

test_that("validate_sampling_weights() returns valid weights invisibly", {
  weights <- c(0.5, 1, 1.5)
  expect_identical(validate_sampling_weights(weights, n = 3), weights)
})

test_that("validate_sampling_weights() rejects a non-numeric vector", {
  expect_balancing_error(
    validate_sampling_weights(c("a", "b"), n = 2)
  )
})

test_that("validate_sampling_weights() rejects a length mismatch", {
  expect_balancing_error(
    validate_sampling_weights(c(1, 2), n = 3)
  )
})

test_that("validate_sampling_weights() rejects missing values", {
  expect_balancing_error(
    validate_sampling_weights(c(1, NA, 1), n = 3)
  )
})

test_that("validate_sampling_weights() rejects negative values", {
  expect_balancing_error(
    validate_sampling_weights(c(1, -1, 1), n = 3)
  )
})
