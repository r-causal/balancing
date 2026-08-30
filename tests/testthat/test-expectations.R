# Unit tests for the custom expectations in helper-expectations.R. An
# expectation that cannot fail reports nothing, so the failing cases are worth
# pinning on their own rather than only through the assertions they stand
# behind.
#
# The failures are caught as conditions rather than through `expect_failure()`,
# which reads exactly one expectation and `expect_finite_column()` makes two.

# `is.finite()` on an absent column returns `logical(0)` and `all(logical(0))`
# is TRUE, so the finiteness half alone would pass on a data frame that does not
# carry the column at all. The membership half is what makes a renamed or
# dropped reporting column a failure instead of silence.
test_that("expect_finite_column() fails when the column is missing", {
  estimates <- data.frame(estimate = 1:3, std.err = c(0.1, 0.2, 0.3))

  expect_error(
    expect_finite_column(estimates, "std.error"),
    class = "expectation_failure"
  )
  expect_finite_column(estimates, "std.err")
})

test_that("expect_finite_column() fails on a non-finite value", {
  expect_error(
    expect_finite_column(data.frame(std.err = c(0.1, NA)), "std.err"),
    class = "expectation_failure"
  )
  expect_error(
    expect_finite_column(data.frame(std.err = c(0.1, Inf)), "std.err"),
    class = "expectation_failure"
  )
})
