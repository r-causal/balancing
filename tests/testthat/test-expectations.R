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

# The same hole reaches every sibling comparison the suite makes on a reported
# column: `all(estimates$std.err > 0)` on a frame that does not carry `std.err`
# compares `NULL` against zero, gets `logical(0)`, and passes. Requiring the
# column to be present closes that, and requiring it to hold at least one value
# closes the case a present-but-empty column would still slip through.
test_that("expect_column_all() fails when the column is missing", {
  estimates <- data.frame(estimate = 1:3, std.err = c(0.1, 0.2, 0.3))

  expect_error(
    expect_column_all(estimates, "std.error", function(x) x > 0),
    class = "expectation_failure"
  )
  expect_column_all(estimates, "std.err", function(x) x > 0)
})

test_that("expect_column_all() fails on a column with no values", {
  expect_error(
    expect_column_all(
      data.frame(std.err = numeric(0)),
      "std.err",
      function(x) x > 0
    ),
    class = "expectation_failure"
  )
})

test_that("expect_column_all() fails when the predicate does not hold", {
  expect_error(
    expect_column_all(
      data.frame(std.err = c(0.1, -0.2)),
      "std.err",
      function(x) x > 0
    ),
    class = "expectation_failure"
  )
})

# Several call sites compare a column against a sibling column, and a missing
# sibling reopens the hole from the other side: `x < NULL` is `logical(0)`
# whatever `x` holds. The predicate is therefore required to answer with one
# value per row rather than only to answer TRUE everywhere it answers at all.
test_that("expect_column_all() fails when the predicate answers a short vector", {
  estimates <- data.frame(std.err = c(0.1, 0.2))

  expect_error(
    expect_column_all(estimates, "std.err", function(x) x < estimates$absent),
    class = "expectation_failure"
  )
})
