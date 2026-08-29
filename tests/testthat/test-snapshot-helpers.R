# Unit tests for the snapshot transform in helper-snapshot.R. The transform runs
# on every cli snapshot in the suite, so what it rewrites and what it leaves
# byte-identical is worth pinning on its own rather than only through the
# snapshots it feeds.

# A balance statistic the fit drove to zero prints as whatever residual the
# solver's floating-point path left behind, and that residual differs on every
# platform down to its exponent. It states nothing about the fit beyond "this
# term balanced", so it is replaced by a placeholder that says exactly that.
test_that("scrub_platform_values() replaces a numerically zero balance value", {
  line <- "1   x1 moment       smd     1  0.552 1.24e-11         0"
  expect_identical(
    scrub_platform_values(line),
    "1   x1 moment       smd     1  0.552 <1e-7         0"
  )
})

# The rule is written on the exponent rather than on the value, so the largest
# number it can match is just under 1e-7 rather than just under 1e-8. The
# placeholder states the cutoff the rule actually applies.
test_that("scrub_platform_values() states the cutoff its exponent rule reaches", {
  line <- "1   x1 moment       smd     1  0.552 9.9e-8         0"
  expect_identical(
    scrub_platform_values(line),
    "1   x1 moment       smd     1  0.552 <1e-7         0"
  )
})

# An imbalance the fit drove to an exact zero renders as a bare "0" rather than
# in exponent form, and it says exactly what a residual of 1e-17 says. The two
# have to record the same placeholder, or which of them a platform reaches
# becomes the difference between a passing and a failing snapshot.
test_that("scrub_platform_values() replaces an exact zero on the largest-imbalance line", {
  line <- "Largest imbalance: 0 (standardized mean difference)"
  expect_identical(
    scrub_platform_values(line),
    "Largest imbalance: <1e-7 (standardized mean difference)"
  )
})

# The rule is on that line alone. A zero anywhere else in a printed fit is a
# count, a tolerance, or a column the fit reports as zero by construction, and
# each of those is portable and worth keeping.
test_that("scrub_platform_values() leaves a zero on any other line alone", {
  line <- "Constraints: 2 terms (tolerance 0)"
  expect_identical(scrub_platform_values(line), line)
})

# A tolerance a condition message reports is a constant the package chose, the
# same on every platform, and it is often the thing the message is about. It sits
# in the magnitude range the near-zero rule scrubs, so the rule has to step
# around it.
test_that("scrub_platform_values() leaves a tolerance a message reports alone", {
  lines <- c(
    "The solve met tolerance 1e-10.",
    "The solve met tolerance = 1e-10.",
    "The solve met tolerance of 1e-10.",
    "The solve met Tolerance 1e-10."
  )
  expect_identical(scrub_platform_values(lines), lines)
})

# The exemption is on the word before the number rather than on the magnitude,
# so a residual of the same size somewhere else is still scrubbed.
test_that("scrub_platform_values() still scrubs a near-zero the word does not cover", {
  expect_identical(
    scrub_platform_values("The solve left 1e-10 behind."),
    "The solve left <1e-7 behind."
  )
})

# A method that balances a term only approximately reports a real distance, and
# that number is the point of the snapshot, so it has to survive untouched.
test_that("scrub_platform_values() leaves a measurable balance value alone", {
  line <- "2   x2 moment       smd     1  0.416 0.00511         0            FALSE"
  expect_identical(scrub_platform_values(line), line)
})
