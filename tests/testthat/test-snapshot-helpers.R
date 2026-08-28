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
    "1   x1 moment       smd     1  0.552 <1e-8         0"
  )
})

# A method that balances a term only approximately reports a real distance, and
# that number is the point of the snapshot, so it has to survive untouched.
test_that("scrub_platform_values() leaves a measurable balance value alone", {
  line <- "2   x2 moment       smd     1  0.416 0.00511         0            FALSE"
  expect_identical(scrub_platform_values(line), line)
})
