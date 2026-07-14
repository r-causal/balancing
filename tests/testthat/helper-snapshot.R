# Snapshot helpers that forward the caller's expression to expect_snapshot() so
# the recorded Code block shows the failing call, not the literal `expr` token.
expect_balancing_error <- function(expr) {
  rlang::inject(testthat::expect_snapshot(
    error = TRUE,
    cnd_class = TRUE,
    !!rlang::enquo(expr)
  ))
}

expect_balancing_warning <- function(expr) {
  rlang::inject(testthat::expect_snapshot(
    cnd_class = TRUE,
    !!rlang::enquo(expr)
  ))
}
