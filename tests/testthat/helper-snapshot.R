expect_balancing_error <- function(expr) {
  testthat::expect_snapshot(
    error = TRUE,
    cnd_class = TRUE,
    expr
  )
}

expect_balancing_warning <- function(expr) {
  testthat::expect_snapshot(
    cnd_class = TRUE,
    expr
  )
}
