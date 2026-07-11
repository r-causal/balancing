test_that("abort() signals a classed balancing error", {
  expect_error(
    abort("Something went wrong."),
    class = "balancing_error"
  )
})

test_that("abort() attaches the requested subclass", {
  expect_error(
    abort("Bad input.", error_class = "balancing_input_error"),
    class = "balancing_input_error"
  )

  cnd <- rlang::catch_cnd(
    abort("Bad input.", error_class = "balancing_input_error")
  )
  expect_s3_class(cnd, c("balancing_input_error", "balancing_error"))
})

test_that("abort() formats messages with cli styling", {
  expect_balancing_error(
    abort(
      c(
        "The weights did not converge.",
        i = "Try increasing the number of iterations."
      ),
      error_class = "balancing_convergence_error"
    )
  )
})

test_that("warn() signals a classed balancing warning", {
  expect_warning(
    warn("Heads up."),
    class = "balancing_warning"
  )
})

test_that("warn() attaches the requested subclass", {
  expect_warning(
    warn("Heads up.", warning_class = "balancing_coercion_warning"),
    class = "balancing_coercion_warning"
  )

  cnd <- rlang::catch_cnd(
    warn("Heads up.", warning_class = "balancing_coercion_warning")
  )
  expect_s3_class(cnd, c("balancing_coercion_warning", "balancing_warning"))
})

test_that("warn() formats messages with cli styling", {
  expect_balancing_warning(
    warn(
      c(
        "Some weights were negative.",
        i = "They were set to zero."
      ),
      warning_class = "balancing_negative_weight_warning"
    )
  )
})

test_that("alert_info() prints when balancing.quiet is FALSE", {
  withr::local_options(balancing.quiet = FALSE)

  expect_message(alert_info("A helpful note."))
})

test_that("alert_info() stays silent when balancing.quiet is TRUE", {
  withr::local_options(balancing.quiet = TRUE)

  expect_silent(alert_info("A helpful note."))
})
