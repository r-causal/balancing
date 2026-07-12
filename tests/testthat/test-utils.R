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

# ---- Thread resolution ----------------------------------------------------

test_that("resolve_threads() honors an explicit thread count", {
  expect_identical(resolve_threads(4), 4L)
  # A count below one is floored to a single thread.
  expect_identical(resolve_threads(0), 1L)
})

test_that("resolve_threads() falls back to the balancing.threads option", {
  withr::local_options(balancing.threads = 3)
  expect_identical(resolve_threads(), 3L)
})

test_that("automatic_threads() returns two under R CMD check", {
  withr::local_envvar(`_R_CHECK_LIMIT_CORES_` = "TRUE")
  expect_identical(automatic_threads(), 2L)
})

test_that("automatic_threads() treats an unknown core count as one", {
  withr::local_envvar(
    `_R_CHECK_LIMIT_CORES_` = NA,
    OMP_THREAD_LIMIT = NA,
    OMP_NUM_THREADS = NA
  )
  testthat::local_mocked_bindings(
    detectCores = function(...) NA_integer_,
    .package = "parallel"
  )
  expect_identical(automatic_threads(), 1L)
})

test_that("env_thread_cap() parses a positive integer and rejects the rest", {
  withr::local_envvar(OMP_THREAD_LIMIT = "4")
  expect_identical(env_thread_cap("OMP_THREAD_LIMIT"), 4L)

  withr::local_envvar(OMP_THREAD_LIMIT = "not-a-number")
  expect_identical(env_thread_cap("OMP_THREAD_LIMIT"), Inf)

  withr::local_envvar(OMP_THREAD_LIMIT = NA)
  expect_identical(env_thread_cap("OMP_THREAD_LIMIT"), Inf)
})
