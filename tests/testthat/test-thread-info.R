test_that("thread_info() round-trips a list from the Rust core", {
  info <- thread_info()

  expect_type(info, "list")
  expect_named(info, c("available", "cap_source"))
})

test_that("thread_info() reports at least one available thread", {
  info <- thread_info()

  expect_type(info$available, "integer")
  expect_length(info$available, 1)
  expect_gte(info$available, 1L)
})

test_that("thread_info() names the source that capped the thread count", {
  info <- thread_info()

  expect_type(info$cap_source, "character")
  expect_length(info$cap_source, 1)
  expect_true(info$cap_source %in% c("system", "OMP_THREAD_LIMIT"))
})
