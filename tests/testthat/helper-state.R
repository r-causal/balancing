# The package records the physical core count for the session, so a test that
# wants to observe the reading has to clear the record on the way in and again
# on the way out. Clearing it on the way out matters as much: a test that mocks
# `parallel::detectCores()` would otherwise leave the mocked answer behind for
# every fit that follows it in the same run.
local_core_count_reset <- function(.env = parent.frame()) {
  reset_physical_cores()
  withr::defer(reset_physical_cores(), envir = .env)
  invisible(NULL)
}
