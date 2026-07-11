# Suppress partial match warnings during tests
op <- options(
  warnPartialMatchDollar = FALSE,
  warnPartialMatchArgs = FALSE,
  warnPartialMatchAttr = FALSE
)

withr::defer(options(op), teardown_env())

# Silence informational alerts by default so snapshots stay stable; individual
# tests that assert on an announcement opt back in with local_options().
quiet_op <- options(balancing.quiet = TRUE)

withr::defer(options(quiet_op), teardown_env())
