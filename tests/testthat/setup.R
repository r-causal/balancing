# Suppress partial match warnings during tests
op <- options(
  warnPartialMatchDollar = FALSE,
  warnPartialMatchArgs = FALSE,
  warnPartialMatchAttr = FALSE
)

withr::defer(options(op), teardown_env())
