# Exposure-type resolution reads a type off the data through causalgenerics,
# honors an explicit type over the detection heuristics, and refuses only a
# declaration the data cannot represent at all or one the method cannot fit.
# The detection specs run through the resolver rather than against
# causalgenerics directly: what balancing depends on is the type a fit resolves
# to, and the resolver is where balancing decides it.

# The type the resolver reads off a vector. `bw_entropy()` fits all three types,
# so the method gate never intervenes and the answer is the detection alone.
detected_type <- function(exposure) {
  resolve_exposure_type("auto", exposure, bw_entropy())
}

# ---- Detected types -------------------------------------------------------

test_that("a two-valued exposure is binary regardless of storage type", {
  expect_identical(detected_type(c(0L, 1L, 0L, 1L)), "binary")
  expect_identical(detected_type(c(TRUE, FALSE, TRUE)), "binary")
  expect_identical(detected_type(factor(c("a", "b", "a"))), "binary")
})

test_that("a many-level factor or character exposure is categorical", {
  expect_identical(detected_type(factor(c("a", "b", "c", "a"))), "categorical")
  expect_identical(detected_type(c("a", "b", "c", "a")), "categorical")
})

test_that("a numeric exposure with few distinct values is categorical", {
  # Five distinct values in 100 observations is a 5 percent unique share, below
  # the 20 percent categorical threshold.
  withr::local_seed(1)
  exposure <- sample(1:5, 100, replace = TRUE)
  expect_identical(detected_type(exposure), "categorical")
})

test_that("a numeric exposure with many distinct values is continuous", {
  withr::local_seed(1)
  expect_identical(detected_type(stats::rnorm(100)), "continuous")
})

test_that("a missing value is not a level of the exposure", {
  # Detection counts observed values, so a two-valued exposure that also
  # carries missingness is binary. Missing values used to count as a level of
  # their own, which read this vector as continuous and a factor with one
  # missing value as categorical. `balance()` refuses a missing exposure before
  # any of this runs, so the rule is reachable only through the resolver.
  expect_identical(detected_type(c(0, 1, NA)), "binary")
  expect_identical(detected_type(factor(c("a", "b", NA))), "binary")
})

test_that("a one-level factor or character exposure is binary", {
  # A factor or character exposure turns categorical once it carries more than
  # two observed values, so a single level reads as binary rather than
  # categorical. `balance()` refuses a one-level exposure either way.
  expect_identical(detected_type(factor(c("a", "a", "a"))), "binary")
  expect_identical(detected_type(c("a", "a", "a")), "binary")
})

# ---- Explicit types -------------------------------------------------------

test_that("an explicit categorical type is honored on a factor exposure", {
  resolved <- resolve_exposure_type(
    "categorical",
    factor(c("a", "b", "c", "a")),
    bw_entropy()
  )
  expect_identical(resolved, "categorical")
})

test_that("an explicit continuous type is honored on a continuous exposure", {
  withr::local_seed(1)
  resolved <- resolve_exposure_type(
    "continuous",
    stats::rnorm(100),
    bw_entropy()
  )
  expect_identical(resolved, "continuous")
})

test_that("an explicit continuous type wins over the categorical heuristic", {
  # Ten distinct doses in 200 observations is a 5 percent unique share, which the
  # heuristic reads as categorical. An explicit type is the caller's declaration
  # of how the exposure is modeled, so it decides the fit.
  exposure <- withr::with_seed(
    1,
    sample(seq(10, 100, by = 10), 200, replace = TRUE)
  )
  expect_identical(detected_type(exposure), "categorical")
  expect_identical(
    resolve_exposure_type("continuous", exposure, bw_entropy()),
    "continuous"
  )
})

test_that("an explicit continuous type wins over two-level detection", {
  # A numeric exposure taking two values is still a dose the caller may model as
  # continuous, so the declaration stands.
  expect_identical(
    resolve_exposure_type("continuous", rep(c(0, 1), 50), bw_entropy()),
    "continuous"
  )
})

test_that("an explicit categorical type wins over continuous detection", {
  exposure <- withr::with_seed(1, stats::rnorm(20))
  expect_identical(
    resolve_exposure_type("categorical", exposure, bw_entropy()),
    "categorical"
  )
})

# ---- The announcement -----------------------------------------------------

test_that("auto resolution announces the detected type once", {
  withr::local_options(balancing.quiet = FALSE)
  announcements <- capture_messages(
    resolve_exposure_type("auto", c(0L, 1L, 0L, 1L), bw_entropy())
  )
  expect_length(announcements, 1L)
  expect_match(announcements, "Treating `.exposure` as binary", fixed = TRUE)
})

test_that("an explicit type announces nothing, whether it stands or not", {
  withr::local_options(balancing.quiet = FALSE)
  expect_length(
    capture_messages(
      resolve_exposure_type("binary", c(0L, 1L, 0L, 1L), bw_entropy())
    ),
    0L
  )
  expect_length(
    capture_messages(
      expect_error(
        resolve_exposure_type(
          "binary",
          withr::with_seed(1, stats::rnorm(100)),
          bw_entropy()
        ),
        class = "causalgenerics_forced_exposure_type"
      )
    ),
    0L
  )
})

# ---- Refusals -------------------------------------------------------------

test_that("a forced type the data contradict raises a classed error", {
  withr::local_seed(1)
  expect_error(
    resolve_exposure_type(
      "binary",
      stats::rnorm(100),
      bw_entropy()
    ),
    class = "causalgenerics_forced_exposure_type"
  )
})

test_that("a continuous type on a non-numeric exposure is a classed error", {
  # Structural impossibility rather than heuristic disagreement: a factor or
  # character exposure carries no dose to correlate the covariates against.
  expect_error(
    resolve_exposure_type(
      "continuous",
      factor(c("a", "b", "c", "a")),
      bw_entropy()
    ),
    class = "causalgenerics_forced_exposure_type"
  )
  expect_error(
    resolve_exposure_type(
      "continuous",
      c("a", "b", "c", "a"),
      bw_entropy()
    ),
    class = "causalgenerics_forced_exposure_type"
  )
})

test_that("a type the method cannot fit is balancing's own refusal", {
  # The structural refusals belong to causalgenerics; which exposure types a
  # method can fit is balancing's own question, so its refusal keeps balancing's
  # condition class.
  expect_error(
    resolve_exposure_type(
      "continuous",
      withr::with_seed(1, stats::rnorm(100)),
      bw_ipt()
    ),
    class = "balancing_exposure_type_error"
  )
})
