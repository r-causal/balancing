# balance() is the orchestrator: it resolves the exposure and covariate
# selections with tidyselect, detects and announces the exposure type, validates
# the estimand and focal level against the method's support, evaluates sampling
# weights, and dispatches to the method's fit. These specs pin the plumbing and
# the classed validation errors; the statistical behavior lives in the
# per-method files and the numerical oracle in expect_balanced().

# ---- Selection ------------------------------------------------------------

test_that("balance() resolves the exposure and covariates with tidyselect", {
  data <- sim_binary(n = 200)
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())

  expect_true(S7::S7_inherits(fit, balancing))
  expect_identical(fit@exposure, "exposure")
  expect_identical(fit@covariates, c("x1", "x2"))
})

test_that("balance() accepts tidyselect helpers for covariates", {
  data <- sim_binary(n = 200)
  fit <- balance(
    data,
    exposure,
    tidyselect::starts_with("x"),
    method = bw_entropy()
  )

  expect_identical(fit@covariates, c("x1", "x2", "x3"))
})

test_that("balance() errors when the exposure selects more than one column", {
  data <- sim_binary(n = 200)
  expect_error(
    balance(data, c(x1, x2), c(x1, x2), method = bw_entropy()),
    class = "balancing_selection_error"
  )
})

test_that("balance() errors when the covariate selection is empty", {
  data <- sim_binary(n = 200)
  expect_error(
    balance(
      data,
      exposure,
      tidyselect::starts_with("nonexistent"),
      method = bw_entropy()
    ),
    class = "balancing_selection_error"
  )
})

# ---- Exposure-type detection ----------------------------------------------

test_that("balance() auto-detects the exposure type and stores it", {
  data <- sim_binary(n = 200)
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())
  expect_identical(fit@exposure_type, "binary")

  fit_cat <- balance(
    sim_categorical(n = 200),
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_identical(fit_cat@exposure_type, "categorical")

  fit_cont <- balance(
    sim_continuous(n = 200),
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_identical(fit_cont@exposure_type, "continuous")
})

test_that("balance() announces the detected exposure type", {
  withr::local_options(balancing.quiet = FALSE)
  data <- sim_binary(n = 200)
  expect_message(
    balance(data, exposure, c(x1, x2), method = bw_entropy()),
    regexp = "binary"
  )
})

test_that("balance() stays silent when the quiet option is set", {
  withr::local_options(balancing.quiet = TRUE)
  data <- sim_binary(n = 200)
  expect_silent(
    balance(data, exposure, c(x1, x2), method = bw_entropy())
  )
})

test_that("a forced exposure type that contradicts the data errors", {
  data <- sim_continuous(n = 200)
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      exposure_type = "binary"
    ),
    class = "balancing_exposure_type_error"
  )
})

# ---- Unused exposure levels -----------------------------------------------

test_that("a binary factor with an unused level fits an estimating-equation method", {
  withr::local_options(balancing.quiet = TRUE)
  data <- sim_binary(n = 200)
  data$exposure <- factor(data$exposure, levels = c(0, 1, 2))
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())
  expect_equal(names(attr(fit@weights, "groups")), c("0", "1"))
  expect_true(all(is.finite(as.numeric(stats::weights(fit)))))
})

test_that("a binary factor with an unused level fits a quadratic-program method", {
  withr::local_options(balancing.quiet = TRUE)
  data <- sim_binary(n = 200)
  data$exposure <- factor(data$exposure, levels = c(0, 1, 2))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    constraints = balance_terms(tolerance = 0.05)
  )
  groups <- attr(fit@weights, "groups")
  w <- as.numeric(weights(fit))
  expect_equal(names(groups), c("0", "1"))
  group_ess <- vapply(
    groups,
    function(idx) sum(w[idx])^2 / sum(w[idx]^2),
    numeric(1)
  )
  expect_true(all(is.finite(group_ess)))
})

test_that("a categorical factor with an unused level fits an estimating-equation method", {
  withr::local_options(balancing.quiet = TRUE)
  data <- sim_categorical(n = 300)
  present <- levels(factor(as.character(data$exposure)))
  data$exposure <- factor(
    as.character(data$exposure),
    levels = c(present, "zzz")
  )
  fit <- balance(data, exposure, c(x1, x2), method = bw_ipt())
  expect_false("zzz" %in% names(attr(fit@weights, "groups")))
  expect_true(all(is.finite(as.numeric(stats::weights(fit)))))
})

test_that("a categorical factor with an unused level keeps a finite effective sample size", {
  withr::local_options(balancing.quiet = TRUE)
  data <- sim_categorical(n = 300)
  present <- levels(factor(as.character(data$exposure)))
  data$exposure <- factor(
    as.character(data$exposure),
    levels = c(present, "zzz")
  )
  fit <- balance(data, exposure, c(x1, x2), method = bw_energy())
  groups <- attr(fit@weights, "groups")
  w <- as.numeric(weights(fit))
  expect_false("zzz" %in% names(groups))
  group_ess <- vapply(
    groups,
    function(idx) sum(w[idx])^2 / sum(w[idx]^2),
    numeric(1)
  )
  expect_true(all(is.finite(group_ess)))
})

test_that("dropping an unused exposure level announces itself", {
  withr::local_options(balancing.quiet = FALSE)
  data <- sim_binary(n = 200)
  data$exposure <- factor(data$exposure, levels = c(0, 1, 2))
  expect_snapshot(
    fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())
  )
})

# ---- Estimand validation and synonyms -------------------------------------

test_that("balance() stores the requested estimand for a binary exposure", {
  data <- sim_binary(n = 200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att"
  )
  expect_identical(fit@estimand, "att")
})

test_that("the untreated-target synonyms resolve to propensity's canonical", {
  data <- sim_binary(n = 200)
  fit_atc <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "atc"
  )
  fit_atu <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "atu"
  )

  # propensity stores "atu" as the canonical string for the untreated target.
  expect_identical(fit_atc@estimand, "atu")
  expect_identical(fit_atu@estimand, "atu")
})

test_that("an estimand unsupported by the method and exposure type errors", {
  data <- sim_binary(n = 200)
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      estimand = "ato"
    ),
    class = "balancing_estimand_error"
  )
})

test_that("a continuous exposure permits only the ate estimand", {
  data <- sim_continuous(n = 200)
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      estimand = "att"
    ),
    class = "balancing_estimand_error"
  )
})

# ---- Focal level ----------------------------------------------------------

test_that("a binary att infers the treated level without focal_level", {
  data <- sim_binary(n = 200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att"
  )
  expect_identical(fit@focal_level, "1")
})

test_that("a categorical att requires focal_level", {
  data <- sim_categorical(n = 200)
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      estimand = "att"
    ),
    class = "balancing_estimand_error"
  )
})

test_that("a categorical att honors a supplied focal_level", {
  data <- sim_categorical(n = 200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att",
    focal_level = "b"
  )
  expect_identical(fit@focal_level, "b")
})

# ---- Sampling weights -----------------------------------------------------

test_that("balance() evaluates sampling_weights given as a bare column", {
  data <- sim_binary(n = 200)
  data$sw <- stats::runif(nrow(data), 0.5, 2)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    sampling_weights = sw
  )
  expect_length(fit@sampling_weights, nrow(data))
  expect_equal(fit@sampling_weights, data$sw)
})

test_that("balance() evaluates sampling_weights given as an external vector", {
  data <- sim_binary(n = 200)
  external <- stats::runif(nrow(data), 0.5, 2)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    sampling_weights = external
  )
  expect_length(fit@sampling_weights, nrow(data))
  expect_equal(fit@sampling_weights, external)
})

test_that("negative sampling weights error", {
  data <- sim_binary(n = 200)
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      sampling_weights = rep(-1, nrow(data))
    ),
    class = "balancing_error"
  )
})

# ---- Dots and missing values ----------------------------------------------

test_that("balance() rejects unnamed arguments through check_dots_empty()", {
  data <- sim_binary(n = 200)
  # A valid call succeeds; an extra unnamed argument trips check_dots_empty().
  expect_true(S7::S7_inherits(
    balance(data, exposure, c(x1, x2), method = bw_entropy()),
    balancing
  ))
  expect_error(
    balance(data, exposure, c(x1, x2), method = bw_entropy(), 5)
  )
})

test_that("missing values in the covariates error", {
  data <- sim_binary(n = 200)
  data$x1[1] <- NA
  expect_error(
    balance(data, exposure, c(x1, x2), method = bw_entropy()),
    class = "balancing_missing_error"
  )
})

test_that("missing values in the exposure error", {
  data <- sim_binary(n = 200)
  data$exposure[1] <- NA
  expect_error(
    balance(data, exposure, c(x1, x2), method = bw_entropy()),
    class = "balancing_missing_error"
  )
})

# ---- Constant covariates --------------------------------------------------

test_that("a constant covariate leaves a continuous-exposure fit intact", {
  data <- sim_continuous(n = 200)
  data$fixed <- 5
  fit <- balance(data, exposure, c(x1, x2, fixed), method = bw_entropy())

  # The constant column is dropped before the fit, so the balance table reports
  # only the covariates that carry information and every verdict resolves.
  expect_false("fixed" %in% fit@balance_table$term)
  expect_true(all(is.finite(fit@balance_table$weighted)))
  expect_true(all(fit@balance_table$within_tolerance))
})

test_that("a single-level factor leaves a continuous-exposure fit intact", {
  data <- sim_continuous(n = 200)
  data$f <- factor(rep("a", nrow(data)))
  fit <- balance(data, exposure, c(x1, x2, f), method = bw_entropy())

  expect_false("f_a" %in% fit@balance_table$term)
  expect_true(all(is.finite(fit@balance_table$weighted)))
  expect_true(all(fit@balance_table$within_tolerance))
})

test_that("a single-level factor leaves a binary-exposure fit intact", {
  data <- sim_binary(n = 200)
  data$f <- factor(rep("a", nrow(data)))
  fit <- balance(data, exposure, c(x1, x2, f), method = bw_entropy())

  expect_false("f_a" %in% fit@balance_table$term)
  expect_true(all(fit@balance_table$within_tolerance))
})

test_that("a constant covariate leaves a binary-exposure fit and ipw() intact", {
  data <- sim_binary(n = 200)
  data$fixed <- 5
  data$y <- withr::with_seed(909, stats::rnorm(nrow(data)))
  fit <- balance(data, exposure, c(x1, x2, fixed), method = bw_entropy())

  expect_false("fixed" %in% fit@balance_table$term)
  expect_true(all(fit@balance_table$within_tolerance))

  data$.wts <- as.numeric(stats::weights(fit))
  outcome_model <- suppressWarnings(stats::glm(
    y ~ exposure,
    data = data,
    weights = .wts
  ))
  result <- ipw(fit, outcome_model)
  expect_true(all(is.finite(result$estimates$estimate)))
  expect_true(all(is.finite(result$estimates$std.err)))
})

test_that("a covariate set with no spread at all is a classed error", {
  data <- sim_continuous(n = 200)
  data$fixed <- 5
  expect_error(
    balance(data, exposure, fixed, method = bw_entropy()),
    class = "balancing_constraints_error"
  )
})

test_that("an unresolved tolerance verdict warns rather than aborting", {
  # Defense in depth for the balance warning: whatever the constraint set, a
  # verdict that does not resolve to TRUE is reported as out of tolerance rather
  # than steering an `if` with a missing value.
  data <- sim_binary(n = 200)
  original <- compute_balance_table
  testthat::local_mocked_bindings(
    compute_balance_table = function(...) {
      table <- original(...)
      table$within_tolerance[[1]] <- NA
      table
    }
  )
  expect_warning(
    balance(data, exposure, c(x1, x2), method = bw_entropy()),
    class = "balancing_balance_warning"
  )
})

# ---- Method argument ------------------------------------------------------

test_that("a bare-string method errors", {
  data <- sim_binary(n = 200)
  expect_error(
    balance(data, exposure, c(x1, x2), method = "entropy"),
    class = "balancing_method_error"
  )
})

test_that("balance() requires a data frame", {
  expect_error(
    balance(list(a = 1), exposure, c(x1, x2), method = bw_entropy()),
    class = "balancing_type_error"
  )
})
