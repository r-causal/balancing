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

test_that("everything() excludes the exposure from its own covariates", {
  # A selection that resolves against the whole data frame reaches the exposure
  # too. Balancing the exposure against itself is infeasible by construction, so
  # the exposure column leaves the covariate selection.
  data <- sim_binary(n = 200)
  fit <- balance(data, exposure, everything(), method = bw_entropy())

  expect_identical(fit@covariates, c("x1", "x2", "x3"))
  expect_false("exposure" %in% fit@balance_table$term)
  expect_true(all(is.finite(as.numeric(stats::weights(fit)))))
})

test_that("a covariate selection naming only the exposure is a classed error", {
  data <- sim_binary(n = 200)
  expect_error(
    balance(data, exposure, exposure, method = bw_entropy()),
    class = "balancing_selection_error"
  )
})

# A rename inside a selection resolves the column by position but supplies a new
# name, and every downstream lookup indexes the data by name. A rename that
# collided with the exposure therefore built the constraints on the exposure
# column and died blaming collinearity; a partial rename died on an unclassed
# dimnames error; a renamed exposure produced a fit naming a column the data do
# not have. Renaming is refused at the selection instead.
test_that("a renamed covariate selection is refused", {
  data <- sim_binary(n = 200)
  expect_error(
    balance(data, exposure, c(exposure = x1), method = bw_entropy()),
    "rename"
  )
  expect_error(
    balance(data, exposure, c(foo = x1, x2), method = bw_entropy()),
    "rename"
  )
})

test_that("a renamed exposure selection is refused", {
  data <- sim_binary(n = 200)
  expect_error(
    balance(data, c(foo = exposure), c(x1, x2), method = bw_entropy()),
    "rename"
  )
})

test_that("an unrenamed selection still resolves", {
  data <- sim_binary(n = 200)
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())
  expect_identical(fit@exposure, "exposure")
  expect_identical(fit@covariates, c("x1", "x2"))
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
    class = "causalgenerics_forced_exposure_type"
  )
})

test_that("an explicit continuous type fits a low-cardinality dose", {
  # Ten distinct doses in 200 observations sit under the categorical heuristic's
  # unique-share threshold, so only the explicit declaration carries the fit to
  # the continuous path, where balance is a weighted correlation.
  data <- withr::with_seed(9, {
    n <- 200
    data.frame(
      dose = sample(seq(10, 100, by = 10), n, replace = TRUE),
      x1 = stats::rnorm(n),
      x2 = stats::rnorm(n)
    )
  })
  fit <- balance(
    data,
    dose,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate",
    exposure_type = "continuous"
  )

  expect_identical(fit@exposure_type, "continuous")
  expect_identical(fit@balance_table$statistic, c("correlation", "correlation"))
  expect_balanced(fit, data)
})

# ---- Unused exposure levels -----------------------------------------------

test_that("a binary factor with an unused level fits an estimating-equation method", {
  withr::local_options(balancing.quiet = TRUE)
  data <- sim_binary(n = 200)
  data$exposure <- factor(data$exposure, levels = c(0, 1, 2))
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())
  expect_equal(fit@exposure_levels, c("0", "1"))
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
  groups <- split(seq_len(nrow(data)), as.character(data$exposure))
  w <- as.numeric(weights(fit))
  expect_equal(fit@exposure_levels, c("0", "1"))
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
  expect_false("zzz" %in% fit@exposure_levels)
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
  groups <- split(seq_len(nrow(data)), as.character(data$exposure))
  w <- as.numeric(weights(fit))
  expect_false("zzz" %in% fit@exposure_levels)
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

# ---- Exposure level order -------------------------------------------------

test_that("numeric exposure levels order numerically, matching factor()", {
  # A character sort puts "10" before "9", which would name the wrong second
  # level. The order follows the data's own type, as base factor() does.
  exposure <- c(9, 10, 9, 10)
  expect_identical(exposure_levels(exposure, "binary"), c("9", "10"))
  expect_identical(
    exposure_levels(exposure, "binary"),
    levels(factor(exposure))
  )
})

test_that("a many-level numeric exposure orders numerically", {
  exposure <- c(2, 10, 1, 20)
  expect_identical(
    exposure_levels(exposure, "categorical"),
    c("1", "2", "10", "20")
  )
  expect_identical(
    exposure_levels(exposure, "categorical"),
    levels(factor(exposure))
  )
})

test_that("character exposure levels keep the character sort", {
  exposure <- c("9", "10", "9", "10")
  expect_identical(exposure_levels(exposure, "binary"), c("10", "9"))
  expect_identical(
    exposure_levels(exposure, "binary"),
    levels(factor(exposure))
  )
})

test_that("factor exposure levels keep the declared order", {
  exposure <- factor(c("hi", "lo", "hi"), levels = c("lo", "hi"))
  expect_identical(exposure_levels(exposure, "binary"), c("lo", "hi"))
})

# ---- The fit's exposure-level property -------------------------------------

# The levels a fit weighted are part of its contract: their first element is the
# reference level every contrast in `ipw()` is measured against. They are recorded
# as a property whose order is `levels(factor(x))`, so a consumer reads them
# without knowing how the fit stored its groups. They used to ride as an
# undocumented attribute on the weight vector.
test_that("the fit records its exposure levels as a property", {
  data <- sim_binary(n = 200)
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())
  expect_identical(fit@exposure_levels, c("0", "1"))
  expect_identical(fit@exposure_levels, levels(factor(data$exposure)))
  expect_null(attr(fit@weights, "groups"))
})

test_that("the recorded exposure levels order numerically, not as characters", {
  data <- sim_binary(n = 200)
  data$exposure <- ifelse(data$exposure == 1, 10, 9)
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())
  expect_identical(fit@exposure_levels, c("9", "10"))
})

test_that("the recorded exposure levels keep a factor's declared order", {
  data <- sim_binary(n = 200)
  data$exposure <- factor(
    ifelse(data$exposure == 1, "hi", "lo"),
    levels = c("hi", "lo")
  )
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())
  expect_identical(fit@exposure_levels, c("hi", "lo"))
})

test_that("the recorded exposure levels list every observed level", {
  data <- sim_categorical(n = 300)
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())
  expect_identical(fit@exposure_levels, levels(factor(data$exposure)))
  expect_gt(length(fit@exposure_levels), 2L)
})

test_that("the recorded exposure levels omit an unused factor level", {
  withr::local_options(balancing.quiet = TRUE)
  data <- sim_binary(n = 200)
  data$exposure <- factor(data$exposure, levels = c(0, 1, 2))
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())
  expect_identical(fit@exposure_levels, c("0", "1"))
})

test_that("a continuous fit records no exposure levels", {
  data <- sim_continuous(n = 150)
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())
  expect_identical(fit@exposure_levels, character(0))
})

test_that("every method family records the same exposure levels", {
  data <- sim_binary(n = 200)
  methods <- list(bw_entropy(), bw_ipt(), bw_cbps(), bw_energy(), bw_cfd())
  for (method in methods) {
    fit <- balance(data, exposure, c(x1, x2), method = method)
    expect_identical(fit@exposure_levels, c("0", "1"))
  }
  # Stable balancing weights need a positive tolerance, their central knob.
  sbw <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    constraints = balance_terms(tolerance = 0.05)
  )
  expect_identical(sbw@exposure_levels, c("0", "1"))
})

# ---- Focal level ----------------------------------------------------------

test_that("a binary att targets the numerically larger exposure level", {
  # With levels 9 and 10 a character sort would infer 9 as the treated level and
  # reweight the wrong group. The treated level is the second level in the data's
  # own order, so it is 10.
  data <- withr::with_seed(11, {
    n <- 300
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    z <- stats::rbinom(n, 1L, stats::plogis(0.6 * x1 - 0.4 * x2))
    data.frame(exposure = ifelse(z == 1L, 10, 9), x1 = x1, x2 = x2)
  })

  att <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att"
  )
  expect_identical(att@focal_level, "10")
  expect_identical(att@exposure_levels, c("9", "10"))

  # Entropy balancing leaves the focal group at its uniform base weights carried
  # to its own total, so a fit that held the wrong group fixed would show level
  # 10 reweighted instead of level 9.
  w <- as.numeric(stats::weights(att))
  focal <- data$exposure == 10
  expect_equal(w[focal], rep(1, sum(focal)))
  expect_false(isTRUE(all.equal(w[!focal], rep(1, sum(!focal)))))

  atc <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "atc"
  )
  expect_identical(atc@focal_level, "9")
})

test_that("a focal estimand needs a level outside the focal group", {
  # A single-level exposure leaves a focal estimand with no group to reweight, so
  # the solve would carry no parameters at all.
  data <- withr::with_seed(5, {
    n <- 60
    data.frame(
      exposure = rep(1L, n),
      x1 = stats::rnorm(n),
      x2 = stats::rnorm(n)
    )
  })
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      estimand = "att",
      focal_level = 1
    ),
    class = "balancing_estimand_error"
  )
})

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

# ---- A focal level the estimand ignores ------------------------------------

# The average treatment effect and the overlap estimand reweight every exposure
# group rather than holding one fixed, so they resolve no focal level and never
# reach the check that the supplied one is an exposure level at all. A level that
# does not exist used to be accepted in silence, which reads as a fit that
# targeted it.
test_that("focal_level with the average treatment effect warns and is ignored", {
  data <- sim_binary(n = 200)
  expect_warning(
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      estimand = "ate",
      focal_level = 1
    ),
    class = "balancing_ignored_argument_warning"
  )
  expect_null(fit@focal_level)
})

test_that("a focal_level that is not an exposure level still warns", {
  data <- sim_binary(n = 200)
  expect_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      estimand = "ate",
      focal_level = "nonesuch"
    ),
    class = "balancing_ignored_argument_warning"
  )
})

test_that("focal_level with the overlap estimand warns and is ignored", {
  data <- sim_binary(n = 200)
  expect_warning(
    fit <- balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(),
      estimand = "ato",
      focal_level = 0
    ),
    class = "balancing_ignored_argument_warning"
  )
  expect_null(fit@focal_level)
})

test_that("a pooled estimand without focal_level is silent", {
  data <- sim_binary(n = 200)
  expect_no_warning(
    balance(data, exposure, c(x1, x2), method = bw_entropy(), estimand = "ate")
  )
  expect_no_warning(
    balance(data, exposure, c(x1, x2), method = bw_cbps(), estimand = "ato")
  )
})

test_that("a focal estimand with focal_level does not warn", {
  data <- sim_categorical(n = 200)
  expect_no_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      estimand = "att",
      focal_level = "b"
    )
  )
})

# ---- An exposure with a single level ---------------------------------------

# Balancing reweights one exposure group toward another, so an exposure that takes
# a single level leaves nothing to balance whatever the estimand. The focal
# estimands refused it already; the pooled estimands fitted the uniform weighting
# they started from and then stopped in the balance table on a maximum over no
# contrasts, an unclassed base error.
single_level_data <- function() {
  withr::with_seed(21, {
    n <- 60
    data.frame(
      exposure = rep(1L, n),
      x1 = stats::rnorm(n),
      x2 = stats::rnorm(n)
    )
  })
}

test_that("a single-level exposure is refused for every estimand", {
  withr::local_options(balancing.quiet = TRUE)
  data <- single_level_data()
  expect_error(
    balance(data, exposure, c(x1, x2), method = bw_entropy(), estimand = "ate"),
    class = "balancing_estimand_error"
  )
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      estimand = "att",
      focal_level = 1
    ),
    class = "balancing_estimand_error"
  )
  expect_error(
    balance(data, exposure, c(x1, x2), method = bw_cbps(), estimand = "ato"),
    class = "balancing_estimand_error"
  )
})

test_that("a single-level exposure is refused across the method families", {
  withr::local_options(balancing.quiet = TRUE)
  data <- single_level_data()
  methods <- list(
    bw_entropy(),
    bw_ipt(),
    bw_cbps(),
    bw_energy(),
    bw_cfd(),
    bw_sbw()
  )
  for (method in methods) {
    expect_error(
      balance(data, exposure, c(x1, x2), method = method, estimand = "ate"),
      class = "balancing_estimand_error"
    )
  }
})

test_that("a single-level factor exposure is refused", {
  withr::local_options(balancing.quiet = TRUE)
  data <- single_level_data()
  data$exposure <- factor(rep("a", nrow(data)))
  expect_error(
    balance(data, exposure, c(x1, x2), method = bw_entropy()),
    class = "balancing_estimand_error"
  )
})

test_that("a declared level no observation takes leaves a single-level exposure", {
  # The unused level is dropped before the count, so a two-level factor with one
  # level unobserved is a single-level exposure and refused as one.
  withr::local_options(balancing.quiet = TRUE)
  data <- single_level_data()
  data$exposure <- factor(rep("a", nrow(data)), levels = c("a", "b"))
  expect_error(
    balance(data, exposure, c(x1, x2), method = bw_entropy()),
    class = "balancing_estimand_error"
  )
})

test_that("a continuous exposure is not measured by the level count", {
  # A continuous exposure carries no levels at all, so the rule does not reach it.
  data <- sim_continuous(n = 150)
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())
  expect_identical(fit@exposure_levels, character(0))
})

# ---- Sampling weights -----------------------------------------------------

test_that("balance() evaluates sampling_weights given as a bare column", {
  data <- sim_binary(n = 200)
  data$sw <- withr::with_seed(414, stats::runif(nrow(data), 0.5, 2))
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
  external <- withr::with_seed(414, stats::runif(nrow(data), 0.5, 2))
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

test_that("infinite sampling weights error", {
  data <- sim_binary(n = 200)
  weights <- rep(1, nrow(data))
  weights[1] <- Inf
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      sampling_weights = weights
    ),
    class = "balancing_range_error"
  )
  # Both arms are refused by `validate_sampling_weights()`, before any method
  # sees the weights, so the second arm pins that the refusal is the fit's own
  # rather than something a particular solver happens to catch.
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(),
      sampling_weights = weights
    ),
    class = "balancing_range_error"
  )
})

test_that("all-zero sampling weights error", {
  data <- sim_binary(n = 200)
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      sampling_weights = rep(0, nrow(data))
    ),
    class = "balancing_range_error"
  )
})

# ---- Zero-measure groups --------------------------------------------------

# A group with no base-measure mass has no weighted mean to target and no total to
# report at, so every quantity built from it is undefined. Depending on the method
# that used to surface as a fit whose balance table was entirely missing values,
# as an infeasibility blamed on the constraint set, or as a convergence error
# blamed on collinearity. The requirement is checked once, before the fit, so
# every method reports the same defect.
test_that("a focal group with no sampling-weight mass is a classed error", {
  data <- sim_binary(n = 200)
  focal_zero <- ifelse(data$exposure == 1L, 0, 1)
  for (method in list(bw_entropy(), bw_ipt(), bw_cbps(), bw_energy())) {
    expect_error(
      balance(
        data,
        exposure,
        c(x1, x2),
        method = method,
        estimand = "att",
        sampling_weights = focal_zero
      ),
      class = "balancing_range_error"
    )
  }
})

test_that("any group with no sampling-weight mass is a classed error", {
  data <- sim_binary(n = 200)
  group_zero <- ifelse(data$exposure == 1L, 0, 1)
  for (method in list(
    bw_entropy(),
    bw_ipt(),
    bw_cbps(),
    bw_energy(),
    bw_cfd()
  )) {
    expect_error(
      balance(
        data,
        exposure,
        c(x1, x2),
        method = method,
        sampling_weights = group_zero
      ),
      class = "balancing_range_error"
    )
  }
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_sbw(),
      sampling_weights = group_zero,
      constraints = balance_terms(tolerance = 0.05)
    ),
    class = "balancing_range_error"
  )
})

test_that("a categorical level with no sampling-weight mass is a classed error", {
  data <- sim_categorical(n = 200)
  level_zero <- ifelse(data$exposure == "a", 0, 1)
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      sampling_weights = level_zero
    ),
    class = "balancing_range_error"
  )
})

# The sampling weights and the base weights multiply, so each can be nonzero
# somewhere while their product is zero everywhere. That left every constraint
# target a ratio of zero totals and died on the group renormalization blaming
# collinearity.
test_that("disjoint sampling and base weight supports are a classed error", {
  data <- sim_binary(n = 200)
  n <- nrow(data)
  sampling <- rep(c(0, 1), length.out = n)
  base <- rep(c(1, 0), length.out = n)
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(base_weights = base),
      sampling_weights = sampling
    ),
    class = "balancing_range_error"
  )
})

test_that("disjoint supports error on the continuous path too", {
  data <- sim_continuous(n = 200)
  n <- nrow(data)
  sampling <- rep(c(0, 1), length.out = n)
  base <- rep(c(1, 0), length.out = n)
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(base_weights = base),
      sampling_weights = sampling
    ),
    class = "balancing_range_error"
  )
})

# Individual zero weights stay legal, and a group that keeps mass under a partly
# zero vector still fits, so the check refuses only a group with nothing left.
test_that("partly zero sampling weights that leave every group mass still fit", {
  data <- sim_binary(n = 200)
  sampling <- rep(c(0, 1), length.out = nrow(data))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    sampling_weights = sampling
  )
  expect_true(all(is.finite(as.numeric(stats::weights(fit)))))
  expect_true(all(is.finite(fit@balance_table$weighted)))
})

# The base-weight length mismatch belongs to the entropy fit, which names the
# argument and the sample size, so the measure check must not intercept it.
test_that("a base weight vector of the wrong length still reports its length", {
  data <- sim_binary(n = 200)
  expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(base_weights = rep(1, 3))
    ),
    "one value per observation"
  )
})

# Sampling weights are naturally supplied as integer counts, and both the solver
# boundary and the fitted object's property take a double, so an integer vector
# fits and is stored coerced.
test_that("integer sampling weights fit and are stored as doubles", {
  data <- sim_binary(n = 200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    sampling_weights = rep(2L, nrow(data))
  )
  expect_type(fit@sampling_weights, "double")
  expect_equal(fit@sampling_weights, rep(2, nrow(data)))
  doubled <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    sampling_weights = rep(2, nrow(data))
  )
  expect_equal(
    as.numeric(stats::weights(fit)),
    as.numeric(stats::weights(doubled))
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

# An infinity survives `anyNA()` and then standardizes to another infinity, so it
# is refused alongside the missing values rather than left to poison the solve.
test_that("infinite values in the covariates error", {
  data <- sim_binary(n = 200)
  data$x1[1] <- Inf
  expect_error(
    balance(data, exposure, c(x1, x2), method = bw_entropy()),
    class = "balancing_range_error"
  )
  data$x1[1] <- -Inf
  expect_error(
    balance(data, exposure, c(x1, x2), method = bw_entropy()),
    class = "balancing_range_error"
  )
})

test_that("infinite values in a continuous exposure error", {
  data <- sim_continuous(n = 200)
  data$exposure[1] <- Inf
  expect_error(
    balance(data, exposure, c(x1, x2), method = bw_entropy()),
    class = "balancing_range_error"
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

# ---- What @covariates records ----------------------------------------------

# The slot names the covariates the fit constrained, not the ones the caller
# selected. A constant or aliased column is dropped at expansion and a `moments`
# request of zero contributes none, so a selected covariate can end up with no
# constraint column at all; listing it anyway claimed balance the fit never
# targeted. The request itself stays visible in the recorded call.
test_that("@covariates omits a covariate whose column was dropped as constant", {
  withr::local_options(balancing.quiet = TRUE)
  data <- sim_binary(n = 200)
  data$fixed <- 5
  fit <- balance(data, exposure, c(x1, x2, fixed), method = bw_entropy())
  expect_identical(fit@covariates, c("x1", "x2"))
  expect_match(paste(deparse(fit@call), collapse = " "), "fixed", fixed = TRUE)
})

test_that("@covariates omits a covariate whose column was dropped as aliased", {
  withr::local_options(balancing.quiet = TRUE)
  data <- sim_binary(n = 200)
  data$copy <- data$x1
  fit <- balance(data, exposure, c(x1, x2, copy), method = bw_entropy())
  expect_identical(fit@covariates, c("x1", "x2"))
})

test_that("@covariates omits a covariate with no moments requested", {
  data <- sim_binary(n = 200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    constraints = balance_terms(moments = c(x1 = 0L, x2 = 1L))
  )
  expect_identical(fit@covariates, "x2")
})

test_that("@covariates keeps a covariate that only partners an interaction", {
  # An interaction column constrains both of its factors, so a covariate with no
  # column of its own is still constrained through the product.
  data <- sim_binary(n = 200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    constraints = balance_terms(
      moments = c(x1 = 1L, x2 = 0L),
      interactions = TRUE
    )
  )
  terms <- fit@balance_table$term
  expect_true(any(grepl("x2", terms, fixed = TRUE)))
  expect_identical(fit@covariates, c("x1", "x2"))
})

test_that("@covariates lists the whole selection when every column is kept", {
  data <- sim_binary(n = 200)
  fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())
  expect_identical(fit@covariates, c("x1", "x2"))
})

test_that("an objective-driven fit with no constraints records no covariates", {
  # Energy and kernel balancing are driven by their objective, which reads every
  # selected covariate, so a fit with no constraint columns balances them without
  # constraining any. The slot records what was constrained, which is nothing.
  data <- sim_binary(n = 200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    constraints = balance_terms(moments = 0L)
  )
  expect_identical(fit@covariates, character(0))
})

test_that("a covariate set with no spread at all is a classed error", {
  data <- sim_continuous(n = 200)
  data$fixed <- 5
  expect_error(
    balance(data, exposure, fixed, method = bw_entropy()),
    class = "balancing_constraints_error"
  )
})

# A balance statistic that is not a number says nothing about how far a fit
# missed, so the warning must not offer it as the largest imbalance. It reports
# that balance could not be assessed instead. No fit reaches this once every group
# is required to carry measure, so the balance table is mocked to reach it.
test_that("an undefined imbalance warns that balance could not be assessed", {
  data <- sim_binary(n = 200)
  original <- compute_balance_table
  testthat::local_mocked_bindings(
    compute_balance_table = function(...) {
      table <- original(...)
      table$weighted[[1]] <- NaN
      table$within_tolerance[[1]] <- FALSE
      table
    }
  )
  expect_warning(
    balance(data, exposure, c(x1, x2), method = bw_entropy()),
    "could not be assessed",
    class = "balancing_balance_warning"
  )
})

test_that("a missing imbalance warns that balance could not be assessed", {
  # A missing statistic makes the maximum NA rather than NaN, which states a
  # distance no better than NaN does, so the same could-not-assess warning covers
  # it. This is the shape an unresolved verdict arrives in alongside a statistic
  # that was never computed.
  data <- sim_binary(n = 200)
  original <- compute_balance_table
  testthat::local_mocked_bindings(
    compute_balance_table = function(...) {
      table <- original(...)
      table$weighted[[1]] <- NA_real_
      table$within_tolerance[[1]] <- NA
      table
    }
  )
  expect_warning(
    balance(data, exposure, c(x1, x2), method = bw_entropy()),
    "could not be assessed",
    class = "balancing_balance_warning"
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

# ---- Solver status routing --------------------------------------------------

# A quadratic-program backend can stop for reasons no iteration cap explains: it
# can break down numerically, stall, or return without having solved at all.
# Advice to raise `max_iterations` is useless in every one of those cases, so each
# raises the convergence error naming the conditioning of the problem instead.
# The statuses come from the compiled backend, so they are exercised through the
# routing directly rather than by provoking a numerical breakdown.
test_that("a solver breakdown reports conditioning rather than the iteration cap", {
  for (status in c("numerical_error", "insufficient_progress", "not_solved")) {
    expect_error(
      check_solver_status(
        list(status = status, converged = FALSE),
        bw_sbw()
      ),
      "stopped without a solution",
      class = "balancing_convergence_error"
    )
  }
})

# The conditioning advice names a knob, so it must name one the method actually
# has. `weight_penalty` is the ridge term the quadratic-program objective carries;
# `bw_energy()` and `bw_cfd()` expose it as an argument, while `bw_sbw()` holds it
# at zero by design and takes no such argument, so advising a caller to raise it
# in `bw_sbw()` sends them after an argument that does not exist. The advice is
# therefore conditional on the method rather than fixed text.
test_that("the conditioning advice names weight_penalty only where there is one", {
  for (status in c("non_convex", "numerical_error")) {
    sbw <- expect_error(
      check_solver_status(
        list(status = status, converged = FALSE),
        bw_sbw()
      ),
      class = "balancing_convergence_error"
    )
    expect_false(grepl("weight_penalty", conditionMessage(sbw), fixed = TRUE))

    energy <- expect_error(
      check_solver_status(
        list(status = status, converged = FALSE),
        bw_energy()
      ),
      class = "balancing_convergence_error"
    )
    expect_true(grepl("weight_penalty", conditionMessage(energy), fixed = TRUE))
  }
})

test_that("a reached iteration cap still advises raising it", {
  expect_warning(
    check_solver_status(
      list(status = "max_iter", converged = FALSE),
      bw_sbw()
    ),
    "max_iterations",
    class = "balancing_convergence_warning"
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

# ---- Sample size ----------------------------------------------------------

# A single observation is not enough to balance. A numeric covariate crosses the
# boundary standardized by its own spread, and the standard deviation of one
# value is a missing value, so the branch that rescues a column with no spread
# used to steer on that missing value and the fit died with the base error
# "missing value where TRUE/FALSE needed". Nothing earlier turns the sample away:
# a lone exposure value is one unique value among one observation, which the
# unique-value heuristic reads as continuous, so the rule requiring two exposure
# levels never measures it. The requirement belongs to the data rather than to a
# method, so it is refused once, before any column is standardized, and every
# method reports the same defect.
test_that("a one-row data frame is a classed error", {
  data <- sim_binary(n = 200)[1, ]

  energy <- expect_error(
    balance(data, exposure, c(x1, x2), method = bw_energy(), estimand = "ate"),
    class = "balancing_empty_error"
  )
  sbw <- expect_error(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_sbw(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0.05)
    ),
    class = "balancing_empty_error"
  )

  # The wording is free to change, but the refusal has to name a size: either the
  # observations that arrived or the minimum they fell short of.
  for (cnd in list(energy, sbw)) {
    expect_match(conditionMessage(cnd), "\\b(1|one|2|two)\\b")
  }
})
