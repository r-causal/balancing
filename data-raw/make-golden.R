# Generate the Rust-side golden fixtures for entropy balancing.
#
# Each fixture captures the entropy solver's inputs at the matrix level (the
# standardized constraint matrix, the targets, the tolerances, the base and
# sampling weights, and the per-group normalization) together with the reference
# weights from a pinned WeightIt. The Rust harness at
# crates/core/tests/golden.rs feeds the captured inputs to solve_discrete() and
# solve_continuous() and checks that our weights reproduce the reference. This
# isolates solver parity from covariate processing, which the R test suite covers
# separately.
#
# This script is the only place outside the tests that names WeightIt. It is
# Rbuildignored and committed. Regenerate the fixtures deliberately with a pinned
# WeightIt and review the resulting JSON diff.
#
# Run from the package root:
#   Rscript data-raw/make-golden.R

library(WeightIt)

# Pin the reference implementation. Bump deliberately and review the fixture
# diff when updating.
weightit_pin <- "1.7.0"
installed <- as.character(utils::packageVersion("WeightIt"))
if (installed != weightit_pin) {
  message(sprintf(
    "WeightIt %s is installed; the fixtures are pinned to %s. Install the pinned version before regenerating.",
    installed,
    weightit_pin
  ))
}

output_dir <- file.path("src", "rust", "crates", "core", "tests", "golden")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Seeded data-generating processes, matching the shapes the test helpers use so
# the fixtures exercise realistic confounding.
sim_binary <- function(n, seed) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  linear_predictor <- 0.8 * x1 - 0.6 * x2
  exposure <- rbinom(n, 1L, plogis(linear_predictor))
  data.frame(exposure = exposure, x1 = x1, x2 = x2)
}

sim_categorical <- function(n, seed) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  eta_b <- 0.7 * x1 - 0.3 * x2
  eta_c <- -0.5 * x1 + 0.6 * x2
  denom <- 1 + exp(eta_b) + exp(eta_c)
  draw <- runif(n)
  level <- ifelse(
    draw < 1 / denom,
    "a",
    ifelse(draw < (1 + exp(eta_b)) / denom, "b", "c")
  )
  data.frame(exposure = factor(level, levels = c("a", "b", "c")), x1, x2)
}

sim_continuous <- function(n, seed) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  exposure <- 0.9 * x1 - 0.6 * x2 + rnorm(n)
  data.frame(exposure = exposure, x1 = x1, x2 = x2)
}

# Standardize a covariate matrix to unweighted mean 0 and unit standard
# deviation, the representation that crosses the FFI boundary.
standardize <- function(m) {
  centers <- colMeans(m)
  scales <- apply(m, 2, sd)
  scales[scales == 0] <- 1
  sweep(sweep(m, 2, centers, "-"), 2, scales, "/")
}

# Renormalize weights so that each group sums to one, matching the n_eff = 1
# normalization the harness solves under. `groups` is an integer group label per
# unit.
normalize_by_group <- function(weights, groups) {
  for (g in unique(groups)) {
    idx <- groups == g
    weights[idx] <- weights[idx] / sum(weights[idx])
  }
  weights
}

write_fixture <- function(name, fixture) {
  path <- file.path(output_dir, paste0(name, ".json"))
  writeLines(jsonlite::toJSON(fixture, auto_unbox = TRUE, digits = 17), path)
  message("wrote ", path)
}

# Build a discrete average-treatment-effect fixture: every group is balanced to
# the pooled covariate means, and each group is normalized to sum to one.
discrete_ate_fixture <- function(name, data, seed) {
  covariate_matrix <- standardize(model.matrix(~ x1 + x2 - 1, data))
  n <- nrow(covariate_matrix)
  p <- ncol(covariate_matrix)
  group_idx <- as.integer(as.factor(data$exposure)) - 1L

  reference <- weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "ebal",
    estimand = "ATE"
  )
  expected <- normalize_by_group(reference$weights, group_idx)

  write_fixture(
    name,
    list(
      kind = "discrete",
      n = n,
      p = p,
      covs = as.numeric(covariate_matrix),
      targets = as.numeric(colSums(covariate_matrix * 0)),
      tols = numeric(p),
      base = rep(1, n),
      s = rep(1, n),
      n_eff = 1,
      group_idx = group_idx,
      expected_weights = as.numeric(expected),
      rel_tol = 1e-6
    )
  )
}

# Build a discrete average-treatment-effect-on-the-treated fixture. The focal
# level keeps its base weights and the other groups are reweighted to the focal
# covariate means. The solver skips the focal rows, so their expected weights are
# zero; the reweighted rows are normalized so their sampling-weighted sum matches
# the focal total, the normalization the solver uses.
discrete_att_fixture <- function(name, data, focal) {
  covariate_matrix <- standardize(model.matrix(~ x1 + x2 - 1, data))
  n <- nrow(covariate_matrix)
  p <- ncol(covariate_matrix)
  exposure <- as.character(data$exposure)
  focal <- as.character(focal)
  is_focal <- exposure == focal
  n_eff <- sum(is_focal)
  targets <- colMeans(covariate_matrix[is_focal, , drop = FALSE])

  # Each non-focal level is an independent block reweighted to the focal means,
  # matching the per-group solve the R layer performs.
  levels_all <- if (is.factor(data$exposure)) {
    levels(data$exposure)
  } else {
    sort(unique(exposure))
  }
  non_focal <- setdiff(levels_all, focal)
  group_idx <- rep(-1L, n)
  for (block in seq_along(non_focal)) {
    group_idx[exposure == non_focal[block]] <- block - 1L
  }

  reference <- weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "ebal",
    estimand = "ATT",
    focal = focal
  )
  expected <- numeric(n)
  for (block in seq_along(non_focal)) {
    idx <- exposure == non_focal[block]
    expected[idx] <- reference$weights[idx] /
      sum(reference$weights[idx]) *
      n_eff
  }

  write_fixture(
    name,
    list(
      kind = "discrete",
      n = n,
      p = p,
      covs = as.numeric(covariate_matrix),
      targets = as.numeric(targets),
      tols = numeric(p),
      base = rep(1, n),
      s = rep(1, n),
      n_eff = n_eff,
      group_idx = as.integer(group_idx),
      expected_weights = as.numeric(expected),
      rel_tol = 1e-6
    )
  )
}

# Build a discrete average-treatment-effect fixture with sampling weights: every
# group is balanced to the sampling-weighted pooled covariate means, and each
# group is normalized so its sampling-weighted weight sum is one, the
# normalization the solver applies.
discrete_ate_sampling_fixture <- function(name, data, s) {
  covariate_matrix <- standardize(model.matrix(~ x1 + x2 - 1, data))
  n <- nrow(covariate_matrix)
  p <- ncol(covariate_matrix)
  group_idx <- as.integer(as.factor(data$exposure)) - 1L
  targets <- apply(covariate_matrix, 2, stats::weighted.mean, w = s)

  reference <- weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "ebal",
    estimand = "ATE",
    s.weights = s
  )
  expected <- reference$weights
  for (g in unique(group_idx)) {
    idx <- group_idx == g
    expected[idx] <- expected[idx] / sum(s[idx] * expected[idx])
  }

  write_fixture(
    name,
    list(
      kind = "discrete",
      n = n,
      p = p,
      covs = as.numeric(covariate_matrix),
      targets = as.numeric(targets),
      tols = numeric(p),
      base = rep(1, n),
      s = as.numeric(s),
      n_eff = 1,
      group_idx = group_idx,
      expected_weights = as.numeric(expected),
      rel_tol = 1e-6
    )
  )
}

# Build a discrete average-treatment-effect fixture with base weights: the base
# weights are the Kullback-Leibler prior the tilt anchors to, and every group is
# balanced to the unweighted pooled covariate means. Each group is normalized so
# its weight sum is one, the normalization the solver applies.
discrete_ate_base_fixture <- function(name, data, base) {
  covariate_matrix <- standardize(model.matrix(~ x1 + x2 - 1, data))
  n <- nrow(covariate_matrix)
  p <- ncol(covariate_matrix)
  group_idx <- as.integer(as.factor(data$exposure)) - 1L

  reference <- weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "ebal",
    estimand = "ATE",
    base.weight = base
  )
  expected <- normalize_by_group(reference$weights, group_idx)

  write_fixture(
    name,
    list(
      kind = "discrete",
      n = n,
      p = p,
      covs = as.numeric(covariate_matrix),
      targets = numeric(p),
      tols = numeric(p),
      base = as.numeric(base),
      s = rep(1, n),
      n_eff = 1,
      group_idx = group_idx,
      expected_weights = as.numeric(expected),
      rel_tol = 1e-6
    )
  )
}

# Build a continuous average-treatment-effect fixture: the exposure and covariate
# marginals are held to the sample and each exposure-covariate product is driven
# to zero, so the weighted correlations vanish.
continuous_ate_fixture <- function(name, data, seed) {
  covariate_matrix <- standardize(model.matrix(~ x1 + x2 - 1, data))
  n <- nrow(covariate_matrix)
  exposure <- data$exposure
  e <- (exposure - mean(exposure)) / sd(exposure)
  covs <- cbind(e, covariate_matrix, covariate_matrix * e)
  p <- ncol(covs)

  reference <- weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "ebal"
  )
  expected <- reference$weights / sum(reference$weights)

  write_fixture(
    name,
    list(
      kind = "continuous",
      n = n,
      p = p,
      covs = as.numeric(covs),
      targets = numeric(p),
      tols = numeric(p),
      base = rep(1, n),
      s = rep(1, n),
      n_eff = 1,
      expected_weights = as.numeric(expected),
      rel_tol = 1e-6
    )
  )
}

for (n in c(500L, 5000L)) {
  binary <- sim_binary(n, seed = 2024)
  categorical <- sim_categorical(n, seed = 2024)
  continuous <- sim_continuous(n, seed = 2024)

  discrete_ate_fixture(
    sprintf("entropy_binary_ate_n%d", n),
    binary,
    seed = 2024
  )
  discrete_ate_fixture(
    sprintf("entropy_categorical_ate_n%d", n),
    categorical,
    seed = 2024
  )
  continuous_ate_fixture(
    sprintf("entropy_continuous_ate_n%d", n),
    continuous,
    seed = 2024
  )

  # Estimand coverage: average treatment effect on the treated and, through a
  # focal untreated level, on the untreated.
  discrete_att_fixture(sprintf("entropy_binary_att_n%d", n), binary, focal = 1)
  discrete_att_fixture(sprintf("entropy_binary_atu_n%d", n), binary, focal = 0)
  discrete_att_fixture(
    sprintf("entropy_categorical_att_n%d", n),
    categorical,
    focal = "b"
  )

  # Sampling-weight coverage for the average treatment effect.
  set.seed(2024)
  sampling_weights <- runif(n, 0.3, 3)
  discrete_ate_sampling_fixture(
    sprintf("entropy_binary_ate_sweights_n%d", n),
    binary,
    sampling_weights
  )

  # Base-weight coverage for the average treatment effect.
  set.seed(99)
  base_weights <- runif(n, 0.5, 2)
  discrete_ate_base_fixture(
    sprintf("entropy_binary_ate_bweights_n%d", n),
    binary,
    base_weights
  )
}

message("Golden fixtures written to ", output_dir)
