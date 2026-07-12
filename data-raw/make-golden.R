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

# Build an inverse probability tilting fixture. The design matrix carries the
# intercept the propensity model needs, so `covs` is the intercept column
# followed by the standardized covariates and `p` counts the intercept. The Rust
# solver returns raw tilt weights, so the reference weights are rescaled to the
# raw convention the solver produces: for the average treatment effect each
# group's sampling-weighted sum matches the whole-sample total, and for a focal
# estimand the focal units carry weight one while the other groups match the
# focal total.
ipt_fixture <- function(name, data, estimand, link, focal = NULL, s = NULL) {
  covariate_matrix <- standardize(model.matrix(~ x1 + x2 - 1, data))
  n <- nrow(covariate_matrix)
  design <- cbind(1, covariate_matrix)
  p <- ncol(design)
  exposure <- as.character(data$exposure)
  levels_all <- if (is.factor(data$exposure)) {
    levels(data$exposure)
  } else {
    sort(unique(exposure))
  }
  treat <- match(exposure, levels_all) - 1L
  if (is.null(s)) {
    s <- rep(1, n)
  }

  args <- list(
    exposure ~ x1 + x2,
    data = data,
    method = "ipt",
    estimand = toupper(estimand),
    link = link
  )
  if (!is.null(focal)) {
    args$focal <- focal
  }
  if (any(s != 1)) {
    args$s.weights <- s
  }
  reference <- do.call(weightit, args)
  wt <- reference$weights

  expected <- numeric(n)
  fixture <- list(
    kind = "ipt",
    n = n,
    p = p,
    covs = as.numeric(design),
    treat = as.integer(treat),
    s = as.numeric(s),
    link = link,
    estimand = if (identical(estimand, "ate")) "ate" else "att",
    rel_tol = 1e-6
  )

  if (identical(estimand, "ate")) {
    total <- sum(s)
    for (g in unique(treat)) {
      idx <- treat == g
      expected[idx] <- wt[idx] / sum(s[idx] * wt[idx]) * total
    }
  } else {
    focal_idx <- match(as.character(focal), levels_all) - 1L
    is_focal <- treat == focal_idx
    n_eff <- sum(s[is_focal])
    expected[is_focal] <- 1
    for (g in setdiff(unique(treat), focal_idx)) {
      idx <- treat == g
      expected[idx] <- wt[idx] / sum(s[idx] * wt[idx]) * n_eff
    }
    fixture$focal <- focal_idx
  }

  fixture$expected_weights <- as.numeric(expected)
  write_fixture(name, fixture)
}

# The balancing weight a unit carries at propensity `p` given its treatment
# indicator `t`, per estimand. These mirror the core's weight forms so the
# reference weights can be rebuilt from a fitted propensity alone.
cbps_weight <- function(estimand, p, t) {
  switch(
    estimand,
    ate = t / p + (1 - t) / (1 - p),
    att = t + (1 - t) * p / (1 - p),
    atc = t * (1 - p) / p + (1 - t),
    ato = t * (1 - p) + (1 - t) * p
  )
}

# Build a binary just-identified covariate balancing propensity score fixture.
# The reference is WeightIt's exactly balancing (over = FALSE) fit; because the
# just-identified moment conditions pin the propensity, the raw solver weights
# are reconstructed from the reference fitted propensity through the estimand's
# weight form. The design matrix carries the intercept the propensity model
# needs, and `p` counts it.
cbps_binary_fixture <- function(name, data, estimand, s = NULL) {
  covariate_matrix <- standardize(model.matrix(~ x1 + x2 - 1, data))
  n <- nrow(covariate_matrix)
  design <- cbind(1, covariate_matrix)
  p <- ncol(design)
  treat <- as.integer(data$exposure)
  if (is.null(s)) {
    s <- rep(1, n)
  }
  focal <- switch(estimand, att = "1", atc = "0", NULL)

  args <- list(
    exposure ~ x1 + x2,
    data = data,
    method = "cbps",
    estimand = toupper(estimand),
    over = FALSE
  )
  if (!is.null(focal)) {
    args$focal <- focal
  }
  if (any(s != 1)) {
    args$s.weights <- s
  }
  reference <- do.call(weightit, args)
  expected <- cbps_weight(estimand, reference$ps, treat)

  write_fixture(
    name,
    list(
      kind = "cbps",
      n = n,
      p = p,
      covs = as.numeric(design),
      treat = as.integer(treat),
      s = as.numeric(s),
      estimand = estimand,
      link = "logit",
      over = FALSE,
      twostep = TRUE,
      expected_weights = as.numeric(expected),
      rel_tol = 1e-6
    )
  )
}

# Build a binary over-identified covariate balancing propensity score fixture.
# WeightIt does not expose the generalized-method-of-moments criterion value, so
# the reference objective is reconstructed here from WeightIt's fitted
# coefficients using the core's own moment and weighting definitions: the
# score-and-balance moment stack, and the two-step weighting matrix formed as the
# eigenvalue-thresholded pseudo-inverse of the moment covariance at the
# maximum-likelihood anchor. The design's tolerance policy requires our solver's
# criterion to sit at or below this reference plus 1e-8, which holds because the
# solver minimizes the same criterion the reconstruction evaluates.
cbps_over_fixture <- function(name, data, twostep = TRUE) {
  covariate_matrix <- standardize(model.matrix(~ x1 + x2 - 1, data))
  n <- nrow(covariate_matrix)
  design <- cbind(1, covariate_matrix)
  p <- ncol(design)
  treat <- as.integer(data$exposure)
  s <- rep(1, n)
  m_total <- 2 * p

  # The propensity is held inside the open unit interval, matching the clamp the
  # criterion applies so a boundary excursion cannot overflow it.
  clamp_prob <- function(pp) pmin(pmax(pp, 1e-8), 1 - 1e-8)

  # The stacked moments: the score residual on the model covariates, then the
  # average-treatment-effect balancing factor times the balance covariates.
  moments <- function(beta) {
    eta <- as.vector(design %*% beta)
    prob <- clamp_prob(stats::plogis(eta))
    factor <- (prob - treat) / (prob * (1 - prob))
    g <- matrix(0, n, m_total)
    for (j in seq_len(p)) {
      g[, j] <- (treat - prob) * design[, j]
      g[, p + j] <- factor * design[, j]
    }
    g
  }

  # The eigenvalue-thresholded symmetric pseudo-inverse the core uses, with the
  # same relative condition floor.
  pseudo_inverse <- function(a) {
    eig <- eigen((a + t(a)) / 2, symmetric = TRUE)
    floor <- 1e-12 * max(abs(eig$values))
    inv <- ifelse(eig$values > floor, 1 / eig$values, 0)
    eig$vectors %*% diag(inv, m_total, m_total) %*% t(eig$vectors)
  }

  anchor <- stats::glm.fit(
    design,
    treat,
    family = stats::binomial()
  )$coefficients
  g_anchor <- moments(anchor)
  covariance <- crossprod(g_anchor * s, g_anchor) / n
  weighting <- pseudo_inverse(covariance)

  reference <- weightit(
    exposure ~ x1 + x2,
    data = data,
    method = "cbps",
    estimand = "ATE",
    over = TRUE,
    twostep = twostep
  )
  # Recover the coefficients that reproduce WeightIt's fitted propensity in the
  # standardized design, so the criterion is evaluated at WeightIt's fit.
  eta_reference <- stats::qlogis(clamp_prob(reference$ps))
  beta_reference <- solve(crossprod(design), crossprod(design, eta_reference))
  m_reference <- colSums(moments(beta_reference) * s) / n
  expected_obj <- as.numeric(t(m_reference) %*% weighting %*% m_reference)

  write_fixture(
    name,
    list(
      kind = "cbps",
      n = n,
      p = p,
      covs = as.numeric(design),
      treat = as.integer(treat),
      s = as.numeric(s),
      estimand = "ate",
      link = "logit",
      over = TRUE,
      twostep = twostep,
      expected_obj = expected_obj
    )
  )
}

# Build a categorical covariate balancing propensity score fixture. The
# just-identified categorical form solves the inverse probability tilting moment
# conditions, so the reference is WeightIt's inverse probability tilting fit and
# the raw weights follow the same convention the tilting fixtures use: the
# average treatment effect matches each group's sampling-weighted sum to the
# whole-sample total, and a focal estimand carries the focal units at weight one
# while the other groups match the focal total.
cbps_multi_fixture <- function(name, data, estimand, focal = NULL) {
  covariate_matrix <- standardize(model.matrix(~ x1 + x2 - 1, data))
  n <- nrow(covariate_matrix)
  design <- cbind(1, covariate_matrix)
  p <- ncol(design)
  exposure <- as.character(data$exposure)
  levels_all <- if (is.factor(data$exposure)) {
    levels(data$exposure)
  } else {
    sort(unique(exposure))
  }
  treat <- match(exposure, levels_all) - 1L
  s <- rep(1, n)

  args <- list(
    exposure ~ x1 + x2,
    data = data,
    method = "ipt",
    estimand = toupper(estimand)
  )
  if (!is.null(focal)) {
    args$focal <- focal
  }
  reference <- do.call(weightit, args)
  wt <- reference$weights

  expected <- numeric(n)
  fixture <- list(
    kind = "cbps_multi",
    n = n,
    p = p,
    covs = as.numeric(design),
    treat = as.integer(treat),
    s = as.numeric(s),
    estimand = if (identical(estimand, "ate")) "ate" else "att",
    link = "logit",
    rel_tol = 1e-6
  )

  if (identical(estimand, "ate")) {
    total <- sum(s)
    for (g in unique(treat)) {
      idx <- treat == g
      expected[idx] <- wt[idx] / sum(s[idx] * wt[idx]) * total
    }
  } else {
    focal_idx <- match(as.character(focal), levels_all) - 1L
    is_focal <- treat == focal_idx
    n_eff <- sum(s[is_focal])
    expected[is_focal] <- 1
    for (g in setdiff(unique(treat), focal_idx)) {
      idx <- treat == g
      expected[idx] <- wt[idx] / sum(s[idx] * wt[idx]) * n_eff
    }
    fixture$focal <- focal_idx
  }

  fixture$expected_weights <- as.numeric(expected)
  write_fixture(name, fixture)
}

# Build a continuous covariate balancing propensity score fixture. WeightIt's
# continuous fits solve a parametric criterion rather than the exact covariance
# balancing this core implements, so no WeightIt method reproduces its weights.
# The reference is instead an independent solve of the same convex dual: the
# minimum-divergence exponential tilt whose intercept column carries the
# exposure-mean condition and whose covariate columns carry the covariance
# conditions. An independent optimizer reaching the same unique optimum is a
# genuine solver-parity check.
cbps_continuous_fixture <- function(name, data) {
  covariate_matrix <- standardize(model.matrix(~ x1 + x2 - 1, data))
  n <- nrow(covariate_matrix)
  design <- cbind(1, covariate_matrix)
  p <- ncol(design)
  s <- rep(1, n)
  expo <- data$exposure

  sbar <- sum(s)
  m_expo <- sum(s * expo) / sbar
  sd_e <- sqrt(sum(s * (expo - m_expo)^2) / sbar)
  if (sd_e == 0) {
    sd_e <- 1e-8
  }
  centered_expo <- (expo - m_expo) / sd_e
  xbar <- colSums(s * design) / sbar
  is_intercept <- apply(design, 2, function(col) {
    first <- col[1]
    first != 0 && all(abs(col - first) < 1e-12)
  })

  feat <- matrix(0, n, p)
  for (j in seq_len(p)) {
    feat[, j] <- if (is_intercept[j]) {
      centered_expo
    } else {
      (design[, j] - xbar[j]) * centered_expo
    }
  }

  objective <- function(gamma) log(sum(s * exp(-(feat %*% gamma))))
  gradient <- function(gamma) {
    w <- s * exp(-(feat %*% gamma))
    -colSums(as.vector(w) * feat) / sum(w)
  }
  fit <- stats::optim(
    rep(0, p),
    objective,
    gradient,
    method = "BFGS",
    control = list(reltol = 1e-15, maxit = 2000)
  )
  weights <- as.vector(exp(-(feat %*% fit$par)))
  weights <- weights * sbar / sum(s * weights)

  write_fixture(
    name,
    list(
      kind = "cbps_cont",
      n = n,
      p = p,
      covs = as.numeric(design),
      expo = as.numeric(expo),
      s = as.numeric(s),
      expected_weights = as.numeric(weights),
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

  # Inverse probability tilting: binary and categorical, the average treatment
  # effect and a focal estimand, plus sampling-weight coverage.
  ipt_fixture(sprintf("ipt_binary_ate_n%d", n), binary, "ate", "logit")
  ipt_fixture(
    sprintf("ipt_binary_att_n%d", n),
    binary,
    "att",
    "logit",
    focal = 1
  )
  ipt_fixture(
    sprintf("ipt_categorical_ate_n%d", n),
    categorical,
    "ate",
    "logit"
  )
  ipt_fixture(
    sprintf("ipt_categorical_att_n%d", n),
    categorical,
    "att",
    "logit",
    focal = "b"
  )
  ipt_fixture(
    sprintf("ipt_binary_ate_sweights_n%d", n),
    binary,
    "ate",
    "logit",
    s = sampling_weights
  )

  # Covariate balancing propensity score. Binary just-identified over the four
  # estimands, the over-identified two-step criterion, the categorical average
  # treatment effect and a focal estimand, the continuous exposure, and
  # sampling-weight coverage.
  cbps_binary_fixture(sprintf("cbps_binary_ate_n%d", n), binary, "ate")
  cbps_binary_fixture(sprintf("cbps_binary_att_n%d", n), binary, "att")
  cbps_binary_fixture(sprintf("cbps_binary_ato_n%d", n), binary, "ato")
  cbps_binary_fixture(
    sprintf("cbps_binary_ate_sweights_n%d", n),
    binary,
    "ate",
    s = sampling_weights
  )
  cbps_over_fixture(sprintf("cbps_binary_over_n%d", n), binary, twostep = TRUE)
  cbps_multi_fixture(
    sprintf("cbps_categorical_ate_n%d", n),
    categorical,
    "ate"
  )
  cbps_multi_fixture(
    sprintf("cbps_categorical_att_n%d", n),
    categorical,
    "att",
    focal = "b"
  )
  cbps_continuous_fixture(sprintf("cbps_continuous_ate_n%d", n), continuous)
}

message("Golden fixtures written to ", output_dir)
