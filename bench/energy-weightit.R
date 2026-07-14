# End-to-end energy balancing benchmark against WeightIt's energy method.
#
# Measures balance(method = bw_energy()) against the equivalent
# WeightIt::weightit(method = "energy") call for a binary exposure and the
# average treatment effect. Both packages default to the same specification:
# a scaled-Euclidean covariate distance, an L2 weight penalty of 1e-4, a
# minimum weight of 1e-8, the improved between-group variant, and an OSQP
# solve. The comparison is therefore as-shipped versus as-shipped, differing
# only in the implementation behind the same objective.
#
# Energy balancing is a quadratic program whose quadratic term is the dense
# n by n covariate distance matrix, so both memory and solve time grow far
# faster than the estimating-equation family. The grid records that scaling.
#
# The improved average-treatment-effect energy value of each weight vector is
# evaluated on a shared scaled-Euclidean distance matrix and the two solutions
# are required to agree, with ours no worse than WeightIt by more than 1e-6,
# before any timing is recorded. A speedup is reported only for a correct fit.
#
# This script is Rbuildignored and is not part of the installed package. Run it
# against an installed release build, not a development load.

suppressMessages({
  library(balancing)
  library(WeightIt)
  library(bench)
})

# Deterministic synthetic data: p standard-normal covariates and a binary
# exposure whose propensity depends on a handful of them, so the groups start
# imbalanced and the reweighting has real work to do.
make_data <- function(n, p, seed = 20240713L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  lp <- 0.4 * X[, 1] - 0.3 * X[, 2] + 0.2 * X[, 3] - 0.2 * X[, 4]
  exposure <- stats::rbinom(n, 1, stats::plogis(lp))
  list(df = data.frame(exposure = exposure, X), X = X, tr = exposure)
}

# Scaled-Euclidean pairwise distance: each covariate divided by its sample
# standard deviation, then the Euclidean distance between rows. This is the
# metric both packages build the energy objective on by default.
scaled_euclidean <- function(X) {
  sds <- apply(X, 2, stats::sd)
  sds[sds == 0] <- 1
  Xs <- sweep(X, 2, sds, "/")
  as.matrix(stats::dist(Xs))
}

# The energy distance of a group to the full sample, with the group weights
# normalized to sum to one. Matches the objective the solver minimizes: the
# cross term to the whole sample minus the within-group term.
energy_to_sample <- function(D, idx, w) {
  n <- nrow(D)
  wn <- w[idx] / sum(w[idx])
  cross <- (2 / n) * sum(wn * rowSums(D[idx, , drop = FALSE]))
  within <- as.numeric(t(wn) %*% D[idx, idx, drop = FALSE] %*% wn)
  cross - within
}

# The between-group energy distance of the improved variant, with each group's
# weights normalized to sum to one.
energy_between <- function(D, idx_a, idx_b, w) {
  wa <- w[idx_a] / sum(w[idx_a])
  wb <- w[idx_b] / sum(w[idx_b])
  cross <- 2 * as.numeric(t(wa) %*% D[idx_a, idx_b, drop = FALSE] %*% wb)
  wa2 <- as.numeric(t(wa) %*% D[idx_a, idx_a, drop = FALSE] %*% wa)
  wb2 <- as.numeric(t(wb) %*% D[idx_b, idx_b, drop = FALSE] %*% wb)
  cross - wa2 - wb2
}

# The improved average-treatment-effect energy objective of a weight vector on
# a shared distance matrix. Lower is better balance. Group-internal weight
# scale cancels, so raw weights on any per-group scale compare directly.
energy_ate_objective <- function(D, tr, w) {
  idx_t <- which(tr == 1)
  idx_c <- which(tr == 0)
  energy_to_sample(D, idx_t, w) +
    energy_to_sample(D, idx_c, w) +
    energy_between(D, idx_t, idx_c, w)
}

fit_ours <- function(df, covs) {
  as.numeric(weights(balance(
    df,
    exposure,
    all_of(covs),
    method = bw_energy(),
    estimand = "ate"
  )))
}

fit_weightit <- function(df, covs) {
  f <- stats::as.formula(paste("exposure ~", paste(covs, collapse = " + ")))
  wt <- weightit(f, data = df, method = "energy", estimand = "ATE")
  wt$weights
}

run <- function(sizes = c(500L, 2000L, 5000L), p = 4L) {
  results <- list()
  parity <- list()

  for (n in sizes) {
    data <- make_data(n, p)
    df <- data$df
    X <- data$X
    tr <- data$tr
    covs <- paste0("x", seq_len(p))

    # As-shipped versus as-shipped: neither package is thread-restricted here,
    # so each runs at its own default. WeightIt's energy solve is single
    # threaded; ours uses its default worker count for distance assembly.
    w_ours <- fit_ours(df, covs)
    w_wi <- fit_weightit(df, covs)

    D <- scaled_euclidean(X)
    obj_ours <- energy_ate_objective(D, tr, w_ours)
    obj_wi <- energy_ate_objective(D, tr, w_wi)
    corr <- stats::cor(
      w_ours / mean(w_ours[tr == 1]),
      w_wi / mean(w_wi[tr == 1])
    )

    parity[[as.character(n)]] <- list(
      obj_ours = obj_ours,
      obj_weightit = obj_wi,
      obj_gap = obj_ours - obj_wi,
      correlation = corr
    )
    cat(sprintf(
      "[n=%d] energy obj ours=%.8e weightit=%.8e gap=%.2e corr=%.5f\n",
      n,
      obj_ours,
      obj_wi,
      obj_ours - obj_wi,
      corr
    ))

    # Ours must be no worse than WeightIt by more than 1e-6 on the shared
    # distance matrix.
    stopifnot(obj_ours <= obj_wi + 1e-6)

    iters <- if (n >= 5000L) 3L else 5L
    mk <- bench::mark(
      ours = fit_ours(df, covs),
      weightit = fit_weightit(df, covs),
      check = FALSE,
      min_iterations = iters,
      max_iterations = iters,
      filter_gc = FALSE
    )
    results[[as.character(n)]] <- mk
  }

  list(results = results, parity = parity, sizes = sizes, p = p)
}

if (sys.nframe() == 0L) {
  out <- run()
  dir.create(
    file.path("scratch", "bench", "raw"),
    showWarnings = FALSE,
    recursive = TRUE
  )
  saveRDS(
    out,
    file.path("scratch", "bench", "raw", "rside-energy-weightit.rds")
  )
  for (nm in names(out$results)) {
    cat("\n== n = ", nm, " ==\n", sep = "")
    print(out$results[[nm]][, c(
      "expression",
      "min",
      "median",
      "itr/sec",
      "mem_alloc",
      "n_itr"
    )])
  }
}
