# End-to-end entropy balancing benchmark against WeightIt's ebal.
#
# Measures balance(method = entropy_balance()) against the equivalent
# WeightIt::weightit(method = "ebal") call at the reference scale used by the
# performance gate: n = 50000 units and 200 first-moment constraints, for the
# ATE and ATT estimands. Weight parity and balance precision are verified before
# any timing is recorded, so a speedup is reported only for correct results.
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
make_data <- function(n, p, seed = 20240711L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  lp <- 0.3 * X[, 1] - 0.3 * X[, 2] + 0.2 * X[, 3] - 0.2 * X[, 4]
  exposure <- stats::rbinom(n, 1, stats::plogis(lp))
  list(df = data.frame(exposure = exposure, X), X = X, tr = exposure)
}

# Largest absolute deviation of a weighted group mean from its target. For the
# ATE both groups target the overall covariate means; for the ATT the control
# group targets the treated covariate means and the treated group is unweighted.
balance_precision <- function(w, X, tr, estimand, focal = 1) {
  if (estimand == "ate") {
    tgt <- colMeans(X)
    m1 <- colSums(w[tr == 1] * X[tr == 1, ]) / sum(w[tr == 1])
    m0 <- colSums(w[tr == 0] * X[tr == 0, ]) / sum(w[tr == 0])
    max(abs(c(m1 - tgt, m0 - tgt)))
  } else {
    tgt <- colMeans(X[tr == focal, ])
    other <- 1 - focal
    mo <- colSums(w[tr == other] * X[tr == other, ]) / sum(w[tr == other])
    max(abs(mo - tgt))
  }
}

# Per-group mean-one normalization, so the two solvers' weights are compared on
# the same scale regardless of each one's internal normalization convention.
norm_groups <- function(w, tr) {
  out <- w
  for (g in unique(tr)) {
    idx <- tr == g
    out[idx] <- w[idx] / mean(w[idx])
  }
  out
}

fit_ours <- function(df, covs, estimand) {
  as.numeric(weights(balance(
    df,
    exposure,
    all_of(covs),
    method = entropy_balance(),
    estimand = estimand
  )))
}

fit_weightit <- function(df, covs, estimand) {
  f <- stats::as.formula(paste("exposure ~", paste(covs, collapse = " + ")))
  wt <- weightit(
    f,
    data = df,
    method = "ebal",
    estimand = toupper(estimand)
  )
  wt$weights
}

run <- function(n = 50000L, p = 200L, iterations = 5L) {
  data <- make_data(n, p)
  df <- data$df
  X <- data$X
  tr <- data$tr
  covs <- paste0("x", seq_len(p))

  results <- list()
  parity <- list()

  for (estimand in c("ate", "att")) {
    # Single-threaded, like-for-like: WeightIt is single-threaded, so cap our
    # core to one worker through the documented thread option.
    options(balancing.threads = 1L)
    w_ours <- fit_ours(df, covs, estimand)
    w_wi <- fit_weightit(df, covs, estimand)

    prec_ours <- balance_precision(w_ours, X, tr, estimand)
    prec_wi <- balance_precision(w_wi, X, tr, estimand)
    a <- norm_groups(w_ours, tr)
    b <- norm_groups(w_wi, tr)
    rel <- max(abs(a - b) / pmax(abs(b), 1e-8))
    parity[[estimand]] <- list(
      prec_ours = prec_ours,
      prec_weightit = prec_wi,
      max_rel_weight_diff = rel,
      correlation = stats::cor(w_ours, w_wi)
    )
    cat(sprintf(
      "[%s] balance precision ours=%.3e weightit=%.3e | max rel weight diff=%.3e\n",
      estimand,
      prec_ours,
      prec_wi,
      rel
    ))

    stopifnot(prec_ours < 1e-6, prec_wi < 1e-6, rel < 1e-6)

    single <- bench::mark(
      ours = fit_ours(df, covs, estimand),
      weightit = fit_weightit(df, covs, estimand),
      check = FALSE,
      min_iterations = iterations,
      max_iterations = iterations,
      filter_gc = FALSE
    )
    results[[paste0(estimand, "_single")]] <- single

    # Threaded row for our method only, reported separately, at 8 workers to
    # match the 8-thread criterion points.
    options(balancing.threads = 8L)
    threaded <- bench::mark(
      ours_threaded = fit_ours(df, covs, estimand),
      check = FALSE,
      min_iterations = iterations,
      max_iterations = iterations,
      filter_gc = FALSE
    )
    results[[paste0(estimand, "_threaded")]] <- threaded
  }

  list(
    results = results,
    parity = parity,
    n = n,
    p = p,
    iterations = iterations
  )
}

if (sys.nframe() == 0L) {
  out <- run()
  saveRDS(
    out,
    file.path("scratch", "bench", "raw", "rside-entropy-weightit.rds")
  )
  for (nm in names(out$results)) {
    cat("\n== ", nm, " ==\n", sep = "")
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
