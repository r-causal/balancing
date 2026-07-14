# End-to-end stable balancing weights benchmark against optweight.
#
# Measures balance(method = sbw()) against the equivalent optweight() call for a
# binary exposure and the average treatment effect, plus one continuous case.
# Both packages solve the same quadratic program: minimize the sum of squared
# weights subject to each reweighted group's covariate means falling inside a
# standardized-mean-difference band. Both default to the L2 norm, a minimum
# weight of 1e-8, and an OSQP solve, so the comparison is as-shipped versus
# as-shipped, differing only in the implementation behind the same objective and
# the same feasible set.
#
# Correctness is checked before any timing: each solution must hold every
# covariate's weighted standardized mean difference inside the requested band,
# and our weight-dispersion objective must be no worse than optweight's by more
# than a small tolerance. A speedup is recorded only for a correct fit.
#
# The continuous fit runs a correlation-refinement loop in R that re-solves the
# quadratic program a few times, tightening the effective tolerance until the
# true weighted correlation sits inside the band. The loop's extra solves are
# counted by tracing the internal solver entry point, and their cost share is
# reported alongside a single-solve time.
#
# This script is Rbuildignored and is not part of the installed package. Run it
# against an installed release build, not a development load.

suppressMessages({
  library(balancing)
  library(optweight)
  library(bench)
})

# Deterministic synthetic data: p standard-normal covariates and a binary
# exposure whose propensity depends on a few of them, so the groups start
# imbalanced but overlap in covariate space and the moment balance is feasible.
make_data <- function(n, p, seed = 20240714L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  lp <- 0.4 * X[, 1] - 0.3 * X[, 2] + 0.2 * X[, 3] - 0.2 * X[, 4]
  exposure <- stats::rbinom(n, 1, stats::plogis(lp))
  list(df = data.frame(exposure = exposure, X), X = X, tr = exposure)
}

# A continuous exposure correlated with the covariates, so uniform weights leave
# a weighted exposure-covariate correlation the constraint rows must pull inside
# the band.
make_data_cont <- function(n, p, seed = 20240714L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  a <- 0.6 * X[, 1] + 0.3 * X[, 2] + stats::rnorm(n)
  list(df = data.frame(exposure = a, X), X = X, a = a)
}

# The worst weighted standardized mean difference of any covariate between the
# reweighted groups and the pooled target, the quantity the balance tolerance
# bounds. Weights are used as returned; standardization is by the pooled sample
# standard deviation.
worst_smd <- function(X, tr, w) {
  sds <- apply(X, 2, stats::sd)
  sds[sds == 0] <- 1
  Xs <- sweep(X, 2, colMeans(X), "-")
  Xs <- sweep(Xs, 2, sds, "/")
  target <- colMeans(Xs)
  worst <- 0
  for (g in unique(tr)) {
    idx <- tr == g
    wm <- colSums(Xs[idx, , drop = FALSE] * w[idx]) / sum(w[idx])
    worst <- max(worst, max(abs(wm - target)))
  }
  worst
}

# The weight-dispersion objective: within each group the weights are normalized
# to mean one, and the objective is the summed squared deviation from one. This
# is the quantity both solvers minimize; the group-internal scale is removed so
# the two packages' normalization conventions compare directly.
dispersion <- function(tr, w) {
  total <- 0
  for (g in unique(tr)) {
    idx <- tr == g
    wg <- w[idx] / mean(w[idx])
    total <- total + sum(wg^2)
  }
  total
}

fit_ours <- function(df, covs, tol) {
  as.numeric(weights(balance(
    df,
    exposure,
    all_of(covs),
    method = sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = tol)
  )))
}

fit_optweight <- function(df, covs, tol) {
  f <- stats::as.formula(paste("exposure ~", paste(covs, collapse = " + ")))
  ow <- optweight(
    f,
    data = df,
    tols = tol,
    estimand = "ATE",
    norm = "l2",
    min.w = 1e-8
  )
  ow$weights
}

run_binary <- function(sizes = c(500L, 2000L, 5000L), p = 4L, tol = 0.05) {
  results <- list()
  parity <- list()

  for (n in sizes) {
    data <- make_data(n, p)
    df <- data$df
    X <- data$X
    tr <- data$tr
    covs <- paste0("x", seq_len(p))

    w_ours <- fit_ours(df, covs, tol)
    w_ow <- fit_optweight(df, covs, tol)

    smd_ours <- worst_smd(X, tr, w_ours)
    smd_ow <- worst_smd(X, tr, w_ow)
    obj_ours <- dispersion(tr, w_ours)
    obj_ow <- dispersion(tr, w_ow)

    parity[[as.character(n)]] <- list(
      smd_ours = smd_ours,
      smd_optweight = smd_ow,
      obj_ours = obj_ours,
      obj_optweight = obj_ow,
      obj_ratio = obj_ours / obj_ow,
      tol = tol
    )
    cat(sprintf(
      "[n=%d tol=%.3f] smd ours=%.4f ow=%.4f | disp ours=%.2f ow=%.2f ratio=%.4f\n",
      n,
      tol,
      smd_ours,
      smd_ow,
      obj_ours,
      obj_ow,
      obj_ours / obj_ow
    ))

    # Both solutions must respect the band, and ours must be no worse on the
    # dispersion objective by more than a small relative slack.
    stopifnot(smd_ours <= tol + 1e-3)
    stopifnot(smd_ow <= tol + 1e-3)
    stopifnot(obj_ours <= obj_ow * (1 + 1e-3) + 1e-6)

    iters <- if (n >= 5000L) 5L else 10L
    mk <- bench::mark(
      ours = fit_ours(df, covs, tol),
      optweight = fit_optweight(df, covs, tol),
      check = FALSE,
      min_iterations = iters,
      max_iterations = iters,
      filter_gc = FALSE
    )
    results[[as.character(n)]] <- mk
  }

  list(results = results, parity = parity, sizes = sizes, p = p, tol = tol)
}

# The continuous fit: measure the full refinement-loop fit, count how many times
# the internal solver runs during one fit by tracing it, and time a single
# solver call so the loop's cost share can be attributed.
run_continuous <- function(n = 2000L, p = 4L, tol = 0.1) {
  data <- make_data_cont(n, p)
  df <- data$df
  covs <- paste0("x", seq_len(p))

  fit_cont <- function() {
    as.numeric(weights(balance(
      df,
      exposure,
      all_of(covs),
      method = sbw(),
      estimand = "ate",
      constraints = balance_terms(tolerance = tol)
    )))
  }

  # Count solver passes in one fit by tracing the internal continuous entry
  # point. Tracing is runtime instrumentation and does not alter the package.
  pass_count <- 0L
  suppressMessages(trace(
    "solve_sbw_cont",
    tracer = quote({
      assign(
        "pass_count",
        get("pass_count", envir = .GlobalEnv) + 1L,
        envir = .GlobalEnv
      )
    }),
    where = asNamespace("balancing"),
    print = FALSE
  ))
  assign("pass_count", 0L, envir = .GlobalEnv)
  w_cont <- fit_cont()
  passes <- get("pass_count", envir = .GlobalEnv)
  suppressMessages(untrace("solve_sbw_cont", where = asNamespace("balancing")))

  # A single solver call at the requested tolerance, timed on its own, so the
  # per-pass cost is isolated from the correlation evaluation and the balance
  # bookkeeping.
  z <- scale(as.matrix(df[covs]))
  s <- rep(1, n)
  solve_sbw_cont <- getFromNamespace("solve_sbw_cont", "balancing")
  single <- bench::mark(
    one = solve_sbw_cont(
      as.numeric(df$exposure),
      z,
      s,
      "l2",
      rep(tol, p),
      1e-8,
      list(threads = 1L)
    ),
    check = FALSE,
    min_iterations = 10L,
    max_iterations = 10L,
    filter_gc = FALSE
  )

  full <- bench::mark(
    full = fit_cont(),
    check = FALSE,
    min_iterations = 10L,
    max_iterations = 10L,
    filter_gc = FALSE
  )

  full_med <- as.numeric(full$median)
  single_med <- as.numeric(single$median)
  cat(sprintf(
    "[continuous n=%d tol=%.2f] passes=%d full_median=%.4fs single_solve_median=%.4fs\n",
    n,
    tol,
    passes,
    full_med,
    single_med
  ))
  cat(sprintf(
    "  solver share of full fit ~ %.1f%% (%d passes x single solve / full)\n",
    100 * passes * single_med / full_med,
    passes
  ))

  list(
    passes = passes,
    full_median = full_med,
    single_median = single_med,
    n = n,
    tol = tol
  )
}

if (sys.nframe() == 0L) {
  options(balancing.threads = 1L)
  bin <- run_binary()
  cont <- run_continuous()
  dir.create(
    file.path("scratch", "bench", "raw"),
    showWarnings = FALSE,
    recursive = TRUE
  )
  saveRDS(
    list(binary = bin, continuous = cont),
    file.path("scratch", "bench", "raw", "rside-sbw-optweight.rds")
  )
  for (nm in names(bin$results)) {
    cat("\n== n = ", nm, " ==\n", sep = "")
    print(bin$results[[nm]][, c(
      "expression",
      "min",
      "median",
      "itr/sec",
      "mem_alloc",
      "n_itr"
    )])
  }
}
