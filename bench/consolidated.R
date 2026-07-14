# Consolidated end-to-end benchmark across every balancing method.
#
# One script covering the six user-facing methods on shared synthetic datasets:
# entropy balancing, inverse probability tilting, the covariate balancing
# propensity score, energy balancing, characteristic function distance
# balancing, and stable balancing weights. It measures each method's own
# end-to-end median through balance() at method-appropriate sizes, re-verifies
# the entropy external reference gate against WeightIt as-shipped, characterizes
# the overhead of the kernel path against energy balancing, compares the
# quadratic-program backends for the kernel method at matched tolerances, and
# spot-checks the stable-balancing-weights certificate fallback at the size that
# triggers it.
#
# Correctness is asserted before any timing: estimating-equation fits must
# converge and hit their balance target, and quadratic-program fits must produce
# a solved status with finite objective. A timing is recorded only for a correct
# fit.
#
# The data generators reuse the exact seeds, covariate counts, and propensity
# coefficients of the per-slice bench scripts (entropy-weightit.R,
# energy-weightit.R, sbw-optweight.R) at the sizes that carry a baseline entry,
# so the medians measure the same workload as the locked baselines and the
# regression comparison is like-for-like.
#
# This script is Rbuildignored and is not part of the installed package. Run it
# against an installed release build, not a development load.

suppressMessages({
  library(balancing)
  library(bench)
})

options(balancing.quiet = TRUE)

RAW_DIR <- file.path("scratch", "bench", "raw")

# ---- Shared synthetic data -------------------------------------------------

# Binary exposure with p standard-normal covariates whose propensity depends on
# a handful of them. The seed and coefficients match the per-slice scripts so a
# given (n, p) reproduces the workload behind the locked baselines.
make_binary <- function(
  n,
  p,
  seed = 20240711L,
  coefs = c(0.3, -0.3, 0.2, -0.2)
) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  k <- min(length(coefs), p)
  lp <- as.numeric(X[, seq_len(k), drop = FALSE] %*% coefs[seq_len(k)])
  exposure <- stats::rbinom(n, 1, stats::plogis(lp))
  list(df = data.frame(exposure = exposure, X), X = X, tr = exposure)
}

# Three-level categorical exposure from a multinomial logit on the covariates.
make_categorical <- function(n, p, seed = 20240715L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  eta2 <- 0.4 * X[, 1] - 0.2 * X[, 2]
  eta3 <- -0.3 * X[, 1] + 0.3 * X[, 3]
  denom <- 1 + exp(eta2) + exp(eta3)
  p1 <- 1 / denom
  p2 <- exp(eta2) / denom
  u <- stats::runif(n)
  g <- ifelse(u < p1, "a", ifelse(u < p1 + p2, "b", "c"))
  list(df = data.frame(exposure = factor(g), X), X = X, tr = g)
}

# Continuous exposure correlated with the covariates.
make_continuous <- function(n, p, seed = 20240714L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  a <- 0.6 * X[, 1] + 0.3 * X[, 2] + stats::rnorm(n)
  list(df = data.frame(exposure = a, X), X = X, a = a)
}

# ---- Correctness checks ----------------------------------------------------

# Largest absolute deviation of a weighted group mean from the ate target
# (overall covariate means). Used to confirm an estimating-equation fit reaches
# exact first-moment balance before it is timed.
ate_precision_binary <- function(w, X, tr) {
  tgt <- colMeans(X)
  m1 <- colSums(w[tr == 1] * X[tr == 1, ]) / sum(w[tr == 1])
  m0 <- colSums(w[tr == 0] * X[tr == 0, ]) / sum(w[tr == 0])
  max(abs(c(m1 - tgt, m0 - tgt)))
}

fit_weights <- function(...) {
  as.numeric(weights(balance(...)))
}

# ---- Per-method internal medians -------------------------------------------

# Each entry fits balance() once, checks correctness, then times it. Sizes are
# method-appropriate: the estimating-equation family scales to tens of
# thousands of units, while the dense quadratic-program family (energy, cfd)
# stays in the low thousands where the n by n solve is tractable.

bench_esteq_method <- function(
  label,
  method_fn,
  data,
  covs,
  estimand = "ate",
  iterations = 5L,
  check_precision = TRUE
) {
  df <- data$df
  fit <- balance(
    df,
    exposure,
    all_of(covs),
    method = method_fn(),
    estimand = estimand
  )
  w <- as.numeric(weights(fit))
  stopifnot(isTRUE(fit@converged))
  if (check_precision && is.numeric(data$tr)) {
    prec <- ate_precision_binary(w, data$X, data$tr)
    stopifnot(prec < 1e-6)
  }
  mk <- bench::mark(
    fit = fit_weights(
      df,
      exposure,
      all_of(covs),
      method = method_fn(),
      estimand = estimand
    ),
    check = FALSE,
    min_iterations = iterations,
    max_iterations = iterations,
    filter_gc = FALSE
  )
  list(label = label, median = as.numeric(mk$median), mark = mk)
}

bench_qp_method <- function(
  label,
  method_call,
  df,
  covs,
  estimand = "ate",
  constraints = NULL,
  iterations = 5L
) {
  do_fit <- function() {
    balance(
      df,
      exposure,
      all_of(covs),
      method = method_call(),
      estimand = estimand,
      constraints = constraints
    )
  }
  fit <- do_fit()
  stopifnot(isTRUE(fit@converged))
  w <- as.numeric(weights(fit))
  stopifnot(all(is.finite(w)))
  mk <- bench::mark(
    fit = as.numeric(weights(do_fit())),
    check = FALSE,
    min_iterations = iterations,
    max_iterations = iterations,
    filter_gc = FALSE
  )
  list(
    label = label,
    median = as.numeric(mk$median),
    mark = mk,
    solver_status = fit@solver_status
  )
}

run_internal_grid <- function() {
  options(balancing.threads = 1L)
  out <- list()

  # Estimating-equation family on binary data, n = 1000 and 10000, p = 10.
  for (n in c(1000L, 10000L)) {
    d <- make_binary(n, 10L)
    covs <- paste0("x", 1:10)
    iters <- if (n >= 10000L) 5L else 10L
    out[[sprintf("entropy_binary_ate_n%d", n)]] <-
      bench_esteq_method("entropy", bw_entropy, d, covs, "ate", iters)
    out[[sprintf("ipt_binary_ate_n%d", n)]] <-
      bench_esteq_method("ipt", bw_ipt, d, covs, "ate", iters)
    # The covariate balancing propensity score satisfies score-weighted moment
    # conditions rather than exact mean balance, so it is checked by convergence
    # only, not by the exact first-moment precision the tilt methods reach.
    out[[sprintf("cbps_binary_ate_n%d", n)]] <-
      bench_esteq_method(
        "cbps",
        bw_cbps,
        d,
        covs,
        "ate",
        iters,
        check_precision = FALSE
      )
  }

  # Categorical entropy and ipt at n = 10000, p = 10 (no exact-precision check;
  # multi-group balance is asserted through convergence).
  dc <- make_categorical(10000L, 10L)
  covsc <- paste0("x", 1:10)
  out[["entropy_categorical_ate_n10000"]] <-
    bench_esteq_method(
      "entropy",
      bw_entropy,
      dc,
      covsc,
      "ate",
      5L,
      check_precision = FALSE
    )
  out[["ipt_categorical_ate_n10000"]] <-
    bench_esteq_method(
      "ipt",
      bw_ipt,
      dc,
      covsc,
      "ate",
      5L,
      check_precision = FALSE
    )

  # Continuous entropy at n = 10000, p = 10.
  dk <- make_continuous(10000L, 10L)
  out[["entropy_continuous_ate_n10000"]] <-
    bench_esteq_method(
      "entropy",
      bw_entropy,
      dk,
      covsc,
      "ate",
      5L,
      check_precision = FALSE
    )

  # Energy balancing, binary ate, n = 500 and 2000, p = 4 (matches the energy
  # baseline workload).
  for (n in c(500L, 2000L)) {
    de <- make_binary(n, 4L, seed = 20240713L, coefs = c(0.4, -0.3, 0.2, -0.2))
    iters <- if (n >= 2000L) 3L else 5L
    out[[sprintf("energy_binary_ate_n%d", n)]] <-
      bench_qp_method(
        "energy",
        bw_energy,
        de$df,
        paste0("x", 1:4),
        "ate",
        iterations = iters
      )
  }

  # Stable balancing weights, binary ate at tol 0.05, n = 500, 2000, 5000, p = 4
  # (matches the sbw baseline workload).
  for (n in c(500L, 2000L, 5000L)) {
    ds <- make_binary(n, 4L, seed = 20240714L, coefs = c(0.4, -0.3, 0.2, -0.2))
    iters <- if (n >= 5000L) 5L else 10L
    out[[sprintf("sbw_binary_ate_n%d", n)]] <-
      bench_qp_method(
        "sbw",
        bw_sbw,
        ds$df,
        paste0("x", 1:4),
        "ate",
        constraints = balance_terms(tolerance = 0.05),
        iterations = iters
      )
  }

  out
}

# ---- Entropy external reference gate ---------------------------------------

# The performance gate: balance(bw_entropy()) must beat the equivalent
# WeightIt ebal call by at least 5x end to end at n = 50000 with 200 first
# moment constraints, at equal balance precision, judged as-shipped versus
# as-shipped. As-shipped, our method resolves its worker count automatically;
# WeightIt is single threaded. A single-threaded row is also recorded so the
# like-for-like solver margin is visible.
run_entropy_gate <- function(n = 50000L, p = 200L, iterations = 3L) {
  if (!requireNamespace("WeightIt", quietly = TRUE)) {
    message("WeightIt not installed; skipping the entropy reference gate")
    return(NULL)
  }
  d <- make_binary(n, p, seed = 20240711L)
  df <- d$df
  X <- d$X
  tr <- d$tr
  covs <- paste0("x", seq_len(p))
  f <- stats::as.formula(paste("exposure ~", paste(covs, collapse = " + ")))

  fit_ours_default <- function() {
    as.numeric(weights(balance(
      df,
      exposure,
      all_of(covs),
      method = bw_entropy(),
      estimand = "ate"
    )))
  }
  fit_ours_1t <- function() {
    old <- getOption("balancing.threads")
    options(balancing.threads = 1L)
    on.exit(options(balancing.threads = old))
    as.numeric(weights(balance(
      df,
      exposure,
      all_of(covs),
      method = bw_entropy(),
      estimand = "ate"
    )))
  }
  fit_wi <- function() {
    WeightIt::weightit(f, data = df, method = "ebal", estimand = "ATE")$weights
  }

  results <- list()
  for (estimand in c("ate", "att")) {
    ours_default <- function() {
      as.numeric(weights(balance(
        df,
        exposure,
        all_of(covs),
        method = bw_entropy(),
        estimand = estimand
      )))
    }
    wi <- function() {
      WeightIt::weightit(
        f,
        data = df,
        method = "ebal",
        estimand = toupper(estimand)
      )$weights
    }
    # As-shipped: automatic threads for ours.
    options(balancing.threads = NULL)
    w_ours <- ours_default()
    w_wi <- wi()
    if (estimand == "ate") {
      prec_ours <- ate_precision_binary(w_ours, X, tr)
      prec_wi <- ate_precision_binary(w_wi, X, tr)
    } else {
      tgt <- colMeans(X[tr == 1, ])
      m0o <- colSums(w_ours[tr == 0] * X[tr == 0, ]) / sum(w_ours[tr == 0])
      m0w <- colSums(w_wi[tr == 0] * X[tr == 0, ]) / sum(w_wi[tr == 0])
      prec_ours <- max(abs(m0o - tgt))
      prec_wi <- max(abs(m0w - tgt))
    }
    stopifnot(prec_ours < 1e-6, prec_wi < 1e-6)

    options(balancing.threads = NULL)
    mk_default <- bench::mark(
      x = ours_default(),
      check = FALSE,
      min_iterations = iterations,
      max_iterations = iterations,
      filter_gc = FALSE
    )
    options(balancing.threads = 1L)
    mk_1t <- bench::mark(
      x = ours_default(),
      check = FALSE,
      min_iterations = iterations,
      max_iterations = iterations,
      filter_gc = FALSE
    )
    options(balancing.threads = NULL)
    mk_wi <- bench::mark(
      x = wi(),
      check = FALSE,
      min_iterations = iterations,
      max_iterations = iterations,
      filter_gc = FALSE
    )

    med_default <- as.numeric(mk_default$median)
    med_1t <- as.numeric(mk_1t$median)
    med_wi <- as.numeric(mk_wi$median)
    results[[estimand]] <- list(
      ours_default = med_default,
      ours_1thread = med_1t,
      weightit = med_wi,
      prec_ours = prec_ours,
      prec_weightit = prec_wi,
      ratio_asshipped = med_wi / med_default,
      ratio_1thread = med_wi / med_1t
    )
    cat(sprintf(
      "[entropy gate %s] ours_default=%.3fs ours_1t=%.3fs weightit=%.3fs | as-shipped %.2fx | 1-thread %.2fx\n",
      estimand,
      med_default,
      med_1t,
      med_wi,
      med_wi / med_default,
      med_wi / med_1t
    ))
  }
  results
}

# ---- cfd overhead versus energy and non-trivial kernel ---------------------

# On the same binary data, the energy kernel of characteristic function distance
# balancing reproduces energy balancing, so the two fits' weights should agree
# and the cfd path's overhead over bw_energy() is the difference in medians. A
# gaussian-kernel fit on the same data times a non-trivial kernel build.
run_cfd_overhead <- function(sizes = c(500L, 2000L), iterations = 3L) {
  options(balancing.threads = 1L)
  options(balancing.qp_backend = "auto")
  out <- list()
  for (n in sizes) {
    d <- make_binary(n, 4L, seed = 20240713L, coefs = c(0.4, -0.3, 0.2, -0.2))
    df <- d$df
    covs <- paste0("x", 1:4)
    w_energy <- as.numeric(weights(balance(
      df,
      exposure,
      all_of(covs),
      method = bw_energy(),
      estimand = "ate"
    )))
    w_cfd_energy <- as.numeric(weights(balance(
      df,
      exposure,
      all_of(covs),
      method = bw_cfd(kernel = "energy"),
      estimand = "ate"
    )))
    rel <- max(abs(w_cfd_energy - w_energy) / pmax(abs(w_energy), 1e-8))
    corr <- stats::cor(w_energy, w_cfd_energy)
    cat(sprintf(
      "[cfd overhead n=%d] energy-kernel vs bw_energy max_rel=%.2e corr=%.6f\n",
      n,
      rel,
      corr
    ))

    iters <- if (n >= 2000L) 3L else 5L
    mk_energy <- bench::mark(
      x = as.numeric(weights(balance(
        df,
        exposure,
        all_of(covs),
        method = bw_energy(),
        estimand = "ate"
      ))),
      check = FALSE,
      min_iterations = iters,
      max_iterations = iters,
      filter_gc = FALSE
    )
    mk_cfd_e <- bench::mark(
      x = as.numeric(weights(balance(
        df,
        exposure,
        all_of(covs),
        method = bw_cfd(kernel = "energy"),
        estimand = "ate"
      ))),
      check = FALSE,
      min_iterations = iters,
      max_iterations = iters,
      filter_gc = FALSE
    )
    mk_cfd_g <- bench::mark(
      x = as.numeric(weights(balance(
        df,
        exposure,
        all_of(covs),
        method = bw_cfd(kernel = "gaussian"),
        estimand = "ate"
      ))),
      check = FALSE,
      min_iterations = iters,
      max_iterations = iters,
      filter_gc = FALSE
    )
    out[[as.character(n)]] <- list(
      energy = as.numeric(mk_energy$median),
      cfd_energy = as.numeric(mk_cfd_e$median),
      cfd_gaussian = as.numeric(mk_cfd_g$median),
      max_rel_weight_diff = rel,
      correlation = corr
    )
    cat(sprintf(
      "  medians: bw_energy=%.4fs cfd_energy=%.4fs cfd_gaussian=%.4fs\n",
      as.numeric(mk_energy$median),
      as.numeric(mk_cfd_e$median),
      as.numeric(mk_cfd_g$median)
    ))
  }
  out
}

# ---- cfd QP backend comparison at matched tolerance ------------------------

# The gaussian kernel is positive semidefinite, so both backends accept the
# spec. Fit the same problem under each backend at the shipped tolerance,
# cross-check that the two solutions agree, then compare medians. This tests the
# design's osqp default for the kernel method against clarabel with current code.
run_cfd_backends <- function(sizes = c(500L, 1000L, 2000L), iterations = 5L) {
  options(balancing.threads = 1L)
  out <- list()
  for (n in sizes) {
    d <- make_binary(n, 4L, seed = 20240713L, coefs = c(0.4, -0.3, 0.2, -0.2))
    df <- d$df
    covs <- paste0("x", 1:4)
    fit_be <- function(backend) {
      old <- getOption("balancing.qp_backend")
      options(balancing.qp_backend = backend)
      on.exit(options(balancing.qp_backend = old))
      balance(
        df,
        exposure,
        all_of(covs),
        method = bw_cfd(kernel = "gaussian"),
        estimand = "ate"
      )
    }
    fo <- fit_be("osqp")
    fc <- fit_be("clarabel")
    wo <- as.numeric(weights(fo))
    wc <- as.numeric(weights(fc))
    stopifnot(isTRUE(fo@converged), isTRUE(fc@converged))
    rel <- max(abs(wo - wc) / pmax(abs(wo), 1e-8))
    obj_gap <- abs(fo@objective - fc@objective) /
      max(abs(fo@objective), 1e-8)
    cat(sprintf(
      "[cfd backend n=%d] osqp vs clarabel weight_rel=%.2e obj_gap=%.2e (osqp obj=%.6e)\n",
      n,
      rel,
      obj_gap,
      fo@objective
    ))
    stopifnot(obj_gap < 1e-4)

    iters <- if (n >= 2000L) 3L else iterations
    mk_o <- bench::mark(
      x = {
        options(balancing.qp_backend = "osqp")
        as.numeric(weights(balance(
          df,
          exposure,
          all_of(covs),
          method = bw_cfd(kernel = "gaussian"),
          estimand = "ate"
        )))
      },
      check = FALSE,
      min_iterations = iters,
      max_iterations = iters,
      filter_gc = FALSE
    )
    mk_c <- bench::mark(
      x = {
        options(balancing.qp_backend = "clarabel")
        as.numeric(weights(balance(
          df,
          exposure,
          all_of(covs),
          method = bw_cfd(kernel = "gaussian"),
          estimand = "ate"
        )))
      },
      check = FALSE,
      min_iterations = iters,
      max_iterations = iters,
      filter_gc = FALSE
    )
    options(balancing.qp_backend = "auto")
    mo <- as.numeric(mk_o$median)
    mc <- as.numeric(mk_c$median)
    out[[as.character(n)]] <- list(
      osqp = mo,
      clarabel = mc,
      ratio_clarabel_over_osqp = mc / mo,
      weight_rel = rel,
      obj_gap = obj_gap,
      osqp_iters = as.integer(fo@iterations),
      clarabel_iters = as.integer(fc@iterations)
    )
    cat(sprintf(
      "  medians: osqp=%.5fs clarabel=%.5fs clarabel/osqp=%.2fx\n",
      mo,
      mc,
      mc / mo
    ))
  }
  out
}

# ---- sbw certificate fallback spot-check -----------------------------------

# The stable-balancing-weights binary ate instance at n = 20000 is the shape on
# which osqp falsely certifies primal infeasibility; the shipped default routing
# re-solves with clarabel on that certificate. Confirm the default fit resolves
# and records the fallback, and that pinning osqp errors while pinning clarabel
# solves.
run_sbw_fallback <- function(n = 20000L) {
  options(balancing.threads = 1L)
  d <- make_binary(n, 4L, seed = 20240714L, coefs = c(0.4, -0.3, 0.2, -0.2))
  df <- d$df
  covs <- paste0("x", 1:4)
  ctrl <- balance_terms(tolerance = 0.01)

  options(balancing.qp_backend = "auto")
  fit_auto <- tryCatch(
    balance(
      df,
      exposure,
      all_of(covs),
      method = bw_sbw(),
      estimand = "ate",
      constraints = ctrl
    ),
    error = function(e) e
  )
  auto_ok <- !inherits(fit_auto, "error") && isTRUE(fit_auto@converged)
  auto_status <- if (inherits(fit_auto, "error")) {
    paste("error:", conditionMessage(fit_auto))
  } else {
    fit_auto@solver_status
  }

  options(balancing.qp_backend = "osqp")
  fit_osqp <- tryCatch(
    balance(
      df,
      exposure,
      all_of(covs),
      method = bw_sbw(),
      estimand = "ate",
      constraints = ctrl
    ),
    error = function(e) e
  )
  osqp_errored <- inherits(fit_osqp, "error") || !isTRUE(fit_osqp@converged)

  options(balancing.qp_backend = "clarabel")
  fit_clar <- tryCatch(
    balance(
      df,
      exposure,
      all_of(covs),
      method = bw_sbw(),
      estimand = "ate",
      constraints = ctrl
    ),
    error = function(e) e
  )
  clar_ok <- !inherits(fit_clar, "error") && isTRUE(fit_clar@converged)

  options(balancing.qp_backend = "auto")
  cat(sprintf(
    "[sbw fallback n=%d] auto_ok=%s status=%s | osqp_errored=%s | clarabel_ok=%s\n",
    n,
    auto_ok,
    auto_status,
    osqp_errored,
    clar_ok
  ))
  list(
    n = n,
    auto_ok = auto_ok,
    auto_status = auto_status,
    osqp_errored = osqp_errored,
    clarabel_ok = clar_ok
  )
}

# ---- Driver ----------------------------------------------------------------

if (sys.nframe() == 0L) {
  dir.create(RAW_DIR, showWarnings = FALSE, recursive = TRUE)

  cat("== internal per-method grid ==\n")
  internal <- run_internal_grid()
  for (nm in names(internal)) {
    cat(sprintf("  %-34s median=%.5fs\n", nm, internal[[nm]]$median))
  }

  cat("\n== cfd overhead versus energy ==\n")
  cfd_overhead <- run_cfd_overhead()

  cat("\n== cfd QP backend comparison ==\n")
  cfd_backends <- run_cfd_backends()

  cat("\n== sbw certificate fallback ==\n")
  sbw_fallback <- run_sbw_fallback()

  cat("\n== entropy external reference gate ==\n")
  entropy_gate <- run_entropy_gate()

  saveRDS(
    list(
      internal = internal,
      cfd_overhead = cfd_overhead,
      cfd_backends = cfd_backends,
      sbw_fallback = sbw_fallback,
      entropy_gate = entropy_gate
    ),
    file.path(RAW_DIR, "consolidated.rds")
  )
  cat("\nsaved: ", file.path(RAW_DIR, "consolidated.rds"), "\n", sep = "")
}
