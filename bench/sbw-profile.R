# Profile one stable balancing fit and attribute the time.
#
# Profiles balance(method = sbw()) on a binary average-treatment problem at
# n = 5000, the size where the quadratic-program solve dominates but the fit
# still runs in well under a second so the profiler collects enough samples. The
# raw profvis object is saved to scratch/bench/raw and summarized with the
# debrief package: the hot lines, the hot call paths, and the R-versus-native
# split, so the report can attribute time to R-side covariate processing versus
# the native solver.
#
# This script is Rbuildignored. Run it against an installed release build.

suppressMessages({
  library(balancing)
  library(profvis)
  library(debrief)
})

options(balancing.threads = 1L)

make_data <- function(n, p, seed = 20240714L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  lp <- 0.4 * X[, 1] - 0.3 * X[, 2] + 0.2 * X[, 3] - 0.2 * X[, 4]
  exposure <- stats::rbinom(n, 1, stats::plogis(lp))
  data.frame(exposure = exposure, X)
}

n <- 5000L
p <- 4L
df <- make_data(n, p)
covs <- paste0("x", seq_len(p))
tol <- 0.05

fit_once <- function() {
  balance(
    df,
    exposure,
    all_of(covs),
    method = sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = tol)
  )
}

# Warm the byte compiler and confirm the fit converges before profiling.
invisible(fit_once())

# Repeat the fit inside the profiler so the solve accumulates enough samples at
# the default 10 ms interval.
pv <- profvis::profvis(
  {
    for (i in seq_len(40)) {
      fit_once()
    }
  },
  interval = 0.005
)

dir.create(
  file.path("scratch", "bench", "raw"),
  showWarnings = FALSE,
  recursive = TRUE
)
saveRDS(pv, file.path("scratch", "bench", "raw", "sbw-profvis-n5000.rds"))

cat("\n===== debrief: full debrief =====\n")
print(debrief::pv_print_debrief(pv))
cat("\n===== debrief: hot lines =====\n")
print(debrief::pv_print_hot_lines(pv))
cat("\n===== debrief: hot paths =====\n")
print(debrief::pv_print_hot_paths(pv))
cat("\n===== debrief: suggestions =====\n")
print(debrief::pv_print_suggestions(pv))
