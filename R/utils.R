be_quiet <- function() {
  getOption("balancing.quiet", default = FALSE)
}

abort <- function(
  ...,
  error_class = NULL,
  call = rlang::caller_env(),
  .envir = parent.frame()
) {
  cli::cli_abort(
    ...,
    class = c(error_class, "balancing_error"),
    call = call,
    .envir = .envir
  )
}

warn <- function(
  ...,
  warning_class = NULL,
  call = rlang::caller_env(),
  .envir = parent.frame()
) {
  cli::cli_warn(
    ...,
    class = c(warning_class, "balancing_warning"),
    call = call,
    .envir = .envir
  )
}

alert_info <- function(.message, .envir = parent.frame()) {
  if (!be_quiet()) {
    cli::cli_alert_info(text = .message, .envir = .envir)
  }
}

# Resolve the worker-thread count for a solver call. The resolution order is the
# explicit `threads` argument, then the `balancing.threads` option, then an
# automatic count. The automatic count is the physical core count, capped by
# OMP_THREAD_LIMIT and OMP_NUM_THREADS and forced to two under R CMD check, so a
# CRAN run never spawns more than the checker permits. The count is computed on
# the R side and passed to the core, which re-checks OMP_THREAD_LIMIT as a
# backstop.
resolve_threads <- function(threads = NULL) {
  if (!is.null(threads)) {
    return(max(1L, as.integer(threads)))
  }
  option <- getOption("balancing.threads")
  if (!is.null(option)) {
    return(max(1L, as.integer(option)))
  }
  automatic_threads()
}

automatic_threads <- function() {
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) {
    return(2L)
  }
  physical <- parallel::detectCores(logical = FALSE)
  if (is.na(physical) || physical < 1L) {
    physical <- 1L
  }
  caps <- c(
    physical,
    env_thread_cap("OMP_THREAD_LIMIT"),
    env_thread_cap("OMP_NUM_THREADS")
  )
  max(1L, as.integer(min(caps)))
}

# Resolve the solver for the exact entropy problem. The shipped default is
# Newton, the only solver that drives the estimating equations to machine
# precision. The default is read from an option so the benchmark promotion
# process can change it in one place without touching the fit path; the
# alternatives are the basin L-BFGS adapter and the L-BFGS-then-Newton hybrid,
# whose Newton polish restores machine-precision estimating equations.
resolve_entropy_solver <- function() {
  choices <- c("newton", "lbfgs", "lbfgs_then_newton")
  solver <- getOption("balancing.entropy_solver", default = "newton")
  if (!is.character(solver) || length(solver) != 1L || !(solver %in% choices)) {
    abort(
      c(
        "The {.code balancing.entropy_solver} option must be one of {.val {choices}}.",
        x = "It is {.val {solver}}."
      ),
      error_class = "balancing_range_error"
    )
  }
  solver
}

# Resolve the quadratic-program backend for the positive-semidefinite methods.
# The shipped default is "auto": the default solver runs first and, on a
# primal-infeasibility certificate, the fit re-solves with the interior-point
# backend, which handles feasible instances the default solver can falsely
# certify infeasible. The value is read from an option so a user can pin a
# backend without a constructor argument, matching the entropy-solver knob; the
# alternatives are "osqp" (no fallback) and "clarabel".
resolve_qp_backend <- function() {
  choices <- c("auto", "osqp", "clarabel")
  backend <- getOption("balancing.qp_backend", default = "auto")
  if (
    !is.character(backend) ||
      length(backend) != 1L ||
      !(backend %in% choices)
  ) {
    abort(
      c(
        "The {.code balancing.qp_backend} option must be one of {.val {choices}}.",
        x = "It is {.val {backend}}."
      ),
      error_class = "balancing_range_error"
    )
  }
  backend
}

# Parse an environment-variable thread cap, returning Inf when the variable is
# unset or not a positive whole number so it does not constrain the minimum.
env_thread_cap <- function(name) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) {
    return(Inf)
  }
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) {
    Inf
  } else {
    parsed
  }
}
