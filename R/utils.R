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

# Sampling-weighted mean of a numeric vector.
weighted_center <- function(x, w) {
  sum(w * x) / sum(w)
}

# Sampling-weighted standard deviation with the reliability (frequency-free)
# denominator `sum(w) - sum(w^2) / sum(w)`. This matches the core's covariate
# transform and reduces exactly to `stats::sd()` when the weights are equal, so a
# uniform-weight fit is standardized identically to the unweighted convention. A
# column with no spread reports zero; callers substitute one to leave it
# unscaled.
weighted_scale <- function(x, w) {
  sw <- sum(w)
  mean <- sum(w * x) / sw
  denom <- sw - sum(w * w) / sw
  variance <- if (denom > 0) sum(w * (x - mean)^2) / denom else 0
  sqrt(max(variance, 0))
}

# ---- The reported weight scale ---------------------------------------------

# The sampling-weighted total each exposure group is reported at. A focal
# estimand reports every group at the focal group's total, which puts the focal
# group at its own base weights carried to that total; those coincide with the
# base weights themselves only when the base weights already sum to it, as
# uniform base weights do. The pooled estimands report each group at its own
# total. A focal level is resolved for exactly the focal estimands, so it carries
# the distinction on its own.
group_target_sums <- function(s, groups, focal_level = NULL) {
  targets <- vapply(groups, function(idx) sum(s[idx]), numeric(1))
  if (!is.null(focal_level)) {
    targets[] <- targets[[focal_level]]
  }
  targets
}

# Move each exposure group's weights onto its reported total. The solvers
# normalize on their own internal convention and the reported convention places
# each group's sampling-weighted total at `targets`, so the correction is one
# constant per group: it changes the reporting scale without disturbing the
# balance the solve achieved. A group whose weights already sum to zero has no
# scale to move to and is left alone. The arguments are the weights together
# with the groups and targets rather than a fitted object, so a caller
# re-evaluating the weights at other parameters can apply the same convention.
#
# A solve that diverged returns weights that are not finite, which leaves the
# group total this divides by as a missing value. That is a failed solve rather
# than a reporting-scale question, so it is refused here with a classed error
# instead of steering the comparison below with a missing value.
renormalize_group_weights <- function(
  w,
  s,
  groups,
  targets,
  call = rlang::caller_env()
) {
  for (level in names(groups)) {
    idx <- groups[[level]]
    current <- sum(s[idx] * w[idx])
    if (!is.finite(current)) {
      abort(
        c(
          "The solver did not produce finite weights.",
          x = "The weights for exposure level {.val {level}} do not sum to a finite total.",
          i = "Check the covariates for collinearity or for a column the exposure determines."
        ),
        error_class = "balancing_convergence_error",
        call = call
      )
    }
    if (current > 0) {
      w[idx] <- w[idx] * (targets[[level]] / current)
    }
  }
  w
}

# Build the container's weight re-evaluation hook from an evaluator returning a
# method's weights on the scale the container stores them at. The reported scale
# differs from the stored scale by one constant per group, and that constant is
# fixed here, at the fit, rather than recomputed at every set of parameters. The
# hook is then a fixed rescaling of the solver's weight path, whose derivative
# is exactly `weight_jacobian` moved to the reported scale. Recomputing the
# constant instead would add the derivative of the group totals, a different
# coupling from the one the container's Jacobian describes.
make_weights_fn <- function(eval_fn, reported, weights_raw) {
  force(eval_fn)
  scale <- reported / weights_raw
  # A unit the solve holds at zero weight has no reporting scale to move to, so
  # it keeps whatever the evaluator returns rather than an undefined ratio.
  scale[!is.finite(scale)] <- 1
  function(theta) {
    as.numeric(eval_fn(theta)) * scale
  }
}

# Guard the trailing dots of a method constructor. Each constructor takes only
# its named tuning parameters, so an unexpected argument, usually a misspelled
# name, raises a classed `balancing_method_error` naming the offending arguments
# rather than the unclassed dots error `rlang::check_dots_empty()` produces.
check_method_dots <- function(..., call = rlang::caller_env()) {
  if (...length() == 0L) {
    return(invisible())
  }
  names <- rlang::names2(rlang::list2(...))
  named <- names[nzchar(names)]
  n_unnamed <- ...length() - length(named)
  if (length(named) > 0) {
    abort(
      c(
        "Unknown tuning argument{?s} {.arg {named}}.",
        i = "Check the argument names against the constructor's help page."
      ),
      error_class = "balancing_method_error",
      call = call
    )
  }
  abort(
    c(
      "This method constructor accepts only named tuning arguments.",
      x = "You passed {n_unnamed} unnamed argument{?s}."
    ),
    error_class = "balancing_method_error",
    call = call
  )
}

# Guard the tolerance a balance specification is constructed with. The property
# validator refuses a missing, infinite, or negative tolerance as well, but an S7
# validator answers with a bare string, so those refusals reach a caller
# unclassed while every other refusal a balance specification raises carries
# `balancing_constraints_error`. Checking the argument here classes them, and
# covers the one shape the validator never sees: `NULL`, which the property type
# check turns away first and reports as a class mismatch rather than as a
# tolerance with no value. The validator keeps its own clauses, since a property
# may also be assigned after construction, and an assignment reaches nothing
# else.
check_tolerance <- function(tolerance, call = rlang::caller_env()) {
  if (is.null(tolerance)) {
    abort(
      c(
        "{.arg tolerance} must be a number.",
        x = "It is {.code NULL}.",
        i = "Use {.code 0} to request exact balance."
      ),
      error_class = "balancing_constraints_error",
      call = call
    )
  }
  if (anyNA(tolerance)) {
    abort(
      c(
        "{.arg tolerance} must not contain missing values.",
        i = "Give every covariate you name a non-negative number."
      ),
      error_class = "balancing_constraints_error",
      call = call
    )
  }
  if (!all(is.finite(tolerance))) {
    abort(
      c(
        "{.arg tolerance} must be finite.",
        x = "An infinite tolerance leaves its constraint unbounded, which asks for no balance at all.",
        i = "Relax a constraint with a large finite tolerance, or drop the covariate."
      ),
      error_class = "balancing_constraints_error",
      call = call
    )
  }
  if (any(tolerance < 0)) {
    abort(
      c(
        "{.arg tolerance} must be non-negative.",
        x = "A tolerance is the largest imbalance permitted, which cannot be less than none.",
        i = "Use {.code 0} to request exact balance."
      ),
      error_class = "balancing_constraints_error",
      call = call
    )
  }
  invisible(tolerance)
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

# Resolve the solver for the exact entropy problem. The default is Newton, the
# only solver that drives the estimating equations to machine precision. The
# default is read from an option so it can change in one place without touching
# the fit path; the alternatives are the basin L-BFGS adapter and the
# L-BFGS-then-Newton hybrid, whose Newton polish restores machine-precision
# estimating equations.
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
# The default is "auto": the default solver runs first and, on a
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
