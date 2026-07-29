# balance() is the orchestrator. It resolves the exposure and covariate
# selections, detects and announces the exposure type, validates the estimand and
# focal level against the method's support, evaluates the sampling weights, builds
# the constraint matrix, dispatches to the method's fit, and assembles the
# balancing result together with its balance table.

#' Estimate balancing weights
#'
#' `balance()` is the entry point to the package. It fits a balancing method to a
#' data frame, returning weights that target covariate balance directly. The
#' exposure and covariates are chosen with tidyselect, the method is one of the
#' method specifications such as [bw_entropy()], and the estimand and
#' constraints control what balance the weights achieve.
#'
#' @details
#' The exposure type is detected automatically and announced through an
#' informational message, which `options(balancing.quiet = TRUE)` suppresses. The
#' estimand vocabulary matches propensity: `"atc"` is accepted as a synonym for
#' the untreated target and stored as `"atu"`. `"att"` and `"atc"` reweight
#' toward a focal exposure level, inferred for a binary exposure and required
#' through `focal_level` for a categorical exposure. Continuous exposures permit
#' only `"ate"`.
#'
#' Constraints default to first-moment balance. Pass a [balance_terms()]
#' specification to balance higher moments, interactions, or quantiles, or to
#' relax exact balance to a tolerance.
#'
#' @param .data A data frame.
#' @param .exposure The exposure column, selected with data-masking. Exactly one
#'   column.
#' @param .covariates The covariate columns, selected with tidyselect. At least
#'   one column, with no default.
#' @param method A [balance_method] specification from one of the method
#'   constructors, such as [bw_entropy()].
#' @param estimand The target estimand: `"ate"`, `"att"`, `"atc"` (stored as
#'   `"atu"`), or `"ato"`. Defaults to `"ate"`.
#' @param ... Reserved; must be empty.
#' @param constraints A [balance_terms()] specification, or `NULL` for the method
#'   default.
#' @param exposure_type One of `"auto"` (the default), `"binary"`,
#'   `"categorical"`, or `"continuous"`.
#' @param focal_level The focal exposure level for `"att"` and `"atc"`. Inferred
#'   for a binary exposure; required for a categorical exposure.
#' @param sampling_weights Sampling weights, given as a bare column name or an
#'   external numeric vector, or `NULL`.
#'
#' @return A [balancing] object.
#'
#' @examples
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' df <- data.frame(
#'   exposure = rbinom(n, 1, plogis(0.5 * x1 - 0.5 * x2)),
#'   x1 = x1,
#'   x2 = x2
#' )
#' fit <- balance(df, exposure, c(x1, x2), method = bw_entropy())
#' fit
#' weights(fit)
#'
#' @export
balance <- function(
  .data,
  .exposure,
  .covariates,
  method = bw_entropy(),
  estimand = c("ate", "att", "atc", "ato"),
  ...,
  constraints = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  focal_level = NULL,
  sampling_weights = NULL
) {
  the_call <- match.call()
  rlang::check_dots_empty()

  validate_data_frame(.data)
  validate_nonempty(.data)

  if (!S7::S7_inherits(method, balance_method)) {
    abort(
      c(
        "{.arg method} must be a balancing method specification.",
        x = "You supplied {.obj_type_friendly {method}}.",
        i = "Construct one with a method constructor, for example {.code bw_entropy()}."
      ),
      error_class = "balancing_method_error"
    )
  }

  exposure_pos <- tidyselect::eval_select(rlang::enquo(.exposure), .data)
  validate_selection(exposure_pos, ".exposure", expected = "one")
  exposure_name <- names(exposure_pos)
  exposure_vec <- .data[[exposure_pos]]

  covariate_pos <- tidyselect::eval_select(rlang::enquo(.covariates), .data)
  validate_selection(covariate_pos, ".covariates", expected = "some")
  covariate_names <- names(covariate_pos)

  n <- nrow(.data)

  sampling_weights_value <- rlang::eval_tidy(
    rlang::enquo(sampling_weights),
    data = .data
  )
  if (!is.null(sampling_weights_value)) {
    validate_sampling_weights(sampling_weights_value, n)
  }

  validate_no_missing(exposure_vec, .data, covariate_names)

  exposure_type <- resolve_exposure_type(exposure_type, exposure_vec, method)

  estimand <- rlang::arg_match0(
    estimand[[1]],
    c("ate", "att", "atc", "atu", "ato"),
    arg_nm = "estimand"
  )
  estimand <- canonical_estimand(estimand)
  supported <- supported_estimands(method, exposure_type)
  if (!estimand %in% supported) {
    abort(
      c(
        "{method_label(method)} does not support the {.val {estimand}} estimand for a {.val {exposure_type}} exposure.",
        i = "Supported estimands are {.val {supported}}."
      ),
      error_class = "balancing_estimand_error"
    )
  }

  exposure_key <- as.character(exposure_vec)
  levels <- exposure_levels(exposure_vec, exposure_type)
  focal_level <- resolve_focal_level(
    estimand,
    exposure_type,
    levels,
    focal_level
  )

  constraints <- constraints %||% default_constraints(method)

  # The quadratic-program family enforces its tolerance box on the standardized
  # columns, so those columns cross the boundary on the sampling-weighted scale the
  # reference measures the box on. The estimating-equation family keeps the
  # unweighted scale, which conditions the Newton step better: its exact problem is
  # invariant to an affine change of the columns, and its inexact problem's box is
  # rescaled to the weighted scale in `solver_box()` regardless of the column
  # scale, so both are unaffected. Either way the balance table reports on the
  # weighted scale.
  constraint_sampling_weights <- if (
    S7::S7_inherits(method, quadratic_program_method)
  ) {
    sampling_weights_value
  } else {
    NULL
  }

  built <- build_constraint_matrix(
    .data,
    covariate_names,
    constraints,
    exposure_type,
    sampling_weights = constraint_sampling_weights
  )

  groups <- if (identical(exposure_type, "continuous")) {
    NULL
  } else {
    stats::setNames(
      lapply(levels, function(level) which(exposure_key == level)),
      levels
    )
  }

  prepared <- list(
    matrix = built$matrix,
    recipe = built$recipe,
    data = .data,
    covariates = covariate_names,
    exposure_vec = exposure_vec,
    exposure_key = exposure_key,
    exposure_type = exposure_type,
    exposure_levels = levels,
    groups = groups,
    estimand = estimand,
    focal_level = focal_level,
    sampling_weights = sampling_weights_value %||% rep(1, n),
    n = n,
    constraints = constraints,
    tolerances = column_tolerances(built$recipe)
  )

  fit <- fit_method(method, prepared)

  # A user interrupt during the solve unwinds the Rust core cleanly and reports
  # itself as a flag rather than a longjmp. Re-signal it here so a cancelled fit
  # raises the interrupt condition instead of returning a partial result behind
  # a convergence warning.
  if (isTRUE(fit$interrupted)) {
    rlang::interrupt()
  }

  check_solver_status(fit, method)

  weights <- new_bw(fit$weights, estimand = estimand, groups = groups)

  # The base measure the solver targets. Only entropy balancing carries base
  # weights; other methods anchor to a uniform measure. It sets the pooled target
  # the average-treatment-effect constraint geometry measures against.
  base_weights <- if ("base_weights" %in% S7::prop_names(method)) {
    method@base_weights %||% rep(1, n)
  } else {
    rep(1, n)
  }
  base_measure <- prepared$sampling_weights * base_weights

  balance_table <- compute_balance_table(
    built$recipe,
    .data,
    exposure_vec,
    exposure_type,
    estimand,
    focal_level,
    groups,
    as.numeric(weights) * prepared$sampling_weights,
    tolerance = 0,
    reference = base_measure,
    constraint_target = fit$constraint_target %||% "pooled",
    sampling_weights = sampling_weights_value
  )

  # A fit warns when a constraint sits outside its tolerance box, judged on the
  # solver's arm-to-target geometry through `within_tolerance`. Methods whose
  # balance is approximate by construction, such as the over-identified covariate
  # balancing propensity score, never consume a tolerance, so they report their
  # criterion rather than warning against a knob that does not reach the fit. A
  # verdict that resolves to anything other than a plain TRUE, which a constraint
  # column with no spread would produce, is reported as outside the box rather
  # than left to steer the branch as a missing value.
  if (
    !isTRUE(fit$approximate) && !all(balance_table$within_tolerance %in% TRUE)
  ) {
    worst <- max(abs(balance_table$weighted))
    warn(
      c(
        "The achieved balance exceeds the requested tolerance.",
        x = "The largest imbalance is {formatC(worst, format = 'f', digits = 4)}.",
        i = "Raise {.arg tolerance} in {.fn balance_terms}, lower the moments, or drop interactions."
      ),
      warning_class = "balancing_balance_warning"
    )
  }

  balancing(
    weights = weights,
    method = method,
    estimand = estimand,
    exposure = exposure_name,
    exposure_type = exposure_type,
    covariates = covariate_names,
    focal_level = focal_level,
    n = as.integer(n),
    constraints = constraints,
    recipe = built$recipe,
    balance_table = balance_table,
    duals = fit$duals,
    coefficients = fit$coefficients,
    converged = fit$converged,
    iterations = fit$iterations,
    objective = fit$objective,
    solver_status = fit$solver_status,
    estimating_equations = fit$estimating_equations,
    sampling_weights = sampling_weights_value,
    call = the_call
  )
}

# Raise or warn on the solver outcome. The quadratic-program family reports a
# terminal status the backend assigns, so an infeasible constraint set raises
# `balancing_infeasible_error` and a hard solver failure raises
# `balancing_convergence_error`, each naming the knob to turn; a reached iteration
# cap with a usable iterate still warns. The estimating-equation family carries no
# status and warns when it did not meet its convergence tolerance.
check_solver_status <- function(fit, method, call = rlang::caller_env()) {
  if (!is.null(fit$status)) {
    if (isTRUE(fit$converged)) {
      return(invisible())
    }
    if (identical(fit$status, "primal_infeasible")) {
      abort(
        c(
          "The balancing problem is infeasible.",
          x = "The solver reported that the constraints cannot be satisfied together.",
          i = "Raise {.arg tolerance} in {.fn balance_terms}, lower the moments, or drop interactions."
        ),
        error_class = "balancing_infeasible_error",
        call = call
      )
    }
    if (fit$status %in% c("non_convex", "dual_infeasible")) {
      abort(
        c(
          "The solver failed to produce a valid solution.",
          x = "It terminated with status {.val {fit$status}}.",
          i = "Check the covariates for collinearity, or raise {.arg weight_penalty} in {.fn {class(method)[1]}}."
        ),
        error_class = "balancing_convergence_error",
        call = call
      )
    }
  }
  if (!isTRUE(fit$converged)) {
    warn(
      c(
        "The solver did not reach its convergence tolerance.",
        i = "Increase {.arg max_iterations} or loosen {.arg convergence_tolerance} in {.fn {class(method)[1]}}."
      ),
      warning_class = "balancing_convergence_warning",
      call = call
    )
  }
  invisible()
}

# The default constraint set for a method with no explicit constraints, dispatched
# on the method so each family states its own default. The estimating-equation
# family balances first moments, its identifying conditions; the objective-driven
# quadratic-program methods let their objective drive balance and add no moment
# constraints unless the caller requests them; stable balancing weights override
# this with first-moment balance, since they need explicit constraints.
default_constraints <- new_generic("default_constraints", "method")

method(default_constraints, balance_method) <- function(method) {
  balance_terms(moments = 1L)
}

method(default_constraints, quadratic_program_method) <- function(method) {
  balance_terms()
}

# The per-column tolerances carried by a recipe, aligned with the constraint
# matrix columns. Empty recipes yield a zero-length vector.
column_tolerances <- function(recipe) {
  vapply(recipe, function(record) record$tolerance, numeric(1))
}

# The distinct exposure levels in a stable order: the factor levels when the
# exposure is a factor, otherwise the unique values sorted in the vector's own
# type and then converted to strings. Sorting the string forms instead would put
# "10" before "9", so a two-level numeric dose of 9 and 10 would name 9 as the
# second level and a focal estimand would reweight the wrong group. Ordering on
# the values reproduces base `factor()`, which is the order propensity and the
# rest of the ecosystem promise. A factor may carry levels no observation takes;
# those empty levels would form groups of size zero that misalign the
# estimating-equation design and yield an undefined effective sample size, so they
# are dropped with an informational alert and only the levels present in the data
# are balanced.
exposure_levels <- function(exposure_vec, exposure_type) {
  if (identical(exposure_type, "continuous")) {
    return(character(0))
  }
  observed <- exposure_vec[!is.na(exposure_vec)]
  if (is.factor(exposure_vec)) {
    present <- unique(as.character(observed))
    all_levels <- levels(exposure_vec)
    unused <- setdiff(all_levels, present)
    if (length(unused) > 0) {
      alert_info("Dropping unused exposure level{?s} {.val {unused}}.")
    }
    all_levels[all_levels %in% present]
  } else {
    as.character(sort(unique(observed)))
  }
}

# The estimand vocabulary accepts "atc" for the untreated target, matching
# propensity, and a fit stores the single canonical spelling "atu". Anything that
# compares a caller's estimand against a stored one canonicalizes first, so the
# spelling that created a fit is also a spelling that names it afterwards.
canonical_estimand <- function(estimand) {
  if (identical(estimand, "atc")) "atu" else estimand
}

# Resolve the focal exposure level for att and atc. A binary exposure infers the
# treated level (the second level) for att and the control level (the first) for
# atc; a categorical exposure requires an explicit focal_level. The average
# treatment effect and the overlap estimand reweight every group rather than hold
# one fixed, so they carry no focal level.
resolve_focal_level <- function(
  estimand,
  exposure_type,
  levels,
  focal_level,
  call = rlang::caller_env()
) {
  if (estimand %in% c("ate", "ato")) {
    return(NULL)
  }

  if (identical(exposure_type, "binary")) {
    if (!is.null(focal_level)) {
      resolved <- as.character(focal_level)
    } else if (identical(estimand, "att")) {
      resolved <- levels[[2]]
    } else {
      resolved <- levels[[1]]
    }
  } else {
    if (is.null(focal_level)) {
      abort(
        c(
          "{.arg focal_level} is required for the {.val {estimand}} estimand with a categorical exposure.",
          i = "Supply the exposure level to target, one of {.val {levels}}."
        ),
        error_class = "balancing_estimand_error",
        call = call
      )
    }
    resolved <- as.character(focal_level)
  }

  if (!resolved %in% levels) {
    abort(
      c(
        "{.arg focal_level} must be an exposure level.",
        x = "{.val {resolved}} is not one of {.val {levels}}."
      ),
      error_class = "balancing_estimand_error",
      call = call
    )
  }

  resolved
}
