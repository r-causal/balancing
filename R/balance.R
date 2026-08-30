# balance() is the orchestrator. It resolves the exposure and covariate
# selections, detects and announces the exposure type, validates the estimand and
# focal level against the method's support, evaluates the sampling weights, builds
# the constraint matrix, dispatches to the method's fit, and assembles the
# balancing result together with its balance table.

#' Estimate balancing weights
#'
#' `balance()` fits a balancing method to a data frame, returning weights that
#' target covariate balance directly. The exposure and covariates are chosen with
#' tidyselect, the method is one of the method specifications such as
#' [bw_entropy()], and the estimand and constraints control what balance the
#' weights achieve.
#'
#' @details
#' The exposure type is detected automatically and announced through an
#' informational message, which `options(balancing.quiet = TRUE)` suppresses. The
#' estimand vocabulary matches propensity: `"atc"` is accepted as a synonym for
#' the untreated target and stored as `"atu"`. `"att"` and `"atc"` reweight
#' toward a focal exposure level, inferred for a binary exposure and required
#' through `.focal_level` for a categorical exposure. Continuous exposures permit
#' only `"ate"`.
#'
#' Constraints default to first-moment balance. Pass a [balance_terms()]
#' specification to balance higher moments, interactions, or quantiles, or to
#' relax exact balance to a tolerance.
#'
#' A factor covariate expands to one indicator per level, and those indicators
#' sum to the constant every balancing method carries. One indicator per factor
#' is therefore redundant with that constant and is dropped, with an
#' informational alert naming the term. The dropped level is the last one, and
#' balancing the levels that remain balances it too. The factor stays in
#' `@covariates`, and `@balance_table` reports the surviving levels rather than
#' the full set.
#'
#' A fit can be interrupted between solver iterations, so a long solve stops at
#' the next iteration rather than at the end of the fit. On Unix the poll reads
#' R's interrupt flag directly and does not service R's event loop, so a
#' [setTimeLimit()] set around the call fires when the call returns rather than
#' partway through the solve.
#'
#' A `difftime` covariate balances as the number it stores, in the unit its own
#' column declares. Nothing rescales it and nothing reinterprets the unit, so its
#' constraints, its recipe, and its balance table match those of the same
#' durations supplied as bare numbers. A `Date` or `POSIXt` covariate balances
#' the same way, as the number `as.numeric()` gives it: days since 1970-01-01 for
#' a date, seconds since then for a date-time. Both date-time representations
#' are read that way, so a `POSIXlt` column balances exactly as the `POSIXct`
#' column holding the same instants does.
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
#' @param .focal_level The focal exposure level for `"att"` and `"atc"`. Inferred
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
  .focal_level = NULL,
  sampling_weights = NULL
) {
  the_call <- match.call()
  rlang::check_dots_empty()

  validate_data_frame(.data)
  validate_row_count(.data)

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

  # A selection may not rename. A rename resolves the column by position but
  # supplies a new name, and every lookup after this indexes the data by the
  # selection's names: a covariate renamed to the exposure would build the
  # constraints on the exposure column, and a renamed exposure would name a column
  # the data do not have. Refusing the rename at the selection reports the real
  # defect where it was written, through tidyselect's own classed error.
  exposure_pos <- tidyselect::eval_select(
    rlang::enquo(.exposure),
    .data,
    allow_rename = FALSE
  )
  validate_selection(exposure_pos, ".exposure", expected = "one")
  exposure_name <- names(exposure_pos)
  exposure_vec <- .data[[exposure_pos]]

  covariate_pos <- tidyselect::eval_select(
    rlang::enquo(.covariates),
    .data,
    allow_rename = FALSE
  )
  covariate_pos <- drop_exposure_covariate(
    covariate_pos,
    exposure_pos,
    exposure_name
  )
  validate_selection(covariate_pos, ".covariates", expected = "some")
  covariate_names <- names(covariate_pos)

  n <- nrow(.data)

  sampling_weights_value <- rlang::eval_tidy(
    rlang::enquo(sampling_weights),
    data = .data
  )
  if (!is.null(sampling_weights_value)) {
    validate_sampling_weights(sampling_weights_value, n)
    # Sampling weights arrive as counts often enough that an integer column is a
    # natural way to supply them, and both the solver boundary and this fit's
    # `sampling_weights` property take a double. The validated vector is coerced
    # once here rather than at each solver call, which also drops any names the
    # source column carried.
    sampling_weights_value <- as.numeric(sampling_weights_value)
  }

  validate_finite_data(exposure_vec, .data, covariate_names)

  exposure_type <- resolve_exposure_type(exposure_type, exposure_vec, method)

  estimand <- rlang::arg_match0(
    estimand[[1]],
    estimand_choices(),
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
  validate_exposure_level_count(levels, exposure_type)
  focal_level <- resolve_focal_level(
    estimand,
    exposure_type,
    levels,
    .focal_level
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

  check_constraint_columns(method, built$matrix, covariate_names)

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

  # The sampling weights are the whole base measure for every method that carries
  # no base weights of its own, so a group left without mass is refused here,
  # before any fit, rather than through whichever failure each method's path
  # happens to reach. Entropy balancing checks the product its base weights form
  # as well, where their length has been validated.
  validate_base_measure(prepared$sampling_weights, groups)

  fit <- fit_method(method, prepared)

  # A user interrupt during the solve unwinds the Rust core cleanly and reports
  # itself as a flag rather than a longjmp. Re-signal it here so a cancelled fit
  # raises the interrupt condition instead of returning a partial result behind
  # a convergence warning.
  if (isTRUE(fit$interrupted)) {
    rlang::interrupt()
  }

  check_solver_status(fit, method)

  weights <- new_bw(fit$weights, estimand = estimand)

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
    sampling_weights = sampling_weights_value,
    matrix = built$matrix,
    enforced_tolerance = fit$enforced_tolerance
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
    warn_balance_exceeded(max(abs(balance_table$weighted)))
  }

  balancing(
    weights = weights,
    method = method,
    estimand = estimand,
    exposure = exposure_name,
    exposure_type = exposure_type,
    exposure_levels = levels,
    covariates = constrained_covariates(built$recipe, covariate_names),
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

# Report a fit whose achieved balance sits outside its tolerance box. The largest
# imbalance names how far the fit missed, so it is only offered when it is a
# number: a maximum that is not finite means at least one constraint's statistic is
# undefined, and "the largest imbalance is NaN" states a distance that was never
# measured. That case reports the assessment as the thing that failed, since no
# tolerance the caller could raise would answer it. Requiring every exposure level
# to carry base-measure mass removes the reachable cause, so this is the guard
# behind that rather than a case a fit reaches. The imbalance prints to significant
# digits rather than to a fixed number of decimals because a tolerance can sit far
# below the fourth decimal, and a fixed-decimal format would then round the value
# that triggered the warning down to a zero that contradicts it.
warn_balance_exceeded <- function(worst, call = rlang::caller_env()) {
  if (!is.finite(worst)) {
    warn(
      c(
        "The achieved balance could not be assessed.",
        x = "At least one constraint's balance statistic is undefined.",
        i = "Check for an exposure level whose weights sum to zero."
      ),
      warning_class = "balancing_balance_warning",
      call = call
    )
    return(invisible())
  }
  warn(
    c(
      "The achieved balance exceeds the requested tolerance.",
      x = "The largest imbalance is {formatC(worst, format = 'g', digits = 3)}.",
      i = "Raise {.arg tolerance} in {.fn balance_terms}, lower the moments, or drop interactions."
    ),
    warning_class = "balancing_balance_warning",
    call = call
  )
  invisible()
}

# The tolerance a quadratic-program method asks of its backend. A method that
# leaves the property NULL takes the core default, which the solver applies as
# both its absolute and its relative tolerance.
qp_default_tolerance <- 1e-8

# A tolerance the quadratic programs reach on problems where the core default
# does not. The energy objective matrix is indefinite, and on a small sample the
# negative curvature it carries puts the alternating-direction residual floor
# above the core default, so a fit asking for more than the iteration can deliver
# spends its whole budget and returns an iterate that has left the optimum. This
# value is the one the sweeps in that regime reach, and it is what the
# non-convergence advice names.
qp_reachable_tolerance <- 1e-6

resolved_qp_tolerance <- function(method) {
  tolerance <- method@convergence_tolerance
  if (is.null(tolerance)) qp_default_tolerance else tolerance
}

# The non-convergence advice for the quadratic-program family, which fails its
# criterion for a reason the estimating-equation family does not share. A descent
# method that spends its iteration cap stopped short of the answer and is helped
# by a larger cap; an alternating-direction iteration on an indefinite form that
# spends its cap has usually passed the residual floor of its problem, past which
# each further iteration moves away from the optimum rather than toward it. So
# the advice leads with the tolerance, names a value the problem can usually
# reach when the fit asked for something tighter, and keeps the cap for last.
#
# Only the indefinite forms carry that floor, so only they call the cap a last
# resort. A positive-semidefinite form keeps descending toward its tolerance for
# as long as the cap allows, and a run that spent the cap there really did stop
# short, so the cap is named as an ordinary lever.
#
# The weights caveat is about a solve that met no tolerance at all rather than
# about one that missed the tolerance asked for. An energy fit that could not
# reach its tolerance reports the iterate of a re-solve at a reachable one and
# still calls itself unconverged, and telling that caller the weights are
# worthless would contradict the advice above it.
quadratic_program_convergence_bullets <- function(method) {
  loosen <- if (resolved_qp_tolerance(method) < qp_reachable_tolerance) {
    "Loosen {.arg convergence_tolerance} in {.fn {class(method)[1]}}, which the problem can usually reach at {.val {qp_reachable_tolerance}}."
  } else {
    "Loosen {.arg convergence_tolerance} in {.fn {class(method)[1]}}."
  }
  cap <- if (has_indefinite_objective(method)) {
    "Raising {.arg max_iterations} is the last resort, and helps only a solve that stopped short of the residual floor rather than past it."
  } else {
    "Raising {.arg max_iterations} is the other lever, since this objective descends toward its tolerance for as long as the cap allows."
  }
  c(
    "The solver did not reach its convergence tolerance.",
    i = loosen,
    x = "The weights of a solve that met no tolerance at all should not be relied on.",
    i = cap
  )
}

# Raise or warn on the solver outcome. The quadratic-program family reports a
# terminal status the backend assigns, so an infeasible constraint set raises
# `balancing_infeasible_error` and a hard solver failure raises
# `balancing_convergence_error`, each naming the knob to turn; a reached iteration
# cap with a usable iterate still warns. The estimating-equation family carries no
# status and warns when it did not meet its convergence tolerance. A fit that
# retried with a second solver carries a record of what it ran, and the warning
# names them, so a caller can see that a different solver has already been tried
# and the remaining levers are the iteration cap and the tolerance.
#
# Which knob a failure names is chosen by the status, so every status a backend
# can assign is routed here. Only a status that genuinely means the solve ran out
# of iterations falls through to the closing warning; a solve that broke down
# numerically or stalled would not be helped by either knob, so it reports the
# conditioning of the problem instead.
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
      advice <- if (tunes_weight_penalty(method)) {
        "Check the covariates for collinearity, or raise {.arg weight_penalty} in {.fn {class(method)[1]}}."
      } else {
        "Check the covariates for collinearity."
      }
      abort(
        c(
          "The solver failed to produce a valid solution.",
          x = "It terminated with status {.val {fit$status}}.",
          i = advice
        ),
        error_class = "balancing_convergence_error",
        call = call
      )
    }
    if (
      fit$status %in%
        c("numerical_error", "insufficient_progress", "not_solved")
    ) {
      advice <- if (tunes_weight_penalty(method)) {
        "Rescale the covariates, loosen {.arg convergence_tolerance}, or raise {.arg weight_penalty} in {.fn {class(method)[1]}}."
      } else {
        "Rescale the covariates, or loosen {.arg convergence_tolerance} in {.fn {class(method)[1]}}."
      }
      abort(
        c(
          "The solver stopped without a solution.",
          x = "It terminated with status {.val {fit$status}}.",
          i = advice
        ),
        error_class = "balancing_convergence_error",
        call = call
      )
    }
  }
  if (!isTRUE(fit$converged)) {
    tried <- solver_labels(fit$solvers_tried)
    bullets <- if (S7::S7_inherits(method, quadratic_program_method)) {
      quadratic_program_convergence_bullets(method)
    } else {
      c(
        "The solver did not reach its convergence tolerance.",
        i = "Increase {.arg max_iterations} or loosen {.arg convergence_tolerance} in {.fn {class(method)[1]}}."
      )
    }
    if (length(tried) > 1L) {
      bullets[[1L]] <- "Neither solver reached its convergence tolerance."
      bullets <- append(bullets, c(x = "The fit tried {tried}."), after = 1L)
    }
    warn(
      bullets,
      warning_class = "balancing_convergence_warning",
      call = call
    )
  }
  invisible()
}

# Drop the exposure column from a resolved covariate selection. A tidyselect
# expression resolves against the whole data frame, so a selection such as
# `everything()` reaches the exposure as well, and balancing the exposure against
# itself is infeasible by construction: no reweighting of a group makes its own
# exposure indicator match the pooled mean. Excluding the response from
# predictors chosen by selection is the established modeling idiom, so the fit
# proceeds rather than erroring, and the exclusion is announced so that a caller
# who named the exposure deliberately can see it did not become a constraint.
drop_exposure_covariate <- function(selection, exposure_pos, exposure_name) {
  keep <- selection != exposure_pos
  if (all(keep)) {
    return(selection)
  }
  alert_info(
    "Dropping the exposure {.val {exposure_name}} from {.arg .covariates}."
  )
  selection[keep]
}

# Refuse a fit whose constraint set is empty for a method the constraints define.
# Entropy balancing, inverse probability tilting, the covariate balancing
# propensity score, and stable balancing weights derive their weights from the
# constraints, so with no constraint columns there is no problem to solve and the
# solver would be handed a design with no parameters. The objective-driven
# quadratic programs, energy and characteristic function distance balancing, have
# an objective that defines the solution on its own and fit as usual, so they
# answer `requires_constraints()` with FALSE.
check_constraint_columns <- function(
  method,
  matrix,
  covariates,
  call = rlang::caller_env()
) {
  if (ncol(matrix) > 0 || !requires_constraints(method)) {
    return(invisible(NULL))
  }
  abort(
    c(
      "{method_label(method)} must have at least one balance constraint.",
      x = "The covariate{?s} {.val {covariates}} contributed no constraint columns.",
      i = "Raise {.arg moments} in {.fn balance_terms}, or balance {.arg quantiles} or {.arg interactions} instead."
    ),
    error_class = "balancing_constraints_error",
    call = call
  )
}

# Whether a method's weights are defined by its balance constraints, so an empty
# constraint set leaves it nothing to solve. The estimating-equation family reads
# its identifying conditions off the constraints; the quadratic-program family is
# driven by its objective, with stable balancing weights the exception that needs
# explicit constraints, mirroring how `default_constraints()` splits the families.
requires_constraints <- new_generic("requires_constraints", "method")

method(requires_constraints, balance_method) <- function(method) {
  TRUE
}

method(requires_constraints, quadratic_program_method) <- function(method) {
  FALSE
}

# Whether a method offers `weight_penalty` as an argument a caller can raise, which
# decides whether the conditioning advice a solver breakdown gives is allowed to
# name it. Energy and characteristic function distance balancing take the penalty
# as a tuning argument; stable balancing weights minimize the weight dispersion
# alone and take no such argument. The property cannot answer this, since stable
# balancing weights hold the inherited penalty at zero rather than dropping it, so
# every quadratic program appears to carry one.
tunes_weight_penalty <- new_generic("tunes_weight_penalty", "method")

method(tunes_weight_penalty, balance_method) <- function(method) {
  FALSE
}

# Whether a method assembles an indefinite quadratic form, which decides whether
# the non-convergence advice may speak of a residual floor. Energy balancing and
# the characteristic function distance energy kernel build their objective from
# the negative pairwise distance, which is conditionally positive semidefinite
# alone and indefinite as a quadratic form; every other objective the package
# assembles is positive semidefinite. The distinction is a property of the
# objective rather than of the family, so it is dispatched on the method and the
# kernel rather than read off the quadratic-program parent.
has_indefinite_objective <- new_generic("has_indefinite_objective", "method")

method(has_indefinite_objective, balance_method) <- function(method) {
  FALSE
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

# The estimand vocabulary, in one place because two entry points match against
# it: the estimand that creates a fit and the estimand that names one afterwards
# in `ipw()`. A spelling either accepts has to be a spelling the other accepts,
# which a second literal list would drift away from.
estimand_choices <- function() {
  c("ate", "att", "atc", "atu", "ato")
}

# The vocabulary accepts "atc" for the untreated target, matching propensity, and
# a fit stores the single canonical spelling "atu". Anything that compares a
# caller's estimand against a stored one canonicalizes first, so the spelling
# that created a fit is also a spelling that names it afterwards.
canonical_estimand <- function(estimand) {
  if (identical(estimand, "atc")) "atu" else estimand
}

# Refuse an exposure that takes a single level. Balancing reweights one exposure
# group toward another, so one level leaves nothing to balance whatever the
# estimand was asked for: a focal estimand would hold the whole sample fixed and
# carry no parameter blocks, and a pooled estimand would return the weighting it
# started from and then have no arm-to-arm contrast to report. A declared factor
# level no observation takes is dropped before the count, so this measures the
# levels the data carry. A continuous exposure has no levels at all and is not
# measured by the rule.
validate_exposure_level_count <- function(
  levels,
  exposure_type,
  call = rlang::caller_env()
) {
  if (identical(exposure_type, "continuous") || length(levels) >= 2) {
    return(invisible())
  }
  abort(
    c(
      "Balancing needs an exposure with at least two levels.",
      x = "The exposure takes the single level {.val {levels}}.",
      i = "Supply an exposure whose values differ across the sample."
    ),
    error_class = "balancing_estimand_error",
    call = call
  )
}

# The covariates that kept at least one constraint column, in the order they were
# selected. The expansion drops constant and aliased columns, and a `moments`
# request of zero contributes none, so a covariate can be selected and constrain
# nothing; listing it on the fit would claim balance the fit never targeted. The
# request stays visible in the recorded call. An interaction column constrains
# both of its factors, so a covariate that appears only as a partner still
# counts.
constrained_covariates <- function(recipe, covariates) {
  sources <- unlist(
    lapply(recipe, function(record) c(record$source, record$partner)),
    use.names = FALSE
  )
  covariates[covariates %in% sources[!is.na(sources)]]
}

# Resolve the focal exposure level for att and atc. A binary exposure infers the
# treated level (the second level) for att and the control level (the first) for
# atc; a categorical exposure requires an explicit `.focal_level`. The average
# treatment effect and the overlap estimand reweight every group rather than hold
# one fixed, so they carry no focal level.
resolve_focal_level <- function(
  estimand,
  exposure_type,
  levels,
  focal_level,
  call = rlang::caller_env()
) {
  # A pooled estimand resolves no focal level, so a supplied one never reaches
  # the check that it names an exposure level and never reaches the fit either. A
  # level that does not exist used to be accepted in silence, which reads as a
  # fit that targeted it, so the request is announced as ignored instead. The
  # convention is the one bw_cbps() uses for a tuning argument its exposure type
  # cannot act on.
  if (estimand %in% c("ate", "ato")) {
    if (!is.null(focal_level)) {
      warn(
        c(
          "{.arg .focal_level} applies to the {.val att} and {.val atc} estimands and is ignored.",
          i = "The {.val {estimand}} estimand reweights every exposure group rather than holding one fixed."
        ),
        warning_class = "balancing_ignored_argument_warning",
        call = call
      )
    }
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
          "{.arg .focal_level} is required for the {.val {estimand}} estimand with a categorical exposure.",
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
        "{.arg .focal_level} must be an exposure level.",
        x = "{.val {resolved}} is not one of {.val {levels}}."
      ),
      error_class = "balancing_estimand_error",
      call = call
    )
  }

  resolved
}
