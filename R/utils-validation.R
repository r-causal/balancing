# Validation helpers. Each raises a classed `balancing_error` on failure and
# returns its input invisibly on success, so validation composes inside the
# orchestrator without altering the value it guards.

#' Validate that an input is a data frame
#'
#' @param .data The object to validate.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `.data`, invisibly, when it is a data frame.
#' @keywords internal
#' @noRd
validate_data_frame <- function(
  .data,
  arg_name = ".data",
  call = rlang::caller_env()
) {
  if (!is.data.frame(.data)) {
    data_class <- class(.data)[1]
    abort(
      "{.arg {arg_name}} must be a data frame, not a {.cls {data_class}}.",
      error_class = "balancing_type_error",
      call = call
    )
  }
  invisible(.data)
}

#' Validate that a data frame has enough rows to balance
#'
#' Two observations are the smallest sample a fit can take, so an absent sample
#' and a sample of one are one refusal rather than two. Every numeric constraint
#' column crosses its boundary standardized by its own spread, and the standard
#' deviation of a single value is a missing value, so the rescue for a column
#' with no spread would steer on that missing value and the fit would die with
#' an unclassed comparison error. Nothing earlier turns a one-row sample away: a
#' lone exposure value is one unique value among one observation, which the
#' unique-value heuristic reads as continuous, so the rule requiring two
#' exposure levels never measures it.
#'
#' Stating the true minimum once matters to the caller who acts on it. Refusing
#' no rows for want of one, then refusing the row they add for want of two,
#' reports a single rule as two, so both counts are measured against the same
#' stated requirement and the bullet names the size that arrived.
#'
#' The requirement belongs to the data rather than to a method, so it is checked
#' ahead of the exposure and constraint machinery. A one-row sample therefore
#' reports its size whatever its columns hold, rather than reporting through the
#' constant-column refusal a factor-only selection would otherwise reach.
#'
#' @param .data The data frame to validate.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `.data`, invisibly, when it has at least two rows.
#' @keywords internal
#' @noRd
validate_row_count <- function(
  .data,
  arg_name = ".data",
  call = rlang::caller_env()
) {
  n <- nrow(.data)
  if (n >= 2L) {
    return(invisible(.data))
  }
  abort(
    c(
      "{.arg {arg_name}} must have at least two rows.",
      x = if (n == 0L) "It has no rows." else "It has {n} row{?s}.",
      i = "Balancing weights reweight a sample toward a target measured on the sample's own spread, which needs more than one observation."
    ),
    error_class = "balancing_empty_error",
    call = call
  )
}

#' Validate a resolved column selection
#'
#' @param selection A named integer vector of resolved column positions.
#' @param arg_name The argument name used in error messages.
#' @param expected One of `"one"` (exactly one column) or `"some"` (at least
#'   one).
#' @param call The calling environment, used to build the error's call.
#'
#' @return `selection`, invisibly, when it selects the required columns.
#' @keywords internal
#' @noRd
validate_selection <- function(
  selection,
  arg_name,
  expected = c("some", "one"),
  call = rlang::caller_env()
) {
  expected <- match.arg(expected)
  if (expected == "one" && length(selection) != 1) {
    abort(
      c(
        "{.arg {arg_name}} must select exactly one column.",
        x = "It selected {length(selection)} column{?s}."
      ),
      error_class = "balancing_selection_error",
      call = call
    )
  }
  if (expected == "some" && length(selection) == 0) {
    abort(
      "{.arg {arg_name}} must select at least one column.",
      error_class = "balancing_selection_error",
      call = call
    )
  }
  invisible(selection)
}

#' Report an infinite value in a numeric input
#'
#' The message and class every infinite numeric input reports, kept in one place
#' so an infinity in a covariate column, in the exposure, or in the sampling
#' weights reads the same wherever it came from. The class is
#' `balancing_range_error`, matching the sign and length failures of the same
#' arguments: an infinity is a value a fit cannot place rather than one it is
#' missing, so the remedy is to correct the value, not to complete it.
#'
#' @param arg_name The argument name used in the error message.
#' @param location Names of the columns at fault, for an argument spanning
#'   several, or `NULL` to report the count of offending values instead.
#' @param count The number of offending values, used when `location` is `NULL`.
#' @param call The calling environment, used to build the error's call.
#'
#' @return Nothing; always raises an error.
#' @keywords internal
#' @noRd
abort_infinite <- function(
  arg_name,
  location = NULL,
  count = NULL,
  call = rlang::caller_env()
) {
  abort(
    c(
      "{.arg {arg_name}} must not contain infinite values.",
      x = if (is.null(location)) {
        "Found {count} infinite value{?s}."
      } else {
        "Infinite values in {.val {location}}."
      },
      i = "A fit centers and scales every numeric input, which an infinity leaves undefined."
    ),
    error_class = "balancing_range_error",
    call = call
  )
}

#' Validate that a numeric input is finite
#'
#' The finiteness gate for a single numeric input. `is.finite()` rejects a missing
#' value, a `NaN`, and an infinity alike, and each of the three poisons the
#' centering and scaling every numeric input passes through on the way to a
#' solver. Both failures are refused here rather than left to surface as the
#' unclassed comparison error a poisoned mean or standard deviation eventually
#' produces. A missing value keeps the `balancing_missing_error` class the
#' exposure and covariate gates use; an infinity reports through
#' `abort_infinite()`. A non-numeric vector still meets the missing-value check,
#' since `anyNA()` reads a missing value of any type; it is the infinity half
#' alone that has nothing to find there.
#'
#' @param x The input to validate.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `x`, invisibly, when every value is finite.
#' @keywords internal
#' @noRd
validate_finite <- function(x, arg_name, call = rlang::caller_env()) {
  if (anyNA(x)) {
    abort(
      "{.arg {arg_name}} must not contain missing values.",
      error_class = "balancing_missing_error",
      call = call
    )
  }
  count <- sum(is.infinite(x))
  if (count > 0) {
    abort_infinite(arg_name, count = count, call = call)
  }
  invisible(x)
}

#' Validate that the exposure and covariates hold usable values
#'
#' The finiteness gate for the fit's data columns. Missingness and infinity are
#' reported separately because their remedies differ, and both name the columns at
#' fault so a wide covariate set points at the offender. The infinity checks read
#' the shared vocabulary in `abort_infinite()`; a non-numeric column can hold no
#' infinity and contributes only to the missingness check.
#'
#' @param exposure_vec The exposure vector.
#' @param data The data frame.
#' @param covariates The covariate column names.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `NULL`, invisibly, when every value is finite.
#' @keywords internal
#' @noRd
validate_finite_data <- function(
  exposure_vec,
  data,
  covariates,
  call = rlang::caller_env()
) {
  if (anyNA(exposure_vec)) {
    abort(
      c(
        "{.arg .exposure} must not contain missing values.",
        i = "balancing has no missingness-indicator machinery; complete the exposure first."
      ),
      error_class = "balancing_missing_error",
      call = call
    )
  }
  exposure_infinite <- sum(is.infinite(exposure_vec))
  if (exposure_infinite > 0) {
    abort_infinite(".exposure", count = exposure_infinite, call = call)
  }
  missing_covariates <- covariates[vapply(
    covariates,
    function(column) anyNA(data[[column]]),
    logical(1)
  )]
  if (length(missing_covariates) > 0) {
    abort(
      c(
        "{.arg .covariates} must not contain missing values.",
        x = "Missing values in {.val {missing_covariates}}.",
        i = "balancing has no missingness-indicator machinery; complete the covariates first."
      ),
      error_class = "balancing_missing_error",
      call = call
    )
  }
  infinite_covariates <- covariates[vapply(
    covariates,
    function(column) any(is.infinite(data[[column]])),
    logical(1)
  )]
  if (length(infinite_covariates) > 0) {
    abort_infinite(".covariates", location = infinite_covariates, call = call)
  }
  invisible(NULL)
}

#' Validate that the base measure carries mass where a fit needs it
#'
#' The base measure is the product of the sampling weights and any base weights,
#' and every group total a fit takes is taken under it: the constraint targets a
#' group is reweighted to, the total each group's reported weights are scaled to,
#' and the reference the balance table measures against. A group whose measure
#' sums to zero leaves each of those a ratio of zero totals, which reaches the
#' solver as a missing target and comes back as weights whose balance cannot be
#' assessed. The check runs before any fit, so every method reports the same
#' defect rather than the collinearity, infeasibility, or missing-value failure its
#' own path happens to reach first.
#'
#' Both vectors are validated non-negative before this, so a positive total is the
#' same condition as a positive value somewhere. A measure with no mass at all
#' needs both vectors to be nonzero somewhere and to be nonzero nowhere in common,
#' which neither vector's own all-zero check can see.
#'
#' [balance()] applies this to the sampling weights, which are the whole measure
#' for every method that carries no base weights. Entropy balancing applies it
#' again to the product its base weights form, once the base-weight length has been
#' checked where that fit reads them.
#'
#' @param measure The base measure, a non-negative numeric vector.
#' @param groups The per-level row indices, or `NULL` for a continuous exposure,
#'   which carries no groups.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `measure`, invisibly, when every required total is positive.
#' @keywords internal
#' @noRd
validate_base_measure <- function(
  measure,
  groups,
  call = rlang::caller_env()
) {
  if (!any(measure > 0)) {
    abort(
      c(
        "The base measure must carry some mass.",
        x = "It is zero for every observation.",
        i = "The base measure is the sampling weights times any base weights, so vectors that are nonzero nowhere in common leave no sample to reweight."
      ),
      error_class = "balancing_range_error",
      call = call
    )
  }
  if (!is.null(groups)) {
    empty <- names(groups)[vapply(
      groups,
      function(idx) !any(measure[idx] > 0),
      logical(1)
    )]
    if (length(empty) > 0) {
      abort(
        c(
          "Every exposure level must carry some base-measure mass.",
          x = "Exposure level{?s} {.val {empty}} {?has/have} a base measure of zero.",
          i = "The base measure is the sampling weights times any base weights; a level with none has no target to balance to and no total to report at."
        ),
        error_class = "balancing_range_error",
        call = call
      )
    }
  }
  invisible(measure)
}

#' Validate sampling weights
#'
#' A unit given no sampling weight is pinned at zero rather than dropped, so an
#' individual zero is a supported input and only a vector with no mass at all is
#' refused: a zero total leaves every weighted mean the fit targets undefined.
#'
#' @param weights The evaluated sampling weights.
#' @param n The number of observations.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `weights`, invisibly, when they are finite, non-negative, length `n`,
#'   and not zero throughout.
#' @keywords internal
#' @noRd
validate_sampling_weights <- function(
  weights,
  n,
  call = rlang::caller_env()
) {
  if (!is.numeric(weights)) {
    weights_class <- class(weights)[1]
    abort(
      "{.arg sampling_weights} must be numeric, not a {.cls {weights_class}}.",
      error_class = "balancing_type_error",
      call = call
    )
  }
  if (length(weights) != n) {
    abort(
      c(
        "{.arg sampling_weights} must have one value per observation.",
        x = "It has length {length(weights)}, but the data have {n} row{?s}."
      ),
      error_class = "balancing_range_error",
      call = call
    )
  }
  validate_finite(weights, "sampling_weights", call = call)
  if (any(weights < 0)) {
    abort(
      c(
        "{.arg sampling_weights} must be non-negative.",
        x = "Found {sum(weights < 0)} negative value{?s}."
      ),
      error_class = "balancing_range_error",
      call = call
    )
  }
  if (all(weights == 0)) {
    abort(
      c(
        "{.arg sampling_weights} must not be zero for every observation.",
        x = "Every weight is zero, which leaves no sample to reweight.",
        i = "Individual zero weights are supported; those units are pinned at zero."
      ),
      error_class = "balancing_range_error",
      call = call
    )
  }
  invisible(weights)
}
