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

#' Validate that a data frame has at least one row
#'
#' @param .data The data frame to validate.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `.data`, invisibly, when it has at least one row.
#' @keywords internal
#' @noRd
validate_nonempty <- function(
  .data,
  arg_name = ".data",
  call = rlang::caller_env()
) {
  if (nrow(.data) == 0L) {
    abort(
      c(
        "{.arg {arg_name}} must have at least one row.",
        x = "It has no rows.",
        i = "Balancing weights require observations to reweight."
      ),
      error_class = "balancing_empty_error",
      call = call
    )
  }
  invisible(.data)
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
#' `abort_infinite()`. A non-numeric vector holds neither, so it passes through.
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
