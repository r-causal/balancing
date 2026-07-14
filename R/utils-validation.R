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

#' Validate that the exposure and covariates have no missing values
#'
#' @param exposure_vec The exposure vector.
#' @param data The data frame.
#' @param covariates The covariate column names.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `NULL`, invisibly, when nothing is missing.
#' @keywords internal
#' @noRd
validate_no_missing <- function(
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
  invisible(NULL)
}

#' Validate sampling weights
#'
#' @param weights The evaluated sampling weights.
#' @param n The number of observations.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `weights`, invisibly, when they are non-negative and length `n`.
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
  if (anyNA(weights)) {
    abort(
      "{.arg sampling_weights} must not contain missing values.",
      error_class = "balancing_missing_error",
      call = call
    )
  }
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
  invisible(weights)
}
