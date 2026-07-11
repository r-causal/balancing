# Exposure-type detection. balancing mirrors the propensity and positively idiom
# so that users experience one r-causal ecosystem: the resolution logic and the
# informational announcement match.

# Does a vector have exactly two distinct values?
has_two_levels <- function(x) {
  length(unique(x)) == 2
}

# Is a numeric exposure categorical by the unique-value heuristic? A numeric
# exposure is treated as categorical when the share of unique values among the
# non-missing observations falls below 20 percent.
is_categorical <- function(.exposure) {
  n_non_na <- sum(!is.na(.exposure))
  if (n_non_na == 0) {
    return(FALSE)
  }
  ratio <- length(unique(.exposure)) / n_non_na
  if (is.nan(ratio)) {
    return(FALSE)
  }
  ratio < 0.2
}

# Classify an exposure as "binary", "categorical", or "continuous". A vector with
# exactly two distinct values is binary; a factor or character vector with more
# than two values is categorical; a numeric vector is categorical when its share
# of unique values is small, and continuous otherwise.
detect_exposure_type <- function(.exposure) {
  if (has_two_levels(.exposure)) {
    "binary"
  } else if (is.factor(.exposure) || is.character(.exposure)) {
    "categorical"
  } else if (is_categorical(.exposure)) {
    "categorical"
  } else {
    "continuous"
  }
}

# Is a forced exposure type consistent with the data? A forced type that the data
# cannot support (for example "binary" on a many-valued continuous exposure) is a
# contradiction the caller must resolve.
forced_type_matches <- function(forced, .exposure) {
  switch(
    forced,
    binary = has_two_levels(.exposure),
    categorical = is.factor(.exposure) ||
      is.character(.exposure) ||
      is_categorical(.exposure),
    continuous = is.numeric(.exposure) &&
      !has_two_levels(.exposure) &&
      !is_categorical(.exposure)
  )
}

#' Resolve the exposure type for a fit
#'
#' Matches `exposure_type` against the permitted values. When it is `"auto"`, the
#' type is inferred from the data and announced through `alert_info()`. An
#' explicit type is honored, but a type the data contradict raises
#' `balancing_exposure_type_error`, as does a type the method does not support.
#'
#' @param exposure_type The `exposure_type` argument.
#' @param exposure_vec The exposure vector.
#' @param method The [balance_method] specification.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A single string: `"binary"`, `"categorical"`, or `"continuous"`.
#' @keywords internal
#' @noRd
resolve_exposure_type <- function(
  exposure_type,
  exposure_vec,
  method,
  call = rlang::caller_env()
) {
  explicit <- rlang::arg_match(
    exposure_type,
    c("auto", "binary", "categorical", "continuous"),
    error_call = call
  )
  detected <- detect_exposure_type(exposure_vec)

  if (explicit == "auto") {
    alert_info("Treating {.arg .exposure} as {detected}.")
    resolved <- detected
  } else {
    if (!forced_type_matches(explicit, exposure_vec)) {
      abort(
        c(
          "{.arg exposure_type} was set to {.val {explicit}}, but the data do not support it.",
          x = "The exposure is detected as {.val {detected}}.",
          i = "Drop {.arg exposure_type} to detect it automatically, or supply an exposure of the forced type."
        ),
        error_class = "balancing_exposure_type_error",
        call = call
      )
    }
    resolved <- explicit
  }

  supported <- supported_exposure_types(method)
  if (!resolved %in% supported) {
    abort(
      c(
        "{method_label(method)} does not support a {.val {resolved}} exposure.",
        i = "Supported exposure types are {.val {supported}}."
      ),
      error_class = "balancing_exposure_type_error",
      call = call
    )
  }

  resolved
}
