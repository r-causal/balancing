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

# Refuse an explicit exposure type the data cannot represent at all. An explicit
# type is the caller's declaration of how the exposure is modeled and wins over
# the detection heuristics, matching propensity's match_exposure_type: a
# ten-level numeric dose declared continuous is fit as a dose, and a numeric
# exposure taking two values declared continuous is the caller's decision. Only a
# declaration the data cannot carry is a contradiction: a binary exposure needs
# exactly two distinct values, and a continuous one needs a numeric vector. Every
# vector can be read as a set of levels, so a categorical declaration is always
# possible.
check_forced_type <- function(
  forced,
  exposure_vec,
  detected,
  call = rlang::caller_env()
) {
  problem <- switch(
    forced,
    binary = if (!has_two_levels(exposure_vec)) {
      "A {.val binary} exposure takes exactly two distinct values, and this one takes {length(unique(exposure_vec))}."
    },
    continuous = if (!is.numeric(exposure_vec)) {
      "A {.val continuous} exposure is numeric, and this one is {.obj_type_friendly {exposure_vec}}."
    },
    categorical = NULL
  )
  if (is.null(problem)) {
    return(invisible(NULL))
  }
  abort(
    c(
      "{.arg exposure_type} was set to {.val {forced}}, but the exposure cannot be treated that way.",
      x = problem,
      i = "Drop {.arg exposure_type} to detect the type from the data, which reads it as {.val {detected}}."
    ),
    error_class = "balancing_exposure_type_error",
    call = call
  )
}

#' Resolve the exposure type for a fit
#'
#' Matches `exposure_type` against the permitted values. When it is `"auto"`, the
#' type is inferred from the data and announced through `alert_info()`. An
#' explicit type wins over the detection heuristics; only a type the data cannot
#' represent raises `balancing_exposure_type_error`, as does a type the method
#' does not support.
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
    check_forced_type(explicit, exposure_vec, detected, call = call)
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
