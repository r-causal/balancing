# Exposure-type resolution. The heuristics that read a type off the data, the
# announcement that reports the one they read, and the refusal of a declared
# type the data cannot carry all belong to causalgenerics, so that a user meets
# one r-causal ecosystem: propensity, positively, and balancing resolve an
# exposure the same way and say the same thing about it. What stays here is the
# question that is balancing's own, which exposure types a given method can fit.

#' Resolve the exposure type for a fit
#'
#' Matches `exposure_type` against the permitted values. When it is `"auto"`,
#' [causalgenerics::detect_exposure_type()] infers the type from the data and
#' announces it. An explicit type wins over the detection heuristics; only a
#' type the data cannot represent at all is refused, by
#' [causalgenerics::check_forced_type()] and with its
#' `causalgenerics_forced_exposure_type` condition. A type the method does not
#' support is balancing's own refusal and raises
#' `balancing_exposure_type_error`.
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

  # Detection is the only path that speaks. causalgenerics announces from
  # inside the detection it performs, so `announce` carries balancing's own
  # quiet option and `options(balancing.quiet = TRUE)` stays the control a
  # balancing user reaches for. The refusal path detects a second time, to name
  # in its message what the data do read as, and asks for silence while it
  # does, so an explicit type announces nothing whether it stands or not.
  if (explicit == "auto") {
    resolved <- causalgenerics::detect_exposure_type(
      exposure_vec,
      announce = !be_quiet(),
      call = call
    )
  } else {
    causalgenerics::check_forced_type(explicit, exposure_vec, call = call)
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
