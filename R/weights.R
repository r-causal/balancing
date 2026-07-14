# The bw class is a sibling of propensity::psw under the shared causal_wts
# parent: balancing weights target covariate balance directly rather than
# through a propensity score, so they are not psw, but they share the
# causal-weight vocabulary (estimand metadata, is_causal_wt()) and the same
# coercion policy of downgrading to a plain double when metadata conflicts.

#' Balancing weight vectors
#'
#' @description
#' `bw` objects are numeric vectors carrying the estimand a set of balancing
#' weights targets. They are the weight vector stored on a [balancing] result
#' and returned by its [weights()] method.
#'
#' `bw` is a sibling of [propensity::psw]: both inherit the `causal_wts` class,
#' so [propensity::is_causal_wt()] and [propensity::estimand()] work on either.
#'
#' @details
#' ## Constructors
#'
#' - `bw()` is the user-facing constructor. It coerces `x` to double before
#'   building the object.
#' - `new_bw()` is the low-level constructor for developers. It assumes `x` is
#'   already a double vector.
#' - `as_bw()` coerces an existing numeric vector to a `bw` object.
#'
#' ## Queries
#'
#' - `is_bw()` tests whether an object is a `bw` vector.
#' - [propensity::estimand()] reads the target estimand.
#'
#' ## Combining
#'
#' Arithmetic preserves the class and estimand, so normalizing weights keeps the
#' metadata. Combining `bw` vectors with matching estimands preserves the class;
#' combining vectors with different estimands, or a `bw` with a
#' [propensity::psw], warns and falls back to a plain double vector.
#'
#' @param x A numeric vector of weights for `bw()` and `new_bw()`, or an object
#'   to test or coerce for `is_bw()` and `as_bw()`.
#' @param estimand A single string naming the target estimand, or `NULL`.
#' @param ... Additional attributes stored on the object (developer use only).
#'
#' @return
#' - `new_bw()`, `bw()`, `as_bw()`: a `bw` vector.
#' - `is_bw()`: a single logical value.
#'
#' @examples
#' w <- bw(c(0.5, 1, 1.5), estimand = "ate")
#' w
#' is_bw(w)
#' propensity::estimand(w)
#'
#' # Arithmetic preserves the class.
#' w / sum(w)
#'
#' @name bw
#' @export
new_bw <- function(x = double(), estimand = NULL, ...) {
  vctrs::vec_assert(x, ptype = double())
  vctrs::new_vctr(
    x,
    estimand = estimand,
    ...,
    class = c("bw", "causal_wts"),
    inherit_base_type = TRUE
  )
}

#' @rdname bw
#' @export
bw <- function(x = double(), estimand = NULL) {
  x <- vctrs::vec_cast(x, to = double())
  attributes(x) <- NULL
  new_bw(x, estimand = estimand)
}

#' @rdname bw
#' @export
as_bw <- function(x, estimand = NULL) {
  x <- vctrs::vec_cast(x, to = double())
  bw(x, estimand = estimand)
}

#' @rdname bw
#' @export
is_bw <- function(x) {
  inherits(x, "bw")
}

# Signal the classed coercion warning when a bw is combined with a vector whose
# metadata cannot be preserved, mirroring propensity's downgrade policy.
warn_bw_downgrade <- function(other) {
  warn(
    c(
      "Cannot combine {.cls bw} weights with {.cls {other}}.",
      i = "The result is a plain {.cls numeric} vector without the estimand.",
      i = "Set a common estimand on both vectors to keep the {.cls bw} class."
    ),
    warning_class = "balancing_class_downgrade_warning"
  )
}

#' @export
vec_ptype_abbr.bw <- function(x, ...) {
  estimand <- propensity::estimand(x)
  if (is.null(estimand)) {
    "bw"
  } else {
    paste0("bw{", estimand, "}")
  }
}

#' @export
vec_ptype_full.bw <- function(x, ...) {
  estimand <- propensity::estimand(x)
  if (is.null(estimand)) {
    "bw{estimand = unknown}"
  } else {
    paste0("bw{estimand = ", estimand, "}")
  }
}

#' @export
vec_restore.bw <- function(x, to, ...) {
  if (inherits(x, "bw")) {
    x <- vctrs::vec_data(x)
  }
  new_bw(x, estimand = propensity::estimand(to))
}

#' @export
vec_ptype2.bw.bw <- function(x, y, ...) {
  if (!identical(propensity::estimand(x), propensity::estimand(y))) {
    warn_bw_downgrade("bw")
    return(double())
  }
  new_bw(estimand = propensity::estimand(x))
}

#' @export
vec_ptype2.bw.double <- function(x, y, ...) {
  double()
}

#' @export
vec_ptype2.double.bw <- function(x, y, ...) {
  double()
}

#' @export
vec_ptype2.bw.integer <- function(x, y, ...) {
  double()
}

#' @export
vec_ptype2.integer.bw <- function(x, y, ...) {
  double()
}

#' @export
vec_cast.bw.bw <- function(x, to, ...) {
  x
}

#' @export
vec_cast.bw.double <- function(x, to, ...) {
  bw(x, estimand = propensity::estimand(to))
}

#' @export
vec_cast.double.bw <- function(x, to, ...) {
  vctrs::vec_data(x)
}

#' @export
vec_cast.bw.integer <- function(x, to, ...) {
  bw(x, estimand = propensity::estimand(to))
}

#' @export
vec_cast.integer.bw <- function(x, to, ...) {
  vctrs::vec_cast(vctrs::vec_data(x), integer(), x_arg = "bw")
}

# Combining a bw with character shares no numeric common type, so the safe
# result is a plain character vector, matching propensity's psw policy.
#' @export
vec_ptype2.bw.character <- function(x, y, ...) {
  warn_bw_downgrade("character")
  character()
}

#' @export
vec_ptype2.character.bw <- function(x, y, ...) {
  warn_bw_downgrade("character")
  character()
}

#' @export
vec_cast.character.bw <- function(x, to, ...) {
  as.character(vctrs::vec_data(x))
}

# A bw combined with a psw shares the causal_wts parent but not the same weight
# semantics, so the safe common type is a plain double.
#' @export
vec_ptype2.bw.psw <- function(x, y, ...) {
  warn_bw_downgrade("psw")
  double()
}

#' @export
vec_ptype2.psw.bw <- function(x, y, ...) {
  warn_bw_downgrade("psw")
  double()
}

#' @export
#' @method vec_arith bw
vec_arith.bw <- function(op, x, y, ...) {
  UseMethod("vec_arith.bw", y)
}

#' @export
#' @method vec_arith.bw default
vec_arith.bw.default <- function(op, x, y, ...) {
  vctrs::stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.bw bw
vec_arith.bw.bw <- function(op, x, y, ...) {
  estimand_x <- propensity::estimand(x)
  estimand_y <- propensity::estimand(y)
  estimand <- if (identical(estimand_x, estimand_y)) {
    estimand_x
  } else {
    paste0(estimand_x, ", ", estimand_y)
  }
  bw(vctrs::vec_arith_base(op, x, y), estimand = estimand)
}

#' @export
#' @method vec_arith.bw numeric
vec_arith.bw.numeric <- function(op, x, y, ...) {
  vctrs::vec_restore(vctrs::vec_arith_base(op, x, y), x)
}

#' @export
#' @method vec_arith.numeric bw
vec_arith.numeric.bw <- function(op, x, y, ...) {
  vctrs::vec_restore(vctrs::vec_arith_base(op, x, y), y)
}

#' @export
#' @method vec_arith.bw integer
vec_arith.bw.integer <- function(op, x, y, ...) {
  vctrs::vec_restore(vctrs::vec_arith_base(op, x, y), x)
}

#' @export
#' @method vec_arith.bw MISSING
vec_arith.bw.MISSING <- function(op, x, y, ...) {
  switch(
    op,
    `-` = vctrs::vec_restore(-vctrs::vec_data(x), x),
    `+` = x,
    vctrs::stop_incompatible_op(op, x, y)
  )
}

#' @export
vec_math.bw <- function(.fn, .x, ...) {
  if (.fn %in% c("cumsum", "cumprod", "cummin", "cummax")) {
    return(vctrs::vec_restore(
      vctrs::vec_math_base(.fn, vctrs::vec_data(.x), ...),
      .x
    ))
  }
  vctrs::vec_math_base(.fn, vctrs::vec_data(.x), ...)
}

#' @export
Summary.bw <- function(..., na.rm = FALSE) {
  args <- lapply(list(...), vctrs::vec_data)
  do.call(.Generic, c(args, list(na.rm = na.rm)))
}

#' @export
min.bw <- function(..., na.rm = FALSE) {
  args <- lapply(list(...), vctrs::vec_data)
  do.call("min", c(args, list(na.rm = na.rm)))
}

#' @export
max.bw <- function(..., na.rm = FALSE) {
  args <- lapply(list(...), vctrs::vec_data)
  do.call("max", c(args, list(na.rm = na.rm)))
}

# median() and quantile() are not part of the Summary group generic, so a bw
# vector needs its own methods to reach them, mirroring propensity's psw. Both
# operate on the underlying double and return a plain numeric summary.
#' @importFrom stats median
#' @export
median.bw <- function(x, na.rm = FALSE, ...) {
  stats::median(vctrs::vec_data(x), na.rm = na.rm, ...)
}

#' @importFrom stats quantile
#' @export
quantile.bw <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE, ...) {
  stats::quantile(vctrs::vec_data(x), probs = probs, na.rm = na.rm, ...)
}

#' @export
`[.bw` <- function(x, i, ...) {
  if (missing(i)) {
    return(NextMethod())
  }
  if (is.matrix(i) || is.array(i)) {
    return(vctrs::vec_data(x)[i, ...])
  }
  NextMethod()
}

# ---- weights() and ess() ---------------------------------------------------

#' Extract balancing weights
#'
#' Returns the [bw] weight vector stored on a [balancing] result. When sampling
#' weights are present they are composed onto the balancing weights by default.
#'
#' @usage NULL
#' @param object A [balancing] result.
#' @param include_sampling_weights Whether to multiply the balancing weights by
#'   the sampling weights. Defaults to `TRUE`.
#' @param ... Ignored.
#'
#' @return A [bw] vector.
#' @export
method(weights, balancing) <- function(
  object,
  ...,
  include_sampling_weights = TRUE
) {
  rlang::check_dots_empty()
  w <- object@weights
  if (include_sampling_weights && !is.null(object@sampling_weights)) {
    w <- w * object@sampling_weights
  }
  w
}

#' Effective sample size for a balancing fit
#'
#' Reports the effective sample size implied by a set of balancing weights. The
#' effective sample size within a group is `sum(w)^2 / sum(w^2)`, which equals
#' the group size when the weights are uniform and shrinks as the weights become
#' more variable. Discrete exposures report one row per exposure level;
#' continuous exposures report a single overall row.
#'
#' This method registers on [causalgenerics::ess()], the shared effective sample
#' size generic. `library(balancing)` re-exports the generic, so `ess(fit)`
#' works without a second attachment.
#'
#' @param x A [balancing] result.
#' @param ... Ignored.
#'
#' @return A tibble with columns `group`, `n`, and `ess`.
#'
#' @examples
#' n <- 200
#' x1 <- rnorm(n)
#' df <- data.frame(exposure = rbinom(n, 1, plogis(0.5 * x1)), x1 = x1)
#' fit <- balance(df, exposure, x1, method = bw_entropy())
#' ess(fit)
#'
#' @name ess.balancing
NULL

causalgenerics_ess <- new_external_generic("causalgenerics", "ess", "x")

method(causalgenerics_ess, balancing) <- function(x, ...) {
  w <- as.numeric(weights(x))
  effective <- function(weights) sum(weights)^2 / sum(weights^2)

  groups <- attr(x@weights, "groups")
  if (is.null(groups)) {
    return(new_balancing_tibble(list(
      group = "overall",
      n = x@n,
      ess = effective(w)
    )))
  }

  levels <- names(groups)
  rows <- lapply(levels, function(level) {
    idx <- groups[[level]]
    list(group = level, n = length(idx), ess = effective(w[idx]))
  })
  new_balancing_tibble(list(
    group = vapply(rows, function(r) r$group, character(1)),
    n = vapply(rows, function(r) as.integer(r$n), integer(1)),
    ess = vapply(rows, function(r) r$ess, numeric(1))
  ))
}
