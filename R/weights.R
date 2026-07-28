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
#' so [causalgenerics::is_causal_wt()] and [causalgenerics::estimand()] work on
#' either.
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
#' - [causalgenerics::estimand()] reads the target estimand.
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
#' estimand(w)
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
  estimand <- estimand(x)
  if (is.null(estimand)) {
    "bw"
  } else {
    paste0("bw{", estimand, "}")
  }
}

#' @export
vec_ptype_full.bw <- function(x, ...) {
  estimand <- estimand(x)
  if (is.null(estimand)) {
    "bw{estimand = unknown}"
  } else {
    paste0("bw{estimand = ", estimand, "}")
  }
}

# `groups` records which rows of the fit belong to each exposure level, so it
# describes one particular vector at one particular length. What the restoration
# keys on is that length: a result that comes back at the size of the object it
# restores from keeps the attribute, and anything else drops it. Arithmetic,
# unary negation, and cumulative math are the operations that hold the size,
# and they rewrite each element where it stands, so the recorded positions still
# describe the rows they name.
#
# Equal size is the available condition rather than the exact one. A reordering
# or a repeat can also arrive at the original size, and those keep positions the
# data no longer matches. Nothing here can tell them apart: `vec_restore()`
# receives the restored data and the object it came from, never the index that
# produced it, so the operation that reordered the rows is not among its
# arguments. Subsetting `groups` alongside the data would take a hook that is
# handed that index.
#' @export
vec_restore.bw <- function(x, to, ...) {
  if (inherits(x, "bw")) {
    x <- vctrs::vec_data(x)
  }
  groups <- if (vctrs::vec_size(x) == vctrs::vec_size(to)) {
    attr(to, "groups")
  } else {
    NULL
  }
  new_bw(x, estimand = estimand(to), groups = groups)
}

#' @export
vec_ptype2.bw.bw <- function(x, y, ...) {
  if (!identical(estimand(x), estimand(y))) {
    warn_bw_downgrade("bw")
    return(double())
  }
  new_bw(estimand = estimand(x))
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
  bw(x, estimand = estimand(to))
}

#' @export
vec_cast.double.bw <- function(x, to, ...) {
  vctrs::vec_data(x)
}

#' @export
vec_cast.bw.integer <- function(x, to, ...) {
  bw(x, estimand = estimand(to))
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
  estimand_x <- estimand(x)
  estimand_y <- estimand(y)
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

# Elementwise and cumulative math, the Summary group generic, min(), max(),
# median(), quantile(), and subsetting all read the underlying double and need
# nothing a bw knows that a causal_wts does not, so causalgenerics supplies them
# on the shared parent. A method registered on bw would shadow the inherited one
# outright, since UseMethod() takes the first match down the class vector.

# ---- weights() -------------------------------------------------------------

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
