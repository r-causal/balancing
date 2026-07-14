# Diagnostic plots for a fitted balancing result. autoplot() draws one of three
# method-diagnostic views and plot() prints it. The generic is ggplot2's, so
# the methods are registered on an external generic wired up at load time by
# S7::methods_register() and only attach once ggplot2 is loaded; ggplot2 stays a
# suggested dependency. These views diagnose the optimization itself: the love
# plot contrasts the balance the weights achieved against the unweighted sample,
# the weight distribution shows how much dispersion the reweighting introduced,
# and the dual-variable chart shows which structural constraints bind. For a
# general covariate-balance assessment workflow, pass the weights to halfmoon.

# ggplot2 is a suggested dependency, so the autoplot method is registered on
# ggplot2's generic through an external generic. S7::methods_register() attaches
# it when ggplot2 loads.
autoplot <- new_external_generic("ggplot2", "autoplot", "object")

#' Diagnostic plots for a balancing fit
#'
#' @description
#' `autoplot()` draws one of three diagnostic views of a [balancing] fit, and
#' `plot()` prints that view. These plots diagnose the optimization: `"balance"`
#' is a love plot contrasting the unweighted and weighted absolute standardized
#' mean differences (exposure-covariate correlations for a continuous exposure)
#' against a dashed tolerance line; `"weights"` shows the distribution of the
#' balancing weights by exposure group; and `"duals"` shows the solver's dual
#' variables as a bar chart, so the constraints that bind the quadratic program
#' stand out.
#'
#' @details
#' The `"duals"` view needs the solver dual variables the quadratic-program
#' family reports. The estimating-equation family (entropy balancing, inverse
#' probability tilting, the covariate balancing propensity score) carries none,
#' so `"duals"` raises an error for those fits.
#'
#' These views target the balancing method itself. To assess covariate balance
#' more broadly, extract the weights with [weights()] and pass them to the
#' halfmoon package, whose love plots and distribution summaries cover a general
#' balance-assessment workflow.
#'
#' @param object,x A [balancing] fit.
#' @param type The view to draw, one of `"balance"`, `"weights"`, or `"duals"`.
#' @param y Not used; present for compatibility with the [plot()] generic.
#' @param ... Passed on. `plot()` forwards these arguments, including `type`, to
#'   `autoplot()`; `autoplot()` itself accepts no further arguments.
#'
#' @return A [ggplot2::ggplot] object. `plot()` returns it invisibly after
#'   printing.
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
#'
#' # Requires ggplot2.
#' if (rlang::is_installed("ggplot2")) {
#'   ggplot2::autoplot(fit, type = "balance")
#'   ggplot2::autoplot(fit, type = "weights")
#' }
#'
#' @name autoplot.balancing
#' @aliases plot.balancing
NULL

method(autoplot, balancing) <- function(
  object,
  type = c("balance", "weights", "duals"),
  ...
) {
  type <- rlang::arg_match(type)
  switch(
    type,
    balance = autoplot_balance(object),
    weights = autoplot_weights(object),
    duals = autoplot_duals(object)
  )
}

# plot() is aliased onto the autoplot page through @aliases in the block above.
# It carries no roxygen of its own so the generated Rd has no \usage section,
# which keeps the documented autoplot arguments valid without an S7-method usage
# line that R CMD check would reject.
method(plot, balancing) <- function(x, y, ...) {
  # plot() draws a ggplot, so ggplot2 must be present; call its generic
  # directly, which dispatches to the method registered above.
  drawing <- ggplot2::autoplot(x, ...)
  print(drawing)
  invisible(drawing)
}

# The love plot. Each constraint term contributes an unweighted and a weighted
# absolute balance statistic, so the point layer carries two rows per term. The
# point layer is added first so that a consumer reading the built layers finds
# the two-per-term balance data where it expects it. The dashed rule marks the
# per-constraint tolerance; an exact fit places it at zero.
autoplot_balance <- function(object) {
  table <- object@balance_table
  statistic <- table$statistic[[1]]
  axis_label <- if (identical(statistic, "correlation")) {
    "Absolute exposure-covariate correlation"
  } else {
    "Absolute standardized mean difference"
  }

  terms <- table$term
  plot_data <- vctrs::data_frame(
    term = factor(rep(terms, 2L), levels = rev(unique(terms))),
    sample = factor(
      rep(c("Unweighted", "Weighted"), each = length(terms)),
      levels = c("Unweighted", "Weighted")
    ),
    value = c(abs(table$unweighted), abs(table$weighted))
  )
  tolerance <- max(table$tolerance)

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$value,
      y = .data$term,
      color = .data$sample,
      shape = .data$sample
    )
  ) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_vline(
      xintercept = tolerance,
      linetype = "dashed",
      color = "grey40"
    ) +
    ggplot2::labs(
      x = axis_label,
      y = NULL,
      color = NULL,
      shape = NULL,
      title = "Balance before and after weighting"
    )
}

# The weight distribution. Weights are shown on the reported scale, faceted by
# exposure group, so the shape of the reweighting is visible per group. A
# continuous exposure has no groups, so it draws a single panel.
autoplot_weights <- function(object) {
  weights <- as.numeric(weights(object))
  groups <- attr(object@weights, "groups")
  if (is.null(groups)) {
    group <- factor(rep("overall", length(weights)))
  } else {
    labels <- character(length(weights))
    for (level in names(groups)) {
      labels[groups[[level]]] <- level
    }
    group <- factor(labels, levels = names(groups))
  }

  plot_data <- vctrs::data_frame(weight = weights, group = group)
  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$weight)) +
    ggplot2::geom_histogram(bins = 30, fill = "grey60", color = "white") +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$group),
      ncol = 1,
      scales = "free_y"
    ) +
    ggplot2::labs(
      x = "Balancing weight",
      y = "Observations",
      title = "Weight distribution by exposure group"
    )
}

# The dual-variable bar chart. Each structural constraint the solver enforced
# contributes one bar, colored by the kind of constraint. The estimating-equation
# family carries no dual variables, so the view is unavailable there.
autoplot_duals <- function(object) {
  duals <- object@duals
  if (is.null(duals)) {
    abort(
      c(
        "This fit has no dual variables to plot.",
        i = "Dual variables are reported by the quadratic-program family (energy, characteristic function distance, and stable balancing weights).",
        i = "The estimating-equation family carries none."
      ),
      error_class = "balancing_autoplot_duals_error"
    )
  }

  plot_data <- vctrs::data_frame(
    constraint = factor(duals$constraint, levels = duals$constraint),
    kind = factor(duals$kind),
    dual = duals$dual
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$constraint, y = .data$dual, fill = .data$kind)
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      x = "Constraint",
      y = "Dual variable",
      fill = NULL,
      title = "Dual variables of the balancing constraints"
    )
}
