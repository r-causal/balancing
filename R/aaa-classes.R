# aaa-classes.R sorts first so the abstract parents and the classes used as
# property types exist before the files that build on them load. S7 methods are
# registered here and wired up at load time by S7::methods_register() in
# .onLoad().

# Render a cli block to a string and write it to stdout. Print methods need their
# output on stdout so that both `print()` at the console and testthat's output
# capture see it; cli would otherwise divert to the message stream whenever
# stdout is redirected.
cat_cli <- function(expr) {
  cat(cli::cli_fmt(expr), sep = "\n")
}

# ---- Method specification hierarchy ---------------------------------------

#' The balancing method specification classes
#'
#' The method constructors ([bw_entropy()] and its siblings) return small
#' S7 objects that carry tuning parameters and never touch data. They share an
#' abstract hierarchy so that [balance()] can query a method's capabilities
#' uniformly. `balance_method` is the abstract root. The estimating-equation
#' family (entropy balancing, inverse probability tilting, the covariate
#' balancing propensity score) subclasses `estimating_equation_method`; the
#' quadratic-program family (energy, characteristic function distance, and
#' stable balancing weights) subclasses `quadratic_program_method`. None of the
#' abstract classes can be constructed directly.
#'
#' @param convergence_tolerance The solver convergence tolerance, or `NULL` for
#'   the core default.
#' @param max_iterations The maximum solver iterations, or `NULL` for the core
#'   default.
#' @param weight_penalty The L2 penalty on the weights.
#' @param min_weight The smallest permitted weight.
#'
#' @return An abstract class object. Constructing a concrete subclass returns a
#'   `balance_method`.
#' @name balance_method
#' @export
balance_method <- new_class(
  "balance_method",
  abstract = TRUE,
  properties = list(
    convergence_tolerance = NULL | class_double,
    max_iterations = NULL | class_integer
  ),
  validator = function(self) {
    if (
      !is.null(self@convergence_tolerance) && self@convergence_tolerance <= 0
    ) {
      "@convergence_tolerance must be a positive number"
    } else if (!is.null(self@max_iterations) && self@max_iterations < 0L) {
      "@max_iterations must be a non-negative whole number"
    }
  }
)

#' @rdname balance_method
#' @export
estimating_equation_method <- new_class(
  "estimating_equation_method",
  parent = balance_method,
  abstract = TRUE
)

#' @rdname balance_method
#' @export
quadratic_program_method <- new_class(
  "quadratic_program_method",
  parent = balance_method,
  abstract = TRUE,
  properties = list(
    weight_penalty = class_double,
    min_weight = class_double
  )
)

# ---- Capability generics ---------------------------------------------------

#' Method capability generics
#'
#' These generics let [balance()] interrogate a method specification without
#' knowing its concrete class. Each balancing method supplies one method per
#' generic.
#'
#' @param method A [balance_method] specification.
#' @param exposure_type One of `"binary"`, `"categorical"`, or `"continuous"`.
#' @param ... Additional context, such as the resolved `constraints`.
#'
#' @return `supported_exposure_types()` and `supported_estimands()` return
#'   character vectors; `supports_estimating_equations()` returns a single
#'   logical; `method_label()` returns a single string; `fit_method()` returns
#'   the internal solve result.
#' @name method_capabilities
#' @keywords internal
#' @export
supported_exposure_types <- new_generic("supported_exposure_types", "method")

#' @rdname method_capabilities
#' @export
supported_estimands <- new_generic(
  "supported_estimands",
  "method",
  function(method, exposure_type) {
    S7_dispatch()
  }
)

#' @rdname method_capabilities
#' @export
supports_estimating_equations <- new_generic(
  "supports_estimating_equations",
  "method",
  function(method, ...) {
    S7_dispatch()
  }
)

#' @rdname method_capabilities
#' @export
method_label <- new_generic("method_label", "method")

#' @rdname method_capabilities
#' @export
fit_method <- new_generic("fit_method", "method", function(method, prepared) {
  S7_dispatch()
})

# ---- balance_terms ---------------------------------------------------------

#' Balance constraints
#'
#' `balance_terms()` records the set of covariate functions a balancing method
#' should equate across exposure groups. It is passed to [balance()] through
#' `constraints`. The specification is data-free: the covariate expansion it
#' describes is applied at fit time.
#'
#' @details
#' The constraint set is built from four ingredients:
#'
#' - `moments`: the highest power of each numeric covariate to balance. A scalar
#'   applies to every covariate; a named integer vector sets powers per
#'   covariate, with unnamed covariates defaulting to `1`. Powers above `1` are
#'   ignored for binary indicator columns.
#' - `interactions`: when `TRUE`, all pairwise products of distinct base columns
#'   are added, excluding products of two indicators of the same factor.
#' - `quantiles`: probabilities in `(0, 1)`. Each probability adds an indicator
#'   column so that mean balance on the indicator is quantile balance on the
#'   covariate. A single probability vector applies to every continuous
#'   covariate; a named list sets probabilities per covariate. Quantile
#'   constraints apply to discrete exposures only.
#' - `tolerance`: the largest absolute standardized mean difference (discrete
#'   exposures) or exposure-covariate correlation (continuous exposures)
#'   permitted per constraint. A scalar applies to every covariate; a named
#'   vector sets tolerances per source covariate, and derived columns inherit
#'   their source covariate's tolerance. `0` requests exact balance. A positive
#'   value selects the inexact problem for entropy balancing and is the central
#'   tuning parameter for stable balancing weights.
#'
#' @param moments The highest covariate power to balance. A single whole number
#'   or a named integer vector; `NULL` (the default) resolves to first moments.
#' @param interactions Whether to add pairwise interactions of the base columns.
#' @param quantiles Quantile probabilities in `(0, 1)`: a numeric vector applied
#'   to every continuous covariate, a named list of probabilities per covariate,
#'   or `NULL` for none.
#' @param tolerance The per-constraint tolerance: a single non-negative number,
#'   or a named vector giving the tolerance per source covariate.
#' @param ... Reserved for future extensions; must be empty.
#'
#' @return A `balance_terms` object.
#'
#' @examples
#' # Balance means and variances of every numeric covariate.
#' balance_terms(moments = 2)
#'
#' # Balance means with a relaxed tolerance.
#' balance_terms(tolerance = 0.05)
#'
#' @export
balance_terms <- new_class(
  "balance_terms",
  properties = list(
    moments = NULL | class_integer,
    interactions = class_logical,
    quantiles = NULL | class_double | class_list,
    tolerance = class_double
  ),
  constructor = function(
    moments = NULL,
    interactions = FALSE,
    quantiles = NULL,
    tolerance = 0,
    ...
  ) {
    rlang::check_dots_empty()
    if (!is.null(moments)) {
      moments <- vctrs::vec_cast(moments, integer(), x_arg = "moments")
    }
    if (!is.null(quantiles)) {
      if (is.list(quantiles)) {
        quantiles <- lapply(quantiles, function(q) {
          vctrs::vec_cast(q, double(), x_arg = "quantiles")
        })
      } else {
        quantiles <- vctrs::vec_cast(quantiles, double(), x_arg = "quantiles")
      }
    }
    tolerance <- vctrs::vec_cast(tolerance, double(), x_arg = "tolerance")
    new_object(
      S7_object(),
      moments = moments,
      interactions = interactions,
      quantiles = quantiles,
      tolerance = tolerance
    )
  },
  validator = function(self) {
    quantile_values <- if (is.list(self@quantiles)) {
      unlist(self@quantiles, use.names = FALSE)
    } else {
      self@quantiles
    }
    if (!is.null(self@moments) && any(self@moments < 0L)) {
      "@moments must be non-negative"
    } else if (
      !is.null(quantile_values) &&
        length(quantile_values) > 0 &&
        any(quantile_values <= 0 | quantile_values >= 1)
    ) {
      "@quantiles must lie strictly between 0 and 1"
    } else if (any(self@tolerance < 0)) {
      "@tolerance must be non-negative"
    }
  }
)

# ---- Estimating-equations container ----------------------------------------

#' The estimating-equations container
#'
#' `balancing_estimating_equations` holds the pieces a stacked sandwich variance
#' needs after balancing: the fitted parameters, the per-unit estimating
#' functions, the analytic Jacobian, and each unit's weight derivatives. It is
#' produced by methods whose weights solve smooth estimating equations and is
#' `NULL` otherwise.
#'
#' @param parameters The fitted parameters, a numeric vector.
#' @param psi The `n` by `p` estimating functions at the solution.
#' @param jacobian The `p` by `p` analytic Jacobian at the solution.
#' @param weight_jacobian The `n` by `p` weight derivatives.
#' @param weights_raw The balancing weights whose derivative is
#'   `weight_jacobian`, a length-`n` numeric vector. The weight derivatives are
#'   stored at whatever per-group reporting scale a method uses internally, so a
#'   consumer that needs the derivative of the reported weights rescales
#'   `weight_jacobian` by the ratio of the reported weights to `weights_raw`.
#' @param psi_fn An optional function re-evaluating `psi` at new parameters.
#' @param weights_fn An optional function returning the reported balancing
#'   weights at new parameters, a plain double vector with the sampling weights
#'   excluded. At the fitted parameters it reproduces
#'   `as.numeric(weights(fit, include_sampling_weights = FALSE))`. The per-group
#'   reporting scale is fixed at the fit rather than recomputed at each set of
#'   parameters, so the function's derivative is `weight_jacobian` rescaled by
#'   the ratio of the reported weights to `weights_raw`, which is the weight
#'   coupling a stacked variance needs.
#'
#' @return A `balancing_estimating_equations` object.
#' @keywords internal
#' @export
balancing_estimating_equations <- new_class(
  "balancing_estimating_equations",
  properties = list(
    parameters = class_double,
    psi = class_double,
    jacobian = class_double,
    weight_jacobian = class_double,
    weights_raw = NULL | class_double,
    psi_fn = NULL | class_function,
    weights_fn = NULL | class_function
  )
)

# ---- Result class ----------------------------------------------------------

#' A fitted balancing result
#'
#' `balancing` is the S7 object [balance()] returns. It carries the balancing
#' weights, the resolved method and estimand, the covariate expansion recipe, a
#' balance table, and the solver diagnostics. The raw data are not stored; the
#' recipe carries what is needed to rebuild the constraint matrix.
#'
#' @usage NULL
#' @param weights The balancing weights, a [bw] vector.
#' @param method The fitted [balance_method] specification.
#' @param estimand The resolved estimand string.
#' @param exposure The exposure column name.
#' @param exposure_type The resolved exposure type.
#' @param covariates The covariate column names.
#' @param focal_level The focal exposure level for `"att"` and `"atc"`, or
#'   `NULL`.
#' @param n The number of observations.
#' @param constraints The resolved [balance_terms] specification, or `NULL`.
#' @param recipe The covariate expansion recipe, a list of per-column records.
#' @param balance_table The achieved balance, one row per constraint term.
#' @param duals Solver dual variables for diagnostics, or `NULL`.
#' @param coefficients The fitted coefficients or dual variables, or `NULL`.
#' @param converged Whether the solver met its convergence criterion.
#' @param iterations The solver iteration count.
#' @param objective The solved objective value.
#' @param solver_status The solver that produced the result.
#' @param estimating_equations The [balancing_estimating_equations] container,
#'   or `NULL`.
#' @param sampling_weights The sampling weights, or `NULL`.
#' @param call The originating call.
#'
#' @return A `balancing` object.
#' @export
balancing <- new_class(
  "balancing",
  properties = list(
    weights = new_S3_class("bw"),
    method = balance_method,
    estimand = class_character,
    exposure = class_character,
    exposure_type = class_character,
    covariates = class_character,
    focal_level = NULL | class_character,
    n = class_integer,
    constraints = NULL | balance_terms,
    recipe = class_list,
    balance_table = class_data.frame,
    duals = NULL | class_data.frame,
    coefficients = NULL | class_double,
    converged = class_logical,
    iterations = class_integer,
    objective = class_double,
    solver_status = class_character,
    estimating_equations = NULL | balancing_estimating_equations,
    sampling_weights = NULL | class_double,
    call = class_call
  )
)

# ---- Shared result methods -------------------------------------------------

method(print, balancing) <- function(x, ...) {
  cat_cli({
    cli::cli_h1("{method_label(x@method)}")
    cli::cli_text("Exposure: {.val {x@exposure}} ({x@exposure_type})")
    if (!is.null(x@focal_level)) {
      cli::cli_text(
        "Estimand: {.val {x@estimand}} (focal level {.val {x@focal_level}})"
      )
    } else {
      cli::cli_text("Estimand: {.val {x@estimand}}")
    }
    cli::cli_text("Observations: {x@n}")

    status <- if (x@converged) "converged" else "did not converge"
    cli::cli_text(
      "Solver: {status} in {x@iterations} iteration{?s}"
    )

    # An objective-driven method may carry no constraint terms at all, and an
    # empty balance table has neither a tolerance nor a largest imbalance to
    # report. Naming the empty set says so, where a maximum over no terms would
    # print negative infinity behind a warning.
    n_constraints <- nrow(x@balance_table)
    if (n_constraints == 0L) {
      cli::cli_text("Constraints: none")
    } else {
      tol <- max(x@balance_table$tolerance)
      cli::cli_text(
        "Constraints: {n_constraints} term{?s} (tolerance {tol})"
      )

      statistic <- x@balance_table$statistic[1]
      largest <- max(abs(x@balance_table$weighted))
      label <- if (identical(statistic, "correlation")) {
        "correlation"
      } else {
        "standardized mean difference"
      }
      cli::cli_text(
        "Largest imbalance: {formatC(largest, format = 'f', digits = 4)} ({label})"
      )
    }
  })
  invisible(x)
}

method(summary, balancing) <- function(object, ...) {
  print(object)
  w <- as.numeric(weights(object))
  mean_w <- mean(w)
  cv <- stats::sd(w) / mean_w
  # The quadratic-program family reports how many weights rest on the
  # minimum-weight floor. The floor is applied to the reported balancing weights
  # after each group is renormalized, so the count is taken there, on the
  # `@weights` scale rather than the sampling-weight-composed scale, with a tight
  # relative tolerance around the floor value.
  report_floor <- S7_inherits(object@method, quadratic_program_method)
  if (report_floor) {
    floor <- object@method@min_weight
    raw <- as.numeric(object@weights)
    at_floor <- sum(raw <= floor * (1 + 1e-6) + 1e-12)
  }
  cat_cli({
    cli::cli_h2("Weights")
    cli::cli_text(
      "Range: {formatC(min(w), format = 'f', digits = 3)} to {formatC(max(w), format = 'f', digits = 3)}"
    )
    cli::cli_text("Mean: {formatC(mean_w, format = 'f', digits = 3)}")
    cli::cli_text(
      "Coefficient of variation: {formatC(cv, format = 'f', digits = 3)}"
    )
    if (report_floor) {
      cli::cli_text(
        "Weights at the minimum-weight floor: {at_floor} of {object@n}"
      )
    }
    cli::cli_h2("Balance")
  })
  # Print the head as a base data frame so the display does not depend on
  # whether the tibble print method is attached.
  print(utils::head(as.data.frame(object@balance_table)))
  invisible(object)
}

#' Extract the estimating-equations container
#'
#' `estimating_equations()` returns the [balancing_estimating_equations] a fit
#' produced, the pieces a stacked sandwich variance needs after balancing. It is
#' available only for fits whose weights solve smooth estimating equations: the
#' estimating-equation family (entropy balancing, inverse probability tilting,
#' just-identified covariate balancing propensity score) with exact balance.
#' Fits without estimating equations, such as any tolerance-relaxed or
#' quadratic-program fit, raise `balancing_ipw_unsupported_error`.
#'
#' @param x A [balancing] result.
#' @param ... Ignored.
#'
#' @return A [balancing_estimating_equations] object.
#'
#' @examples
#' n <- 200
#' x1 <- rnorm(n)
#' df <- data.frame(exposure = rbinom(n, 1, plogis(0.5 * x1)), x1 = x1)
#' fit <- balance(df, exposure, x1, method = bw_entropy())
#' estimating_equations(fit)
#'
#' @export
estimating_equations <- new_generic("estimating_equations", "x")

method(estimating_equations, balancing) <- function(x, ...) {
  rlang::check_dots_empty()
  if (is.null(x@estimating_equations)) {
    abort(
      c(
        "This fit has no estimating equations.",
        i = "Estimating equations are produced by the estimating-equation family with exact balance.",
        i = "Use the bootstrap workflow in the inference vignette for variance instead."
      ),
      error_class = "balancing_ipw_unsupported_error"
    )
  }
  x@estimating_equations
}
