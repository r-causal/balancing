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
#' @param convergence_tolerance The solver convergence tolerance, or `NULL` to
#'   leave it to the solver. The value resolved for `NULL` differs by family:
#'   `1e-10` on the gradient for the estimating-equation methods, and `1e-8` as
#'   both the absolute and the relative tolerance for the quadratic-program
#'   methods.
#' @param max_iterations The maximum solver iterations, or `NULL` to leave the
#'   cap to the solver. The value resolved for `NULL` is 1000 for the
#'   estimating-equation methods and 200000 for the quadratic-program methods.
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
  # Each solver knob is optional, so a supplied value is checked for being a
  # single usable number rather than only for its sign. A missing value or a
  # vector of the wrong length would otherwise decide the sign comparison and stop
  # the validator itself with a base error, leaving the property unnamed. The
  # tolerance divides and multiplies inside the solve, so an infinity is refused
  # alongside the missing value; the iteration cap is an integer and cannot hold
  # one. Every concrete method inherits these checks rather than repeating them.
  validator = function(self) {
    if (!is.null(self@convergence_tolerance)) {
      if (
        length(self@convergence_tolerance) != 1 ||
          !is.finite(self@convergence_tolerance) ||
          self@convergence_tolerance <= 0
      ) {
        return("@convergence_tolerance must be a single finite positive number")
      }
    }
    if (!is.null(self@max_iterations)) {
      if (
        length(self@max_iterations) != 1 ||
          is.na(self@max_iterations) ||
          self@max_iterations < 0L
      ) {
        return("@max_iterations must be a single non-negative whole number")
      }
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
  ),
  # The penalty and the minimum-weight floor belong to every quadratic program, so
  # they are validated where they are declared rather than in each concrete
  # constructor. The floor has a ceiling as well as a sign: the constraint set pins
  # each reweighted arm's mean weight at one, so a floor at one leaves the uniform
  # weighting as the only feasible point and a floor above one leaves none at all.
  # Both used to reach the solver and return as an infeasibility blamed on the
  # constraints. An infinite penalty is refused for the same reason a missing one
  # is: neither names a quadratic the solver can form.
  validator = function(self) {
    if (
      length(self@weight_penalty) != 1 ||
        !is.finite(self@weight_penalty) ||
        self@weight_penalty < 0
    ) {
      return("@weight_penalty must be a single finite non-negative number")
    }
    if (
      length(self@min_weight) != 1 ||
        is.na(self@min_weight) ||
        self@min_weight < 0
    ) {
      return("@min_weight must be a single non-negative number")
    }
    if (self@min_weight >= 1) {
      return(
        "@min_weight must be less than one, the mean weight the constraint set pins each reweighted arm at"
      )
    }
  }
)

# The validator clause for a `distribution_moments` property, shared by the two
# methods that hold marginal distribution moments for a continuous exposure.
# Both read the property the same way, as a single count of one or more, and both
# used to reach the comparison in their fit path with whatever was supplied, so
# the clause is written once. Returns `NULL` when the value is usable, which is
# what an S7 validator returns to accept an object.
validate_distribution_moments <- function(moments) {
  if (length(moments) != 1 || is.na(moments)) {
    return("@distribution_moments must be a single whole number")
  }
  if (moments < 1L) {
    return("@distribution_moments must be a positive whole number")
  }
  NULL
}

# ---- Capability generics ---------------------------------------------------

#' Method capability generics
#'
#' These generics let [balance()] interrogate a method specification without
#' knowing its concrete class. Each balancing method supplies one method per
#' generic.
#'
#' @param method A [balance_method] specification.
#' @param exposure_type One of `"binary"`, `"categorical"`, or `"continuous"`.
#' @param ... Additional context for `supports_estimating_equations()`, passed by
#'   name. Every method accepts both `exposure_type` and `constraints` and
#'   ignores whichever of the two its answer does not depend on, so one call
#'   shape puts the question to any method.
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
#' A factor covariate contributes one indicator per level rather than the
#' reference coding a model formula would use. Those indicators sum to the
#' constant every balancing method carries, so one of them is redundant and the
#' expansion drops it, naming the term in an informational alert. The level
#' dropped is the last one, and the constraints that remain balance it as well:
#' with the other level proportions equated across exposure groups, the omitted
#' one follows. The balance table reports the surviving levels.
#'
#' @param moments The highest covariate power to balance. A single whole number
#'   or a named integer vector; `NULL` (the default) resolves to first moments.
#' @param interactions Whether to add pairwise interactions of the base columns.
#'   These expand the constraint set the weights must balance, adding the
#'   pairwise products of the base columns to the covariate functions a fit
#'   constrains: equated across the exposure groups of a discrete exposure, and
#'   driven to zero correlation with a continuous one. They say nothing about
#'   causal interaction between two exposures, which is an effect rather than a
#'   constraint and is reported for a joint exposure by
#'   [`ipw()`][ipw.balancing].
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
    check_tolerance(tolerance)
    new_object(
      S7_object(),
      moments = moments,
      interactions = interactions,
      quantiles = quantiles,
      tolerance = tolerance
    )
  },
  # Each range check reads a whole vector, since every one of these properties may
  # carry a value per covariate. A missing value anywhere in one would decide the
  # comparison it reaches and stop the validator with a base error, so missingness
  # is refused first, per property, and the range checks then run on values they
  # can compare.
  validator = function(self) {
    quantile_values <- if (is.list(self@quantiles)) {
      unlist(self@quantiles, use.names = FALSE)
    } else {
      self@quantiles
    }
    if (anyNA(self@moments)) {
      return("@moments must not contain missing values")
    }
    if (anyNA(quantile_values)) {
      return("@quantiles must not contain missing values")
    }
    if (anyNA(self@tolerance)) {
      return("@tolerance must not contain missing values")
    }
    if (!is.null(self@moments) && any(self@moments < 0L)) {
      "@moments must be non-negative"
    } else if (
      !is.null(quantile_values) &&
        length(quantile_values) > 0 &&
        any(quantile_values <= 0 | quantile_values >= 1)
    ) {
      "@quantiles must lie strictly between 0 and 1"
    } else if (!all(is.finite(self@tolerance))) {
      "@tolerance must be finite"
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
#' @param parts_fn An optional function returning both of the above at one set
#'   of parameters, a list with elements `weights` and `psi` holding exactly
#'   what `weights_fn` and `psi_fn` return there. A method whose estimating
#'   functions are a transformation of its own weights computes the pair
#'   together for the price of one, and a consumer that needs both at every
#'   parameter vector, as a stacked sandwich does, halves its work by asking for
#'   them together. It is an optimization rather than a contract: a consumer
#'   reads it when it is present and falls back to the two functions when it is
#'   not, so a method supplies it only when the saving is real.
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
    weights_fn = NULL | class_function,
    parts_fn = NULL | class_function
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
#' @param exposure_levels The exposure levels the fit weighted, as character
#'   strings in the order `levels(factor(exposure))` gives: a factor's declared
#'   order for a factor exposure, and the values sorted in their own type
#'   otherwise, so a numeric dose of 9 and 10 orders 9 first. Levels no
#'   observation takes are dropped and do not appear. The first element is the
#'   reference level every contrast in [`ipw()`][ipw.balancing] is measured
#'   against. A continuous exposure carries no levels, so the vector is empty.
#' @param covariates The covariate column names that contributed at least one
#'   retained constraint column, in selection order. The constraint expansion
#'   drops constant and aliased columns, and a covariate can be selected without
#'   contributing any column at all, so this records what the fit constrained
#'   rather than what was requested; the request itself stays in `call`. A fit
#'   whose balance comes from its objective rather than from constraints, such as
#'   energy or kernel balancing with no moment constraints, therefore records no
#'   covariates even though its objective reads every selected one.
#' @param focal_level The focal exposure level for `"att"` and `"atc"`, or
#'   `NULL`. This is the fitted object's property, set from the `.focal_level`
#'   argument of [balance()].
#' @param n The number of observations.
#' @param constraints The resolved [balance_terms] specification, or `NULL`.
#' @param recipe The covariate expansion recipe, a list of per-column records.
#' @param balance_table The achieved balance, one row per constraint term.
#' @param duals Solver dual variables for diagnostics, or `NULL`.
#' @param coefficients The fitted coefficients or dual variables, or `NULL`.
#' @param converged Whether the solver met its convergence criterion.
#' @param iterations The solver iteration count. An energy fit that could not
#'   reach its tolerance re-solves at a reachable one, and when that re-solve
#'   converges this sums the original and the fallback solve, so it can exceed
#'   the requested `max_iterations`. When the re-solve does not converge the
#'   fit reports the original solve alone, so the count stays within the cap.
#'   See [bw_energy()] for the fuller account.
#' @param objective The solved objective value.
#' @param solver_status The solver that produced the result.
#' @param estimating_equations The [balancing_estimating_equations] container,
#'   or `NULL`.
#' @param vcov The covariance of the fit's own weight parameters, a `p` by `p`
#'   matrix named for them, or `NULL`. A covariance for those parameters comes
#'   from the stacked system an [`ipw()`][ipw.balancing] result assembles, so
#'   [balance()] leaves this empty and the result fills it in on the copy of the
#'   fit it stores, where [stats::vcov()] reads it back.
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
    exposure_levels = class_character,
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
    vcov = NULL | class_double,
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

      # The largest imbalance names how far the fit missed, so it is only offered
      # when it is a number. A maximum that is not finite means at least one
      # constraint's statistic is undefined, and printing "NaN" as the largest
      # imbalance states a distance that was never measured. Saying the
      # assessment is what failed is the same ruling the balance warning follows.
      # The figure is rendered to three significant digits, matching the balance
      # warning, because a well-solved fit can leave an imbalance far below the
      # fourth decimal and a fixed format prints that as a zero, contradicting
      # the warning that names the same number.
      statistic <- x@balance_table$statistic[1]
      largest <- max(abs(x@balance_table$weighted))
      if (is.finite(largest)) {
        label <- if (identical(statistic, "correlation")) {
          "correlation"
        } else {
          "standardized mean difference"
        }
        cli::cli_text(
          "Largest imbalance: {formatC(largest, format = 'g', digits = 3)} ({label})"
        )
      } else {
        cli::cli_text("Largest imbalance: could not be assessed")
      }
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
