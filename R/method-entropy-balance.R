# Entropy balancing: the maximum-entropy reweighting whose weighted covariate
# moments match a set of targets. The method specification carries only tuning
# parameters; fit_method() prepares the constraint targets, calls the Rust
# solver, applies the estimand's group-sum normalization, and assembles the
# smooth estimating equations when the problem is exact.

#' Entropy balancing
#'
#' `entropy_balance()` specifies entropy balancing for [balance()]. The weights
#' minimize the Kullback-Leibler divergence from a set of base weights subject to
#' the covariate constraints, so among all reweightings that achieve balance the
#' solution stays as close as possible to the base weights. Entropy balancing
#' supports binary, categorical, and continuous exposures.
#'
#' @details
#' For a binary exposure the average treatment effect reweights each exposure
#' group to the pooled covariate means, and the average treatment effect on the
#' treated reweights the control group to the treated covariate means while the
#' treated group keeps its base weights. When every requested tolerance is zero
#' the constraints hold exactly and the weights solve smooth estimating
#' equations, which [balance()] records for later linearized inference. A
#' positive `tolerance` in [balance_terms()] selects the inexact problem, which
#' balances each constraint to within the tolerance and does not produce
#' estimating equations.
#'
#' @param base_weights A numeric vector of base weights, one per observation, or
#'   `NULL` for uniform base weights. The estimated weights minimize
#'   `sum(w * log(w / base_weights))`.
#' @param distribution_moments For continuous exposures, the number of exposure
#'   and covariate marginal moments held equal to the unweighted sample, or
#'   `NULL` for the constraint moments. Raised automatically when smaller than
#'   the constraint moments.
#' @param convergence_tolerance The solver convergence tolerance on the gradient.
#' @param max_iterations The maximum solver iterations, or `NULL` for the core
#'   default.
#' @param ... Reserved for future extensions; must be empty. Tuning parameters
#'   must be passed by name.
#'
#' @return An `entropy_balance` specification, a [balance_method].
#'
#' @references
#' Hainmueller, J. (2012). Entropy balancing for causal effects: A multivariate
#' reweighting method to produce balanced samples in observational studies.
#' *Political Analysis*, 20(1), 25-46.
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
#' fit <- balance(df, exposure, c(x1, x2), method = entropy_balance())
#' fit
#'
#' @export
entropy_balance <- new_class(
  "entropy_balance",
  parent = estimating_equation_method,
  properties = list(
    base_weights = NULL | class_double,
    distribution_moments = NULL | class_integer
  ),
  constructor = function(
    ...,
    base_weights = NULL,
    distribution_moments = NULL,
    convergence_tolerance = 1e-10,
    max_iterations = NULL
  ) {
    rlang::check_dots_empty()
    if (!is.null(base_weights)) {
      base_weights <- vctrs::vec_cast(
        base_weights,
        double(),
        x_arg = "base_weights"
      )
    }
    if (!is.null(distribution_moments)) {
      distribution_moments <- vctrs::vec_cast(
        distribution_moments,
        integer(),
        x_arg = "distribution_moments"
      )
    }
    new_object(
      S7_object(),
      base_weights = base_weights,
      distribution_moments = distribution_moments,
      convergence_tolerance = convergence_tolerance,
      max_iterations = max_iterations
    )
  },
  validator = function(self) {
    if (!is.null(self@base_weights) && any(self@base_weights < 0)) {
      "@base_weights must be non-negative"
    }
  }
)

method(method_label, entropy_balance) <- function(method) {
  "Entropy balancing"
}

method(supported_exposure_types, entropy_balance) <- function(method) {
  c("binary", "categorical", "continuous")
}

method(supported_estimands, entropy_balance) <- function(
  method,
  exposure_type
) {
  switch(
    exposure_type,
    binary = c("ate", "att", "atu"),
    categorical = c("ate", "att"),
    continuous = "ate"
  )
}

method(supports_estimating_equations, entropy_balance) <- function(
  method,
  ...,
  constraints = NULL
) {
  rlang::check_dots_empty()
  if (is.null(constraints)) {
    return(TRUE)
  }
  all(constraints@tolerance <= 0)
}

# Assemble the Rust option list, dropping the tuning parameters left at the core
# default so the solver applies its own. The worker-thread count is resolved on
# the R side and passed down on every call. The inexact problem is solved by
# FISTA, whose termination measures the relative change in the loss rather than
# the gradient; that criterion is weaker than the Newton gradient tolerance, so
# the convergence tolerance is tightened for the inexact path to keep the
# achieved balance inside the requested box.
entropy_options <- function(method, inexact = FALSE) {
  options <- list(threads = resolve_threads())
  tolerance <- method@convergence_tolerance
  if (inexact) {
    tolerance <- min(tolerance %||% 1e-10, 1e-14)
  }
  if (!is.null(tolerance)) {
    options$convergence_tolerance <- tolerance
  }
  if (!is.null(method@max_iterations)) {
    options$max_iterations <- as.integer(method@max_iterations)
  }
  options
}

# Scale each column's SMD-scale tolerance to the raw scale the solver's box
# constrains. Numeric columns cross standardized to unit scale, so their box
# equals the SMD tolerance; indicator and quantile columns cross raw, so their
# box must be multiplied by the column's standard deviation for the achieved
# standardized mean difference to bind at the requested tolerance.
solver_box <- function(z, tolerances) {
  column_sd <- apply(z, 2, stats::sd)
  column_sd[column_sd == 0] <- 1
  tolerances * column_sd
}

# Weighted column means of Z over the rows in `idx`, the constraint targets.
target_means <- function(z, idx, s) {
  weights <- s[idx] / sum(s[idx])
  as.numeric(crossprod(z[idx, , drop = FALSE], weights))
}

method(fit_method, entropy_balance) <- function(method, prepared) {
  if (identical(prepared$exposure_type, "continuous")) {
    fit_entropy_continuous(method, prepared)
  } else {
    fit_entropy_discrete(method, prepared)
  }
}

fit_entropy_discrete <- function(method, prepared) {
  z <- prepared$matrix
  n <- prepared$n
  p <- ncol(z)
  s <- prepared$sampling_weights
  base <- method@base_weights %||% rep(1, n)
  if (length(base) != n) {
    abort(
      c(
        "{.arg base_weights} must have one value per observation.",
        x = "It has length {length(base)}, but the data have {n} row{?s}."
      ),
      error_class = "balancing_range_error"
    )
  }
  tolerances <- prepared$tolerances
  inexact <- any(tolerances > 0)
  tols <- solver_box(z, tolerances)
  levels <- prepared$exposure_levels
  groups <- prepared$groups
  estimand <- prepared$estimand
  focal <- prepared$focal_level

  # The reference distribution the tilt anchors to is the base measure, so that
  # units already balanced under their base weights are left proportional to
  # them.
  measure <- s * base

  if (identical(estimand, "ate")) {
    targets <- target_means(z, seq_len(n), measure)
    solved_levels <- levels
    group_idx <- match(prepared$exposure_key, solved_levels) - 1L
    n_eff <- 1
  } else {
    focal_idx <- groups[[focal]]
    targets <- target_means(z, focal_idx, measure)
    solved_levels <- setdiff(levels, focal)
    group_idx <- rep(-1L, n)
    for (block in seq_along(solved_levels)) {
      group_idx[groups[[solved_levels[block]]]] <- block - 1L
    }
    n_eff <- sum(s[focal_idx])
  }

  result <- solve_entropy(
    z,
    as.integer(group_idx),
    targets,
    base,
    s,
    tols,
    n_eff,
    entropy_options(method, inexact = inexact)
  )

  weights_rust <- result$weights
  w <- weights_rust
  if (!identical(estimand, "ate")) {
    w[focal_idx] <- base[focal_idx]
  }

  for (level in levels) {
    idx <- groups[[level]]
    target_sum <- if (identical(estimand, "ate")) sum(s[idx]) else n_eff
    current <- sum(s[idx] * w[idx])
    if (current > 0) {
      w[idx] <- w[idx] * (target_sum / current)
    }
  }

  solved_groups <- lapply(solved_levels, function(level) groups[[level]])
  block_scales <- group_scale_factors(solved_groups, s, w, weights_rust)

  list(
    weights = w,
    coefficients = as.numeric(result$duals),
    converged = isTRUE(result$converged),
    interrupted = isTRUE(result$interrupted),
    iterations = as.integer(result$iterations),
    objective = entropy_objective(s, w, base),
    solver_status = result$solver,
    estimating_equations = rescale_estimating_equations(
      result,
      ncol(z),
      block_scales
    ),
    groups = groups
  )
}

fit_entropy_continuous <- function(method, prepared) {
  z <- prepared$matrix
  n <- prepared$n
  s <- prepared$sampling_weights
  base <- method@base_weights %||% rep(1, n)
  tolerances <- prepared$tolerances

  exposure <- as.numeric(prepared$exposure_vec)
  exposure_scale <- stats::sd(exposure)
  if (exposure_scale == 0) {
    exposure_scale <- 1
  }
  e <- (exposure - mean(exposure)) / exposure_scale

  # The marginal distribution constraints hold the exposure and covariate
  # marginals equal to the unweighted sample. The correlation constraints drive
  # each weighted exposure-covariate product to zero. The distribution moments
  # extend the marginals: they are raised to at least the constraint moments,
  # with an alert when the requested value is smaller.
  covariate_moments <- covariate_constraint_moments(
    prepared$recipe,
    prepared$covariates
  )
  constraint_moments <- max(1L, max(covariate_moments, 0L))
  moments <- resolve_distribution_moments(
    method@distribution_moments,
    constraint_moments
  )

  exposure_marginals <- moment_columns(exposure, moments)
  extra_marginals <- higher_covariate_marginals(
    prepared$data,
    covariate_moments,
    moments
  )
  marginals <- do.call(cbind, c(list(exposure_marginals, z), extra_marginals))
  products <- z * e

  covs <- cbind(marginals, products)
  n_marginal <- ncol(marginals)
  n_product <- ncol(products)
  dist_ind <- c(rep(1L, n_marginal), rep(0L, n_product))
  targets <- rep(0, ncol(covs))

  # The exposure and each covariate column cross standardized and their
  # marginals are held to unit variance, so the weighted mean of their product
  # is the weighted exposure-covariate correlation. The tolerance is therefore
  # applied to the product columns directly, on the correlation scale the design
  # specifies, without the standard-deviation rescaling the discrete indicators
  # need.
  tols <- c(rep(0, n_marginal), tolerances)
  inexact <- any(tolerances > 0)
  n_eff <- sum(s)

  result <- solve_entropy_cont(
    covs,
    targets,
    tols,
    as.integer(dist_ind),
    base,
    s,
    n_eff,
    entropy_options(method, inexact = inexact)
  )

  weights_rust <- result$weights
  w <- weights_rust
  current <- sum(s * w)
  if (current > 0) {
    w <- w * (n_eff / current)
  }

  block_scales <- group_scale_factors(list(seq_len(n)), s, w, weights_rust)

  list(
    weights = w,
    coefficients = as.numeric(result$duals),
    converged = isTRUE(result$converged),
    interrupted = isTRUE(result$interrupted),
    iterations = as.integer(result$iterations),
    objective = entropy_objective(s, w, base),
    solver_status = result$solver,
    estimating_equations = rescale_estimating_equations(
      result,
      ncol(covs),
      block_scales
    ),
    groups = NULL
  )
}

# The entropy objective is the achieved Kullback-Leibler divergence from the base
# weights, summed over the reweighted units.
entropy_objective <- function(s, w, base) {
  positive <- w > 0 & base > 0
  sum(s[positive] * w[positive] * log(w[positive] / base[positive]))
}

# The per-group scalar relating the reported weights to the weights the solver
# returned. The R layer renormalizes each group's weights to the estimand's
# target sum after the solve; that scalar multiplies the estimating-equation
# rows and Jacobian blocks the core computed at its own normalization.
group_scale_factors <- function(solved_groups, s, w, weights_rust) {
  vapply(
    solved_groups,
    function(idx) {
      denominator <- sum(s[idx] * weights_rust[idx])
      if (denominator == 0) {
        1
      } else {
        sum(s[idx] * w[idx]) / denominator
      }
    },
    numeric(1)
  )
}

# Build the estimating-equations container from the matrices the core returned,
# applying the per-group rescaling. The core computes the per-unit estimating
# functions, the Jacobian, and the weight derivatives at the solution; the R
# layer performs only the generic per-block scalar multiply and never
# re-derives the method's moment conditions. Returns `NULL` for the inexact
# problem, whose core fields are absent.
rescale_estimating_equations <- function(result, p, block_scales) {
  psi <- result$psi
  jacobian <- result$jac
  weight_jacobian <- result$dw_dbeta
  if (is.null(psi) || is.null(jacobian) || is.null(weight_jacobian)) {
    return(NULL)
  }
  for (block in seq_along(block_scales)) {
    scale <- block_scales[[block]]
    if (scale == 1) {
      next
    }
    columns <- (block - 1L) * p + seq_len(p)
    psi[, columns] <- psi[, columns] * scale
    weight_jacobian[, columns] <- weight_jacobian[, columns] * scale
    jacobian[columns, columns] <- jacobian[columns, columns] * scale
  }
  balancing_estimating_equations(
    parameters = as.numeric(result$duals),
    psi = psi,
    jacobian = jacobian,
    weight_jacobian = weight_jacobian,
    psi_fn = NULL
  )
}

# The highest balanced power of each numeric covariate, read from the recipe.
# Covariates that contribute only indicator columns report zero.
covariate_constraint_moments <- function(recipe, covariates) {
  moments <- stats::setNames(rep(0L, length(covariates)), covariates)
  for (record in recipe) {
    if (identical(record$type, "numeric")) {
      source <- record$source
      moments[[source]] <- max(moments[[source]], as.integer(record$power))
    }
  }
  moments
}

# Resolve the distribution moments to at least the constraint moments, alerting
# when a smaller value was requested.
resolve_distribution_moments <- function(requested, constraint_moments) {
  if (is.null(requested)) {
    return(constraint_moments)
  }
  requested <- as.integer(requested)
  if (requested < constraint_moments) {
    alert_info(
      "Raising {.arg distribution_moments} to the constraint moments ({constraint_moments})."
    )
    return(constraint_moments)
  }
  requested
}

# Standardized centered powers 1..moments of a numeric vector, the marginal
# moment columns for a distribution constraint. Each column is centered to mean
# zero and scaled to unit standard deviation, so the target is the unweighted
# sample value of zero.
moment_columns <- function(x, moments) {
  centered <- x - mean(x)
  columns <- lapply(seq_len(moments), function(power) {
    standardize_vector(centered^power)
  })
  do.call(cbind, columns)
}

# Marginal moment columns for the covariate powers above those already in the
# constraint matrix, up to the requested moments. Only numeric covariates
# contribute; a covariate whose constraint moment already reaches the requested
# moments adds nothing.
higher_covariate_marginals <- function(data, covariate_moments, moments) {
  columns <- list()
  for (cov in names(covariate_moments)) {
    start <- covariate_moments[[cov]] + 1L
    if (covariate_moments[[cov]] < 1L || start > moments) {
      next
    }
    values <- as.numeric(data[[cov]])
    centered <- values - mean(values)
    for (power in start:moments) {
      columns[[length(columns) + 1]] <- standardize_vector(centered^power)
    }
  }
  columns
}

# Center a vector to mean zero and scale to unit standard deviation, guarding a
# degenerate column.
standardize_vector <- function(x) {
  scale <- stats::sd(x)
  if (scale == 0) {
    scale <- 1
  }
  (x - mean(x)) / scale
}
