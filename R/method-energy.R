# Energy balancing: weights that minimize an energy statistic of covariate
# balance. For a binary or categorical exposure the objective is the energy
# distance between each reweighted exposure group and the target sample, with the
# improved variant adding the between-group energy distance for the average
# treatment effect; for a continuous exposure the objective is the weighted
# distance covariance between the exposure and the covariates. Every form is a
# quadratic program with a simplex-type constraint set, so the method spec
# carries only tuning parameters and fit_method() assembles the covariate
# distance inputs, calls the Rust solver, renormalizes each group to its estimand
# target total, and reports the solver's dual variables for diagnostics. The
# quadratic-program family has no estimating equations.

#' Energy balancing
#'
#' `bw_energy()` specifies energy balancing for [balance()]. The weights
#' minimize the energy distance between the reweighted exposure groups and a
#' target sample, subject to a simplex-type constraint set, so the reweighting
#' improves multivariate covariate balance without positing a propensity model.
#' Energy balancing supports binary, categorical, and continuous exposures.
#'
#' @details
#' For a binary or categorical exposure the objective is the sum of each group's
#' energy distance to the target sample. The improved variant for the average
#' treatment effect adds the between-group energy distance, which balances the
#' groups against one another as well as against the sample. A focal estimand
#' reweights the non-focal groups toward the focal group, whose units keep their
#' base weight. The energy distance is built from a pairwise covariate distance
#' matrix; `distance` selects how that matrix is formed.
#'
#' For a continuous exposure the objective is the weighted distance covariance
#' between the exposure and the covariates, following Huling, Greifer, and Chen,
#' plus the marginal energy distances of the weighted exposure and covariate
#' distributions. `distribution_moments` sets how many exposure and covariate
#' marginal moments are held equal to the unweighted sample, and
#' `dimension_adjustment` reweights the covariate energy distance by the
#' covariate dimensionality.
#'
#' Energy balancing belongs to the quadratic-program family, which has no
#' estimating equations, so a fit produces no linearized-inference container and
#' the tolerance in [balance_terms()] relaxes any added moment constraints rather
#' than selecting an inexact solver. A tolerance supplied without moment
#' constraints has nothing to relax, so it is warned and ignored.
#'
#' @param distance The covariate distance definition the energy objective is
#'   built on, one of `"scaled_euclidean"` (each covariate divided by its
#'   standard deviation), `"mahalanobis"`, or `"euclidean"`.
#' @param improved Whether to add the between-group energy distance of the
#'   improved variant for the average treatment effect with a discrete exposure.
#' @param weight_penalty The L2 penalty on the weights, which stabilizes the
#'   quadratic program.
#' @param min_weight The smallest permitted weight.
#' @param distribution_moments For a continuous exposure, the number of exposure
#'   and covariate marginal moments held equal to the unweighted sample, or
#'   `NULL` for the constraint moments. Raised automatically when smaller than
#'   the constraint moments.
#' @param dimension_adjustment For a continuous exposure, whether to weight the
#'   covariate energy distance by the covariate dimensionality adjustment.
#' @param convergence_tolerance The quadratic-program solver tolerance, or `NULL`
#'   for the core default.
#' @param max_iterations The maximum solver iterations, or `NULL` for the core
#'   default.
#' @param ... Reserved for future extensions; must be empty. Tuning parameters
#'   must be passed by name.
#'
#' @return An `bw_energy` specification, a [balance_method].
#'
#' @references
#' Huling, J. D. and Mak, S. (2024). Energy balancing of covariate distributions.
#' *Journal of Causal Inference*, 12(1), 20220029.
#'
#' Huling, J. D., Greifer, N., and Chen, G. (2024). Independence weights for
#' causal inference with continuous treatments. *Journal of the American
#' Statistical Association*, 119(546), 1657-1670.
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
#' fit <- balance(df, exposure, c(x1, x2), method = bw_energy())
#' fit
#'
#' @export
bw_energy <- new_class(
  "bw_energy",
  parent = quadratic_program_method,
  properties = list(
    distance = class_character,
    improved = class_logical,
    distribution_moments = NULL | class_integer,
    dimension_adjustment = class_logical
  ),
  constructor = function(
    distance = c("scaled_euclidean", "mahalanobis", "euclidean"),
    improved = TRUE,
    weight_penalty = 1e-4,
    min_weight = 1e-8,
    distribution_moments = NULL,
    dimension_adjustment = TRUE,
    convergence_tolerance = NULL,
    max_iterations = NULL,
    ...
  ) {
    rlang::check_dots_empty()
    distance <- rlang::arg_match(distance)
    weight_penalty <- vctrs::vec_cast(
      weight_penalty,
      double(),
      x_arg = "weight_penalty"
    )
    min_weight <- vctrs::vec_cast(min_weight, double(), x_arg = "min_weight")
    if (!is.null(distribution_moments)) {
      distribution_moments <- vctrs::vec_cast(
        distribution_moments,
        integer(),
        x_arg = "distribution_moments"
      )
    }
    if (!is.null(convergence_tolerance)) {
      convergence_tolerance <- vctrs::vec_cast(
        convergence_tolerance,
        double(),
        x_arg = "convergence_tolerance"
      )
    }
    if (!is.null(max_iterations)) {
      max_iterations <- vctrs::vec_cast(
        max_iterations,
        integer(),
        x_arg = "max_iterations"
      )
    }
    new_object(
      S7_object(),
      distance = distance,
      improved = improved,
      weight_penalty = weight_penalty,
      min_weight = min_weight,
      distribution_moments = distribution_moments,
      dimension_adjustment = dimension_adjustment,
      convergence_tolerance = convergence_tolerance,
      max_iterations = max_iterations
    )
  },
  validator = function(self) {
    if (
      length(self@weight_penalty) != 1 ||
        is.na(self@weight_penalty) ||
        self@weight_penalty < 0
    ) {
      return("@weight_penalty must be a single non-negative number")
    }
    if (
      length(self@min_weight) != 1 ||
        is.na(self@min_weight) ||
        self@min_weight < 0
    ) {
      return("@min_weight must be a single non-negative number")
    }
    if (length(self@improved) != 1 || is.na(self@improved)) {
      return("@improved must be a single logical value")
    }
    if (
      length(self@dimension_adjustment) != 1 ||
        is.na(self@dimension_adjustment)
    ) {
      return("@dimension_adjustment must be a single logical value")
    }
    if (!is.null(self@distribution_moments) && self@distribution_moments < 1L) {
      return("@distribution_moments must be a positive whole number")
    }
  }
)

method(method_label, bw_energy) <- function(method) {
  "Energy balancing"
}

method(supported_exposure_types, bw_energy) <- function(method) {
  c("binary", "categorical", "continuous")
}

method(supported_estimands, bw_energy) <- function(method, exposure_type) {
  switch(
    exposure_type,
    # The overlap estimand is legal only for the covariate balancing propensity
    # score, so energy balancing offers the same set as entropy balancing.
    binary = c("ate", "att", "atu"),
    categorical = c("ate", "att"),
    continuous = "ate"
  )
}

# The quadratic-program family has no estimating equations for any exposure type
# or constraint set, so the answer is always FALSE. The context arguments are
# accepted so the generic call shape matches the estimating-equation family.
method(supports_estimating_equations, bw_energy) <- function(
  method,
  ...,
  exposure_type = NULL,
  constraints = NULL
) {
  rlang::check_dots_empty()
  FALSE
}

# Assemble the Rust option list, dropping the tuning parameters left at the core
# default so the quadratic-program solver applies its own. The worker-thread
# count is resolved on the R side and passed down on every call.
energy_options <- function(method) {
  options <- list(threads = resolve_threads())
  if (!is.null(method@convergence_tolerance)) {
    options$convergence_tolerance <- method@convergence_tolerance
  }
  if (!is.null(method@max_iterations)) {
    options$max_iterations <- as.integer(method@max_iterations)
  }
  options
}

# Whether the constraint set requests moment, quantile, or interaction balance,
# which energy balancing enforces as quadratic-program moment rows. For the
# quadratic-program family a bare tolerance requests no constraints of its own,
# and the default constraint set carries none, so the energy objective alone
# drives balance unless one of these is present.
requests_moments <- function(constraints) {
  !is.null(constraints) &&
    (!is.null(constraints@moments) ||
      !is.null(constraints@quantiles) ||
      isTRUE(constraints@interactions))
}

# Whether any positive tolerance is present, which the family relaxes moment
# constraints with when constraints are present and warns about otherwise.
has_positive_tolerance <- function(constraints) {
  !is.null(constraints) && any(constraints@tolerance > 0)
}

# The numeric covariate matrix the distance is built from. Numeric and logical
# covariates cross as themselves; a factor or character covariate contributes one
# indicator column per level, since each level is a coordinate of the covariate
# distance. The Rust core forms the pairwise distance from these columns under
# the named distance definition.
distance_covariates <- function(data, covariates) {
  columns <- lapply(covariates, function(covariate) {
    values <- data[[covariate]]
    if (is.numeric(values) || is.logical(values)) {
      matrix(as.numeric(values), ncol = 1)
    } else {
      levels <- sort(unique(as.character(values)))
      do.call(
        cbind,
        lapply(levels, function(level) {
          as.numeric(as.character(values) == level)
        })
      )
    }
  })
  do.call(cbind, columns)
}

# The tolerance in a balance_terms() specification relaxes added moment
# constraints; with none present it has nothing to act on, so warn and proceed
# with the pure energy objective. A continuous fit holds its distribution moments
# exactly as identifying conditions and never adds relaxable constraints, so any
# positive tolerance is ignored there as well.
warn_ignored_tolerance <- function(call = rlang::caller_env()) {
  warn(
    c(
      "{.arg tolerance} relaxes added moment constraints, but this fit has none to relax.",
      i = "Drop {.arg tolerance} from {.fn balance_terms}, or add moment constraints with {.arg moments}, {.arg quantiles}, or {.arg interactions} for a discrete exposure."
    ),
    warning_class = "balancing_ignored_argument_warning",
    call = call
  )
}

method(fit_method, bw_energy) <- function(method, prepared) {
  if (identical(prepared$exposure_type, "continuous")) {
    if (has_positive_tolerance(prepared$constraints)) {
      warn_ignored_tolerance()
    }
    return(fit_energy_continuous(method, prepared))
  }

  enforce <- requests_moments(prepared$constraints)
  if (!enforce && has_positive_tolerance(prepared$constraints)) {
    warn_ignored_tolerance()
  }

  fit_energy_discrete(method, prepared, enforce)
}

fit_energy_discrete <- function(method, prepared, enforce) {
  n <- prepared$n
  s <- prepared$sampling_weights
  covs <- distance_covariates(prepared$data, prepared$covariates)
  z <- prepared$matrix
  levels <- prepared$exposure_levels
  groups <- prepared$groups
  estimand <- prepared$estimand
  focal <- prepared$focal_level
  key <- prepared$exposure_key

  # The moment-constraint rows hold each reweighted group's weighted mean of the
  # standardized covariate columns within a tolerance band of a target. For the
  # average treatment effect the target is the pooled sample mean and the band is
  # split across the groups so their difference stays within the tolerance; for a
  # focal estimand the target is the focal group's mean and the non-focal groups
  # are pulled to it within the full tolerance. Numeric columns cross at unit
  # scale, so their band is the tolerance directly; indicator columns cross raw,
  # so their band scales by the column standard deviation. The energy objective
  # alone drives balance when no moment constraints are requested.
  measure <- s
  if (enforce) {
    if (identical(estimand, "ate")) {
      targets <- target_means(z, seq_len(n), measure)
    } else {
      targets <- target_means(z, groups[[focal]], measure)
    }
    tols <- solver_box(z, prepared$tolerances)
    moment_covs <- z
  } else {
    targets <- numeric(0)
    tols <- numeric(0)
    moment_covs <- matrix(numeric(0), nrow = n, ncol = 0)
  }

  options <- energy_options(method)

  if (identical(prepared$exposure_type, "binary")) {
    # The focal group is coded one and held fixed; the average treatment effect
    # codes the second level as the reweighting reference. Both focal estimands
    # map to the core's focal solve, which holds the coded-one group.
    if (identical(estimand, "ate")) {
      treat <- as.integer(key == levels[[2]])
      core_estimand <- "ate"
      nvar <- n
      n_group_levels <- 2L
    } else {
      treat <- as.integer(key == focal)
      core_estimand <- "att"
      nvar <- n - length(groups[[focal]])
      n_group_levels <- 1L
    }
    result <- solve_energy(
      covs,
      treat,
      s,
      method@distance,
      core_estimand,
      method@improved,
      moment_covs,
      targets,
      tols,
      method@min_weight,
      method@weight_penalty,
      options
    )
  } else {
    treat_idx <- match(key, levels) - 1L
    if (identical(estimand, "ate")) {
      focal_idx <- 0L
      core_estimand <- "ate"
      nvar <- n
      n_group_levels <- length(levels)
    } else {
      focal_idx <- match(focal, levels) - 1L
      core_estimand <- "att"
      nvar <- n - length(groups[[focal]])
      n_group_levels <- length(levels) - 1L
    }
    result <- solve_energy_multi(
      covs,
      as.integer(treat_idx),
      as.integer(focal_idx),
      s,
      method@distance,
      core_estimand,
      method@improved,
      moment_covs,
      targets,
      tols,
      method@min_weight,
      method@weight_penalty,
      options
    )
  }

  duals <- energy_duals_frame(result$duals, nvar, n_group_levels)
  assemble_energy(result, method, prepared, duals, approximate = !enforce)
}

fit_energy_continuous <- function(method, prepared) {
  n <- prepared$n
  s <- prepared$sampling_weights
  covs <- distance_covariates(prepared$data, prepared$covariates)
  z <- prepared$matrix
  exposure <- as.numeric(prepared$exposure_vec)

  # The distribution-moment constraints hold the weighted exposure and covariate
  # marginals equal to the unweighted sample. They are raised to at least the
  # constraint moments, with an alert when the requested value is smaller. The
  # weighted distance covariance the objective minimizes is what drives the
  # exposure-covariate association toward zero, so no separate correlation
  # constraint is added in the default fit.
  covariate_moments <- covariate_constraint_moments(
    prepared$recipe,
    prepared$covariates
  )
  constraint_moments <- max(1L, max(covariate_moments, 0L))
  moments <- resolve_distribution_moments(
    method@distribution_moments,
    constraint_moments
  )

  d_treat <- moment_columns(exposure, moments)
  extra_covariate_marginals <- higher_covariate_marginals(
    prepared$data,
    covariate_moments,
    moments
  )
  d_covs <- do.call(cbind, c(list(z), extra_covariate_marginals))

  bal_covs <- matrix(numeric(0), nrow = n, ncol = 0)
  bal_tols <- numeric(0)

  options <- energy_options(method)

  result <- solve_energy_cont(
    covs,
    exposure,
    s,
    method@distance,
    method@dimension_adjustment,
    method@min_weight,
    method@weight_penalty,
    d_covs,
    d_treat,
    bal_covs,
    bal_tols,
    options
  )

  # The continuous solve carries one total-sum row followed by the distribution
  # rows; the box rows bound each of the n units.
  n_structural_leading <- 1L
  duals <- energy_duals_frame(result$duals, n, n_structural_leading)
  assemble_energy(result, method, prepared, duals, approximate = TRUE)
}

# The solver's dual variables for the structural constraint rows, dropping the
# per-unit box duals. The quadratic program stacks the box rows first, then the
# group-sum rows, then the moment rows, so the structural duals follow the `nvar`
# box rows. Returns `NULL` when no structural row resolved.
energy_duals_frame <- function(duals, nvar, n_group_rows) {
  duals <- as.numeric(duals)
  m <- length(duals)
  if (m <= nvar) {
    return(NULL)
  }
  structural <- duals[(nvar + 1):m]
  n_structural <- length(structural)
  n_group_rows <- min(n_group_rows, n_structural)
  kind <- c(
    rep("group", n_group_rows),
    rep("moment", n_structural - n_group_rows)
  )
  new_balancing_tibble(list(
    constraint = seq_len(n_structural),
    kind = kind,
    dual = structural
  ))
}

# Renormalize each exposure group to its estimand target total and pack the fit
# result. The core returns weights on the per-group quadratic-program scale and
# applies the minimum-weight floor there; the reported weights renormalize each
# group to its estimand target total, a change of reporting convention, and the
# floor is re-applied on that reported scale so the documented minimum holds
# where the specs assert it. For the average treatment effect each group is
# scaled to its own sampling-weighted total; for a focal estimand every group is
# scaled to the focal total, which leaves the focal group at its base weight. A
# continuous fit is a single group scaled to the sampling-weight total. The
# quadratic-program family carries no estimating equations.
assemble_energy <- function(result, method, prepared, duals, approximate) {
  s <- prepared$sampling_weights
  estimand <- prepared$estimand
  focal <- prepared$focal_level
  groups <- prepared$groups
  levels <- prepared$exposure_levels

  w <- as.numeric(result$weights)
  if (is.null(groups)) {
    current <- sum(s * w)
    if (current > 0) {
      w <- w * (sum(s) / current)
    }
  } else {
    focal_estimand <- estimand %in% c("att", "atu")
    for (level in levels) {
      idx <- groups[[level]]
      target_sum <- if (focal_estimand) {
        sum(s[groups[[focal]]])
      } else {
        sum(s[idx])
      }
      current <- sum(s[idx] * w[idx])
      if (current > 0) {
        w[idx] <- w[idx] * (target_sum / current)
      }
    }
  }
  w[w < method@min_weight] <- method@min_weight

  list(
    weights = w,
    coefficients = NULL,
    duals = duals,
    converged = isTRUE(result$converged),
    interrupted = isTRUE(result$interrupted),
    iterations = as.integer(result$iterations),
    objective = as.numeric(result$objective),
    solver_status = result$solver_status,
    # The backend's terminal status name, carried so the orchestrator maps an
    # infeasible constraint set and a hard solver failure to their designed
    # conditions rather than the generic convergence warning.
    status = result$status,
    estimating_equations = NULL,
    # The energy objective drives approximate balance when no moment constraints
    # are enforced, so the balance warning must not fire against a tolerance the
    # fit does not target.
    approximate = approximate,
    groups = groups
  )
}
