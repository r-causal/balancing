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
#' marginal moments are held equal to the sample under the base measure, and
#' `dimension_adjustment` reweights the covariate energy distance by the
#' covariate dimensionality.
#'
#' The two knobs a continuous fit carries are separate, and the WeightIt package
#' names them the same way: `moments` in [balance_terms()] is WeightIt's
#' `moments`, adding a constraint that holds the weighted correlation of the
#' exposure with each covariate power within its tolerance, which defaults to
#' zero, and `distribution_moments` here is
#' WeightIt's `d.moments`, pinning the marginal moments of the exposure and of
#' the covariates. Neither sets the other here. WeightIt does couple them in one
#' direction: `weightit()` raises its own `d.moments` to its `moments`, so a fit
#' matching a WeightIt call with `moments = k` for `k` above one sets
#' `distribution_moments = k` here as well as `moments = k` in
#' [balance_terms()]. The correlation rows are held within the tolerance
#' [balance_terms()] carries, and reaching that band takes more than one solve.
#' The quadratic program bounds a linearized correlation whose exposure and
#' covariate scales are fixed at the sample, and the spread of energy weights
#' shrinks both weighted standard deviations, so a single solve at the requested
#' band overshoots it: a band of `0.05` lands between 0.070 and 0.086 at 200 to
#' 1000 observations. The fit therefore tightens the bound it hands the program
#' and re-solves, up to eight passes, until the reported correlation sits inside
#' the band. A band of `0.05` took two passes at 350 and at 1000 observations, so
#' it costs about two solves against the one the same fit at exact balance takes,
#' exact balance having nothing to tighten. A band the passes cannot reach is
#' reported at its last iterate, and the balance warning judges it as it judges
#' any other fit.
#'
#' Without those rows the continuous objective targets distributional
#' independence between the exposure and the covariates rather than zero
#' correlations, and it does not drive the correlations to zero. A residual
#' weighted correlation of roughly 0.1 to 0.3 is ordinary at a few hundred to a
#' few thousand observations, and WeightIt's continuous energy method leaves the
#' same residual. What holds it up is `weight_penalty`, which trades that
#' residual against effective sample size: at its default of `1e-4` the penalty
#' term is about three quarters of the objective at 1000 observations, leaving a
#' largest correlation near 0.22 to 0.25 at an effective sample size near 71
#' percent, while a penalty of zero brings the correlation down to 0.05 to 0.07
#' and the effective sample size down to about 20 percent. Ask for
#' `balance_terms(moments = 1)` to remove the correlation outright, at a cost in
#' effective sample size of its own.
#'
#' Energy balancing belongs to the quadratic-program family, which has no
#' estimating equations, so a fit produces no estimating-equations container and
#' the tolerance in [balance_terms()] relaxes the constraints a fit added rather
#' than selecting an inexact solver. A tolerance supplied without those
#' constraints has nothing to relax, so it is warned and ignored, and the balance
#' table reports the tolerance the fit enforced, which is zero.
#'
#' @param distance The covariate distance definition the energy objective is
#'   built on, one of `"scaled_euclidean"` (each covariate centered at its
#'   weighted mean and divided by its weighted standard deviation),
#'   `"mahalanobis"`, or `"euclidean"`.
#' @param improved Whether to add the between-group energy distance of the
#'   improved variant for the average treatment effect with a discrete exposure.
#' @param weight_penalty The L2 penalty on the weights, which stabilizes the
#'   quadratic program. For a continuous exposure it is also what sets the
#'   residual exposure-covariate correlation, as the details section explains.
#' @param min_weight The smallest permitted weight. The reported weights average
#'   one within each exposure group, so a floor approaching one leaves almost no
#'   room above it: the weight spread shrinks in proportion to the headroom
#'   `1 - min_weight`, and the fit degenerates smoothly into uniform weights and
#'   reports the balance uniform weights achieve. Nothing warns at that boundary,
#'   because the problem stays feasible and the solution is a real one.
#'   [bw_sbw()], whose tolerances are hard constraints rather than an objective,
#'   refuses the same floor as infeasible instead.
#' @param distribution_moments For a continuous exposure, the number of exposure
#'   and covariate marginal moments held equal to the sample under the base
#'   measure, or `NULL` for the first moments. This is WeightIt's `d.moments`,
#'   and it is the only route to those rows: the `moments` in [balance_terms()]
#'   asks for exposure-covariate correlation constraints instead and leaves the
#'   marginals here. Energy balancing carries no base weights, so the base
#'   measure is the sampling weights, and without them the marginals are held
#'   equal to the unweighted sample.
#' @param dimension_adjustment For a continuous exposure, whether to weight the
#'   covariate energy distance by the covariate dimensionality adjustment.
#' @param convergence_tolerance The quadratic-program solver tolerance, which
#'   the solver applies as both its absolute and its relative tolerance, or
#'   `NULL` for the family default of `1e-8`. Energy balancing defaults to
#'   `1e-6` rather than to that family default because its quadratic form is
#'   indefinite: on a small sample the alternating-direction residual floors
#'   above `1e-8`, and a run that keeps going past that floor walks away from
#'   the optimum instead of stalling at it. The energy objective always solves
#'   through the alternating-direction backend, whatever `balancing.qp_backend`
#'   names, so there is no backend to choose here: a tolerance below what the
#'   problem can reach spends the full iteration cap, then warns and reports the
#'   iterate of a re-solve at a tolerance the problem does reach, provided that
#'   re-solve converges within the same `max_iterations`. When it does not, the
#'   fit reports the iterate of the original solve.
#' @param max_iterations The maximum solver iterations, or `NULL` for the
#'   resolved default of 200000. The re-solve above is given the same cap, and
#'   when it converges the reported `@iterations` sums the two solves, so an
#'   energy fit that could not reach its tolerance can report more iterations
#'   than this. The refinement passes of a continuous fit with a positive
#'   tolerance are summed the same way, each pass being a solve of its own.
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
    ...,
    distance = c("scaled_euclidean", "mahalanobis", "euclidean"),
    improved = TRUE,
    weight_penalty = 1e-4,
    min_weight = 1e-8,
    distribution_moments = NULL,
    dimension_adjustment = TRUE,
    convergence_tolerance = 1e-6,
    max_iterations = NULL
  ) {
    check_method_dots(...)
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
  # The weight penalty and the minimum-weight floor are validated by the
  # quadratic-program parent, which declares them. The distribution moments reach a
  # comparison in the continuous fit path, so a missing value or a vector is
  # refused here rather than left to stop that fit with a base error.
  validator = function(self) {
    if (length(self@improved) != 1 || is.na(self@improved)) {
      return("@improved must be a single logical value")
    }
    if (
      length(self@dimension_adjustment) != 1 ||
        is.na(self@dimension_adjustment)
    ) {
      return("@dimension_adjustment must be a single logical value")
    }
    if (!is.null(self@distribution_moments)) {
      return(validate_distribution_moments(self@distribution_moments))
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

# The weight penalty is a tuning argument here, so a solver breakdown may advise
# raising it.
method(tunes_weight_penalty, bw_energy) <- function(method) {
  TRUE
}

# The energy objective is built from the negative pairwise distance, so its
# quadratic form is indefinite.
method(has_indefinite_objective, bw_energy) <- function(method) {
  TRUE
}

# Assemble the Rust option list, dropping the tuning parameters left at the core
# default so the quadratic-program solver applies its own. The worker-thread
# count and the quadratic-program backend are resolved on the R side and passed
# down on every call.
energy_options <- function(method, backend) {
  options <- list(threads = resolve_threads(), backend = backend)
  if (!is.null(method@convergence_tolerance)) {
    options$convergence_tolerance <- method@convergence_tolerance
  }
  if (!is.null(method@max_iterations)) {
    options$max_iterations <- as.integer(method@max_iterations)
  }
  options
}

# Solve, and on a run that spent its iteration cap solve once more at a
# tolerance the objective can reach. The energy quadratic form is indefinite, so
# the alternating-direction iteration is a contraction only until its residuals
# reach the floor of the problem; a tolerance below that floor keeps the run
# going, and the iterate it carries at the cap has left the optimum rather than
# stopped short of it. Renormalizing that iterate per group gives it the shape of
# a real answer, so it must not be what the fit reports. The retry costs one
# extra solve, and only on a fit that already failed. The fit still reports
# itself as unconverged, because the tolerance the caller asked for was not met,
# and the convergence warning that follows names the tolerance to ask for
# instead. A run that ends on any other terminal status is diagnosed by
# check_solver_status() and is not retried, and neither is a fit whose tolerance
# is already at or above the reachable one, where the retry would repeat the
# solve that just failed.
solve_energy_with_fallback <- function(method, options, solve) {
  result <- solve(options)
  reached_cap <- !isTRUE(result$converged) &&
    identical(result$status, "max_iter")
  if (!reached_cap || resolved_qp_tolerance(method) >= qp_reachable_tolerance) {
    return(result)
  }
  options$convergence_tolerance <- qp_reachable_tolerance
  retry <- solve(options)
  if (!isTRUE(retry$converged)) {
    return(result)
  }
  retry$converged <- FALSE
  retry$status <- result$status
  retry$iterations <- result$iterations + retry$iterations
  retry
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
#
# The column is read through the same accessor the constraint builder uses, so a
# duration reaches the distance as the number it stores rather than falling to
# the categorical branch and becoming one indicator per distinct duration.
distance_covariates <- function(data, covariates) {
  columns <- lapply(covariates, function(covariate) {
    values <- covariate_values(data, covariate)
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

# The tolerance in a balance_terms() specification relaxes added constraints;
# with none present it has nothing to act on, so warn and proceed with the pure
# objective. Both `moments` and `interactions` add constraints for either
# exposure type: moment rows for a discrete exposure, exposure-covariate
# correlation rows on the continuous energy path. Only `quantiles` is confined to
# a discrete exposure, so it is the only one the advice qualifies.
warn_ignored_tolerance <- function(call = rlang::caller_env()) {
  warn(
    c(
      "{.arg tolerance} relaxes added constraints, but this fit has none to relax.",
      i = "Drop {.arg tolerance} from {.fn balance_terms}, or add constraints with {.arg moments} or {.arg interactions}, or with {.arg quantiles} for a discrete exposure."
    ),
    warning_class = "balancing_ignored_argument_warning",
    call = call
  )
}

# Resolve the quadratic-program backend for an energy fit, and announce a pinned
# interior-point backend the objective cannot use. The energy Gram matrix is only
# conditionally positive semidefinite, so the quadratic term the fit assembles is
# indefinite and the interior-point backend refuses an indefinite form up front.
# Energy balancing therefore routes to the ADMM backend whatever was requested,
# which is the correct route rather than an optional one, so the fit proceeds;
# what it owes the caller is to say the request was dropped, as a method
# constructor does for a tuning argument its exposure type cannot use. The
# automatic policy asks for no particular backend and has nothing to report, and
# the recorded solver status names what ran either way. Routing the option
# through the shared resolver is also what makes an unknown value an error here
# rather than a silent default.
resolve_energy_backend <- function(call = rlang::caller_env()) {
  backend <- resolve_qp_backend()
  if (!identical(backend, "clarabel")) {
    return(backend)
  }
  warn(
    c(
      "The {.code balancing.qp_backend} option is {.val clarabel}, which energy balancing cannot use, and is ignored.",
      x = "The energy objective's quadratic form is indefinite, and the interior-point backend solves only positive-semidefinite forms.",
      i = "The fit used {.val osqp} instead, which {.code @solver_status} records."
    ),
    warning_class = "balancing_ignored_argument_warning",
    call = call
  )
  "osqp"
}

method(fit_method, bw_energy) <- function(method, prepared) {
  backend <- resolve_energy_backend()
  enforce <- requests_moments(prepared$constraints)

  if (!enforce && has_positive_tolerance(prepared$constraints)) {
    warn_ignored_tolerance()
  }

  if (identical(prepared$exposure_type, "continuous")) {
    return(fit_energy_continuous(method, prepared, enforce, backend))
  }

  fit_energy_discrete(method, prepared, enforce, backend)
}

fit_energy_discrete <- function(method, prepared, enforce, backend) {
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
    tols <- solver_box(z, prepared$tolerances, s)
    moment_covs <- z
  } else {
    targets <- numeric(0)
    tols <- numeric(0)
    moment_covs <- matrix(numeric(0), nrow = n, ncol = 0)
  }

  options <- energy_options(method, backend)

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
    result <- solve_energy_with_fallback(method, options, function(opts) {
      solve_energy(
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
        opts
      )
    })
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
    result <- solve_energy_with_fallback(method, options, function(opts) {
      solve_energy_multi(
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
        opts
      )
    })
  }

  duals <- energy_duals_frame(result$duals, nvar, n_group_levels)
  assemble_energy(
    result,
    method,
    prepared,
    duals,
    approximate = !enforce,
    enforced_tolerance = if (enforce) NULL else 0
  )
}

# Shift each distribution-moment column to its mean under the base measure. The
# quadratic program pins every distribution row at a weighted mean of zero, so
# where a column is centered decides which population its row targets: centering
# on the base measure pins the row at that measure's mean of the incoming
# column. Numeric covariate columns arrive centered on the sampling weights and
# do not move. The exposure powers and the higher covariate marginals arrive
# centered on the unweighted sample, so leaving them there would hold the
# exposure to one population while the covariates are held to another. An
# indicator column arrives raw, so a row pinned at its zero would drive the
# indicated stratum to no weight at all rather than to its proportion.
center_on_measure <- function(columns, measure) {
  centers <- as.numeric(crossprod(columns, measure / sum(measure)))
  sweep(columns, 2, centers, "-")
}

# The first-moment marginal columns of a continuous fit: one indicator per level
# of a factor covariate, the covariate itself where it is already an indicator,
# and the standardized first power of a numeric covariate. Returns the columns
# alongside the highest marginal power each covariate reached, which
# higher_covariate_marginals() continues from.
#
# They are built from the covariates rather than read off the constraint recipe
# because a covariate's marginal distribution is not what the constraint set
# selects. Reading them off the recipe left `moments` a second route to the
# marginals: `balance_terms(moments = c(x1 = 0))` drops x1's constraint record,
# and with it x1's marginal row, so a fit asked to leave x1 out of the
# correlation rows stopped holding x1's own distribution as well.
#
# The columns cross the boundary on the same scale the constraint matrix uses,
# so they are built through the same records and the same rebuild, and the
# constant and aliased columns are dropped exactly as the constraint build drops
# them. A factor's indicators sum to the constant every method carries, so one
# of them is redundant against the fit's own total-sum row. The drops are silent
# here: the constraint build has already reported whatever it dropped, and these
# rows are the fit's own bookkeeping rather than a set the caller named.
marginal_distribution_columns <- function(data, covariates, sampling_weights) {
  center_fn <- if (is.null(sampling_weights)) {
    mean
  } else {
    function(x) weighted_center(x, sampling_weights)
  }
  scale_fn <- if (is.null(sampling_weights)) {
    stats::sd
  } else {
    function(x) weighted_scale(x, sampling_weights)
  }

  records <- list()
  for (cov in covariates) {
    v <- covariate_values(data, cov)
    if (is.factor(v) || is.character(v)) {
      levels <- if (is.factor(v)) {
        levels(v)
      } else {
        sort(unique(as.character(v)))
      }
      for (level in levels) {
        records[[length(records) + 1L]] <- new_recipe_record(
          term = paste0(cov, "_", level),
          kind = "moment",
          type = "indicator",
          source = cov,
          level = level
        )
      }
    } else if (is.logical(v) || is_binary_numeric(v)) {
      records[[length(records) + 1L]] <- new_recipe_record(
        term = cov,
        kind = "moment",
        type = "indicator",
        source = cov,
        level = NA_character_
      )
    } else {
      base_center <- center_fn(v)
      raw <- v - base_center
      scale <- scale_fn(raw)
      if (scale == 0) {
        scale <- 1
      }
      records[[length(records) + 1L]] <- new_recipe_record(
        term = cov,
        kind = "moment",
        type = "numeric",
        source = cov,
        power = 1L,
        base_center = base_center,
        center = center_fn(raw),
        scale = scale
      )
    }
  }

  columns <- rebuild_constraint_matrix(records, data)
  for (drop in list(constant_columns, aliased_columns)) {
    dropped <- drop(columns)
    if (length(dropped) > 0) {
      keep <- setdiff(seq_along(records), dropped)
      records <- records[keep]
      columns <- columns[, keep, drop = FALSE]
    }
  }

  list(
    columns = columns,
    moments = covariate_constraint_moments(records, covariates)
  )
}

fit_energy_continuous <- function(method, prepared, enforce, backend) {
  n <- prepared$n
  s <- prepared$sampling_weights
  covs <- distance_covariates(prepared$data, prepared$covariates)
  z <- prepared$matrix
  exposure <- as.numeric(prepared$exposure_vec)

  # The reference distribution the marginals are held to is the base measure, as
  # it is for entropy balancing. The quadratic-program family carries no base
  # weights of its own, so that measure is the sampling weights, and without them
  # it is uniform and the marginals reduce to the unweighted sample.
  measure <- s

  # The distribution-moment constraints hold the weighted exposure and covariate
  # marginals equal to the sample under the base measure. Every one of those rows
  # takes the same measure, so a fit under informative sampling holds the
  # exposure and the covariates to one population rather than two.
  # `distribution_moments` is the only argument that sets how many of them there
  # are: the constraint set a caller passes to balance_terms() asks for
  # exposure-covariate correlation rows here, and the marginal rows are built
  # from the covariates and `distribution_moments` alone.
  moments <- method@distribution_moments %||% 1L
  marginals <- marginal_distribution_columns(
    prepared$data,
    prepared$covariates,
    prepared$constraint_sampling_weights
  )

  d_treat <- center_on_measure(moment_columns(exposure, moments), measure)
  extra_covariate_marginals <- higher_covariate_marginals(
    prepared$data,
    marginals$moments,
    moments
  )
  d_covs <- center_on_measure(
    do.call(
      cbind,
      c(list(marginals$columns), extra_covariate_marginals)
    ),
    measure
  )

  # The correlation rows hold the weighted correlation of the exposure with each
  # constraint column inside that column's tolerance, which is what `moments` and
  # `interactions` in balance_terms() ask for on a continuous exposure and what
  # WeightIt's `moments` argument means for a continuous treatment. The core
  # standardizes the exposure on the base measure itself and the first
  # distribution row pins its weighted mean there, so a row driven to zero is a
  # weighted covariance of zero rather than one offset by the gap between the two
  # exposure means, and a column whose own weighted mean is not pinned is covered
  # as well. Without a requested constraint set the weighted distance covariance
  # the objective minimizes is what drives the association toward zero, and no
  # correlation row is added.
  target <- if (enforce) {
    prepared$tolerances
  } else {
    numeric(0)
  }
  bal_covs <- if (enforce) {
    z
  } else {
    matrix(numeric(0), nrow = n, ncol = 0)
  }

  options <- energy_options(method, backend)

  # The row the quadratic program bounds is a linearized correlation whose
  # exposure and covariate scales are fixed at the sampling-weight sample.
  # Reweighting to meet the bound shrinks both weighted standard deviations, so
  # the true weighted Pearson correlation the fit is judged on runs above the
  # bound by the product of the two shrinkage ratios: a single solve at a
  # tolerance of 0.05 lands between 0.070 and 0.086 at 200 to 1000 observations.
  # The bound handed to the program is therefore tightened over a few passes
  # until the reported correlation sits inside the requested band, the same
  # refinement fit_sbw_continuous() runs against the same statistic. Each pass
  # rescales a binding column's bound toward its target, never above it, so the
  # loop tightens monotonically, and it stops once every column is inside the
  # band the balance table judges it against. A band the passes cannot reach is
  # kept at the last iterate that converged and reported: the table then judges
  # it out of balance and the fit warns through the ordinary balance warning.
  #
  # Exact balance and a fit with no correlation rows have nothing to tighten and
  # take a single pass. Each pass costs a whole solve, and the reported
  # iterations sum every one of them.
  effective <- target
  iterations <- 0L
  result <- NULL
  last_converged <- NULL
  for (pass in seq_len(correlation_refinement_passes)) {
    result <- solve_energy_with_fallback(method, options, function(opts) {
      solve_energy_cont(
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
        effective,
        opts
      )
    })
    iterations <- iterations + as.integer(result$iterations)
    if (!isTRUE(result$converged)) {
      # A tightened pass that certifies infeasibility means the requested
      # correlation band is unreachable, which surfaces honestly as the
      # infeasible condition. A pass that merely ran out of iterations falls
      # back to the last converged iterate, which the balance warning then
      # judges, rather than reporting the unsettled iterate a tightened bound
      # left behind.
      if (
        !identical(result$status, "primal_infeasible") &&
          !is.null(last_converged)
      ) {
        result <- last_converged
      }
      break
    }
    # A fit with nothing to tighten leaves before the measurement as well as
    # before the second solve, so the default fit pays for no correlation it
    # would not have computed.
    if (!any(target > 0)) {
      break
    }
    last_converged <- result
    # cov.wt normalizes internally, so the composed sampling weights, not their
    # renormalized copy, carry the reweighting the reported statistic reflects.
    composed <- as.numeric(result$weights) * s
    achieved <- weighted_exposure_correlations(exposure, z, composed)
    binding <- target > 0 & achieved > target + balance_margin(target)
    if (!any(binding)) {
      break
    }
    ratio <- ifelse(achieved > 0, target / achieved, 1)
    effective[binding] <- pmin(
      target[binding],
      effective[binding] * ratio[binding] * correlation_refinement_safety
    )
  }
  result$iterations <- iterations

  # The continuous solve carries one total-sum row followed by the distribution
  # rows; the box rows bound each of the n units.
  n_structural_leading <- 1L
  duals <- energy_duals_frame(result$duals, n, n_structural_leading)
  assemble_energy(
    result,
    method,
    prepared,
    duals,
    approximate = !enforce,
    enforced_tolerance = if (enforce) NULL else 0
  )
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
assemble_energy <- function(
  result,
  method,
  prepared,
  duals,
  approximate,
  enforced_tolerance = NULL
) {
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
    # The tolerance the fit held its rows at, where that is the method's own
    # value rather than the one the specification asked for, so the balance
    # table reports what the program enforced. A fit that added no constraint
    # rows is the case: a tolerance reached no row of the program, so what it
    # enforced is zero, and reporting the requested band instead would print a
    # box nothing was placed in. A fit that did add rows holds them inside the
    # requested band, so it leaves this `NULL` and the table reads the
    # per-column tolerances the specification named.
    enforced_tolerance = enforced_tolerance,
    groups = groups
  )
}
