# Stable balancing weights: minimum-dispersion weights under approximate
# covariate balance. Among all weightings that hold each reweighted group's
# covariate means inside a tolerance band, the method selects the one of least
# weight dispersion, so the balance tolerance is the central tuning parameter
# rather than an optional relaxation. The problem is a strictly convex quadratic
# program, so the method spec carries only tuning parameters and fit_method()
# assembles the standardized balance columns, calls the Rust solver, renormalizes
# each group to its estimand target total, and reports the solver's dual
# variables for diagnostics. The quadratic-program family has no estimating
# equations.

#' Stable balancing weights
#'
#' `bw_sbw()` specifies stable balancing weights for [balance()]. Among all
#' weightings that hold each reweighted exposure group's covariate means inside a
#' tolerance band, stable balancing weights select the one of least dispersion,
#' following Zubizarreta. The default `"l2"` norm minimizes the sum of squared
#' weights, so for a fixed per-group total it minimizes the weight variance.
#' Stable balancing weights support binary, categorical, and continuous
#' exposures.
#'
#' @details
#' The balance tolerance is the method's central tuning parameter. It is set
#' through `constraints = balance_terms(tolerance = ...)` in [balance()], not on
#' the method spec, and it must be positive: with an exact (zero) tolerance the
#' problem reduces to exact moment balance, which abandons the minimum-variance
#' rationale and is prone to infeasibility, so a fit without a positive tolerance
#' is refused. Given a positive tolerance the weights hold each group's weighted
#' covariate means within the band while minimizing the weight dispersion. For a
#' discrete exposure no feasible reweighting has smaller dispersion than the fit
#' returns.
#'
#' A single scalar tolerance applies to every covariate; a named vector sets a
#' tolerance per covariate. A covariate left at zero in a named vector demands
#' exact balance on that covariate while the others are relaxed, the
#' infeasibility-prone case, so a mixed specification is an explicit choice rather
#' than a convenience.
#'
#' For a discrete exposure the constraints default to first-moment balance; pass
#' [balance_terms()] to balance higher moments, interactions, or quantiles, each
#' inside its tolerance band. For a focal estimand the focal group keeps unit
#' weight and the other groups are pulled to the focal group's covariate means.
#'
#' For a continuous exposure the weighted exposure-covariate correlations are held
#' within the tolerance. The quadratic program bounds a linearized correlation
#' whose scales are fixed at the sample, so the fit tightens that internal bound
#' over a few passes until the reported weighted correlation sits inside the
#' requested band. The returned weights therefore minimize dispersion over the
#' tightened internal band rather than over every weighting that meets the
#' reported band, so the strict minimum-dispersion guarantee is stated for
#' discrete exposures only.
#'
#' Stable balancing weights belong to the quadratic-program family, which has no
#' estimating equations, so a fit produces no estimating-equations container.
#'
#' The `norm` argument selects how the weight dispersion is measured, always
#' against the uniform baseline of one within each reweighted group. `"l2"`
#' minimizes the sum of squared weights, so for a fixed per-group total it
#' minimizes the weight variance. `"l1"` minimizes the sum of absolute deviations
#' from one, which tends to leave many weights untouched and concentrate the
#' reweighting on a few units. `"linf"` minimizes the single largest absolute
#' deviation from one, which spreads the reweighting as evenly as the balance
#' constraints allow. The `"l1"` and `"linf"` problems are linear programs solved
#' through the same quadratic-program backends as `"l2"`; their solutions can be
#' non-unique, so a fit reports the achieved dispersion rather than promising a
#' unique weighting.
#'
#' @param norm The weight-dispersion norm to minimize, one of `"l2"` (the sum of
#'   squared weights, minimum variance), `"l1"` (the sum of absolute deviations
#'   from one), or `"linf"` (the largest absolute deviation from one).
#' @param min_weight The smallest permitted weight.
#' @param convergence_tolerance The quadratic-program solver tolerance, or `NULL`
#'   for the core default.
#' @param max_iterations The maximum solver iterations, or `NULL` for the core
#'   default.
#' @param ... Reserved for future extensions; must be empty. Tuning parameters
#'   must be passed by name.
#'
#' @return An `bw_sbw` specification, a [balance_method].
#'
#' @references
#' Zubizarreta, J. R. (2015). Stable weights that balance covariates for
#' estimation with incomplete outcome data. *Journal of the American Statistical
#' Association*, 110(511), 910-922.
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
#' fit <- balance(
#'   df,
#'   exposure,
#'   c(x1, x2),
#'   method = bw_sbw(),
#'   constraints = balance_terms(tolerance = 0.05)
#' )
#' fit
#'
#' @export
bw_sbw <- new_class(
  "bw_sbw",
  parent = quadratic_program_method,
  properties = list(
    norm = class_character
  ),
  constructor = function(
    ...,
    norm = c("l2", "l1", "linf"),
    min_weight = 1e-8,
    convergence_tolerance = NULL,
    max_iterations = NULL
  ) {
    check_method_dots(...)
    norm <- rlang::arg_match(norm)
    min_weight <- vctrs::vec_cast(min_weight, double(), x_arg = "min_weight")
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
      norm = norm,
      # Stable balancing weights minimize the weight dispersion alone; there is no
      # separate penalty, so the parent's weight penalty is held at zero.
      weight_penalty = 0,
      min_weight = min_weight,
      convergence_tolerance = convergence_tolerance,
      max_iterations = max_iterations
    )
  },
  # The minimum-weight floor is validated by the quadratic-program parent, which
  # declares it, together with the weight penalty this method holds at zero.
  validator = function(self) {
    if (!self@norm %in% c("l2", "l1", "linf")) {
      return("@norm must be one of \"l2\", \"l1\", or \"linf\"")
    }
  }
)

method(method_label, bw_sbw) <- function(method) {
  "Stable balancing weights"
}

method(supported_exposure_types, bw_sbw) <- function(method) {
  c("binary", "categorical", "continuous")
}

method(supported_estimands, bw_sbw) <- function(method, exposure_type) {
  switch(
    exposure_type,
    # The overlap estimand is legal only for the covariate balancing propensity
    # score, so stable balancing weights offer the same set as energy balancing.
    binary = c("ate", "att", "atu"),
    categorical = c("ate", "att"),
    continuous = "ate"
  )
}

# The quadratic-program family has no estimating equations for any exposure type
# or constraint set, so the answer is always FALSE. The context arguments are
# accepted so the generic call shape matches the estimating-equation family.
method(supports_estimating_equations, bw_sbw) <- function(
  method,
  ...,
  exposure_type = NULL,
  constraints = NULL
) {
  rlang::check_dots_empty()
  FALSE
}

# Stable balancing weights need explicit balance constraints, so their default is
# first-moment balance rather than the empty set the objective-driven
# quadratic-program methods use.
method(default_constraints, bw_sbw) <- function(method) {
  balance_terms(moments = 1L)
}

# The minimum-dispersion objective is minimized by the uniform weighting on its
# own, so the balance constraints are what makes the problem a balancing problem.
# Stable balancing weights therefore refuse an empty constraint set rather than
# returning uniform weights, unlike the other quadratic-program methods whose
# objective measures balance itself.
method(requires_constraints, bw_sbw) <- function(method) {
  TRUE
}

# Assemble the Rust option list, dropping the tuning parameters left at the core
# default so the quadratic-program solver applies its own. The worker-thread
# count is resolved on the R side and passed down on every call.
sbw_options <- function(method) {
  options <- list(
    threads = resolve_threads(),
    backend = resolve_qp_backend()
  )
  if (!is.null(method@convergence_tolerance)) {
    options$convergence_tolerance <- method@convergence_tolerance
  }
  if (!is.null(method@max_iterations)) {
    options$max_iterations <- as.integer(method@max_iterations)
  }
  options
}

# Announce that the default backend certified the problem infeasible and the fit
# re-solved with the interior-point backend. Carried as a classed message so it is
# snapshot-tested and a caller can mute or catch it; gated by the quiet option
# like the covariate-expansion alerts.
alert_backend_fallback <- function() {
  if (be_quiet()) {
    return(invisible())
  }
  cli::cli_inform(
    c(
      i = "The default solver certified this problem infeasible; re-solved with the {.val clarabel} backend."
    ),
    class = c("balancing_fallback_message", "balancing_message")
  )
}

# The core estimand string. The orchestrator stores the untreated target as
# "atu"; the binary core names it "atc".
sbw_core_estimand <- function(estimand) {
  switch(estimand, ate = "ate", att = "att", atu = "atc")
}

# Refuse a fit without a positive balance tolerance. The tolerance is the
# method's central tuning parameter, so an absent or zero tolerance is a
# constraints error naming the knob to set rather than a solver failure.
enforce_positive_tolerance <- function(prepared, call = rlang::caller_env()) {
  if (!any(prepared$tolerances > 0)) {
    abort(
      c(
        "Stable balancing weights require a positive balance tolerance.",
        x = "No constraint carries a tolerance above zero.",
        i = "Set {.arg tolerance} in {.fn balance_terms} to a positive value, the central tuning parameter for {.fn bw_sbw}."
      ),
      error_class = "balancing_constraints_error",
      call = call
    )
  }
}

# OSQP's C core writes its setup-validation errors, such as a non-positive
# variable count, straight to stderr, past the verbosity setting the backend
# turns off. Nothing on the R side can capture that, so the guarantee has to be
# that such a spec is never assembled. It is not: `balance()` refuses an empty
# data frame, a single-level exposure, sampling weights that are zero throughout,
# and any exposure level with no base-measure mass, all before a method is
# dispatched at all. Every quadratic program therefore reaches the backend with
# at least one unit and at least one constraint row, which is the whole of what
# the setup validation asks for. The same holds for energy and characteristic
# function distance balancing, which share those refusals and the same assembler.
method(fit_method, bw_sbw) <- function(method, prepared) {
  enforce_positive_tolerance(prepared)

  if (identical(prepared$exposure_type, "continuous")) {
    fit_sbw_continuous(method, prepared)
  } else {
    fit_sbw_discrete(method, prepared)
  }
}

fit_sbw_discrete <- function(method, prepared) {
  z <- prepared$matrix
  n <- prepared$n
  s <- prepared$sampling_weights
  levels <- prepared$exposure_levels
  groups <- prepared$groups
  estimand <- prepared$estimand
  focal <- prepared$focal_level
  key <- prepared$exposure_key

  # The moment rows hold each reweighted group's weighted mean of the
  # standardized covariate columns within a tolerance band of a target. For the
  # average treatment effect the target is the pooled sample mean; for a focal
  # estimand it is the focal group's mean. The band is applied at its full width
  # against the target for every estimand, following Zubizarreta's formulation, so
  # the achieved arm-to-target standardized mean difference matches the tolerance
  # the method is asked for. For a focal estimand this is the reference's own
  # single held-fixed band; for the average treatment effect it is a strict
  # superset of the reference's feasible set, which additionally pins the pair
  # average and so binds each arm at half width. In both cases the reference
  # weights are feasible in this band, so their dispersion bounds ours from above.
  # Numeric columns cross standardized to unit scale, so their band is the
  # tolerance directly; indicator columns cross raw, so their band scales by the
  # column standard deviation.
  if (identical(estimand, "ate")) {
    targets <- target_means(z, seq_len(n), s)
  } else {
    targets <- target_means(z, groups[[focal]], s)
  }
  tols <- solver_box(z, prepared$tolerances, s)
  options <- sbw_options(method)

  if (identical(prepared$exposure_type, "binary")) {
    core_estimand <- sbw_core_estimand(estimand)
    if (identical(estimand, "ate")) {
      # For the average treatment effect every group is pulled to the shared
      # target, so the treated coding is only the second level's indicator.
      treat <- as.integer(key == levels[[2]])
      nvar <- n
      n_group_rows <- 2L
    } else {
      # A focal estimand holds the focal group at unit weight. The core reads the
      # focal group from the treatment coding: the treated target codes the focal
      # group one, the untreated target codes it zero.
      treat <- if (identical(estimand, "att")) {
        as.integer(key == focal)
      } else {
        as.integer(key != focal)
      }
      nvar <- n - length(groups[[focal]])
      n_group_rows <- 1L
    }
    result <- solve_sbw(
      treat,
      s,
      core_estimand,
      method@norm,
      z,
      targets,
      tols,
      method@min_weight,
      options
    )
  } else {
    treat_idx <- match(key, levels) - 1L
    if (identical(estimand, "ate")) {
      focal_idx <- 0L
      core_estimand <- "ate"
      nvar <- n
      n_group_rows <- length(levels)
    } else {
      focal_idx <- match(focal, levels) - 1L
      core_estimand <- "att"
      nvar <- n - length(groups[[focal]])
      n_group_rows <- length(levels) - 1L
    }
    result <- solve_sbw_multi(
      as.integer(treat_idx),
      as.integer(focal_idx),
      s,
      core_estimand,
      method@norm,
      z,
      targets,
      tols,
      method@min_weight,
      options
    )
  }

  # Each reweighted group carries one sum row and one moment row per covariate, so
  # the sum and balance rows number `n_group_rows * (1 + ncol(z))`; the
  # absolute-deviation norms add auxiliary rows beyond these that the dual report
  # excludes.
  n_keep <- n_group_rows * (1L + ncol(z))
  duals <- sbw_duals_frame(result$duals, nvar, n_group_rows, n_keep)
  assemble_sbw(result, method, prepared, duals = duals)
}

# The largest number of correlation-refinement passes and the fraction of the
# room to the target the effective tolerance is tightened to on each pass, held
# just under one so a converged fit sits inside the band rather than on its edge.
sbw_cont_max_passes <- 8L
sbw_cont_safety <- 0.98

# Absolute weighted exposure-covariate Pearson correlations under weights `w`, the
# statistic the continuous fit is judged on and the quantity the balance table
# reports, so the refinement loop measures the same thing the specs assert. A
# column with no weighted spread has no correlation to report: it is met by every
# weighting, so it reads as zero rather than carrying an undefined value into the
# comparison that decides which tolerances still bind.
sbw_weighted_correlations <- function(exposure, z, w) {
  vapply(
    seq_len(ncol(z)),
    function(j) {
      correlation <- stats::cov.wt(
        cbind(exposure, z[, j]),
        wt = w,
        cor = TRUE
      )$cor[1, 2]
      if (is.finite(correlation)) abs(correlation) else 0
    },
    numeric(1)
  )
}

fit_sbw_continuous <- function(method, prepared) {
  z <- prepared$matrix
  s <- prepared$sampling_weights
  exposure <- as.numeric(prepared$exposure_vec)
  options <- sbw_options(method)

  # The quadratic program bounds a linearized weighted correlation whose exposure
  # and covariate scales are fixed at the sampling-weight sample. Reweighting to
  # meet the bound shrinks both weighted standard deviations, so the true weighted
  # Pearson correlation the fit is judged on runs above the linearized bound by the
  # product of the two shrinkage ratios. A one-shot solve at the requested
  # tolerance therefore overshoots the correlation band. The effective tolerance
  # passed to the solver is tightened over a few passes until the achieved true
  # correlation sits inside the requested band; each pass measures the same
  # weighted Pearson correlation the specs check, on the reported weights (the
  # balancing weights composed with the sampling weights), and rescales each
  # column's effective tolerance toward its target, never above it, so the loop
  # tightens monotonically and stops once every column is inside its band.
  target <- prepared$tolerances
  effective <- target
  result <- NULL
  last_converged <- NULL
  for (pass in seq_len(sbw_cont_max_passes)) {
    result <- solve_sbw_cont(
      exposure,
      z,
      s,
      method@norm,
      effective,
      method@min_weight,
      options
    )
    if (!isTRUE(result$converged)) {
      # A tightened pass that certifies infeasibility means the requested
      # correlation band is unreachable, which surfaces honestly as the infeasible
      # condition. A pass that merely ran out of iterations falls back to the last
      # converged iterate, which the balance warning then judges, rather than
      # returning floor-level weights behind a convergence warning.
      if (
        !identical(result$status, "primal_infeasible") &&
          !is.null(last_converged)
      ) {
        result <- last_converged
      }
      break
    }
    last_converged <- result
    # cov.wt normalizes internally, so the composed sampling weights, not their
    # renormalized copy, carry the reweighting the reported statistic reflects.
    composed <- as.numeric(result$weights) * s
    achieved <- sbw_weighted_correlations(exposure, z, composed)
    binding <- target > 0 & achieved > target + 1e-8
    if (!any(binding)) {
      break
    }
    ratio <- ifelse(achieved > 0, target / achieved, 1)
    effective[binding] <- pmin(
      target[binding],
      effective[binding] * ratio[binding] * sbw_cont_safety
    )
  }

  # The continuous solve carries one total-sum row followed by the correlation
  # rows; the box rows bound each of the n units, and the absolute-deviation norms
  # add auxiliary rows the dual report excludes.
  n_keep <- 1L + ncol(z)
  duals <- sbw_duals_frame(
    result$duals,
    prepared$n,
    1L,
    n_keep,
    group_kind = "total",
    moment_kind = "correlation"
  )
  assemble_sbw(result, method, prepared, duals = duals)
}

# The solver's dual variables for the structural constraint rows, dropping the
# per-unit box duals. The quadratic program stacks the box rows first, then the
# sum rows, then the balance rows, so the structural duals follow the `nvar` box
# rows. A discrete fit labels the sum rows "group" and the balance rows "moment"; a
# continuous fit has a single total-sum row and correlation balance rows, so its
# labels differ. The absolute-deviation norms append auxiliary and conversion rows
# after the balance rows, so `n_keep` caps the reported duals at the sum and
# balance rows and excludes those auxiliary rows; for the squared norm it equals
# the full structural count, so nothing is dropped. Returns `NULL` when no
# structural row resolved.
sbw_duals_frame <- function(
  duals,
  nvar,
  n_group_rows,
  n_keep,
  group_kind = "group",
  moment_kind = "moment"
) {
  duals <- as.numeric(duals)
  m <- length(duals)
  if (m <= nvar) {
    return(NULL)
  }
  last <- min(nvar + n_keep, m)
  structural <- duals[(nvar + 1):last]
  n_structural <- length(structural)
  n_group_rows <- min(n_group_rows, n_structural)
  kind <- c(
    rep(group_kind, n_group_rows),
    rep(moment_kind, n_structural - n_group_rows)
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
# floor is re-applied on that reported scale so the documented minimum holds where
# the specs assert it. For the average treatment effect each group is scaled to
# its own sampling-weighted total; for a focal estimand every group is scaled to
# the focal total, which leaves the focal group at unit weight. A continuous fit
# is a single group scaled to the sampling-weight total. Balance is enforced
# through the tolerance box rather than approximated, so the fit is not flagged
# approximate and the balance warning judges it against its tolerance.
assemble_sbw <- function(result, method, prepared, duals = NULL) {
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

  # The backend that produced the weights, and whether the automatic fallback
  # engaged. A fallback is recorded as "clarabel_fallback" so the behavior
  # difference is never silent, and it is announced with an informational alert.
  fell_back <- isTRUE(result$fell_back)
  backend <- result$solver_status
  solver_status <- if (fell_back) "clarabel_fallback" else backend
  # Announce a fallback only when it produced a usable solution; a fallback that
  # also certifies infeasibility surfaces as the infeasible error, not a rescue.
  if (fell_back && isTRUE(result$converged)) {
    alert_backend_fallback()
  }

  list(
    weights = w,
    coefficients = NULL,
    duals = duals,
    converged = isTRUE(result$converged),
    interrupted = isTRUE(result$interrupted),
    iterations = as.integer(result$iterations),
    objective = as.numeric(result$objective),
    solver_status = solver_status,
    # The backend's terminal status name, carried so the orchestrator maps an
    # infeasible constraint set and a hard solver failure to their designed
    # conditions rather than the generic convergence warning.
    status = result$status,
    estimating_equations = NULL,
    approximate = FALSE,
    groups = groups
  )
}
