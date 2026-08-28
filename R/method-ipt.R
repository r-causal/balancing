# Inverse probability tilting: a propensity model whose tilted score equations
# force the weighted covariate means to their estimand targets, so the achieved
# balance is exact on the requested moments. The method specification carries
# only tuning parameters; fit_method() prepends the intercept the propensity
# model needs, calls the Rust solver, applies the estimand's group-sum
# normalization, and passes the raw estimating equations through unchanged.

#' Inverse probability tilting
#'
#' `bw_ipt()` specifies inverse probability tilting for [balance()]. A propensity
#' model is fit not by maximum likelihood but by a tilted moment condition that
#' forces each treatment group's weighted covariate means to their estimand
#' targets, so balance on the requested moments is exact by construction.
#' Inverse probability tilting supports binary and categorical exposures.
#'
#' @details
#' For each treatment level the propensity `p_i = G(x_i' beta)` is estimated so
#' that the weighted covariate total matches the target population total. The
#' average treatment effect tilts every level to the whole sample, weighting a
#' unit by the inverse of its modeled propensity. A focal estimand tilts each
#' non-focal level to the focal level, weighting a unit by `(1 - p_i) / p_i`,
#' and leaves the focal units at weight one. With mean balance and the logit
#' link the tilt solves the same treated-target problem as entropy balancing, so
#' the two methods produce the same average-treatment-effect-on-the-treated
#' weights for a binary exposure.
#'
#' The weights solve smooth estimating equations regardless of the requested
#' tolerance, which [balance()] records for the M-estimation variance in
#' [`ipw()`][ipw.balancing].
#'
#' @param link The propensity link, one of `"logit"`, `"probit"`, or
#'   `"cloglog"`.
#' @param convergence_tolerance The solver convergence tolerance on the tilting
#'   moment. `1e-10` is both this argument's default and the value the solver
#'   resolves for `NULL`.
#' @param max_iterations The maximum solver iterations, or `NULL` for the
#'   resolved default of 1000.
#' @param ... Reserved for future extensions; must be empty. Tuning parameters
#'   must be passed by name.
#'
#' @return An `bw_ipt` specification, a [balance_method].
#'
#' @references
#' Graham, B. S., Pinto, C. C. de X., and Egel, D. (2012). Inverse probability
#' tilting for moment condition models with missing data. *The Review of
#' Economic Studies*, 79(3), 1053-1079.
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
#' fit <- balance(df, exposure, c(x1, x2), method = bw_ipt())
#' fit
#'
#' @export
bw_ipt <- new_class(
  "bw_ipt",
  parent = estimating_equation_method,
  properties = list(
    link = class_character
  ),
  constructor = function(
    ...,
    link = c("logit", "probit", "cloglog"),
    convergence_tolerance = 1e-10,
    max_iterations = NULL
  ) {
    check_method_dots(...)
    link <- rlang::arg_match(link)
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
      link = link,
      convergence_tolerance = convergence_tolerance,
      max_iterations = max_iterations
    )
  }
)

method(method_label, bw_ipt) <- function(method) {
  "Inverse probability tilting"
}

method(supported_exposure_types, bw_ipt) <- function(method) {
  c("binary", "categorical")
}

method(supported_estimands, bw_ipt) <- function(method, exposure_type) {
  switch(
    exposure_type,
    binary = c("ate", "att", "atu"),
    categorical = c("ate", "att"),
    continuous = character(0)
  )
}

# Inverse probability tilting always solves smooth estimating equations. Unlike
# entropy balancing, the tilt keeps them regardless of the requested tolerance,
# so the constraints do not change the answer. The tilt is derived for a discrete
# exposure and supports no continuous fit, so it has no continuous estimating
# equations to offer either.
method(supports_estimating_equations, bw_ipt) <- function(
  method,
  ...,
  exposure_type = NULL,
  constraints = NULL
) {
  rlang::check_dots_empty()
  !identical(exposure_type, "continuous")
}

# Assemble the Rust option list, dropping the tuning parameters left at the core
# default so the solver applies its own. The worker-thread count is resolved on
# the R side and passed down on every call.
ipt_options <- function(method) {
  options <- list(threads = resolve_threads())
  if (!is.null(method@convergence_tolerance)) {
    options$convergence_tolerance <- method@convergence_tolerance
  }
  if (!is.null(method@max_iterations)) {
    options$max_iterations <- as.integer(method@max_iterations)
  }
  options
}

method(fit_method, bw_ipt) <- function(method, prepared) {
  if (identical(prepared$exposure_type, "binary")) {
    fit_ipt_binary(method, prepared)
  } else {
    fit_ipt_categorical(method, prepared)
  }
}

fit_ipt_binary <- function(method, prepared) {
  n <- prepared$n
  s <- prepared$sampling_weights
  # The propensity model carries an intercept, whose tilting moment fixes each
  # group's weighted total and so drives the group-sum normalization.
  covs <- cbind(1, prepared$matrix)
  levels <- prepared$exposure_levels
  key <- prepared$exposure_key
  estimand <- prepared$estimand
  focal <- prepared$focal_level

  if (identical(estimand, "ate")) {
    treat <- as.integer(key == levels[[2]])
    core_estimand <- "ate"
  } else {
    # A focal estimand keeps the focal group at weight one, which the core does
    # for its level one. Encoding the focal level as one and solving the
    # treated-focal problem covers both the treated and the untreated targets.
    treat <- as.integer(key == focal)
    core_estimand <- "att"
  }

  result <- solve_ipt(
    covs,
    treat,
    s,
    core_estimand,
    method@link,
    ipt_options(method)
  )

  # The focal path encodes the focal level as one and solves the treated-focal
  # problem, so the re-evaluation hook tilts to that level; the average treatment
  # effect ignores the focal index.
  focal_idx <- if (identical(estimand, "ate")) 0L else 1L
  psi_fn <- make_ipt_psi_fn(
    covs,
    treat,
    focal_idx,
    s,
    core_estimand,
    method@link
  )
  weights_eval <- make_ipt_weights_eval(
    covs,
    treat,
    focal_idx,
    s,
    core_estimand,
    method@link
  )

  assemble_ipt(result, prepared, psi_fn, weights_eval)
}

fit_ipt_categorical <- function(method, prepared) {
  s <- prepared$sampling_weights
  covs <- cbind(1, prepared$matrix)
  levels <- prepared$exposure_levels
  key <- prepared$exposure_key
  estimand <- prepared$estimand
  focal <- prepared$focal_level

  treat_idx <- match(key, levels) - 1L

  if (identical(estimand, "ate")) {
    # The core requires a focal index in range even though the average treatment
    # effect ignores it.
    focal_idx <- 0L
    core_estimand <- "ate"
  } else {
    focal_idx <- match(focal, levels) - 1L
    core_estimand <- "att"
  }

  result <- solve_ipt_multi(
    covs,
    as.integer(treat_idx),
    as.integer(focal_idx),
    s,
    core_estimand,
    method@link,
    ipt_options(method)
  )

  psi_fn <- make_ipt_psi_fn(
    covs,
    treat_idx,
    focal_idx,
    s,
    core_estimand,
    method@link
  )
  weights_eval <- make_ipt_weights_eval(
    covs,
    treat_idx,
    focal_idx,
    s,
    core_estimand,
    method@link
  )

  assemble_ipt(result, prepared, psi_fn, weights_eval)
}

# A closure re-evaluating the tilting estimating functions at new coefficients,
# over the Rust eval entrypoint and the solve inputs captured here. The captured
# design matrix is the memory cost the optional container accepts.
make_ipt_psi_fn <- function(covs, treat_idx, focal, s, estimand, link) {
  force(covs)
  force(treat_idx)
  force(focal)
  force(s)
  force(estimand)
  force(link)
  function(theta) {
    eval_psi_ipt(
      as.numeric(theta),
      covs,
      as.integer(treat_idx),
      as.integer(focal),
      s,
      estimand,
      link
    )
  }
}

# A closure re-evaluating the tilting weights at new coefficients, over the Rust
# eval entrypoint and the solve inputs captured here. The entrypoint returns the
# raw M-estimator weights, the scale the container stores, so the reported
# per-group normalization is applied on top of it rather than inside it.
make_ipt_weights_eval <- function(covs, treat_idx, focal, s, estimand, link) {
  force(covs)
  force(treat_idx)
  force(focal)
  force(s)
  force(estimand)
  force(link)
  function(theta) {
    eval_weights_ipt(
      as.numeric(theta),
      covs,
      as.integer(treat_idx),
      as.integer(focal),
      s,
      estimand,
      link
    )
  }
}

# Normalize each exposure group to its estimand target sum and pack the fit
# result. For the average treatment effect the core's intercept moment fixes
# each group's sampling-weighted total at the whole-sample total, so rescaling
# to the group's own total is a deliberate change of reporting convention, a
# real per-group factor of roughly n_k / n, not drift removal. A focal estimand
# already places each group at the focal total, with the focal group at weight
# one, so there its rescaling only removes numerical drift. Either way the
# reported weights are not the raw M-estimator weights: the estimating-equations
# container is stored separately at the raw solution (see
# `ipt_estimating_equations`), and downstream inference reads the container, not
# these weights.
assemble_ipt <- function(result, prepared, psi_fn = NULL, weights_eval = NULL) {
  s <- prepared$sampling_weights
  focal <- prepared$focal_level
  groups <- prepared$groups

  weights_raw <- as.numeric(result$weights)
  w <- renormalize_group_weights(
    result$weights,
    s,
    groups,
    group_target_sums(s, groups, focal)
  )

  weights_fn <- if (is.null(weights_eval)) {
    NULL
  } else {
    make_weights_fn(weights_eval, reported = w, weights_raw = weights_raw)
  }

  list(
    weights = w,
    coefficients = as.numeric(result$coefs),
    converged = isTRUE(result$converged),
    interrupted = isTRUE(result$interrupted),
    iterations = as.integer(result$iterations),
    objective = result$grad_norm,
    solver_status = "newton",
    estimating_equations = ipt_estimating_equations(
      result,
      weights_raw,
      psi_fn,
      weights_fn
    ),
    groups = groups
  )
}

# Build the estimating-equations container from the matrices the core returned.
# Inverse probability tilting's per-unit estimating functions carry a
# weight-independent target term, so the container is stored at the raw
# M-estimator solution rather than at the renormalized reporting weights: the
# per-group normalization applied for output is not a linear scaling of the
# estimating functions, so scaling them would break the M-estimation
# representation. The stored functions therefore sum to zero column by column at
# the fitted parameters. `weight_jacobian` is the derivative of the raw weights,
# recorded on `weights_raw` so a consumer can rescale to the reported convention.
ipt_estimating_equations <- function(
  result,
  weights_raw,
  psi_fn = NULL,
  weights_fn = NULL
) {
  balancing_estimating_equations(
    parameters = as.numeric(result$coefs),
    psi = result$psi,
    jacobian = result$jac,
    weight_jacobian = result$dw_dbeta,
    weights_raw = weights_raw,
    psi_fn = psi_fn,
    weights_fn = weights_fn
  )
}
