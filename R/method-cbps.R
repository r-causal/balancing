# Covariate balancing propensity score: a propensity model whose parameters are
# chosen to satisfy covariate balancing moment conditions. In the just-identified
# form the number of moment conditions equals the number of parameters, so the
# balancing conditions hold exactly and the achieved balance matches the
# requested moments. In the over-identified form the model score equations are
# stacked onto the balancing conditions and a generalized-method-of-moments
# criterion is minimized, so balance is approximate and the fit reports the
# criterion value rather than estimating equations. The method specification
# carries only tuning parameters; fit_method() prepares the design, calls the
# Rust solver, applies the estimand's group-sum normalization, and passes the raw
# estimating equations through unchanged for the just-identified discrete form.

#' Covariate balancing propensity score
#'
#' `bw_cbps()` specifies the covariate balancing propensity score for [balance()]. A
#' propensity model is fit so that its parameters satisfy covariate balancing
#' moment conditions rather than the maximum-likelihood score alone. In the
#' just-identified form the moment conditions equal the parameter count, so
#' balance on the requested moments is exact by construction and the weights
#' solve smooth estimating equations. The covariate balancing propensity score
#' supports binary, categorical, and continuous exposures, and is the only method
#' that supports the overlap estimand `"ato"` for a binary exposure.
#'
#' @details
#' For a binary exposure the just-identified fit weights a unit by a function of
#' its modeled propensity determined by the estimand: the inverse propensity for
#' the average treatment effect, the inverse odds for the average treatment
#' effect on the treated, and the overlap factor for the overlap estimand. With
#' mean balance and the logit link the just-identified average-treatment-
#' effect-on-the-treated fit solves the same treated-target moment conditions as
#' entropy balancing and inverse probability tilting, so the three methods
#' produce the same weights.
#'
#' Setting `over_identified = TRUE` stacks the propensity model's own score
#' equations onto the balancing conditions and minimizes a generalized-method-of-
#' moments criterion. Balance is then approximate, the fit records the criterion
#' value on its objective, and it supplies no estimating equations. `two_step`
#' selects the two-step weighting matrix for that criterion; it has no effect on
#' the just-identified fit and is warned and ignored there. Every just-identified
#' discrete fit, the overlap estimand included, supplies estimating equations;
#' only the over-identified form and a continuous exposure do not.
#'
#' For a continuous exposure the covariate balancing conditions require the
#' weighted exposure mean to match the sample mean and the weighted covariance
#' between the exposure and every covariate to vanish. The weights that meet these
#' conditions with the least departure from uniformity are the minimum-divergence
#' exponential tilt, the same reweighting the continuous form of entropy balancing
#' uses. This is the nonparametric reading of covariate balancing for a continuous
#' exposure and departs from the parametric generalized propensity score of the
#' Fong, Hazlett, and Imai reference, which derives the weights from a Gaussian
#' density ratio and can be unstable; the two share the balancing conditions but
#' not the weight family.
#'
#' @param over_identified Whether to add the propensity model's score equations
#'   and minimize the generalized-method-of-moments criterion. `FALSE` fits the
#'   just-identified form, whose balance is exact.
#' @param two_step Whether to use the two-step weighting matrix for the
#'   over-identified criterion, rather than the continuously updating criterion.
#'   Ignored, with a warning, when `over_identified` is `FALSE`.
#' @param link The propensity link, one of `"logit"`, `"probit"`, or
#'   `"cloglog"`. Binary exposures only.
#' @param convergence_tolerance The solver convergence tolerance.
#' @param max_iterations The maximum solver iterations, or `NULL` for the core
#'   default.
#' @param ... Reserved for future extensions; must be empty. Tuning parameters
#'   must be passed by name.
#'
#' @return A `bw_cbps` specification, a [balance_method].
#'
#' @references
#' Imai, K. and Ratkovic, M. (2014). Covariate balancing propensity score.
#' *Journal of the Royal Statistical Society: Series B (Statistical
#' Methodology)*, 76(1), 243-263.
#'
#' Fong, C., Hazlett, C., and Imai, K. (2018). Covariate balancing propensity
#' score for a continuous treatment: Application to the efficacy of political
#' advertisements. *The Annals of Applied Statistics*, 12(1), 156-177.
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
#' fit <- balance(df, exposure, c(x1, x2), method = bw_cbps())
#' fit
#'
#' @export
bw_cbps <- new_class(
  "bw_cbps",
  parent = estimating_equation_method,
  properties = list(
    over_identified = class_logical,
    two_step = class_logical,
    link = class_character
  ),
  constructor = function(
    over_identified = FALSE,
    two_step = TRUE,
    link = c("logit", "probit", "cloglog"),
    convergence_tolerance = 1e-10,
    max_iterations = NULL,
    ...
  ) {
    rlang::check_dots_empty()
    link <- rlang::arg_match(link)
    if (!is.null(max_iterations)) {
      max_iterations <- vctrs::vec_cast(
        max_iterations,
        integer(),
        x_arg = "max_iterations"
      )
    }
    new_object(
      S7_object(),
      over_identified = over_identified,
      two_step = two_step,
      link = link,
      convergence_tolerance = convergence_tolerance,
      max_iterations = max_iterations
    )
  },
  validator = function(self) {
    if (length(self@over_identified) != 1 || is.na(self@over_identified)) {
      return("@over_identified must be a single logical value")
    }
    if (length(self@two_step) != 1 || is.na(self@two_step)) {
      return("@two_step must be a single logical value")
    }
    if (
      !is.null(self@convergence_tolerance) && self@convergence_tolerance <= 0
    ) {
      return("@convergence_tolerance must be positive")
    }
    if (!is.null(self@max_iterations) && self@max_iterations < 0) {
      return("@max_iterations must be non-negative")
    }
  }
)

method(method_label, bw_cbps) <- function(method) {
  "Covariate balancing propensity score"
}

method(supported_exposure_types, bw_cbps) <- function(method) {
  c("binary", "categorical", "continuous")
}

method(supported_estimands, bw_cbps) <- function(method, exposure_type) {
  switch(
    exposure_type,
    # The overlap estimand is legal only for a binary exposure.
    binary = c("ate", "att", "atu", "ato"),
    categorical = c("ate", "att"),
    continuous = "ate"
  )
}

# The just-identified discrete form supplies smooth estimating equations. The
# over-identified form minimizes a generalized-method-of-moments criterion and a
# continuous exposure balances the exposure-covariate covariance through an
# exponential tilt, so neither carries estimating equations.
method(supports_estimating_equations, bw_cbps) <- function(
  method,
  ...,
  exposure_type = NULL,
  constraints = NULL
) {
  rlang::check_dots_empty()
  if (method@over_identified) {
    return(FALSE)
  }
  if (identical(exposure_type, "continuous")) {
    return(FALSE)
  }
  TRUE
}

# Assemble the Rust option list, dropping the tuning parameters left at the core
# default so the solver applies its own. The worker-thread count is resolved on
# the R side and passed down on every call.
cbps_options <- function(method) {
  options <- list(threads = resolve_threads())
  if (!is.null(method@convergence_tolerance)) {
    options$convergence_tolerance <- method@convergence_tolerance
  }
  if (!is.null(method@max_iterations)) {
    options$max_iterations <- as.integer(method@max_iterations)
  }
  options
}

method(fit_method, bw_cbps) <- function(method, prepared) {
  # The two-step weighting matrix belongs to the over-identified criterion.
  # Setting it while the fit is just-identified has no effect, so warn and
  # proceed rather than fail validation.
  if (!method@over_identified && !method@two_step) {
    warn(
      c(
        "{.arg two_step} applies only to the over-identified fit and is ignored.",
        i = "Set {.code over_identified = TRUE} in {.fn bw_cbps} to use the two-step weighting matrix."
      ),
      warning_class = "balancing_ignored_argument_warning"
    )
  }

  switch(
    prepared$exposure_type,
    binary = fit_cbps_binary(method, prepared),
    categorical = fit_cbps_categorical(method, prepared),
    continuous = fit_cbps_continuous(method, prepared)
  )
}

# The estimand string the core expects. The orchestrator stores the untreated
# target as "atu"; the core names it "atc".
cbps_core_estimand <- function(estimand) {
  switch(estimand, ate = "ate", att = "att", atu = "atc", ato = "ato")
}

fit_cbps_binary <- function(method, prepared) {
  s <- prepared$sampling_weights
  # The propensity model carries an intercept, so the design is the intercept
  # column followed by the standardized covariates. The just-identified fit
  # solves the same design for the model and the balance conditions.
  covs <- cbind(1, prepared$matrix)
  levels <- prepared$exposure_levels
  key <- prepared$exposure_key

  # The second level is the treated level; the estimand determines which group
  # keeps its base weights and which is reweighted toward it.
  treat <- as.integer(key == levels[[2]])
  core_estimand <- cbps_core_estimand(prepared$estimand)

  result <- solve_cbps(
    covs,
    covs,
    treat,
    s,
    core_estimand,
    method@link,
    method@over_identified,
    method@two_step,
    cbps_options(method)
  )

  # Only the just-identified fit produces a smooth estimating-equations
  # container, so the re-evaluation hook is built for that path alone. The
  # over-identified generalized-method-of-moments fit has no container.
  psi_fn <- if (isTRUE(method@over_identified)) {
    NULL
  } else {
    make_cbps_psi_fn(covs, treat, s, core_estimand, method@link)
  }

  assemble_cbps(result, method, prepared, psi_fn)
}

# A closure re-evaluating the binary just-identified estimating functions at new
# coefficients, over the Rust eval entrypoint and the solve inputs captured
# here. The captured design matrix is the memory cost the optional container
# accepts.
make_cbps_psi_fn <- function(covs, treat, s, estimand, link) {
  force(covs)
  force(treat)
  force(s)
  force(estimand)
  force(link)
  function(theta) {
    eval_psi_cbps(
      as.numeric(theta),
      covs,
      as.integer(treat),
      s,
      estimand,
      link
    )
  }
}

fit_cbps_categorical <- function(method, prepared) {
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

  result <- solve_cbps_multi(
    covs,
    as.integer(treat_idx),
    as.integer(focal_idx),
    s,
    core_estimand,
    method@link,
    cbps_options(method)
  )

  # The categorical just-identified fit coincides with per-level covariate
  # balancing, so its estimating functions share the tilting form; the hook
  # re-evaluates them through the shared tilting entrypoint.
  psi_fn <- make_ipt_psi_fn(
    covs,
    treat_idx,
    focal_idx,
    s,
    core_estimand,
    method@link
  )

  assemble_cbps(result, method, prepared, psi_fn)
}

fit_cbps_continuous <- function(method, prepared) {
  s <- prepared$sampling_weights
  # The continuous form balances the weighted exposure-covariate covariance
  # through an exponential tilt. The intercept column carries the exposure-mean
  # condition, which fixes the weighted exposure mean at the sample mean, and the
  # standardized covariate columns carry the covariance conditions; the core
  # centers the exposure internally.
  covs <- cbind(1, prepared$matrix)
  exposure <- as.numeric(prepared$exposure_vec)

  result <- solve_cbps_cont(
    covs,
    exposure,
    s,
    cbps_options(method)
  )

  # The continuous form supplies no smooth estimating equations, so it carries no
  # container and no re-evaluation hook.
  assemble_cbps(result, method, prepared, NULL)
}

# Normalize each exposure group to its estimand target sum and pack the fit
# result. The core returns the raw M-estimator weights; the reported weights
# renormalize each group to its estimand target total, a deliberate change of
# reporting convention rather than drift removal. For the average treatment
# effect and the overlap estimand each group is scaled to its own
# sampling-weighted total; for a focal estimand every group is scaled to the
# focal total, which leaves the focal group at its base weight. The
# estimating-equations container is stored separately at the raw solution (see
# `cbps_estimating_equations`), so downstream inference reads the container, not
# these renormalized weights.
assemble_cbps <- function(result, method, prepared, psi_fn = NULL) {
  s <- prepared$sampling_weights
  estimand <- prepared$estimand
  focal <- prepared$focal_level
  groups <- prepared$groups
  levels <- prepared$exposure_levels

  w <- result$weights
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

  # The over-identified fit reports its generalized-method-of-moments criterion;
  # the just-identified and continuous fits report the solver's merit value.
  objective <- if (isTRUE(method@over_identified) && !is.null(result$gmm_obj)) {
    result$gmm_obj
  } else {
    result$obj_value
  }

  # Only the binary over-identified generalized-method-of-moments fit is used
  # here; the categorical and continuous paths do not consume `over_identified`.
  over_identified <- isTRUE(method@over_identified) &&
    identical(prepared$exposure_type, "binary")

  solver_status <- if (over_identified) "lbfgs" else "newton"

  list(
    weights = w,
    coefficients = as.numeric(result$coefs),
    converged = isTRUE(result$converged),
    interrupted = isTRUE(result$interrupted),
    iterations = as.integer(result$iterations),
    objective = objective,
    solver_status = solver_status,
    estimating_equations = cbps_estimating_equations(result, psi_fn),
    # The over-identified criterion targets no per-constraint tolerance, so its
    # approximate balance must not raise the tolerance warning.
    approximate = over_identified,
    # The covariate balancing conditions equate the exposure arms directly, so the
    # tolerance geometry is arm-to-arm rather than the box family's pooled target.
    constraint_target = "arms",
    groups = groups
  )
}

# Build the estimating-equations container from the matrices the core returned.
# The just-identified discrete form returns the per-unit estimating functions,
# the Jacobian, and the weight derivatives at the raw M-estimator solution, where
# the per-unit functions sum to zero column by column. The over-identified and
# continuous forms return `NULL` for these fields, so the container is absent.
cbps_estimating_equations <- function(result, psi_fn = NULL) {
  psi <- result$psi
  jacobian <- result$jac
  weight_jacobian <- result$dw_dbeta
  if (is.null(psi) || is.null(jacobian) || is.null(weight_jacobian)) {
    return(NULL)
  }
  balancing_estimating_equations(
    parameters = as.numeric(result$coefs),
    psi = psi,
    jacobian = jacobian,
    weight_jacobian = weight_jacobian,
    weights_raw = as.numeric(result$weights),
    psi_fn = psi_fn
  )
}
