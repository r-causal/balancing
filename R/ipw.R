# balancing registers a method on propensity's ipw() generic so that a fitted
# balancing object can drive the same bring-your-own-model workflow as a
# propensity score fit. The variance is a stacked M-estimator: the weight
# parameters solve the estimating equations the fit carries, the outcome model
# supplies its score equations, and the marginal means supply their own, all
# sandwiched together so the standard errors account for having estimated the
# weights. R performs only generic matrix algebra with the container pieces and
# the outcome-model family functions; it never re-derives the balancing math.

#' Inverse probability weighting for a balancing fit
#'
#' @description
#' A [balancing] fit registers a method on [propensity::ipw()], so a set of
#' balancing weights drives the same bring-your-own-model workflow as a
#' propensity score fit. You supply the fit and a weighted outcome model, and
#' `ipw()` returns causal effect estimates with standard errors that account for
#' having estimated the weights.
#'
#' @details
#' The point estimates are the g-computation marginal means: the outcome model
#' is predicted with the exposure fixed to each level and averaged. For a binary
#' outcome the method returns the risk difference (`rd`), the log risk ratio
#' (`log(rr)`), and the log odds ratio (`log(or)`); for a continuous outcome it
#' returns the difference in means (`diff`).
#'
#' The standard errors come from a stacked M-estimator. The stacked parameter
#' vector holds the weight parameters, the outcome-model coefficients, and the
#' marginal means. Its estimating functions are the weight-parameter estimating
#' equations the fit carries, the outcome-model score equations, and the
#' marginal-mean equations. The sandwich variance of this system propagates the
#' uncertainty from estimating the weights into the effect standard errors, which
#' a variance that treats the weights as fixed would understate.
#'
#' The method is available only for fits whose weights solve smooth estimating
#' equations with a binary exposure: the estimating-equation family (entropy
#' balancing, inverse probability tilting, and the just-identified covariate
#' balancing propensity score) with exact balance. Any other fit, including a
#' tolerance-relaxed fit, an over-identified or quadratic-program fit, or a
#' categorical or continuous exposure, raises `balancing_ipw_unsupported_error`
#' and points to the bootstrap workflow described in the inference vignette.
#'
#' The outcome model must be the marginal model whose only predictor is the
#' exposure, such as `y ~ exposure`. The stacked variance is derived for that
#' form, where the g-computation means reduce to the weighted group means. A
#' covariate-adjusted outcome model raises `balancing_ipw_input_error`; use the
#' bootstrap workflow in the inference vignette for those models.
#'
#' @references
#' Kostouraki A, Hajage D, Rachet B, et al. On variance estimation of the
#' inverse probability-of-treatment weighting estimator: A tutorial for
#' different types of propensity score weights. *Statistics in Medicine*.
#' 2024;43(13):2672-2694. \doi{10.1002/sim.10078}
#'
#' @param ps_mod A [balancing] fit that produced the weights.
#' @param outcome_mod A weighted marginal outcome model of class [stats::glm()]
#'   or [stats::lm()], fitted with the balancing weights and with the exposure as
#'   its only predictor.
#' @param .data The data frame holding the exposure and outcome. If `NULL`, the
#'   values are taken from the outcome model frame.
#' @param estimand The causal estimand. If `NULL`, the fit's estimand is used.
#'   Supplying an estimand that disagrees with the fit raises
#'   `balancing_estimand_error`.
#' @param conf_level The confidence level for the intervals. Default `0.95`.
#' @param ... Ignored, for compatibility with the generic.
#'
#' @return An object of class `ipw`, the shared return contract of
#'   [propensity::ipw()].
#'
#' @examples
#' n <- 200
#' x1 <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.5 * x1))
#' y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.3 * x1))
#' df <- data.frame(exposure = z, x1 = x1, y = y)
#'
#' fit <- balance(df, exposure, x1, method = bw_entropy(), estimand = "ate")
#' df$.wts <- as.numeric(weights(fit))
#' outcome_mod <- glm(y ~ exposure, data = df, family = binomial(), weights = .wts)
#'
#' propensity::ipw(fit, outcome_mod)
#'
#' @name ipw.balancing
#' @importFrom propensity ipw
#' @importFrom stats getCall
NULL

# Expose the originating call through the standard model accessor so tools that
# summarize a fit, including the ipw() print method, can label it. Without a
# method the S7 object is not subsettable and the accessor would fail.
getCall_generic <- new_external_generic("stats", "getCall", "x")

method(getCall_generic, balancing) <- function(x, ...) {
  x@call
}

propensity_ipw <- new_external_generic("propensity", "ipw", "ps_mod")

method(propensity_ipw, balancing) <- function(
  ps_mod,
  outcome_mod,
  .data = NULL,
  estimand = NULL,
  conf_level = 0.95,
  ...
) {
  container <- ps_mod@estimating_equations
  if (is.null(container)) {
    abort_ipw_unsupported(reason = "no_equations")
  }
  if (!identical(ps_mod@exposure_type, "binary")) {
    abort_ipw_unsupported(
      reason = "exposure_type",
      exposure_type = ps_mod@exposure_type
    )
  }

  estimand <- resolve_ipw_estimand(estimand, ps_mod@estimand)

  exposure_name <- ps_mod@exposure
  validate_ipw_outcome_model(outcome_mod, exposure_name)

  frame <- if (is.null(.data)) stats::model.frame(outcome_mod) else .data
  if (!exposure_name %in% names(frame)) {
    abort(
      c(
        "The exposure {.val {exposure_name}} is not in the outcome model frame.",
        i = "Supply {.arg .data} containing the exposure column when the outcome formula transforms it."
      ),
      error_class = "balancing_ipw_input_error"
    )
  }
  levels <- sort(unique(frame[[exposure_name]]))
  if (length(levels) != 2L) {
    abort(
      c(
        "{.fun ipw} supports binary exposures only.",
        x = "The exposure {.val {exposure_name}} has {length(levels)} observed level{?s} in the data."
      ),
      error_class = "balancing_ipw_input_error"
    )
  }

  weights <- as.numeric(weights(ps_mod))
  outcome <- as.numeric(stats::model.response(stats::model.frame(outcome_mod)))
  family <- stats::family(outcome_mod)
  continuous <- is_gaussian_outcome(outcome_mod)

  # The fitted outcome model supplies its own design and working residuals; the
  # fixed-exposure designs supply the marginal-mean equations. All of these come
  # from generic model-matrix and family machinery.
  fitted <- outcome_model_pieces(outcome_mod, family)
  design0 <- fixed_exposure_pieces(outcome_mod, frame, exposure_name, levels[1])
  design1 <- fixed_exposure_pieces(outcome_mod, frame, exposure_name, levels[2])

  mu0 <- mean(design0$mu)
  mu1 <- mean(design1$mu)

  covariance <- stacked_sandwich(
    container = container,
    weights = weights,
    outcome = outcome,
    fitted = fitted,
    design0 = design0,
    design1 = design1,
    mu0 = mu0,
    mu1 = mu1
  )

  estimates <- ipw_estimates(
    mu0 = mu0,
    mu1 = mu1,
    covariance = covariance,
    conf_level = conf_level,
    continuous = continuous
  )

  structure(
    list(
      estimand = estimand,
      ps_mod = ps_mod,
      outcome_mod = outcome_mod,
      estimates = estimates
    ),
    class = "ipw"
  )
}

# The unsupported condition is shared by every configuration ipw() cannot
# handle, so both the missing-equations path and the non-binary path route
# through one message that names the reason and points to the bootstrap
# workflow.
abort_ipw_unsupported <- function(
  reason = c("no_equations", "exposure_type"),
  exposure_type = NULL,
  call = rlang::caller_env()
) {
  reason <- rlang::arg_match(reason)
  detail <- if (identical(reason, "no_equations")) {
    c(
      x = "This fit's weights do not solve smooth estimating equations, so the stacked variance is unavailable.",
      i = "Estimating equations come from the estimating-equation family (entropy balancing, inverse probability tilting, just-identified covariate balancing propensity score) with exact balance."
    )
  } else {
    c(
      x = "This fit has a {exposure_type} exposure, and only binary exposures are supported.",
      i = "The stacked variance is derived for a binary exposure."
    )
  }
  abort(
    c(
      "{.fun ipw} cannot compute a stacked variance for this balancing fit.",
      detail,
      i = "See the inference vignette for a bootstrap workflow."
    ),
    error_class = "balancing_ipw_unsupported_error",
    call = call,
    .envir = environment()
  )
}

# The estimand either comes from the fit or must agree with it. A supplied
# estimand that contradicts the fit is a specification error naming the knob.
resolve_ipw_estimand <- function(
  estimand,
  fit_estimand,
  call = rlang::caller_env()
) {
  if (is.null(estimand)) {
    return(fit_estimand)
  }
  if (!identical(estimand, fit_estimand)) {
    abort(
      c(
        "The requested {.arg estimand} does not match the fit.",
        x = "The fit targets {.val {fit_estimand}}.",
        x = "You requested {.val {estimand}}."
      ),
      error_class = "balancing_estimand_error",
      call = call
    )
  }
  estimand
}

is_gaussian_outcome <- function(outcome_mod) {
  if (inherits(outcome_mod, "glm")) {
    return(identical(stats::family(outcome_mod)$family, "gaussian"))
  }
  # A plain lm is a linear model.
  TRUE
}

# The stacked variance is derived for the marginal outcome model, whose only
# predictor is the exposure. That form makes the g-computation means equal the
# weighted group means and makes the per-group score sums vanish, which the
# variance relies on. A covariate-adjusted model would return a silently wrong
# standard error, so the contract is validated rather than trusted.
validate_ipw_outcome_model <- function(
  outcome_mod,
  exposure_name,
  call = rlang::caller_env()
) {
  if (!inherits(outcome_mod, c("glm", "lm"))) {
    abort(
      c(
        "{.arg outcome_mod} must be a fitted outcome model.",
        i = "Supply a model of class {.cls glm} or {.cls lm}.",
        x = "{.arg outcome_mod} has class {.cls {class(outcome_mod)}}."
      ),
      error_class = "balancing_ipw_input_error",
      call = call
    )
  }
  term_labels <- attr(stats::terms(outcome_mod), "term.labels")
  if (!identical(term_labels, exposure_name)) {
    abort(
      c(
        "{.arg outcome_mod} must be the marginal outcome model.",
        i = "The exposure {.val {exposure_name}} must be its only predictor.",
        i = "See the inference vignette for a bootstrap workflow with covariate-adjusted outcome models."
      ),
      error_class = "balancing_ipw_input_error",
      call = call
    )
  }
}

# The fitted outcome model's design, working residual, and working weight. The
# working residual is the score contribution per unit and the working weight is
# the Fisher-scoring weight, both from the family's mean and variance functions,
# so the assembly stays generic across families.
outcome_model_pieces <- function(outcome_mod, family) {
  design <- stats::model.matrix(outcome_mod)
  eta <- as.numeric(design %*% stats::coef(outcome_mod))
  mu <- family$linkinv(eta)
  mu_eta <- family$mu.eta(eta)
  variance <- family$variance(mu)
  list(
    design = design,
    mu = mu,
    mu_eta = mu_eta,
    variance = variance
  )
}

# The outcome design and predictions with the exposure fixed to one level, for
# the marginal-mean equations. The model's terms, contrasts, and factor levels
# are reused so the columns line up with the fitted coefficients.
fixed_exposure_pieces <- function(outcome_mod, frame, exposure_name, level) {
  frame[[exposure_name]] <- level
  terms <- stats::delete.response(stats::terms(outcome_mod))
  model_frame <- stats::model.frame(terms, frame, xlev = outcome_mod$xlevels)
  design <- stats::model.matrix(
    terms,
    model_frame,
    contrasts.arg = outcome_mod$contrasts
  )
  eta <- as.numeric(design %*% stats::coef(outcome_mod))
  family <- stats::family(outcome_mod)
  list(
    design = design,
    eta = eta,
    mu = family$linkinv(eta),
    mu_eta = family$mu.eta(eta)
  )
}

# The stacked sandwich covariance of (theta, beta, mu0, mu1). The bread is the
# derivative of the per-unit estimating functions with respect to the stacked
# parameters; the meat is their average outer product. The weight block's rows
# come from the container, the outcome-score rows carry the weight derivatives
# so the weight uncertainty propagates, and the marginal-mean rows close the
# system. Only the marginal-mean covariance block is returned, since the effect
# contrasts touch only those parameters.
stacked_sandwich <- function(
  container,
  weights,
  outcome,
  fitted,
  design0,
  design1,
  mu0,
  mu1
) {
  psi <- container@psi
  jacobian <- container@jacobian
  weight_jacobian <- container@weight_jacobian

  design <- fitted$design
  n <- nrow(design)
  p <- ncol(psi)
  q <- ncol(design)

  # The outcome-score derivative uses the expected (Fisher-scoring) information,
  # which equals the observed information for a canonical link such as the logit
  # for a binomial outcome or the identity for a gaussian one.
  residual <- (outcome - fitted$mu) * fitted$mu_eta / fitted$variance
  working_weight <- weights * fitted$mu_eta^2 / fitted$variance
  score <- (weights * residual) * design

  moment0 <- design0$mu - mu0
  moment1 <- design1$mu - mu1

  stacked <- cbind(psi, score, moment0, moment1)
  meat <- crossprod(stacked) / n

  # The score, meat, and outcome-information blocks use the weights the outcome
  # model was fitted with, which compose the sampling weights onto the reported
  # balancing weights. The container's weight derivatives are stored at each
  # method's own weight scale, so the coupling block rescales each unit's row by
  # the ratio of the model weight to that scale. The ratio is the sampling weight
  # for methods whose container is already at the reported scale, and the
  # sampling weight times the per-group reporting factor otherwise. The
  # per-group factor's own parameter derivative cancels for the marginal outcome
  # model, so the level ratio is the whole correction.
  coupling <- weight_jacobian * (weights / container@weights_raw)

  size <- p + q + 2L
  theta <- seq_len(p)
  beta <- p + seq_len(q)
  index0 <- p + q + 1L
  index1 <- p + q + 2L

  bread <- matrix(0, size, size)
  bread[theta, theta] <- jacobian / n
  bread[beta, theta] <- crossprod(design * residual, coupling) / n
  bread[beta, beta] <- -crossprod(design * working_weight, design) / n
  bread[index0, beta] <- colSums(design0$mu_eta * design0$design) / n
  bread[index1, beta] <- colSums(design1$mu_eta * design1$design) / n
  bread[index0, index0] <- -1
  bread[index1, index1] <- -1

  bread_inv <- solve(bread)
  covariance <- bread_inv %*% meat %*% t(bread_inv) / n
  covariance[c(index0, index1), c(index0, index1), drop = FALSE]
}

# The effect rows and their delta-method standard errors. Each effect is a
# smooth contrast of the two marginal means, so its variance is the gradient
# quadratic form against the marginal-mean covariance block.
ipw_estimates <- function(mu0, mu1, covariance, conf_level, continuous) {
  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)

  effect_row <- function(effect, estimate, gradient) {
    std_err <- sqrt(as.numeric(t(gradient) %*% covariance %*% gradient))
    z <- estimate / std_err
    data.frame(
      effect = effect,
      estimate = estimate,
      std.err = std_err,
      z = z,
      ci.lower = estimate - z_value * std_err,
      ci.upper = estimate + z_value * std_err,
      conf.level = conf_level,
      p.value = 2 * (1 - stats::pnorm(abs(z)))
    )
  }

  difference <- effect_row(
    if (continuous) "diff" else "rd",
    mu1 - mu0,
    c(-1, 1)
  )
  if (continuous) {
    return(difference)
  }

  log_rr <- effect_row(
    "log(rr)",
    log(mu1) - log(mu0),
    c(-1 / mu0, 1 / mu1)
  )
  log_or <- effect_row(
    "log(or)",
    log(mu1 / (1 - mu1)) - log(mu0 / (1 - mu0)),
    c(-1 / (mu0 * (1 - mu0)), 1 / (mu1 * (1 - mu1)))
  )

  rbind(difference, log_rr, log_or)
}
