# The variance for a balancing fit is a stacked M-estimator, and this file
# assembles that stack as a single estimating-function closure that deli
# differentiates and sandwiches. Where the hand-assembled path writes each block
# of the bread out analytically, here the whole system is written once as the
# estimating functions themselves and `deli::compute_sandwich()` finite
# differences the bread from them. Nothing is re-solved: every parameter enters
# at the value its own fit already found, and the closure only re-evaluates the
# estimating functions around that point.
#
# The stack is ordered [theta_w | beta | mu0 | mu1 | contrasts]. The weight
# parameters come first because everything downstream depends on them and
# nothing upstream does. The outcome-model coefficients follow, coupled to the
# weight parameters through the weights the score carries. The marginal means
# read off the outcome model with the exposure fixed to each level, and the
# effect contrasts close the system. Making each contrast a parameter of the
# stack is what removes the delta method from the caller: the contrast's
# standard error is already on the diagonal of the returned covariance.

#' The deli-backed stacked sandwich for a balancing fit
#'
#' Assembles the stacked estimating functions for the weight parameters, the
#' outcome model, the marginal means, and the effect contrasts, then returns
#' the sandwich covariance of the whole system.
#'
#' @param container The fit's [balancing_estimating_equations].
#' @param outcome_mod The fitted weighted marginal outcome model.
#' @param frame The data frame holding the exposure.
#' @param exposure_name The exposure column name.
#' @param sampling_weights The fit's sampling weights, or `NULL`.
#'
#' @return A list with `theta`, the stacked parameter vector, and `vcov`, its
#'   covariance on the standard-error scale, both named by stacked block order.
#'
#' @noRd
ipw_deli_sandwich <- function(
  container,
  outcome_mod,
  frame,
  exposure_name,
  sampling_weights = NULL
) {
  n <- nrow(frame)
  family <- stats::family(outcome_mod)
  continuous <- is_gaussian_outcome(outcome_mod)
  distribution <- deli_distribution(family)
  outcome <- resolve_outcome_response(outcome_mod)
  design <- stats::model.matrix(outcome_mod)
  offset <- outcome_mod$offset

  # The marginal-mean equations predict the outcome model with the exposure
  # fixed to each level, so they reuse the same fixed-exposure designs the
  # hand-assembled path builds. An offset is part of the linear predictor rather
  # than of the design, so it is carried alongside and added to eta.
  levels <- sort(unique(frame[[exposure_name]]))
  design0 <- fixed_exposure_pieces(
    outcome_mod,
    frame,
    exposure_name,
    levels[1],
    offset = offset
  )
  design1 <- fixed_exposure_pieces(
    outcome_mod,
    frame,
    exposure_name,
    levels[2],
    offset = offset
  )

  # Every parameter enters at the value its own fit already produced: the weight
  # parameters from the container, the coefficients from the outcome model, and
  # the means and contrasts as plug-in values of those two.
  weight_parameters <- container@parameters
  p <- length(weight_parameters)
  coefficients <- stats::coef(outcome_mod)
  q <- length(coefficients)
  mu0 <- mean(design0$mu)
  mu1 <- mean(design1$mu)
  contrasts <- ipw_contrast_values(mu0, mu1, continuous)
  effects <- ipw_contrast_names(continuous)
  k <- length(effects)

  theta <- c(weight_parameters, coefficients, mu0, mu1, contrasts)
  names(theta) <- c(
    paste0("theta_w", seq_len(p)),
    paste0("beta_", colnames(design)),
    "mu0",
    "mu1",
    effects
  )

  # The sampling weights compose multiplicatively onto the reported balancing
  # weights, which is the scale the outcome model was fitted at. Holding them
  # fixed here is correct: they are a design quantity, not an estimate.
  sampling <- sampling_weights %||% rep(1, n)

  # The reported weights and the weight-parameter estimating functions both come
  # from the container, so no method's math is restated. Differentiating the
  # score block through `weights_fn()` is what propagates the uncertainty in the
  # weights into the effect standard errors.
  stacked_equations <- function(theta) {
    weight_theta <- theta[seq_len(p)]
    beta <- theta[p + seq_len(q)]
    mean0 <- theta[[p + q + 1L]]
    mean1 <- theta[[p + q + 2L]]
    contrast_theta <- theta[p + q + 2L + seq_len(k)]

    weights <- as.numeric(container@weights_fn(weight_theta)) * sampling
    score <- deli::ee_glm(
      beta,
      X = design,
      y = outcome,
      distribution = distribution,
      link = family$link,
      weights = weights,
      offset = offset
    )

    eta0 <- as.numeric(design0$design %*% beta)
    eta1 <- as.numeric(design1$design %*% beta)
    if (!is.null(offset)) {
      eta0 <- eta0 + offset
      eta1 <- eta1 + offset
    }

    # The contrasts are deterministic functions of the two means, so their rows
    # are the same value for every unit. They contribute nothing to the meat at
    # the solution, where that value is zero, and everything to the bread, which
    # is where the delta method they replace used to be applied.
    contrast_rows <- matrix(
      ipw_contrast_values(mean0, mean1, continuous) - contrast_theta,
      nrow = k,
      ncol = n
    )

    rbind(
      t(container@psi_fn(weight_theta)),
      score,
      family$linkinv(eta0) - mean0,
      family$linkinv(eta1) - mean1,
      contrast_rows
    )
  }

  # A central difference trades truncation error, of order the step squared,
  # against cancellation error, of order the double epsilon over the step; deli's
  # 1e-9 default sits far into the cancellation regime, where agreement with the
  # analytic bread is near 2e-7 rather than the 2e-10 a 1e-6 step reaches.
  covariance <- deli::compute_sandwich(
    stacked_equations,
    theta,
    deriv_method = "capprox",
    dx = 1e-6,
    allow_pinv = FALSE
  ) /
    n

  # A bread holding missing values is warned about and answered with `NULL`
  # rather than an error, which would otherwise surface much later as a
  # complaint about dimnames applied to a non-array. Name the real cause here
  # instead, at the point where it is still legible.
  if (!is.matrix(covariance)) {
    abort(
      c(
        "The stacked variance could not be computed for this outcome model.",
        x = "The stacked estimating functions are not finite at the fitted parameters.",
        i = "See the inference vignette for a bootstrap workflow."
      ),
      error_class = "balancing_ipw_unsupported_error"
    )
  }
  dimnames(covariance) <- list(names(theta), names(theta))

  list(theta = theta, vcov = covariance)
}

# The effect contrasts of the two marginal means, and the labels the estimates
# table reports them under. They are parameters of the stack rather than a
# post-hoc transformation, so the values and the names are needed separately:
# the values to seed the stack, the names to label its blocks.
ipw_contrast_values <- function(mu0, mu1, continuous) {
  if (continuous) {
    return(mu1 - mu0)
  }
  c(
    mu1 - mu0,
    log(mu1) - log(mu0),
    log(mu1 / (1 - mu1)) - log(mu0 / (1 - mu0))
  )
}

ipw_contrast_names <- function(continuous) {
  if (continuous) "diff" else c("rd", "log(rr)", "log(or)")
}

# deli names distributions the way Python delicatessen does, which agrees with
# R's family names for the families that matter here once case is normalized;
# only the inverse gaussian is spelled differently. A family deli does not know
# raises deli's own error rather than being silently mapped onto a different
# variance function.
#
# The gamma and negative binomial estimating equations are the exception that
# has to be caught here rather than there. Both estimate a dispersion parameter
# alongside the coefficients, so they read the last element of the coefficient
# vector as a log dispersion and return an extra row. Passed a plain coefficient
# vector they would return a wrong-shaped block built from a misread parameter,
# which no downstream check would notice.
deli_distribution <- function(family, call = rlang::caller_env()) {
  distribution <- switch(
    family$family,
    inverse.gaussian = "inverse_normal",
    tolower(family$family)
  )
  if (distribution %in% c("gamma", "negative_binomial", "nb")) {
    abort(
      c(
        "{.fun ipw} cannot compute a stacked variance for a {.val {family$family}} outcome model.",
        x = "Its estimating equation estimates a dispersion parameter alongside the coefficients, which the stacked system does not carry.",
        i = "See the inference vignette for a bootstrap workflow."
      ),
      error_class = "balancing_ipw_unsupported_error",
      call = call
    )
  }
  distribution
}
