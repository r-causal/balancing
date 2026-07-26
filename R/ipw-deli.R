# The variance for a balancing fit is a stacked M-estimator, and this file
# assembles that stack as a single estimating-function closure that deli
# differentiates and sandwiches. Rather than writing each block of the bread out
# analytically, the whole system is written once as the estimating functions
# themselves and `deli::compute_sandwich()` finite differences the bread from
# them. Nothing is re-solved: every parameter enters at the value its own fit
# already found, and the closure only re-evaluates the estimating functions
# around that point.
#
# The stack is ordered [theta_w | beta | mu0 | mu1 | contrasts]. The weight
# parameters come first because everything downstream depends on them and
# nothing upstream does. The outcome-model coefficients follow, coupled to the
# weight parameters through the weights the score carries. The marginal means
# read off the outcome model with the exposure fixed to each level, and the
# effect contrasts close the system. Making each contrast a parameter of the
# stack is what removes the delta method from the caller: the contrast's
# standard error is already on the diagonal of the returned covariance.
#
# Two details of that system are there for the outcome models that adjust for
# covariates, and both reduce to what a marginal model already did.
#
# The marginal means standardize over the estimand's target population rather
# than over every unit. A marginal model predicts one value per exposure level,
# so which units are averaged over cannot matter; an adjusted model predicts a
# value per unit, and the population averaged over is then part of the estimand.
# The mean rows therefore carry the target population's indicator, which is
# every unit for a pooled estimand and the focal group for a focal one, times
# the sampling weights. This is the tilted g-computation the propensity package
# performs with the tilting function of its fitted score; a balancing fit
# carries no propensity model, so the population is read off the data and the
# indicator contributes nothing to the derivative.
#
# The weights the outcome score carries are the reported weights re-derived at
# the perturbed weight parameters, renormalization included. The reported scale
# carries each exposure group to a target total, and holding that per-group
# factor at the value the fit found would drop the part of the derivative that
# comes from the total itself moving. For a marginal model the dropped part
# contributes nothing, since the per-group weighted score sums vanish in every
# design direction and a per-group rescale therefore leaves the coefficients
# alone. For an adjusted model the covariate directions do not vanish, so the
# renormalization is applied again at each set of weight parameters.

#' The deli-backed stacked sandwich for a balancing fit
#'
#' Assembles the stacked estimating functions for the weight parameters, the
#' outcome model, the marginal means, and the effect contrasts, then returns
#' the sandwich covariance of the whole system.
#'
#' @param container The fit's [balancing_estimating_equations].
#' @param outcome_mod The fitted weighted outcome model.
#' @param frame The data frame holding the exposure.
#' @param exposure_name The exposure column name.
#' @param sampling_weights The fit's sampling weights, or `NULL`.
#' @param focal_level The fit's focal exposure level, or `NULL` for a pooled
#'   estimand. It names the target population the marginal means standardize
#'   over and the group total the reported weights are carried to.
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
  sampling_weights = NULL,
  focal_level = NULL
) {
  n <- nrow(frame)
  family <- stats::family(outcome_mod)
  continuous <- is_gaussian_outcome(outcome_mod)
  distribution <- deli_distribution(family)
  outcome <- resolve_outcome_response(outcome_mod)
  design <- stats::model.matrix(outcome_mod)
  offset <- outcome_mod$offset

  # The marginal-mean equations predict the outcome model with the exposure
  # fixed to each level, so they build a design per level from the model's own
  # terms. An offset is part of the linear predictor rather than of the design,
  # so it is carried alongside and added to eta.
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

  # The sampling weights compose multiplicatively onto the reported balancing
  # weights, which is the scale the outcome model was fitted at. Holding them
  # fixed here is correct: they are a design quantity, not an estimate.
  sampling <- sampling_weights %||% rep(1, n)

  # The standardization weight, one per unit: the target population's indicator
  # times the sampling weights. A pooled estimand targets every unit, so the
  # indicator is one throughout and the means are the sampling-weighted averages
  # of the fixed-exposure predictions; a focal estimand targets the focal group,
  # which is the treated group for the average effect on the treated and the
  # untreated group for the average effect on the untreated, since the fit
  # resolves the focal level to the group it holds fixed. Nothing here depends
  # on the parameters, so the standardization contributes no derivative of its
  # own.
  #
  # The exposure groups and the totals the reported weights are carried to come
  # from the same two facts, so they are built once here for the weight map
  # below.
  key <- as.character(frame[[exposure_name]])
  groups <- split(seq_len(n), key)
  targets <- group_target_sums(sampling, groups, focal_level)
  tilt <- if (is.null(focal_level)) {
    sampling
  } else {
    sampling * (key == focal_level)
  }

  # Every parameter enters at the value its own fit already produced: the weight
  # parameters from the container, the coefficients from the outcome model, and
  # the means and contrasts as plug-in values of those two. The means are the
  # standardized ones, which is the root of the mean rows the closure returns.
  weight_parameters <- container@parameters
  p <- length(weight_parameters)
  coefficients <- stats::coef(outcome_mod)
  q <- length(coefficients)
  mu0 <- sum(tilt * design0$mu) / sum(tilt)
  mu1 <- sum(tilt * design1$mu) / sum(tilt)
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

  # The reported weights and the weight-parameter estimating functions both come
  # from the container, so no method's math is restated. Differentiating the
  # score block through `weights_fn()` is what propagates the uncertainty in the
  # weights into the effect standard errors.
  #
  # The container's hook carries the reporting scale as the fixed per-group
  # factor the fit found, which is the scale its `weight_jacobian` describes and
  # the contract its own tests pin. Renormalizing the hook's result to the group
  # totals recovers the reported weight map itself, the one an outcome model
  # refitted at other weight parameters would have been given: the hook's fixed
  # factor is constant within a group, so it cancels between the numerator and
  # the group sum, leaving the raw weights carried to the target total. The
  # renormalization is a per-group rescale of the hook's output rather than
  # another crossing into the method, so it is cheap and it caches with the
  # values it is derived from.
  #
  # Those two hooks are the expensive part of the closure: each crosses into the
  # method's own evaluation entrypoint over the whole data set. Both are pure
  # functions of the weight block, and the finite difference presents the same
  # weight sub-vector many times over, because perturbing a coordinate outside
  # that block leaves the sub-vector exactly at its fitted value. Two cached
  # entries cover every repeat. The fitted sub-vector is pinned, since each of
  # the central difference's two sweeps returns to it in a long run and a single
  # most-recent entry would lose it in between; one further entry holds the most
  # recent perturbation, which a sweep asks for twice in succession. The keys
  # are short numeric vectors, so comparing them outright is cheaper than
  # hashing them.
  evaluate_hooks <- function(weight_theta) {
    list(
      key = weight_theta,
      weights = renormalize_group_weights(
        as.numeric(container@weights_fn(weight_theta)),
        sampling,
        groups,
        targets
      ),
      psi = t(container@psi_fn(weight_theta))
    )
  }
  base_hooks <- evaluate_hooks(as.numeric(weight_parameters))
  recent_hooks <- base_hooks
  hooks_at <- function(weight_theta) {
    if (identical(weight_theta, base_hooks$key)) {
      return(base_hooks)
    }
    if (!identical(weight_theta, recent_hooks$key)) {
      recent_hooks <<- evaluate_hooks(weight_theta)
    }
    recent_hooks
  }

  stacked_equations <- function(theta) {
    beta <- theta[p + seq_len(q)]
    mean0 <- theta[[p + q + 1L]]
    mean1 <- theta[[p + q + 2L]]
    contrast_theta <- theta[p + q + 2L + seq_len(k)]

    hooks <- hooks_at(as.numeric(theta[seq_len(p)]))
    weights <- hooks$weights * sampling
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
    # is what carries their standard errors without a delta method.
    contrast_rows <- matrix(
      ipw_contrast_values(mean0, mean1, continuous) - contrast_theta,
      nrow = k,
      ncol = n
    )

    # The mean rows are weighted by the standardization weight, so their root is
    # the mean of the fixed-exposure predictions over the target population
    # rather than over every unit. A marginal model predicts one value per
    # exposure level, which makes the weighted row a constant multiple of the
    # unweighted one and leaves the sandwich exactly where it was.
    rbind(
      hooks$psi,
      score,
      tilt * (family$linkinv(eta0) - mean0),
      tilt * (family$linkinv(eta1) - mean1),
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
# The quasibinomial family is the one deliberate rename. deli carries no quasi
# families and needs none: the quasibinomial variance function is the binomial
# one, so the estimating equations are the same equations. The dispersion the
# quasi family estimates never reaches the sandwich, which is built from the
# score alone, and scaling a score by a constant scales the bread by that
# constant and the meat by its square, leaving the sandwich unchanged. Mapping
# the family onto deli's binomial therefore reproduces the binomial answer
# rather than approximating it.
#
# The gamma and negative binomial estimating equations are the exception that
# has to be caught here rather than there. Both estimate a dispersion parameter
# alongside the coefficients, so they read the last element of the coefficient
# vector as a log dispersion and return an extra row. Passed a plain coefficient
# vector they would return a wrong-shaped block built from a misread parameter,
# which no downstream check would notice. A negative binomial fit spells its
# estimated dispersion into the family name itself, as "Negative Binomial(2)",
# so the refusal matches on the prefix; an equality test against a fixed spelling
# would never fire.
deli_distribution <- function(family, call = rlang::caller_env()) {
  distribution <- switch(
    family$family,
    inverse.gaussian = "inverse_normal",
    quasibinomial = "binomial",
    tolower(family$family)
  )
  estimates_dispersion <- distribution %in%
    c("gamma", "negative_binomial", "nb") ||
    startsWith(distribution, "negative binomial")
  if (estimates_dispersion) {
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
