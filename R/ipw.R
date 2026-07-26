# balancing registers a method on propensity's ipw() generic so that a fitted
# balancing object can drive the same bring-your-own-model workflow as a
# propensity score fit. The variance is a stacked M-estimator: the weight
# parameters solve the estimating equations the fit carries, the outcome model
# supplies its score equations, and the marginal means and effect contrasts
# supply their own, all sandwiched together so the standard errors account for
# having estimated the weights. This file validates the inputs and reads the
# effect table off that system; the system itself is assembled in ipw-deli.R,
# where the container's own hooks and the outcome-model family functions are the
# only sources of method math.

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
#' is predicted with the exposure fixed to each level and averaged over the
#' estimand's target population. For a binary
#' outcome the method returns the risk difference (`rd`), the log risk ratio
#' (`log(rr)`), and the log odds ratio (`log(or)`); for a continuous outcome it
#' returns the difference in means (`diff`).
#'
#' The standard errors come from a stacked M-estimator that the deli package
#' differentiates and sandwiches. The stacked parameter vector holds four
#' blocks: the
#' weight parameters, the outcome-model coefficients, the two marginal means,
#' and the effect contrasts. Their estimating functions are, in the same order,
#' the weight-parameter estimating equations the fit carries, re-evaluated at
#' new weight parameters through the hooks the fit's
#' [balancing_estimating_equations] container supplies; the outcome-model score,
#' from [deli::ee_glm()], carrying the balancing weights as the fit would have
#' reported them at those parameters, per-group renormalization included; the
#' marginal-mean equations, which predict the outcome model with the exposure
#' fixed to each level and standardize over the estimand's target population;
#' and one deterministic row per contrast, setting the contrast parameter
#' equal to its formula in the two means.
#' [deli::compute_sandwich()] differentiates that system at the fitted values
#' and returns its empirical sandwich covariance.
#'
#' Nothing is re-solved along the way. Every parameter enters at the value its
#' own fit already found, and the stacked estimating functions are only
#' re-evaluated around that point. Because each contrast is a parameter of the
#' stack rather than a transformation applied afterward, its standard error is
#' already on the diagonal of the joint covariance and no delta-method step
#' stands between the sandwich and the reported effects. The joint covariance
#' propagates the uncertainty from estimating the weights into the effect
#' standard errors, which a variance that treats the weights as fixed would
#' understate.
#'
#' The method is available only for fits whose weights solve smooth estimating
#' equations with a binary exposure: the estimating-equation family (entropy
#' balancing, inverse probability tilting, and the just-identified covariate
#' balancing propensity score) with exact balance. Any other fit, including an
#' entropy fit at a positive tolerance, an over-identified or quadratic-program
#' fit, or a categorical or continuous exposure, raises
#' `balancing_ipw_unsupported_error` and points to the bootstrap workflow
#' described in the inference vignette.
#'
#' The outcome model must carry the exposure among its predictors, and may
#' adjust for covariates alongside it, including in interactions with the
#' exposure: `y ~ exposure`, `y ~ exposure + x1 + x2`, and `y ~ exposure * x1`
#' are all supported. A model without an exposure term raises
#' `balancing_ipw_input_error`, since its two fixed-exposure predictions would
#' be the same prediction and every contrast it reported would be zero.
#'
#' Which population the marginal means are averaged over is part of the
#' estimand, and matters as soon as the outcome model adjusts for anything. The
#' marginal model of a binary exposure is saturated, one free parameter per
#' exposure level, so absent an offset it predicts a single value per level, its
#' marginal means are the weighted group means whatever link the family carries,
#' and no choice of population can change them. An adjusted model predicts a
#' value per unit, so the means are standardized over the estimand's target
#' population: every unit for a pooled estimand, and the focal group's units for
#' `"att"` or `"atc"`. Sampling weights, where the fit has them, weight that
#' average as well.
#'
#' The link enters the outcome-model score, and it enters exactly: the bread is
#' differentiated from the estimating functions themselves rather than read off
#' an information-matrix formula, so a non-canonical link such as probit or
#' cloglog is handled exactly rather than approximately, for an adjusted model
#' as much as for a marginal one.
#'
#' Two further conditions on the outcome model raise the same
#' `balancing_ipw_input_error`. Its
#' family must be binomial, quasibinomial, or gaussian, which includes a plain
#' [stats::lm()], since the reported effects are the contrasts derived for those
#' families' marginal means. And it must have been fitted with the weights the
#' fit produced, since the stacked variance differentiates the outcome-model
#' score through those weights: the model's weights are compared against the
#' fit's, per unit at a relative tolerance of 1e-6. Those are the weights
#' `weights(fit)` returns, which already carry the fit's sampling weights if it
#' has any. Sampling weights compose multiplicatively onto the balancing weights
#' and the stack holds them fixed, since they are a design quantity rather than
#' an estimate.
#'
#' An offset is supported, written either as an `offset()` term in the outcome
#' formula or passed through the model's `offset` argument. It is carried
#' through both the outcome-model score and the fixed-exposure linear
#' predictors, so the marginal means are the g-computation means with each
#' unit's offset held at its observed value.
#'
#' @references
#' Kostouraki A, Hajage D, Rachet B, et al. On variance estimation of the
#' inverse probability-of-treatment weighting estimator: A tutorial for
#' different types of propensity score weights. *Statistics in Medicine*.
#' 2024;43(13):2672-2694. \doi{10.1002/sim.10078}
#'
#' @param ps_mod A [balancing] fit that produced the weights.
#' @param outcome_mod A weighted outcome model of class [stats::glm()] or
#'   [stats::lm()], fitted with the balancing weights and carrying the exposure
#'   among its predictors. It may adjust for covariates alongside the exposure.
#' @param .data The data frame holding the exposure and outcome. If `NULL`, the
#'   values are taken from the outcome model frame.
#' @param estimand The causal estimand. If `NULL`, the fit's estimand is used.
#'   Supplying an estimand that disagrees with the fit raises
#'   `balancing_estimand_error`.
#' @param conf_level The confidence level for the intervals. Default `0.95`.
#' @param ... Ignored, for compatibility with the generic.
#'
#' @return An object of class `ipw`, the shared return contract of
#'   [propensity::ipw()]. Alongside `estimand`, `ps_mod`, `outcome_mod`, and the
#'   `estimates` table, the result carries two fields describing the variance:
#'
#'   * `se_method`, the string `"mestimation"`, naming how the standard errors
#'     were computed.
#'   * `fit`, the fitted variance system, a list of `theta`, the stacked
#'     parameter vector, and `vcov`, its sandwich covariance. Both are named by
#'     stacked block: `theta_w1` onward for the weight parameters, `beta_`
#'     followed by the design column name for the outcome-model coefficients,
#'     then `mu0`, `mu1`, and one name per effect. The standard errors in
#'     `estimates` are `sqrt(diag(fit$vcov))` read at the effect names.
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
#' ipw(fit, outcome_mod)
#'
#' # The outcome model may also adjust for covariates, in which case the
#' # marginal means are standardized over the estimand's target population.
#' adjusted_mod <- glm(
#'   y ~ exposure + x1,
#'   data = df,
#'   family = binomial(),
#'   weights = .wts
#' )
#'
#' ipw(fit, adjusted_mod)
#'
#' @name ipw.balancing
#' @importFrom causalgenerics ipw
#' @importFrom stats getCall
NULL

# Expose the originating call through the standard model accessor so tools that
# summarize a fit, including the ipw() print method, can label it. Without a
# method the S7 object is not subsettable and the accessor would fail.
getCall_generic <- new_external_generic("stats", "getCall", "x")

method(getCall_generic, balancing) <- function(x, ...) {
  x@call
}

causalgenerics_ipw <- new_external_generic("causalgenerics", "ipw", "ps_mod")

method(causalgenerics_ipw, balancing) <- function(
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
  if (is.null(container@psi_fn) || is.null(container@weights_fn)) {
    abort_ipw_unsupported(reason = "no_hooks")
  }

  estimand <- resolve_ipw_estimand(estimand, ps_mod@estimand)

  exposure_name <- ps_mod@exposure
  weights <- as.numeric(weights(ps_mod))
  validate_ipw_outcome_model(outcome_mod, exposure_name, weights)

  frame <- if (is.null(.data)) stats::model.frame(outcome_mod) else .data
  if (!is.null(.data) && nrow(frame) != ps_mod@n) {
    abort(
      c(
        "{.arg .data} must have one row per observation in the fit.",
        x = "It has {nrow(frame)} row{?s}, but the fit used {ps_mod@n}."
      ),
      error_class = "balancing_ipw_input_error"
    )
  }
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

  # The variance engine composes the sampling weights onto the weights the
  # container's own hook returns, so it takes the fit's sampling weights raw.
  # The preflight above compares against the composed weights instead, because
  # those are the weights the outcome model was fitted with. The focal level
  # goes with them: it names both the population the marginal means standardize
  # over and the group total the reported weights are carried to, and the
  # container records neither.
  variance_system <- ipw_deli_sandwich(
    container = container,
    outcome_mod = outcome_mod,
    frame = frame,
    exposure_name = exposure_name,
    sampling_weights = ps_mod@sampling_weights,
    focal_level = ps_mod@focal_level
  )

  estimates <- ipw_estimates(
    theta = variance_system$theta,
    vcov = variance_system$vcov,
    conf_level = conf_level,
    continuous = is_gaussian_outcome(outcome_mod)
  )

  # The result carries the same fields propensity's own method returns. The
  # stacked parameter vector and its covariance are the whole of the fitted
  # variance system here, so they stand in for the solver object propensity
  # reports: nothing was solved, since every parameter entered at the value its
  # own fit had already found.
  structure(
    list(
      estimand = estimand,
      ps_mod = ps_mod,
      outcome_mod = outcome_mod,
      estimates = estimates,
      se_method = "mestimation",
      fit = variance_system
    ),
    class = "ipw"
  )
}

# The unsupported condition is shared by every configuration ipw() cannot
# handle, so the missing-equations path, the non-binary path, and the
# missing-hooks path all route through one message that names the reason and
# points to the bootstrap workflow.
abort_ipw_unsupported <- function(
  reason = c("no_equations", "exposure_type", "no_hooks"),
  exposure_type = NULL,
  call = rlang::caller_env()
) {
  reason <- rlang::arg_match(reason)
  detail <- switch(
    reason,
    no_equations = c(
      x = "This fit's weights do not solve smooth estimating equations, so the stacked variance is unavailable.",
      i = "Estimating equations come from the estimating-equation family (entropy balancing, inverse probability tilting, just-identified covariate balancing propensity score) with exact balance."
    ),
    exposure_type = c(
      x = "This fit has a {exposure_type} exposure, and only binary exposures are supported.",
      i = "The stacked variance is derived for a binary exposure."
    ),
    no_hooks = c(
      x = "This fit's container does not carry re-evaluation hooks, which the stacked variance differentiates the weight path through.",
      i = "The hooks re-evaluate the estimating functions and the reported weights at new weight parameters."
    )
  )
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

# The outcome model may adjust for covariates, and may interact them with the
# exposure, but it must carry the exposure itself. The marginal means are
# computed by fixing the exposure to each level and predicting, so a model
# without an exposure term has two identical fixed-exposure designs: every
# contrast it reports would be zero, and the table would look like an estimate
# of no effect rather than the absence of an estimator.
#
# The shape checks come first, since a model of the wrong class or the wrong
# form cannot be interrogated for anything else. The three that follow all guard
# against the same failure mode as the exposure check: a model that runs and
# returns an effect table nobody could tell was wrong.
validate_ipw_outcome_model <- function(
  outcome_mod,
  exposure_name,
  expected_weights,
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
  if (!exposure_name %in% term_labels) {
    abort(
      c(
        "{.arg outcome_mod} must include the exposure among its predictors.",
        x = "The exposure {.val {exposure_name}} is not one of its terms.",
        i = "The model may adjust for covariates alongside the exposure."
      ),
      error_class = "balancing_ipw_input_error",
      call = call
    )
  }
  validate_ipw_outcome_family(outcome_mod, call = call)
  validate_ipw_weight_consistency(outcome_mod, expected_weights, call = call)

  # Resolving the response is itself a check: it refuses a factor response the
  # model discarded. The variance engine resolves it again for its own use, but
  # doing it here as well is what attributes the refusal to ipw() rather than to
  # the engine the caller never named.
  resolve_outcome_response(outcome_mod, call = call)
  invisible(NULL)
}

# The effects the method reports are contrasts of two marginal means, and each
# contrast is derived for a particular reading of those means. A binary outcome
# gives probabilities, whose contrasts are the risk difference, the log risk
# ratio, and the log odds ratio; a gaussian outcome gives conditional means,
# whose contrast is their difference. A family outside that set has marginal
# means neither reading describes. A count outcome is the clearest case: its
# means are rates, so the odds ratio is undefined and its row would be reported
# as a missing value beside a risk difference the label does not fit.
#
# The quasibinomial family belongs with the binomial one. It shares the binomial
# variance function, so the marginal means are probabilities and every contrast
# reads the same; only the dispersion differs, and the dispersion does not enter
# the sandwich. A plain lm reports a gaussian family through the same accessor,
# so it needs no separate branch.
validate_ipw_outcome_family <- function(
  outcome_mod,
  call = rlang::caller_env()
) {
  supported <- c("binomial", "quasibinomial", "gaussian")
  family <- stats::family(outcome_mod)$family
  if (family %in% supported) {
    return(invisible(NULL))
  }
  abort(
    c(
      "{.arg outcome_mod} must come from a supported outcome family.",
      x = "Its family is {.val {family}}.",
      i = "The supported families are {.val {supported}}.",
      i = "See the inference vignette for a bootstrap workflow with other families."
    ),
    error_class = "balancing_ipw_input_error",
    call = call
  )
}

# The stacked variance differentiates the outcome-model score through the
# weights the fit produced, so it describes the system that was actually solved
# only when the outcome model was fitted at those weights. Fitted at any other
# weights, or at none, the model's coefficients and its variance describe two
# different estimators, and the result is still an ordinary-looking effect
# table. The weights compared against are the fit's composed weights, which
# already carry the sampling weights, since those are the weights the outcome
# model is meant to have been fitted with.
#
# The comparison is per unit and relative: every unit's model weight must agree
# with the weight the fit reports for it to within 1e-6 of that weight's own
# magnitude. A mean relative difference would let one badly wrong unit hide
# behind the rest, and an absolute difference would not travel across weight
# scales, since a fit reporting weights that sum to the sample size and one
# reporting weights that average one differ by a factor of that size. Weights
# below one are compared against one, which makes the tolerance absolute in that
# range rather than demanding relative agreement near zero that floating point
# arithmetic cannot deliver.
#
# A model fitted without weights records none at all rather than a vector of
# ones, so a missing vector is read as ones. That is what such a model actually
# fitted, and reading it that way is what makes an unweighted model on a
# weighted fit an error rather than a case that quietly skips the check.
validate_ipw_weight_consistency <- function(
  outcome_mod,
  expected_weights,
  call = rlang::caller_env()
) {
  model_weights <- stats::weights(outcome_mod)
  if (is.null(model_weights)) {
    model_weights <- rep(1, length(expected_weights))
  }
  model_weights <- as.numeric(model_weights)

  # A weight vector of a different length cannot belong to this fit, so it is
  # reported as a mismatch rather than compared through a recycled subtraction.
  magnitude <- pmax(abs(expected_weights), 1)
  deviation <- if (length(model_weights) == length(expected_weights)) {
    max(abs(model_weights - expected_weights) / magnitude)
  } else {
    Inf
  }
  if (deviation <= 1e-6) {
    return(invisible(NULL))
  }
  abort(
    c(
      "{.arg outcome_mod} must be fitted with the weights from {.arg ps_mod}.",
      x = "Its weights differ from the fit's, compared per unit at relative tolerance 1e-6.",
      i = "Refit it with {.code weights = weights(fit)}, where {.arg fit} is the balancing fit."
    ),
    error_class = "balancing_ipw_input_error",
    call = call
  )
}

# The sandwich needs the response on the scale the outcome model actually
# modeled, which is not always the scale the model frame stores. A binomial glm
# fitted on a two-level factor models zero for the first level and one for the
# other, while the model frame keeps the factor, whose integer codes are one and
# two. Coercing the model-frame response would put such an outcome on the wrong
# scale and corrupt every standard error, silently, because the point estimates
# read off the coefficients and would not move. A glm records the scale it
# modeled in `$y`, so that is the reliable source for every family.
#
# A glm fitted with `y = FALSE` keeps no stored response. A numeric model-frame
# response is unambiguous and is used directly, but a factor one leaves the
# modeled scale a guess, so it is refused instead. An lm does not store `$y` by
# default and its model-frame response is already numeric, so it takes the same
# direct path.
resolve_outcome_response <- function(outcome_mod, call = rlang::caller_env()) {
  if (inherits(outcome_mod, "glm") && !is.null(outcome_mod$y)) {
    return(as.numeric(outcome_mod$y))
  }
  response <- stats::model.response(stats::model.frame(outcome_mod))
  if (is.factor(response)) {
    abort(
      c(
        "{.arg outcome_mod} must carry the response it modeled.",
        x = "It has a factor response but was fitted with {.code y = FALSE}, which discards that response.",
        i = "Refit it with {.code y = TRUE}, or on a numeric 0/1 response."
      ),
      error_class = "balancing_ipw_input_error",
      call = call
    )
  }
  as.numeric(response)
}

# The outcome model's terms with the response and any offset removed, so that
# they describe the design alone. The offset has to go for two reasons. It is
# not a column of the design, and the fitted model already carries its per-unit
# value, so re-deriving it is redundant. And an offset written into the formula
# cannot be re-derived from the model frame at all: the frame stores it under
# the deparsed call, `offset(v)`, while the terms ask for the variable `v`,
# which is not a column of the frame. Rebuilding the terms from the term labels,
# which never include an offset, keeps a formula offset and an `offset` argument
# on the same path.
design_terms <- function(outcome_mod) {
  terms <- stats::delete.response(stats::terms(outcome_mod))
  if (is.null(attr(terms, "offset"))) {
    return(terms)
  }
  stats::terms(
    stats::reformulate(
      attr(terms, "term.labels"),
      intercept = attr(terms, "intercept") == 1L,
      env = environment(terms)
    )
  )
}

# The outcome design and predictions with the exposure fixed to one level, for
# the marginal-mean equations. The model's terms, contrasts, and factor levels
# are reused so the columns line up with the fitted coefficients.
#
# An offset is part of the linear predictor rather than of the design, so a
# model that carries one supplies it separately and it is added to eta. Fixing
# the exposure does not change it, since it is a known per-unit quantity.
fixed_exposure_pieces <- function(
  outcome_mod,
  frame,
  exposure_name,
  level,
  offset = NULL
) {
  frame[[exposure_name]] <- level
  terms <- design_terms(outcome_mod)
  model_frame <- stats::model.frame(terms, frame, xlev = outcome_mod$xlevels)
  design <- stats::model.matrix(
    terms,
    model_frame,
    contrasts.arg = outcome_mod$contrasts
  )
  eta <- as.numeric(design %*% stats::coef(outcome_mod))
  if (!is.null(offset)) {
    eta <- eta + offset
  }
  family <- stats::family(outcome_mod)
  list(
    design = design,
    eta = eta,
    mu = family$linkinv(eta),
    mu_eta = family$mu.eta(eta)
  )
}

# The effect rows. Each effect is a parameter of the stacked system, so its
# estimate is that parameter's entry in the stacked parameter vector and its
# standard error is the square root of the matching diagonal entry of the
# covariance. Nothing is contrasted or differentiated here, which is what keeps
# the contrast formulas stated once, in the stack itself.
ipw_estimates <- function(theta, vcov, conf_level, continuous) {
  effects <- ipw_contrast_names(continuous)
  estimate <- unname(theta[effects])
  std_err <- unname(sqrt(diag(vcov)[effects]))
  z <- estimate / std_err
  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)

  data.frame(
    effect = effects,
    estimate = estimate,
    std.err = std_err,
    z = z,
    ci.lower = estimate - z_value * std_err,
    ci.upper = estimate + z_value * std_err,
    conf.level = conf_level,
    p.value = 2 * (1 - stats::pnorm(abs(z)))
  )
}
