# balancing registers a method on the shared ipw() generic so that a fitted
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
#' A [balancing] fit registers a method on [causalgenerics::ipw()], so a set of
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
#' A categorical exposure reports those same measures for each non-reference
#' level against the reference level, which is the first of the fit's own levels:
#' a factor's declared level order for a factor exposure, and the sorted values
#' otherwise. A K-level exposure therefore contributes K marginal means and one
#' block of measures per non-reference level, and the estimates table gains a
#' `comparison` column, placed after `effect`, naming each contrast as
#' `"<level> vs <reference>"`. A binary exposure keeps the table it has always
#' returned, with no `comparison` column.
#'
#' A continuous exposure has no levels to contrast, so there is no pair of
#' marginal means to difference. What the method reports instead is the
#' dose-response coefficient of a weighted marginal structural model: the
#' balancing weights break the exposure-covariate association, the outcome model
#' carries exactly one term in the exposure, and that term's coefficient is the
#' effect of a one-unit change in the exposure on the model's own link scale. The
#' estimates table holds a single row, keeping the columns it holds for every
#' other exposure and gaining no `comparison` column, and that row is named for
#' the link: `slope` for an identity link, whether the model arrives as a
#' [stats::lm()] or as a gaussian [stats::glm()]; `log(or)` for a logit; and
#' `log(rr)` for a log link. Another link raises `balancing_ipw_input_error`,
#' since its coefficient is none of those three. A continuous fit targets the
#' average treatment effect and nothing else, so an `estimand` supplied
#' alongside it either agrees or raises `balancing_estimand_error`, as it does
#' for any other fit.
#'
#' The standard errors come from a stacked M-estimator that the deli package
#' differentiates and sandwiches. The stacked parameter vector holds four
#' blocks: the
#' weight parameters, the outcome-model coefficients, the marginal means, one per
#' exposure level, and the effect contrasts. Their estimating functions are, in
#' the same order,
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
#' A continuous exposure stacks the same system with the g-computation half
#' removed. There are no fixed-exposure predictions to standardize and no
#' contrasts to form, so the stack holds the weight parameters and the marginal
#' structural model's coefficients alone, and the effect is already one of those
#' coefficients rather than a parameter derived from them.
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
#' That covariance is a large-sample one, and how large a sample it takes
#' differs by exposure. The binary risk-difference standard error is calibrated
#' at a few hundred observations; the continuous slope's is anticonservative
#' there. Over 500 draws its ratio of mean standard error to the standard
#' deviation of the estimates was 0.836 at 300 observations and 0.942 at 1200.
#' Read a continuous interval as the asymptotic statement it is, and prefer the
#' bootstrap of the inference vignette at small sample sizes.
#'
#' The method is available only for fits whose weights solve smooth estimating
#' equations: the estimating-equation family (entropy balancing, inverse
#' probability tilting, and the just-identified covariate balancing propensity
#' score) with exact balance at a binary or categorical exposure, and entropy
#' balancing with exact balance at a continuous one. Any other fit, including an
#' entropy fit at a positive tolerance, an over-identified or quadratic-program
#' fit, or a continuous fit from a method that solves no estimating equations,
#' raises `balancing_ipw_unsupported_error` and points to the bootstrap workflow
#' described in the inference vignette.
#'
#' The outcome model must carry the same exposure levels the fit weighted.
#' Fitting it on data that drop a level, or that carry one the fit never saw,
#' raises `balancing_ipw_input_error`: the counterfactual predictions and the
#' weights would then describe different exposures, and the resulting table would
#' name a contrast it had not computed.
#'
#' The outcome model must carry the exposure among its predictors, and may
#' adjust for covariates alongside it, including in interactions with the
#' exposure: `y ~ exposure`, `y ~ exposure + x1 + x2`, and `y ~ exposure * x1`
#' are all supported. A categorical exposure enters as a factor, so its
#' fixed-exposure designs come from the model's own contrasts. A model without an
#' exposure term raises `balancing_ipw_input_error`, since its fixed-exposure
#' predictions would all be the same prediction and every contrast it reported
#' would be zero.
#'
#' A continuous exposure narrows that to exactly one term, the exposure itself,
#' since what it reports is one coefficient of the model rather than a contrast
#' of predictions from it. `y ~ exposure` and `y ~ exposure + x1` are supported,
#' while `y ~ exposure + I(exposure^2)`, `y ~ poly(exposure, 2)`, and
#' `y ~ exposure * x1` raise `balancing_ipw_input_error`: each carries a second
#' design column in the exposure, so the effect of a one-unit change depends on
#' where it is read and no single coefficient is it. The check reads the model's
#' terms rather than the text of its formula, so a transformed or interacted
#' exposure term is caught however it is written.
#'
#' An offset reaches the linear predictor without being a term, so it is checked
#' on its own, and for every exposure type. An offset expression naming the
#' exposure, written either into the formula or passed through the model's
#' `offset` argument, raises `balancing_ipw_input_error`. An offset is held at
#' its observed value while the exposure is fixed to each level, so a discrete
#' exposure's marginal means would read one exposure in the design and another in
#' the offset, and a continuous exposure's coefficient would be something other
#' than the effect of a one-unit change. An exposure-free offset stays supported.
#' That check is static, so it accepts any offset whose stored expression does
#' not name the exposure: a precomputed vector under some other symbol, and
#' equally a wrapper that forwards the offset through its dots, which records
#' `..1` in the fitted call. Keeping the exposure out of such an offset is the
#' caller's to honor, and it matters most for a discrete exposure, where a
#' laundered offset corrupts the fixed-exposure marginal means and can reverse
#' the sign of the reported contrast rather than merely shifting a coefficient.
#'
#' The exposure may be a factor, a character column, or an integer code, and a
#' term that transforms it counts as carrying it. A character column becomes a
#' factor in a model formula on its own, while an integer code is written
#' `factor(exposure)` for a model with one parameter per level; left
#' untransformed, an integer code enters as a slope in the codes and the marginal
#' means are the g-computation means of that model rather than of a saturated
#' one. A transformed exposure is not a column of the model frame, since the
#' frame stores the transformation, so such a model is passed with `.data`.
#'
#' Which population the marginal means are averaged over is part of the
#' estimand, and matters as soon as the outcome model adjusts for anything. A
#' marginal model is saturated in the exposure, one free parameter per exposure
#' level, so absent an offset it predicts a single value per level, its
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
#' Of the two binomial families, [stats::quasibinomial()] is the one to fit a
#' binary outcome with here, and the examples below use it: balancing weights are
#' not counts, so [stats::binomial()] warns about non-integer successes at every
#' fit. The two solve the same estimating equation, since they share the binomial
#' variance function, and the dispersion the quasi family estimates never reaches
#' the sandwich, which is built from the score alone. The results are therefore
#' identical rather than merely close.
#'
#' An offset is supported, written either as an `offset()` term in the outcome
#' formula or passed through the model's `offset` argument, so long as it does
#' not read the exposure. It is carried through both the outcome-model score and
#' the fixed-exposure linear predictors, so the marginal means are the
#' g-computation means with each unit's offset held at its observed value. That
#' is the right treatment of a quantity the exposure does not move and the wrong
#' treatment of one it does, which is why the exposure-reading case is refused
#' above.
#'
#' @references
#' Kostouraki A, Hajage D, Rachet B, et al. On variance estimation of the
#' inverse probability-of-treatment weighting estimator: A tutorial for
#' different types of propensity score weights. *Statistics in Medicine*.
#' 2024;43(13):2672-2694. \doi{10.1002/sim.10078}
#'
#' @param wt_mod A [balancing] fit that produced the weights.
#' @param outcome_mod A weighted outcome model of class [stats::glm()] or
#'   [stats::lm()], fitted with the balancing weights and carrying the exposure
#'   among its predictors. It may adjust for covariates alongside the exposure.
#'   For a continuous exposure it is a marginal structural model carrying
#'   exactly one term in the exposure.
#' @param .data The data frame holding the exposure and outcome. If `NULL`, the
#'   values are taken from the outcome model frame. It carries the exposure
#'   column the fixed-exposure predictions are built from, so it has nothing to
#'   supply for a continuous exposure, which makes no such predictions, and is
#'   ignored there.
#' @param estimand The causal estimand. If `NULL`, the fit's estimand is used.
#'   As in [balance()], `"atc"` is accepted as a synonym for `"atu"`. Supplying
#'   an estimand that disagrees with the fit raises `balancing_estimand_error`.
#' @param conf_level The confidence level for the intervals. Default `0.95`.
#' @param ... Ignored, for compatibility with the generic.
#'
#' @return An object of class `ipw`, an implementation of
#'   [causalgenerics::ipw()]. Alongside `estimand`, `wt_mod`, `outcome_mod`, and
#'   the `estimates` table, the result carries two fields describing the
#'   variance:
#'
#'   * `se_method`, the string `"mestimation"`, naming how the standard errors
#'     were computed.
#'   * `fit`, the fitted variance system, a list of `theta`, the stacked
#'     parameter vector, and `vcov`, its sandwich covariance. Both are named by
#'     stacked block: `theta_w1` onward for the weight parameters, `beta_`
#'     followed by the design column name for the outcome-model coefficients,
#'     then the marginal means and one name per contrast. A binary exposure names
#'     its means `mu0` and `mu1` and its contrasts by measure alone; a
#'     categorical exposure names each mean `mu_` followed by its level and each
#'     contrast by measure and level, as `rd_b`. A continuous exposure carries
#'     neither block, since the effect is one of the outcome-model coefficients;
#'     that coefficient is named for the effect, as `slope`, in place of the
#'     `beta_` name the others carry. The standard errors in `estimates` are
#'     `sqrt(diag(fit$vcov))` read at those effect names.
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
#'
#' # quasibinomial() solves the same estimating equation as binomial() and does
#' # not warn that weights are not counts.
#' outcome_mod <- glm(
#'   y ~ exposure,
#'   data = df,
#'   family = quasibinomial(),
#'   weights = .wts
#' )
#'
#' ipw(fit, outcome_mod)
#'
#' # The outcome model may also adjust for covariates, in which case the
#' # marginal means are standardized over the estimand's target population.
#' adjusted_mod <- glm(
#'   y ~ exposure + x1,
#'   data = df,
#'   family = quasibinomial(),
#'   weights = .wts
#' )
#'
#' ipw(fit, adjusted_mod)
#'
#' # A categorical exposure reports each level against the reference level, and
#' # the estimates table names the comparison.
#' odds_b <- exp(0.6 * x1)
#' odds_c <- exp(-0.5 * x1)
#' denominator <- 1 + odds_b + odds_c
#' draw <- runif(n)
#' df$arm <- factor(
#'   ifelse(
#'     draw < 1 / denominator,
#'     "a",
#'     ifelse(draw < (1 + odds_b) / denominator, "b", "c")
#'   ),
#'   levels = c("a", "b", "c")
#' )
#' df$relapse <- rbinom(
#'   n,
#'   1,
#'   plogis(-0.4 + 0.5 * (df$arm == "b") + 0.9 * (df$arm == "c") + 0.3 * x1)
#' )
#'
#' arm_fit <- balance(df, arm, x1, method = bw_ipt(), estimand = "ate")
#' df$.arm_wts <- as.numeric(weights(arm_fit))
#' arm_mod <- glm(
#'   relapse ~ arm,
#'   data = df,
#'   family = quasibinomial(),
#'   weights = .arm_wts
#' )
#'
#' ipw(arm_fit, arm_mod)
#'
#' # A continuous exposure reports one effect, the exposure coefficient of a
#' # weighted marginal structural model, named for that model's link.
#' df$dose <- 0.7 * x1 + rnorm(n)
#' df$score <- 2 + 0.5 * df$dose + 0.4 * x1 + rnorm(n)
#'
#' dose_fit <- balance(df, dose, x1, method = bw_entropy(), estimand = "ate")
#' df$.dose_wts <- as.numeric(weights(dose_fit))
#' dose_mod <- lm(score ~ dose, data = df, weights = .dose_wts)
#'
#' ipw(dose_fit, dose_mod)
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

causalgenerics_ipw <- new_external_generic("causalgenerics", "ipw", "wt_mod")

method(causalgenerics_ipw, balancing) <- function(
  wt_mod,
  outcome_mod,
  .data = NULL,
  estimand = NULL,
  conf_level = 0.95,
  ...
) {
  container <- wt_mod@estimating_equations
  if (is.null(container)) {
    abort_ipw_unsupported(reason = "no_equations")
  }
  if (is.null(container@psi_fn) || is.null(container@weights_fn)) {
    abort_ipw_unsupported(reason = "no_hooks")
  }
  # Which exposure the fit weighted decides which stacked system is assembled,
  # not whether one can be. A method that solves smooth estimating equations at
  # a continuous fit carries the same container with the same hooks, so the two
  # refusals above are the whole gate: a continuous fit from a method that
  # solves none, or an entropy fit at a positive tolerance, arrives here without
  # a container and is turned away for that reason rather than for its exposure.
  categorical <- identical(wt_mod@exposure_type, "categorical")
  continuous_exposure <- identical(wt_mod@exposure_type, "continuous")

  estimand <- resolve_ipw_estimand(estimand, wt_mod@estimand)

  exposure_name <- wt_mod@exposure
  weights <- as.numeric(weights(wt_mod))
  validate_ipw_outcome_model(
    outcome_mod,
    exposure_name,
    weights,
    continuous_exposure = continuous_exposure
  )

  # A continuous exposure has no counterfactual designs to build, so it needs
  # neither the model frame nor the exposure levels: its whole stack is the
  # weight parameters and the marginal structural model's coefficients, and the
  # effect it reports is one of those coefficients rather than a contrast of
  # marginal means.
  if (continuous_exposure) {
    variance_system <- ipw_deli_msm_sandwich(
      container = container,
      outcome_mod = outcome_mod,
      exposure_name = exposure_name,
      sampling_weights = wt_mod@sampling_weights,
      call = rlang::current_env()
    )
    effect <- msm_effect_name(outcome_mod)
    estimates <- ipw_estimate_rows(
      theta = variance_system$theta,
      vcov = variance_system$vcov,
      conf_level = conf_level,
      keys = effect,
      effects = effect
    )
  } else {
    frame <- resolve_ipw_frame(outcome_mod, .data, exposure_name, wt_mod@n)

    # The exposure levels come off the fit rather than being sorted out of the
    # frame again. The fit recorded the ordering its solve used, whose first
    # element is the reference level every contrast below is measured against.
    # Sorting the data's own values here would agree with that by coincidence and
    # disagree silently whenever a factor declares its levels out of alphabetical
    # order, which would report every contrast against the wrong level.
    levels <- wt_mod@exposure_levels
    validate_ipw_exposure_levels(frame[[exposure_name]], levels, exposure_name)

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
      levels = levels,
      categorical = categorical,
      sampling_weights = wt_mod@sampling_weights,
      focal_level = wt_mod@focal_level,
      call = rlang::current_env()
    )

    estimates <- ipw_estimates(
      theta = variance_system$theta,
      vcov = variance_system$vcov,
      conf_level = conf_level,
      continuous = is_gaussian_outcome(outcome_mod),
      levels = if (categorical) levels else NULL
    )
  }

  # The result is built by the shared constructor, so it carries the fields the
  # common print and as.data.frame methods read. The stacked parameter vector
  # and its covariance are the whole of the fitted variance system here, so they
  # stand in for the solver object a propensity score fit reports: nothing was
  # solved, since every parameter entered at the value its own fit had already
  # found.
  causalgenerics::new_ipw(
    estimand = estimand,
    wt_mod = wt_mod,
    outcome_mod = outcome_mod,
    estimates = estimates,
    se_method = "mestimation",
    fit = variance_system
  )
}

# The data the counterfactual designs of a discrete exposure are built from,
# which is the outcome model's own frame unless the caller supplied one. A frame
# of the wrong length cannot belong to this fit, and a frame without the exposure
# column cannot have the exposure fixed in it, which is what a formula that
# transforms the exposure leaves behind: the frame stores the transformation
# rather than the column, so such a model needs `.data`.
resolve_ipw_frame <- function(
  outcome_mod,
  .data,
  exposure_name,
  n,
  call = rlang::caller_env()
) {
  frame <- if (is.null(.data)) stats::model.frame(outcome_mod) else .data
  if (!is.null(.data) && nrow(frame) != n) {
    abort(
      c(
        "{.arg .data} must have one row per observation in the fit.",
        x = "It has {nrow(frame)} row{?s}, but the fit used {n}."
      ),
      error_class = "balancing_ipw_input_error",
      call = call,
      .envir = environment()
    )
  }
  if (!exposure_name %in% names(frame)) {
    abort(
      c(
        "The exposure {.val {exposure_name}} is not in the outcome model frame.",
        i = "Supply {.arg .data} containing the exposure column when the outcome formula transforms it."
      ),
      error_class = "balancing_ipw_input_error",
      call = call,
      .envir = environment()
    )
  }
  frame
}

# The unsupported condition is shared by every configuration ipw() cannot
# handle, so the missing-equations path and the missing-hooks path both route
# through one message that names the reason and points to the bootstrap
# workflow.
abort_ipw_unsupported <- function(
  reason = c("no_equations", "no_hooks"),
  call = rlang::caller_env()
) {
  reason <- rlang::arg_match(reason)
  detail <- switch(
    reason,
    no_equations = c(
      x = "This fit's weights do not solve smooth estimating equations, so the stacked variance is unavailable.",
      i = "Estimating equations come from the estimating-equation family (entropy balancing, inverse probability tilting, just-identified covariate balancing propensity score) with exact balance."
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

# The outcome model describes the same exposure the fit weighted only when it
# carries the same levels. A level the fit weighted but the data no longer
# contain has no counterfactual design to predict, and a level the data carry but
# the fit never saw was never balanced and has no weights behind it. Either way
# the effect table would look ordinary while describing a different contrast
# than the one it names, so the mismatch is refused with both sides named.
validate_ipw_exposure_levels <- function(
  exposure,
  levels,
  exposure_name,
  call = rlang::caller_env()
) {
  observed <- unique(as.character(exposure[!is.na(exposure)]))
  absent <- setdiff(levels, observed)
  unexpected <- setdiff(observed, levels)
  if (length(absent) == 0 && length(unexpected) == 0) {
    return(invisible(NULL))
  }

  detail <- character(0)
  if (length(absent) > 0) {
    detail <- c(
      detail,
      x = "The fit weighted {.val {absent}}, which the data do not contain."
    )
  }
  if (length(unexpected) > 0) {
    detail <- c(
      detail,
      x = "The data contain {.val {unexpected}}, which the fit did not weight."
    )
  }
  abort(
    c(
      "The exposure {.val {exposure_name}} must carry the same levels in the outcome model as in the fit.",
      detail,
      i = "Fit the outcome model on the data the weights were fitted from."
    ),
    error_class = "balancing_ipw_input_error",
    call = call,
    .envir = environment()
  )
}

# The estimand either comes from the fit or must agree with it. A supplied
# estimand that contradicts the fit is a specification error naming the knob.
# Agreement is judged on the canonical spelling, since the fit stores one name
# for the untreated target and accepts two, and the request is echoed back as the
# caller spelled it.
#
# The vocabulary is matched first, against the same choices `balance()` matches
# against. A name no estimand carries is a different failure from one the fit
# does not target, and reporting a misspelling as a disagreement with the fit
# names the fit's estimand while leaving the caller to notice that theirs is not
# an estimand at all.
resolve_ipw_estimand <- function(
  estimand,
  fit_estimand,
  call = rlang::caller_env()
) {
  if (is.null(estimand)) {
    return(fit_estimand)
  }
  estimand <- rlang::arg_match0(
    estimand,
    estimand_choices(),
    arg_nm = "estimand",
    error_call = call
  )
  requested <- canonical_estimand(estimand)
  if (!identical(requested, fit_estimand)) {
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
  requested
}

is_gaussian_outcome <- function(outcome_mod) {
  if (inherits(outcome_mod, "glm")) {
    return(identical(stats::family(outcome_mod)$family, "gaussian"))
  }
  # A plain lm is a linear model.
  TRUE
}

# The variables the outcome model's terms are built from, which is what the
# exposure is looked for among. The exposure counts as present whenever a term
# reads it, not only when a term is spelled exactly like it: an exposure stored
# as a character column or as an integer code is written `factor(exposure)` to
# give the outcome model a level per group, and that model carries the exposure
# as surely as one whose formula names the column outright. The fixed-exposure
# designs are rebuilt from the model's own terms with the exposure column set to
# each level, so a transformation is applied again at each of them; what the
# transformation costs is the model frame, which stores the transformed column
# rather than the exposure, so such a model needs `.data`.
model_term_variables <- function(outcome_mod) {
  labels <- attr(stats::terms(outcome_mod), "term.labels")
  unique(unlist(lapply(labels, function(label) all.vars(str2lang(label)))))
}

# The outcome model may adjust for covariates, and may interact them with the
# exposure, but it must carry the exposure itself. The marginal means are
# computed by fixing the exposure to each level and predicting, so a model
# without an exposure term has two identical fixed-exposure designs: every
# contrast it reports would be zero, and the table would look like an estimate
# of no effect rather than the absence of an estimator.
#
# The shape checks come first, since a model of the wrong class or the wrong
# form cannot be interrogated for anything else. Those that follow all guard
# against the same failure mode as the exposure check: a model that runs and
# returns an effect table nobody could tell was wrong. Some of them apply to a
# continuous exposure alone, whose reported effect is a coefficient of this model
# rather than a contrast of predictions from it.
validate_ipw_outcome_model <- function(
  outcome_mod,
  exposure_name,
  expected_weights,
  continuous_exposure = FALSE,
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
  validate_ipw_response_shape(outcome_mod, call = call)
  if (!exposure_name %in% model_term_variables(outcome_mod)) {
    # What the model may do with the exposure depends on the exposure type. A
    # discrete exposure is read through predictions with the exposure fixed to
    # each level, so a transformation of it is fine and `factor()` is the usual
    # one. A continuous exposure reports the exposure's own coefficient, and the
    # slope check below refuses every transformation, so offering one here would
    # send that caller to a second refusal.
    latitude <- if (continuous_exposure) {
      "The model may adjust for covariates alongside the exposure, which must enter as a term of its own."
    } else {
      "The model may adjust for covariates alongside the exposure, and may carry the exposure inside a transformation such as {.fun factor}."
    }
    # A model whose only mention of the exposure is an offset reaches here, since
    # an offset is not a term, and would otherwise be told the exposure is absent
    # while the caller can see it written in the formula. The offset check below
    # would refuse such a model anyway once a real exposure term were added, so
    # the pointer saves a round trip as well as the confusion.
    offset_note <- if (
      any(offsets_read_exposure(offset_expressions(outcome_mod), exposure_name))
    ) {
      "An offset is not a term, so an exposure that reaches the model only through one is not carried by it."
    }
    abort(
      c(
        "{.arg outcome_mod} must include the exposure among its predictors.",
        x = "The exposure {.val {exposure_name}} appears in none of its terms.",
        i = latitude,
        i = offset_note
      ),
      error_class = "balancing_ipw_input_error",
      call = call
    )
  }
  validate_ipw_outcome_family(outcome_mod, call = call)
  # An offset carrying the exposure spoils both exposure types, so the check runs
  # before the branch rather than inside it.
  validate_ipw_exposure_offset(outcome_mod, exposure_name, call = call)
  if (continuous_exposure) {
    validate_ipw_exposure_slope(outcome_mod, exposure_name, call = call)
    # The link names the reported effect, so a link no effect name describes is
    # refused here, where the frame the message describes is still the ipw()
    # call the caller made rather than the variance engine they never named.
    msm_effect_name(outcome_mod, call = call)
  }
  validate_ipw_weight_consistency(outcome_mod, expected_weights, call = call)

  # Resolving the response is itself a check: it refuses a factor response the
  # model discarded. The variance engine resolves it again for its own use, but
  # doing it here as well is what attributes the refusal to ipw() rather than to
  # the engine the caller never named.
  resolve_outcome_response(outcome_mod, call = call)
  invisible(NULL)
}

# A response of more than one column is two different models, and the stack
# handles neither. Through `glm()` it is the grouped binomial form: each row
# carries a count of successes and a count of failures rather than one Bernoulli
# draw, and the model's prior weights are the weights it was given times each
# row's trial count, a scale the fit knows nothing about. Through `lm()` it is a
# multivariate fit, whose weights are the weights it was given but whose
# coefficients are one block per response column and whose score is therefore not
# a single per-unit vector. Either way the stacked variance rebuilds an
# outcome-model score from the weights the fit reports at each set of weight
# parameters, and that score is not the one the model fitted.
#
# The shape is therefore refused, and refused here rather than left to the weight
# preflight, which sees the grouped binomial form's scaled prior weights and
# would report a mismatch to a caller who supplied exactly the fit's weights, and
# sees nothing at all wrong with the multivariate one.
validate_ipw_response_shape <- function(
  outcome_mod,
  call = rlang::caller_env()
) {
  response <- stats::model.response(stats::model.frame(outcome_mod))
  columns <- NCOL(response)
  if (columns <= 1L) {
    return(invisible(NULL))
  }
  abort(
    c(
      "{.fun ipw} cannot compute a stacked variance for a multi-column response.",
      x = "Its response is a matrix of {columns} columns, which is the grouped binomial form for {.fun glm} and a multivariate fit for {.fun lm}.",
      i = "A grouped binomial fit scales the weights it was given by each row's trial count; a multivariate fit carries one coefficient block per response column. Neither leaves a single per-unit score the stack can rebuild.",
      i = "Fit the weights and the outcome model on data with one row per unit and one response column, or see the inference vignette for a bootstrap workflow."
    ),
    error_class = "balancing_ipw_unsupported_error",
    call = call,
    .envir = environment()
  )
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

# A continuous exposure reports one coefficient of the outcome model, so the
# model has to have one coefficient to report: exactly one term reading the
# exposure, and that term the exposure itself. A second term in the exposure, a
# term that expands to several columns, or an interaction with a covariate all
# describe a dose-response surface whose effect depends on where it is read, and
# no single coefficient of such a model is the effect the table would name.
#
# The check reads the model's terms rather than its formula text, so an exposure
# that arrives transformed or interacted is caught however it is written:
# `y ~ a + I(a^2)`, `y ~ poly(a, 2)`, and `y ~ a * x1` are refused while `y ~ a`
# and `y ~ a + v` are accepted. Covariates that do not interact with the exposure
# stay free, since they leave the exposure's own column and its coefficient in
# place. An offset is not a term and so is not this check's to see; the check
# below reads it.
validate_ipw_exposure_slope <- function(
  outcome_mod,
  exposure_name,
  call = rlang::caller_env()
) {
  labels <- attr(stats::terms(outcome_mod), "term.labels")
  reads_exposure <- vapply(
    labels,
    function(label) exposure_name %in% all.vars(str2lang(label)),
    logical(1),
    USE.NAMES = FALSE
  )
  exposure_terms <- labels[reads_exposure]

  # The coefficient the term contributes has to be there under the exposure's own
  # name, which is what the reported effect is read from. A lone exposure term
  # entering as anything but a numeric column would carry the right label and the
  # wrong design. Both comparisons run against the quoted spelling, since a term
  # label and a coefficient name carry the back-quotes a non-syntactic column
  # needs while the fit records the bare name.
  quoted <- quoted_name(exposure_name)
  linear <- identical(exposure_terms, quoted) &&
    quoted %in% names(stats::coef(outcome_mod))
  if (!linear) {
    abort(
      c(
        "{.arg outcome_mod} must carry exactly one coefficient in the exposure.",
        x = "The exposure {.val {exposure_name}} enters through {.val {exposure_terms}}.",
        i = "A continuous exposure reports the coefficient of {.code {quoted}} itself, which a second exposure term, a transformation, or an interaction leaves undefined.",
        i = "The model may adjust for covariates that do not interact with the exposure."
      ),
      error_class = "balancing_ipw_input_error",
      call = call,
      .envir = environment()
    )
  }
  invisible(NULL)
}

# Both places an offset can enter a fitted model, since neither records the
# other: a formula offset lives in the terms object, whose `offset` attribute
# indexes it among the variables, and the `offset` argument lives in the fitted
# call. A model carrying no offset yields the one `NULL` the argument left
# behind, which names nothing.
offset_expressions <- function(outcome_mod) {
  terms <- stats::terms(outcome_mod)
  variables <- attr(terms, "variables")
  c(
    lapply(
      attr(terms, "offset"),
      function(position) variables[[position + 1L]]
    ),
    list(stats::getCall(outcome_mod)$offset)
  )
}

# Which of a model's offset expressions name the exposure. Two checks read this:
# the refusal below, which reports the offending expressions, and the
# exposure-presence check, which points a caller at an offset when that is the
# only place the exposure is written.
offsets_read_exposure <- function(expressions, exposure_name) {
  vapply(
    expressions,
    function(expression) exposure_name %in% all.vars(expression),
    logical(1)
  )
}

# An offset is the second way the exposure can reach the linear predictor, and
# it reaches it past the term labels the check above reads: an offset is not a
# term, so a model whose offset is the exposure carries exactly one exposure
# term and is accepted by that check while what it reports is no longer the
# effect of the exposure. `y ~ a + offset(a)` estimates a coefficient one below
# the slope it would report, and a table naming that number the slope would be
# wrong in a way nothing about it shows.
#
# A discrete exposure is spoiled the same way through a different route. Its
# marginal means come from predictions with the exposure fixed to each level,
# and an offset is held at its observed value across them, since an offset is a
# known per-unit quantity the counterfactual does not move. An offset computed
# from the exposure is not known that way, so each prediction reads one exposure
# in the design and another in the offset, and the means are neither factual nor
# counterfactual. On the package's own binary fixture the contrast of those
# means comes back with the opposite sign to the g-computation the same model
# implies, which is why the check runs for every exposure type.
#
# The inspection is static, so what it refuses is an offset expression that
# names the exposure. An offset that does not, including a precomputed vector
# whose symbol carries some other name, is beyond reach: nothing in the fitted
# model distinguishes such a vector from any other per-unit quantity. The
# one-term contract stands on the documentation there.
validate_ipw_exposure_offset <- function(
  outcome_mod,
  exposure_name,
  call = rlang::caller_env()
) {
  expressions <- offset_expressions(outcome_mod)
  reads_exposure <- offsets_read_exposure(expressions, exposure_name)
  if (!any(reads_exposure)) {
    return(invisible(NULL))
  }
  # The offending expressions are spelled the way the model writes them, which
  # back-quotes a name R cannot parse as a symbol. An offset supplied through the
  # argument is often a bare column name, so this is the spelling the slope
  # message already uses for the exposure and the two agree.
  offending <- vapply(
    expressions[reads_exposure],
    function(expression) deparse1(expression, backtick = TRUE),
    character(1)
  )
  abort(
    c(
      "{.arg outcome_mod} must not carry an offset that reads the exposure.",
      x = "{cli::qty(offending)}The offset{?s} {.val {offending}} {?reads/read} the exposure {.val {exposure_name}}.",
      i = "An offset is held at its observed value while the exposure is fixed to each level, so the marginal means would read one exposure in the design and another in the offset.",
      i = "For a continuous exposure the same offset leaves the exposure coefficient something other than the effect of a one-unit change.",
      i = "An offset that does not read the exposure, such as the person-time offset of a rate model, is supported."
    ),
    error_class = "balancing_ipw_input_error",
    call = call,
    .envir = environment()
  )
}

# A column name R cannot parse as a symbol is back-quoted wherever the model
# writes it down, so a term label, a coefficient name, and a design column all
# carry the quotes while the fit records the bare name. Comparisons against any
# of those run through here, which leaves a syntactic name exactly as it is.
quoted_name <- function(name) {
  deparse1(as.name(name), backtick = TRUE)
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
      "{.arg outcome_mod} must be fitted with the weights from {.arg wt_mod}.",
      x = "Its weights differ from the fit's, compared per unit at relative tolerance 1e-6.",
      i = "Refit it with {.code weights = weights(fit)}, where {.arg fit} is the balancing fit.",
      i = "A fit with sampling weights composes them into {.code weights(fit)}, so the outcome model takes that composed vector rather than either factor on its own."
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

# The effect rows, and the shared column contract every exposure's estimates
# table keeps. Each effect is a parameter of the stacked system, so its estimate
# is that parameter's entry in the stacked parameter vector and its standard
# error is the square root of the matching diagonal entry of the covariance.
# Nothing is contrasted or differentiated here, which is what keeps the contrast
# formulas stated once, in the stack itself.
#
# `keys` names the stacked entries the rows are read from and `effects` names
# what each row is called, which are the same strings for a continuous exposure,
# whose one effect is the exposure coefficient under the name that coefficient
# already carries in the stack, and differ for a categorical one, whose measures
# repeat once per comparison under distinct stacked names.
ipw_estimate_rows <- function(theta, vcov, conf_level, keys, effects) {
  estimate <- unname(theta[keys])
  std_err <- unname(sqrt(diag(vcov)[keys]))
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

# A categorical exposure reports one block of measures per non-reference level,
# so the `effect` column alone no longer identifies a row: the same three
# measures appear once per comparison. The table therefore gains a `comparison`
# column naming the two levels, placed immediately after `effect`, which is where
# propensity puts it and where its own print method looks for it. A binary
# exposure has a single comparison and keeps the eight-column table, since a
# column repeating one label on every row identifies nothing.
ipw_estimates <- function(theta, vcov, conf_level, continuous, levels = NULL) {
  keys <- ipw_contrast_names(continuous, levels)
  measures <- ipw_contrast_names(continuous)
  estimates <- ipw_estimate_rows(
    theta = theta,
    vcov = vcov,
    conf_level = conf_level,
    keys = keys,
    effects = rep(measures, times = length(keys) / length(measures))
  )
  if (is.null(levels)) {
    return(estimates)
  }

  comparison <- rep(
    paste(levels[-1], "vs", levels[[1]]),
    each = length(measures)
  )
  cbind(
    estimates["effect"],
    comparison = comparison,
    estimates[setdiff(names(estimates), "effect")]
  )
}
