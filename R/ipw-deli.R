# The variance for a balancing fit is a stacked M-estimator, and this file
# assembles that stack as a single estimating-function closure that deli
# differentiates and sandwiches. Rather than writing each block of the bread out
# analytically, the whole system is written once as the estimating functions
# themselves and `deli::compute_sandwich()` finite differences the bread from
# them. Nothing is re-solved: every parameter enters at the value its own fit
# already found, and the closure only re-evaluates the estimating functions
# around that point.
#
# The stack is ordered [theta_w | beta | means | contrasts]. The weight
# parameters come first because everything downstream depends on them and
# nothing upstream does. The outcome-model coefficients follow, coupled to the
# weight parameters through the weights the score carries. The marginal means
# read off the outcome model with the exposure fixed to each level, one mean per
# level, and the effect contrasts close the system. Making each contrast a
# parameter of the stack is what removes the delta method from the caller: the
# contrast's standard error is already on the diagonal of the returned
# covariance.
#
# An exposure with K levels contributes K means and one block of contrasts per
# non-reference level, each measured against the reference level, which is the
# first of the fit's own levels. A binary exposure is that system at K equal to
# two, and the only thing its two levels change is the naming: its means are
# `mu0` and `mu1` and its contrasts are unsuffixed, since one contrast needs no
# label to tell it from another.
#
# A continuous exposure has no levels to fix, so the last two blocks are absent
# and the stack stops at the outcome-model coefficients. What it reports is the
# exposure coefficients of a weighted marginal structural model, each of which is
# already a parameter of the stack, so the same reading applies: every effect's
# standard error is on the diagonal of the returned covariance under the name its
# own coefficient carries.
#
# A `.by` request adds two more blocks at the end: the marginal means within
# each stratum of the modifier, then the contrasts of those means within each
# stratum followed by the contrasts of one stratum's contrasts against the
# reference stratum's. They come last on purpose. Nothing already in the stack
# reads a parameter of theirs, so the bread stays block lower triangular, the
# leading block of its inverse is the leading block the ungrouped system
# produces, and the whole-sample rows a grouped result reports are the rows the
# same fit reports without a request rather than a second computation of them.
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
#' @param levels The exposure levels the fit weighted, as strings in the fit's
#'   own order. The first is the reference level every contrast is measured
#'   against.
#' @param categorical Whether the fit's exposure is categorical, which decides
#'   how the mean and contrast blocks are named.
#' @param by The strata a `.by` request named, as `ipw_resolve_by()` resolves
#'   them, or `NULL` when no request was made. A request appends one mean and
#'   one contrast block per stratum, and one contrast block per non-reference
#'   stratum against the reference one, after every block the ungrouped system
#'   carries.
#' @param joint The surface a declared crossing is reported under, as
#'   `ipw_joint_plan()` resolves it, or `NULL` when the exposure declares none.
#'   A declared crossing keeps the mean block and replaces the contrast block
#'   with the simple effects and their interaction, rather than adding to it.
#' @param sampling_weights The fit's sampling weights, or `NULL`.
#' @param focal_level The fit's focal exposure level, or `NULL` for a pooled
#'   estimand. It names the target population the marginal means standardize
#'   over and the group total the reported weights are carried to.
#' @param call The frame a refusal reports as the failing call. The engine is
#'   internal, so a caller who reached it through [ipw()] must be sent to
#'   `ipw()` rather than to a function they cannot go and read; the default
#'   names this function, which is what a direct call deserves.
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
  levels,
  categorical = FALSE,
  by = NULL,
  joint = NULL,
  sampling_weights = NULL,
  focal_level = NULL,
  call = rlang::current_env()
) {
  n <- nrow(frame)
  family <- stats::family(outcome_mod)
  continuous <- is_gaussian_outcome(outcome_mod)
  distribution <- deli_distribution(family)
  outcome <- resolve_outcome_response(outcome_mod)
  design <- stats::model.matrix(outcome_mod)
  offset <- outcome_mod$offset

  # The marginal-mean equations predict the outcome model with the exposure
  # fixed to each level, so they build one design per level from the model's own
  # terms. An offset is part of the linear predictor rather than of the design,
  # so it is carried alongside and added to eta.
  pieces <- lapply(
    resolve_level_values(frame[[exposure_name]], levels),
    function(value) {
      fixed_exposure_pieces(
        outcome_mod,
        frame,
        exposure_name,
        value,
        offset = offset
      )
    }
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
  # below. The groups are built in the fit's level order rather than by
  # splitting, which would sort them, so that the per-group scale this applies
  # is the one the fit itself reported at.
  key <- as.character(frame[[exposure_name]])
  groups <- stats::setNames(
    lapply(levels, function(level) which(key == level)),
    levels
  )
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
  means <- vapply(
    pieces,
    function(piece) {
      sum(tilt * piece$mu) / sum(tilt)
    },
    numeric(1)
  )
  m <- length(means)
  # A declared crossing reports the same means under contrasts written in the
  # two treatments, so it takes the contrast block over rather than sitting
  # beside it. Everything before that block is what it always was, which is what
  # makes the two surfaces agree on the rows both report.
  contrasts <- ipw_contrast_values(means, continuous)
  effects <- ipw_contrast_names(continuous, if (categorical) levels else NULL)
  if (!is.null(joint)) {
    contrasts <- ipw_joint_values(joint, means, continuous)
    effects <- ipw_joint_names(joint)
  }
  k <- length(effects)

  # The stratum blocks are seeded from the same pieces and the same tilt the
  # whole-sample blocks are, restricted to one stratum at a time, and they are
  # appended after every block above rather than interleaved with them. Nothing
  # already in the stack reads a parameter of theirs, so the bread stays block
  # lower triangular and the leading block of the covariance is the ungrouped
  # system's own.
  by_stack <- ipw_by_stack(by, pieces, tilt, continuous, levels, categorical)
  m_by <- length(by_stack$means)
  k_by <- length(by_stack$contrasts)

  theta <- c(
    weight_parameters,
    coefficients,
    means,
    contrasts,
    by_stack$means,
    by_stack$contrasts
  )
  names(theta) <- c(
    paste0("theta_w", seq_len(p)),
    paste0("beta_", colnames(design)),
    ipw_mean_names(levels, categorical),
    effects,
    by_stack$mean_names,
    by_stack$contrast_names
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
  # The renormalization is named on its own because the cache applies it to
  # whichever hook supplied the weights, and the reported weight map is named on
  # its own because the rank check below differences the same map the closure
  # carries, at parameters the closure never asks for.
  rescale <- function(weights) {
    renormalize_group_weights(weights, sampling, groups, targets)
  }
  weights_at <- function(weight_theta) {
    rescale(as.numeric(container@weights_fn(weight_theta)))
  }
  hooks_at <- make_hooks_cache(
    container,
    rescale,
    as.numeric(weight_parameters)
  )

  stacked_equations <- function(theta) {
    beta <- theta[p + seq_len(q)]
    mean_theta <- theta[p + q + seq_len(m)]
    contrast_theta <- theta[p + q + m + seq_len(k)]

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

    # The fixed-exposure predictions at these coefficients, one vector per
    # exposure level. The whole-sample mean rows and every stratum's mean rows
    # standardize the same predictions over different populations, so they are
    # computed once here and weighted twice rather than derived twice.
    fixed <- lapply(seq_len(m), function(j) {
      eta <- as.numeric(pieces[[j]]$design %*% beta)
      if (!is.null(offset)) {
        eta <- eta + offset
      }
      family$linkinv(eta)
    })

    # The mean rows are weighted by the standardization weight, so their root is
    # the mean of the fixed-exposure predictions over the target population
    # rather than over every unit. A marginal model predicts one value per
    # exposure level, which makes the weighted row a constant multiple of the
    # unweighted one and leaves the sandwich exactly where it was.
    mean_rows <- do.call(
      rbind,
      lapply(seq_len(m), function(j) tilt * (fixed[[j]] - mean_theta[[j]]))
    )

    # The contrasts are deterministic functions of the means, so their rows are
    # the same value for every unit. They contribute nothing to the meat at the
    # solution, where that value is zero, and everything to the bread, which is
    # what carries their standard errors without a delta method.
    contrast_rows <- if (is.null(joint)) {
      matrix(
        ipw_contrast_values(mean_theta, continuous) - contrast_theta,
        nrow = k,
        ncol = n
      )
    } else {
      ipw_joint_rows(joint, mean_theta, contrast_theta, continuous, n)
    }

    by_rows <- ipw_by_rows(
      by_stack = by_stack,
      fixed = fixed,
      mean_theta = theta[p + q + m + k + seq_len(m_by)],
      contrast_theta = theta[p + q + m + k + m_by + seq_len(k_by)],
      continuous = continuous,
      n = n
    )

    rbind(
      hooks$psi,
      score,
      mean_rows,
      contrast_rows,
      by_rows$mean,
      by_rows$contrast
    )
  }

  validate_stacked_bread(
    container@jacobian,
    weights_at,
    as.numeric(weight_parameters),
    call = call
  )

  list(
    theta = theta,
    vcov = stacked_covariance(
      stacked_equations,
      theta,
      n,
      container@jacobian,
      call = call
    )
  )
}

#' The deli-backed stacked sandwich for a continuous exposure
#'
#' The same M-estimator with the g-computation half removed. There are no
#' exposure levels to fix and no pair of marginal means to contrast, so the
#' stack is `[theta_w | beta]` alone: the weight parameters, then the
#' coefficients of the weighted marginal structural model whose score they
#' enter. The reported effects are the exposure's own coefficients, so they need
#' no parameters of their own and no contrast rows to carry them. Naming those
#' entries for the labels the estimates table reports them under is what puts
#' them on the same footing as the discrete path's contrasts, whose standard
#' errors are read off the same diagonal under the same names.
#'
#' The weights the score carries come straight from the container's hook. The
#' reported scale for a continuous fit carries the whole sample to one total,
#' and the hook already normalizes to it at every set of parameters, so there is
#' no per-group renormalization to apply on top.
#'
#' @param container The fit's [balancing_estimating_equations].
#' @param outcome_mod The fitted weighted marginal structural model.
#' @param exposure_name The exposure column name, which decides which design
#'   columns the reported effects are read from: those the model's
#'   exposure-reading terms expanded to.
#' @param sampling_weights The fit's sampling weights, or `NULL`.
#' @param call The frame a refusal reports as the failing call, on the same
#'   terms as the discrete engine's.
#'
#' @return A list with `theta`, the stacked parameter vector, and `vcov`, its
#'   covariance on the standard-error scale, both named by stacked block order.
#'
#' @noRd
ipw_deli_msm_sandwich <- function(
  container,
  outcome_mod,
  exposure_name,
  sampling_weights = NULL,
  call = rlang::current_env()
) {
  family <- stats::family(outcome_mod)
  distribution <- deli_distribution(family)
  outcome <- resolve_outcome_response(outcome_mod)
  design <- stats::model.matrix(outcome_mod)
  offset <- outcome_mod$offset
  n <- nrow(design)

  # The sampling weights compose multiplicatively onto the balancing weights and
  # the stack holds them fixed, exactly as it does for a discrete exposure.
  sampling <- sampling_weights %||% rep(1, n)

  weight_parameters <- container@parameters
  p <- length(weight_parameters)
  coefficients <- stats::coef(outcome_mod)
  q <- length(coefficients)

  theta <- c(weight_parameters, coefficients)
  names(theta) <- c(
    paste0("theta_w", seq_len(p)),
    paste0("beta_", colnames(design))
  )
  # The reported effects are coefficients the stack already carries, so naming
  # them is the whole of what this route does with the surface: each
  # exposure-reading column takes the label its estimates row is read under, in
  # place of the `beta_` name a covariate's column keeps.
  surface <- msm_coefficient_identity(outcome_mod, exposure_name, call = call)
  names(theta)[p + surface$columns] <- surface$keys

  weights_at <- function(weight_theta) {
    as.numeric(container@weights_fn(weight_theta))
  }
  hooks_at <- make_hooks_cache(
    container,
    identity,
    as.numeric(weight_parameters)
  )

  stacked_equations <- function(theta) {
    beta <- theta[p + seq_len(q)]
    hooks <- hooks_at(as.numeric(theta[seq_len(p)]))
    score <- deli::ee_glm(
      beta,
      X = design,
      y = outcome,
      distribution = distribution,
      link = family$link,
      weights = hooks$weights * sampling,
      offset = offset
    )
    rbind(hooks$psi, score)
  }

  validate_stacked_bread(
    container@jacobian,
    weights_at,
    as.numeric(weight_parameters),
    call = call
  )

  list(
    theta = theta,
    vcov = stacked_covariance(
      stacked_equations,
      theta,
      n,
      container@jacobian,
      call = call
    )
  )
}

# The container's two hooks are the expensive part of either stack: each crosses
# into the method's own evaluation entrypoint over the whole data set. Both are
# pure functions of the weight block, and the finite difference presents the
# same weight sub-vector many times over, because perturbing a coordinate
# outside that block leaves the sub-vector exactly at its fitted value. Two
# cached entries cover every repeat. The fitted sub-vector is pinned, since each
# of the central difference's two sweeps returns to it in a long run and a
# single most-recent entry would lose it in between; one further entry holds the
# most recent perturbation, which a sweep asks for twice in succession. The keys
# are short numeric vectors, so comparing them outright is cheaper than hashing
# them.
#
# A container may carry the two hooks as one, which is what a method whose
# estimating functions are a transformation of its own weights can offer: the
# pair costs what one of them costs. Every cache miss wants both, so the combined
# hook is read when it is there and the two separate ones when it is not. The
# two paths return the same values to the bit, so which one a container takes is
# a matter of how much work reaching them costs.
#
# The reported scale is the caller's to decide: a discrete exposure renormalizes
# the hook's output per group and a continuous one takes it as it comes. So
# `rescale` is passed in and applied to whichever hook supplied the weights.
make_hooks_cache <- function(container, rescale, parameters) {
  evaluate_hooks <- if (is.null(container@parts_fn)) {
    function(weight_theta) {
      list(
        key = weight_theta,
        weights = rescale(as.numeric(container@weights_fn(weight_theta))),
        psi = t(container@psi_fn(weight_theta))
      )
    }
  } else {
    function(weight_theta) {
      parts <- container@parts_fn(weight_theta)
      list(
        key = weight_theta,
        weights = rescale(as.numeric(parts$weights)),
        psi = t(parts$psi)
      )
    }
  }
  base_hooks <- evaluate_hooks(parameters)
  recent_hooks <- base_hooks
  function(weight_theta) {
    if (identical(weight_theta, base_hooks$key)) {
      return(base_hooks)
    }
    if (!identical(weight_theta, recent_hooks$key)) {
      recent_hooks <<- evaluate_hooks(weight_theta)
    }
    recent_hooks
  }
}

# The empirical sandwich covariance of a stacked system at its root, named by
# stacked block.
#
# A central difference trades truncation error, of order the step squared,
# against cancellation error, of order the double epsilon over the step; deli's
# 1e-9 default sits far into the cancellation regime, where agreement with the
# analytic bread is near 2e-7 rather than the 2e-10 a 1e-6 step reaches.
#
# A bread the engine cannot invert is refused by the engine, under a class of its
# own, from a frame no caller ever wrote. Translating that refusal here is what
# keeps it readable: the caller meets the package's own classed refusal, reported
# at the call they made, carrying the bootstrap route out that every
# unsupported-variance refusal offers. The engine's own account of why it refused
# is left behind rather than chained onto it, because following that account
# means reading `compute_sandwich()`, which is not code the caller can go and act
# on. Both of the engine's reasons, a bread that is not finite and a bread it
# reads as singular, are named in the refusal instead, unless the fit itself
# settles which of them it was.
#
# It often does. The stacked bread is block lower triangular in the fit's own
# Jacobian, so a rank-deficient fit block makes the whole stack singular, and the
# container carries that Jacobian analytically. Reading its rank by the same rule
# `validate_stacked_bread()` applies leaves nothing ambiguous: the refusal names
# the rank it found and the constraint columns to go and look at, in place of the
# reading that says the estimating functions may not be finite. Reaching here
# with a deficient block is the tolerated case, a direction the reported weights
# are flat along, which nothing reported is read from and which the engine may
# still refuse to invert on one platform and accept on another. A full-rank
# block leaves both readings open, since the fit is then not what went wrong,
# and the generic account stands unchanged.
#
# A bread holding missing values may also come back as a value rather than a
# condition, which would otherwise surface much later as a complaint about
# dimnames applied to a non-array. Name the real cause here instead, at the point
# where it is still legible.
stacked_covariance <- function(
  stacked_equations,
  theta,
  n,
  jacobian,
  call = rlang::caller_env()
) {
  covariance <- rlang::try_fetch(
    deli::compute_sandwich(
      stacked_equations,
      theta,
      deriv_method = "capprox",
      dx = 1e-6,
      allow_pinv = FALSE
    ) /
      n,
    deli_bread_not_invertible = function(cnd) {
      deficiency <- measure_jacobian_rank(jacobian)
      rank <- deficiency$rank
      parameter_count <- deficiency$count
      reading <- if (rank < parameter_count) {
        c(
          x = "The balancing fit's estimating equations have rank {rank} of {parameter_count}, so the stacked bread is singular.",
          i = "Refit the weights on covariates whose constraint columns are independent."
        )
      } else {
        c(
          i = "Either the stacked estimating functions are not finite around the fit, or the stacked bread is singular there."
        )
      }
      abort(
        c(
          "The stacked variance could not be computed for this outcome model.",
          x = "The stacked bread has no inverse at the fitted parameters.",
          reading,
          i = "See the inference vignette for a bootstrap workflow."
        ),
        error_class = "balancing_ipw_unsupported_error",
        call = call,
        .envir = environment()
      )
    }
  )

  if (!is.matrix(covariance)) {
    abort(
      c(
        "The stacked variance could not be computed for this outcome model.",
        x = "The stacked estimating functions are not finite at the fitted parameters.",
        i = "See the inference vignette for a bootstrap workflow."
      ),
      error_class = "balancing_ipw_unsupported_error",
      call = call
    )
  }
  dimnames(covariance) <- list(names(theta), names(theta))
  covariance
}

# The surface a continuous exposure reports: which columns of the outcome design
# carry the dose response, what each of their rows is called, and the stacked
# name each row is read under. Nothing is standardized on this route, so a row is
# exactly a coefficient of the weighted fit and the description is the whole of
# what separates one shape of marginal structural model from another.
#
# The columns are found by variable membership rather than by name. Every term
# reading the exposure reads it alone, which the validator has already settled,
# so the columns those terms expanded to are the exposure's however the terms are
# written: one for a bare or singly transformed exposure, several for a curve
# written out term by term or handed to a basis constructor. The design's own
# `assign` attribute is what maps a term back to its columns, which is why a
# basis needs no frame and no rebuilding of the design.
#
# The rows are named on the same footing as every other surface in the package.
# An exposure entering through one column is the whole of the dose response, so
# its coefficient is that response's slope everywhere, the row keeps the word the
# link gives it, and nothing further is named: a contrast column repeating one
# value down a one-row table would read as a contrast that was named. An exposure
# entering through several has no such row, since a curve has a different slope
# at every dose, so each row is named after the coefficient it reports and the
# scale word steps back to `coef` at an identity link.
msm_coefficient_identity <- function(
  outcome_mod,
  exposure_name,
  call = rlang::caller_env()
) {
  design <- stats::model.matrix(outcome_mod)
  reads_exposure <- vapply(
    model_term_variable_sets(outcome_mod),
    function(variables) exposure_name %in% variables,
    logical(1)
  )
  columns <- which(attr(design, "assign") %in% which(reads_exposure))

  named <- length(columns) > 1L
  effect <- msm_effect_name(outcome_mod, named = named, call = call)
  contrast <- if (named) colnames(design)[columns]
  list(
    columns = columns,
    keys = if (named) paste(effect, contrast) else effect,
    effect = rep(effect, length(columns)),
    contrast = contrast,
    group = NULL
  )
}

# The name a reported continuous effect is measured under, which is set by the
# outcome model's link, since the link is what a coefficient of that model is a
# one-unit effect on. A logit moves the log odds, so the coefficient is a log
# odds ratio; a log link moves the log mean, so it is a log risk ratio. Both
# words are honest of a coefficient whatever column it multiplies, so neither
# depends on how the exposure entered.
#
# An identity link is the one that does. It moves the mean itself, so a lone
# exposure column's coefficient is the slope of the dose response, while one of
# several is no slope at all and the row is named after its own coefficient
# instead, under the word `coef`, which claims nothing beyond that. Keeping the
# two apart is what stops `slope` from meaning an evaluated-at claim in one place
# and a bare coefficient in another.
#
# Another link leaves the coefficient without a name of any of those kinds: a
# probit coefficient is a shift in a latent standard normal scale, which is
# neither a slope on the response nor the log of any ratio, and reporting it
# under one of those labels would name an effect the model does not estimate.
# Such a model is refused rather than labeled.
msm_effect_name <- function(
  outcome_mod,
  named = FALSE,
  call = rlang::caller_env()
) {
  link <- stats::family(outcome_mod)$link
  supported <- c("identity", "logit", "log")
  if (!link %in% supported) {
    abort(
      c(
        "{.arg outcome_mod} must use a link the continuous effect can be named for.",
        x = "Its link is {.val {link}}.",
        i = "The supported links are {.val {supported}}, whose exposure coefficients are a slope, a log odds ratio, and a log risk ratio.",
        i = "See the inference vignette for a bootstrap workflow with other links."
      ),
      error_class = "balancing_ipw_input_error",
      call = call
    )
  }
  switch(
    link,
    identity = if (named) "coef" else "slope",
    logit = "log(or)",
    log = "log(rr)"
  )
}

# The rank of an estimating-equation Jacobian, and the directions its deficiency
# lies along.
#
# A rank read off a decomposition needs a tolerance, since a deficiency arrived
# at in floating point leaves a singular value of rounding size rather than an
# exact zero. The cutoff is relative to the largest singular value, which is what
# makes it independent of the scale the estimating equations happen to be written
# at, and it is set far below any singular value a well-conditioned fit block
# produces and far above the rounding a deficient one leaves. Both readers of a
# deficiency, the check made before the stack is differenced and the refusal
# raised when the engine cannot invert it, measure it here so the two never
# disagree about whether a given fit is deficient.
measure_jacobian_rank <- function(jacobian) {
  decomposition <- svd(jacobian)
  deficient <- decomposition$d <= decomposition$d[[1]] * 1e-8
  list(
    rank = sum(!deficient),
    count = ncol(jacobian),
    directions = decomposition$v[, deficient, drop = FALSE]
  )
}

# Refuse a stacked system whose weight parameters are not identified in a
# direction the reported effects can see.
#
# The stacked bread is block lower triangular: the weight-parameter equations
# depend on the weight parameters alone, so the system's determinant is that
# block's times the rest, and the stack is singular whenever the fit's own
# Jacobian is. The converse needs the stack's other diagonal blocks to be
# nonsingular, which holds for a converged glm with estimable coefficients and,
# on the discrete path, for the mean and contrast blocks, whose diagonals are
# minus the standardization total and minus one. With those blocks nonsingular a
# singular stack means a singular weight block and nothing else, so the
# deficiency this check has to detect is exactly the one the fit's own Jacobian
# carries. Asking the finite difference to report that is asking too much of it.
# A deficiency that is a second-order cancellation comes back as a pivot of
# rounding size rather than as a zero, and `solve()` accepts it, so
# `allow_pinv = FALSE` refuses only the deficiencies that survive to the last bit
# and answers the rest with standard errors resting on rounding error. The check
# is therefore made on the analytic Jacobian the container carries, which is the
# matrix the difference is approximating.
#
# A deficiency is not automatically fatal, and the distinction is not a judgment
# call. Everything downstream of the weight block reads the weight parameters
# only through the reported weight map, so a direction the map is flat along
# cannot reach the outcome-model score, the marginal means, or the contrasts: it
# contaminates the weight block of the covariance, which nothing reported is read
# from. That is the factor-covariate case, where a factor's level indicators sum
# to the constant function and the level-sum direction rescales each group's
# weights by a constant the renormalization divides out again. A direction the
# map is not flat along is the opposite case: the effects inherit the inverse of
# a pivot that is rounding error, and their standard errors come back finite and
# wrong by whatever factor that pivot happened to take.
#
# So flatness is measured rather than assumed. Each deficient direction is
# differenced through the same reported weight map the stack carries, at the step
# the bread is differenced at, and the movement is read against the scale of the
# weights themselves. The two populations are far apart: on the package's own
# fixtures a flat direction moves the weights by around 1e-10 relative and a live
# direction by order one, so no threshold between them is delicate.
validate_stacked_bread <- function(
  jacobian,
  weights_at,
  parameters,
  call = rlang::caller_env()
) {
  deficiency <- measure_jacobian_rank(jacobian)
  if (deficiency$rank == deficiency$count) {
    return(invisible(NULL))
  }

  step <- 1e-6
  magnitude <- max(1, max(abs(weights_at(parameters))))
  directions <- deficiency$directions
  movement <- vapply(
    seq_len(ncol(directions)),
    function(j) {
      shift <- step * directions[, j]
      derivative <- (weights_at(parameters + shift) -
        weights_at(parameters - shift)) /
        (2 * step)
      max(abs(derivative)) / magnitude
    },
    numeric(1)
  )
  # A movement that is not a number counts as movement. The weight map then has
  # no derivative along that direction at all, which is not the flatness the
  # tolerance measures, and a comparison against a missing value must not be left
  # to decide the branch below: it would stop the check with a base error naming
  # neither the fit nor the deficiency, and one such direction would carry the
  # classification of every other direction with it.
  moving <- sum(is.na(movement) | movement > 1e-6)
  if (moving == 0L) {
    return(invisible(NULL))
  }

  rank <- deficiency$rank
  parameter_count <- deficiency$count
  abort(
    c(
      "{.fun ipw} cannot compute a stacked variance for this balancing fit.",
      x = "Its estimating equations have rank {rank} of {parameter_count}, so the stacked bread is singular.",
      x = "The reported weights move along {moving} unidentified direction{?s}, which carries the deficiency into the effect standard errors.",
      i = "Refit the weights on covariates whose constraint columns are independent, or see the inference vignette for a bootstrap workflow."
    ),
    error_class = "balancing_ipw_unsupported_error",
    call = call,
    .envir = environment()
  )
}

# The effect contrasts of the marginal means, and the labels the estimates table
# reports them under. They are parameters of the stack rather than a post-hoc
# transformation, so the values and the names are needed separately: the values
# to seed the stack, the names to label its blocks.
#
# Every contrast is measured against the reference level, the first of the means,
# and the blocks run level-major: each non-reference level contributes all of its
# effect measures before the next level begins. A binary exposure has one such
# block and no need to say which contrast it belongs to, so its labels are the
# bare measure names; a categorical exposure suffixes each label with the level
# it compares, which is what keeps the names unique across blocks.
#
# `collapsible_only` drops the log odds ratio, which is the set of measures a
# stratum and a contrast of strata report. An odds ratio is noncollapsible: the
# odds ratio over a sample is not an average of the odds ratios within its
# subgroups, and the difference of two of them is not the difference in effect
# it reads as. The whole-sample rows keep it, since nothing there averages
# anything over subgroups.
ipw_contrast_values <- function(means, continuous, collapsible_only = FALSE) {
  reference <- means[[1]]
  values <- lapply(means[-1], function(mu) {
    if (continuous) {
      return(mu - reference)
    }
    odds <- if (collapsible_only) {
      numeric(0)
    } else {
      log(mu / (1 - mu)) - log(reference / (1 - reference))
    }
    c(mu - reference, log(mu) - log(reference), odds)
  })
  unlist(values, use.names = FALSE)
}

ipw_contrast_names <- function(
  continuous,
  levels = NULL,
  collapsible_only = FALSE
) {
  measures <- if (continuous) "diff" else c("rd", "log(rr)", "log(or)")
  if (collapsible_only) {
    measures <- setdiff(measures, "log(or)")
  }
  if (is.null(levels)) {
    return(measures)
  }
  unlist(
    lapply(levels[-1], function(level) paste0(measures, "_", level)),
    use.names = FALSE
  )
}

# The seed values and the stacked names of the blocks a `.by` request appends,
# or `NULL` when no request was made. The means come first, one per exposure
# level per stratum, stratum-major so a stratum's tuple sits together; then the
# contrasts of those means, once per stratum, in the order the whole-sample
# contrast block runs in and over the measures a stratum reports; then those
# same contrasts once more for each non-reference stratum against the reference
# one.
#
# The standardization weight is the whole-sample tilt restricted to a stratum,
# which is what makes a stratum's means the g-computation means over that
# stratum's share of the estimand's target population: every unit of it for a
# pooled estimand and its focal units for a focal one. Reading the tilt from the
# whole-sample block rather than rebuilding it is what keeps the two readings of
# a focal estimand from drifting apart.
ipw_by_stack <- function(by, pieces, tilt, continuous, levels, categorical) {
  if (is.null(by)) {
    return(NULL)
  }

  strata <- length(by$labels)
  tilts <- lapply(
    seq_len(strata),
    function(s) tilt * by$indicators[, s]
  )
  means <- lapply(tilts, function(weight) {
    vapply(
      pieces,
      function(piece) sum(weight * piece$mu) / sum(weight),
      numeric(1)
    )
  })
  stratum_contrasts <- lapply(
    means,
    function(mu) ipw_contrast_values(mu, continuous, collapsible_only = TRUE)
  )
  contrast_names <- ipw_contrast_names(
    continuous,
    if (categorical) levels else NULL,
    collapsible_only = TRUE
  )

  list(
    tilts = tilts,
    strata = strata,
    per_stratum = length(contrast_names),
    means = unlist(means, use.names = FALSE),
    contrasts = c(
      unlist(stratum_contrasts, use.names = FALSE),
      unlist(
        lapply(
          seq_len(strata - 1L),
          function(s) stratum_contrasts[[s + 1L]] - stratum_contrasts[[1L]]
        ),
        use.names = FALSE
      )
    ),
    mean_names = ipw_by_names(
      ipw_mean_names(levels, categorical),
      by$labels
    ),
    contrast_names = ipw_by_names(
      contrast_names,
      c(by$labels, by$em_labels)
    )
  )
}

# The stratum blocks of one evaluation of the stacked estimating functions.
#
# The mean rows are the whole-sample mean rows restricted to a stratum: the root
# of the row weighted by that stratum's tilt is the tilt-weighted mean of the
# fixed-exposure predictions over the stratum's units. They read the predictions
# the whole-sample rows were built from rather than predicting again, so the two
# sets of means cannot describe different counterfactuals.
#
# The contrast rows come in two groups, both deterministic and so constant
# across units. A stratum's contrasts transform that stratum's means the way the
# whole-sample contrasts transform theirs. The rows contrasting two strata are
# read off the stratum contrast parameters rather than recomputed from the
# means, so each is the difference of two parameters the system already carries
# and its derivative is exact.
ipw_by_rows <- function(
  by_stack,
  fixed,
  mean_theta,
  contrast_theta,
  continuous,
  n
) {
  if (is.null(by_stack)) {
    return(list(mean = NULL, contrast = NULL))
  }

  n_levels <- length(fixed)
  per_stratum <- by_stack$per_stratum
  mean_rows <- do.call(
    rbind,
    lapply(seq_along(mean_theta), function(row) {
      stratum <- (row - 1L) %/% n_levels + 1L
      level <- (row - 1L) %% n_levels + 1L
      by_stack$tilts[[stratum]] * (fixed[[level]] - mean_theta[[row]])
    })
  )

  stratum_values <- unlist(
    lapply(seq_len(by_stack$strata), function(s) {
      ipw_contrast_values(
        mean_theta[(s - 1L) * n_levels + seq_len(n_levels)],
        continuous,
        collapsible_only = TRUE
      )
    }),
    use.names = FALSE
  )
  em_values <- unlist(
    lapply(seq_len(by_stack$strata - 1L), function(s) {
      contrast_theta[s * per_stratum + seq_len(per_stratum)] -
        contrast_theta[seq_len(per_stratum)]
    }),
    use.names = FALSE
  )

  list(
    mean = mean_rows,
    contrast = matrix(
      c(stratum_values, em_values) - contrast_theta,
      nrow = length(contrast_theta),
      ncol = n
    )
  )
}

# The names of the marginal-mean block. A categorical exposure names each mean
# after its level, which is what the estimates table's contrast labels and the
# contrast names are keyed on. A binary exposure keeps the positional `mu0` and
# `mu1`, the names its results have always carried and the ones propensity uses
# for the same block.
ipw_mean_names <- function(levels, categorical) {
  if (categorical) {
    return(paste0("mu_", levels))
  }
  paste0("mu", seq_along(levels) - 1L)
}

# The exposure values that stand for the fit's levels, taken from the data
# itself. The fit records its levels as strings, while the exposure column may be
# a factor, a character vector, or numeric, and fixing the exposure to a level
# has to leave that column the type the outcome model was fitted on: assigning
# the string "1" into a numeric exposure would turn the design column into a
# factor contrast the fitted coefficients do not describe. Taking the value from
# the first observation at each level preserves the type exactly, factor levels
# included.
resolve_level_values <- function(exposure, levels) {
  key <- as.character(exposure)
  lapply(levels, function(level) exposure[[match(level, key)]])
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
