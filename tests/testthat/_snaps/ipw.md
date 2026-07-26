# an ipw() result prints for a balancing fit

    Code
      print(result)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      
      Propensity Score Model:
        Call: balance(.data = data, .exposure = exposure, .covariates = c(x1, 
          x2), method = bw_entropy(), estimand = "ate") 
      
      Outcome Model:
        Call: stats::glm(formula = formula, family = family, data = data, weights = .wts) 
      
      Estimates:
              estimate  std.err      z   ci.lower ci.upper conf.level p.value  
      rd      0.148050 0.076088 1.9458 -0.0010793  0.29718       0.95 0.05168 .
      log(rr) 0.313825 0.162944 1.9260 -0.0055384  0.63319       0.95 0.05411 .
      log(or) 0.598059 0.311491 1.9200 -0.0124530  1.20857       0.95 0.05486 .
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# ipw() rejects an estimand that contradicts the fit

    Code
      stop(cnd)
    Condition <balancing_estimand_error>
      Error in `propensity::ipw()`:
      ! The requested `estimand` does not match the fit.
      x The fit targets "ate".
      x You requested "att".

# ipw() rejects an outcome model that is not a glm or lm

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `propensity::ipw()`:
      ! `outcome_mod` must be a fitted outcome model.
      i Supply a model of class <glm> or <lm>.
      x `outcome_mod` has class <list>.

# ipw() rejects a covariate-adjusted outcome model

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `propensity::ipw()`:
      ! `outcome_mod` must be the marginal outcome model.
      i The exposure "exposure" must be its only predictor.
      i See the inference vignette for a bootstrap workflow with covariate-adjusted outcome models.

# ipw() rejects an outcome model fitted without weights

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `propensity::ipw()`:
      ! `outcome_mod` must be fitted with the weights from `ps_mod`.
      x Its weights differ from the fit's, compared per unit at relative tolerance 1e-6.
      i Refit it with `weights = weights(fit)`, where `fit` is the balancing fit.

# ipw() rejects an outcome model with an offset term

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `propensity::ipw()`:
      ! `outcome_mod` must not carry an offset.
      x Its linear predictor includes an offset, which the stacked variance does not yet carry.
      i See the inference vignette for a bootstrap workflow.

# ipw() rejects a poisson outcome model

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `propensity::ipw()`:
      ! `outcome_mod` must come from a supported outcome family.
      x Its family is "poisson".
      i The supported families are "binomial", "quasibinomial", and "gaussian".
      i See the inference vignette for a bootstrap workflow with other families.

# ipw() rejects a factor outcome model fitted without its response

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `propensity::ipw()`:
      ! `outcome_mod` must carry the response it modeled.
      x It has a factor response but was fitted with `y = FALSE`, which discards that response.
      i Refit it with `y = TRUE`, or on a numeric 0/1 response.

# ipw() rejects a categorical-exposure fit that has a container

    Code
      stop(cnd)
    Condition <balancing_ipw_unsupported_error>
      Error in `propensity::ipw()`:
      ! `ipw()` cannot compute a stacked variance for this balancing fit.
      x This fit has a categorical exposure, and only binary exposures are supported.
      i The stacked variance is derived for a binary exposure.
      i See the inference vignette for a bootstrap workflow.

# the unsupported-weights ipw error carries the bootstrap pointer

    Code
      stop(cnd)
    Condition <balancing_ipw_unsupported_error>
      Error in `propensity::ipw()`:
      ! `ipw()` cannot compute a stacked variance for this balancing fit.
      x This fit's weights do not solve smooth estimating equations, so the stacked variance is unavailable.
      i Estimating equations come from the estimating-equation family (entropy balancing, inverse probability tilting, just-identified covariate balancing propensity score) with exact balance.
      i See the inference vignette for a bootstrap workflow.

