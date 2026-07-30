# an ipw() result prints for a balancing fit

    Code
      print(result)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      
      Weight Estimator:
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

# ipw() requires the exposure among the outcome model's predictors

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `outcome_mod` must include the exposure among its predictors.
      x The exposure "exposure" appears in none of its terms.
      i The model may adjust for covariates alongside the exposure, and may carry the exposure inside a transformation such as `factor()`.

# ipw() rejects an estimand that contradicts the fit

    Code
      stop(cnd)
    Condition <balancing_estimand_error>
      Error in `ipw()`:
      ! The requested `estimand` does not match the fit.
      x The fit targets "ate".
      x You requested "att".

# ipw() rejects an estimand outside the vocabulary

    Code
      stop(cnd)
    Condition <rlang_error>
      Error in `ipw()`:
      ! `estimand` must be one of "ate", "att", "atc", "atu", or "ato", not "bogus".

# ipw() rejects an outcome model that is not a glm or lm

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `outcome_mod` must be a fitted outcome model.
      i Supply a model of class <glm> or <lm>.
      x `outcome_mod` has class <list>.

# ipw() rejects an outcome model fitted without weights

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `outcome_mod` must be fitted with the weights from `wt_mod`.
      x Its weights differ from the fit's, compared per unit at relative tolerance 1e-6.
      i Refit it with `weights = weights(fit)`, where `fit` is the balancing fit.
      i A fit with sampling weights composes them into `weights(fit)`, so the outcome model takes that composed vector rather than either factor on its own.

# the weight-mismatch message points at the composed weights

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `outcome_mod` must be fitted with the weights from `wt_mod`.
      x Its weights differ from the fit's, compared per unit at relative tolerance 1e-6.
      i Refit it with `weights = weights(fit)`, where `fit` is the balancing fit.
      i A fit with sampling weights composes them into `weights(fit)`, so the outcome model takes that composed vector rather than either factor on its own.

# ipw() rejects a poisson outcome model

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `outcome_mod` must come from a supported outcome family.
      x Its family is "poisson".
      i The supported families are "binomial", "quasibinomial", and "gaussian".
      i See the inference vignette for a bootstrap workflow with other families.

# ipw() rejects a factor outcome model fitted without its response

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `outcome_mod` must carry the response it modeled.
      x It has a factor response but was fitted with `y = FALSE`, which discards that response.
      i Refit it with `y = TRUE`, or on a numeric 0/1 response.

# ipw() refuses a grouped binomial outcome model

    Code
      stop(cnd)
    Condition <balancing_ipw_unsupported_error>
      Error in `ipw()`:
      ! `ipw()` cannot compute a stacked variance for a grouped binomial outcome model.
      x Its response is a matrix of 2 columns, so `glm()` scaled the weights it was given by each row's trial count.
      i Fit the weights and the outcome model on data with one row per trial, or see the inference vignette for a bootstrap workflow.

# ipw() rejects a fit whose container carries no re-evaluation hooks

    Code
      stop(cnd)
    Condition <balancing_ipw_unsupported_error>
      Error in `ipw()`:
      ! `ipw()` cannot compute a stacked variance for this balancing fit.
      x This fit's container does not carry re-evaluation hooks, which the stacked variance differentiates the weight path through.
      i The hooks re-evaluate the estimating functions and the reported weights at new weight parameters.
      i See the inference vignette for a bootstrap workflow.

# the unsupported-weights ipw error carries the bootstrap pointer

    Code
      stop(cnd)
    Condition <balancing_ipw_unsupported_error>
      Error in `ipw()`:
      ! `ipw()` cannot compute a stacked variance for this balancing fit.
      x This fit's weights do not solve smooth estimating equations, so the stacked variance is unavailable.
      i Estimating equations come from the estimating-equation family (entropy balancing, inverse probability tilting, just-identified covariate balancing propensity score) with exact balance.
      i See the inference vignette for a bootstrap workflow.

# the stacked variance refuses a deficiency that moves the weights

    Code
      stop(cnd)
    Condition <balancing_ipw_unsupported_error>
      Error in `ipw_deli_sandwich()`:
      ! `ipw()` cannot compute a stacked variance for this balancing fit.
      x Its estimating equations have rank 4 of 5, so the stacked bread is singular.
      x The reported weights move along 1 unidentified direction, which carries the deficiency into the effect standard errors.
      i Refit the weights on covariates whose constraint columns are independent, or see the inference vignette for a bootstrap workflow.

