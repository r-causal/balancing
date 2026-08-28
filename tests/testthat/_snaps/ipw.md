# an ipw() result prints for a balancing fit

    Code
      print(result)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: marginal (population-averaged) 
      
      Weight Estimator:
        Call: balance(.data = data, .exposure = exposure, .covariates = c(x1, 
          x2), method = bw_entropy(), estimand = "ate") 
      
      Outcome Model:
        Call: stats::glm(formula = formula, family = family, data = data, weights = .wts) 
      
      Marginal estimates:
                     estimate  std.err      z   ci.lower ci.upper conf.level
      mean 0         0.401600 0.050558 7.9434  0.3025082  0.50069       0.95
      mean 1         0.549649 0.056862 9.6664  0.4382025  0.66110       0.95
      rd 1 vs 0      0.148050 0.076088 1.9458 -0.0010793  0.29718       0.95
      log(rr) 1 vs 0 0.313825 0.162944 1.9260 -0.0055384  0.63319       0.95
      log(or) 1 vs 0 0.598059 0.311491 1.9200 -0.0124530  1.20857       0.95
                       p.value    
      mean 0         1.967e-15 ***
      mean 1         < 2.2e-16 ***
      rd 1 vs 0        0.05168 .  
      log(rr) 1 vs 0   0.05411 .  
      log(or) 1 vs 0   0.05486 .  
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

# ipw() names the offset when it carries the only mention of the exposure

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `outcome_mod` must include the exposure among its predictors.
      x The exposure "exposure" appears in none of its terms.
      i The model may adjust for covariates alongside the exposure, and may carry the exposure inside a transformation such as `factor()`.
      i An offset is not a term, so an exposure that reaches the model only through one is not carried by it.

# ipw() asks a continuous outcome model for a term that reads the exposure

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `outcome_mod` must include the exposure among its predictors.
      x The exposure "exposure" appears in none of its terms.
      i The model may adjust for covariates alongside the exposure, and may carry the exposure inside a transformation such as `poly()`, so long as no term reads a covariate alongside it.

# ipw() refuses a continuous outcome model with an unnamed link

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `outcome_mod` must use a link the continuous effect can be named for.
      x Its link is "probit".
      i The supported links are "identity", "logit", and "log", whose exposure coefficients are a slope, a log odds ratio, and a log risk ratio.
      i See the inference vignette for a bootstrap workflow with other links.

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

# ipw() rejects an outcome model with an aliased exposure coefficient

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `outcome_mod` must have an estimate for every coefficient.
      x It is rank deficient, so the coefficient "I(2 * exposure)" is not estimable.
      i Drop the aliased term from `outcome_mod` and fit it again before calling `ipw()`.

# ipw() rejects an outcome model with an aliased covariate coefficient

    Code
      stop(cnd)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `outcome_mod` must have an estimate for every coefficient.
      x It is rank deficient, so the coefficient "I(2 * x1)" is not estimable.
      i Drop the aliased term from `outcome_mod` and fit it again before calling `ipw()`.

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
      ! `ipw()` cannot compute a stacked variance for a multi-column response.
      x Its response is a matrix of 2 columns, which is the grouped binomial form for `glm()` and a multivariate fit for `lm()`.
      i A grouped binomial fit scales the weights it was given by each row's trial count; a multivariate fit carries one coefficient block per response column. Neither leaves a single per-unit score the stack can rebuild.
      i Fit the weights and the outcome model on data with one row per unit and one response column, or see the inference vignette for a bootstrap workflow.

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

# ipw() names the rank of a deficient fit block deli refuses

    Code
      stop(cnd)
    Condition <balancing_ipw_unsupported_error>
      Error in `ipw()`:
      ! The stacked variance could not be computed for this outcome model.
      x The stacked bread has no inverse at the fitted parameters.
      x The balancing fit's estimating equations have rank 4 of 5, so the stacked bread is singular.
      i Refit the weights on covariates whose constraint columns are independent.
      i See the inference vignette for a bootstrap workflow.

