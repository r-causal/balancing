# Inverse probability weighting for a balancing fit

A
[balancing](https://r-causal.github.io/balancing/reference/balancing.md)
fit registers a method on
[`causalgenerics::ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html),
so a set of balancing weights drives the same bring-your-own-model
workflow as a propensity score fit. You supply the fit and a weighted
outcome model, and
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
returns causal effect estimates with standard errors that account for
having estimated the weights.

## Arguments

- wt_mod:

  A
  [balancing](https://r-causal.github.io/balancing/reference/balancing.md)
  fit that produced the weights.

- outcome_mod:

  A weighted outcome model of class
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) or
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html), fitted with the
  balancing weights and carrying the exposure among its predictors. It
  may adjust for covariates alongside the exposure. For a continuous
  exposure it is a marginal structural model carrying exactly one term
  in the exposure.

- .data:

  The data frame holding the exposure and outcome. If `NULL`, the values
  are taken from the outcome model frame. It carries the exposure column
  the fixed-exposure predictions are built from, so it has nothing to
  supply for a continuous exposure, which makes no such predictions, and
  is ignored there.

- estimand:

  The causal estimand. If `NULL`, the fit's estimand is used. As in
  [`balance()`](https://r-causal.github.io/balancing/reference/balance.md),
  `"atc"` is accepted as a synonym for `"atu"`. Supplying an estimand
  that disagrees with the fit raises `balancing_estimand_error`.

- conf_level:

  The confidence level for the intervals. Default `0.95`.

- ...:

  Ignored, for compatibility with the generic.

## Value

An object of class `ipw`, an implementation of
[`causalgenerics::ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html).
Alongside `estimand`, `wt_mod`, `outcome_mod`, and the `estimates`
table, the result carries two fields describing the variance:

- `se_method`, the string `"mestimation"`, naming how the standard
  errors were computed.

- `fit`, the fitted variance system, a list of `theta`, the stacked
  parameter vector, and `vcov`, its sandwich covariance. Both are named
  by stacked block: `theta_w1` onward for the weight parameters, `beta_`
  followed by the design column name for the outcome-model coefficients,
  then the marginal means and one name per contrast. A binary exposure
  names its means `mu0` and `mu1` and its contrasts by measure alone; a
  categorical exposure names each mean `mu_` followed by its level and
  each contrast by measure and level, as `rd_b`. A continuous exposure
  carries neither block, since the effect is one of the outcome-model
  coefficients; that coefficient is named for the effect, as `slope`, in
  place of the `beta_` name the others carry. The standard errors in
  `estimates` are `sqrt(diag(fit$vcov))` read at those effect names.

## Details

The point estimates are the g-computation marginal means: the outcome
model is predicted with the exposure fixed to each level and averaged
over the estimand's target population. For a binary outcome the method
returns the risk difference (`rd`), the log risk ratio (`log(rr)`), and
the log odds ratio (`log(or)`); for a continuous outcome it returns the
difference in means (`diff`).

A categorical exposure reports those same measures for each
non-reference level against the reference level, which is the first of
the fit's own levels: a factor's declared level order for a factor
exposure, and the sorted values otherwise. A K-level exposure therefore
contributes K marginal means and one block of measures per non-reference
level, and the estimates table gains a `comparison` column, placed after
`effect`, naming each contrast as `"<level> vs <reference>"`. A binary
exposure keeps the table it has always returned, with no `comparison`
column.

A continuous exposure has no levels to contrast, so there is no pair of
marginal means to difference. What the method reports instead is the
dose-response coefficient of a weighted marginal structural model: the
balancing weights break the exposure-covariate association, the outcome
model carries exactly one term in the exposure, and that term's
coefficient is the effect of a one-unit change in the exposure on the
model's own link scale. The estimates table holds a single row, keeping
the columns it holds for every other exposure and gaining no
`comparison` column, and that row is named for the link: `slope` for an
identity link, whether the model arrives as a
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) or as a gaussian
[`stats::glm()`](https://rdrr.io/r/stats/glm.html); `log(or)` for a
logit; and `log(rr)` for a log link. Another link raises
`balancing_ipw_input_error`, since its coefficient is none of those
three. A continuous fit targets the average treatment effect and nothing
else, so an `estimand` supplied alongside it either agrees or raises
`balancing_estimand_error`, as it does for any other fit.

The standard errors come from a stacked M-estimator that the deli
package differentiates and sandwiches. The stacked parameter vector
holds four blocks: the weight parameters, the outcome-model
coefficients, the marginal means, one per exposure level, and the effect
contrasts. Their estimating functions are, in the same order, the
weight-parameter estimating equations the fit carries, re-evaluated at
new weight parameters through the hooks the fit's
[balancing_estimating_equations](https://r-causal.github.io/balancing/reference/balancing_estimating_equations.md)
container supplies; the outcome-model score, from
[`deli::ee_glm()`](https://r-causal.github.io/deli/reference/ee_glm.html),
carrying the balancing weights as the fit would have reported them at
those parameters, per-group renormalization included; the marginal-mean
equations, which predict the outcome model with the exposure fixed to
each level and standardize over the estimand's target population; and
one deterministic row per contrast, setting the contrast parameter equal
to its formula in the two means.
[`deli::compute_sandwich()`](https://r-causal.github.io/deli/reference/compute_sandwich.html)
differentiates that system at the fitted values and returns its
empirical sandwich covariance.

A continuous exposure stacks the same system with the g-computation half
removed. There are no fixed-exposure predictions to standardize and no
contrasts to form, so the stack holds the weight parameters and the
marginal structural model's coefficients alone, and the effect is
already one of those coefficients rather than a parameter derived from
them.

Nothing is re-solved along the way. Every parameter enters at the value
its own fit already found, and the stacked estimating functions are only
re-evaluated around that point. Because each contrast is a parameter of
the stack rather than a transformation applied afterward, its standard
error is already on the diagonal of the joint covariance and no
delta-method step stands between the sandwich and the reported effects.
The joint covariance propagates the uncertainty from estimating the
weights into the effect standard errors, which a variance that treats
the weights as fixed would understate.

That covariance is a large-sample one, and how large a sample it takes
differs by exposure. The binary risk-difference standard error is
calibrated at a few hundred observations; the continuous slope's is
anticonservative there. Over 500 draws its ratio of mean standard error
to the standard deviation of the estimates was 0.836 at 300 observations
and 0.942 at 1200. Read a continuous interval as the asymptotic
statement it is, and prefer the bootstrap of the inference vignette at
small sample sizes.

The method is available only for fits whose weights solve smooth
estimating equations: the estimating-equation family (entropy balancing,
inverse probability tilting, and the just-identified covariate balancing
propensity score) with exact balance at a binary or categorical
exposure, and entropy balancing with exact balance at a continuous one.
Any other fit, including an entropy fit at a positive tolerance, an
over-identified or quadratic-program fit, or a continuous fit from a
method that solves no estimating equations, raises
`balancing_ipw_unsupported_error` and points to the bootstrap workflow
described in the inference vignette.

The outcome model must carry the same exposure levels the fit weighted.
Fitting it on data that drop a level, or that carry one the fit never
saw, raises `balancing_ipw_input_error`: the counterfactual predictions
and the weights would then describe different exposures, and the
resulting table would name a contrast it had not computed.

The outcome model must carry the exposure among its predictors, and may
adjust for covariates alongside it, including in interactions with the
exposure: `y ~ exposure`, `y ~ exposure + x1 + x2`, and
`y ~ exposure * x1` are all supported. A categorical exposure enters as
a factor, so its fixed-exposure designs come from the model's own
contrasts. A model without an exposure term raises
`balancing_ipw_input_error`, since its fixed-exposure predictions would
all be the same prediction and every contrast it reported would be zero.

A continuous exposure narrows that to exactly one term, the exposure
itself, since what it reports is one coefficient of the model rather
than a contrast of predictions from it. `y ~ exposure` and
`y ~ exposure + x1` are supported, while `y ~ exposure + I(exposure^2)`,
`y ~ poly(exposure, 2)`, and `y ~ exposure * x1` raise
`balancing_ipw_input_error`: each carries a second design column in the
exposure, so the effect of a one-unit change depends on where it is read
and no single coefficient is it. The check reads the model's terms
rather than the text of its formula, so a transformed or interacted
exposure term is caught however it is written.

An offset reaches the linear predictor without being a term, so it is
checked on its own, and for every exposure type. An offset expression
naming the exposure, written either into the formula or passed through
the model's `offset` argument, raises `balancing_ipw_input_error`. An
offset is held at its observed value while the exposure is fixed to each
level, so a discrete exposure's marginal means would read one exposure
in the design and another in the offset, and a continuous exposure's
coefficient would be something other than the effect of a one-unit
change. An exposure-free offset stays supported. That check is static,
so it accepts any offset whose stored expression does not name the
exposure: a precomputed vector under some other symbol, and equally a
wrapper that forwards the offset through its dots, which records `..1`
in the fitted call. Keeping the exposure out of such an offset is the
caller's to honor, and it matters most for a discrete exposure, where a
laundered offset corrupts the fixed-exposure marginal means and can
reverse the sign of the reported contrast rather than merely shifting a
coefficient.

The exposure may be a factor, a character column, or an integer code,
and a term that transforms it counts as carrying it. A character column
becomes a factor in a model formula on its own, while an integer code is
written `factor(exposure)` for a model with one parameter per level;
left untransformed, an integer code enters as a slope in the codes and
the marginal means are the g-computation means of that model rather than
of a saturated one. A transformed exposure is not a column of the model
frame, since the frame stores the transformation, so such a model is
passed with `.data`.

Which population the marginal means are averaged over is part of the
estimand, and matters as soon as the outcome model adjusts for anything.
A marginal model is saturated in the exposure, one free parameter per
exposure level, so absent an offset it predicts a single value per
level, its marginal means are the weighted group means whatever link the
family carries, and no choice of population can change them. An adjusted
model predicts a value per unit, so the means are standardized over the
estimand's target population: every unit for a pooled estimand, and the
focal group's units for `"att"` or `"atc"`. Sampling weights, where the
fit has them, weight that average as well.

The link enters the outcome-model score, and it enters exactly: the
bread is differentiated from the estimating functions themselves rather
than read off an information-matrix formula, so a non-canonical link
such as probit or cloglog is handled exactly rather than approximately,
for an adjusted model as much as for a marginal one.

Two further conditions on the outcome model raise the same
`balancing_ipw_input_error`. Its family must be binomial, quasibinomial,
or gaussian, which includes a plain
[`stats::lm()`](https://rdrr.io/r/stats/lm.html), since the reported
effects are the contrasts derived for those families' marginal means.
And it must have been fitted with the weights the fit produced, since
the stacked variance differentiates the outcome-model score through
those weights: the model's weights are compared against the fit's, per
unit at a relative tolerance of 1e-6. Those are the weights
`weights(fit)` returns, which already carry the fit's sampling weights
if it has any. Sampling weights compose multiplicatively onto the
balancing weights and the stack holds them fixed, since they are a
design quantity rather than an estimate.

Of the two binomial families,
[`stats::quasibinomial()`](https://rdrr.io/r/stats/family.html) is the
one to fit a binary outcome with here, and the examples below use it:
balancing weights are not counts, so
[`stats::binomial()`](https://rdrr.io/r/stats/family.html) warns about
non-integer successes at every fit. The two solve the same estimating
equation, since they share the binomial variance function, and the
dispersion the quasi family estimates never reaches the sandwich, which
is built from the score alone. The results are therefore identical
rather than merely close.

An offset is supported, written either as an
[`offset()`](https://rdrr.io/r/stats/offset.html) term in the outcome
formula or passed through the model's `offset` argument, so long as it
does not read the exposure. It is carried through both the outcome-model
score and the fixed-exposure linear predictors, so the marginal means
are the g-computation means with each unit's offset held at its observed
value. That is the right treatment of a quantity the exposure does not
move and the wrong treatment of one it does, which is why the
exposure-reading case is refused above.

## References

Kostouraki A, Hajage D, Rachet B, et al. On variance estimation of the
inverse probability-of-treatment weighting estimator: A tutorial for
different types of propensity score weights. *Statistics in Medicine*.
2024;43(13):2672-2694.
[doi:10.1002/sim.10078](https://doi.org/10.1002/sim.10078)

## Examples

``` r
n <- 200
x1 <- rnorm(n)
z <- rbinom(n, 1, plogis(0.5 * x1))
y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.3 * x1))
df <- data.frame(exposure = z, x1 = x1, y = y)

fit <- balance(df, exposure, x1, method = bw_entropy(), estimand = "ate")
#> ℹ Treating `.exposure` as binary.
df$.wts <- weights(fit)

# quasibinomial() solves the same estimating equation as binomial() and does
# not warn that weights are not counts.
outcome_mod <- glm(
  y ~ exposure,
  data = df,
  family = quasibinomial(),
  weights = .wts
)

ipw(fit, outcome_mod)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> 
#> Propensity Score Model:
#>   Call: NULL 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ exposure, family = quasibinomial(), data = df, 
#>     weights = .wts) 
#> 
#> Estimates:
#>         estimate  std.err        z ci.lower ci.upper conf.level   p.value    
#> rd       0.28118 0.066447 4.231587   0.1509  0.41141       0.95 2.320e-05 ***
#> log(rr)  0.59317 0.154680 3.834844   0.2900  0.89634       0.95 0.0001256 ***
#> log(or)  1.15663 0.288897 4.003589   0.5904  1.72285       0.95 6.239e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# The outcome model may also adjust for covariates, in which case the
# marginal means are standardized over the estimand's target population.
adjusted_mod <- glm(
  y ~ exposure + x1,
  data = df,
  family = quasibinomial(),
  weights = .wts
)

ipw(fit, adjusted_mod)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> 
#> Propensity Score Model:
#>   Call: NULL 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ exposure + x1, family = quasibinomial(), data = df, 
#>     weights = .wts) 
#> 
#> Estimates:
#>         estimate std.err       z ci.lower ci.upper conf.level   p.value    
#> rd       0.28010 0.06674 4.19690   0.1493  0.41091       0.95 2.706e-05 ***
#> log(rr)  0.59109 0.15589 3.79170   0.2856  0.89664       0.95 0.0001496 ***
#> log(or)  1.15198 0.29001 3.97213   0.5836  1.72039       0.95 7.123e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# A categorical exposure reports each level against the reference level, and
# the estimates table names the comparison.
odds_b <- exp(0.6 * x1)
odds_c <- exp(-0.5 * x1)
denominator <- 1 + odds_b + odds_c
draw <- runif(n)
df$arm <- factor(
  ifelse(
    draw < 1 / denominator,
    "a",
    ifelse(draw < (1 + odds_b) / denominator, "b", "c")
  ),
  levels = c("a", "b", "c")
)
df$relapse <- rbinom(
  n,
  1,
  plogis(-0.4 + 0.5 * (df$arm == "b") + 0.9 * (df$arm == "c") + 0.3 * x1)
)

arm_fit <- balance(df, arm, x1, method = bw_ipt(), estimand = "ate")
#> ℹ Treating `.exposure` as categorical.
df$.arm_wts <- weights(arm_fit)
arm_mod <- glm(
  relapse ~ arm,
  data = df,
  family = quasibinomial(),
  weights = .arm_wts
)

ipw(arm_fit, arm_mod)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> 
#> Propensity Score Model:
#>   Call: NULL 
#> 
#> Outcome Model:
#>   Call: glm(formula = relapse ~ arm, family = quasibinomial(), data = df, 
#>     weights = .arm_wts) 
#> 
#> Estimates:
#> Warning: non-unique values when setting 'row.names': ‘log(or)’, ‘log(rr)’, ‘rd’
#> Error in `.rowNamesDF<-`(x, value = value): duplicate 'row.names' are not allowed

# A continuous exposure reports one effect, the exposure coefficient of a
# weighted marginal structural model, named for that model's link.
df$dose <- 0.7 * x1 + rnorm(n)
df$score <- 2 + 0.5 * df$dose + 0.4 * x1 + rnorm(n)

dose_fit <- balance(df, dose, x1, method = bw_entropy(), estimand = "ate")
#> ℹ Treating `.exposure` as continuous.
df$.dose_wts <- weights(dose_fit)
dose_mod <- lm(score ~ dose, data = df, weights = .dose_wts)

ipw(dose_fit, dose_mod)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> 
#> Propensity Score Model:
#>   Call: NULL 
#> 
#> Outcome Model:
#>   Call: lm(formula = score ~ dose, data = df, weights = .dose_wts) 
#> 
#> Estimates:
#>       estimate  std.err        z ci.lower ci.upper conf.level   p.value    
#> slope  0.48987 0.076449 6.407895     0.34  0.63971       0.95 1.475e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
