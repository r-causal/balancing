# Inference after balancing

Balancing weights are estimated, not fixed. A weighted outcome model
treats its weights as known constants, so the standard errors it reports
do not account for the uncertainty in the weights themselves. This
article explains why that matters and shows two ways to get honest
standard errors: M-estimation through
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
for the estimating-equation methods, and a bootstrap for everything
else.

``` r

library(balancing)
```

## Why the weights’ uncertainty matters

Consider a binary exposure with a binary outcome and two confounders.

``` r

n <- 1500
age <- rnorm(n)
score <- rnorm(n)
exposure <- rbinom(n, 1, plogis(0.5 * age - 0.5 * score))
event <- rbinom(n, 1, plogis(-0.4 + 0.7 * exposure + 0.4 * age - 0.3 * score))

study <- data.frame(exposure, age, score, event)
```

We fit entropy balancing weights for the ATE and a weighted outcome
model.

``` r

fit <- balance(
  study,
  exposure,
  c(age, score),
  method = bw_entropy(),
  estimand = "ate"
)
study$w <- weights(fit)

outcome_mod <- glm(
  event ~ exposure,
  data = study,
  family = quasibinomial(),
  weights = w
)
```

The outcome model reports a standard error for the exposure coefficient,
but it is conditional on the weights. It treats the reweighting as given
and so leaves out one of the sources of variation in the estimate.

``` r

naive <- summary(outcome_mod)$coefficients["exposure", ]
naive
#>     Estimate   Std. Error      t value     Pr(>|t|) 
#> 6.472828e-01 1.047788e-01 6.177611e+00 8.368374e-10
```

The right standard error propagates the uncertainty from estimating the
weights into the effect estimate. Whether it is larger or smaller than
the naive value depends on the estimand and the design, but it is
generally not the same, and using the naive value can give intervals
with the wrong coverage.

## M-estimation with `ipw()`

The estimating-equation methods
([`bw_entropy()`](https://r-causal.github.io/balancing/reference/bw_entropy.md),
[`bw_ipt()`](https://r-causal.github.io/balancing/reference/bw_ipt.md),
and the just-identified
[`bw_cbps()`](https://r-causal.github.io/balancing/reference/bw_cbps.md))
fit weights that solve smooth estimating equations. That structure lets
the uncertainty be propagated analytically, without resampling.
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
writes one stacked estimating-equation system whose parameters are the
weight parameters, the outcome-model coefficients, the two marginal
means, and the effect contrasts, and the deli package differentiates
that system at the fitted values and returns its empirical sandwich
covariance. Nothing is re-solved: every parameter enters at the value
its own fit already found. Because each contrast is a parameter of the
stack rather than a transformation applied afterward, its standard error
comes straight off the diagonal of the joint covariance, with no
delta-method step in between. The result is a standard error that
accounts for having estimated the weights.

``` r

result <- ipw(fit, outcome_mod)
result
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: balance(.data = study, .exposure = exposure, .covariates = c(age, 
#>     score), method = bw_entropy(), estimand = "ate") 
#> 
#> Outcome Model:
#>   Call: glm(formula = event ~ exposure, family = quasibinomial(), data = study, 
#>     weights = w) 
#> 
#> Marginal estimates:
#>                estimate  std.err       z ci.lower ci.upper conf.level   p.value
#> mean 0         0.407009 0.019144 21.2600  0.36949  0.44453       0.95 < 2.2e-16
#> mean 1         0.567323 0.018809 30.1628  0.53046  0.60419       0.95 < 2.2e-16
#> rd 1 vs 0      0.160314 0.026838  5.9734  0.10771  0.21292       0.95 2.323e-09
#> log(rr) 1 vs 0 0.332094 0.057547  5.7709  0.21931  0.44488       0.95 7.886e-09
#> log(or) 1 vs 0 0.647283 0.110286  5.8691  0.43113  0.86344       0.95 4.381e-09
#>                   
#> mean 0         ***
#> mean 1         ***
#> rd 1 vs 0      ***
#> log(rr) 1 vs 0 ***
#> log(or) 1 vs 0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The `rd` row is the risk difference on the same scale as the outcome
model’s marginal effect. Comparing its standard error to the naive one
shows the correction.

``` r

rd <- as.data.frame(result$estimates)
rd <- rd[rd$effect == "rd", ]

data.frame(
  source = c("naive (weights fixed)", "ipw (weights estimated)"),
  std.err = c(naive[["Std. Error"]], rd$std.err)
)
#>                    source    std.err
#> 1   naive (weights fixed) 0.10477882
#> 2 ipw (weights estimated) 0.02683792
```

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
computes this variance for an estimating-equation fit at exact balance:
any of the three methods at a binary or categorical exposure, and
entropy balancing at a continuous one. The outcome model may be a
[`glm()`](https://rdrr.io/r/stats/glm.html) with a binomial,
quasibinomial, or gaussian family, or a plain
[`lm()`](https://rdrr.io/r/stats/lm.html). For a discrete exposure any
link those families carry is handled exactly, including a non-canonical
one such as probit, because the sandwich differentiates the score it is
given rather than assuming an information matrix. Model offsets are
supported natively, written either as an
[`offset()`](https://rdrr.io/r/stats/offset.html) term in the formula or
passed through the model’s `offset` argument, and are carried through
both the outcome-model score and the marginal means. A continuous
exposure reports coefficients rather than contrasts of predictions, so
an offset that reads the exposure moves them off the effects they name
and is refused there; an exposure-free offset, such as the person-time
offset of a rate model, is supported as it is everywhere else. That
refusal reads the offset expression, so an offset arriving as a
precomputed vector is beyond it and keeping the exposure out of one is
yours to honor.

### Categorical exposures

A categorical exposure works the same way, with one marginal mean per
level. The stacked system carries all of them, and the table reports
each of those means and then, for a plain categorical exposure, the
contrasts of each non-reference level against the reference level, which
is the first of the exposure’s own levels. A crossing of two treatments
declared with
[`causalgenerics::joint_exposure()`](https://r-causal.github.io/causalgenerics/reference/joint_exposure.html)
is reported in the treatments rather than in the cells, and the section
on joint exposures below describes those rows. The outcome model enters
the exposure as a factor.

``` r

odds_medium <- exp(0.6 * age - 0.3 * score)
odds_high <- exp(-0.4 * age + 0.7 * score)
denominator <- 1 + odds_medium + odds_high
draw <- runif(n)

study$arm <- factor(
  ifelse(
    draw < 1 / denominator,
    "low",
    ifelse(draw < (1 + odds_medium) / denominator, "medium", "high")
  ),
  levels = c("low", "medium", "high")
)
study$relapse <- rbinom(
  n,
  1,
  plogis(
    -0.4 +
      0.5 * (study$arm == "medium") +
      0.9 * (study$arm == "high") +
      0.4 * age -
      0.3 * score
  )
)

arm_fit <- balance(
  study,
  arm,
  c(age, score),
  method = bw_ipt(),
  estimand = "ate"
)
study$arm_w <- weights(arm_fit)

arm_mod <- glm(
  relapse ~ arm,
  data = study,
  family = quasibinomial(),
  weights = arm_w
)

ipw(arm_fit, arm_mod)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: balance(.data = study, .exposure = arm, .covariates = c(age, 
#>     score), method = bw_ipt(), estimand = "ate") 
#> 
#> Outcome Model:
#>   Call: glm(formula = relapse ~ arm, family = quasibinomial(), data = study, 
#>     weights = arm_w) 
#> 
#> Marginal estimates:
#>                       estimate  std.err       z ci.lower ci.upper conf.level
#> mean low              0.386610 0.022665 17.0577 0.342188  0.43103       0.95
#> mean medium           0.493063 0.026132 18.8685 0.441846  0.54428       0.95
#> mean high             0.606883 0.025079 24.1986 0.557729  0.65604       0.95
#> rd medium vs low      0.106453 0.034223  3.1105 0.039377  0.17353       0.95
#> log(rr) medium vs low 0.243221 0.078185  3.1108 0.089982  0.39646       0.95
#> log(or) medium vs low 0.433836 0.140134  3.0959 0.159179  0.70849       0.95
#> rd high vs low        0.220273 0.033521  6.5713 0.154574  0.28597       0.95
#> log(rr) high vs low   0.450921 0.071158  6.3369 0.311454  0.59039       0.95
#> log(or) high vs low   0.895815 0.140884  6.3585 0.619686  1.17194       0.95
#>                         p.value    
#> mean low              < 2.2e-16 ***
#> mean medium           < 2.2e-16 ***
#> mean high             < 2.2e-16 ***
#> rd medium vs low       0.001867 ** 
#> log(rr) medium vs low  0.001866 ** 
#> log(or) medium vs low  0.001962 ** 
#> rd high vs low        4.988e-11 ***
#> log(rr) high vs low   2.344e-10 ***
#> log(or) high vs low   2.037e-10 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The table leads with the marginal mean of each arm and then reports each
effect measure once per contrast, and the `contrast` column names what
each row belongs to: the level itself for a mean, and the two levels
compared for an effect. A binary exposure is read by the same rule, its
means named for the two levels and its single comparison named as a
contrast of them. Everything else carries over unchanged: the outcome
model may adjust for covariates, the marginal means standardize over the
estimand’s target population, and the standard errors account for having
estimated the weights.

### Effect modification with `.by`

`.by` names a modifier, and the result reports the rows it would report
without one, then the marginal means and effects within each level of
the modifier, then each non-reference subgroup against the reference
one. The table gains a `group` column naming the subgroup a row was
estimated in, placed after `contrast`, since a subgroup qualifies a
whole comparison rather than one side of it.

The data below carry a treatment whose effect is larger among men than
among women by construction.

``` r

set.seed(7)
m <- 500
severity <- rnorm(m)
comorbidity <- rnorm(m)
sex <- factor(
  sample(c("female", "male"), m, replace = TRUE),
  levels = c("female", "male")
)
treated <- rbinom(m, 1, plogis(0.5 * severity - 0.4 * comorbidity))
recovery <- rbinom(
  m,
  1,
  plogis(
    -0.4 +
      0.45 * treated +
      0.9 * treated * (sex == "male") +
      0.4 * severity -
      0.3 * comorbidity
  )
)

trial <- data.frame(treated, severity, comorbidity, sex, recovery)
```

The weights are fit once, over the whole sample, and the outcome model
interacts the exposure with the modifier so that the fitted response is
free to differ between the two groups.

``` r

trial_fit <- balance(
  trial,
  treated,
  c(severity, comorbidity),
  method = bw_entropy(),
  estimand = "ate"
)
trial$trial_w <- weights(trial_fit)

trial_mod <- glm(
  recovery ~ treated * sex,
  data = trial,
  family = quasibinomial(),
  weights = trial_w
)

ipw(trial_fit, trial_mod, .by = sex)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: balance(.data = trial, .exposure = treated, .covariates = c(severity, 
#>     comorbidity), method = bw_entropy(), estimand = "ate") 
#> 
#> Outcome Model:
#>   Call: glm(formula = recovery ~ treated * sex, family = quasibinomial(), 
#>     data = trial, weights = trial_w) 
#> 
#> Marginal estimates:
#>                                           estimate  std.err       z ci.lower
#> mean 0 overall                            0.317188 0.029794 10.6460 0.258793
#> mean 1 overall                            0.607507 0.032381 18.7611 0.544041
#> rd 1 vs 0 overall                         0.290319 0.044181  6.5711 0.203725
#> log(rr) 1 vs 0 overall                    0.649869 0.108379  5.9963 0.437449
#> log(or) 1 vs 0 overall                    1.203569 0.194093  6.2010 0.823154
#> mean 0 sex = female                       0.349733 0.042586  8.2125 0.266267
#> mean 1 sex = female                       0.509935 0.046754 10.9068 0.418299
#> mean 0 sex = male                         0.283042 0.042281  6.6944 0.200173
#> mean 1 sex = male                         0.709877 0.045528 15.5921 0.620644
#> rd 1 vs 0 sex = female                    0.160201 0.063241  2.5332 0.036251
#> log(rr) 1 vs 0 sex = female               0.377111 0.152424  2.4741 0.078365
#> rd 1 vs 0 sex = male                      0.426836 0.062132  6.8698 0.305058
#> log(rr) 1 vs 0 sex = male                 0.919497 0.162565  5.6562 0.600875
#> rd 1 vs 0 sex = male vs sex = female      0.266634 0.089849  2.9676 0.090534
#> log(rr) 1 vs 0 sex = male vs sex = female 0.542386 0.225250  2.4079 0.100904
#>                                           ci.upper conf.level   p.value    
#> mean 0 overall                             0.37558       0.95 < 2.2e-16 ***
#> mean 1 overall                             0.67097       0.95 < 2.2e-16 ***
#> rd 1 vs 0 overall                          0.37691       0.95 4.995e-11 ***
#> log(rr) 1 vs 0 overall                     0.86229       0.95 2.019e-09 ***
#> log(or) 1 vs 0 overall                     1.58398       0.95 5.611e-10 ***
#> mean 0 sex = female                        0.43320       0.95 < 2.2e-16 ***
#> mean 1 sex = female                        0.60157       0.95 < 2.2e-16 ***
#> mean 0 sex = male                          0.36591       0.95 2.166e-11 ***
#> mean 1 sex = male                          0.79911       0.95 < 2.2e-16 ***
#> rd 1 vs 0 sex = female                     0.28415       0.95  0.011303 *  
#> log(rr) 1 vs 0 sex = female                0.67586       0.95  0.013358 *  
#> rd 1 vs 0 sex = male                       0.54861       0.95 6.431e-12 ***
#> log(rr) 1 vs 0 sex = male                  1.23812       0.95 1.548e-08 ***
#> rd 1 vs 0 sex = male vs sex = female       0.44273       0.95  0.003001 ** 
#> log(rr) 1 vs 0 sex = male vs sex = female  0.98387       0.95  0.016043 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The `overall` block is the result an ungrouped call reports, and the
blocks after it are the effect within each level of `sex`. A subgroup’s
marginal means are the g-computation means over that subgroup’s units
alone, and its risk difference is the contrast of those two means.
Nothing is refit inside a subgroup: the weights are the ones the single
fit produced, and no subgroup gets balance conditions of its own. The
last two rows compare the subgroups, and since every row here is a
parameter of the same stacked system, that comparison is a parameter as
well rather than a difference taken between two separately fitted
effects. Its standard error therefore carries the covariance between the
two subgroups, which share the weight parameters and the outcome model’s
coefficients; adding two independently estimated variances would miss
that covariance, in a direction the design decides.

The log odds ratio stays among the whole-sample rows. It is
noncollapsible, so the odds ratio over a sample is not an average of the
odds ratios within its subgroups, and a difference of two of them moves
with each subgroup’s baseline risk whether or not the effect is modified
at all. The collapsible measures, `rd` and `log(rr)`, are the ones
reported for the subgroups and for their contrast.

The subgroup contrast is not what the interaction term in
`recovery ~ treated * sex` reports. That coefficient belongs to the
outcome model’s conditional surface: it is the change in the model’s log
odds ratio between the two levels of `sex`, on the link scale the model
was fitted with. The `.by` rows are an estimand instead, the marginal
effect of treating standardized over each subgroup’s own covariate
distribution, and their contrast is the difference between two such
effects on each collapsible scale. The two are different quantities and
generally different numbers, so a term in the model is not a reading of
the request and the request is not a reading of the term. What `.by`
does need from the model is room for the effect to differ. The subgroup
rows are g-computation on the model as it was specified, so a model
carrying no term that reads both the exposure and the modifier is warned
about rather than refused: the modification may reach it through a
column derived from the modifier instead.

### Joint exposures

Intervening on two treatments at once is a third question again, and it
is not the one `.by` answers.
[`causalgenerics::joint_exposure()`](https://r-causal.github.io/causalgenerics/reference/joint_exposure.html)
crosses two discrete treatments into one categorical exposure and
records the crossing on the column.
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
weights that column as it weights any factor over the same cells, so
what the declaration changes is which effects are reported rather than
which population is balanced.

``` r

set.seed(21)
k <- 500
frailty <- rnorm(k)
drug <- rbinom(k, 1, plogis(0.4 * frailty))
diet <- rbinom(k, 1, plogis(-0.2 + 0.5 * frailty + 0.6 * drug))
remission <- rbinom(
  k,
  1,
  plogis(-1.4 + 0.3 * drug + 0.3 * diet + 1.6 * drug * diet + 0.5 * frailty)
)

regimens <- data.frame(
  frailty = frailty,
  remission = remission,
  drug = factor(drug),
  diet = factor(diet)
)
# Built after the frame so the crossing reads the same two columns the frame
# carries.
regimens$regimen <- causalgenerics::joint_exposure(
  drug = regimens$drug,
  diet = regimens$diet
)
```

``` r

regimen_fit <- balance(
  regimens,
  regimen,
  c(frailty),
  method = bw_ipt(),
  estimand = "ate"
)
regimens$regimen_w <- weights(regimen_fit)

regimen_mod <- glm(
  remission ~ regimen + frailty,
  data = regimens,
  family = quasibinomial(),
  weights = regimen_w
)

ipw(regimen_fit, regimen_mod)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: balance(.data = regimens, .exposure = regimen, .covariates = c(frailty), 
#>     method = bw_ipt(), estimand = "ate") 
#> 
#> Outcome Model:
#>   Call: glm(formula = remission ~ regimen + frailty, family = quasibinomial(), 
#>     data = regimens, weights = regimen_w) 
#> 
#> Marginal estimates:
#>                                           estimate  std.err       z   ci.lower
#> mean drug = 0, diet = 0 overall           0.183250 0.036263  5.0534  0.1121765
#> mean drug = 1, diet = 0 overall           0.297125 0.043441  6.8397  0.2119823
#> mean drug = 0, diet = 1 overall           0.236567 0.037595  6.2926  0.1628826
#> mean drug = 1, diet = 1 overall           0.643617 0.037633 17.1023  0.5698571
#> rd drug: 1 vs 0 diet = 0                  0.113875 0.056162  2.0276  0.0037984
#> log(rr) drug: 1 vs 0 diet = 0             0.483300 0.244244  1.9788  0.0045901
#> rd drug: 1 vs 0 diet = 1                  0.407050 0.052230  7.7934  0.3046814
#> log(rr) drug: 1 vs 0 diet = 1             1.000873 0.167350  5.9807  0.6728723
#> rd diet: 1 vs 0 drug = 0                  0.053317 0.051672  1.0318 -0.0479579
#> log(rr) diet: 1 vs 0 drug = 0             0.255378 0.251131  1.0169 -0.2368303
#> rd diet: 1 vs 0 drug = 1                  0.346492 0.056676  6.1135  0.2354089
#> log(rr) diet: 1 vs 0 drug = 1             0.772951 0.155941  4.9567  0.4673111
#> rd drug: 1 vs 0 diet = 1 vs diet = 0      0.293176 0.076862  3.8143  0.1425290
#> log(rr) drug: 1 vs 0 diet = 1 vs diet = 0 0.517573 0.296887  1.7433 -0.0643152
#>                                           ci.upper conf.level   p.value    
#> mean drug = 0, diet = 0 overall            0.25432       0.95 4.340e-07 ***
#> mean drug = 1, diet = 0 overall            0.38227       0.95 7.934e-12 ***
#> mean drug = 0, diet = 1 overall            0.31025       0.95 3.123e-10 ***
#> mean drug = 1, diet = 1 overall            0.71738       0.95 < 2.2e-16 ***
#> rd drug: 1 vs 0 diet = 0                   0.22395       0.95 0.0426015 *  
#> log(rr) drug: 1 vs 0 diet = 0              0.96201       0.95 0.0478434 *  
#> rd drug: 1 vs 0 diet = 1                   0.50942       0.95 6.522e-15 ***
#> log(rr) drug: 1 vs 0 diet = 1              1.32887       0.95 2.222e-09 ***
#> rd diet: 1 vs 0 drug = 0                   0.15459       0.95 0.3021495    
#> log(rr) diet: 1 vs 0 drug = 0              0.74759       0.95 0.3091965    
#> rd diet: 1 vs 0 drug = 1                   0.45758       0.95 9.744e-10 ***
#> log(rr) diet: 1 vs 0 drug = 1              1.07859       0.95 7.171e-07 ***
#> rd drug: 1 vs 0 diet = 1 vs diet = 0       0.44382       0.95 0.0001366 ***
#> log(rr) drug: 1 vs 0 diet = 1 vs diet = 0  1.09946       0.95 0.0812756 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The table leads with the counterfactual mean of each cell, one per
combination of the two treatments, each standardized over the whole
sample. The simple effects follow: each treatment against its own
reference level, within a level of the other that the `group` column
names. The last pair of rows is the interaction, the difference between
the drug’s effect at `diet = 1` and its effect at `diet = 0`, which is
large here because each treatment does much less on its own than the two
do together. Reported as a plain categorical exposure, these four cells
would give each cell against the reference cell, under contrasts like
`"drug = 1, diet = 0 vs drug = 0, diet = 0"`; the declaration replaces
those rows rather than adding to them, so no row here is named that way.
The interaction appears once, under the first treatment’s framing,
because it is symmetric in the two: the change in the drug’s effect
across the levels of the diet is the change in the diet’s effect across
the levels of the drug, and reporting it twice would put one quantity in
the table under two names. No contrast row carries a log odds ratio, for
the reason no subgroup row does.

The declaration is read off the exposure column of the frame
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
resolves, which is the outcome model’s own frame here, so dropping it
with [`factor()`](https://rdrr.io/r/base/factor.html) returns the
cell-against-cell rows. A declared crossing is reported for the ATE
alone, and `.by` on one is refused: effect modification of a joint
intervention is a three-way question, and this surface reports neither
it nor a projection of it.
[`?ipw.balancing`](https://r-causal.github.io/balancing/reference/ipw.balancing.md)
gives the full row set and both restrictions.

### Continuous exposures

A continuous exposure has no levels to contrast, so there is no pair of
marginal means to difference.
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
reports the dose response of a weighted marginal structural model
instead: entropy balancing removes the association between the exposure
and the covariates, and every coefficient of the outcome model that
reads the exposure is an effect on the model’s own link scale.

``` r

study$dose <- 0.8 * age - 0.5 * score + rnorm(n)
study$response <- 1 + 0.4 * study$dose + 0.5 * age - 0.3 * score + rnorm(n)

dose_fit <- balance(
  study,
  dose,
  c(age, score),
  method = bw_entropy(),
  estimand = "ate"
)
study$dose_w <- weights(dose_fit)

dose_mod <- lm(response ~ dose, data = study, weights = dose_w)

ipw(dose_fit, dose_mod)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: balance(.data = study, .exposure = dose, .covariates = c(age, 
#>     score), method = bw_entropy(), estimand = "ate") 
#> 
#> Outcome Model:
#>   Call: lm(formula = response ~ dose, data = study, weights = dose_w) 
#> 
#> Marginal estimates:
#>       estimate std.err      z ci.lower ci.upper conf.level   p.value    
#> slope  0.39271 0.03888 10.101  0.31651  0.46892       0.95 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

An exposure entering through one design column is the whole of the dose
response, so its coefficient is that response’s slope everywhere. The
table holds one row, named for the link: `slope` for an identity link,
`log(or)` for a logit, and `log(rr)` for a log link.

The exposure may also enter through several columns. What the model has
to keep is variable membership: every term reading the exposure must
read the exposure alone, however many columns it expands to, so a curve
written out term by term and one handed to a basis constructor are both
reported.

``` r

curve_mod <- lm(response ~ splines::ns(dose, 3), data = study, weights = dose_w)

ipw(dose_fit, curve_mod)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: conditional (outcome model) 
#> 
#> Weight Estimator:
#>   Call: balance(.data = study, .exposure = dose, .covariates = c(age, 
#>     score), method = bw_entropy(), estimand = "ate") 
#> 
#> Outcome Model:
#>   Call: lm(formula = response ~ splines::ns(dose, 3), data = study, weights = dose_w) 
#> 
#> Conditional estimates (outcome model):
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept)           -0.86446    0.84010 -1.0290  0.303480    
#> splines::ns(dose, 3)1  1.99573    0.37897  5.2662 1.393e-07 ***
#> splines::ns(dose, 3)2  4.24376    1.84443  2.3009  0.021400 *  
#> splines::ns(dose, 3)3  2.79960    0.89688  3.1215  0.001799 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The table then holds one row per coefficient and gains a `contrast`
column naming each row after the coefficient the fit names. The scale
word steps back to `coef` at an identity link, since a curve has a
different slope at every dose and no one of its coefficients is that
slope; a logit still reports `log(or)` and a log link `log(rr)`, because
a coefficient of those models is a log ratio whatever column it
multiplies. Nothing is standardized on this path either way, so each row
is a coefficient of the weighted fit and what the stack adds is the
standard error.

A term reading a covariate alongside the exposure is refused.
`response ~ dose * age` contributes a coefficient that is a change in
the dose response per unit of `age`, so there is no one effect for a row
to report and no value of `age` a row could name it at. Only entropy
balancing at exact balance reaches this path, and only for the ATE,
which is the only estimand a continuous fit targets.

The continuous standard error reaches its nominal coverage more slowly
than the binary risk difference does. Over 500 draws, its ratio of mean
standard error to the standard deviation of the estimates was 0.836 at
300 observations and 0.942 at 1200, so an interval at a few hundred
observations is too narrow. Bootstrap it at that size.

### Covariate-adjusted outcome models

The outcome model must carry the exposure among its predictors, but it
may adjust for covariates alongside it, and may interact them with the
exposure. Adjusting for the same covariates the weights balance is a
common way to reduce residual confounding and tighten the estimate.

``` r

adjusted_mod <- glm(
  event ~ exposure + age + score,
  data = study,
  family = quasibinomial(),
  weights = w
)

ipw(fit, adjusted_mod)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: balance(.data = study, .exposure = exposure, .covariates = c(age, 
#>     score), method = bw_entropy(), estimand = "ate") 
#> 
#> Outcome Model:
#>   Call: glm(formula = event ~ exposure + age + score, family = quasibinomial(), 
#>     data = study, weights = w) 
#> 
#> Marginal estimates:
#>                estimate  std.err       z ci.lower ci.upper conf.level   p.value
#> mean 0         0.407496 0.019320 21.0920  0.36963  0.44536       0.95 < 2.2e-16
#> mean 1         0.566987 0.018973 29.8846  0.52980  0.60417       0.95 < 2.2e-16
#> rd 1 vs 0      0.159491 0.026781  5.9555  0.10700  0.21198       0.95 2.593e-09
#> log(rr) 1 vs 0 0.330306 0.057430  5.7515  0.21775  0.44287       0.95 8.848e-09
#> log(or) 1 vs 0 0.643897 0.110020  5.8525  0.42826  0.85953       0.95 4.842e-09
#>                   
#> mean 0         ***
#> mean 1         ***
#> rd 1 vs 0      ***
#> log(rr) 1 vs 0 ***
#> log(or) 1 vs 0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Adjustment changes what the marginal means average over. A marginal
model predicts one value per exposure level, so it makes no difference
which units you average across; an adjusted model predicts a value per
unit, and the population you average across becomes part of the
estimand.
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
standardizes over the estimand’s target population: every unit for the
ATE, and the focal group’s units for the ATT or the ATC. Sampling
weights, where the fit has them, weight that average as well. The
standard error accounts for the adjustment the same way it accounts for
everything else, by carrying the adjusted model’s own score equations in
the stack.

A model that leaves the exposure out has nothing to contrast, so it is
refused rather than reported as a null effect.

``` r

ipw(fit, glm(event ~ age + score, data = study, family = quasibinomial(), weights = w))
#> Error in `ipw()`:
#> ! `outcome_mod` must include the exposure among its predictors.
#> ✖ The exposure "exposure" appears in none of its terms.
#> ℹ The model may adjust for covariates alongside the exposure, and may carry the
#>   exposure inside a transformation such as `factor()`.
```

Everything else requires the bootstrap. The quadratic-program methods
([`bw_energy()`](https://r-causal.github.io/balancing/reference/bw_energy.md),
[`bw_cfd()`](https://r-causal.github.io/balancing/reference/bw_cfd.md),
and
[`bw_sbw()`](https://r-causal.github.io/balancing/reference/bw_sbw.md))
and an entropy fit at a positive tolerance solve inequality-constrained
problems, whose solutions are characterized by which constraints bind at
the optimum rather than by a smooth system set to zero. An
over-identified
[`bw_cbps()`](https://r-causal.github.io/balancing/reference/bw_cbps.md)
fit minimizes a generalized-method-of-moments criterion in more moment
conditions than it has parameters, so its balancing conditions are not
solved to zero either. None of these has an estimating-equation
representation to stack, so the fit carries none and
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
raises an error pointing here. Inverse probability tilting is the
exception to the tolerance rule: it solves the same smooth estimating
equations whatever tolerance is requested, so
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
handles it at any tolerance.

``` r

sbw_fit <- balance(
  study,
  exposure,
  c(age, score),
  method = bw_sbw(),
  constraints = balance_terms(tolerance = 0.02)
)
ipw(sbw_fit, outcome_mod)
#> Error in `ipw()`:
#> ! `ipw()` cannot compute a stacked variance for this balancing fit.
#> ✖ This fit's weights do not solve smooth estimating equations, so the stacked
#>   variance is unavailable.
#> ℹ Estimating equations come from the estimating-equation family (entropy
#>   balancing, inverse probability tilting, just-identified covariate balancing
#>   propensity score) with exact balance.
#> ℹ See the inference vignette for a bootstrap workflow.
```

## Bootstrapping the other fits

Where a fit carries no estimating equations, inference uses resampling
instead. The recipe is a standard nonparametric bootstrap: resample the
rows with replacement, refit the weights and the outcome model on each
resample, and take the spread of the effect estimates across resamples
as the standard error. The important detail is that the reweighting is
refit inside the loop, so every source of variation, including the
estimation of the weights, is captured. The example below uses
[`bw_sbw()`](https://r-causal.github.io/balancing/reference/bw_sbw.md),
but the recipe is the same for any fit and it also covers the cases
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
declines for reasons other than the weights, such as an outcome model
outside the supported families.

We’ll show a hand-rolled approach, but there are several tools for
bootstrapping in R. See in particular the `{boot}` package and the
[rsample](https://rsample.tidymodels.org) package. See also the appendix
on the bootstrap in [Causal Inference in
R](https://www.r-causal.org/appendices/a-bootstrap).

We keep the example small so it runs quickly: a few hundred resamples on
a modest sample.

``` r

boot_study <- study[seq_len(400), ]

fit_effect <- function(data) {
  fit <- balance(
    data,
    exposure,
    c(age, score),
    method = bw_sbw(),
    constraints = balance_terms(tolerance = 0.02)
  )
  data$w <- weights(fit)
  model <- glm(
    event ~ exposure,
    data = data,
    family = quasibinomial(),
    weights = w
  )
  # The marginal risk difference by g-computation
  p1 <- predict(model, transform(data, exposure = 1), type = "response")
  p0 <- predict(model, transform(data, exposure = 0), type = "response")
  mean(p1) - mean(p0)
}

set.seed(1)
point_estimate <- fit_effect(boot_study)

n_boot <- 200
boot_estimates <- vapply(seq_len(n_boot), function(b) {
  rows <- sample(nrow(boot_study), replace = TRUE)
  fit_effect(boot_study[rows, ])
}, numeric(1))
```

The bootstrap standard error is the standard deviation of the resampled
estimates, and a percentile interval comes from their quantiles.

``` r

data.frame(
  estimate = point_estimate,
  std.err = sd(boot_estimates),
  ci.lower = quantile(boot_estimates, 0.025),
  ci.upper = quantile(boot_estimates, 0.975),
  row.names = NULL
)
#>    estimate    std.err   ci.lower  ci.upper
#> 1 0.1446761 0.05401617 0.03356418 0.2587607
```

A real analysis would use more resamples, typically at least one or two
thousand, and might use the bias-corrected and accelerated interval
rather than the percentile one. The structure stays the same: refit the
weights inside every resample.

## Choosing an approach

Use
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
when the method carries estimating equations, which is the
estimating-equation family with a binary or categorical exposure and
exact balance, and entropy balancing with exact balance at a continuous
exposure. It is faster and needs no tuning. Use the bootstrap for the
quadratic-program methods, for an entropy fit at a positive tolerance,
for an over-identified
[`bw_cbps()`](https://r-causal.github.io/balancing/reference/bw_cbps.md)
fit, and for a continuous fit from a method that solves no estimating
equations. A continuous exposure is worth bootstrapping at a small
sample size in any case: its stacked standard error is a large-sample
one and runs anticonservative at a few hundred observations.
Covariate-adjusted outcome models need neither route in particular:
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
stacks them directly. The bootstrap is more general but more expensive,
and its precision improves with the number of resamples.
