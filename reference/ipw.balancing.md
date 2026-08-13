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
  exposure it is a marginal structural model whose every
  exposure-reading term reads the exposure alone.

- .data:

  The data frame holding the exposure and outcome. If `NULL`, the values
  are taken from the outcome model frame. It carries the exposure column
  the fixed-exposure predictions are built from, so it has nothing to
  supply for a continuous exposure, which makes no such predictions, and
  is ignored there.

  It must hold the rows both models were fitted on, in the same order.
  Half of the stacked system reads it and half reads the fit: the
  counterfactual predictions, the subgroup indicators, and a focal
  estimand's standardization come from `.data`, while the weight
  equations and the outcome-model score stay in the fit's order. A frame
  holding the right rows in another order therefore leaves every effect
  estimate unchanged, since a weighted mean does not care in which order
  it is summed, and makes every standard error wrong. Each column
  `.data` and the outcome model frame name in common is compared value
  by value, and a disagreement raises `balancing_ipw_input_error` naming
  the column and the first row it disagrees at. Columns the outcome
  model never saw, the modifier `.by` names among them, are free.

- estimand:

  The causal estimand. If `NULL`, the fit's estimand is used. As in
  [`balance()`](https://r-causal.github.io/balancing/reference/balance.md),
  `"atc"` is accepted as a synonym for `"atu"`. Supplying an estimand
  that disagrees with the fit raises `balancing_estimand_error`.

- conf_level:

  The confidence level for the intervals. Default `0.95`.

- effects:

  The presentation mode the result records, either `"marginal"` (the
  default) or `"conditional"`. The marginal reading reports the
  population-averaged causal contrasts described above; the conditional
  reading reports the outcome model's coefficient surface. Both surfaces
  are computed whichever mode is named, since the stacked system is
  solved either way, so the argument settles which one the result
  presents and nothing else.
  [`causalgenerics::as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
  and
  [`causalgenerics::as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
  move a result between the two readings afterwards, a pooled result as
  much as an unpooled one, and the accessors take an `effects` argument
  of their own for a single call.

  The conditional reading reports the coefficients of the stored
  `outcome_mod` against the outcome block of the stacked sandwich, which
  is the block the wrapper around that model carries, so its standard
  errors account for having estimated the weights rather than treating
  them as fixed. The stored `wt_mod` is the balancing fit itself rather
  than a wrapper of one, since the weight block of the same stack is
  already its own: the stored copy carries that block, the covariance of
  the fit's weight parameters, and
  [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html) on it reads the
  block back. The fit that went into the call is untouched, and reports
  no covariance of its own.

- .by:

  A modifier to report the effects within the levels of, given unquoted
  and selected with tidyselect out of `.data` where one was supplied and
  out of the outcome model frame otherwise. The default, `NULL`, is the
  absence of a request: the result then reports the whole-sample effects
  alone and names no subgroups at all. A selection reaching any number
  of columns other than one is refused, since the effects are reported
  within the levels of a single variable. The effect modification
  section below describes the rows a request adds and the configurations
  it refuses, and the joint exposure section describes why a declared
  crossing takes none.

- ...:

  Ignored, for compatibility with the generic.

## Value

An object of class `ipw`, an implementation of
[`causalgenerics::ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html).
Alongside `estimand`, `wt_mod`, `outcome_mod`, the `estimates` table,
and the `effects` field recording the presentation mode described above,
the result carries two fields describing the variance:

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
  carries neither block, since its effects are outcome-model
  coefficients; each of those coefficients is named for the label its
  estimates row carries, as `slope` or `coef I(exposure^2)`, in place of
  the `beta_` name the others keep. A `.by` request appends, after all
  of those, every subgroup's mean block in subgroup order, then every
  subgroup's contrast block in that same order, then one contrast block
  per non-reference subgroup against the reference one. Each of their
  names is the name of the block it repeats, suffixed with its group, as
  `mu0_sex = female` and `rd_sex = female vs sex = male`. A declared
  joint exposure keeps the mean block and replaces the contrast block,
  naming each of its own contrasts for the measure and the row, as
  `rd_a: 1 vs 0 e = 0`. The standard errors in `estimates` are
  `sqrt(diag(fit$vcov))` read at those effect names.

The `estimates` table carries the covariance of the reported effects as
its `ipw_vcov` attribute, which is what
[`stats::vcov()`](https://rdrr.io/r/stats/vcov.html) returns in the
marginal reading. Both its dimnames are the display labels described
above, as `"rd b vs a sex = female"`. The stored `outcome_mod` is
wrapped by
[`causalgenerics::new_ipw_model()`](https://r-causal.github.io/causalgenerics/reference/new_ipw_model.html),
which carries the outcome-model block of `fit$vcov` under the model's
own coefficient names, so [`vcov()`](https://rdrr.io/r/stats/vcov.html)
on it reports the joint-estimation variance. The stored `wt_mod` carries
the leading block of the same covariance under the `theta_w` names
above, which name the fit's own parameters on either route, so
[`vcov()`](https://rdrr.io/r/stats/vcov.html) on it reports the
covariance of the weight parameters.
[`stats::df.residual()`](https://rdrr.io/r/stats/df.residual.html)
returns `NA_integer_`, since the stacked system is not a fit with
residual degrees of freedom of its own.

[`stats::nobs()`](https://rdrr.io/r/stats/nobs.html) delegates to the
stored outcome model, which counts the rows it was fitted on that carry
a nonzero weight. A unit given no sampling weight is pinned at zero
rather than dropped, so the weight vector stays the length of the data
the fit saw while the outcome model counts one row fewer for each pinned
unit, and [`nobs()`](https://rdrr.io/r/stats/nobs.html) on the result is
then smaller than `length(weights(result))`.

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
level, and the estimates table gains a `contrast` column, placed after
`effect`, naming each contrast as `"<level> vs <reference>"`. A binary
exposure keeps the table it has always returned, with no `contrast`
column. A categorical exposure declared as a crossing of two treatments
by
[`causalgenerics::joint_exposure()`](https://r-causal.github.io/causalgenerics/reference/joint_exposure.html)
replaces those level-against-reference rows with the surface described
under Joint exposures below.

A continuous exposure has no levels to contrast, so there is no pair of
marginal means to difference. What the method reports instead is the
dose response of a weighted marginal structural model: the balancing
weights break the exposure-covariate association, and every coefficient
of the outcome model that reads the exposure is an effect on the model's
own link scale. The estimates table holds one row per such coefficient,
read straight off the weighted fit, since nothing is standardized here.

An exposure entering through one design column, whether as a bare term
or as a transformation of one, is the whole of the dose response, so its
coefficient is that response's slope everywhere. Such a model keeps the
single-row table it has always returned, with the columns every other
exposure's table holds and no `contrast` column, and the row is named
for the link: `slope` for an identity link, whether the model arrives as
a [`stats::lm()`](https://rdrr.io/r/stats/lm.html) or as a gaussian
[`stats::glm()`](https://rdrr.io/r/stats/glm.html); `log(or)` for a
logit; and `log(rr)` for a log link.

An exposure entering through several columns, as in
`y ~ exposure + I(exposure^2)` or `y ~ splines::ns(exposure, 3)`, has no
such row, since a curve has a different slope at every dose. The table
then holds one row per column, gains a `contrast` column naming each row
after the coefficient the fit names, and the scale word steps back to
`coef` at an identity link. A logit still reports `log(or)` and a log
link `log(rr)`, because a coefficient of those models is a log ratio
whatever column it multiplies. Another link raises
`balancing_ipw_input_error` either way, since its coefficients are none
of those things. A continuous fit targets the average treatment effect
and nothing else, so an `estimand` supplied alongside it either agrees
or raises `balancing_estimand_error`, as it does for any other fit.

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
marginal structural model's coefficients alone, and each reported effect
is already one of those coefficients rather than a parameter derived
from them. A basis widens that stack by its own columns and changes
nothing else.

Nothing is re-solved along the way. Every parameter enters at the value
its own fit already found, and the stacked estimating functions are only
re-evaluated around that point. Because each contrast is a parameter of
the stack rather than a transformation applied afterward, its standard
error is already on the diagonal of the joint covariance and no
delta-method step stands between the sandwich and the reported effects.
The joint covariance propagates the uncertainty from estimating the
weights into the effect standard errors, which a variance that treats
the weights as fixed would understate.

The result reads through the accessors causalgenerics registers on the
class, which every fitting package's results share:
[`stats::coef()`](https://rdrr.io/r/stats/coef.html) for the reported
effects under their display labels,
[`stats::vcov()`](https://rdrr.io/r/stats/vcov.html) for their
covariance, [`stats::confint()`](https://rdrr.io/r/stats/confint.html)
for their intervals,
[`stats::nobs()`](https://rdrr.io/r/stats/nobs.html) for the number of
observations, and
[`stats::weights()`](https://rdrr.io/r/stats/weights.html) for the
weights the outcome model was fitted with. The covariance those
accessors read is the block of the stacked one belonging to the reported
effects, which is on the estimates table rather than derived from it:
the effect measures are transformations of the same pair of marginal
means, so they covary, and the off-diagonals a caller combining two of
them needs are not recoverable from the standard errors alone. The
stored outcome model carries its own block of the same covariance, so
[`vcov()`](https://rdrr.io/r/stats/vcov.html) on it reports the variance
of its coefficients with the uncertainty from estimating the weights
included, where a bare refit of the same weighted model treats the
weights as fixed and understates it.

A row's display label is built from the identity columns of the
estimates table, in the order the table carries them: the effect
measure, then the contrast where the surface names one, then the group
where it names one. A contrast names a categorical exposure's pair of
levels, a continuous surface's basis coefficient, or a joint exposure's
treatment; a group names the subgroup a `.by` request reports the row
within, or the level a joint exposure holds the other treatment at. That
one label names the row wherever a caller reads it, so the printed
table, the names of
[`stats::coef()`](https://rdrr.io/r/stats/coef.html), the dimnames of
[`stats::vcov()`](https://rdrr.io/r/stats/vcov.html), and the rownames
of [`stats::confint()`](https://rdrr.io/r/stats/confint.html) agree row
for row with the estimates table. The stacked `fit$theta` and `fit$vcov`
carry names of their own, which name blocks of the estimating-equation
system rather than reported rows and are described under Value.

The sandwich covariance is a large-sample one, and how large a sample it
takes differs by exposure. The binary risk-difference standard error is
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

A continuous exposure narrows that by variable membership, since what it
reports are coefficients of the model rather than contrasts of
predictions from it. Every term reading the exposure must read the
exposure alone, however many design columns it expands to, so
`y ~ exposure`, `y ~ exposure + x1`, `y ~ exposure + I(exposure^2)`,
`y ~ sin(exposure)`, `y ~ poly(exposure, 2)` and
`y ~ splines::ns(exposure, 3)` are all supported. A term reading a
covariate alongside the exposure raises `balancing_ipw_input_error`, so
`y ~ exposure * x1`, `y ~ exposure + exposure:x1` and
`y ~ I(exposure * x1)` are refused: each contributes a coefficient that
is a change in the dose response per unit of that covariate, so there is
no one effect for a row to report and no covariate value a row could
name it at. The check reads the model's terms rather than the text of
its formula, so the mixing is caught however it is written.

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

## Effect modification

`.by` names a modifier, and a result carrying one reports the effects it
reports without a request, then those same effects within each of the
modifier's levels, then each non-reference level against the reference
one. The estimates table gains a `group` column naming the subgroup each
row was estimated in, placed after `contrast` where a categorical
exposure names one and after `effect` where it does not, since a
subgroup qualifies the whole comparison rather than one side of it. The
whole-sample rows come first, under the group `"overall"`; a subgroup's
rows are named `"var = value"`, as `"sex = female"`; and a contrast of
subgroups joins the two, as `"sex = female vs sex = male"`. The
reference subgroup is the modifier's first level, which for a factor is
the first of its declared levels rather than the first in sorted order.
A character modifier declares no levels, so it is read as a factor on
the way in and its reference subgroup is its alphabetically first value,
whatever order the values appear in. Declare the column a factor to
measure the contrasts against some other subgroup.

A subgroup reports the collapsible measures alone. For a binary outcome
that is `rd` and `log(rr)`; a continuous outcome reports `diff`, the
only measure it has. The log odds ratio stays among the whole-sample
rows. An odds ratio is noncollapsible, so the odds ratio over a sample
is not an average of the odds ratios within its subgroups and the
difference of two of them is not the difference in effect it reads as;
nothing in the whole-sample rows averages anything over subgroups, which
is why they keep it. A categorical exposure crosses its contrasts with
the subgroups, so each subgroup block reports each contrast in the order
the whole-sample block reports it.

A subgroup's marginal means are the g-computation means over that
subgroup alone, standardized the way the whole-sample means are: over
every unit of the subgroup for a pooled estimand, and over its focal
units for `"att"` or `"atc"`. Sampling weights weight that average as
well.

Every row a request adds is a parameter of the same stacked system the
whole-sample rows are parameters of. The blocks are appended after
everything the ungrouped stack carries and nothing earlier reads a
parameter of theirs, so the leading blocks are the ones an ungrouped fit
produces and the whole-sample rows of a grouped result are the rows it
reported. What that buys is the covariance: the subgroups share the
weight parameters and the outcome model's coefficients, so their effects
covary, and each contrast of subgroups is a parameter of the joint
system rather than a difference of two separate fits, whose variance
would have to be read as a sum.

A level of the modifier that no unit carries names an empty subgroup and
is dropped rather than refused, which is what lets a modifier be subset
without being recoded first. Four configurations are refused with
`balancing_ipw_input_error`: a selection naming any number of columns
other than one, a modifier that is neither a factor nor a character
column, a modifier carrying missing values, and a modifier one of whose
subgroups holds fewer than all of the exposure levels. That last one has
no contrast to report in the subgroup that is short a level, and
refitting either model does not supply a comparison the data do not
hold. `.by` with a continuous exposure raises
`balancing_ipw_unsupported_error`: what such a fit reports is the
marginal structural model's own exposure coefficients rather than
contrasts of standardized means, so there is no effect within a subgroup
for the argument to name, and `.by` with a declared joint exposure
raises it for the reason given under Joint exposures below. A modifier
that reaches the argument through `.data` carries that frame's row-order
requirement with it, described under `.data` above: the subgroups are
built from the rows it holds while the weights stay in the fit's order.

An outcome model with no term reading both the exposure and the modifier
raises `balancing_ipw_by_interaction_warning` and the result is still
built. The subgroup effects are g-computation on the model as it was
specified, so the effect differs across subgroups only where a term
reads both columns, and that may be the modeling choice a caller meant
to make. The warning reports which terms were read rather than
announcing that the effect is the same in every subgroup, because a
model may carry the modification through a column derived from the
modifier instead: `y ~ exposure * sex_female` reported by `.by = sex`
names no term reading `sex`, and its subgroup effects differ all the
same.

## Joint exposures

[`causalgenerics::joint_exposure()`](https://r-causal.github.io/causalgenerics/reference/joint_exposure.html)
crosses two discrete treatments into one categorical exposure and
records the crossing on the vector. The result is a factor, so
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
weights it as it weights any factor over those cells, and the crossing
changes which effects are reported rather than which population is
balanced.

Reported as a plain categorical exposure, such a fit gives each cell
against the reference cell, under contrasts like
`"a = 1, e = 0 vs a = 0, e = 0"`. Those rows are arithmetically right
and they answer a question nobody asked. A declared crossing is reported
in the two treatments instead, and the cell-against-cell rows are
replaced rather than supplemented, so their labels appear nowhere a row
is named: in no estimates column, no coefficient name, no covariance
dimname, no interval rowname, and no printed row. The surface holds
three kinds of row:

- the counterfactual mean of each cell, under the effect label `"mean"`,
  with the cell as its contrast and `"overall"` as its group;

- the simple effects, each treatment's effect within a fixed level of
  the other, with the treatment as the contrast, written `"a: 1 vs 0"`,
  and the level the other is held at as the group, written `"e = 0"`.
  These include the comparisons cell-against-cell reporting cannot
  express at all, such as the first treatment's effect with the second
  set to one;

- the interaction, the difference between two of the first treatment's
  simple effects, with the two compared levels of the second as its
  group, written `"e = 1 vs e = 0"`. It is reported once. Interaction is
  symmetric in the two treatments, so the difference between the first's
  simple effects is the difference between the second's, and reporting
  both would put one quantity in the table under two names.

A two-by-two crossing therefore reports fourteen rows for a binary
outcome: four means, four simple effects on each of two scales, and the
interaction on each of them. A continuous outcome reports nine, since it
has one scale.

The group column names something different here than it names under
`.by`. A group naming a level of the other treatment names the value
that treatment is set to in the two cell means the row contrasts, which
is a setting of the intervention rather than a subgroup of units. Every
row on this surface, the cell means included, standardizes over the
whole sample, so the contrasts are differences and double differences of
cell means over one population.

No contrast row carries a log odds ratio, for the reason no subgroup row
does: an odds ratio is noncollapsible, so neither a simple effect
reported beside one nor a difference of two of them says what it appears
to. The `"mean"` rows are means and carry no scale of their own.

Every row is a parameter of the same stacked system. The cell means are
the block the categorical path already carries, and the contrast block
is written over those same means, so a simple effect and the
cell-against-cell contrast that happens to equal it report the same
estimate and the same standard error. Each interaction row is the
difference of two simple-effect parameters, which makes it the double
difference of the four means by construction rather than by two
arithmetic routes that have to agree.

The declaration is read off the exposure column of the frame the method
resolves, which is `.data` where the caller supplied one and the outcome
model's own frame otherwise. The fit records its levels as plain strings
and keeps no memory of the crossing, so it is the frame and not the fit
that decides which surface is reported, and dropping the declaration
with `factor(x)` returns the cell-against-cell rows.

A crossing whose two treatments carry one name, as in
`joint_exposure(a = x, a = y)`, raises `balancing_ipw_input_error`.
Every row is keyed by the treatment it contrasts and the level the other
is held at, both written from the component names, so one name over two
treatments names two different effects the same way. Declare each
treatment under a name of its own.

Two further configurations raise `balancing_ipw_unsupported_error`. A
declared crossing with `.by` is refused: effect modification of a joint
intervention is a three-way question, the interaction between two
treatments within the levels of a third variable, and this surface
reports neither that nor a projection of it. A declared crossing
weighted for anything but `"ate"` is refused as well: every cell mean
here standardizes to one population, and a tilted estimand standardizes
each of them to a population the simple effects and the interaction are
not defined over.

## Multiple imputation

Missing covariate data is handled by imputing first and analyzing within
each completed dataset: impute, then balance, weight, and fit the
outcome model once per imputation, then pool the per-imputation results
with
[`causalgenerics::pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html).
That function combines them by Rubin's rules with a Barnard-Rubin
degrees-of-freedom adjustment, and documents the rules and the
agreements the results have to satisfy.

The per-imputation analysis is written as an
[`lapply()`](https://rdrr.io/r/base/lapply.html) over the completed
datasets rather than as `with(imp, ...)`.
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
takes its data as the first argument and selects covariates out of it
with tidyselect, so it needs that frame as a named object; inside
[`with()`](https://rdrr.io/r/base/with.html) the completed frame is the
evaluation environment instead and has no name to give. Walking the
imputations directly also keeps the frame in hand for the weight column
and the outcome model, which the same step needs.

Only the estimating-equation methods at exact balance reach this
workflow at all, for the reason they are the only ones
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
accepts anywhere: the pooled standard error is built from per-imputation
standard errors, and a fit carrying no estimating equations produces no
result to pool. The refusal comes at the per-imputation
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
call rather than at the pooling step.

The complete-data degrees of freedom are found without being told. A
balancing result reports none of its own, since the stacked system is
not a fit with residual degrees of freedom, so
[`causalgenerics::pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html)
falls through to the outcome models and takes the smallest residual
degrees of freedom among them. Pooling the same results with
[`mice::pool()`](https://amices.org/mice/reference/pool.html) instead
needs two things this route supplies on its own: the propensity package
loaded, since the `tidy()` method it reaches an `ipw` result through
lives there rather than in balancing, and that count passed explicitly
as its `dfcom` argument, since it reads only what the results themselves
report.

A pooled balancing result carries both readings.
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
hands every outcome model over already wrapped with its block of the
corrected covariance, so a set of balancing results always has a
conditional surface to pool beside the marginal one, and
[`causalgenerics::pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html)
pools both whichever reading it is asked for.
[`causalgenerics::as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
and
[`causalgenerics::as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
therefore move the pooled result between the readings after pooling, as
they move an unpooled one, and the pooled
[`stats::coef()`](https://rdrr.io/r/stats/coef.html),
[`stats::vcov()`](https://rdrr.io/r/stats/vcov.html),
[`stats::confint()`](https://rdrr.io/r/stats/confint.html), and
[`base::as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
take an `effects` argument naming a reading for one call.

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
#> ℹ Treating `.exposure` as binary
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

# `.by` reports the effects again within the levels of a modifier, then
# contrasts each level against the first of them.
df$grp <- factor(ifelse(x1 > 0, "high", "low"), levels = c("low", "high"))
by_mod <- glm(
  y ~ exposure * grp,
  data = df,
  family = quasibinomial(),
  weights = .wts
)

ipw(fit, by_mod, .by = grp)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> 
#> Propensity Score Model:
#>   Call: NULL 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ exposure * grp, family = quasibinomial(), data = df, 
#>     weights = .wts) 
#> 
#> Estimates:
#> Warning: non-unique values when setting 'row.names': ‘log(rr)’, ‘rd’
#> Error in `.rowNamesDF<-`(x, value = value): duplicate 'row.names' are not allowed

# A categorical exposure reports each level against the reference level, and
# the estimates table names the contrast.
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
#> ℹ Treating `.exposure` as categorical
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

# A continuous exposure reports the dose response of a weighted marginal
# structural model. An exposure entering through one design column is that
# response's slope, reported as one row named for the model's link.
df$dose <- 0.7 * x1 + rnorm(n)
df$score <- 2 + 0.5 * df$dose + 0.4 * x1 + rnorm(n)

dose_fit <- balance(df, dose, x1, method = bw_entropy(), estimand = "ate")
#> ℹ Treating `.exposure` as continuous
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

# An exposure entering through several columns reports one row per
# coefficient, named after the coefficient the fit names.
curve_mod <- lm(score ~ poly(dose, 2), data = df, weights = .dose_wts)

ipw(dose_fit, curve_mod)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> 
#> Propensity Score Model:
#>   Call: NULL 
#> 
#> Outcome Model:
#>   Call: lm(formula = score ~ poly(dose, 2), data = df, weights = .dose_wts) 
#> 
#> Estimates:
#> Warning: non-unique value when setting 'row.names': ‘coef’
#> Error in `.rowNamesDF<-`(x, value = value): duplicate 'row.names' are not allowed

# With missing covariate data, analyze within each completed dataset and
# pool the results afterward.
set.seed(2024)
n <- 150
mx1 <- rnorm(n)
mx2 <- rnorm(n)
mz <- rbinom(n, 1, plogis(0.6 * mx1 - 0.4 * mx2))
my <- rbinom(n, 1, plogis(-0.3 + 0.5 * mz + 0.4 * mx1))
incomplete <- data.frame(exposure = mz, x1 = mx1, x2 = mx2, y = my)
incomplete$x1[rbinom(n, 1, plogis(-1.2 + 0.5 * mx2)) == 1] <- NA

imp <- mice::mice(incomplete, m = 2, print = FALSE, seed = 4321)

fits <- lapply(mice::complete(imp, "all"), function(completed) {
  completed_fit <- balance(
    completed,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  completed$.wts <- weights(completed_fit)
  ipw(
    completed_fit,
    glm(
      y ~ exposure,
      data = completed,
      family = quasibinomial(),
      weights = .wts
    )
  )
})
#> ℹ Treating `.exposure` as binary
#> ℹ Treating `.exposure` as binary

pool_ipw(fits)
#> Pooled Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> Imputations: 2 
#> Complete-data df: 148 
#> 
#> Pooled marginal estimates:
#>         estimate  std.err      t     df  ci.lower ci.upper conf.level p.value
#> rd      0.093662 0.086512 1.0827 64.596 -0.079134  0.26646       0.95  0.2830
#> log(rr) 0.178972 0.165070 1.0842 67.103 -0.150501  0.50844       0.95  0.2821
#> log(or) 0.376786 0.350486 1.0750 64.419 -0.323302  1.07687       0.95  0.2864

# The pooled result carries both readings, so it moves to the outcome
# models' coefficients after pooling.
as_conditional(pool_ipw(fits))
#> Pooled Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: conditional (outcome model) 
#> Imputations: 2 
#> Complete-data df: 148 
#> 
#> Pooled conditional estimates (outcome model):
#>              estimate   std.err       t      df ci.lower ci.upper conf.level
#> (Intercept) -0.089361  0.223575 -0.3997 127.868 -0.53175  0.35302       0.95
#> exposure     0.376786  0.350486  1.0750  64.419 -0.32330  1.07687       0.95
#>             p.value
#> (Intercept)  0.6901
#> exposure     0.2864
```
