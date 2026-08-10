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
#' `contrast` column, placed after `effect`, naming each contrast as
#' `"<level> vs <reference>"`. A binary exposure keeps the table it has always
#' returned, with no `contrast` column. A categorical exposure declared as a
#' crossing of two treatments by [causalgenerics::joint_exposure()] replaces
#' those level-against-reference rows with the surface described under Joint
#' exposures below.
#'
#' A continuous exposure has no levels to contrast, so there is no pair of
#' marginal means to difference. What the method reports instead is the
#' dose response of a weighted marginal structural model: the balancing weights
#' break the exposure-covariate association, and every coefficient of the outcome
#' model that reads the exposure is an effect on the model's own link scale. The
#' estimates table holds one row per such coefficient, read straight off the
#' weighted fit, since nothing is standardized here.
#'
#' An exposure entering through one design column, whether as a bare term or as a
#' transformation of one, is the whole of the dose response, so its coefficient
#' is that response's slope everywhere. Such a model keeps the single-row table
#' it has always returned, with the columns every other exposure's table holds
#' and no `contrast` column, and the row is named for the link: `slope` for an
#' identity link, whether the model arrives as a [stats::lm()] or as a gaussian
#' [stats::glm()]; `log(or)` for a logit; and `log(rr)` for a log link.
#'
#' An exposure entering through several columns, as in `y ~ exposure +
#' I(exposure^2)` or `y ~ splines::ns(exposure, 3)`, has no such row, since a
#' curve has a different slope at every dose. The table then holds one row per
#' column, gains a `contrast` column naming each row after the coefficient the
#' fit names, and the scale word steps back to `coef` at an identity link. A
#' logit still reports `log(or)` and a log link `log(rr)`, because a coefficient
#' of those models is a log ratio whatever column it multiplies. Another link
#' raises `balancing_ipw_input_error` either way, since its coefficients are none
#' of those things. A continuous fit targets the average treatment effect and
#' nothing else, so an `estimand` supplied alongside it either agrees or raises
#' `balancing_estimand_error`, as it does for any other fit.
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
#' structural model's coefficients alone, and each reported effect is already one
#' of those coefficients rather than a parameter derived from them. A basis
#' widens that stack by its own columns and changes nothing else.
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
#' The result reads through the accessors causalgenerics registers on the class,
#' which every fitting package's results share: [stats::coef()] for the reported
#' effects under their display labels, [stats::vcov()] for their covariance,
#' [stats::confint()] for their intervals, [stats::nobs()] for the number of
#' observations, and [stats::weights()] for the weights the outcome model was
#' fitted with. The covariance those accessors read is the block of the stacked
#' one belonging to the reported effects, which is on the estimates table rather
#' than derived from it: the effect measures are transformations of the same pair
#' of marginal means, so they covary, and the off-diagonals a caller combining
#' two of them needs are not recoverable from the standard errors alone. The
#' stored outcome model carries its own block of the same covariance, so
#' `vcov()` on it reports the variance of its coefficients with the uncertainty
#' from estimating the weights included, where a bare refit of the same weighted
#' model treats the weights as fixed and understates it.
#'
#' A row's display label is built from the identity columns of the estimates
#' table, in the order the table carries them: the effect measure, then the
#' contrast where the surface names one, then the group where it names one. A
#' contrast names a categorical exposure's pair of levels, a continuous
#' surface's basis coefficient, or a joint exposure's treatment; a group names
#' the subgroup a `.by` request reports the row within, or the level a joint
#' exposure holds the other treatment at. That one label names the row wherever
#' a caller reads it, so the printed table, the names of [stats::coef()], the
#' dimnames of [stats::vcov()], and the rownames of [stats::confint()] agree row
#' for row with the estimates table. The stacked `fit$theta` and `fit$vcov`
#' carry names of their own, which name blocks of the estimating-equation
#' system rather than reported rows and are described under Value.
#'
#' The sandwich covariance is a large-sample one, and how large a sample it
#' takes differs by exposure. The binary risk-difference standard error is
#' calibrated at a few hundred observations; the continuous slope's is
#' anticonservative there. Over 500 draws its ratio of mean standard error to
#' the standard deviation of the estimates was 0.836 at 300 observations and
#' 0.942 at 1200. Read a continuous interval as the asymptotic statement it is,
#' and prefer the bootstrap of the inference vignette at small sample sizes.
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
#' A continuous exposure narrows that by variable membership, since what it
#' reports are coefficients of the model rather than contrasts of predictions
#' from it. Every term reading the exposure must read the exposure alone, however
#' many design columns it expands to, so `y ~ exposure`, `y ~ exposure + x1`,
#' `y ~ exposure + I(exposure^2)`, `y ~ sin(exposure)`, `y ~ poly(exposure, 2)`
#' and `y ~ splines::ns(exposure, 3)` are all supported. A term reading a
#' covariate alongside the exposure raises `balancing_ipw_input_error`, so
#' `y ~ exposure * x1`, `y ~ exposure + exposure:x1` and `y ~ I(exposure * x1)`
#' are refused: each contributes a coefficient that is a change in the dose
#' response per unit of that covariate, so there is no one effect for a row to
#' report and no covariate value a row could name it at. The check reads the
#' model's terms rather than the text of its formula, so the mixing is caught
#' however it is written.
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
#' # Effect modification
#'
#' `.by` names a modifier, and a result carrying one reports the effects it
#' reports without a request, then those same effects within each of the
#' modifier's levels, then each non-reference level against the reference one.
#' The estimates table gains a `group` column naming the subgroup each row was
#' estimated in, placed after `contrast` where a categorical exposure names one
#' and after `effect` where it does not, since a subgroup qualifies the whole
#' comparison rather than one side of it. The whole-sample rows come first,
#' under the group `"overall"`; a subgroup's rows are named `"var = value"`, as
#' `"sex = female"`; and a contrast of subgroups joins the two, as
#' `"sex = female vs sex = male"`. The reference subgroup is the modifier's first
#' level, which for a factor is the first of its declared levels rather than the
#' first in sorted order. A character modifier declares no levels, so it is read
#' as a factor on the way in and its reference subgroup is its alphabetically
#' first value, whatever order the values appear in. Declare the column a factor
#' to measure the contrasts against some other subgroup.
#'
#' A subgroup reports the collapsible measures alone. For a binary outcome that
#' is `rd` and `log(rr)`; a continuous outcome reports `diff`, the only measure
#' it has. The log odds ratio stays among the whole-sample rows. An odds ratio
#' is noncollapsible, so the odds ratio over a sample is not an average of the
#' odds ratios within its subgroups and the difference of two of them is not
#' the difference in effect it reads as; nothing in the whole-sample rows
#' averages anything over subgroups, which is why they keep it. A categorical
#' exposure crosses its
#' contrasts with the subgroups, so each subgroup block reports each contrast in
#' the order the whole-sample block reports it.
#'
#' A subgroup's marginal means are the g-computation means over that subgroup
#' alone, standardized the way the whole-sample means are: over every unit of
#' the subgroup for a pooled estimand, and over its focal units for `"att"` or
#' `"atc"`. Sampling weights weight that average as well.
#'
#' Every row a request adds is a parameter of the same stacked system the
#' whole-sample rows are parameters of. The blocks are appended after everything
#' the ungrouped stack carries and nothing earlier reads a parameter of theirs,
#' so the leading blocks are the ones an ungrouped fit produces and the
#' whole-sample rows of a grouped result are the rows it reported. What that
#' buys is the covariance: the subgroups share the weight parameters and the
#' outcome model's coefficients, so their effects covary, and each contrast of
#' subgroups is a parameter of the joint system rather than a difference of two
#' separate fits, whose variance would have to be read as a sum.
#'
#' A level of the modifier that no unit carries names an empty subgroup and is
#' dropped rather than refused, which is what lets a modifier be subset without
#' being recoded first. Four configurations are refused with
#' `balancing_ipw_input_error`: a selection naming any number of columns other
#' than one, a modifier that is neither a factor nor a character column, a
#' modifier carrying missing values, and a modifier one of whose subgroups holds
#' fewer than all of the exposure levels. That last one has no contrast to
#' report in the subgroup that is short a level, and refitting either model does
#' not supply a comparison the data do not hold. `.by` with a continuous
#' exposure raises `balancing_ipw_unsupported_error`: what such a fit reports is
#' the marginal structural model's own exposure coefficients rather than
#' contrasts of standardized means, so there is no effect within a subgroup for
#' the argument to name, and `.by` with a declared joint exposure raises it for
#' the reason given under Joint exposures below. A modifier that reaches the
#' argument through `.data` carries that frame's row-order requirement with it,
#' described under `.data` above: the subgroups are built from the rows it holds
#' while the weights stay in the fit's order.
#'
#' An outcome model with no term reading both the exposure and the modifier
#' raises `balancing_ipw_by_interaction_warning` and the result is still built.
#' The subgroup effects are g-computation on the model as it was specified, so
#' the effect differs across subgroups only where a term reads both columns, and
#' that may be the modeling choice a caller meant to make. The warning reports
#' which terms were read rather than announcing that the effect is the same in
#' every subgroup, because a model may carry the modification through a column
#' derived from the modifier instead: `y ~ exposure * sex_female` reported by
#' `.by = sex` names no term reading `sex`, and its subgroup effects differ all
#' the same.
#'
#' # Joint exposures
#'
#' [causalgenerics::joint_exposure()] crosses two discrete treatments into one
#' categorical exposure and records the crossing on the vector. The result is a
#' factor, so [balance()] weights it as it weights any factor over those cells,
#' and the crossing changes which effects are reported rather than which
#' population is balanced.
#'
#' Reported as a plain categorical exposure, such a fit gives each cell against
#' the reference cell, under contrasts like `"a = 1, e = 0 vs a = 0, e = 0"`.
#' Those rows are arithmetically right and they answer a question nobody asked.
#' A declared crossing is reported in the two treatments instead, and the
#' cell-against-cell rows are replaced rather than supplemented, so their labels
#' appear nowhere a row is named: in no estimates column, no coefficient name,
#' no covariance dimname, no interval rowname, and no printed row. The surface
#' holds three kinds of row:
#'
#' * the counterfactual mean of each cell, under the effect label `"mean"`, with
#'   the cell as its contrast and `"overall"` as its group;
#' * the simple effects, each treatment's effect within a fixed level of the
#'   other, with the treatment as the contrast, written `"a: 1 vs 0"`, and the
#'   level the other is held at as the group, written `"e = 0"`. These include
#'   the comparisons cell-against-cell reporting cannot express at all, such as
#'   the first treatment's effect with the second set to one;
#' * the interaction, the difference between two of the first treatment's simple
#'   effects, with the two compared levels of the second as its group, written
#'   `"e = 1 vs e = 0"`. It is reported once. Interaction is symmetric in the two
#'   treatments, so the difference between the first's simple effects is the
#'   difference between the second's, and reporting both would put one quantity
#'   in the table under two names.
#'
#' A two-by-two crossing therefore reports fourteen rows for a binary outcome:
#' four means, four simple effects on each of two scales, and the interaction on
#' each of them. A continuous outcome reports nine, since it has one scale.
#'
#' The group column names something different here than it names under `.by`. A
#' group naming a level of the other treatment names the value that treatment is
#' set to in the two cell means the row contrasts, which is a setting of the
#' intervention rather than a subgroup of units. Every row on this surface, the
#' cell means included, standardizes over the whole sample, so the contrasts are
#' differences and double differences of cell means over one population.
#'
#' No contrast row carries a log odds ratio, for the reason no subgroup row does:
#' an odds ratio is noncollapsible, so neither a simple effect reported beside
#' one nor a difference of two of them says what it appears to. The `"mean"`
#' rows are means and carry no scale of their own.
#'
#' Every row is a parameter of the same stacked system. The cell means are the
#' block the categorical path already carries, and the contrast block is written
#' over those same means, so a simple effect and the cell-against-cell contrast
#' that happens to equal it report the same estimate and the same standard
#' error. Each interaction row is the difference of two simple-effect
#' parameters, which makes it the double difference of the four means by
#' construction rather than by two arithmetic routes that have to agree.
#'
#' The declaration is read off the exposure column of the frame the method
#' resolves, which is `.data` where the caller supplied one and the outcome
#' model's own frame otherwise. The fit records its levels as plain strings and
#' keeps no memory of the crossing, so it is the frame and not the fit that
#' decides which surface is reported, and dropping the declaration with
#' `factor(x)` returns the cell-against-cell rows.
#'
#' A crossing whose two treatments carry one name, as in
#' `joint_exposure(a = x, a = y)`, raises `balancing_ipw_input_error`. Every row
#' is keyed by the treatment it contrasts and the level the other is held at,
#' both written from the component names, so one name over two treatments names
#' two different effects the same way. Declare each treatment under a name of its
#' own.
#'
#' Two further configurations raise `balancing_ipw_unsupported_error`. A declared
#' crossing with `.by` is refused: effect modification of a joint intervention
#' is a three-way question, the interaction between two treatments within the
#' levels of a third variable, and this surface reports neither that nor a
#' projection of it. A declared crossing weighted for anything but `"ate"` is
#' refused as well: every cell mean here standardizes to one population, and a
#' tilted estimand standardizes each of them to a population the simple effects
#' and the interaction are not defined over.
#'
#' # Multiple imputation
#'
#' Missing covariate data is handled by imputing first and analyzing within each
#' completed dataset: impute, then balance, weight, and fit the outcome model
#' once per imputation, then pool the per-imputation results with
#' [causalgenerics::pool_ipw()]. That function combines them by Rubin's rules
#' with a Barnard-Rubin degrees-of-freedom adjustment, and documents the rules
#' and the agreements the results have to satisfy.
#'
#' The per-imputation analysis is written as an `lapply()` over the completed
#' datasets rather than as `with(imp, ...)`. `balance()` takes its data as the
#' first argument and selects covariates out of it with tidyselect, so it needs
#' that frame as a named object; inside `with()` the completed frame is the
#' evaluation environment instead and has no name to give. Walking the
#' imputations directly also keeps the frame in hand for the weight column and
#' the outcome model, which the same step needs.
#'
#' Only the estimating-equation methods at exact balance reach this workflow at
#' all, for the reason they are the only ones `ipw()` accepts anywhere: the
#' pooled standard error is built from per-imputation standard errors, and a fit
#' carrying no estimating equations produces no result to pool. The refusal
#' comes at the per-imputation `ipw()` call rather than at the pooling step.
#'
#' The complete-data degrees of freedom are found without being told. A
#' balancing result reports none of its own, since the stacked system is not a
#' fit with residual degrees of freedom, so [causalgenerics::pool_ipw()] falls
#' through to the outcome models and takes the smallest residual degrees of
#' freedom among them. Pooling the same results with [mice::pool()] instead
#' needs two things this route supplies on its own: the propensity package
#' loaded, since the `tidy()` method it reaches an `ipw` result through lives
#' there rather than in balancing, and that count passed explicitly as its
#' `dfcom` argument, since it reads only what the results themselves report.
#'
#' A pooled balancing result carries both readings. `ipw()` hands every outcome
#' model over already wrapped with its block of the corrected covariance, so a
#' set of balancing results always has a conditional surface to pool beside the
#' marginal one, and [causalgenerics::pool_ipw()] pools both whichever reading
#' it is asked for. [causalgenerics::as_marginal()] and
#' [causalgenerics::as_conditional()] therefore move the pooled result between
#' the readings after pooling, as they move an unpooled one, and the pooled
#' [stats::coef()], [stats::vcov()], [stats::confint()], and
#' [base::as.data.frame()] take an `effects` argument naming a reading for one
#' call.
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
#'   For a continuous exposure it is a marginal structural model whose every
#'   exposure-reading term reads the exposure alone.
#' @param .data The data frame holding the exposure and outcome. If `NULL`, the
#'   values are taken from the outcome model frame. It carries the exposure
#'   column the fixed-exposure predictions are built from, so it has nothing to
#'   supply for a continuous exposure, which makes no such predictions, and is
#'   ignored there.
#'
#'   It must hold the rows both models were fitted on, in the same order.
#'   Half of the stacked system reads it and half reads the fit: the
#'   counterfactual predictions, the subgroup indicators, and a focal estimand's
#'   standardization come from `.data`, while the weight equations and the
#'   outcome-model score stay in the fit's order. A frame holding the right rows
#'   in another order therefore leaves every effect estimate unchanged, since a
#'   weighted mean does not care in which order it is summed, and makes every
#'   standard error wrong. Each column `.data` and the outcome model frame name
#'   in common is compared value by value, and a disagreement raises
#'   `balancing_ipw_input_error` naming the column and the first row it
#'   disagrees at. Columns the outcome model never saw, the modifier `.by` names
#'   among them, are free.
#' @param estimand The causal estimand. If `NULL`, the fit's estimand is used.
#'   As in [balance()], `"atc"` is accepted as a synonym for `"atu"`. Supplying
#'   an estimand that disagrees with the fit raises `balancing_estimand_error`.
#' @param conf_level The confidence level for the intervals. Default `0.95`.
#' @param effects The presentation mode the result records, either `"marginal"`
#'   (the default) or `"conditional"`. The marginal reading reports the
#'   population-averaged causal contrasts described above; the conditional
#'   reading reports the outcome model's coefficient surface. Both surfaces are
#'   computed whichever mode is named, since the stacked system is solved either
#'   way, so the argument settles which one the result presents and nothing
#'   else. [causalgenerics::as_marginal()] and [causalgenerics::as_conditional()]
#'   move a result between the two readings afterwards, a pooled result as much
#'   as an unpooled one, and the accessors take an `effects` argument of their
#'   own for a single call.
#'
#'   The conditional reading reports the coefficients of the stored
#'   `outcome_mod` against the outcome block of the stacked sandwich, which is
#'   the block the wrapper around that model carries, so its standard errors
#'   account for having estimated the weights rather than treating them as
#'   fixed. The stored `wt_mod` is the balancing fit itself rather than a
#'   wrapper of one, since the weight block of the same stack is already its
#'   own: the stored copy carries that block, the covariance of the fit's
#'   weight parameters, and [stats::vcov()] on it reads the block back. The fit
#'   that went into the call is untouched, and reports no covariance of its own.
#' @param .by A modifier to report the effects within the levels of, given
#'   unquoted and selected with tidyselect out of `.data` where one was supplied
#'   and out of the outcome model frame otherwise. The default, `NULL`, is the
#'   absence of a request: the result then reports the whole-sample effects
#'   alone and names no subgroups at all. A selection reaching any number of
#'   columns other than one is refused, since the effects are reported within
#'   the levels of a single variable. The effect modification section below
#'   describes the rows a request adds and the configurations it refuses, and
#'   the joint exposure section describes why a declared crossing takes none.
#' @param ... Ignored, for compatibility with the generic.
#'
#' @return An object of class `ipw`, an implementation of
#'   [causalgenerics::ipw()]. Alongside `estimand`, `wt_mod`, `outcome_mod`, the
#'   `estimates` table, and the `effects` field recording the presentation mode
#'   described above, the result carries two fields describing the variance:
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
#'     neither block, since its effects are outcome-model coefficients; each of
#'     those coefficients is named for the label its estimates row carries, as
#'     `slope` or `coef I(exposure^2)`, in place of the `beta_` name the others
#'     keep. A `.by` request appends, after all of those, every subgroup's mean
#'     block in subgroup order, then every subgroup's contrast block in that
#'     same order, then one contrast block per non-reference subgroup against
#'     the reference one. Each of their names is the name of the block it
#'     repeats, suffixed with its group, as
#'     `mu0_sex = female` and `rd_sex = female vs sex = male`. A declared joint
#'     exposure keeps the mean block and replaces the contrast block, naming
#'     each of its own contrasts for the measure and the row, as
#'     `rd_a: 1 vs 0 e = 0`. The standard errors in `estimates` are
#'     `sqrt(diag(fit$vcov))` read at those effect names.
#'
#'   The `estimates` table carries the covariance of the reported effects as its
#'   `ipw_vcov` attribute, which is what [stats::vcov()] returns in the marginal
#'   reading. Both its dimnames are the display labels described above, as
#'   `"rd b vs a sex = female"`. The stored `outcome_mod` is wrapped by
#'   [causalgenerics::new_ipw_model()], which carries the outcome-model block of
#'   `fit$vcov` under the model's own coefficient names, so `vcov()` on it
#'   reports the joint-estimation variance. The stored `wt_mod` carries the
#'   leading block of the same covariance under the `theta_w` names above, which
#'   name the fit's own parameters on either route, so `vcov()` on it reports
#'   the covariance of the weight parameters. [stats::df.residual()] returns
#'   `NA_integer_`, since the stacked system is not a fit with residual degrees
#'   of freedom of its own.
#'
#'   [stats::nobs()] delegates to the stored outcome model, which counts the
#'   rows it was fitted on that carry a nonzero weight. A unit given no sampling
#'   weight is pinned at zero rather than dropped, so the weight vector stays the
#'   length of the data the fit saw while the outcome model counts one row fewer
#'   for each pinned unit, and `nobs()` on the result is then smaller than
#'   `length(weights(result))`.
#'
#' @examples
#' n <- 200
#' x1 <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.5 * x1))
#' y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.3 * x1))
#' df <- data.frame(exposure = z, x1 = x1, y = y)
#'
#' fit <- balance(df, exposure, x1, method = bw_entropy(), estimand = "ate")
#' df$.wts <- weights(fit)
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
#' # `.by` reports the effects again within the levels of a modifier, then
#' # contrasts each level against the first of them.
#' df$grp <- factor(ifelse(x1 > 0, "high", "low"), levels = c("low", "high"))
#' by_mod <- glm(
#'   y ~ exposure * grp,
#'   data = df,
#'   family = quasibinomial(),
#'   weights = .wts
#' )
#'
#' ipw(fit, by_mod, .by = grp)
#'
#' # A categorical exposure reports each level against the reference level, and
#' # the estimates table names the contrast.
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
#' df$.arm_wts <- weights(arm_fit)
#' arm_mod <- glm(
#'   relapse ~ arm,
#'   data = df,
#'   family = quasibinomial(),
#'   weights = .arm_wts
#' )
#'
#' ipw(arm_fit, arm_mod)
#'
#' # A continuous exposure reports the dose response of a weighted marginal
#' # structural model. An exposure entering through one design column is that
#' # response's slope, reported as one row named for the model's link.
#' df$dose <- 0.7 * x1 + rnorm(n)
#' df$score <- 2 + 0.5 * df$dose + 0.4 * x1 + rnorm(n)
#'
#' dose_fit <- balance(df, dose, x1, method = bw_entropy(), estimand = "ate")
#' df$.dose_wts <- weights(dose_fit)
#' dose_mod <- lm(score ~ dose, data = df, weights = .dose_wts)
#'
#' ipw(dose_fit, dose_mod)
#'
#' # An exposure entering through several columns reports one row per
#' # coefficient, named after the coefficient the fit names.
#' curve_mod <- lm(score ~ poly(dose, 2), data = df, weights = .dose_wts)
#'
#' ipw(dose_fit, curve_mod)
#'
#' @examplesIf requireNamespace("mice", quietly = TRUE)
#' # With missing covariate data, analyze within each completed dataset and
#' # pool the results afterward.
#' set.seed(2024)
#' n <- 150
#' mx1 <- rnorm(n)
#' mx2 <- rnorm(n)
#' mz <- rbinom(n, 1, plogis(0.6 * mx1 - 0.4 * mx2))
#' my <- rbinom(n, 1, plogis(-0.3 + 0.5 * mz + 0.4 * mx1))
#' incomplete <- data.frame(exposure = mz, x1 = mx1, x2 = mx2, y = my)
#' incomplete$x1[rbinom(n, 1, plogis(-1.2 + 0.5 * mx2)) == 1] <- NA
#'
#' imp <- mice::mice(incomplete, m = 2, print = FALSE, seed = 4321)
#'
#' fits <- lapply(mice::complete(imp, "all"), function(completed) {
#'   completed_fit <- balance(
#'     completed,
#'     exposure,
#'     c(x1, x2),
#'     method = bw_entropy(),
#'     estimand = "ate"
#'   )
#'   completed$.wts <- weights(completed_fit)
#'   ipw(
#'     completed_fit,
#'     glm(
#'       y ~ exposure,
#'       data = completed,
#'       family = quasibinomial(),
#'       weights = .wts
#'     )
#'   )
#' })
#'
#' pool_ipw(fits)
#'
#' # The pooled result carries both readings, so it moves to the outcome
#' # models' coefficients after pooling.
#' as_conditional(pool_ipw(fits))
#'
#' @name ipw.balancing
#' @importFrom causalgenerics ipw
#' @importFrom stats getCall
#' @importFrom stats vcov
NULL

# Expose the originating call through the standard model accessor so tools that
# summarize a fit, including the ipw() print method, can label it. Without a
# method the S7 object is not subsettable and the accessor would fail.
getCall_generic <- new_external_generic("stats", "getCall", "x")

method(getCall_generic, balancing) <- function(x, ...) {
  x@call
}

# Report the covariance of the weight parameters through the standard model
# accessor. Nothing in this package computes one for a fit on its own. The block
# arrives only as a by-product of the stacked assembly an ipw() result performs,
# which fills it in on the copy of the fit it stores, so a fit that has not been
# through that assembly has no block to report and refuses rather than returning
# an empty matrix.
#
# One class covers both ways of arriving at that refusal, since the caller's
# position is the same either way, and the guidance is what separates them. A
# fit whose weights solve estimating equations has the ipw() route open to it,
# so the guidance names it. A fit from a method that solves none does not:
# ipw() turns such a fit away as well, so naming the route would walk its holder
# into a second refusal, and the guidance names the bootstrap workflow instead.
vcov_generic <- new_external_generic("stats", "vcov", "object")

method(vcov_generic, balancing) <- function(object, ...) {
  if (is.null(object@vcov)) {
    guidance <- if (is.null(object@estimating_equations)) {
      c(
        i = "This fit's weights solve no estimating equations, so there is no such result to read one from.",
        i = "Use the bootstrap workflow in the inference vignette for variance instead."
      )
    } else {
      c(
        i = "Build an {.fun ipw} result from this fit and read {.fun vcov} off the fit it stores as {.arg wt_mod}."
      )
    }
    abort(
      c(
        "This fit carries no covariance for its weight parameters.",
        x = "A covariance for them comes from the stacked system an {.fun ipw} result assembles, and this fit has not been through one.",
        guidance
      ),
      error_class = "balancing_vcov_error"
    )
  }
  object@vcov
}

causalgenerics_ipw <- new_external_generic("causalgenerics", "ipw", "wt_mod")

method(causalgenerics_ipw, balancing) <- function(
  wt_mod,
  outcome_mod,
  .data = NULL,
  estimand = NULL,
  conf_level = 0.95,
  effects = c("marginal", "conditional"),
  .by = NULL,
  ...
) {
  # The reading reaches only the constructor: it names no part of the stacked
  # system and nothing below branches on it. So it is settled ahead of the
  # checks on the two models, and a call that is wrong in the reading and in a
  # model reports the reading rather than making the caller fix the model and
  # meet this refusal on the next attempt.
  effects <- rlang::arg_match(effects)

  # The modifier is selected out of a frame this function has not resolved yet,
  # so the request is defused here and evaluated below, once the frame the
  # counterfactual designs are built from is in hand.
  .by <- rlang::enquo(.by)

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

  # A continuous exposure has no effect within a subgroup for a request to name,
  # and that is settled by the fit alone. Refusing here rather than after the
  # outcome model is inspected is what keeps a caller from being sent to fix a
  # model whose result the request could not have been answered from anyway.
  if (continuous_exposure) {
    check_ipw_by_exposure(.by)
  }

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
  # effects it reports are among those coefficients rather than contrasts of
  # marginal means. A basis is reported from the fitted objects alone for the
  # same reason, which is what lets one whose model frame records the basis and
  # carries no exposure column be reported without `.data`.
  if (continuous_exposure) {
    variance_system <- ipw_deli_msm_sandwich(
      container = container,
      outcome_mod = outcome_mod,
      exposure_name = exposure_name,
      sampling_weights = wt_mod@sampling_weights,
      call = rlang::current_env()
    )
    estimates <- ipw_estimates_from_identity(
      msm_coefficient_identity(outcome_mod, exposure_name),
      theta = variance_system$theta,
      vcov = variance_system$vcov,
      conf_level = conf_level
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

    # A supplied frame has to be the fit's rows in the fit's order, which the
    # row count alone does not say. The exposure check runs first because a
    # frame describing a different set of exposure levels is a more specific
    # complaint than a misaligned one, and reporting it that way keeps the
    # remedy pointed at the column the caller changed.
    validate_ipw_frame_alignment(.data, outcome_mod)

    # A crossing is declared on the exposure column rather than recorded by the
    # fit, which stores its levels as plain strings and keeps no memory of what
    # they were built from. So the frame resolved above is the only place the
    # declaration can be read, and it is read once, whichever arm filled the
    # frame in.
    #
    # All three refusals come before anything is reported, since each of them is
    # about whether the surface can be written at all rather than about what it
    # would say. The shared-name check comes first among them: it asks whether
    # the declaration can be named at all, where the other two ask whether this
    # surface can be reported for the estimand the fit carries and beside a
    # modifier.
    joint <- ipw_joint_plan(
      frame[[exposure_name]],
      levels,
      is_gaussian_outcome(outcome_mod)
    )
    check_ipw_joint_components(joint, frame[[exposure_name]])
    check_ipw_joint_estimand(joint, estimand)
    check_ipw_joint_by(joint, .by)

    # The counterfactual designs fix the exposure to one cell at a time, and a
    # declared crossing asked to give up its other cells gives up its
    # declaration and says so. The plan is read off the declaration above and
    # everything below works from the cells alone, so the column goes on as a
    # plain factor over them.
    if (!is.null(joint)) {
      frame[[exposure_name]] <- ipw_joint_bare(frame[[exposure_name]])
    }

    # The modifier is read out of the same frame the counterfactual designs are
    # built from, so the strata and the predictions they standardize describe
    # one set of rows. The exposure goes with it, since a stratum holding only
    # some of the exposure levels identifies no contrast there.
    by <- ipw_resolve_by(
      .by,
      frame = frame,
      exposure = frame[[exposure_name]],
      exposure_levels = levels,
      exposure_name = exposure_name,
      outcome_mod = outcome_mod
    )

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
      by = by,
      joint = joint,
      sampling_weights = wt_mod@sampling_weights,
      focal_level = wt_mod@focal_level,
      call = rlang::current_env()
    )

    estimates <- ipw_estimates(
      theta = variance_system$theta,
      vcov = variance_system$vcov,
      conf_level = conf_level,
      continuous = is_gaussian_outcome(outcome_mod),
      levels = if (categorical) levels else NULL,
      by = by,
      joint = joint
    )
  }

  # The weight parameters lead the stack on either route, so their own
  # covariance is its leading block and one slice serves both exposure paths.
  # The fit is where that block belongs, since it describes the parameters the
  # fit solved for rather than anything the result reports, and filling it in
  # here reaches only this function's copy: an S7 object is a value, so the fit
  # the caller passed still carries no covariance and still refuses to report
  # one. The block keeps the stacked `theta_w` names, which name the fit's own
  # parameters on either route.
  weight_parameters <- length(container@parameters)
  wt_mod@vcov <- variance_system$vcov[
    seq_len(weight_parameters),
    seq_len(weight_parameters),
    drop = FALSE
  ]

  # The result is built by the shared constructor, so it carries the fields the
  # common print and as.data.frame methods read. The stacked parameter vector
  # and its covariance are the whole of the fitted variance system here, so they
  # stand in for the solver object a propensity score fit reports: nothing was
  # solved, since every parameter entered at the value its own fit had already
  # found.
  causalgenerics::new_ipw(
    estimand = estimand,
    wt_mod = wt_mod,
    outcome_mod = wrap_outcome_model(
      outcome_mod,
      variance_system$vcov,
      weight_parameters
    ),
    estimates = estimates,
    se_method = "mestimation",
    fit = variance_system,
    effects = effects
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

# Refuse a supplied frame whose rows are not the rows the models were fitted on,
# in the same order.
#
# Half of the stacked system reads the fit's row order and the other half reads
# `.data`. The weight-parameter equations, the outcome-model score, and every
# cross term of the meat come from the fit and the fitted model; the
# counterfactual designs, the stratum indicators, and a focal estimand's tilt
# are built from `.data`. A frame carrying the right rows in the wrong order
# pairs each unit's prediction with another unit's weight. The marginal means
# survive that, since a weighted mean does not care in which order it is summed,
# so the point estimates come back unchanged while every standard error moves,
# which is the shape of mistake nothing downstream can notice.
#
# The comparison is against the outcome model's own frame, which is the one
# object known to be in the fit's order: the weight preflight has already
# compared its prior weights against the fit's weights unit by unit, so a model
# fitted on reordered rows never reaches here. Every column the two frames name
# in common is compared element by element, which catches a reordering and an
# altered value alike.
#
# Values are compared rather than objects. A factor supplied with a level the
# model frame's copy does not declare is the same column read under a wider
# declaration, which is what the dropped-level path relies on, so factors are
# read as their labels and two missing values count as agreeing.
#
# What the comparison reaches is the shared column names and no more, which is
# worth stating because one workflow shares fewer of them than the rest. A model
# whose formula transforms the exposure stores `factor(arm)` rather than `arm`,
# so a frame supplied for that model is checked on whatever else it names in
# common, which is usually the response alone, and a permutation that leaves the
# response fixed would pass. Closing that means rebuilding the model's own
# variables from the supplied frame and comparing those, which reaches every
# model at the cost of a second design build.
validate_ipw_frame_alignment <- function(
  .data,
  outcome_mod,
  call = rlang::caller_env()
) {
  if (is.null(.data)) {
    return(invisible(NULL))
  }
  model_frame <- stats::model.frame(outcome_mod)

  for (column in intersect(names(model_frame), names(.data))) {
    supplied <- alignment_values(.data[[column]])
    fitted <- alignment_values(model_frame[[column]])

    # Two missing values agree and one missing value does not, which takes
    # writing out: a comparison against a missing value is itself missing, and
    # `which()` drops a missing element rather than reporting it, so a row the
    # two frames disagree about by one of them being absent would otherwise
    # read as a row they agree about.
    equal <- supplied == fitted
    equal[is.na(equal)] <- FALSE
    differs <- which(!equal & !(is.na(supplied) & is.na(fitted)))
    if (length(differs) == 0L) {
      next
    }

    position <- differs[[1L]]
    abort(
      c(
        "{.arg .data} must hold the rows the models were fitted on, in the same order.",
        x = "Its {.val {column}} column disagrees with the outcome model frame's at row {position}, and at {length(differs)} row{?s} in all.",
        i = "The counterfactual predictions, the subgroup indicators, and a focal estimand's standardization are built from {.arg .data}, while the weight equations and the outcome-model score stay in the fit's order.",
        i = "A frame holding the right rows in the wrong order therefore leaves the effect estimates unchanged and every standard error wrong.",
        i = "Supply the frame the outcome model was fitted on, or omit {.arg .data} when the model frame already carries the exposure."
      ),
      error_class = "balancing_ipw_input_error",
      call = call,
      .envir = environment()
    )
  }
  invisible(NULL)
}

# One column on the scale the alignment comparison reads it. A factor becomes
# its labels, so that two columns holding the same values under different level
# declarations agree; everything else is compared as it stands, which lets an
# integer column and the double the model frame stored agree as numbers rather
# than disagreeing as types.
alignment_values <- function(x) {
  if (is.factor(x)) as.character(x) else x
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
  unique(unlist(model_term_variable_sets(outcome_mod)))
}

# The same reading term by term, which is what a check about two variables
# entering one term needs: the union above says only that both are somewhere in
# the model, and a model carrying each of them alone would answer that question
# the same way as one carrying their interaction.
model_term_variable_sets <- function(outcome_mod) {
  labels <- attr(stats::terms(outcome_mod), "term.labels")
  lapply(labels, function(label) all.vars(str2lang(label)))
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
    # one. A continuous exposure reports the coefficients of the terms reading
    # the exposure, so a transformation is fine there as well while a term
    # reading a covariate alongside it is not, and the check below refuses that
    # one; offering the covariate latitude here would send that caller to a
    # second refusal.
    latitude <- if (continuous_exposure) {
      "The model may adjust for covariates alongside the exposure, and may carry the exposure inside a transformation such as {.fun poly}, so long as no term reads a covariate alongside it."
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
    validate_ipw_exposure_terms(outcome_mod, exposure_name, call = call)
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

# A continuous exposure reports coefficients of the outcome model rather than
# contrasts of predictions from it, so what the model may say about the exposure
# is set by which variables each of its terms reads. A term reading the exposure
# alone contributes coefficients of the dose response itself, however it is
# written and however many design columns it expands to, and each of them is a
# row of the reported surface. A term reading a covariate as well contributes a
# coefficient that is a change in the dose response per unit of that covariate,
# so the effect it describes is defined only at a value of the covariate no row
# could name, and such a model is refused.
#
# The check reads the model's terms rather than its formula text, so the mixing
# is caught however it is written: `y ~ a * x1`, `y ~ a + a:x1`, `y ~ I(a * x1)`
# and `y ~ poly(a, 2):x1` are refused while `y ~ a`, `y ~ a + v`,
# `y ~ a + I(a^2)`, `y ~ sin(a)` and `y ~ splines::ns(a, 3)` are accepted.
# Covariates that do not enter a term with the exposure stay free, since they
# leave the exposure's own columns and their coefficients in place. An offset is
# not a term and so is not this check's to see; the check below reads it.
validate_ipw_exposure_terms <- function(
  outcome_mod,
  exposure_name,
  call = rlang::caller_env()
) {
  labels <- attr(stats::terms(outcome_mod), "term.labels")
  term_variables <- model_term_variable_sets(outcome_mod)
  reads_exposure <- vapply(
    term_variables,
    function(variables) exposure_name %in% variables,
    logical(1)
  )
  reads_more <- vapply(
    term_variables,
    function(variables) length(setdiff(variables, exposure_name)) > 0,
    logical(1)
  )
  mixed <- labels[reads_exposure & reads_more]
  if (length(mixed) == 0) {
    return(invisible(NULL))
  }

  abort(
    c(
      "{.arg outcome_mod} must read the exposure in terms that read nothing else.",
      x = "{cli::qty(mixed)}The term{?s} {.val {mixed}} {?reads/read} the exposure {.val {exposure_name}} alongside another variable.",
      i = "Such a term contributes a coefficient that is a change in the dose response per unit of what it reads, so there is no one effect for a row to report and no value a row could name it at.",
      i = "A term reading the exposure alone is admitted however it is written, so {.code {exposure_name} + I({exposure_name}^2)} and {.code poly({exposure_name}, 2)} each report one row per coefficient.",
      i = "The model may adjust for covariates that do not enter a term with the exposure."
    ),
    error_class = "balancing_ipw_input_error",
    call = call,
    .envir = environment()
  )
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
  # back-quotes a name R cannot parse as a symbol, so an offset supplied through
  # the argument as a bare column name is reported as the model records it.
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
# what each row is called. The two are the same strings for a continuous exposure
# entering through one column, whose single effect is that column's coefficient
# under the name it already carries in the stack, and differ wherever a measure
# repeats: across the contrasts of a categorical exposure, across the subgroups
# of a `.by` request, and across the coefficients of a basis dose response.
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
# measures appear once per contrast. The table therefore gains a `contrast`
# column naming the two levels, placed immediately after `effect`, since the
# column qualifies the measure it follows. That is the column the causalgenerics
# contract names, and every surface built from the table reads it from there. A
# binary exposure has a single contrast and keeps the eight-column table, since a
# column repeating one label on every row identifies nothing.
#
# A `.by` request repeats the measures again across subgroups, so the table
# gains a `group` column on the same terms, after the contrast column where
# there is one and after `effect` where there is not, since a subgroup qualifies
# the whole comparison rather than one side of it. The subgroup blocks follow
# the whole-sample block and run subgroup-major, each of them repeating the
# whole-sample block's own contrast-major order over the measures a subgroup
# reports.
# A declared crossing replaces the vs-reference block rather than adding to it,
# so its rows are described on their own terms in R/ipw-joint.R and the two
# grammars meet here, where each of them is read the same way: keys into the
# stack, and the identity columns the rows are named by.
ipw_estimates <- function(
  theta,
  vcov,
  conf_level,
  continuous,
  levels = NULL,
  by = NULL,
  joint = NULL
) {
  identity <- if (is.null(joint)) {
    ipw_contrast_identity(continuous, levels, by)
  } else {
    ipw_joint_identity(joint)
  }
  ipw_estimates_from_identity(
    identity,
    theta = theta,
    vcov = vcov,
    conf_level = conf_level
  )
}

# The estimates table of any surface, given the description of its rows: which
# entries of the stack they are read from, and which identity columns name them.
# Both exposure types meet here, since a continuous exposure describes its own
# rows rather than deriving them from levels and a modifier, and what follows is
# the same for either: read the rows, add the identity columns that name
# something, attach the covariance.
ipw_estimates_from_identity <- function(identity, theta, vcov, conf_level) {
  estimates <- ipw_estimate_rows(
    theta = theta,
    vcov = vcov,
    conf_level = conf_level,
    keys = identity$keys,
    effects = identity$effect
  )
  identity_columns <- c(
    if (!is.null(identity$contrast)) list(contrast = identity$contrast),
    if (!is.null(identity$group)) list(group = identity$group)
  )
  if (length(identity_columns) > 0) {
    estimates <- cbind(
      estimates["effect"],
      as.data.frame(identity_columns, stringsAsFactors = FALSE),
      estimates[setdiff(names(estimates), "effect")]
    )
  }
  # The covariance is attached last because `cbind()` rebuilds the table, which
  # would drop an attribute attached before it, and because the labels it is
  # named by are read off the finished table.
  attach_effect_covariance(estimates, vcov = vcov, keys = identity$keys)
}

# The stacked keys and the identity columns of the vs-reference surface, which
# is what every exposure reports when nothing else is declared or requested.
ipw_contrast_identity <- function(continuous, levels, by) {
  keys <- ipw_contrast_names(continuous, levels)
  measures <- ipw_contrast_names(continuous)
  contrast_labels <- if (is.null(levels)) {
    NULL
  } else {
    paste(levels[-1], "vs", levels[[1]])
  }

  effect <- rep(measures, times = length(keys) / length(measures))
  contrast <- if (is.null(contrast_labels)) {
    NULL
  } else {
    rep(contrast_labels, each = length(measures))
  }
  group <- NULL

  if (!is.null(by)) {
    groups <- c(by$labels, by$em_labels)
    by_measures <- ipw_contrast_names(continuous, collapsible_only = TRUE)
    by_keys <- ipw_contrast_names(continuous, levels, collapsible_only = TRUE)
    group <- c(
      rep(ipw_overall_group, length(effect)),
      rep(groups, each = length(by_keys))
    )
    keys <- c(keys, ipw_by_names(by_keys, groups))
    effect <- c(
      effect,
      rep(
        rep(by_measures, times = length(by_keys) / length(by_measures)),
        times = length(groups)
      )
    )
    if (!is.null(contrast_labels)) {
      contrast <- c(
        contrast,
        rep(
          rep(contrast_labels, each = length(by_measures)),
          times = length(groups)
        )
      )
    }
  }

  list(keys = keys, effect = effect, contrast = contrast, group = group)
}

# The covariance of the reported effects, carried on the estimates table as the
# attribute the shared accessors read. Each effect is a parameter of the stacked
# system, so that covariance is already a block of the stacked one: the block is
# taken across as the sandwich produced it, at the stacked names, and relabeled
# to the labels the table's own rows carry. Rebuilding it from the reported
# standard errors would lose what it is attached for, since the effect measures
# are transformations of one and the same pair of marginal means and covary,
# which a table of standard errors cannot say.
#
# The rows of the block follow the rows of the table. A categorical exposure
# writes its rows contrast-major, all of one contrast's measures before the next
# contrast begins, and names its stacked contrasts level-major, which is the
# same order under two spellings.
attach_effect_covariance <- function(estimates, vcov, keys) {
  labels <- ipw_effect_labels(estimates)
  covariance <- vcov[keys, keys, drop = FALSE]
  dimnames(covariance) <- list(labels, labels)
  attr(estimates, "ipw_vcov") <- covariance
  estimates
}

# The display label of each estimates row: the effect measure, then the contrast
# where a categorical exposure has left the measure repeating across contrasts,
# then the subgroup where a `.by` request has left it repeating across those.
# This is the rule causalgenerics labels the printed rows and the accessor
# output by, restated here so the covariance's dimnames name the effects the
# same way every other surface of the result does. The two rules have to grow
# together: a label built from fewer columns than the table is keyed by repeats
# itself, and a covariance carries that as duplicated dimnames rather than as an
# error.
#
# The subgroup comes last because it qualifies the whole contrast rather than
# one side of it, which is the order causalgenerics reads the identity columns
# in.
#
# Only `attach_effect_covariance()` calls this, and the frame it is handed was
# built moments earlier by `ipw_estimates_from_identity()`, which adds each
# identity column itself. The canonical names are therefore the only ones that
# frame can carry, so no alias needs reading.
ipw_effect_labels <- function(estimates) {
  columns <- intersect(c("effect", "contrast", "group"), names(estimates))
  do.call(paste, unname(lapply(columns, function(column) estimates[[column]])))
}

# The outcome model with its own block of the stacked covariance carried
# alongside it, which is what makes `vcov()` on the stored model report a
# variance that accounts for having estimated the weights. A bare refit of the
# same weighted model treats them as fixed and understates it.
#
# The block is addressed by position rather than by name, since the model's
# coefficient names are not the stack's: the stack prefixes them with `beta_`,
# and a continuous exposure renames the entries its reported effects are read
# from. Position says the same thing for either stack, because the weight
# parameters lead and the outcome-model coefficients follow in both.
wrap_outcome_model <- function(outcome_mod, vcov, weight_parameters) {
  coefficients <- names(stats::coef(outcome_mod))
  block <- weight_parameters + seq_along(coefficients)
  covariance <- vcov[block, block, drop = FALSE]
  dimnames(covariance) <- list(coefficients, coefficients)
  causalgenerics::new_ipw_model(outcome_mod, covariance)
}
