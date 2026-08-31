# balancing_type_error: non-data-frame input

    Code
      balance(list(a = 1), exposure, c(x1, x2), method = bw_entropy())
    Condition <balancing_type_error>
      Error in `balance()`:
      ! `.data` must be a data frame, not a <list>.

# balancing_selection_error: exposure selects two columns

    Code
      balance(data, c(x1, x2), c(x1, x2), method = bw_entropy())
    Condition <balancing_selection_error>
      Error in `balance()`:
      ! `.exposure` must select exactly one column.
      x It selected 2 columns.

# balancing_missing_error: missing covariate values

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy())
    Condition <balancing_missing_error>
      Error in `balance()`:
      ! `.covariates` must not contain missing values.
      x Missing values in "x1".
      i balancing has no missingness-indicator machinery; complete the covariates first.

# balancing_method_error: a bare-string method

    Code
      balance(data, exposure, c(x1, x2), method = "entropy")
    Condition <balancing_method_error>
      Error in `balance()`:
      ! `method` must be a balancing method specification.
      x You supplied a string.
      i Construct one with a method constructor, for example `bw_entropy()`.

# balancing_method_error: an unknown tuning argument

    Code
      bw_entropy(bogus = 1)
    Condition <balancing_method_error>
      Error in `bw_entropy()`:
      ! Unknown tuning argument `bogus`.
      i Check the argument names against the constructor's help page.

# balancing_empty_error: a zero-row data frame

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy())
    Condition <balancing_empty_error>
      Error in `balance()`:
      ! `.data` must have at least two rows.
      x It has no rows.
      i Balancing weights reweight a sample toward a target measured on the sample's own spread, which needs more than one observation.

# balancing_empty_error: a one-row data frame

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy())
    Condition <balancing_empty_error>
      Error in `balance()`:
      ! `.data` must have at least two rows.
      x It has 1 row.
      i Balancing weights reweight a sample toward a target measured on the sample's own spread, which needs more than one observation.

# balancing_estimand_error: an unsupported estimand

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), estimand = "ato")
    Condition <balancing_estimand_error>
      Error in `balance()`:
      ! Entropy balancing does not support the "ato" estimand for a "binary" exposure.
      i Supported estimands are "ate", "att", and "atu".

# balancing_exposure_type_error: a method rejects an exposure type

    Code
      balance(data, exposure, c(x1, x2), method = bw_ipt(), estimand = "ate")
    Condition <balancing_exposure_type_error>
      Error in `balance()`:
      ! Inverse probability tilting does not support a "continuous" exposure.
      i Supported exposure types are "binary" and "categorical".

# balancing_estimand_error: a categorical att without .focal_level

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), estimand = "att")
    Condition <balancing_estimand_error>
      Error in `balance()`:
      ! `.focal_level` is required for the "att" estimand with a categorical exposure.
      i Supply the exposure level to target, one of "a", "b", and "c".

# causalgenerics_forced_exposure_type: forced type contradicts data

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), exposure_type = "binary")
    Condition <causalgenerics_forced_exposure_type>
      Error in `balance()`:
      ! `exposure_type` was set to "binary", but `.exposure` cannot be treated that way.
      x A "binary" exposure takes exactly two observed values, and `.exposure` takes 150.
      i Drop `exposure_type` to detect the type from the data, which reads `.exposure` as "continuous".

# balancing_constraints_error: quantiles with a continuous exposure

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), estimand = "ate",
      constraints = balance_terms(quantiles = 0.5))
    Condition <balancing_constraints_error>
      Error in `balance()`:
      ! Quantile constraints require a discrete exposure.
      x The exposure is "continuous".
      i Drop `quantiles` from `balance_terms()`, or balance moments instead.

# balancing_constraints_error: an unnamed multi-element tolerance

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), constraints = balance_terms(
        tolerance = c(0.1, 0.2)))
    Condition <balancing_constraints_error>
      Error in `balance()`:
      ! `tolerance` must be a single number or a named vector.
      x It has length 2 and no names.
      i Supply one value for every covariate, or name each element with a covariate.

# balancing_constraints_error: a tolerance named for a non-covariate

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), constraints = balance_terms(
        tolerance = c(nonesuch = 0.1)))
    Condition <balancing_constraints_error>
      Error in `balance()`:
      ! `tolerance` names must be covariates.
      x Not a covariate: "nonesuch".
      i Name each element with one of "x1" and "x2".

# balancing_constraints_error: an unnamed multi-element moments vector

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), constraints = balance_terms(
        moments = c(2L, 3L)))
    Condition <balancing_constraints_error>
      Error in `balance()`:
      ! `moments` must be a single whole number or a named vector.
      x It has length 2 and no names.
      i Supply one value for every covariate, or name each element with a covariate.

# balancing_constraints_error: moments named for a non-covariate

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), constraints = balance_terms(
        moments = c(nonesuch = 2L)))
    Condition <balancing_constraints_error>
      Error in `balance()`:
      ! `moments` names must be covariates.
      x Not a covariate: "nonesuch".
      i Name each element with one of "x1" and "x2".

# balancing_constraints_error: an empty constraint set

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), constraints = balance_terms(
        moments = 0L))
    Condition <balancing_constraints_error>
      Error in `balance()`:
      ! Entropy balancing must have at least one balance constraint.
      x The covariates "x1" and "x2" contributed no constraint columns.
      i Raise `moments` in `balance_terms()`, or balance `quantiles` or `interactions` instead.

# balancing_estimand_error: a focal estimand with one exposure level

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), estimand = "att",
      .focal_level = 1)
    Condition <balancing_estimand_error>
      Error in `balance()`:
      ! Balancing needs an exposure with at least two levels.
      x The exposure takes the single level "1".
      i Supply an exposure whose values differ across the sample.

# balancing_estimand_error: a pooled estimand with one exposure level

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), estimand = "ate")
    Condition <balancing_estimand_error>
      Error in `balance()`:
      ! Balancing needs an exposure with at least two levels.
      x The exposure takes the single level "1".
      i Supply an exposure whose values differ across the sample.

# balancing_range_error: base weights of the wrong length

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(base_weights = rep(1, 3)))
    Condition <balancing_range_error>
      Error in `fit_entropy_discrete()`:
      ! `base_weights` must have one value per observation.
      x It has length 3, but the data have 100 rows.

# balancing_range_error: base weights of the wrong length, continuous

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(base_weights = rep(1, 3)))
    Condition <balancing_range_error>
      Error in `fit_entropy_continuous()`:
      ! `base_weights` must have one value per observation.
      x It has length 3, but the data have 100 rows.

# balancing_range_error: infinite covariate values

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy())
    Condition <balancing_range_error>
      Error in `balance()`:
      ! `.covariates` must not contain infinite values.
      x Infinite values in "x1".
      i A fit centers and scales every numeric input, which an infinity leaves undefined.

# balancing_range_error: an infinite continuous exposure

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy())
    Condition <balancing_range_error>
      Error in `balance()`:
      ! `.exposure` must not contain infinite values.
      x Found 1 infinite value.
      i A fit centers and scales every numeric input, which an infinity leaves undefined.

# balancing_range_error: infinite sampling weights

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), sampling_weights = weights)
    Condition <balancing_range_error>
      Error in `balance()`:
      ! `sampling_weights` must not contain infinite values.
      x Found 1 infinite value.
      i A fit centers and scales every numeric input, which an infinity leaves undefined.

# balancing_range_error: sampling weights that are all zero

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), sampling_weights = rep(
        0, nrow(data)))
    Condition <balancing_range_error>
      Error in `balance()`:
      ! `sampling_weights` must not be zero for every observation.
      x Every weight is zero, which leaves no sample to reweight.
      i Individual zero weights are supported; those units are pinned at zero.

# balancing_range_error: an exposure level with no measure

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), sampling_weights = ifelse(
        data$exposure == 1L, 0, 1))
    Condition <balancing_range_error>
      Error in `balance()`:
      ! Every exposure level must carry some base-measure mass.
      x Exposure level "1" has a base measure of zero.
      i The base measure is the sampling weights times any base weights; a level with none has no target to balance to and no total to report at.

# balancing_range_error: disjoint sampling and base weight supports

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(base_weights = 1 -
        alternating), sampling_weights = alternating)
    Condition <balancing_range_error>
      Error in `fit_entropy_discrete()`:
      ! The base measure must carry some mass.
      x It is zero for every observation.
      i The base measure is the sampling weights times any base weights, so vectors that are nonzero nowhere in common leave no sample to reweight.

# balancing_constraints_error: a duplicated moments name

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), constraints = balance_terms(
        moments = c(x1 = 2L, x1 = 3L)))
    Condition <balancing_constraints_error>
      Error in `balance()`:
      ! `moments` names must be unique.
      x "x1" is named more than once.
      i Give each covariate one value.

# balancing_constraints_error: a partially named tolerance

    Code
      balance(data, exposure, c(x1, x2), method = bw_sbw(), constraints = balance_terms(
        tolerance = c(x1 = 0.1, 0.2)))
    Condition <balancing_constraints_error>
      Error in `balance()`:
      ! Every element of `tolerance` must be named when any element is.
      x 1 element carries no name.
      i Name each element with one of "x1" and "x2", or supply a single unnamed value for every covariate.

# balancing_range_error: an invalid entropy solver option

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy())
    Condition <balancing_range_error>
      Error in `resolve_entropy_solver()`:
      ! The `balancing.entropy_solver` option must be one of "newton", "lbfgs", and "lbfgs_then_newton".
      x It is "nope".

# balancing_range_error: an invalid quadratic-program backend option

    Code
      balance(data, exposure, c(x1, x2), method = bw_sbw(), constraints = balance_terms(
        tolerance = 0.05))
    Condition <balancing_range_error>
      Error in `resolve_qp_backend()`:
      ! The `balancing.qp_backend` option must be one of "auto", "osqp", and "clarabel".
      x It is "nope".

# balancing_ipw_unsupported_error: estimating_equations() when absent

    Code
      estimating_equations(fit)
    Condition <balancing_ipw_unsupported_error>
      Error in `method(estimating_equations, balancing::balancing)`:
      ! This fit has no estimating equations.
      i Estimating equations are produced by the estimating-equation family with exact balance.
      i Use the bootstrap workflow in the inference vignette for variance instead.

# balancing_vcov_error: vcov() on a fit with estimating equations

    Code
      stats::vcov(fit)
    Condition <balancing_vcov_error>
      Error in `stats::vcov()`:
      ! This fit carries no covariance for its weight parameters.
      x A covariance for them comes from the stacked system an `ipw()` result assembles, and this fit has not been through one.
      i Build an `ipw()` result from this fit and read `vcov()` off the fit it stores as `wt_mod`.

# balancing_vcov_error: vcov() on a fit without estimating equations

    Code
      stats::vcov(fit)
    Condition <balancing_vcov_error>
      Error in `stats::vcov()`:
      ! This fit carries no covariance for its weight parameters.
      x A covariance for them comes from the stacked system an `ipw()` result assembles, and this fit has not been through one.
      i This fit's weights solve no estimating equations, so there is no such result to read one from.
      i Use the bootstrap workflow in the inference vignette for variance instead.

