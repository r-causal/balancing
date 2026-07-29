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
      ! `.data` must have at least one row.
      x It has no rows.
      i Balancing weights require observations to reweight.

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

# balancing_estimand_error: a categorical att without focal_level

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), estimand = "att")
    Condition <balancing_estimand_error>
      Error in `balance()`:
      ! `focal_level` is required for the "att" estimand with a categorical exposure.
      i Supply the exposure level to target, one of "a", "b", and "c".

# balancing_exposure_type_error: a forced type contradicts the data

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(), exposure_type = "binary")
    Condition <balancing_exposure_type_error>
      Error in `balance()`:
      ! `exposure_type` was set to "binary", but the exposure cannot be treated that way.
      x A "binary" exposure takes exactly two distinct values, and this one takes 150.
      i Drop `exposure_type` to detect the type from the data, which reads it as "continuous".

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
      focal_level = 1)
    Condition <balancing_estimand_error>
      Error in `balance()`:
      ! The "att" estimand needs an exposure level outside the focal group.
      x The exposure takes the single level "1".
      i Supply an exposure with at least two levels, or use the "ate" estimand.

# balancing_range_error: base weights of the wrong length

    Code
      balance(data, exposure, c(x1, x2), method = bw_entropy(base_weights = rep(1, 3)))
    Condition <balancing_range_error>
      Error in `fit_entropy_discrete()`:
      ! `base_weights` must have one value per observation.
      x It has length 3, but the data have 100 rows.

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

