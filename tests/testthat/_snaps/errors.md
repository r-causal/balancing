# balancing_type_error: non-data-frame input

    Code
      expr
    Condition <balancing_type_error>
      Error in `balance()`:
      ! `.data` must be a data frame, not a <list>.

# balancing_selection_error: exposure selects two columns

    Code
      expr
    Condition <balancing_selection_error>
      Error in `balance()`:
      ! `.exposure` must select exactly one column.
      x It selected 2 columns.

# balancing_missing_error: missing covariate values

    Code
      expr
    Condition <balancing_missing_error>
      Error in `balance()`:
      ! `.covariates` must not contain missing values.
      x Missing values in "x1".
      i balancing has no missingness-indicator machinery; complete the covariates first.

# balancing_method_error: a bare-string method

    Code
      expr
    Condition <balancing_method_error>
      Error in `balance()`:
      ! `method` must be a balancing method specification.
      x You supplied a string.
      i Construct one with a method constructor, for example `entropy_balance()`.

# balancing_estimand_error: an unsupported estimand

    Code
      expr
    Condition <balancing_estimand_error>
      Error in `balance()`:
      ! Entropy balancing does not support the "ato" estimand for a "binary" exposure.
      i Supported estimands are "ate", "att", and "atu".

# balancing_estimand_error: a categorical att without focal_level

    Code
      expr
    Condition <balancing_estimand_error>
      Error in `balance()`:
      ! `focal_level` is required for the "att" estimand with a categorical exposure.
      i Supply the exposure level to target, one of "a", "b", and "c".

# balancing_exposure_type_error: a forced type contradicts the data

    Code
      expr
    Condition <balancing_exposure_type_error>
      Error in `balance()`:
      ! `exposure_type` was set to "binary", but the data do not support it.
      x The exposure is detected as "continuous".
      i Drop `exposure_type` to detect it automatically, or supply an exposure of the forced type.

# balancing_constraints_error: quantiles with a continuous exposure

    Code
      expr
    Condition <balancing_constraints_error>
      Error in `build_constraint_matrix()`:
      ! Quantile constraints require a discrete exposure.
      x The exposure is "continuous".
      i Drop `quantiles` from `balance_terms()`, or balance moments instead.

# balancing_constraints_error: an unnamed multi-element tolerance

    Code
      expr
    Condition <balancing_constraints_error>
      Error in `build_constraint_matrix()`:
      ! `tolerance` must be a single number or a named vector.
      x It has length 2 and no names.
      i Supply one value for every covariate, or name each element with a covariate.

# balancing_constraints_error: a tolerance named for a non-covariate

    Code
      expr
    Condition <balancing_constraints_error>
      Error in `build_constraint_matrix()`:
      ! `tolerance` names must be covariates.
      x Not a covariate: "nonesuch".
      i Name each element with one of "x1" and "x2".

# balancing_range_error: base weights of the wrong length

    Code
      expr
    Condition <balancing_range_error>
      Error in `fit_entropy_discrete()`:
      ! `base_weights` must have one value per observation.
      x It has length 3, but the data have 100 rows.

# balancing_range_error: an invalid entropy solver option

    Code
      expr
    Condition <balancing_range_error>
      Error in `resolve_entropy_solver()`:
      ! The `balancing.entropy_solver` option must be one of "newton", "lbfgs", and "lbfgs_then_newton".
      x It is "nope".

# balancing_ipw_unsupported_error: estimating_equations() when absent

    Code
      expr
    Condition <balancing_ipw_unsupported_error>
      Error in `method(estimating_equations, balancing::balancing)`:
      ! This fit has no estimating equations.
      i Estimating equations are produced by the estimating-equation family with exact balance.
      i Use the bootstrap workflow in the inference vignette for variance instead.

