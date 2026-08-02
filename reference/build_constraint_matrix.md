# Build the constraint matrix and recipe

Build the constraint matrix and recipe

## Usage

``` r
build_constraint_matrix(
  .data,
  .covariates,
  constraints,
  exposure_type,
  sampling_weights = NULL,
  call = rlang::caller_env()
)
```

## Arguments

- .data:

  The data frame.

- .covariates:

  A character vector of covariate column names.

- constraints:

  A
  [balance_terms](https://r-causal.github.io/balancing/reference/balance_terms.md)
  specification, or `NULL` for the default first-moment constraints.

- exposure_type:

  One of `"binary"`, `"categorical"`, or `"continuous"`.

- sampling_weights:

  Optional sampling weights. When supplied, numeric columns are
  standardized to weighted mean zero and unit weighted standard
  deviation rather than the unweighted sample scale.

- call:

  The calling environment, used to build the error's call so a
  constraint error names the user-facing function.

## Value

A list with `matrix` and `recipe`.
