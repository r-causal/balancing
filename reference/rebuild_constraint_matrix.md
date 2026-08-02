# Rebuild the constraint matrix from a recipe

Rebuild the constraint matrix from a recipe

## Usage

``` r
rebuild_constraint_matrix(recipe, .data)
```

## Arguments

- recipe:

  The covariate expansion recipe stored on a
  [balancing](https://r-causal.github.io/balancing/reference/balancing.md)
  result.

- .data:

  The data frame the recipe was built from.

## Value

The numeric constraint matrix.
