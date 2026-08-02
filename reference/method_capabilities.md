# Method capability generics

These generics let
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
interrogate a method specification without knowing its concrete class.
Each balancing method supplies one method per generic.

## Usage

``` r
supported_exposure_types(method, ...)

supported_estimands(method, exposure_type)

supports_estimating_equations(method, ...)

method_label(method, ...)

fit_method(method, prepared)
```

## Arguments

- method:

  A
  [balance_method](https://r-causal.github.io/balancing/reference/balance_method.md)
  specification.

- ...:

  Additional context for `supports_estimating_equations()`, passed by
  name. Every method accepts both `exposure_type` and `constraints` and
  ignores whichever of the two its answer does not depend on, so one
  call shape puts the question to any method.

- exposure_type:

  One of `"binary"`, `"categorical"`, or `"continuous"`.

## Value

`supported_exposure_types()` and `supported_estimands()` return
character vectors; `supports_estimating_equations()` returns a single
logical; `method_label()` returns a single string; `fit_method()`
returns the internal solve result.
