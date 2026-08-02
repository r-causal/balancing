# The estimating-equations container

`balancing_estimating_equations` holds the pieces a stacked sandwich
variance needs after balancing: the fitted parameters, the per-unit
estimating functions, the analytic Jacobian, and each unit's weight
derivatives. It is produced by methods whose weights solve smooth
estimating equations and is `NULL` otherwise.

## Usage

``` r
balancing_estimating_equations(
  parameters = numeric(0),
  psi = numeric(0),
  jacobian = numeric(0),
  weight_jacobian = numeric(0),
  weights_raw = NULL,
  psi_fn = NULL,
  weights_fn = NULL,
  parts_fn = NULL
)
```

## Arguments

- parameters:

  The fitted parameters, a numeric vector.

- psi:

  The `n` by `p` estimating functions at the solution.

- jacobian:

  The `p` by `p` analytic Jacobian at the solution.

- weight_jacobian:

  The `n` by `p` weight derivatives.

- weights_raw:

  The balancing weights whose derivative is `weight_jacobian`, a
  length-`n` numeric vector. The weight derivatives are stored at
  whatever per-group reporting scale a method uses internally, so a
  consumer that needs the derivative of the reported weights rescales
  `weight_jacobian` by the ratio of the reported weights to
  `weights_raw`.

- psi_fn:

  An optional function re-evaluating `psi` at new parameters.

- weights_fn:

  An optional function returning the reported balancing weights at new
  parameters, a plain double vector with the sampling weights excluded.
  At the fitted parameters it reproduces
  `as.numeric(weights(fit, include_sampling_weights = FALSE))`. The
  per-group reporting scale is fixed at the fit rather than recomputed
  at each set of parameters, so the function's derivative is
  `weight_jacobian` rescaled by the ratio of the reported weights to
  `weights_raw`, which is the weight coupling a stacked variance needs.

- parts_fn:

  An optional function returning both of the above at one set of
  parameters, a list with elements `weights` and `psi` holding exactly
  what `weights_fn` and `psi_fn` return there. A method whose estimating
  functions are a transformation of its own weights computes the pair
  together for the price of one, and a consumer that needs both at every
  parameter vector, as a stacked sandwich does, halves its work by
  asking for them together. It is an optimization rather than a
  contract: a consumer reads it when it is present and falls back to the
  two functions when it is not, so a method supplies it only when the
  saving is real.

## Value

A `balancing_estimating_equations` object.
