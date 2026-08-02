# The balancing method specification classes

The method constructors
([`bw_entropy()`](https://r-causal.github.io/balancing/reference/bw_entropy.md)
and its siblings) return small S7 objects that carry tuning parameters
and never touch data. They share an abstract hierarchy so that
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
can query a method's capabilities uniformly. `balance_method` is the
abstract root. The estimating-equation family (entropy balancing,
inverse probability tilting, the covariate balancing propensity score)
subclasses `estimating_equation_method`; the quadratic-program family
(energy, characteristic function distance, and stable balancing weights)
subclasses `quadratic_program_method`. None of the abstract classes can
be constructed directly.

## Usage

``` r
balance_method(convergence_tolerance = NULL, max_iterations = NULL)

estimating_equation_method(convergence_tolerance = NULL, max_iterations = NULL)

quadratic_program_method(
  convergence_tolerance = NULL,
  max_iterations = NULL,
  weight_penalty = numeric(0),
  min_weight = numeric(0)
)
```

## Arguments

- convergence_tolerance:

  The solver convergence tolerance, or `NULL` for the core default.

- max_iterations:

  The maximum solver iterations, or `NULL` for the core default.

- weight_penalty:

  The L2 penalty on the weights.

- min_weight:

  The smallest permitted weight.

## Value

An abstract class object. Constructing a concrete subclass returns a
`balance_method`.
