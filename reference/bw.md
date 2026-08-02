# Balancing weight vectors

`bw` objects are numeric vectors carrying the estimand a set of
balancing weights targets. They are the weight vector stored on a
[balancing](https://r-causal.github.io/balancing/reference/balancing.md)
result and returned by its
[`weights()`](https://rdrr.io/r/stats/weights.html) method.

`bw` is a sibling of
[propensity::psw](https://r-causal.github.io/propensity/reference/psw.html):
both inherit the `causal_wts` class, so
[`causalgenerics::is_causal_wt()`](https://r-causal.github.io/causalgenerics/reference/causal-weights.html)
and
[`causalgenerics::estimand()`](https://r-causal.github.io/causalgenerics/reference/causal-weights.html)
work on either.

## Usage

``` r
new_bw(x = double(), estimand = NULL, ...)

bw(x = double(), estimand = NULL)

as_bw(x, estimand = NULL)

is_bw(x)
```

## Arguments

- x:

  A numeric vector of weights for `bw()` and `new_bw()`, or an object to
  test or coerce for `is_bw()` and `as_bw()`.

- estimand:

  A single string naming the target estimand, or `NULL`.

- ...:

  Additional attributes stored on the object (developer use only).

## Value

- `new_bw()`, `bw()`, `as_bw()`: a `bw` vector.

- `is_bw()`: a single logical value.

## Details

### Constructors

- `bw()` is the user-facing constructor. It coerces `x` to double before
  building the object.

- `new_bw()` is the low-level constructor for developers. It assumes `x`
  is already a double vector.

- `as_bw()` coerces an existing numeric vector to a `bw` object.

### Queries

- `is_bw()` tests whether an object is a `bw` vector.

- [`causalgenerics::estimand()`](https://r-causal.github.io/causalgenerics/reference/causal-weights.html)
  reads the target estimand.

### Combining

Arithmetic preserves the class and estimand, so normalizing weights
keeps the metadata. Combining `bw` vectors with matching estimands
preserves the class; combining vectors with different estimands, or a
`bw` with a
[propensity::psw](https://r-causal.github.io/propensity/reference/psw.html),
warns and falls back to a plain double vector.

## Examples

``` r
w <- bw(c(0.5, 1, 1.5), estimand = "ate")
w
#> <bw{estimand = ate}[3]>
#> [1] 0.5 1.0 1.5
is_bw(w)
#> [1] TRUE
estimand(w)
#> [1] "ate"

# Arithmetic preserves the class.
w / sum(w)
#> <bw{estimand = ate}[3]>
#> [1] 0.1666667 0.3333333 0.5000000
```
