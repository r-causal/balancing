# Extract balancing weights

Returns the weights a
[balancing](https://r-causal.github.io/balancing/reference/balancing.md)
result fitted, as a
[bw](https://r-causal.github.io/balancing/reference/bw.md) vector
carrying the estimand they target. When sampling weights are present
they are composed onto the balancing weights by default.

## Arguments

- object:

  A
  [balancing](https://r-causal.github.io/balancing/reference/balancing.md)
  result.

- ...:

  Ignored.

- include_sampling_weights:

  Whether to multiply the balancing weights by the sampling weights.
  Defaults to `TRUE`.

## Value

A [bw](https://r-causal.github.io/balancing/reference/bw.md) vector of
the same length as the data the fit was built from, in the same row
order. It carries the estimand and no other metadata about the fit.
