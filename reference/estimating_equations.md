# Extract the estimating-equations container

`estimating_equations()` returns the
[balancing_estimating_equations](https://r-causal.github.io/balancing/reference/balancing_estimating_equations.md)
a fit produced, the pieces a stacked sandwich variance needs after
balancing. It is available only for fits whose weights solve smooth
estimating equations: the estimating-equation family (entropy balancing,
inverse probability tilting, just-identified covariate balancing
propensity score) with exact balance. Fits without estimating equations,
such as any tolerance-relaxed or quadratic-program fit, raise
`balancing_ipw_unsupported_error`.

## Usage

``` r
estimating_equations(x, ...)
```

## Arguments

- x:

  A
  [balancing](https://r-causal.github.io/balancing/reference/balancing.md)
  result.

- ...:

  Ignored.

## Value

A
[balancing_estimating_equations](https://r-causal.github.io/balancing/reference/balancing_estimating_equations.md)
object.

## Examples

``` r
n <- 200
x1 <- rnorm(n)
df <- data.frame(exposure = rbinom(n, 1, plogis(0.5 * x1)), x1 = x1)
fit <- balance(df, exposure, x1, method = bw_entropy())
#> ℹ Treating `.exposure` as binary.
estimating_equations(fit)
#> <balancing::balancing_estimating_equations>
#>  @ parameters     : num [1:2] -0.0821 0.0922
#>  @ psi            : num [1:200, 1:2] 0 0.64 -1.78 0 0 ...
#>  @ jacobian       : num [1:2, 1:2] -106.3 0 0 -95.8
#>  @ weight_jacobian: num [1:200, 1:2] 0 -0.64 1.78 0 0 ...
#>  @ weights_raw    : num [1:200] 1.181 1.055 0.843 1.013 1.002 ...
#>  @ psi_fn         : function (theta)  
#>  @ weights_fn     : function (theta)  
#>  @ parts_fn       : function (theta)  
```
