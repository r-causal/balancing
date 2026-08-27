# Estimate balancing weights

`balance()` fits a balancing method to a data frame, returning weights
that target covariate balance directly. The exposure and covariates are
chosen with tidyselect, the method is one of the method specifications
such as
[`bw_entropy()`](https://r-causal.github.io/balancing/reference/bw_entropy.md),
and the estimand and constraints control what balance the weights
achieve.

## Usage

``` r
balance(
  .data,
  .exposure,
  .covariates,
  method = bw_entropy(),
  estimand = c("ate", "att", "atc", "ato"),
  ...,
  constraints = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  focal_level = NULL,
  sampling_weights = NULL
)
```

## Arguments

- .data:

  A data frame.

- .exposure:

  The exposure column, selected with data-masking. Exactly one column.

- .covariates:

  The covariate columns, selected with tidyselect. At least one column,
  with no default.

- method:

  A
  [balance_method](https://r-causal.github.io/balancing/reference/balance_method.md)
  specification from one of the method constructors, such as
  [`bw_entropy()`](https://r-causal.github.io/balancing/reference/bw_entropy.md).

- estimand:

  The target estimand: `"ate"`, `"att"`, `"atc"` (stored as `"atu"`), or
  `"ato"`. Defaults to `"ate"`.

- ...:

  Reserved; must be empty.

- constraints:

  A
  [`balance_terms()`](https://r-causal.github.io/balancing/reference/balance_terms.md)
  specification, or `NULL` for the method default.

- exposure_type:

  One of `"auto"` (the default), `"binary"`, `"categorical"`, or
  `"continuous"`.

- focal_level:

  The focal exposure level for `"att"` and `"atc"`. Inferred for a
  binary exposure; required for a categorical exposure.

- sampling_weights:

  Sampling weights, given as a bare column name or an external numeric
  vector, or `NULL`.

## Value

A
[balancing](https://r-causal.github.io/balancing/reference/balancing.md)
object.

## Details

The exposure type is detected automatically and announced through an
informational message, which `options(balancing.quiet = TRUE)`
suppresses. The estimand vocabulary matches propensity: `"atc"` is
accepted as a synonym for the untreated target and stored as `"atu"`.
`"att"` and `"atc"` reweight toward a focal exposure level, inferred for
a binary exposure and required through `focal_level` for a categorical
exposure. Continuous exposures permit only `"ate"`.

Constraints default to first-moment balance. Pass a
[`balance_terms()`](https://r-causal.github.io/balancing/reference/balance_terms.md)
specification to balance higher moments, interactions, or quantiles, or
to relax exact balance to a tolerance.

A factor covariate expands to one indicator per level, and those
indicators sum to the constant every balancing method carries. One
indicator per factor is therefore redundant with that constant and is
dropped, with an informational alert naming the term. The dropped level
is the last one, and balancing the levels that remain balances it too.
The factor stays in `@covariates`, and `@balance_table` reports the
surviving levels rather than the full set.

## Examples

``` r
n <- 200
x1 <- rnorm(n)
x2 <- rnorm(n)
df <- data.frame(
  exposure = rbinom(n, 1, plogis(0.5 * x1 - 0.5 * x2)),
  x1 = x1,
  x2 = x2
)
fit <- balance(df, exposure, c(x1, x2), method = bw_entropy())
#> ℹ Treating `.exposure` as binary
fit
#> 
#> ── Entropy balancing ───────────────────────────────────────────────────────────
#> Exposure: "exposure" (binary)
#> Estimand: "ate"
#> Observations: 200
#> Solver: converged in 4 iterations
#> Constraints: 2 terms (tolerance 0)
#> Largest imbalance: 0.0000 (standardized mean difference)
weights(fit)
#> <bw{estimand = ate}[200]>
#>   [1] 1.4822934 0.7994941 0.5217720 1.1479092 0.4533513 0.7719240 1.4035649
#>   [8] 0.9802867 0.6591079 0.4800601 0.9333154 0.8102705 1.6120219 1.2614811
#>  [15] 1.0726080 0.5965857 0.7849825 0.9401063 0.6117172 0.6289439 0.8520625
#>  [22] 0.9468626 0.7128332 0.8910249 0.4571927 0.7437911 0.7347726 1.2079409
#>  [29] 1.7007547 1.1320231 1.3291447 0.7122425 1.2337525 0.4245500 1.0657892
#>  [36] 0.9448055 1.5472549 0.8184924 0.7120407 2.3243778 1.3922519 1.8993094
#>  [43] 0.7821393 0.7156895 1.7979486 0.6757129 1.0673768 0.8853535 0.9829299
#>  [50] 0.5794024 1.7340581 1.1372705 0.8287521 1.3162716 0.7278441 0.7254614
#>  [57] 0.4972972 0.4126609 1.1278622 1.5171804 1.4193092 0.8940481 0.7640096
#>  [64] 1.0504507 1.2124434 0.8180128 1.3138355 0.7945054 0.9039191 0.5995180
#>  [71] 0.6741742 0.5422249 0.6920180 1.1133355 0.7207668 1.5469144 0.7167294
#>  [78] 1.3624360 0.9150708 1.1298465 0.7888205 0.5610904 1.1505776 0.9957281
#>  [85] 0.8201578 1.1014723 0.6259533 1.2059720 1.3995685 1.0559207 0.9040765
#>  [92] 0.7320147 0.5912157 1.2506114 0.8025971 0.4236445 1.4067153 1.1222699
#>  [99] 1.3436226 0.5525230 1.1053023 1.0277711 0.5807392 1.9199148 0.5959668
#> [106] 0.9663253 0.7643925 0.5334122 0.6881712 2.2038752 0.9940144 1.2074641
#> [113] 0.9934255 0.7497418 1.0115677 0.5664137 1.3148216 1.6978816 0.7475722
#> [120] 0.8199103 1.6592912 1.1152774 1.0917257 0.9671046 1.8852798 1.6019299
#> [127] 1.0046930 0.6775991 0.8265485 0.9489581 0.5690831 2.3531584 0.7904469
#> [134] 0.7541519 0.7417943 0.8708424 1.1211099 0.8646106 1.3589596 1.3417948
#> [141] 1.0663867 0.4739887 0.5866959 0.5802799 0.5821500 1.1516238 0.7786035
#> [148] 1.4043011 0.6995948 0.8366562 0.8069833 1.1529196 0.6574264 0.5655263
#> [155] 0.9931426 1.3352181 1.1989259 2.0035946 0.7287538 1.7205435 0.6605466
#> [162] 1.0639135 1.4911738 0.7110839 1.8300625 0.6460642 1.2048766 0.7678042
#> [169] 1.6031276 0.8783361 0.9869689 0.9613194 0.5594809 0.6132884 1.4764363
#> [176] 1.0704189 1.4017109 0.8254835 1.0099264 0.7658067 0.7131065 1.4019112
#> [183] 0.6616522 0.6809120 1.0140668 0.5869015 1.0177748 0.6325730 0.6441757
#> [190] 1.3773450 0.6814893 1.5846447 0.5071859 1.0649916 0.6279650 1.8233348
#> [197] 0.6641230 1.1021803 1.6417696 1.4233126
```
