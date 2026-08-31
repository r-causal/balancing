# Balance constraints

`balance_terms()` records the set of covariate functions a balancing
method should equate across exposure groups. It is passed to
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md)
through `constraints`. The specification is data-free: the covariate
expansion it describes is applied at fit time.

## Usage

``` r
balance_terms(
  moments = NULL,
  interactions = FALSE,
  quantiles = NULL,
  tolerance = 0,
  ...
)
```

## Arguments

- moments:

  The highest covariate power to balance. A single whole number or a
  named integer vector; `NULL` (the default) resolves to first moments.
  For a continuous exposure each power is held at zero weighted
  correlation with the exposure instead.

- interactions:

  Whether to add pairwise interactions of the base columns. These expand
  the constraint set the weights must balance, adding the pairwise
  products of the base columns to the covariate functions a fit
  constrains: equated across the exposure groups of a discrete exposure,
  and driven to zero correlation with a continuous one. They say nothing
  about causal interaction between two exposures, which is an effect
  rather than a constraint and is reported for a joint exposure by
  [`ipw()`](https://r-causal.github.io/balancing/reference/ipw.balancing.md).

- quantiles:

  Quantile probabilities in `(0, 1)`: a numeric vector applied to every
  continuous covariate, a named list of probabilities per covariate, or
  `NULL` for none.

- tolerance:

  The per-constraint tolerance: a single non-negative number, or a named
  vector giving the tolerance per source covariate.

- ...:

  Reserved for future extensions; must be empty.

## Value

A `balance_terms` object.

## Details

The constraint set is built from four ingredients:

- `moments`: the highest power of each numeric covariate to balance. A
  scalar applies to every covariate; a named integer vector sets powers
  per covariate, with unnamed covariates defaulting to `1`. Powers above
  `1` are ignored for binary indicator columns.

- `interactions`: when `TRUE`, all pairwise products of distinct base
  columns are added, excluding products of two indicators of the same
  factor.

- `quantiles`: probabilities in `(0, 1)`. Each probability adds an
  indicator column so that mean balance on the indicator is quantile
  balance on the covariate. A single probability vector applies to every
  continuous covariate; a named list sets probabilities per covariate.
  Quantile constraints apply to discrete exposures only.

- `tolerance`: the largest absolute standardized mean difference
  (discrete exposures) or exposure-covariate correlation (continuous
  exposures) permitted per constraint. A scalar applies to every
  covariate; a named vector sets tolerances per source covariate, and
  derived columns inherit their source covariate's tolerance. `0`
  requests exact balance. A positive value selects the inexact problem
  for entropy balancing and is the central tuning parameter for stable
  balancing weights.

For a continuous exposure there are no groups to equate, so a constraint
column is instead held within `tolerance` of zero weighted correlation
with the exposure. On the continuous energy path this is what `moments`
and `interactions` request, and there the correlation is held exactly
whatever `tolerance` says, for the reason
[`bw_energy()`](https://r-causal.github.io/balancing/reference/bw_energy.md)
records. The marginal distribution of the exposure and of the covariates
is a separate matter, held by `distribution_moments` in
[`bw_energy()`](https://r-causal.github.io/balancing/reference/bw_energy.md)
and
[`bw_entropy()`](https://r-causal.github.io/balancing/reference/bw_entropy.md);
asking for correlation constraints does not add marginal rows, and
raising `distribution_moments` adds no correlation constraint.

A factor covariate contributes one indicator per level rather than the
reference coding a model formula would use. Those indicators sum to the
constant every balancing method carries, so one of them is redundant and
the expansion drops it, naming the term in an informational alert. The
level dropped is the last one, and the constraints that remain balance
it as well: with the other level proportions equated across exposure
groups, the omitted one follows. The balance table reports the surviving
levels.

## Examples

``` r
# Balance means and variances of every numeric covariate.
balance_terms(moments = 2)
#> <balancing::balance_terms>
#>  @ moments     : int 2
#>  @ interactions: logi FALSE
#>  @ quantiles   : NULL
#>  @ tolerance   : num 0

# Balance means with a relaxed tolerance.
balance_terms(tolerance = 0.05)
#> <balancing::balance_terms>
#>  @ moments     : NULL
#>  @ interactions: logi FALSE
#>  @ quantiles   : NULL
#>  @ tolerance   : num 0.05
```
