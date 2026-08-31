# Characteristic function distance balancing

`bw_cfd()` specifies characteristic function distance balancing, also
called kernel balancing, for
[`balance()`](https://r-causal.github.io/balancing/reference/balance.md).
The weights minimize a kernel measure of the distance between the
reweighted exposure groups and a target sample, following Wong and Chan.
Every kernel but `"energy"` is positive semidefinite, so the objective
is a convex quadratic program; the energy kernel is conditionally
positive semidefinite alone, so its quadratic form is indefinite and the
solve routes to the ADMM backend. Either way the reweighting improves
multivariate covariate balance without positing a propensity model.
Characteristic function distance balancing supports binary and
categorical exposures.

## Usage

``` r
bw_cfd(
  ...,
  kernel = c("gaussian", "matern", "laplace", "t", "energy"),
  smoothness = 1.5,
  degrees_of_freedom = 5,
  simulation_draws = 5000,
  improved = TRUE,
  weight_penalty = 1e-04,
  min_weight = 1e-08,
  convergence_tolerance = NULL,
  max_iterations = NULL
)
```

## Arguments

- ...:

  Reserved for future extensions; must be empty. Tuning parameters must
  be passed by name.

- kernel:

  The covariate kernel the objective is built on, one of `"gaussian"`,
  `"matern"`, `"laplace"`, `"t"`, or `"energy"`.

- smoothness:

  The Matern smoothness order, one of `0.5`, `1.5`, or `2.5`. Used only
  by the Matern kernel.

- degrees_of_freedom:

  The degrees of freedom of the t kernel's frequency distribution, a
  finite number greater than two. Used only by the `"t"` kernel. The t
  kernel approaches the Gaussian kernel as the degrees of freedom grow,
  so that limit is requested as `kernel = "gaussian"` rather than as an
  infinite degrees of freedom, which names no frequency distribution to
  draw from.

- simulation_draws:

  The number of Monte Carlo frequency projections for the `"t"` kernel.

- improved:

  Whether to add the between-group term of the improved variant for the
  average treatment effect with a discrete exposure.

- weight_penalty:

  The L2 penalty on the weights, which stabilizes the quadratic program.

- min_weight:

  The smallest permitted weight.

- convergence_tolerance:

  The quadratic-program solver tolerance, or `NULL` for the resolved
  default of `1e-8`, which the solver applies as both its absolute and
  its relative tolerance. Under the alternating-direction backend, which
  the other kernels take by default and which the `"energy"` kernel
  always takes because its quadratic form is indefinite, a tolerance
  below what the problem can reach spends the full iteration cap and
  then warns.

- max_iterations:

  The maximum solver iterations, or `NULL` for the resolved default of
  200000.

## Value

A `bw_cfd` specification, a
[balance_method](https://r-causal.github.io/balancing/reference/balance_method.md).

## Details

The objective is the sum of each group's kernel distance to the target
sample, built from a covariate kernel matrix rather than a pairwise
distance matrix. `kernel` selects the kernel. The `"gaussian"`,
`"laplace"`, and `"matern"` kernels are distance based: their bandwidth
is the median of the pairwise covariate distances, so the fit is
invariant to a uniform rescaling of the covariates. `smoothness` sets
the Matern order, one of `0.5`, `1.5`, or `2.5`. The `"t"` kernel
approximates a heavy-tailed spectral kernel by Monte Carlo:
`simulation_draws` frequency vectors are drawn from a multivariate t
distribution with `degrees_of_freedom` degrees of freedom, on the R side
under R's random number generator, so a fixed seed reproduces the
weights. Raising the degrees of freedom takes the t kernel toward the
Gaussian kernel, which is available exactly as `kernel = "gaussian"` and
needs no Monte Carlo draws. The improved variant for the average
treatment effect adds the between-group term of the kernel mean
embedding, balancing the groups against one another as well as against
the sample.

Setting `kernel = "energy"` uses the negative pairwise distance, which
reproduces
[`bw_energy()`](https://r-causal.github.io/balancing/reference/bw_energy.md)
with its `"scaled_euclidean"` distance on the same data and constraints.
That equivalence is the reference point for the method, since the energy
kernel has an external implementation while the other kernels do not.

Characteristic function distance balancing belongs to the
quadratic-program family, which has no estimating equations, so a fit
produces no estimating-equations container and the guarantee is
objective-level rather than exact moment balance: without moment
constraints the kernel objective drives balance, and with them the
constraint rows hold within tolerance. The tolerance in
[`balance_terms()`](https://r-causal.github.io/balancing/reference/balance_terms.md)
relaxes any added moment constraints rather than selecting an inexact
solver. A tolerance supplied without moment constraints has nothing to
relax, so it is warned and ignored.

What the optimization pins is the objective and the balance it achieves,
not the individual weights. A kernel matrix built on a smooth kernel is
ill-conditioned by construction, since nearby units contribute nearly
proportional rows, and a quadratic program whose objective is flat along
such directions has a set of solutions rather than one. Two
quadratic-program backends given the same ill-conditioned problem agree
on the objective to seven significant figures while individual unit
weights differ materially; a Gaussian kernel at its median bandwidth is
the case this was measured on. Treat the per-unit weights of a
`bw_cfd()` fit as one member of an equivalence set rather than as a
uniquely determined quantity: what the method determines is the
reweighted distribution, so estimands and balance statistics computed
from the weights are stable where a claim about a particular unit's
weight is not.

## References

Wong, R. K. W. and Chan, K. C. G. (2018). Kernel-based covariate
functional balancing for observational studies. *Biometrika*, 105(1),
199-213.

Huling, J. D. and Mak, S. (2024). Energy balancing of covariate
distributions. *Journal of Causal Inference*, 12(1), 20220029.

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
fit <- balance(df, exposure, c(x1, x2), method = bw_cfd())
#> ℹ Treating `.exposure` as binary
fit
#> 
#> ── Characteristic function distance balancing ──────────────────────────────────
#> Exposure: "exposure" (binary)
#> Estimand: "ate"
#> Observations: 200
#> Solver: converged in 1825 iterations
#> Constraints: 2 terms (tolerance 0)
#> Largest imbalance: 0.0117 (standardized mean difference)
```
