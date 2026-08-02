# balancing: Optimization-Based Balancing Weights for Causal Inference

Calculates optimization-based balancing weights for causal inference,
including entropy balancing, inverse probability tilting, the covariate
balancing propensity score, energy balancing, characteristic function
distance balancing, and stable balancing weights. Weights target
covariate balance directly by solving a convex optimization problem,
with support for binary, categorical, and continuous exposures across a
range of causal estimands. A Rust core provides the numerical solvers.

## See also

Useful links:

- <https://r-causal.github.io/balancing/>

- <https://github.com/r-causal/balancing>

- Report bugs at <https://github.com/r-causal/balancing/issues>

## Author

**Maintainer**: Malcolm Barrett <malcolmbarrett@gmail.com>
([ORCID](https://orcid.org/0000-0003-0299-5825)) \[copyright holder\]

Authors:

- Malcolm Barrett <malcolmbarrett@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-0299-5825)) \[copyright holder\]

Other contributors:

- The authors of the vendored Rust crates (see src/rust/Cargo.lock and
  the bundled crate sources for the full list) \[copyright holder\]
