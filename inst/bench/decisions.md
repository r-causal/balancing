# Solver and backend decisions

This file records the solver and quadratic-program backend chosen for each
method, with the measured evidence behind the choice and the date it was taken.
Each method exposes its non-default alternatives at runtime (`solver` and
`backend` options), so a user can override any default. A default changes only
through a benchmark showing the alternative solves every instance the current
default solves, meets the same accuracy, and is at least 1.25x faster in the
median with no instance more than 1.5x slower.

Measurement context for every entry below unless noted otherwise:

- Machine: Apple M4 (Mac16,12), 10 cores, 32 GB.
- Toolchain: rustc 1.93.0, R 4.6.0; release build (`lto = "thin"`,
  `codegen-units = 1`).
- Date: 2026-07-14. Commit: 7658e61.

## Estimating-equation family (entropy, inverse probability tilting, just-identified covariate balancing propensity score)

Default: damped Newton. Optional: an L-BFGS-then-Newton hybrid. Demoted: the
basin L-BFGS adapter as a general solver.

Evidence. On identical entropy problems across four regimes, Newton converged in
3 to 4 iterations everywhere, at gradient sup norms at or below 1e-12. The basin
L-BFGS adapter was competitive only on a well-conditioned continuous problem
(0.77x Newton) and failed on a near-collinear instance, running the full 1000
iteration ceiling without satisfying its convergence criterion (final gradient
norm 1.3e-9, against Newton's 3.5e-16 in 3 iterations) at roughly 530x Newton's
time. The hybrid recovered that ill-conditioned instance to machine precision
through its Newton polish (12 iterations) at a small premium over pure Newton,
so it is kept as an option for users who want an L-BFGS approach on well-behaved
problems without giving up the machine-precision estimating equations the
sandwich variance needs. Basin is retained for warm starts and the
over-identified covariate balancing propensity score objective, not as a general
default.

External reference. End to end against WeightIt 1.7.0 at n = 50000 with 200
first-moment constraints, as-shipped versus as-shipped, at equal first-moment
precision (both at or below 3e-11): entropy balancing for the average treatment
effect is 6.73x faster (3.25 s versus 21.9 s); for the average treatment effect
on the treated it is 4.26x faster (2.42 s versus 10.3 s), the smaller margin
reflecting that WeightIt's treated-group-unweighted solve is itself cheaper.

## Energy balancing

Default: osqp. Required, not merely preferred.

Evidence. Energy balancing's quadratic form is indefinite by construction. The
interior-point backend (clarabel) rejects every energy specification up front on
its convexity tag. Forced past the tag on an identical indefinite form, clarabel
terminates at its iteration cap with no solution while osqp solves the same form
to optimality (objective -1195.600385). Interior-point methods require a convex
quadratic; osqp's operator-splitting method does not, and it is the only backend
that solves this family. On the separate convex simplex projection that both
backends accept, clarabel is about 1.1x osqp, so it is not faster even where it
is eligible. Energy solve time grows near-cubically in n through the dense n by n
distance matrix; this is inherent to the method.

## Stable balancing weights

Default: osqp, with an automatic clarabel fallback triggered by an osqp
primal-infeasibility certificate. The backend that produced the weights is
recorded on the fit.

Evidence. Stable balancing weights minimize a diagonal positive quadratic form,
which both backends accept, and on most instances osqp is faster (the binary
average-treatment-effect family runs in 75 to 100 osqp iterations). But
osqp falsely certifies primal infeasibility on some feasible instances: binary
average-treatment-effect problems at n = 20000 and ill-scaled constraint
matrices, where clarabel solves to a constraint violation at or below 1e-15 in 7
to 9 interior-point iterations. The false certificate is driven by problem shape
and conditioning, not size alone, so it cannot be routed around by a size
threshold; the only reliable trigger is the certificate itself. Blanket
clarabel is not a promotable default (median 1.35x osqp over the feasible-both
grid, nine of 27 instances above 1.5x). The fallback keeps osqp's speed on the
happy path and clarabel's correctness on the instances osqp gets wrong. Verified
at n = 20000, binary average treatment effect: the default path
resolves and records `clarabel_fallback`, an explicit osqp pin errors, and an
explicit clarabel pin solves.

## Characteristic function distance balancing (kernel balancing)

Default: osqp, retained provisionally. The routing is under review: current
measurements do not support an osqp speed advantage on the gaussian kernel and
in fact favor clarabel at small sizes. A change is deferred until the full kernel
family is measured.

Evidence. On the gaussian kernel at matched tolerance, with objectives
cross-checked before timing, clarabel was faster than osqp at n = 500 (0.50x)
and n = 1000 (0.57x), converging in 11 to 12 interior-point iterations against
1550 to 1725 osqp iterations, and reached parity at n = 2000
(1.04x) as its dense factorization cost caught up. The two backends reached the
same objective (relative gaps at or below 6e-8) but produced very different
weight vectors: the kernel-only program with a small weight ridge is weakly
convex near the minimum-weight floor and admits many near-optimal weight
vectors, so the two backends select different ones at equal balance quality. A
routing change therefore needs the full kernel family (gaussian, matern,
laplace, t, and energy) across sizes and estimands, and a decision on whether
the weight non-uniqueness is acceptable, before osqp is displaced. Until then
osqp remains the default and clarabel is available through the `backend` option.

The energy kernel reproduces energy balancing exactly (identical weights,
correlation 1.000000) with at or below 2.5 percent overhead over the direct
energy method, confirming the reference equivalence the kernel method documents.
