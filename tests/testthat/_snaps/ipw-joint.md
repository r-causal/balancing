# balancing_ipw_unsupported_error: .by on a declared crossing

    Code
      ipw(fit, outcome_mod, .by = modifier)
    Condition <balancing_ipw_unsupported_error>
      Error in `ipw()`:
      ! `ipw()` does not support `.by` for a joint exposure.
      x A joint exposure already reports an interaction between two treatments, and reporting it again within the levels of a modifier is a three-way question this surface does not answer.
      i Drop the declaration with `factor(x)` to report each cell against the reference cell within the levels of `.by`.

# balancing_ipw_unsupported_error: a focal estimand on a crossing

    Code
      ipw(fit, outcome_mod)
    Condition <balancing_ipw_unsupported_error>
      Error in `ipw()`:
      ! `ipw()` reports a joint exposure for the "ate" estimand only.
      x The weights were built for the "att" estimand.
      i Every cell mean on the joint surface standardizes to one population, and a tilted estimand standardizes each of them to a population the simple effects and the interaction are not defined over.
      i Refit the weights for "ate", or drop the declaration with `factor(x)` to report each cell against the reference cell under the estimand you have.

