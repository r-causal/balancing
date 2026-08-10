# balancing_ipw_input_error: .by selects no column

    Code
      ipw(fit, outcome_mod, .by = tidyselect::starts_with("no_such_column"))
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `.by` must name exactly one modifier.
      x It names 0 columns.
      i The effects are reported within the levels of a single variable. Cross two variables into one column, with `interaction()`, and name that column instead.

# balancing_ipw_input_error: .by selects two columns

    Code
      ipw(fit, outcome_mod, .by = c(modifier, x1))
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `.by` must name exactly one modifier.
      x It names 2 columns.
      i The effects are reported within the levels of a single variable. Cross two variables into one column, with `interaction()`, and name that column instead.

# balancing_ipw_input_error: .by names a modifier with missing values

    Code
      ipw(fit, outcome_mod, .data = supplied, .by = patchy)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `.by` must name a modifier with no missing values.
      x "patchy" has 3 missing values.
      i A missing value names no subgroup, so the units carrying one belong to none of the strata the effects would be reported within.
      i Drop those rows and refit both models, or recode the missing values as a level of their own.

# balancing_ipw_input_error: .by names a numeric modifier

    Code
      ipw(fit, outcome_mod, .data = data, .by = modifier_hi)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `.by` must name a factor or a character modifier.
      x "modifier_hi" has class <numeric>.
      i The effects are reported within the levels of the modifier, so it has to name a fixed set of subgroups.
      i Cut a continuous column into groups with `cut()`, or convert a logical column or a numeric code to a factor, and name that column instead.

# balancing_ipw_input_error: a stratum holds one exposure level

    Code
      ipw(fit, outcome_mod, .data = supplied, .by = thin)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `.by` must name a modifier whose subgroups each hold every exposure level.
      x "thin = narrow" holds no unit with "exposure" set to "0".
      i An effect within a subgroup contrasts the exposure levels inside it, so a subgroup missing one of them has no contrast to report there.
      i Use a coarser modifier, one whose subgroups each hold every exposure level. Refitting either model does not help, since the data hold no comparison there.

# balancing_ipw_input_error: .data holds the fit's rows out of order

    Code
      ipw(fit, outcome_mod, .data = sorted, .by = modifier)
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `.data` must hold the rows the models were fitted on, in the same order.
      x Its "y" column disagrees with the outcome model frame's at row 1, and at 208 rows in all.
      i The counterfactual predictions, the subgroup indicators, and a focal estimand's standardization are built from `.data`, while the weight equations and the outcome-model score stay in the fit's order.
      i A frame holding the right rows in the wrong order therefore leaves the effect estimates unchanged and every standard error wrong.
      i Supply the frame the outcome model was fitted on, or omit `.data` when the model frame already carries the exposure.

# balancing_ipw_unsupported_error: .by on a continuous exposure

    Code
      ipw(fit, outcome_mod, .by = modifier)
    Condition <balancing_ipw_unsupported_error>
      Error in `ipw()`:
      ! `ipw()` does not support `.by` for a continuous exposure.
      x A continuous exposure reports the marginal structural model's own exposure coefficient rather than a contrast of standardized means, so there is no effect within a subgroup to report.
      i Omit `.by` to report the whole-sample effect.
      i Fitting each subgroup on its own subset reports a coefficient per subgroup and no covariance between them, so the difference between two subgroups cannot be tested from those fits.

# balancing_ipw_by_interaction_warning: no term reads both columns

    Code
      invisible(ipw(fit, outcome_mod, .by = modifier))
    Condition <balancing_ipw_by_interaction_warning>
      Warning in `ipw()`:
      `outcome_mod` has no term reading both "exposure" and "modifier".
      i The subgroup effects are g-computation on `outcome_mod` as it was specified, so the effect differs across subgroups only where a term reads the exposure and the modifier together.
      i Add `exposure:modifier` to `outcome_mod` and refit it to let the effect differ across the levels of "modifier".
      i A model that carries the modification through a column derived from "modifier", such as an indicator, reads that column rather than this one, and its subgroup effects differ even so.

