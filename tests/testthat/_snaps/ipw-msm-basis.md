# a basis fit announces the reading it records

    Code
      invisible(ipw(fit, outcome_mod))
    Message
      i `ipw()` reports only the conditional reading because the exposure enters `outcome_mod` through more than one term, such as a spline or polynomial.
      i With a nonlinear dose-response, no single coefficient is the effect of the exposure, so there is no marginal effect to report.
      i Use the marginaleffects package to marginalize over the dose: `avg_slopes()` for slopes, `avg_comparisons()` for contrasts, and `avg_predictions()` for causal dose-response functions. See <https://marginaleffects.com/chapters/interactions.html>.
      i Set `effects = "conditional"` to silence this message.

# a basis fit refuses the marginal reading at construction

    Code
      ipw(fit, outcome_mod, effects = "marginal")
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `ipw()` reports only the conditional reading because the exposure enters `outcome_mod` through more than one term, such as a spline or polynomial.
      x `effects = "marginal"` is not available for this model.
      i With a nonlinear dose-response, no single coefficient is the effect of the exposure, so there is no marginal effect to report.
      i Use the marginaleffects package to marginalize over the dose: `avg_slopes()` for slopes, `avg_comparisons()` for contrasts, and `avg_predictions()` for causal dose-response functions. See <https://marginaleffects.com/chapters/interactions.html>.
      i Use `effects = "conditional"` or omit `effects`.

