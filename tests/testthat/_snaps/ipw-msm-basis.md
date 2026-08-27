# a basis fit announces the reading it records

    Code
      invisible(ipw(fit, outcome_mod))
    Message
      i `ipw()` reports only the conditional reading when the exposure enters `outcome_mod` through several columns.
      i The coefficient surface is the outcome model's own, and no single row of it is a causal effect.
      i Marginalizing over the dose is left to the marginaleffects package: call `avg_slopes()` or `avg_comparisons()` on the conditional result. See <https://marginaleffects.com/chapters/interactions.html>.
      i Set `effects = "conditional"` to silence this message.

# a basis fit refuses the marginal reading at construction

    Code
      ipw(fit, outcome_mod, effects = "marginal")
    Condition <balancing_ipw_input_error>
      Error in `ipw()`:
      ! `ipw()` reports only the conditional reading when the exposure enters `outcome_mod` through several columns.
      x `effects = "marginal"` asks for a reading this model has none of.
      i The coefficient surface is the outcome model's own, and no single row of it is a causal effect.
      i Marginalizing over the dose is left to the marginaleffects package: call `avg_slopes()` or `avg_comparisons()` on the conditional result. See <https://marginaleffects.com/chapters/interactions.html>.
      i Use `effects = "conditional"` or omit `effects`.

