# abort() formats messages with cli styling

    Code
      abort(c("The weights did not converge.", i = "Try increasing the number of iterations."),
      error_class = "balancing_convergence_error")
    Condition <balancing_convergence_error>
      Error:
      ! The weights did not converge.
      i Try increasing the number of iterations.

# warn() formats messages with cli styling

    Code
      warn(c("Some weights were negative.", i = "They were set to zero."),
      warning_class = "balancing_negative_weight_warning")
    Condition <balancing_negative_weight_warning>
      Warning:
      Some weights were negative.
      i They were set to zero.

