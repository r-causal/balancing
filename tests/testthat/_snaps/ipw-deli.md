# balancing_internal_error: a matrix block of the wrong width

    Code
      stack_psi_blocks(list(wide_enough, too_narrow), n)
    Condition <balancing_internal_error>
      Error in `stack_psi_blocks()`:
      ! Every block of the stacked estimating function must cover every observation.
      x Block 2 carries 2 columns for a sample of 4.

# balancing_internal_error: a bare-vector block of the wrong length

    Code
      stack_psi_blocks(list(wide_enough, c(1, 2)), n)
    Condition <balancing_internal_error>
      Error in `stack_psi_blocks()`:
      ! Every block of the stacked estimating function must cover every observation.
      x Block 2 holds 2 values for a sample of 4.

# balancing_internal_error: a bare-vector block the reduction refuses

    Code
      sum_psi_blocks(list(wide_enough, c(1, 2)), n)
    Condition <balancing_internal_error>
      Error in `sum_psi_blocks()`:
      ! Every block of the stacked estimating function must cover every observation.
      x Block 2 holds 2 values for a sample of 4.

