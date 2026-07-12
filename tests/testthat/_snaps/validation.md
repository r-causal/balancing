# validate_sampling_weights() rejects a non-numeric vector

    Code
      expr
    Condition <balancing_type_error>
      Error:
      ! `sampling_weights` must be numeric, not a <character>.

# validate_sampling_weights() rejects a length mismatch

    Code
      expr
    Condition <balancing_range_error>
      Error:
      ! `sampling_weights` must have one value per observation.
      x It has length 2, but the data have 3 rows.

# validate_sampling_weights() rejects missing values

    Code
      expr
    Condition <balancing_missing_error>
      Error:
      ! `sampling_weights` must not contain missing values.

# validate_sampling_weights() rejects negative values

    Code
      expr
    Condition <balancing_range_error>
      Error:
      ! `sampling_weights` must be non-negative.
      x Found 1 negative value.

