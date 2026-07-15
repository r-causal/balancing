# expect_balanced() is the primary oracle for the property tests. It rebuilds
# the constraint matrix from the recipe stored on the fitted object, recomputes
# the achieved balance from the extracted weights independently of the object's
# own balance table, and asserts that every constraint term sits within the
# requested tolerance plus a 1e-6 numerical margin.
#
# The result object does not store the raw data by design (only the recipe
# needed to rebuild the matrix), so the source data frame is supplied through
# `.data`. Balance is measured on the standardized scale: for discrete
# exposures each term is a standardized mean difference against the estimand's
# target group, and for continuous exposures each term is a weighted
# exposure-covariate correlation.
#
# Sampling and base weights move the reference the target group defines. The
# balance target is the base-measure-weighted mean, where the base measure is
# the product of the sampling weights and the method's base weights; the
# achieved side uses the extracted weights, which already compose the sampling
# weights. Without either, the base measure is uniform and the targets reduce to
# the unweighted means.
expect_balanced <- function(x, .data, tolerance = 0) {
  constraint_matrix <- rebuild_constraint_matrix(x@recipe, .data)
  w <- as.numeric(stats::weights(x))
  n <- nrow(.data)

  # Standardize on the same scale the constraint columns cross the boundary on:
  # weighted mean zero and unit weighted standard deviation under the sampling
  # weights, reducing to the unweighted scale when none are present.
  z <- standardize_columns(constraint_matrix, x@sampling_weights)

  sampling <- x@sampling_weights %||% rep(1, n)
  base <- tryCatch(x@method@base_weights, error = function(e) NULL) %||%
    rep(1, n)
  measure <- sampling * base

  exposure <- .data[[x@exposure]]

  if (identical(x@exposure_type, "continuous")) {
    achieved <- apply(z, 2, function(column) {
      abs(stats::cov.wt(cbind(exposure, column), wt = w, cor = TRUE)$cor[1, 2])
    })
  } else {
    groups <- as.character(exposure)
    if (identical(x@estimand, "ate")) {
      target <- apply(z, 2, stats::weighted.mean, w = measure)
      achieved <- vapply(
        unique(groups),
        function(g) {
          idx <- groups == g
          weighted <- apply(
            z[idx, , drop = FALSE],
            2,
            stats::weighted.mean,
            w = w[idx]
          )
          max(abs(weighted - target))
        },
        numeric(1)
      )
      achieved <- max(achieved)
    } else {
      focal <- x@focal_level
      focal_idx <- groups == focal
      target <- apply(
        z[focal_idx, , drop = FALSE],
        2,
        stats::weighted.mean,
        w = measure[focal_idx]
      )
      others <- setdiff(unique(groups), focal)
      achieved <- vapply(
        others,
        function(g) {
          idx <- groups == g
          weighted <- apply(
            z[idx, , drop = FALSE],
            2,
            stats::weighted.mean,
            w = w[idx]
          )
          max(abs(weighted - target))
        },
        numeric(1)
      )
    }
  }

  testthat::expect_lte(max(achieved), tolerance + 1e-6)
}
