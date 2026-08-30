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

# expect_stacked_psi_matches_rbind() pins the assembly of the stacked estimating
# function against the `rbind()` of its blocks. The stacked sandwich builds the
# S-by-n matrix block by block through `stack_psi_blocks()`, and the claim that
# helper makes is that filling a preallocated matrix produces exactly what
# stacking the same blocks with `rbind()` produces, to the bit.
#
# The comparison is made on the matrix itself rather than on the variance the
# fit reports, because the variance would not see every difference:
# `stacked_covariance()` overwrites the covariance dimnames with the stacked
# parameter names, so an assembly that dropped or invented row names would still
# return the reference variance. It is also made on the blocks a real fit builds
# rather than on a hand-written list, because the shapes vary with the exposure,
# the outcome model, and whether subgroups or a crossing were declared, and no
# fixed list reaches all of them.
#
# The stand-in delegates to the real helper and returns its value, so the fit
# running under it is the fit the package performs and every evaluation the
# finite difference asks for is compared, not only the one at the root.
# Mismatches are collected and reported once at the end: the finite difference
# calls the closure twice per stacked coordinate, and an expectation inside the
# stand-in would turn a single defect into hundreds of failures. The call count
# is asserted too, so a route that stopped assembling its psi matrix through the
# helper would fail here rather than pass vacuously.
expect_stacked_psi_matches_rbind <- function(expr) {
  assemble <- stack_psi_blocks
  mismatches <- character()
  calls <- 0L

  testthat::local_mocked_bindings(
    stack_psi_blocks = function(blocks, n) {
      calls <<- calls + 1L
      stacked <- assemble(blocks, n)
      expected <- do.call(rbind, blocks)
      if (!identical(stacked, expected)) {
        mismatches <<- c(
          mismatches,
          paste0(
            "call ",
            calls,
            ": assembled ",
            paste(dim(stacked), collapse = " by "),
            ", rbind gives ",
            paste(dim(expected), collapse = " by ")
          )
        )
      }
      stacked
    }
  )

  value <- force(expr)

  testthat::expect_identical(mismatches, character())
  testthat::expect_gt(calls, 0L)

  invisible(value)
}

# expect_finite_column() asserts that a data frame carries the named column and
# that every value in it is finite.
#
# Both halves are load-bearing, and the first is the reason the helper exists.
# `is.finite()` on a column a data frame does not have returns `logical(0)`,
# `all(logical(0))` is TRUE, and so a bare `all(is.finite(df$column))` passes
# without reading anything whenever the column it names is absent. The suite
# asserts finiteness on reported columns dozens of times, and every one of those
# assertions would go quiet under a rename of the reported schema rather than
# reporting it. Testing membership first turns that case into a failure that
# names the missing column.
expect_finite_column <- function(object, column) {
  testthat::expect_true(
    column %in% names(object),
    info = paste0("expected a column named ", column)
  )
  testthat::expect_true(
    all(is.finite(object[[column]])),
    info = paste0("expected every value of ", column, " to be finite")
  )
  invisible(object)
}
