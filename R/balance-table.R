# Balance-table construction. The achieved balance is recomputed from the
# rebuilt constraint matrix and the extracted weights so that the reported table
# matches what expect_balanced() verifies independently: discrete exposures
# report standardized mean differences against the estimand's target group, and
# continuous exposures report weighted exposure-covariate correlations.

# The margin a weighted statistic may exceed its tolerance by before the fit
# reports it as out of balance. The absolute floor absorbs floating-point noise
# in the exactly balanced case; the small relative term absorbs the structural
# gap between the weighted association the continuous solver constrains and the
# weighted Pearson correlation reported here, so a fit that met its box does not
# warn against its own tolerance.
balance_margin <- function(tolerance) {
  1e-6 + 0.02 * tolerance
}

# Build a bare tibble without a tibble dependency, matching how positively and
# the tidyverse store display tables on result objects.
new_balancing_tibble <- function(cols) {
  vctrs::new_data_frame(
    cols,
    class = c("tbl_df", "tbl", "data.frame")
  )
}

# Standardize a matrix to unweighted mean 0 and unit standard deviation, matching
# the re-standardization expect_balanced() applies before measuring balance.
standardize_columns <- function(m) {
  centers <- colMeans(m)
  scales <- apply(m, 2, stats::sd)
  scales[scales == 0] <- 1
  sweep(sweep(m, 2, centers, "-"), 2, scales, "/")
}

# Weighted mean of every column of `z` over the rows in `idx`.
weighted_column_means <- function(z, idx, w) {
  apply(z[idx, , drop = FALSE], 2, stats::weighted.mean, w = w[idx])
}

compute_balance_table <- function(
  recipe,
  data,
  exposure_vec,
  exposure_type,
  estimand,
  focal_level,
  groups,
  weights,
  tolerance,
  reference = NULL
) {
  if (is.null(reference)) {
    reference <- rep(1, length(weights))
  }
  matrix <- rebuild_constraint_matrix(recipe, data)
  z <- standardize_columns(matrix)
  p <- ncol(z)
  terms <- vapply(recipe, function(term) term$term, character(1))
  kinds <- vapply(recipe, function(term) term$kind, character(1))
  tolerances <- vapply(
    recipe,
    function(term) term$tolerance %||% tolerance,
    numeric(1)
  )

  if (identical(exposure_type, "continuous")) {
    unweighted <- vapply(
      seq_len(p),
      function(j) {
        abs(stats::cor(exposure_vec, z[, j]))
      },
      numeric(1)
    )
    weighted <- vapply(
      seq_len(p),
      function(j) {
        abs(stats::cov.wt(
          cbind(exposure_vec, z[, j]),
          wt = weights,
          cor = TRUE
        )$cor[1, 2])
      },
      numeric(1)
    )
    group_label <- rep("overall", p)
    statistic <- rep("correlation", p)
  } else {
    unit_weights <- rep(1, length(weights))
    if (identical(estimand, "ate")) {
      # For the average treatment effect every group is reweighted toward the
      # same reference distribution, so the imbalance is how far each group's
      # weighted mean sits from that shared reference, not from zero. The
      # reference is the base measure the solver targets, which reduces to the
      # unweighted pooled mean when there are no base or sampling weights.
      level_names <- names(groups)
      pooled_unweighted <- colMeans(z)
      pooled_weighted <- apply(z, 2, stats::weighted.mean, w = reference)
      unweighted <- imbalance_over_groups(
        z,
        groups,
        unit_weights,
        pooled_unweighted
      )
      weighted_by_group <- lapply(level_names, function(level) {
        abs(
          weighted_column_means(z, groups[[level]], weights) - pooled_weighted
        )
      })
      weighted <- do.call(pmax, weighted_by_group)
      argmax <- max.col(
        t(do.call(rbind, weighted_by_group)),
        ties.method = "first"
      )
      group_label <- level_names[argmax]
    } else {
      focal_idx <- groups[[focal_level]]
      target_u <- colMeans(z[focal_idx, , drop = FALSE])
      target_w <- weighted_column_means(z, focal_idx, weights)
      others <- setdiff(names(groups), focal_level)
      unweighted_by_group <- lapply(others, function(level) {
        abs(colMeans(z[groups[[level]], , drop = FALSE]) - target_u)
      })
      weighted_by_group <- lapply(others, function(level) {
        abs(weighted_column_means(z, groups[[level]], weights) - target_w)
      })
      unweighted <- do.call(pmax, unweighted_by_group)
      weighted <- do.call(pmax, weighted_by_group)
      argmax <- max.col(
        t(do.call(rbind, weighted_by_group)),
        ties.method = "first"
      )
      group_label <- others[argmax]
    }
    statistic <- rep("smd", p)
  }

  new_balancing_tibble(list(
    term = terms,
    kind = kinds,
    statistic = statistic,
    group = group_label,
    unweighted = unweighted,
    weighted = weighted,
    tolerance = tolerances,
    within_tolerance = weighted <= tolerances + balance_margin(tolerances)
  ))
}

# Largest absolute weighted column mean across the exposure groups, used for the
# unweighted (unit-weight) side of an ate table.
imbalance_over_groups <- function(z, groups, w, target) {
  per_group <- lapply(names(groups), function(level) {
    means <- weighted_column_means(z, groups[[level]], w)
    if (is.null(target)) abs(means) else abs(means - target)
  })
  do.call(pmax, per_group)
}
