# Balance-table construction. The achieved balance is recomputed from the
# rebuilt constraint matrix and the extracted weights. Two geometries are in
# play, and they are reported on separate columns. The headline `unweighted` and
# `weighted` columns carry the arm-to-arm standardized mean difference, the
# cobalt and halfmoon love-plot convention. The `within_tolerance` column and the
# balance warning instead compare against the solver's constraint geometry, the
# arm-to-target residual the tolerance box actually bounds, so a fit that meets
# its documented per-constraint tolerance is never reported as violating it.
# Continuous exposures report weighted exposure-covariate correlations, where the
# two geometries coincide.

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

# Standardize a matrix to weighted mean 0 and unit weighted standard deviation
# under the sampling weights, matching the scale the constraint columns cross the
# boundary on and the denominator the reference implementations report standardized
# mean differences against. Without sampling weights the weighted statistics reduce
# to the unweighted ones. The re-standardization expect_balanced() applies uses the
# same convention.
standardize_columns <- function(m, sampling_weights = NULL) {
  if (is.null(sampling_weights)) {
    centers <- colMeans(m)
    scales <- apply(m, 2, stats::sd)
  } else {
    centers <- apply(m, 2, weighted_center, w = sampling_weights)
    scales <- apply(m, 2, weighted_scale, w = sampling_weights)
  }
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
  reference = NULL,
  constraint_target = c("pooled", "arms"),
  sampling_weights = NULL
) {
  if (is.null(reference)) {
    reference <- rep(1, length(weights))
  }
  constraint_target <- match.arg(constraint_target)
  matrix <- rebuild_constraint_matrix(recipe, data)
  z <- standardize_columns(matrix, sampling_weights)
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
    # For a continuous exposure the reported correlation is the quantity the
    # solver constrains, so the two geometries coincide.
    constraint_residual <- weighted
  } else {
    # Headline statistic: the arm-to-arm standardized mean difference. Each
    # non-target level's weighted mean is compared to a single target level's
    # weighted mean, standardized by the pooled unweighted standard deviation,
    # which `z` already applies, and the largest such contrast is reported. The
    # average treatment effect and the overlap estimand reweight every group, so
    # the target is the first (reference) level: a binary exposure reports the
    # treated-versus-control difference and a multi-level exposure the largest
    # contrast against the reference level. A focal estimand holds the focal group
    # fixed, so its target is the focal level.
    unit_weights <- rep(1, length(weights))
    target_level <- if (estimand %in% c("att", "atu")) {
      focal_level
    } else {
      names(groups)[[1]]
    }
    target_idx <- groups[[target_level]]
    target_u <- weighted_column_means(z, target_idx, unit_weights)
    target_w <- weighted_column_means(z, target_idx, weights)
    other_levels <- setdiff(names(groups), target_level)

    unweighted_by_group <- lapply(other_levels, function(level) {
      abs(weighted_column_means(z, groups[[level]], unit_weights) - target_u)
    })
    weighted_by_group <- lapply(other_levels, function(level) {
      abs(weighted_column_means(z, groups[[level]], weights) - target_w)
    })
    unweighted <- do.call(pmax, unweighted_by_group)
    weighted <- do.call(pmax, weighted_by_group)
    argmax <- max.col(
      t(do.call(rbind, weighted_by_group)),
      ties.method = "first"
    )
    group_label <- other_levels[argmax]
    statistic <- rep("smd", p)

    # Constraint geometry for `within_tolerance` and the balance warning. The
    # residual is each constraint's distance from the target the solver actually
    # enforces, which the arm-to-arm headline may exceed. The geometry depends on
    # the method, not only the estimand: the estimating-equation box family
    # (entropy balancing, inverse probability tilting) holds each arm within
    # tolerance of a shared pooled target for the average treatment effect, so the
    # residual is each arm's distance to it and the headline can reach twice as
    # much when two arms sit at opposite edges of the box. The covariate balancing
    # propensity score instead equates the arms directly, so its residual is the
    # arm-to-arm headline. The focal and overlap estimands always compare against a
    # single held-fixed arm, where the two geometries coincide.
    constraint_residual <- if (
      identical(estimand, "ate") && identical(constraint_target, "pooled")
    ) {
      pooled_target <- apply(z, 2, stats::weighted.mean, w = reference)
      per_arm <- lapply(names(groups), function(level) {
        abs(weighted_column_means(z, groups[[level]], weights) - pooled_target)
      })
      do.call(pmax, per_arm)
    } else {
      weighted
    }
  }

  new_balancing_tibble(list(
    term = terms,
    kind = kinds,
    statistic = statistic,
    group = group_label,
    unweighted = unweighted,
    weighted = weighted,
    tolerance = tolerances,
    within_tolerance = constraint_residual <=
      tolerances + balance_margin(tolerances)
  ))
}
