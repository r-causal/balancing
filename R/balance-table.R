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
#
# The centers and scales are taken as column arithmetic over the whole matrix
# rather than column by column, since the balance table is assembled once per fit
# over as many columns as the fit has constraints. `colSums()` accumulates in the
# same extended precision `sum()` does, so the weighted branch reproduces the
# per-column `weighted_center()` and `weighted_scale()` values (R/utils.R) bit for
# bit rather than merely approximating them: the reliability denominator and the
# zero-scale guard are both carried over unchanged, the former as a single scalar
# because it depends only on the weights.
#
# The unweighted scale takes its sum of squares about a corrected two-pass
# center, the first column mean plus the mean of the residuals from it.
# `stats::sd()` applies the same correction, but carries it entirely in long
# double, whereas this reproduction rounds to double between the two passes. The
# two therefore agree to rounding rather than exactly, occasionally differing by
# a final unit in the last place. The correction still earns its pass: without it
# the plain two-pass form matches `stats::sd()` only while a column's values are
# comparable to their spread, and its error grows with the ratio of the column's
# offset to that spread, whereas the corrected form tracks `stats::sd()` to
# rounding at any offset. The reported centering stays on the plain column mean,
# as before.
standardize_columns <- function(m, sampling_weights = NULL) {
  if (is.null(sampling_weights)) {
    centers <- colMeans(m)
    centered <- sweep(m, 2, centers, "-")
    corrected <- sweep(m, 2, centers + colMeans(centered), "-")
    scales <- sqrt(colSums(corrected^2) / (nrow(m) - 1))
  } else {
    centers <- column_weighted_means(m, sampling_weights)
    centered <- sweep(m, 2, centers, "-")
    scales <- centered_column_scales(centered, sampling_weights)
  }
  constant <- column_is_constant(m)
  centered[, constant] <- 0
  scales[constant | scales == 0] <- 1
  sweep(centered, 2, scales, "/")
}

# Which columns hold one value repeated, read from the values rather than from
# the scale computed above.
#
# A column with no spread is supposed to come out of the centering as zeros and
# be left unscaled, and on the unweighted path with a well-behaved value it
# does. It need not. Both centers round: the weighted one divides a sum of
# products by a sum of weights, the unweighted one a sum by a count, and unless
# the repeated value survives that arithmetic exactly the centered column holds a
# rounding residual instead of a zero. The scale then reports the residual as the
# column's spread and the column standardizes to arbitrary order-one values: a
# constant 0.98 under non-uniform sampling weights comes out at 0.97 in every
# row. The reported balance for such a column is then a report on the rounding.
#
# The reading is exact equality rather than a floor on the computed scale, and
# that is the point. A floor has to be calibrated against a residual whose size
# depends on the column's magnitude, on the sample size, and on whether the
# platform's `long double` is wider than its `double`, and any floor wide enough
# to cover the residual at five thousand rows is wide enough to flatten a column
# offset far from its own spread, which is a column the corrected two-pass center
# above exists to standardize correctly. Equality needs no calibration and
# cannot reach a column that varies at all.
#
# The comparison is against the first row broadcast down the matrix, which reads
# every column in one vectorized pass. The covariate columns carry no missing
# values, which the fit validates before a constraint matrix is built, so the
# missing-value branch below governs only the direct callers.
#
# A column carrying one is read as not constant. The comparison answers a
# missing value with a missing value, and both consumers of this reading index
# with it: `standardize_columns()` writes `centered[, constant] <- 0` and
# `solver_box()` (R/method-entropy.R) writes `column_sd[...] <- 1`. A missing
# subscript makes each of those a silent no-op, so which columns the guards
# reached would depend on a value the guards say nothing about. Reading such a
# column as varying makes that explicit and leaves the standardization behaving
# as it did: the column keeps its computed center and scale, and the missing
# value carries through them into the standardized column.
column_is_constant <- function(m) {
  differences <- colSums(m != rep(m[1L, ], each = nrow(m)))
  !is.na(differences) & differences == 0L
}

# Weighted mean of every column of a matrix, as one pass of column arithmetic.
column_weighted_means <- function(z, w) {
  colSums(z * w) / sum(w)
}

# Sampling-weighted standard deviation of every column of an already-centered
# matrix, with the reliability denominator `weighted_scale()` (R/utils.R) uses.
# The caller centers because both callers hold the centered matrix already: the
# standardization sweeps by the centers it just took, and the solver's tolerance
# box takes them for this alone.
#
# This is the same arithmetic in the same order as the per-column form, so it
# reproduces it bit for bit rather than approximating it. `colSums()` accumulates
# the way `sum()` does, and the products it accumulates are the same products.
# What it saves is the per-column dispatch: on a 20000 by 30 matrix it takes 5.0
# milliseconds where `apply(z, 2, weighted_scale, w = )` takes 7.0.
centered_column_scales <- function(centered, w) {
  total <- sum(w)
  denominator <- total - sum(w * w) / total
  variances <- if (denominator > 0) {
    colSums(centered^2 * w) / denominator
  } else {
    rep(0, ncol(centered))
  }
  sqrt(pmax(variances, 0))
}

# Weighted mean of every column of `z` over the rows in `idx`.
weighted_column_means <- function(z, idx, w) {
  column_weighted_means(z[idx, , drop = FALSE], w[idx])
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
  sampling_weights = NULL,
  matrix = NULL
) {
  if (is.null(reference)) {
    reference <- rep(1, length(weights))
  }
  constraint_target <- match.arg(constraint_target)
  # A caller that already holds the constraint matrix passes it rather than
  # letting the table rebuild one. Rebuilding costs about as much as the fit
  # itself on a wide constraint set, and the recipe reproduces the matrix
  # exactly, so the second build only repeats work. A caller holding nothing but
  # the recipe, such as a diagnostic run against a stored result, leaves this
  # NULL and gets the rebuild.
  if (is.null(matrix)) {
    matrix <- rebuild_constraint_matrix(recipe, data)
  }
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
      pooled_target <- column_weighted_means(z, reference)
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
