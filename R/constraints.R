# Covariate expansion. build_constraint_matrix() turns a covariate selection and
# a balance_terms() specification into the numeric constraint matrix and a
# serializable recipe; rebuild_constraint_matrix() reconstructs the matrix from
# the recipe alone. Every column is assembled by rebuild_column() from the
# recipe so that a rebuild reproduces the original matrix exactly. Numeric
# columns and their powers cross the boundary standardized to unit scale;
# indicator and quantile columns cross as raw zero/one contrasts, since every
# level is itself a balance target.

# A single per-column recipe record. Fields the balance table and rebuild both
# read are always present; fields that do not apply to a column type are filled
# with NA so that vapply() over the recipe stays type-stable.
new_recipe_record <- function(
  term,
  kind,
  type,
  source,
  power = NA_integer_,
  base_center = NA_real_,
  center = 0,
  scale = 1,
  level = NA_character_,
  partner = NA_character_,
  partner_level = NA_character_,
  probability = NA_real_,
  cutpoint = NA_real_,
  tolerance = 0
) {
  list(
    term = term,
    kind = kind,
    type = type,
    source = source,
    power = power,
    base_center = base_center,
    center = center,
    scale = scale,
    level = level,
    partner = partner,
    partner_level = partner_level,
    probability = probability,
    cutpoint = cutpoint,
    tolerance = tolerance
  )
}

# Raw representation of one base column, used to form interactions and to rebuild
# a column from the data.
base_values <- function(source, level, data) {
  column <- data[[source]]
  if (is.na(level)) {
    as.numeric(column)
  } else {
    as.numeric(as.character(column) == level)
  }
}

# Reconstruct one column from its recipe record and the data.
rebuild_column <- function(record, data) {
  switch(
    record$type,
    numeric = {
      x <- as.numeric(data[[record$source]])
      ((x - record$base_center)^record$power - record$center) / record$scale
    },
    indicator = base_values(record$source, record$level, data),
    interaction = {
      left <- base_values(record$source, record$level, data)
      right <- base_values(record$partner, record$partner_level, data)
      (left * right - record$center) / record$scale
    },
    quantile = as.numeric(as.numeric(data[[record$source]]) <= record$cutpoint)
  )
}

#' Rebuild the constraint matrix from a recipe
#'
#' @param recipe The covariate expansion recipe stored on a [balancing] result.
#' @param .data The data frame the recipe was built from.
#'
#' @return The numeric constraint matrix.
#' @keywords internal
#' @export
rebuild_constraint_matrix <- function(recipe, .data) {
  if (length(recipe) == 0) {
    return(matrix(numeric(0), nrow = nrow(.data), ncol = 0))
  }
  columns <- lapply(recipe, rebuild_column, data = .data)
  matrix <- do.call(cbind, columns)
  colnames(matrix) <- vapply(recipe, function(record) record$term, character(1))
  matrix
}

# Check the names of a per-covariate specification before it is expanded, the
# rules `moments` and `tolerance` share. The callers reach this only for a vector
# that carries names, since an unnamed one is a scalar applied to every covariate.
# A vector that carries names for some elements and not others leaves the unnamed
# ones with no covariate to apply to; the empty name used to be reported as a
# covariate that does not exist, which named the wrong defect. A name given twice
# silently kept the first value and dropped the rest, so the second request never
# reached the fit. Both join the unknown-name check here so one set of rules covers
# both arguments.
#
# A missing name counts as no name. `nzchar()` reads a missing string as a name of
# some length, so a missing name would walk past the check below and be reported
# as a covariate the data does not have, which names the wrong defect again:
# nothing was misspelled, an element was left unnamed.
check_covariate_names <- function(
  values,
  covariates,
  arg,
  call = rlang::caller_env()
) {
  element_names <- names(values)
  unnamed <- sum(is.na(element_names) | !nzchar(element_names))
  if (unnamed > 0) {
    abort(
      c(
        "Every element of {.arg {arg}} must be named when any element is.",
        x = "{unnamed} element{?s} carr{?ies/y} no name.",
        i = "Name each element with one of {.val {covariates}}, or supply a single unnamed value for every covariate."
      ),
      error_class = "balancing_constraints_error",
      call = call
    )
  }
  repeated <- unique(element_names[duplicated(element_names)])
  if (length(repeated) > 0) {
    abort(
      c(
        "{.arg {arg}} names must be unique.",
        x = "{.val {repeated}} {?is/are} named more than once.",
        i = "Give each covariate one value."
      ),
      error_class = "balancing_constraints_error",
      call = call
    )
  }
  unknown <- setdiff(element_names, covariates)
  if (length(unknown) > 0) {
    abort(
      c(
        "{.arg {arg}} names must be covariates.",
        x = "Not {cli::qty(unknown)} {?a covariate/covariates}: {.val {unknown}}.",
        i = "Name each element with one of {.val {covariates}}."
      ),
      error_class = "balancing_constraints_error",
      call = call
    )
  }
  invisible(values)
}

# Expand `moments` into a per-covariate named integer vector. A scalar applies to
# every covariate; a named vector overrides the default of one moment. An unnamed
# vector of length other than one is a classed error rather than a silent recycle,
# and `check_covariate_names()` refuses a partly named, repeated, or unknown name,
# matching how `resolve_tolerance()` treats the same shapes.
resolve_moments <- function(moments, covariates, call = rlang::caller_env()) {
  resolved <- stats::setNames(rep(1L, length(covariates)), covariates)
  if (is.null(moments) || length(moments) == 0) {
    return(resolved)
  }
  if (is.null(names(moments))) {
    if (length(moments) != 1) {
      abort(
        c(
          "{.arg moments} must be a single whole number or a named vector.",
          x = "It has length {length(moments)} and no names.",
          i = "Supply one value for every covariate, or name each element with a covariate."
        ),
        error_class = "balancing_constraints_error",
        call = call
      )
    }
    resolved[] <- as.integer(moments[[1]])
    return(resolved)
  }
  check_covariate_names(moments, covariates, "moments", call = call)
  resolved[names(moments)] <- as.integer(moments[names(moments)])
  resolved
}

# Expand `tolerance` into a per-covariate named numeric vector. A scalar applies
# to every covariate; a named vector sets tolerances per covariate, with unnamed
# covariates left at exact balance (0). Derived columns inherit their source
# covariate's tolerance. An unnamed vector of length other than one is a classed
# error rather than a silent misrecycle, and `check_covariate_names()` refuses a
# partly named, repeated, or unknown name.
resolve_tolerance <- function(
  tolerance,
  covariates,
  call = rlang::caller_env()
) {
  resolved <- stats::setNames(rep(0, length(covariates)), covariates)
  if (is.null(tolerance) || length(tolerance) == 0) {
    return(resolved)
  }
  if (is.null(names(tolerance))) {
    if (length(tolerance) != 1) {
      abort(
        c(
          "{.arg tolerance} must be a single number or a named vector.",
          x = "It has length {length(tolerance)} and no names.",
          i = "Supply one value for every covariate, or name each element with a covariate."
        ),
        error_class = "balancing_constraints_error",
        call = call
      )
    }
    resolved[] <- tolerance
    return(resolved)
  }
  check_covariate_names(tolerance, covariates, "tolerance", call = call)
  resolved[names(tolerance)] <- tolerance
  resolved
}

# Is a numeric covariate a zero/one indicator, for which powers above one repeat
# the column?
is_binary_numeric <- function(v) {
  values <- unique(v[!is.na(v)])
  length(values) <= 2 && all(values %in% c(0, 1))
}

# Which measure the numeric columns are standardized under is the caller's to
# choose, and the two method families choose differently on purpose. The
# estimating-equation methods standardize on the unweighted sample: their targets
# are computed from these same columns, so a per-column scaling cancels between a
# constraint and its target and the fit is invariant to it, while the unweighted
# scale conditions the Newton step better. Where the choice does bite, the
# continuous entropy path, that path centers its own targets on the base measure
# rather than asking for a different column scale here. The quadratic-program
# family standardizes under the sampling weights, because its constraint rows are
# pinned at zero and nothing downstream re-centers them: such a row says the
# weighted mean matches the sample only if the column was centered under the same
# measure the fit reports on, so the sampling weights have to sit in the column
# scale itself.

#' Build the constraint matrix and recipe
#'
#' @param .data The data frame.
#' @param .covariates A character vector of covariate column names.
#' @param constraints A [balance_terms] specification, or `NULL` for the default
#'   first-moment constraints.
#' @param exposure_type One of `"binary"`, `"categorical"`, or `"continuous"`.
#' @param sampling_weights Optional sampling weights. When supplied, numeric
#'   columns are standardized to weighted mean zero and unit weighted standard
#'   deviation rather than the unweighted sample scale.
#' @param call The calling environment, used to build the error's call so a
#'   constraint error names the user-facing function.
#'
#' @return A list with `matrix` and `recipe`.
#' @keywords internal
#' @export
build_constraint_matrix <- function(
  .data,
  .covariates,
  constraints,
  exposure_type,
  sampling_weights = NULL,
  call = rlang::caller_env()
) {
  if (is.null(constraints)) {
    constraints <- balance_terms(moments = 1L)
  }
  quantiles <- constraints@quantiles
  if (!is.null(quantiles) && identical(exposure_type, "continuous")) {
    abort(
      c(
        "Quantile constraints require a discrete exposure.",
        x = "The exposure is {.val continuous}.",
        i = "Drop {.arg quantiles} from {.fn balance_terms}, or balance moments instead."
      ),
      error_class = "balancing_constraints_error",
      call = call
    )
  }

  moments <- resolve_moments(constraints@moments, .covariates, call = call)
  tolerances <- resolve_tolerance(
    constraints@tolerance,
    .covariates,
    call = call
  )

  # Numeric columns cross the boundary standardized to weighted mean zero and unit
  # weighted standard deviation, matching the scale the reference implementations
  # and the core's covariate transform use. Without sampling weights the weighted
  # statistics reduce to the unweighted ones, so an ordinary fit is unaffected.
  center_fn <- if (is.null(sampling_weights)) {
    mean
  } else {
    function(x) weighted_center(x, sampling_weights)
  }
  scale_fn <- if (is.null(sampling_weights)) {
    stats::sd
  } else {
    function(x) weighted_scale(x, sampling_weights)
  }

  records <- list()
  interaction_bases <- list()

  for (cov in .covariates) {
    v <- .data[[cov]]
    if (is.factor(v) || is.character(v)) {
      levels <- if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
      for (level in levels) {
        records[[length(records) + 1]] <- new_recipe_record(
          term = paste0(cov, "_", level),
          kind = "moment",
          type = "indicator",
          source = cov,
          level = level,
          tolerance = tolerances[[cov]]
        )
        interaction_bases[[length(interaction_bases) + 1]] <- list(
          source = cov,
          level = level,
          is_factor = TRUE
        )
      }
    } else if (is.logical(v) || is_binary_numeric(v)) {
      if (moments[[cov]] > 1L) {
        alert_info(
          "Ignoring moments above one for the binary covariate {.val {cov}}."
        )
      }
      records[[length(records) + 1]] <- new_recipe_record(
        term = cov,
        kind = "moment",
        type = "indicator",
        source = cov,
        level = NA_character_,
        tolerance = tolerances[[cov]]
      )
      interaction_bases[[length(interaction_bases) + 1]] <- list(
        source = cov,
        level = NA_character_,
        is_factor = FALSE
      )
    } else {
      base_center <- center_fn(v)
      for (power in seq_len(moments[[cov]])) {
        raw <- (v - base_center)^power
        center <- center_fn(raw)
        scale <- scale_fn(raw)
        if (scale == 0) {
          scale <- 1
        }
        records[[length(records) + 1]] <- new_recipe_record(
          term = if (power == 1L) cov else paste0(cov, "^", power),
          kind = if (power == 1L) "moment" else "power",
          type = "numeric",
          source = cov,
          power = as.integer(power),
          base_center = base_center,
          center = center,
          scale = scale,
          tolerance = tolerances[[cov]]
        )
      }
      interaction_bases[[length(interaction_bases) + 1]] <- list(
        source = cov,
        level = NA_character_,
        is_factor = FALSE
      )
    }
  }

  if (isTRUE(constraints@interactions)) {
    records <- c(
      records,
      interaction_records(
        interaction_bases,
        .data,
        tolerances,
        center_fn,
        scale_fn
      )
    )
  }

  if (!is.null(quantiles)) {
    records <- c(
      records,
      quantile_records(.covariates, .data, quantiles, tolerances)
    )
  }

  matrix <- rebuild_constraint_matrix(records, .data)

  # Constant columns come out first. They are not aliased with anything, so the
  # rank check does not always reach them, and a constant column that survives
  # leaves the achieved balance undefined rather than met.
  constant <- constant_columns(matrix)
  if (length(constant) > 0) {
    constant_terms <- record_terms(records[constant])
    if (length(constant) == ncol(matrix)) {
      abort(
        c(
          "Every constraint column takes a single value.",
          x = "The constraint{?s} {.val {constant_terms}} {?is/are} constant.",
          i = "Supply at least one covariate that varies across the sample."
        ),
        error_class = "balancing_constraints_error",
        call = call
      )
    }
    alert_info(
      "Dropping constant constraint{?s} {.val {constant_terms}}."
    )
    keep <- setdiff(seq_along(records), constant)
    records <- records[keep]
    matrix <- matrix[, keep, drop = FALSE]
  }

  dropped <- aliased_columns(matrix)
  if (length(dropped) > 0) {
    dropped_terms <- record_terms(records[dropped])
    alert_info(
      "Dropping aliased constraint{?s} {.val {dropped_terms}}."
    )
    keep <- setdiff(seq_along(records), dropped)
    records <- records[keep]
    matrix <- matrix[, keep, drop = FALSE]
  }

  list(matrix = matrix, recipe = records)
}

record_terms <- function(records) {
  vapply(records, function(record) record$term, character(1))
}

# Build the interaction records: pairwise products of distinct base columns,
# skipping products of two indicators of the same factor.
interaction_records <- function(bases, data, tolerances, center_fn, scale_fn) {
  records <- list()
  n_bases <- length(bases)
  for (i in seq_len(n_bases)) {
    for (j in seq_len(n_bases)) {
      if (j <= i) {
        next
      }
      left <- bases[[i]]
      right <- bases[[j]]
      same_factor <- left$is_factor &&
        right$is_factor &&
        identical(left$source, right$source)
      if (same_factor) {
        next
      }
      left_values <- base_values(left$source, left$level, data)
      right_values <- base_values(right$source, right$level, data)
      product <- left_values * right_values
      center <- center_fn(product)
      scale <- scale_fn(product)
      if (scale == 0) {
        scale <- 1
      }
      records[[length(records) + 1]] <- new_recipe_record(
        term = interaction_term(left, right),
        kind = "interaction",
        type = "interaction",
        source = left$source,
        level = left$level,
        partner = right$source,
        partner_level = right$level,
        center = center,
        scale = scale,
        tolerance = tolerances[[left$source]]
      )
    }
  }
  records
}

interaction_term <- function(left, right) {
  left_label <- if (is.na(left$level)) {
    left$source
  } else {
    paste0(left$source, "_", left$level)
  }
  right_label <- if (is.na(right$level)) {
    right$source
  } else {
    paste0(right$source, "_", right$level)
  }
  paste0(left_label, ":", right_label)
}

# Build the quantile indicator records for the numeric covariates. `quantiles`
# is either one probability vector applied to every continuous covariate or a
# named list giving the probabilities per covariate; a covariate absent from the
# list contributes no quantile columns.
quantile_records <- function(covariates, data, quantiles, tolerances) {
  records <- list()
  for (cov in covariates) {
    v <- data[[cov]]
    if (!is.numeric(v) || is_binary_numeric(v)) {
      next
    }
    probabilities <- if (is.list(quantiles)) quantiles[[cov]] else quantiles
    if (is.null(probabilities)) {
      next
    }
    for (probability in probabilities) {
      cutpoint <- stats::quantile(v, probs = probability, names = FALSE)
      records[[length(records) + 1]] <- new_recipe_record(
        term = paste0(cov, "_q", probability),
        kind = "quantile",
        type = "quantile",
        source = cov,
        probability = probability,
        cutpoint = cutpoint,
        tolerance = tolerances[[cov]]
      )
    }
  }
  records
}

# Positions of columns that take a single value across the sample. A constant
# column carries no balance information: every weighting meets it exactly, and
# the statistics reported on it, a standardized mean difference or an
# exposure-covariate correlation, divide by its zero spread. The check reads the
# assembled columns, so one rule covers a constant numeric covariate, the
# indicator of a factor with a single level, and the product of two constant
# bases alike. Every element of such a column is the same arithmetic applied to
# the same inputs, so the extremes agree to the bit and the comparison needs no
# tolerance. Standardization does not always leave a constant column at zero:
# when the center it subtracts differs from the column by a rounding step, the
# division by a rounding-scale spread returns a column of unit magnitude that is
# still constant, which this rule catches and a zero test would not.
constant_columns <- function(columns) {
  which(vapply(
    seq_len(ncol(columns)),
    function(j) {
      values <- columns[, j]
      values <- values[is.finite(values)]
      length(values) == 0L || max(values) == min(values)
    },
    logical(1)
  ))
}

# Positions of columns a rank-revealing QR identifies as aliased. The check runs
# on the assembled constraint columns together with a constant, which is the
# geometry the solver actually sees: the entropy dual normalizes within each
# exposure group, inverse probability tilting and the covariate balancing
# propensity score bind an explicit intercept column, and the stable balancing
# weights, energy, and characteristic function distance quadratic programs each
# carry a group-sum row. A constraint set that is affinely dependent on that
# constant is rank deficient in the solver's geometry even when the columns on
# their own are independent, and a full set of factor level indicators is exactly
# that shape. Counting the constant also makes the drop independent of which
# affine representative the data happen to carry: a covariate that is another
# covariate plus a shift, and a zero/one indicator paired with its complement,
# reduce to the same surviving set as their centered counterparts.
#
# The constant goes in first and is never a candidate for removal, since the
# solver carries it whatever the constraints do. R's default LINPACK `dqrdc2`
# keeps column order and moves only deficient columns to the end, so column `j`
# of the constraints is dropped exactly when its pivot position, `j + 1` in the
# augmented matrix, falls beyond the rank. The later member of an affine set is
# therefore the one dropped, which for a factor is its last level, and the choice
# is deterministic rather than a function of column ordering within the pivot.
#
# Rank is a tolerance question rather than an exact one, so the tolerance is
# named here rather than left to the default: a column whose residual, after the
# constant and the columns ahead of it are projected out, falls below `tol` times
# that column's own norm is moved beyond the rank and dropped as aliased. The
# default 1e-07 is the right order for this matrix because every column arrives
# at unit scale, the numeric and interaction columns standardized and the
# indicator and quantile columns zero or one, so the relative test reads against
# the same magnitude column by column. It sits far above the residual that
# rounding leaves on a set that is dependent in exact arithmetic and far below
# the residual a column that varies on its own keeps.
aliased_columns <- function(columns) {
  if (ncol(columns) == 0L) {
    return(integer(0))
  }
  augmented <- cbind(1, columns)
  decomposition <- qr(augmented, tol = 1e-07)
  if (decomposition$rank == ncol(augmented)) {
    return(integer(0))
  }
  kept <- decomposition$pivot[seq_len(decomposition$rank)]
  sort(setdiff(seq_len(ncol(columns)), kept - 1L))
}
