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

# Expand `moments` into a per-covariate named integer vector. A scalar applies to
# every covariate; a named vector overrides the default of one moment.
resolve_moments <- function(moments, covariates) {
  resolved <- stats::setNames(rep(1L, length(covariates)), covariates)
  if (is.null(moments)) {
    return(resolved)
  }
  if (is.null(names(moments))) {
    resolved[] <- as.integer(moments[[1]])
    return(resolved)
  }
  named <- intersect(names(moments), covariates)
  resolved[named] <- as.integer(moments[named])
  resolved
}

# Expand `tolerance` into a per-covariate named numeric vector. A scalar applies
# to every covariate; a named vector sets tolerances per covariate, with unnamed
# covariates left at exact balance (0). Derived columns inherit their source
# covariate's tolerance. An unnamed vector of length other than one, or a name
# that is not a covariate, is a classed error rather than a silent misrecycle.
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
  unknown <- setdiff(names(tolerance), covariates)
  if (length(unknown) > 0) {
    abort(
      c(
        "{.arg tolerance} names must be covariates.",
        x = "Not {cli::qty(unknown)} {?a covariate/covariates}: {.val {unknown}}.",
        i = "Name each element with one of {.val {covariates}}."
      ),
      error_class = "balancing_constraints_error",
      call = call
    )
  }
  resolved[names(tolerance)] <- tolerance
  resolved
}

# Is a numeric covariate a zero/one indicator, for which powers above one repeat
# the column?
is_binary_numeric <- function(v) {
  values <- unique(v[!is.na(v)])
  length(values) <= 2 && all(values %in% c(0, 1))
}

#' Build the constraint matrix and recipe
#'
#' @param .data The data frame.
#' @param .covariates A character vector of covariate column names.
#' @param constraints A [balance_terms] specification, or `NULL` for the default
#'   first-moment constraints.
#' @param exposure_type One of `"binary"`, `"categorical"`, or `"continuous"`.
#' @param sampling_weights Optional sampling weights, reserved for weighted
#'   standardization.
#'
#' @return A list with `matrix` and `recipe`.
#' @keywords internal
#' @export
build_constraint_matrix <- function(
  .data,
  .covariates,
  constraints,
  exposure_type,
  sampling_weights = NULL
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
      error_class = "balancing_constraints_error"
    )
  }

  moments <- resolve_moments(constraints@moments, .covariates)
  tolerances <- resolve_tolerance(constraints@tolerance, .covariates)

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
      base_center <- mean(v)
      for (power in seq_len(moments[[cov]])) {
        raw <- (v - base_center)^power
        center <- mean(raw)
        scale <- stats::sd(raw)
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
      interaction_records(interaction_bases, .data, tolerances)
    )
  }

  if (!is.null(quantiles)) {
    records <- c(
      records,
      quantile_records(.covariates, .data, quantiles, tolerances)
    )
  }

  matrix <- rebuild_constraint_matrix(records, .data)
  dropped <- aliased_columns(records, .data)
  if (length(dropped) > 0) {
    dropped_terms <- vapply(
      records[dropped],
      function(record) record$term,
      character(1)
    )
    alert_info(
      "Dropping aliased constraint{?s} {.val {dropped_terms}}."
    )
    keep <- setdiff(seq_along(records), dropped)
    records <- records[keep]
    matrix <- matrix[, keep, drop = FALSE]
  }

  list(matrix = matrix, recipe = records)
}

# Build the interaction records: pairwise products of distinct base columns,
# skipping products of two indicators of the same factor.
interaction_records <- function(bases, data, tolerances) {
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
      center <- mean(product)
      scale <- stats::sd(product)
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

# The raw, uncentered value of one column, used only for the rank check. Genuine
# duplication (an identical covariate) is detected here, while covariates related
# by an affine shift, which become collinear only after centering, stay distinct.
raw_check_column <- function(record, data) {
  switch(
    record$type,
    numeric = as.numeric(data[[record$source]])^record$power,
    indicator = base_values(record$source, record$level, data),
    interaction = base_values(record$source, record$level, data) *
      base_values(record$partner, record$partner_level, data),
    quantile = as.numeric(as.numeric(data[[record$source]]) <= record$cutpoint)
  )
}

# Positions of columns a rank-revealing QR identifies as aliased.
aliased_columns <- function(records, data) {
  if (length(records) <= 1) {
    return(integer(0))
  }
  raw <- do.call(cbind, lapply(records, raw_check_column, data = data))
  decomposition <- qr(raw)
  if (decomposition$rank == ncol(raw)) {
    return(integer(0))
  }
  kept <- decomposition$pivot[seq_len(decomposition$rank)]
  sort(setdiff(seq_len(ncol(raw)), kept))
}
