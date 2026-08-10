# Effect modification on the stacked variance an `ipw()` result assembles. `.by`
# names a modifier, and a result carrying one reports the effects it reports
# without a request, then those same effects within each of the modifier's
# strata, then each non-reference stratum against the reference one.
#
# What lives here is the part of a request that is not the stack: whether a
# request was made at all, the refusals that decide whether one can be answered,
# and the strata the stacked blocks are built from. The blocks themselves live
# in R/ipw-deli.R, appended to the whole-sample blocks whose arithmetic they
# repeat, which is what keeps a grouped system's leading blocks the ungrouped
# system's own rather than a copy of them.

# The group label the whole-sample rows carry. Every other row names a stratum,
# so the rows estimated over everyone need a name of their own for the column to
# name a subgroup in every row of the table rather than in some of them.
ipw_overall_group <- "overall"

# Whether a captured `.by` names nothing, which is the argument's default. A
# result built without a request is the result the method returned before the
# argument existed, down to the columns of its estimates table, so absence is
# asked about here rather than inferred further down from a resolved modifier
# that came back empty.
ipw_by_absent <- function(.by) {
  is.null(.by) || rlang::quo_is_null(.by)
}

# Refuse `.by` for a continuous exposure, the one exposure type whose result
# reports no contrast of standardized means to modify. What the method reports
# there is the marginal structural model's own exposure coefficient, a single
# number for the whole sample with no per-stratum counterpart the stacked blocks
# could build. Fitting the model within each stratum instead reports a
# coefficient apiece and no covariance between them, so the difference between
# two strata could not be tested from those fits; the refusal names that route
# rather than quietly taking it.
check_ipw_by_exposure <- function(.by, call = rlang::caller_env()) {
  if (ipw_by_absent(.by)) {
    return(invisible(NULL))
  }

  abort(
    c(
      "{.fun ipw} does not support {.arg .by} for a continuous exposure.",
      x = "A continuous exposure reports the marginal structural model's own exposure coefficient rather than a contrast of standardized means, so there is no effect within a subgroup to report.",
      i = "Omit {.arg .by} to report the whole-sample effect.",
      i = "Fitting each subgroup on its own subset reports a coefficient per subgroup and no covariance between them, so the difference between two subgroups cannot be tested from those fits."
    ),
    error_class = "balancing_ipw_unsupported_error",
    call = call
  )
}

# The strata a `.by` request names, or `NULL` when it names none. `frame` is the
# frame the modifier is read out of, which is `.data` where the caller supplied
# one and the outcome model's own model frame otherwise, so the modifier is
# selected from the same rows the counterfactual designs are built from.
#
# The checks run in the order a request fails in rather than in the order the
# pieces are computed: which column the request names has to be settled before
# what that column holds, and both before the outcome model is inspected for a
# term letting the effect differ across the strata. That last one is a warning,
# so a request that will not be answered is refused rather than diagnosed first.
ipw_resolve_by <- function(
  .by,
  frame,
  exposure,
  exposure_levels,
  exposure_name,
  outcome_mod,
  call = rlang::caller_env()
) {
  if (ipw_by_absent(.by)) {
    return(NULL)
  }

  name <- ipw_by_column(.by, frame, call = call)
  values <- frame[[name]]

  check_ipw_by_missing(values, name, call = call)
  check_ipw_by_type(values, name, call = call)

  # A level no unit carries names an empty stratum, which has no mean to
  # standardize and no contrast to report. Dropping it is what lets a modifier
  # be subset without being recoded first, and it is the reading both fitted
  # models were built under anyway, since a model frame drops a factor's unused
  # levels on the way in.
  values <- droplevels(as.factor(values))
  levels <- levels(values)
  labels <- paste0(name, " = ", levels)
  indicators <- matrix(
    vapply(
      levels,
      function(level) as.numeric(values == level),
      numeric(length(values))
    ),
    nrow = length(values),
    dimnames = list(NULL, levels)
  )

  check_ipw_by_levels(
    indicators,
    exposure,
    exposure_levels,
    exposure_name,
    labels,
    call = call
  )
  check_ipw_by_interaction(outcome_mod, exposure_name, name, call = call)

  list(
    name = name,
    labels = labels,
    # Each non-reference stratum against the reference one, which is the
    # modifier's first level rather than the first level in sorted order.
    em_labels = if (length(labels) > 1L) {
      paste(labels[-1], "vs", labels[[1]])
    } else {
      character(0)
    },
    indicators = indicators
  )
}

# The one column of `frame` a request names. Selection is tidyselect's, the way
# `balance()` selects its covariates, so a bare name, a string, and a name held
# in a variable all reach the same column and a column that is not there is
# reported by tidyselect in the terms the caller wrote.
#
# Any other number of columns is refused at both ends. Nothing selected is not
# the same request as `.by = NULL`, which is the argument's default and asks for
# no subgroups at all, and two columns describe a crossing whose labels and
# reference stratum the caller has to settle for themselves.
ipw_by_column <- function(.by, frame, call = rlang::caller_env()) {
  selection <- tidyselect::eval_select(
    .by,
    frame,
    allow_rename = FALSE,
    error_call = call
  )

  if (!identical(length(selection), 1L)) {
    abort(
      c(
        "{.arg .by} must name exactly one modifier.",
        x = "It names {length(selection)} column{?s}.",
        i = "The effects are reported within the levels of a single variable. Cross two variables into one column, with {.fun interaction}, and name that column instead."
      ),
      error_class = "balancing_ipw_input_error",
      call = call,
      .envir = environment()
    )
  }

  names(selection)
}

# A missing value names no subgroup, so the units carrying one belong to none of
# the strata the effects would be reported within. Dropping them here would
# report every stratum over a different sample than the whole-sample rows use,
# and the two sets of rows would then describe two populations while sitting in
# one table.
check_ipw_by_missing <- function(values, name, call = rlang::caller_env()) {
  missing <- sum(is.na(values))
  if (missing == 0L) {
    return(invisible(NULL))
  }

  abort(
    c(
      "{.arg .by} must name a modifier with no missing values.",
      x = "{.val {name}} has {missing} missing value{?s}.",
      i = "A missing value names no subgroup, so the units carrying one belong to none of the strata the effects would be reported within.",
      i = "Drop those rows and refit both models, or recode the missing values as a level of their own."
    ),
    error_class = "balancing_ipw_input_error",
    call = call,
    .envir = environment()
  )
}

# The subgroups are the levels of the modifier, so the modifier has to name a
# fixed set of them. A numeric column names one subgroup per distinct value, and
# a logical one names subgroups whose labels would read as the condition that
# produced them rather than as levels a caller declared. Either recoding is the
# caller's to make, since a cut point and a level order are modeling choices
# this method has no way to guess.
check_ipw_by_type <- function(values, name, call = rlang::caller_env()) {
  if (is.factor(values) || is.character(values)) {
    return(invisible(NULL))
  }

  abort(
    c(
      "{.arg .by} must name a factor or a character modifier.",
      x = "{.val {name}} has class {.cls {class(values)}}.",
      i = "The effects are reported within the levels of the modifier, so it has to name a fixed set of subgroups.",
      i = "Cut a continuous column into groups with {.fun cut}, or convert a logical column or a numeric code to a factor, and name that column instead."
    ),
    error_class = "balancing_ipw_input_error",
    call = call,
    .envir = environment()
  )
}

# Require every stratum to hold every exposure level. The effects within a
# stratum are contrasts of the counterfactual means taken over that stratum's
# units, and a stratum in which nobody took some level identifies neither the
# mean there nor any contrast against it. The outcome model still predicts at
# that level, so the fit would return a number: an extrapolation from the strata
# that do hold it, reported as though it had been estimated there.
#
# Membership is what is required, level by level. Counting the distinct values a
# stratum holds is not enough once the exposure has more than two levels, since
# a stratum can hold two of three and still be missing a contrast.
#
# Only the first incomplete stratum is named. Every fix is a coarser modifier,
# after which the check runs again over the new strata, and neither model can be
# refit to supply a comparison the data do not hold.
check_ipw_by_levels <- function(
  indicators,
  exposure,
  exposure_levels,
  exposure_name,
  labels,
  call = rlang::caller_env()
) {
  key <- as.character(exposure)
  absent <- lapply(
    seq_along(labels),
    function(s) {
      setdiff(as.character(exposure_levels), key[indicators[, s] == 1])
    }
  )
  incomplete <- lengths(absent) > 0L
  if (!any(incomplete)) {
    return(invisible(NULL))
  }

  first <- which(incomplete)[[1L]]
  stratum <- labels[[first]]
  missing_levels <- absent[[first]]

  abort(
    c(
      "{.arg .by} must name a modifier whose subgroups each hold every exposure level.",
      x = "{.val {stratum}} holds no unit with {.val {exposure_name}} set to {.val {missing_levels}}.",
      i = "An effect within a subgroup contrasts the exposure levels inside it, so a subgroup missing one of them has no contrast to report there.",
      i = "Use a coarser modifier, one whose subgroups each hold every exposure level. Refitting either model does not help, since the data hold no comparison there."
    ),
    error_class = "balancing_ipw_input_error",
    call = call,
    .envir = environment()
  )
}

# Report an outcome model with no term reading the exposure and the modifier
# together. The subgroup effects are g-computation on the model as it was
# specified, so the effect differs across subgroups only where a term reads both
# columns, and a model carrying no such term is worth saying so about.
#
# What the message says is what was checked, and no more. A model may carry the
# modification through a column derived from the modifier rather than through
# the modifier itself, `y ~ exposure * modifier_hi` reported by `.by = modifier`
# being the case that arises, and its subgroup effects genuinely differ. A
# message announcing that the effect is the same in every subgroup would be
# false there, while the term report stays true. It is a warning rather than a
# refusal either way, since the fit is an honest reading of the model supplied.
check_ipw_by_interaction <- function(
  outcome_mod,
  exposure_name,
  name,
  call = rlang::caller_env()
) {
  reads_both <- vapply(
    model_term_variable_sets(outcome_mod),
    function(variables) {
      exposure_name %in% variables && name %in% variables
    },
    logical(1)
  )
  if (any(reads_both)) {
    return(invisible(NULL))
  }

  term <- paste0(exposure_name, ":", name)
  warn(
    c(
      "{.arg outcome_mod} has no term reading both {.val {exposure_name}} and {.val {name}}.",
      i = "The subgroup effects are g-computation on {.arg outcome_mod} as it was specified, so the effect differs across subgroups only where a term reads the exposure and the modifier together.",
      i = "Add {.code {term}} to {.arg outcome_mod} and refit it to let the effect differ across the levels of {.val {name}}.",
      i = "A model that carries the modification through a column derived from {.val {name}}, such as an indicator, reads that column rather than this one, and its subgroup effects differ even so."
    ),
    warning_class = "balancing_ipw_by_interaction_warning",
    call = call,
    .envir = environment()
  )

  invisible(NULL)
}

# The stacked name of every block a stratum, or a contrast of strata,
# contributes: the name the whole-sample block carries, suffixed with the group
# it belongs to. The means and the contrasts are named by the one rule, so a
# name in the stack and the key the estimates table reads a row back by cannot
# come apart.
ipw_by_names <- function(names, groups) {
  unlist(
    lapply(groups, function(group) paste0(names, "_", group)),
    use.names = FALSE
  )
}
