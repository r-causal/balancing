# Reporting for a joint exposure: two treatments intervened on at once.
#
# `causalgenerics::joint_exposure()` crosses two discrete treatments into one
# categorical exposure and records the crossing on the vector. The result is a
# factor, so the categorical path already fits it and already reports what it
# reports for any exposure over those cells: each cell against the reference
# cell. Those rows are arithmetically right and they answer a question nobody
# asked. What the declaration buys is a surface written in the two treatments
# rather than in the cells:
#
#   * the counterfactual mean of each cell as a row of its own, under the effect
#     label `"mean"`, keyed by the cell it belongs to;
#   * the simple effects, each treatment's effect within a fixed level of the
#     other, including the comparisons vs-reference reporting cannot express at
#     all;
#   * the interaction, the difference between two of the first treatment's
#     simple effects, on each collapsible scale.
#
# The cell means are already parameters of the stacked system, as the
# categorical mean block, so the surface adds no means: it reports that block and
# replaces the vs-reference contrast block with contrasts written over the same
# means. That replacement is what makes the cell-against-cell labels unreachable
# rather than merely unreported. A simple effect is a contrast of two mean
# parameters, and an interaction row is the difference of two simple-effect
# parameters, so the double difference is exact by construction rather than by
# two arithmetic routes that have to agree to the last bits.
#
# The odds ratio is off every contrast row here for the reason it is off a
# stratum row: it is noncollapsible, so neither a simple effect reported beside
# one nor a difference of two of them says what it appears to.

# The plan a declared exposure is reported under, or `NULL` when the exposure
# declares no crossing, which is what leaves every undeclared exposure on the
# path it was already on.
#
# `levels` is the fit's own level order, and the plan's cell positions index it,
# since that is the order the stacked mean block is built in. The labels come
# from the declaration instead, which varies its first component fastest and so
# puts the reference cell first. Reading the positions through the labels rather
# than assuming the two orders agree is what keeps a simple effect naming the
# pair of means it actually contrasts.
#
# The plan is positions and labels and nothing numeric, so the seed values, the
# psi rows, and the estimates table read one description of the surface rather
# than three that have to agree.
ipw_joint_plan <- function(exposure, levels, continuous) {
  if (!causalgenerics::is_joint_exposure(exposure)) {
    return(NULL)
  }

  components <- causalgenerics::joint_components(exposure)
  component_names <- names(components)
  component_levels <- unname(components)
  sizes <- lengths(component_levels)
  cells <- levels(exposure)
  cell_at <- function(first, second) {
    match(cells[[first + (second - 1L) * sizes[[1L]]]], levels)
  }

  simple <- list()
  for (component in 1:2) {
    other <- 3L - component
    for (held in seq_len(sizes[[other]])) {
      for (level in seq_len(sizes[[component]])[-1]) {
        pair <- if (component == 1L) {
          c(cell_at(level, held), cell_at(1L, held))
        } else {
          c(cell_at(held, level), cell_at(held, 1L))
        }
        simple[[length(simple) + 1L]] <- list(
          contrast = ipw_joint_contrast_label(
            component_names[[component]],
            component_levels[[component]],
            level
          ),
          group = paste0(
            component_names[[other]],
            " = ",
            component_levels[[other]][[held]]
          ),
          hi = pair[[1L]],
          lo = pair[[2L]]
        )
      }
    }
  }

  # The interaction is reported once, under the first component's framing,
  # rather than twice. It is symmetric in the two treatments, so the difference
  # between the first's simple effects is the difference between the second's,
  # and reporting both would put one quantity in the table under two names.
  #
  # Each row differences two entries of the first component's block, which runs
  # with the held level outer and the compared level inner.
  per_held <- sizes[[1L]] - 1L
  interaction <- list()
  for (held in seq_len(sizes[[2L]])[-1]) {
    for (level in seq_len(sizes[[1L]])[-1]) {
      interaction[[length(interaction) + 1L]] <- list(
        contrast = ipw_joint_contrast_label(
          component_names[[1L]],
          component_levels[[1L]],
          level
        ),
        group = paste0(
          component_names[[2L]],
          " = ",
          component_levels[[2L]][[held]],
          " vs ",
          component_names[[2L]],
          " = ",
          component_levels[[2L]][[1L]]
        ),
        hi = (held - 1L) * per_held + (level - 1L),
        lo = level - 1L
      )
    }
  }

  list(
    cells = cells,
    forms = ipw_contrast_names(continuous, collapsible_only = TRUE),
    simple = simple,
    interaction = interaction
  )
}

# How a treatment's own contrast is written: the treatment, then the level it is
# set to against its reference level. The cells carry `"name = level"` labels
# already, so the colon is what keeps a contrast over one treatment from reading
# as a cell of the crossing.
ipw_joint_contrast_label <- function(name, levels, level) {
  paste0(name, ": ", levels[[level]], " vs ", levels[[1L]])
}

# The stacked name of each joint contrast: the measure, then the row's own
# label. It is the same rule the subgroup blocks are named under, the measure
# suffixed with the identity of the row, so a name in the stack and the key the
# estimates table reads a row back by cannot come apart.
ipw_joint_names <- function(joint) {
  unlist(
    lapply(c(joint$simple, joint$interaction), function(entry) {
      paste0(joint$forms, "_", entry$contrast, " ", entry$group)
    }),
    use.names = FALSE
  )
}

# Each simple effect at the measures its pair of means implies. The pair goes
# through the same contrast arithmetic the vs-reference block uses, handed the
# held mean first and the compared one second, so the measures are defined in
# one place for the whole package rather than restated for this surface.
ipw_joint_simple_values <- function(joint, means, continuous) {
  unlist(
    lapply(joint$simple, function(entry) {
      ipw_contrast_values(
        c(means[[entry$lo]], means[[entry$hi]]),
        continuous,
        collapsible_only = TRUE
      )
    }),
    use.names = FALSE
  )
}

# The two entries of the joint block an interaction row differences, read off
# whichever vector carries the simple effects.
ipw_joint_interaction_values <- function(joint, simple) {
  forms <- length(joint$forms)
  unlist(
    lapply(joint$interaction, function(entry) {
      simple[(entry$hi - 1L) * forms + seq_len(forms)] -
        simple[(entry$lo - 1L) * forms + seq_len(forms)]
    }),
    use.names = FALSE
  )
}

# The values the joint contrast block is seeded at: the simple effects the
# fitted means imply, then the interactions those simple effects imply, which
# is the exact root of each row below.
ipw_joint_values <- function(joint, means, continuous) {
  simple <- ipw_joint_simple_values(joint, means, continuous)
  c(simple, ipw_joint_interaction_values(joint, simple))
}

# The joint contrast rows of one evaluation of the stacked estimating functions.
# They are deterministic functions of the means and of each other, so each row
# is the same value for every unit: nothing at the solution, where that value is
# zero, and everything to the bread, which is what carries their standard errors
# without a delta method.
#
# The simple effects are read off the mean parameters and the interactions off
# the simple-effect parameters rather than off the means again. That is what
# makes an interaction row the difference of two parameters the system already
# carries: its derivative is exact, and the row equals the double difference by
# construction rather than to within the accuracy of a finite difference.
ipw_joint_rows <- function(joint, mean_theta, contrast_theta, continuous, n) {
  simple <- ipw_joint_simple_values(joint, mean_theta, continuous)
  matrix(
    c(simple, ipw_joint_interaction_values(joint, contrast_theta)) -
      contrast_theta,
    nrow = length(contrast_theta),
    ncol = n
  )
}

# The identity columns of the rows a declared crossing reports, and the stacked
# keys each of them is read from: the cell means first, each keyed by its cell
# and belonging to no subgroup, then the simple effects keyed by the treatment
# they contrast and the level the other is held at, then the interaction keyed
# by the same treatment and the two levels of the other being compared.
#
# The mean rows are read out of the stack by name rather than by position, so
# they are reported in the crossing's own cell order whatever order the fit
# solved its levels in.
ipw_joint_identity <- function(joint) {
  entries <- c(joint$simple, joint$interaction)
  forms <- length(joint$forms)

  list(
    keys = c(
      ipw_mean_names(joint$cells, categorical = TRUE),
      ipw_joint_names(joint)
    ),
    effect = c(
      rep("mean", length(joint$cells)),
      rep(joint$forms, times = length(entries))
    ),
    contrast = c(
      joint$cells,
      rep(vapply(entries, function(e) e$contrast, character(1)), each = forms)
    ),
    group = c(
      rep(ipw_overall_group, length(joint$cells)),
      rep(vapply(entries, function(e) e$group, character(1)), each = forms)
    )
  )
}

# The exposure column with its declaration set aside, which is what the
# counterfactual designs are built from. Fixing the exposure to one cell leaves
# the other three unpopulated, and a joint exposure asked to give up cells gives
# up its declaration and says so, so a design built by writing into the declared
# column warns once per cell. The cells and their order are all the designs need,
# and a plain factor over the same levels carries both.
ipw_joint_bare <- function(exposure) {
  factor(as.character(exposure), levels = levels(exposure))
}

# ---- Refusals --------------------------------------------------------------

# Refuse a crossing whose two components carry one name. Every row of the
# surface is keyed by the treatment it contrasts and the level the other
# treatment is held at, and both of those are written from the component names,
# so one name over two treatments writes one key over two different effects.
# The rows are read back out of the stacked system by name, so what a caller
# would get is the first of each colliding pair reported twice in place of the
# two effects the crossing holds.
#
# causalgenerics builds such a crossing: it checks each component in turn, and
# the cells it writes stay distinct, so the pair of names is a defect only once
# something reports in the two treatments rather than in the cells. That is why
# the refusal belongs here, and why it is an input error rather than an
# unsupported one: nothing about the surface is missing, and the crossing means
# what it was declared to mean as soon as each treatment carries its own name.
#
# The names come off the declaration through the same accessor the plan reads
# them with. Recovering them from the labels the plan built would ask which of
# two treatments a colliding label was written from, which is the question the
# collision makes unanswerable.
check_ipw_joint_components <- function(
  joint,
  exposure,
  call = rlang::caller_env()
) {
  if (is.null(joint)) {
    return(invisible(NULL))
  }

  component_names <- names(causalgenerics::joint_components(exposure))
  if (!identical(component_names[[1L]], component_names[[2L]])) {
    return(invisible(NULL))
  }

  abort(
    c(
      "{.fun ipw} cannot report a joint exposure whose two treatments share a name.",
      x = "Both components of the crossing are named {.val {component_names[[1L]]}}.",
      i = "Every row is keyed by the treatment it contrasts and the level the other treatment is held at, so one name over two treatments names two different effects the same way.",
      i = "Declare the crossing with a name of its own for each treatment, as in {.code joint_exposure(a = x, b = y)}."
    ),
    error_class = "balancing_ipw_input_error",
    call = call,
    .envir = environment()
  )
}

# Refuse a focal estimand on a declared crossing. The surface standardizes every
# cell mean to one population, and a focal estimand names one cell as that
# population. The cell means would then be the outcomes among the units treated
# at one corner of the crossing, which is a coherent quantity and not the one
# the simple effects and the interaction rows combine, so reporting it under
# these labels would say something the numbers do not.
check_ipw_joint_estimand <- function(
  joint,
  estimand,
  call = rlang::caller_env()
) {
  if (is.null(joint) || identical(estimand, "ate")) {
    return(invisible(NULL))
  }

  abort(
    c(
      "{.fun ipw} reports a joint exposure for the {.val ate} estimand only.",
      x = "The weights were built for the {.val {estimand}} estimand.",
      i = "Every cell mean on the joint surface standardizes to one population, and a tilted estimand standardizes each of them to a population the simple effects and the interaction are not defined over.",
      i = "Refit the weights for {.val ate}, or drop the declaration with {.code factor(x)} to report each cell against the reference cell under the estimand you have."
    ),
    error_class = "balancing_ipw_unsupported_error",
    call = call,
    .envir = environment()
  )
}

# Refuse `.by` on a declared crossing. Effect modification of a joint
# intervention is a three-way question, the interaction between two treatments
# within the levels of a third variable, and this surface reports neither that
# nor a projection of it that could stand in for it.
check_ipw_joint_by <- function(joint, .by, call = rlang::caller_env()) {
  if (is.null(joint) || ipw_by_absent(.by)) {
    return(invisible(NULL))
  }

  abort(
    c(
      "{.fun ipw} does not support {.arg .by} for a joint exposure.",
      x = "A joint exposure already reports an interaction between two treatments, and reporting it again within the levels of a modifier is a three-way question this surface does not answer.",
      i = "Drop the declaration with {.code factor(x)} to report each cell against the reference cell within the levels of {.arg .by}."
    ),
    error_class = "balancing_ipw_unsupported_error",
    call = call
  )
}
