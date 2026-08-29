# Snapshot helpers that forward the caller's expression to expect_snapshot() so
# the recorded Code block shows the failing call, not the literal `expr` token.
expect_balancing_error <- function(expr) {
  rlang::inject(testthat::expect_snapshot(
    error = TRUE,
    cnd_class = TRUE,
    !!rlang::enquo(expr)
  ))
}

expect_balancing_warning <- function(expr) {
  rlang::inject(testthat::expect_snapshot(
    cnd_class = TRUE,
    transform = scrub_platform_values,
    !!rlang::enquo(expr)
  ))
}

# A fit's print or summary block. The warning helper above carries the same
# transform because a fit-time warning records the fit it warned about, so its
# snapshot holds a print block too. The error helper does not: an error aborts
# the fit, so nothing is printed.
expect_balancing_snapshot <- function(expr) {
  rlang::inject(testthat::expect_snapshot(
    transform = scrub_platform_values,
    !!rlang::enquo(expr)
  ))
}

# Everything in a printed fit that reports where and how the suite ran rather
# than what the fit did. The iterative solvers take a different floating-point
# path on each platform, so the iteration count, the number of weights resting on
# the minimum-weight floor, and the low-order digits of a balance statistic all
# move between the machine a snapshot was recorded on and the machines the
# package is checked on, while the display around them does not. Scrubbing them
# keeps what these snapshots are for, the shape and wording of the output, and
# gives up the values that cannot be pinned portably.
#
# A balance value below 1e-7 goes further and loses its digits entirely. A
# term the fit drove to zero leaves behind whatever residual its arithmetic
# happened to accumulate, and at that magnitude the residual is a report on the
# platform's floating-point path rather than on the fit: the mantissa and the
# exponent both differ from one machine to the next, so rounding cannot make two
# platforms agree the way it can for a number the fit actually resolved. The one
# thing such a value states, that the term balanced, is what the placeholder
# keeps. The rule is on magnitude rather than on a column, because which columns
# hold a driven-to-zero value depends on the method and the constraint set.
#
# testthat passes a transform to the Output and Condition blocks but not to the
# recorded Code block, so a call in a snapshot is never rewritten. The patterns
# are still written to leave a short number alone, which keeps them safe for the
# tolerances and counts that appear in a condition message.
scrub_platform_values <- function(lines) {
  # cli opens a header with a blank line whose presence depends on what reached
  # the console immediately before it, anywhere in the process. Another test
  # failing earlier in the same run prints a diff and flips it, so the line
  # reports the run rather than the fit and is dropped. Only a leading one: the
  # transform receives a fit's whole block in one piece, so the blank lines
  # between its sections are interior and stay.
  if (length(lines) > 0L && lines[[1L]] == "") {
    lines <- lines[-1L]
  }
  # Both the count and cli's pluralization of "iteration" move with it, so the
  # whole tail goes and the status is what survives.
  lines <- sub(
    "^(\\s*Solver: (?:converged|did not converge)) in [0-9]+ iterations?$",
    "\\1 in <n> iterations",
    lines,
    perl = TRUE
  )
  # The total is the observation count, which is data-determined and stable, so
  # only the count of weights at the floor is replaced.
  lines <- sub(
    "^(\\s*Weights at the minimum-weight floor: )[0-9]+( of [0-9]+)$",
    "\\1<n>\\2",
    lines,
    perl = TRUE
  )
  # Ahead of the rounding: a value like 6.438292e-11 is wide enough to round,
  # and rounding it to 6.44e-11 would leave an exponent that still says nothing
  # portable. Matching the exponent at or below -8 catches both forms, single
  # digit and padded, and the guards on either side keep the pattern off a
  # number that merely ends in something exponent-shaped.
  #
  # The rule reads the exponent rather than the value, so the largest number it
  # can match is just under 1e-7 rather than just under 1e-8, and the
  # placeholder states the cutoff the rule reaches rather than one it does not.
  # Narrowing the rule to the values below 1e-8 alone would mean parsing every
  # match, which buys nothing: a value between the two is as unportable as one
  # below both.
  #
  # A number the text introduces as a tolerance is exempt. The package's own
  # tolerances live in this magnitude range, and a message that reports one is
  # reporting a constant the package chose rather than a residual a platform
  # arrived at, so scrubbing it would hide the value the message is about. The
  # exemption is written as three fixed-width lookbehinds because PCRE takes no
  # variable-width one, and they cover the forms a message uses: the bare word,
  # an equals sign, and "of".
  lines <- gsub(
    paste0(
      "(?<!(?i:tolerance) )(?<!(?i:tolerance) = )(?<!(?i:tolerance) of )",
      "(?<![0-9.])[-+]?[0-9]+(?:[.][0-9]+)?e-(?:0*[89]|0*[1-9][0-9]+)(?![0-9])"
    ),
    "<1e-7",
    lines,
    perl = TRUE
  )
  # The same value at the other end of the formatter. A fit that drove its
  # largest imbalance to an exact zero prints a bare "0" there, which the
  # exponent rule above cannot see, and whether a platform lands on that zero or
  # on a residual of 1e-17 is the floating-point path rather than the fit. The
  # two therefore record the same placeholder. The rewrite is confined to that
  # line, because a zero anywhere else in a printed fit is a count, a tolerance,
  # or a column that is zero by construction.
  lines <- sub(
    "^(\\s*Largest imbalance: )0(?![0-9.])",
    "\\1<1e-7",
    lines,
    perl = TRUE
  )
  collapse_balance_table_spacing(round_wide_decimals(lines))
}

# The columns a printed balance table carries, which is what its header is
# recognized by.
balance_table_column_names <- c(
  "term",
  "kind",
  "statistic",
  "group",
  "unweighted",
  "weighted",
  "tolerance",
  "within_tolerance"
)

# Give up the balance table's alignment and keep its values.
#
# `print.data.frame()` sizes each column to the widest thing in it and
# right-aligns the rest under a header padded to match, so the width of a column
# is a function of the values it holds. A value whose last digit differs by a
# platform therefore moves the header even when the value itself survives the
# rounding above: 0.00511 and 0.0051 print one character apart, and a snapshot
# recorded on one machine fails on another with no visible change in what it
# reports. Collapsing every run of spaces inside the table leaves the values and
# their order, which is what the snapshot is for, and gives up the padding, which
# is not.
#
# The rewrite is confined to the table. Everything else in a printed fit is cli
# output whose spacing is written rather than computed, so it is stable and worth
# pinning as it stands.
#
# The block runs from a header to its last row. A header is a line made only of
# the table's own column names, which covers both the first one and the
# continuation headers `print.data.frame()` emits when the table is too wide for
# the console; a row is a line opening with its row number. The block ends at the
# first line that is neither.
collapse_balance_table_spacing <- function(lines) {
  header <- vapply(
    lines,
    is_balance_table_header,
    logical(1),
    USE.NAMES = FALSE
  )
  row <- grepl("^\\s*[0-9]+\\s", lines, perl = TRUE)

  inside <- FALSE
  for (i in seq_along(lines)) {
    inside <- header[[i]] || (inside && row[[i]])
    if (inside) {
      lines[[i]] <- gsub("\\s{2,}", " ", trimws(lines[[i]]), perl = TRUE)
    }
  }
  lines
}

is_balance_table_header <- function(line) {
  words <- strsplit(trimws(line), "\\s+", perl = TRUE)[[1L]]
  length(words) > 0L &&
    nzchar(words[[1L]]) &&
    all(words %in% balance_table_column_names)
}

# Round a printed decimal to three significant digits once it carries more than
# four, which is the width at which the platform difference starts. Anything
# shorter is left byte-identical: a tolerance such as `0.05`, an imbalance such
# as `0.00515` from formatC(format = "g", digits = 3), a weight range such as
# `0.220`, and any integer all survive untouched. Leading zeros do not count
# toward the width, so `0.005149743` is seven digits wide rather than ten. That
# three-significant-digit rendering also puts an imbalance the fit drove to zero
# in exponent form, where the near-zero rule above replaces it with the `<1e-7`
# placeholder rather than leaving a platform-specific mantissa.
#
# Rounding the parsed value rather than truncating the text is what makes two
# platforms agree: 6.438292e-11 and 6.438290e-11 are different doubles that both
# carry 6.44e-11.
round_wide_decimals <- function(lines) {
  found <- gregexpr("[0-9]*[.][0-9]+(e[-+][0-9]+)?", lines, perl = TRUE)
  regmatches(lines, found) <- lapply(
    regmatches(lines, found),
    function(number) {
      mantissa <- sub("e[-+][0-9]+$", "", number, perl = TRUE)
      digits <- gsub("[.]", "", mantissa, perl = TRUE)
      wide <- nchar(sub("^0+", "", digits, perl = TRUE)) > 4
      number[wide] <- formatC(
        as.numeric(number[wide]),
        format = "g",
        digits = 3
      )
      number
    }
  )
  lines
}
