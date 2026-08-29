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
# A balance value at or below 1e-8 goes further and loses its digits entirely. A
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
  lines <- gsub(
    "(?<![0-9.])[-+]?[0-9]+(?:[.][0-9]+)?e-(?:0*[89]|0*[1-9][0-9]+)(?![0-9])",
    "<1e-8",
    lines,
    perl = TRUE
  )
  round_wide_decimals(lines)
}

# Round a printed decimal to three significant digits once it carries more than
# four, which is the width at which the platform difference starts. Anything
# shorter is left byte-identical: a tolerance such as `0.05`, an imbalance such
# as `0.00515` from formatC(format = "g", digits = 3), a weight range such as
# `0.220`, and any integer all survive untouched. Leading zeros do not count
# toward the width, so `0.005149743` is seven digits wide rather than ten. That
# three-significant-digit rendering also puts an imbalance the fit drove to zero
# in exponent form, where the near-zero rule above replaces it with the `<1e-8`
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
