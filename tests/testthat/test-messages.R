# Snapshot the user-facing warnings and informational alerts the entropy slice
# raises. The classed fit-time warnings, the class-downgrade warning, and the
# three covariate-expansion alerts are recorded with their message text and
# condition class so a change in wording or class is caught.

# ---- Fit-time warnings ----------------------------------------------------

test_that("balancing_convergence_warning: the iteration cap is reached", {
  # Three Newton steps drive the balance essentially to zero but do not meet the
  # gradient tolerance, so the fit warns about convergence alone. The solver is
  # pinned so that the failed solve is reported rather than retried with the
  # hybrid, which clears this cap: what the snapshot records is the wording a
  # single solver's failure carries, and the two-solver wording is asserted
  # against directly in test-method-entropy.R.
  withr::local_options(balancing.entropy_solver = "newton")
  data <- sim_binary()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(max_iterations = 3L),
      estimand = "ate"
    )
  )
})

# The two families fail the same criterion for opposite reasons, so they are
# advised differently and the difference is asserted rather than snapshotted: a
# snapshot would record the sentences but not that one leads with the tolerance
# and the other with the cap. cli wraps a bullet at the console width, so each
# message is compared with its line breaks folded into single spaces.
flatten_message <- function(condition) {
  gsub("[[:space:]]+", " ", conditionMessage(condition))
}

test_that("balancing_convergence_warning: a quadratic program leads with the tolerance", {
  # The energy quadratic form is indefinite, and once the alternating-direction
  # iteration passes its residual floor it walks away from the optimum instead of
  # stalling at it. A run that spent its cap therefore did not stop short of the
  # answer, and raising the cap makes the iterate worse rather than better. The
  # advice leads with loosening the tolerance, says the weights are not to be
  # relied on, and mentions the cap last.
  data <- sim_binary()
  condition <- expect_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(max_iterations = 5L),
      estimand = "ate"
    ),
    class = "balancing_convergence_warning"
  )
  message <- flatten_message(condition)
  expect_match(message, "convergence_tolerance", fixed = TRUE)
  expect_match(message, "max_iterations", fixed = TRUE)
  expect_lt(
    regexpr("convergence_tolerance", message, fixed = TRUE),
    regexpr("max_iterations", message, fixed = TRUE)
  )
  expect_match(message, "weights", fixed = TRUE)
  expect_match(message, "rely on|relied on")
})

test_that("balancing_convergence_warning: a quadratic program names a reachable tolerance", {
  # Telling a caller to loosen a tolerance is no help without a value to loosen
  # it to, so a fit that asked for more than the solver can deliver is given one
  # the objective reaches.
  data <- sim_binary()
  condition <- expect_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(convergence_tolerance = 1e-10, max_iterations = 5L),
      estimand = "ate"
    ),
    class = "balancing_convergence_warning"
  )
  expect_match(flatten_message(condition), "1e-0?6")
})

test_that("balancing_convergence_warning: an indefinite objective keeps the residual-floor caveat", {
  # Energy balancing and the characteristic function distance energy kernel are
  # the two indefinite quadratic forms, and the caveat about the residual floor
  # is theirs: past that floor the iteration walks away from the optimum, so a
  # larger cap makes the iterate worse.
  data <- sim_binary()
  for (method in list(
    bw_energy(max_iterations = 5L),
    bw_cfd(kernel = "energy", max_iterations = 5L)
  )) {
    condition <- expect_warning(
      balance(data, exposure, c(x1, x2), method = method, estimand = "ate"),
      class = "balancing_convergence_warning"
    )
    expect_match(flatten_message(condition), "residual floor", fixed = TRUE)
  }
})

test_that("balancing_convergence_warning: a positive-semidefinite objective is not given the residual-floor caveat", {
  # Every kernel but energy assembles a positive-semidefinite quadratic form,
  # whose alternating-direction iteration descends toward the tolerance for as
  # long as the cap allows. A run that spent its cap there really did stop short,
  # so the cap is an ordinary lever rather than a last resort past a floor.
  data <- sim_binary()
  condition <- expect_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cfd(max_iterations = 5L),
      estimand = "ate"
    ),
    class = "balancing_convergence_warning"
  )
  message <- flatten_message(condition)
  expect_match(message, "max_iterations", fixed = TRUE)
  expect_no_match(message, "residual floor", fixed = TRUE)
  expect_no_match(message, "last resort", fixed = TRUE)
})

test_that("balancing_convergence_warning: the weights caveat names a solve that met no tolerance", {
  # The energy fallback re-solves at a tolerance the problem does reach, and the
  # fit still reports itself unconverged because the requested tolerance was not
  # met. The caveat therefore has to be about a solve that met no tolerance at
  # all rather than about one that missed the tolerance asked for, which would
  # contradict the advice above it.
  data <- sim_binary()
  condition <- expect_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(max_iterations = 5L),
      estimand = "ate"
    ),
    class = "balancing_convergence_warning"
  )
  expect_match(
    flatten_message(condition),
    "met no tolerance",
    fixed = TRUE
  )
})

test_that("balancing_convergence_warning: the estimating-equation wording is unchanged", {
  # The entropy and tilting solvers descend monotonically, so a run that spent its
  # cap really did stop short and more iterations really do help. Their advice
  # keeps leading with the cap, which is what separates the two families, and this
  # spec holds it fixed while the quadratic-program wording moves.
  withr::local_options(balancing.entropy_solver = "newton")
  data <- sim_binary()
  condition <- expect_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(max_iterations = 3L),
      estimand = "ate"
    ),
    class = "balancing_convergence_warning"
  )
  expect_match(
    flatten_message(condition),
    "Increase `max_iterations` or loosen `convergence_tolerance`",
    fixed = TRUE
  )
})

test_that("balancing_balance_warning: achieved balance exceeds the tolerance", {
  # A continuous tolerance without the second distribution moment leaves the
  # exposure variance free, so the weighted correlation exceeds the requested
  # bound and the fit warns.
  data <- sim_continuous()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      estimand = "ate",
      constraints = balance_terms(tolerance = 0.1)
    )
  )
})

# The two cases below assert the rendered magnitude directly rather than through
# a snapshot, because the value is the whole point of the sentence and a snapshot
# would record whatever the formatter happens to produce.

test_that("balancing_balance_warning: a residual imbalance keeps its magnitude", {
  # The warning names how far the fit missed, so an imbalance smaller than the
  # display precision still has to read as a number. Printed at a fixed number of
  # decimal places, 4.9e-06 becomes "0.0000", which states that the fit balanced
  # exactly and contradicts the sentence above it.
  condition <- expect_warning(
    warn_balance_exceeded(4.9e-06),
    class = "balancing_balance_warning"
  )
  expect_match(
    conditionMessage(condition),
    "The largest imbalance is 4.9e-06.",
    fixed = TRUE
  )
})

test_that("balancing_balance_warning: an ordinary imbalance stays legible", {
  # The same format has to leave an imbalance at the scale a caller acts on
  # readable, which is three significant digits rather than four decimal places.
  condition <- expect_warning(
    warn_balance_exceeded(0.1751),
    class = "balancing_balance_warning"
  )
  expect_match(
    conditionMessage(condition),
    "The largest imbalance is 0.175.",
    fixed = TRUE
  )
})

test_that("balancing_ignored_argument_warning: two_step without over_identified", {
  # The two-step weighting matrix belongs to the over-identified criterion, so
  # requesting it on a just-identified fit has no effect; the fit warns that the
  # argument is ignored and proceeds.
  data <- sim_binary()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(two_step = FALSE, over_identified = FALSE),
      estimand = "ate"
    )
  )
})

test_that("balancing_ignored_argument_warning: over_identified for a categorical exposure", {
  # The over-identified criterion is defined for a binary exposure alone, so a
  # categorical fit warns that the request is ignored and returns the exactly
  # balancing solution. The message names the exposure type it was raised for.
  data <- sim_categorical()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(over_identified = TRUE),
      estimand = "ate"
    )
  )
})

test_that("balancing_ignored_argument_warning: over_identified for a continuous exposure", {
  data <- sim_continuous()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(over_identified = TRUE),
      estimand = "ate"
    )
  )
})

test_that("balancing_ignored_argument_warning: link for a continuous exposure", {
  # The link names a propensity model, which the continuous form does not fit:
  # it balances the exposure-covariate covariance through an exponential tilt
  # instead. The message names the argument it dropped and why the exposure type
  # has nothing to apply it to.
  data <- sim_continuous()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(link = "probit"),
      estimand = "ate"
    )
  )
})

test_that("balancing_ignored_argument_warning: every argument a categorical fit ignores", {
  # Once the over-identified request is ignored the fit is not over-identified,
  # so the two-step weighting matrix has no criterion to weight either. Each
  # ignored argument carries its own warning rather than the first standing in
  # for both.
  data <- sim_categorical()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cbps(over_identified = TRUE, two_step = FALSE),
      estimand = "ate"
    )
  )
})

test_that("balancing_ignored_argument_warning: a clarabel pin the energy kernel cannot honor", {
  # The energy kernel's quadratic term is indefinite, which the interior-point
  # backend refuses, so a pinned clarabel request cannot be honored for it. The
  # fit names the request it dropped and the backend that ran instead.
  withr::local_options(balancing.qp_backend = "clarabel")
  data <- sim_binary()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_cfd(kernel = "energy"),
      estimand = "ate"
    )
  )
})

test_that("balancing_ignored_argument_warning: a clarabel pin energy balancing cannot honor", {
  # Energy balancing assembles the same indefinite quadratic form as the energy
  # kernel, so the interior-point backend refuses it and a pinned clarabel
  # request cannot be honored. The fit names the request it dropped and the
  # backend that ran instead.
  withr::local_options(balancing.qp_backend = "clarabel")
  data <- sim_binary()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_energy(),
      estimand = "ate"
    )
  )
})

test_that("balancing_ignored_argument_warning: focal_level with a pooled estimand", {
  # The average treatment effect reweights every exposure group rather than
  # holding one fixed, so it has no focal level to resolve and a supplied one is
  # never validated against the data. The fit names the estimand that ignores it.
  data <- sim_binary()
  expect_balancing_warning(
    balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      estimand = "ate",
      focal_level = 1
    )
  )
})

# ---- Class-downgrade warning ----------------------------------------------

test_that("balancing_class_downgrade_warning: mismatched estimands", {
  x <- bw(c(1, 2), estimand = "ate")
  y <- bw(c(3, 4), estimand = "att")
  expect_balancing_warning(vctrs::vec_c(x, y))
})

# ---- Covariate-expansion alerts -------------------------------------------

test_that("alert: the detected exposure type is announced", {
  withr::local_options(balancing.quiet = FALSE)
  data <- sim_binary()
  expect_snapshot(
    invisible(balance(
      data,
      exposure,
      c(x1, x2),
      method = bw_entropy(),
      estimand = "ate"
    ))
  )
})

test_that("alert: the exposure is excluded from a covariate selection", {
  withr::local_options(balancing.quiet = FALSE)
  data <- sim_binary()
  expect_snapshot(
    invisible(balance(
      data,
      exposure,
      everything(),
      method = bw_entropy(),
      estimand = "ate"
    ))
  )
})

test_that("alert: aliased constraint columns are dropped", {
  withr::local_options(balancing.quiet = FALSE)
  data <- sim_binary()
  data$x1_copy <- data$x1
  expect_snapshot(
    invisible(balance(
      data,
      exposure,
      c(x1, x1_copy, x2),
      method = bw_entropy(),
      estimand = "ate",
      exposure_type = "binary"
    ))
  )
})

test_that("alert: moments above one on a binary covariate are ignored", {
  withr::local_options(balancing.quiet = FALSE)
  data <- sim_binary()
  data$flag <- as.integer(data$x1 > 0)
  expect_snapshot(
    invisible(balance(
      data,
      exposure,
      c(x2, flag),
      method = bw_entropy(),
      estimand = "ate",
      exposure_type = "binary",
      constraints = balance_terms(moments = 2L)
    ))
  )
})
