# glance() reports a one-row fit summary. The imbalance column is named for the
# balance statistic the exposure type uses: a standardized mean difference for
# discrete exposures and a correlation for continuous ones.

test_that("glance() of a binary fit reports the max absolute SMD", {
  data <- sim_binary()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bal_entropy(),
    estimand = "ate"
  )
  glanced <- generics::glance(fit)
  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 1L)
  expect_true("max_absolute_smd" %in% names(glanced))
  expect_false("max_absolute_correlation" %in% names(glanced))
  expect_identical(glanced$method, "Entropy balancing")
  expect_identical(glanced$estimand, "ate")
  expect_identical(glanced$exposure_type, "binary")
  expect_gt(glanced$ess, 0)
})

test_that("glance() of a continuous fit reports the max absolute correlation", {
  data <- sim_continuous()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bal_entropy(),
    estimand = "ate"
  )
  glanced <- generics::glance(fit)
  expect_true("max_absolute_correlation" %in% names(glanced))
  expect_false("max_absolute_smd" %in% names(glanced))
  expect_identical(glanced$exposure_type, "continuous")
})

# The full column contract: a discrete fit and a continuous fit differ only in
# the name of the imbalance column.
test_that("glance() reports exactly the documented columns", {
  discrete <- generics::glance(balance(
    sim_binary(200),
    exposure,
    c(x1, x2),
    method = bal_entropy(),
    estimand = "ate"
  ))
  expect_identical(
    names(discrete),
    c(
      "method",
      "estimand",
      "exposure_type",
      "n",
      "ess",
      "n_constraints",
      "max_absolute_smd",
      "converged",
      "iterations",
      "objective"
    )
  )

  continuous <- generics::glance(balance(
    sim_continuous(200),
    exposure,
    c(x1, x2),
    method = bal_entropy(),
    estimand = "ate"
  ))
  expect_identical(
    names(continuous),
    c(
      "method",
      "estimand",
      "exposure_type",
      "n",
      "ess",
      "n_constraints",
      "max_absolute_correlation",
      "converged",
      "iterations",
      "objective"
    )
  )
})

# glance() summarizes every method family with a single row and a live method
# label, an effective sample size bounded by n, and a finite objective.
glance_specs <- list(
  list(
    label = "entropy",
    method = quote(bal_entropy()),
    constraints = quote(balance_terms())
  ),
  list(
    label = "ipt",
    method = quote(bal_ipt()),
    constraints = quote(balance_terms())
  ),
  list(
    label = "cbps",
    method = quote(bal_cbps()),
    constraints = quote(balance_terms())
  ),
  list(
    label = "energy",
    method = quote(bal_energy()),
    constraints = quote(balance_terms())
  ),
  list(
    label = "cfd",
    method = quote(bal_cfd()),
    constraints = quote(balance_terms())
  ),
  list(
    label = "sbw",
    method = quote(bal_sbw()),
    constraints = quote(balance_terms(tolerance = 0.05))
  )
)

for (spec in glance_specs) {
  local({
    spec <- spec
    test_that(paste0("glance() reports one row for ", spec$label), {
      fit <- balance(
        sim_binary(200),
        exposure,
        c(x1, x2),
        method = eval(spec$method),
        estimand = "ate",
        constraints = eval(spec$constraints)
      )
      glanced <- generics::glance(fit)

      expect_s3_class(glanced, "tbl_df")
      expect_identical(nrow(glanced), 1L)
      expect_identical(glanced$exposure_type, "binary")
      expect_true("max_absolute_smd" %in% names(glanced))
      expect_type(glanced$method, "character")
      expect_gt(glanced$ess, 0)
      expect_lte(glanced$ess, glanced$n + 1e-8)
      expect_identical(glanced$n_constraints, nrow(tidy(fit)))
      expect_true(is.finite(glanced$objective))
    })
  })
}
