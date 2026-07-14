# tidy() reports the balance table as a one-row-per-constraint tibble. These
# specs pin the column contract across every method and supported exposure type,
# and confirm the covariate expansion is labeled with the documented `kind`
# values. The balance statistic is a standardized mean difference for discrete
# exposures and an exposure-covariate correlation for continuous ones.

tidy_columns <- c(
  "term",
  "kind",
  "statistic",
  "group",
  "unweighted",
  "weighted",
  "tolerance",
  "within_tolerance"
)

# One entry per method and supported exposure type. Stable balancing weights
# tune through a positive tolerance, so its constraints carry one; the rest use
# exact balance. Each expects the balance statistic its exposure type reports.
tidy_specs <- list(
  list(
    label = "entropy binary",
    data = quote(sim_binary(200)),
    method = quote(bw_entropy()),
    estimand = "ate",
    statistic = "smd"
  ),
  list(
    label = "entropy categorical",
    data = quote(sim_categorical(200)),
    method = quote(bw_entropy()),
    estimand = "ate",
    statistic = "smd"
  ),
  list(
    label = "entropy continuous",
    data = quote(sim_continuous(200)),
    method = quote(bw_entropy()),
    estimand = "ate",
    statistic = "correlation"
  ),
  list(
    label = "ipt binary",
    data = quote(sim_binary(200)),
    method = quote(bw_ipt()),
    estimand = "ate",
    statistic = "smd"
  ),
  list(
    label = "ipt categorical",
    data = quote(sim_categorical(200)),
    method = quote(bw_ipt()),
    estimand = "ate",
    statistic = "smd"
  ),
  list(
    label = "cbps binary",
    data = quote(sim_binary(200)),
    method = quote(bw_cbps()),
    estimand = "ate",
    statistic = "smd"
  ),
  list(
    label = "cbps categorical",
    data = quote(sim_categorical(200)),
    method = quote(bw_cbps()),
    estimand = "ate",
    statistic = "smd"
  ),
  list(
    label = "cbps continuous",
    data = quote(sim_continuous(200)),
    method = quote(bw_cbps()),
    estimand = "ate",
    statistic = "correlation"
  ),
  list(
    label = "energy binary",
    data = quote(sim_binary(200)),
    method = quote(bw_energy()),
    estimand = "ate",
    statistic = "smd"
  ),
  list(
    label = "energy categorical",
    data = quote(sim_categorical(200)),
    method = quote(bw_energy()),
    estimand = "ate",
    statistic = "smd"
  ),
  list(
    label = "energy continuous",
    data = quote(sim_continuous(200)),
    method = quote(bw_energy()),
    estimand = "ate",
    statistic = "correlation"
  ),
  list(
    label = "cfd binary",
    data = quote(sim_binary(200)),
    method = quote(bw_cfd()),
    estimand = "ate",
    statistic = "smd"
  ),
  list(
    label = "cfd categorical",
    data = quote(sim_categorical(200)),
    method = quote(bw_cfd()),
    estimand = "ate",
    statistic = "smd"
  ),
  list(
    label = "sbw binary",
    data = quote(sim_binary(200)),
    method = quote(bw_sbw()),
    estimand = "ate",
    statistic = "smd",
    tolerance = 0.05
  ),
  list(
    label = "sbw categorical",
    data = quote(sim_categorical(200)),
    method = quote(bw_sbw()),
    estimand = "ate",
    statistic = "smd",
    tolerance = 0.05
  ),
  list(
    label = "sbw continuous",
    data = quote(sim_continuous(200)),
    method = quote(bw_sbw()),
    estimand = "ate",
    statistic = "correlation",
    tolerance = 0.05
  )
)

for (spec in tidy_specs) {
  local({
    spec <- spec
    test_that(
      paste0("tidy() reports the constraint table for ", spec$label),
      {
        constraints <- if (is.null(spec$tolerance)) {
          balance_terms()
        } else {
          balance_terms(tolerance = spec$tolerance)
        }
        fit <- balance(
          eval(spec$data),
          exposure,
          c(x1, x2),
          method = eval(spec$method),
          estimand = spec$estimand,
          constraints = constraints
        )
        tidied <- tidy(fit)

        expect_s3_class(tidied, "tbl_df")
        expect_identical(names(tidied), tidy_columns)
        # One row per expanded constraint term.
        expect_identical(nrow(tidied), length(fit@recipe))
        expect_true(all(
          tidied$kind %in% c("moment", "power", "interaction", "quantile")
        ))
        expect_true(all(tidied$statistic == spec$statistic))
        expect_type(tidied$unweighted, "double")
        expect_type(tidied$weighted, "double")
        expect_type(tidied$within_tolerance, "logical")
      }
    )
  })
}

# The `kind` labels track the covariate expansion the constraints request. Raw
# centered powers above the first are "power"; a first moment stays "moment";
# pairwise products are "interaction"; a quantile indicator is "quantile".

test_that("tidy() labels higher moments as power terms", {
  fit <- balance(
    sim_binary(200),
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate",
    constraints = balance_terms(moments = 2)
  )
  kinds <- tidy(fit)$kind
  expect_true("moment" %in% kinds)
  expect_true("power" %in% kinds)
})

test_that("tidy() labels pairwise products as interaction terms", {
  fit <- balance(
    sim_binary(200),
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate",
    constraints = balance_terms(interactions = TRUE)
  )
  expect_true("interaction" %in% tidy(fit)$kind)
})

test_that("tidy() labels quantile indicators as quantile terms", {
  fit <- balance(
    sim_binary(200),
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate",
    constraints = balance_terms(quantiles = c(0.25, 0.75))
  )
  expect_true("quantile" %in% tidy(fit)$kind)
})

# A continuous exposure has no exposure groups to contrast, so every row reports
# the single overall correlation.

test_that("tidy() of a continuous fit labels the group as overall", {
  fit <- balance(
    sim_continuous(200),
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  tidied <- tidy(fit)
  expect_true(all(tidied$group == "overall"))
  expect_true(all(tidied$statistic == "correlation"))
})
