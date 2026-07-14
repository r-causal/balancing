# autoplot() draws the three diagnostic views a fitted result supports: a love
# plot of unweighted versus weighted balance, a weight distribution by group,
# and, for the quadratic-program family, a dual-variable bar chart. It is
# registered on ggplot2's generic, so these specs guard behind ggplot2 and
# assert on the object's layers, aesthetics, and data rather than a rendered
# image. plot() prints the same object. The estimating-equation family carries
# no dual variables, so the "duals" view raises a classed error there.

geom_classes <- function(plot) {
  vapply(plot$layers, function(layer) class(layer$geom)[[1]], character(1))
}

# ---- The balance love plot ------------------------------------------------

test_that("autoplot() draws a love plot with points and a tolerance line", {
  skip_if_not_installed("ggplot2")
  fit <- balance(
    sim_binary(200),
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  plot <- ggplot2::autoplot(fit, type = "balance")

  expect_s3_class(plot, "ggplot")
  geoms <- geom_classes(plot)
  expect_true("GeomPoint" %in% geoms)
  # The dashed tolerance line is a vertical rule at the requested tolerance.
  expect_true("GeomVline" %in% geoms)
})

test_that("autoplot() defaults to the balance view", {
  skip_if_not_installed("ggplot2")
  fit <- balance(
    sim_binary(200),
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  plot <- ggplot2::autoplot(fit)

  expect_s3_class(plot, "ggplot")
  expect_true("GeomVline" %in% geom_classes(plot))
})

test_that("the love plot contrasts unweighted and weighted balance", {
  skip_if_not_installed("ggplot2")
  fit <- balance(
    sim_binary(200),
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  plot <- ggplot2::autoplot(fit, type = "balance")
  built <- ggplot2::ggplot_build(plot)$data[[1]]
  table <- tidy(fit)

  # Two points per constraint term, the unweighted and weighted balance, so the
  # point layer carries twice as many rows as there are terms.
  expect_identical(nrow(built), 2L * nrow(table))

  # The plotted values are the absolute balance statistics keyed by series: the
  # unweighted points carry the unweighted column and the weighted points the
  # weighted column, both in balance-table order. A transposed mapping, a wrong
  # column, or a scaling error would fail here.
  plot_data <- plot$data
  expect_setequal(
    as.character(unique(plot_data$sample)),
    c("Unweighted", "Weighted")
  )
  expect_equal(
    plot_data$value[plot_data$sample == "Unweighted"],
    abs(table$unweighted)
  )
  expect_equal(
    plot_data$value[plot_data$sample == "Weighted"],
    abs(table$weighted)
  )
  # The point layer draws exactly those values on the x axis.
  expect_equal(
    sort(built$x),
    sort(unname(c(abs(table$unweighted), abs(table$weighted))))
  )
})

test_that("the love plot tolerance line sits at the requested tolerance", {
  skip_if_not_installed("ggplot2")
  # A positive tolerance places the dashed rule away from the origin, so a
  # mispositioned line cannot masquerade as the zero-tolerance default.
  fit <- balance(
    sim_binary(200),
    exposure,
    c(x1, x2),
    method = bw_sbw(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  plot <- ggplot2::autoplot(fit, type = "balance")
  vline_index <- which(geom_classes(plot) == "GeomVline")
  expect_length(vline_index, 1L)

  built <- ggplot2::ggplot_build(plot)$data[[vline_index]]
  expect_equal(unique(built$xintercept), max(tidy(fit)$tolerance))
  expect_gt(unique(built$xintercept), 0)
})

# ---- The weight distribution ----------------------------------------------

test_that("autoplot() draws a weight distribution grouped by exposure", {
  skip_if_not_installed("ggplot2")
  data <- sim_binary(200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  plot <- ggplot2::autoplot(fit, type = "weights")

  expect_s3_class(plot, "ggplot")
  expect_gte(length(plot$layers), 1L)

  # Each observation's panel group is its own exposure level and the plotted
  # value is its extracted weight, both in observation order. A mislabeled facet
  # grouping would break this identity.
  plot_data <- plot$data
  expect_identical(as.character(plot_data$group), as.character(data$exposure))
  expect_equal(plot_data$weight, as.numeric(weights(fit)))
})

# ---- The dual-variable bar chart ------------------------------------------

for (spec in list(
  list(
    label = "energy",
    method = quote(bw_energy()),
    constraints = quote(balance_terms())
  ),
  list(
    label = "cfd",
    method = quote(bw_cfd()),
    constraints = quote(balance_terms())
  ),
  list(
    label = "sbw",
    method = quote(bw_sbw()),
    constraints = quote(balance_terms(tolerance = 0.05))
  )
)) {
  local({
    spec <- spec
    test_that(
      paste0("autoplot() draws a dual-variable bar chart for ", spec$label),
      {
        skip_if_not_installed("ggplot2")
        fit <- balance(
          sim_binary(200),
          exposure,
          c(x1, x2),
          method = eval(spec$method),
          estimand = "ate",
          constraints = eval(spec$constraints)
        )
        plot <- ggplot2::autoplot(fit, type = "duals")

        expect_s3_class(plot, "ggplot")
        expect_true(any(geom_classes(plot) %in% c("GeomCol", "GeomBar")))
      }
    )
  })
}

test_that("the dual-variable view raises a classed error without duals", {
  skip_if_not_installed("ggplot2")
  # The estimating-equation family produces no solver dual variables.
  fit <- balance(
    sim_binary(200),
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_null(fit@duals)
  expect_error(
    ggplot2::autoplot(fit, type = "duals"),
    class = "balancing_autoplot_duals_error"
  )
})

test_that("autoplot() rejects an unknown view", {
  skip_if_not_installed("ggplot2")
  fit <- balance(
    sim_binary(200),
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_error(ggplot2::autoplot(fit, type = "nonsense"))
})

# ---- plot() prints the autoplot -------------------------------------------

test_that("plot() draws and returns the autoplot invisibly", {
  skip_if_not_installed("ggplot2")
  fit <- balance(
    sim_binary(200),
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  # Draw to a null device so the test leaves no Rplots.pdf behind.
  withr::local_pdf(NULL)
  drawn <- withVisible(plot(fit, type = "weights"))

  expect_s3_class(drawn$value, "ggplot")
  expect_false(drawn$visible)
})

# ---- Dual-variable table --------------------------------------------------

# The dual view is fed by the labeled dual table the quadratic-program family
# stores. Each row is one structural constraint, labeled by the kind of
# constraint it enforces and carrying its dual value. The estimating-equation
# family carries no dual table.

for (spec in list(
  list(
    label = "energy",
    method = quote(bw_energy()),
    constraints = quote(balance_terms())
  ),
  list(
    label = "cfd",
    method = quote(bw_cfd()),
    constraints = quote(balance_terms())
  ),
  list(
    label = "sbw",
    method = quote(bw_sbw()),
    constraints = quote(balance_terms(tolerance = 0.05))
  )
)) {
  local({
    spec <- spec
    test_that(paste0("the ", spec$label, " fit exposes a labeled dual table"), {
      fit <- balance(
        sim_binary(200),
        exposure,
        c(x1, x2),
        method = eval(spec$method),
        estimand = "ate",
        constraints = eval(spec$constraints)
      )
      duals <- fit@duals

      expect_s3_class(duals, "data.frame")
      expect_identical(names(duals), c("constraint", "kind", "dual"))
      expect_gt(nrow(duals), 0)
      expect_type(duals$dual, "double")
      expect_true(all(duals$kind %in% c("group", "moment")))
    })
  })
}

test_that("the estimating-equation family carries no dual table", {
  fit <- balance(
    sim_binary(200),
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  expect_null(fit@duals)
})
