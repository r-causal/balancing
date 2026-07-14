# Specs for the ipw() integration: the estimating-equations contract that
# variance estimation depends on, and balancing's method on propensity's ipw()
# generic. propensity::ipw() dispatches on the balancing fit and computes a
# stacked-sandwich variance from the fit's estimating-equations container.

# A small, fully seeded binary-exposure data set with a binary and a continuous
# outcome, used for the numerical-oracle tests so the arithmetic is fixed.
ipw_fixture <- function(n = 200) {
  withr::with_seed(101, {
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    z <- stats::rbinom(n, 1L, stats::plogis(0.7 * x1 - 0.5 * x2))
    y_bin <- stats::rbinom(n, 1L, stats::plogis(-0.3 + 0.5 * z + 0.4 * x1))
    y_cont <- 1 + 0.6 * z + 0.5 * x1 - 0.3 * x2 + stats::rnorm(n)
    data.frame(exposure = z, x1 = x1, x2 = x2, y = y_bin, y_cont = y_cont)
  })
}

# A weighted outcome model of the marginal form ipw() expects: the exposure is
# the only predictor, so the fitted marginal means are the weighted group means.
# The weights ride along as a column so the model frame resolves them cleanly.
fit_outcome <- function(formula, data, wts, family) {
  data[[".wts"]] <- wts
  suppressWarnings(
    stats::glm(formula, data = data, family = family, weights = .wts)
  )
}

# The marginal means ipw() reports: predict the outcome model with the exposure
# fixed to each level and average, matching propensity's g-computation.
marginal_means <- function(outcome_mod, data, exposure_name = "exposure") {
  levels <- sort(unique(data[[exposure_name]]))
  d1 <- d0 <- data
  d0[[exposure_name]] <- levels[[1]]
  d1[[exposure_name]] <- levels[[2]]
  list(
    mu0 = mean(stats::predict(outcome_mod, newdata = d0, type = "response")),
    mu1 = mean(stats::predict(outcome_mod, newdata = d1, type = "response"))
  )
}

# A naive risk-difference standard error that treats the weights as fixed: the
# stacked system drops the weight-parameter block, so its only estimating
# functions are the marginal-mean equations w_i z_i (y_i - mu1) and
# w_i (1 - z_i) (y_i - mu0). It is the sandwich a weighted regression reports
# when it ignores that the weights were estimated, used to show the method's
# correction moves the standard error.
naive_rd_se <- function(w, z, y, mu1, mu0) {
  n <- length(z)
  g1 <- w * z * (y - mu1)
  g0 <- w * (1 - z) * (y - mu0)
  bread <- matrix(0, 2L, 2L)
  bread[1L, 1L] <- -sum(w * z) / n
  bread[2L, 2L] <- -sum(w * (1 - z)) / n
  meat <- crossprod(cbind(g1, g0)) / n
  bread_inv <- solve(bread)
  cov <- bread_inv %*% meat %*% t(bread_inv) / n
  contrast <- c(1, -1)
  sqrt(as.numeric(t(contrast) %*% cov %*% contrast))
}

# An independent, scale-coherent risk-difference standard error. The
# implementation composes the sampling weights onto the reported balancing
# weights and rescales the container's weight derivatives to that scale; this
# oracle instead builds the whole stacked M-estimator at the scale the container
# stores natively, so it shares none of the implementation's rescaling algebra.
# `estimating_equations()@weights_raw` is the balancing weight whose derivative
# is `weight_jacobian`, so refitting the outcome model at those weights and
# coupling through `weight_jacobian` directly gives a coherent stack. Coherent
# stacks at any consistent per-group weight scale return the same risk-difference
# variance, so agreement with the implementation is exact up to floating point.
coherent_rd_se <- function(fit, data, sampling = NULL) {
  ee <- estimating_equations(fit)
  n <- nrow(data)
  s <- sampling %||% rep(1, n)
  composed <- s * ee@weights_raw
  data[[".coherent_w"]] <- composed
  outcome_mod <- suppressWarnings(stats::glm(
    y ~ exposure,
    data = data,
    family = stats::binomial(),
    weights = .coherent_w
  ))

  family <- outcome_mod$family
  design <- stats::model.matrix(outcome_mod)
  q <- ncol(design)
  eta <- as.numeric(design %*% stats::coef(outcome_mod))
  mu <- family$linkinv(eta)
  mu_eta <- family$mu.eta(eta)
  variance <- family$variance(mu)
  residual <- (data$y - mu) * mu_eta / variance
  working_weight <- composed * mu_eta^2 / variance
  score <- (composed * residual) * design

  levels <- sort(unique(data$exposure))
  d0 <- d1 <- data
  d0$exposure <- levels[[1]]
  d1$exposure <- levels[[2]]
  terms <- stats::delete.response(stats::terms(outcome_mod))
  x0 <- stats::model.matrix(terms, stats::model.frame(terms, d0))
  x1 <- stats::model.matrix(terms, stats::model.frame(terms, d1))
  eta0 <- as.numeric(x0 %*% stats::coef(outcome_mod))
  eta1 <- as.numeric(x1 %*% stats::coef(outcome_mod))
  h0 <- family$linkinv(eta0)
  h1 <- family$linkinv(eta1)
  mu0 <- mean(h0)
  mu1 <- mean(h1)

  psi <- ee@psi
  p <- ncol(psi)
  stacked <- cbind(psi, score, h0 - mu0, h1 - mu1)
  meat <- crossprod(stacked) / n

  size <- p + q + 2L
  theta <- seq_len(p)
  beta <- p + seq_len(q)
  index0 <- p + q + 1L
  index1 <- p + q + 2L
  bread <- matrix(0, size, size)
  bread[theta, theta] <- ee@jacobian / n
  bread[beta, theta] <- crossprod(design * residual, s * ee@weight_jacobian) / n
  bread[beta, beta] <- -crossprod(design * working_weight, design) / n
  bread[index0, beta] <- colSums(family$mu.eta(eta0) * x0) / n
  bread[index1, beta] <- colSums(family$mu.eta(eta1) * x1) / n
  bread[index0, index0] <- -1
  bread[index1, index1] <- -1

  bread_inv <- solve(bread)
  cov <- bread_inv %*% meat %*% t(bread_inv) / n
  contrast <- numeric(size)
  contrast[index1] <- 1
  contrast[index0] <- -1
  sqrt(as.numeric(t(contrast) %*% cov %*% contrast))
}

# ---- Estimating-equations container contract ------------------------------

# These pin the container ipw() consumes. The dimension and column-sum
# assertions hold as soon as a fit populates the container, so they pass today;
# the finite-difference Jacobian assertion additionally requires the psi
# re-evaluation hook and is the contract variance estimation leans on.

test_that("an entropy binary ate fit exposes a consistent container", {
  data <- sim_binary(200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  ee <- estimating_equations(fit)
  n <- nrow(data)
  p <- ncol(ee@psi)
  expect_identical(nrow(ee@psi), n)
  expect_identical(dim(ee@jacobian), c(p, p))
  expect_identical(dim(ee@weight_jacobian), dim(ee@psi))
  expect_length(ee@parameters, p)
  expect_lt(max(abs(colSums(ee@psi))), 1e-6)
})

test_that("an bw_ipt binary ate fit exposes a consistent container", {
  data <- sim_binary(200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  ee <- estimating_equations(fit)
  p <- ncol(ee@psi)
  expect_identical(nrow(ee@psi), nrow(data))
  expect_identical(dim(ee@jacobian), c(p, p))
  expect_identical(dim(ee@weight_jacobian), dim(ee@psi))
  expect_lt(max(abs(colSums(ee@psi))), 1e-6)
})

test_that("a just-identified bw_cbps binary ate fit exposes a consistent container", {
  data <- sim_binary(200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ate"
  )
  ee <- estimating_equations(fit)
  p <- ncol(ee@psi)
  expect_identical(nrow(ee@psi), nrow(data))
  expect_identical(dim(ee@jacobian), c(p, p))
  expect_identical(dim(ee@weight_jacobian), dim(ee@psi))
  expect_lt(max(abs(colSums(ee@psi))), 1e-6)
})

# The analytic Jacobian is the derivative of the estimating-function column sums
# at the solution. A central finite difference through the psi re-evaluation
# hook must reproduce it. The hook re-evaluates psi at new parameters; without
# it the container cannot be finite-difference checked, so these fail until it
# is populated.
for (spec in list(
  list(label = "entropy", method = quote(bw_entropy())),
  list(label = "ipt", method = quote(bw_ipt())),
  list(label = "cbps just-identified", method = quote(bw_cbps()))
)) {
  local({
    spec <- spec
    test_that(
      paste0(
        "the ",
        spec$label,
        " Jacobian matches a finite difference of psi column sums"
      ),
      {
        data <- sim_binary(200)
        fit <- balance(
          data,
          exposure,
          c(x1, x2),
          method = eval(spec$method),
          estimand = "ate"
        )
        ee <- estimating_equations(fit)
        psi_fn <- ee@psi_fn
        theta <- ee@parameters
        eps <- 1e-6
        finite_diff <- vapply(
          seq_along(theta),
          function(j) {
            up <- theta
            down <- theta
            up[j] <- up[j] + eps
            down[j] <- down[j] - eps
            (colSums(psi_fn(up)) - colSums(psi_fn(down))) / (2 * eps)
          },
          numeric(length(theta))
        )
        expect_equal(finite_diff, ee@jacobian, tolerance = 1e-4)
      }
    )
  })
}

# ---- The ipw() method: effect rows and point estimates --------------------

test_that("ipw() returns the binary-outcome effect rows for an entropy fit", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att",
    focal_level = "1"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  result <- propensity::ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)

  expect_s3_class(result, "ipw")
  expect_identical(estimates$effect, c("rd", "log(rr)", "log(or)"))
})

test_that("ipw() returns a single difference row for a continuous outcome", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "att",
    focal_level = "1"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y_cont ~ exposure, data, w, stats::gaussian())

  result <- propensity::ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)

  expect_identical(estimates$effect, "diff")
})

test_that("an ipw() result prints for a balancing fit", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  result <- propensity::ipw(fit, outcome_mod)

  expect_snapshot(print(result))
})

test_that("ipw() point estimates match the plain weighted-glm computation", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())
  means <- marginal_means(outcome_mod, data)

  result <- propensity::ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)
  rd <- estimates$estimate[estimates$effect == "rd"]
  log_rr <- estimates$estimate[estimates$effect == "log(rr)"]

  expect_equal(rd, means$mu1 - means$mu0, tolerance = 1e-8)
  expect_equal(log_rr, log(means$mu1 / means$mu0), tolerance = 1e-8)
})

# ---- Standard errors: qualitative properties and a numerical oracle -------

test_that("ipw() standard errors are finite and positive", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  result <- propensity::ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)

  expect_true(all(is.finite(estimates$std.err)))
  expect_true(all(estimates$std.err > 0))
})

test_that("ipw() standard errors differ from the naive weights-fixed sandwich", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())
  means <- marginal_means(outcome_mod, data)
  naive_se <- naive_rd_se(
    w,
    data$exposure,
    data$y,
    means$mu1,
    means$mu0
  )

  result <- propensity::ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)
  rd_se <- estimates$std.err[estimates$effect == "rd"]

  # Accounting for weight estimation moves the standard error; it does not equal
  # the sandwich that treats the weights as fixed.
  expect_false(isTRUE(all.equal(rd_se, naive_se)))
})

# The stacked sandwich covers the whole estimating-equation family, not entropy
# alone, and across estimands. Inverse probability tilting and the
# just-identified covariate balancing propensity score store their containers at
# the raw M-estimator solution while the reported weights are renormalized per
# group; the average treatment effect renormalizes each group by a different
# per-group factor, so a coherent stack must carry that factor into the weight
# coupling. Each case compares the method against the independent, scale-coherent
# oracle at a tight tolerance.
for (spec in list(
  list(
    label = "an entropy ate fit",
    method = quote(bw_entropy()),
    estimand = "ate",
    focal = NULL
  ),
  list(
    label = "an bw_ipt ate fit",
    method = quote(bw_ipt()),
    estimand = "ate",
    focal = NULL
  ),
  list(
    label = "an bw_ipt att fit",
    method = quote(bw_ipt()),
    estimand = "att",
    focal = "1"
  ),
  list(
    label = "a just-identified bw_cbps ate fit",
    method = quote(bw_cbps()),
    estimand = "ate",
    focal = NULL
  ),
  list(
    label = "a just-identified bw_cbps att fit",
    method = quote(bw_cbps()),
    estimand = "att",
    focal = "1"
  )
)) {
  local({
    spec <- spec
    test_that(
      paste0(
        "ipw() risk-difference standard error matches the coherent oracle for ",
        spec$label
      ),
      {
        data <- ipw_fixture()
        fit <- balance(
          data,
          exposure,
          c(x1, x2),
          method = eval(spec$method),
          estimand = spec$estimand,
          focal_level = spec$focal
        )
        w <- as.numeric(stats::weights(fit))
        outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())
        oracle_se <- coherent_rd_se(fit, data)

        result <- propensity::ipw(fit, outcome_mod)
        estimates <- as.data.frame(result)
        rd_se <- estimates$std.err[estimates$effect == "rd"]

        expect_identical(estimates$effect, c("rd", "log(rr)", "log(or)"))
        expect_true(all(is.finite(estimates$std.err)))
        expect_true(all(estimates$std.err > 0))
        expect_equal(rd_se, oracle_se, tolerance = 1e-8)
      }
    )
  })
}

test_that("ipw() risk-difference standard error is coherent with sampling weights", {
  data <- ipw_fixture()
  data$sw <- withr::with_seed(7, stats::runif(nrow(data), 0.5, 2))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate",
    sampling_weights = sw
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())
  oracle_se <- coherent_rd_se(fit, data, sampling = data$sw)

  result <- propensity::ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)
  rd_se <- estimates$std.err[estimates$effect == "rd"]

  expect_equal(rd_se, oracle_se, tolerance = 1e-8)
})

test_that("ipw() standard errors track a nonparametric bootstrap", {
  skip_on_cran()
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  result <- propensity::ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)
  rd_se <- estimates$std.err[estimates$effect == "rd"]

  n <- nrow(data)
  boot_rd <- withr::with_seed(2024, {
    vapply(
      seq_len(200),
      function(b) {
        idx <- sample.int(n, n, replace = TRUE)
        resampled <- data[idx, , drop = FALSE]
        # A resampled data set can legitimately fail to converge; that replicate
        # drops out through the error handler, and its convergence warning is
        # suppressed so it does not leak into the suite output.
        out <- tryCatch(
          suppressWarnings({
            boot_fit <- balance(
              resampled,
              exposure,
              c(x1, x2),
              method = bw_entropy(),
              estimand = "ate"
            )
            boot_w <- as.numeric(stats::weights(boot_fit))
            boot_mod <- fit_outcome(
              y ~ exposure,
              resampled,
              boot_w,
              stats::binomial()
            )
            means <- marginal_means(boot_mod, resampled)
            means$mu1 - means$mu0
          }),
          error = function(e) NA_real_
        )
        out
      },
      numeric(1)
    )
  })
  boot_se <- stats::sd(boot_rd, na.rm = TRUE)

  # The bootstrap is noisy at this replicate count, so the agreement is loose.
  expect_equal(rd_se, boot_se, tolerance = 0.15)
})

# ---- Arguments: conf_level and estimand -----------------------------------

test_that("ipw() respects conf_level", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  wide <- as.data.frame(propensity::ipw(fit, outcome_mod, conf_level = 0.95))
  narrow <- as.data.frame(propensity::ipw(fit, outcome_mod, conf_level = 0.80))

  wide_width <- wide$ci.upper - wide$ci.lower
  narrow_width <- narrow$ci.upper - narrow$ci.lower

  expect_true(all(narrow$conf.level == 0.80))
  expect_true(all(narrow_width < wide_width))
})

test_that("ipw() rejects an estimand that contradicts the fit", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  expect_error(
    propensity::ipw(fit, outcome_mod, estimand = "att"),
    class = "balancing_estimand_error"
  )

  cnd <- rlang::catch_cnd(
    propensity::ipw(fit, outcome_mod, estimand = "att"),
    classes = "balancing_estimand_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
})

# ---- Input validation -----------------------------------------------------

test_that("ipw() rejects an outcome model that is not a glm or lm", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )

  cnd <- rlang::catch_cnd(
    propensity::ipw(fit, list(coefficients = 1)),
    classes = "balancing_ipw_input_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
})

test_that("ipw() rejects a covariate-adjusted outcome model", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure + x1, data, w, stats::binomial())

  cnd <- rlang::catch_cnd(
    propensity::ipw(fit, outcome_mod),
    classes = "balancing_ipw_input_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
})

test_that("ipw() rejects a supplied data frame without two exposure levels", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())
  one_level <- data
  one_level$exposure <- 1L

  expect_error(
    propensity::ipw(fit, outcome_mod, .data = one_level),
    class = "balancing_ipw_input_error"
  )
})

# ---- Re-evaluation hook validation ----------------------------------------

# The container's psi re-evaluation hook crosses into Rust; a wrong-length
# parameter vector must surface as an R condition rather than a panic.

test_that("the bw_ipt psi_fn rejects a wrong-length parameter vector", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  psi_fn <- estimating_equations(fit)@psi_fn
  expect_error(psi_fn(c(1, 2)))
})

test_that("the entropy psi_fn rejects a wrong-length parameter vector", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  psi_fn <- estimating_equations(fit)@psi_fn
  expect_error(psi_fn(c(1, 2)))
})

# The categorical covariate balancing propensity score reuses the tilt's
# estimating-function entrypoint, so pin that the container's psi_fn reproduces
# the stored psi at the fitted parameters for a categorical fit and a focal
# estimand, guarding the reuse against drift.
test_that("the categorical cbps psi_fn reproduces the stored psi", {
  data <- sim_categorical(200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ate"
  )
  ee <- estimating_equations(fit)
  expect_equal(ee@psi_fn(ee@parameters), ee@psi, tolerance = 1e-10)
})

test_that("the categorical att cbps psi_fn reproduces the stored psi", {
  data <- sim_categorical(200)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "att",
    focal_level = "b"
  )
  ee <- estimating_equations(fit)
  expect_equal(ee@psi_fn(ee@parameters), ee@psi, tolerance = 1e-10)
})

# ---- Unsupported configurations -------------------------------------------

# Fits whose weights do not solve smooth estimating equations cannot supply the
# stacked sandwich, so ipw() raises the shared unsupported condition pointing to
# the bootstrap workflow. The tolerance-relaxed entropy fit, the over-identified
# bw_cbps fit, and the continuous-exposure bw_cbps fit each lack a container.

test_that("ipw() rejects a tolerance-relaxed entropy fit", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  expect_error(
    propensity::ipw(fit, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  )
})

test_that("ipw() rejects an over-identified bw_cbps fit", {
  data <- ipw_fixture()
  fit <- suppressWarnings(balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(over_identified = TRUE),
    estimand = "ate"
  ))
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  expect_error(
    propensity::ipw(fit, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  )
})

test_that("ipw() rejects a continuous-exposure bw_cbps fit", {
  data <- sim_continuous(200)
  data$y <- stats::rbinom(nrow(data), 1L, 0.5)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  expect_error(
    propensity::ipw(fit, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  )
})

# A container alone is not sufficient: v1 supports binary exposures only. A
# categorical or continuous fit carries estimating equations yet still raises
# the unsupported condition with the bootstrap pointer.

test_that("ipw() rejects a categorical-exposure fit that has a container", {
  data <- sim_categorical(200)
  data$y <- stats::rbinom(nrow(data), 1L, 0.5)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_false(is.null(fit@estimating_equations))
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  expect_error(
    propensity::ipw(fit, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  )

  cnd <- rlang::catch_cnd(
    propensity::ipw(fit, outcome_mod),
    classes = "balancing_ipw_unsupported_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
})

test_that("ipw() rejects a categorical-exposure bw_ipt fit", {
  data <- sim_categorical(200)
  data$y <- stats::rbinom(nrow(data), 1L, 0.5)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  expect_false(is.null(fit@estimating_equations))
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  expect_error(
    propensity::ipw(fit, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  )
})

test_that("ipw() rejects a continuous-exposure fit that has a container", {
  data <- sim_continuous(200)
  data$y <- stats::rbinom(nrow(data), 1L, 0.5)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_false(is.null(fit@estimating_equations))
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::gaussian())

  expect_error(
    propensity::ipw(fit, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  )
})

test_that("the unsupported-weights ipw error carries the bootstrap pointer", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  cnd <- rlang::catch_cnd(
    propensity::ipw(fit, outcome_mod),
    classes = "balancing_ipw_unsupported_error"
  )
  expect_snapshot(
    error = TRUE,
    cnd_class = TRUE,
    stop(cnd)
  )
})
