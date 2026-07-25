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

# The container also exposes the weights as a function of the parameters.
# `weights_fn(theta)` returns the weights the fit reports, with the sampling
# weights excluded, so a variance estimator can differentiate the weight path
# without the R layer reimplementing any method's math. The contract has three
# parts: at the fitted parameters it reproduces the reported weights; it carries
# the reporting scale the container records on `weights_raw`; and its derivative
# is `weight_jacobian` moved to that same reporting scale. The third part is the
# one downstream inference leans on, because a weight coupling left at the
# container's storage scale is wrong by the per-group reporting factor.

# Assert the whole `weights_fn` contract for a fitted result. The central
# difference uses a step of 1e-6 scaled by each coordinate's magnitude, which
# balances a truncation error of order the step squared against a cancellation
# error of order the double epsilon over the step; both sit near 1e-10 relative
# to the Jacobian's scale. The 1e-6 relative tolerance therefore leaves four
# orders of margin over what a correct implementation reaches, while still
# separating it from a derivative left at the wrong per-group scale, whose error
# is of the order of the scale factor itself.
expect_weights_fn_contract <- function(fit, data) {
  ee <- estimating_equations(fit)
  weights_fn <- ee@weights_fn
  expect_true(rlang::is_function(weights_fn))
  evaluate <- function(theta) as.numeric(weights_fn(theta))

  theta <- ee@parameters
  reported <- as.numeric(
    stats::weights(fit, include_sampling_weights = FALSE)
  )
  expect_equal(evaluate(theta), reported, tolerance = 1e-10)

  # The move from the container's storage scale to the reporting scale is a
  # per-unit ratio, so it is defined only where the stored weight is nonzero.
  # Every method that populates the container stores a nonzero weight for every
  # unit, including the focal units an entropy or tilting fit reports at their
  # renormalized base weights, so the ratio is well defined here rather than
  # merely guarded.
  raw <- ee@weights_raw
  expect_true(all(raw != 0))
  rescaled <- (reported / raw) * ee@weight_jacobian

  finite_diff <- vapply(
    seq_along(theta),
    function(j) {
      step <- 1e-6 * max(1, abs(theta[[j]]))
      up <- theta
      down <- theta
      up[[j]] <- up[[j]] + step
      down[[j]] <- down[[j]] - step
      (evaluate(up) - evaluate(down)) / (2 * step)
    },
    numeric(length(reported))
  )
  expect_equal(finite_diff, rescaled, tolerance = 1e-6)

  # A focal group is reported at weights that do not move with the parameters,
  # so its rows of the difference vanish identically rather than to a tolerance.
  if (!is.null(fit@focal_level)) {
    focal <- which(as.character(data[[fit@exposure]]) == fit@focal_level)
    expect_identical(max(abs(finite_diff[focal, ])), 0)
  }
}

# Base weights away from one for the entropy case that varies them. The focal
# group is reported at its base weights carried to the focal total, so a base
# vector whose own total misses that target separates the reported focal weights
# from the base weights themselves, which a uniform base leaves indistinguishable.
contract_base_weights <- function(n) {
  withr::with_seed(303, stats::runif(n, 0.5, 2))
}

for (spec in list(
  list(
    label = "an entropy ate fit",
    data = quote(sim_binary(200)),
    method = quote(bw_entropy()),
    estimand = "ate",
    focal = NULL
  ),
  list(
    label = "an entropy att fit",
    data = quote(sim_binary(200)),
    method = quote(bw_entropy()),
    estimand = "att",
    focal = "1"
  ),
  list(
    label = "an entropy att fit with base weights",
    data = quote(sim_binary(200)),
    method = quote(bw_entropy(
      base_weights = contract_base_weights(nrow(data))
    )),
    estimand = "att",
    focal = "1"
  ),
  list(
    label = "an bw_ipt ate fit",
    data = quote(sim_binary(200)),
    method = quote(bw_ipt()),
    estimand = "ate",
    focal = NULL
  ),
  list(
    label = "an bw_ipt att fit",
    data = quote(sim_binary(200)),
    method = quote(bw_ipt()),
    estimand = "att",
    focal = "1"
  ),
  list(
    label = "an bw_ipt categorical ate fit",
    data = quote(sim_categorical(200)),
    method = quote(bw_ipt()),
    estimand = "ate",
    focal = NULL
  ),
  list(
    label = "an bw_ipt categorical att fit",
    data = quote(sim_categorical(200)),
    method = quote(bw_ipt()),
    estimand = "att",
    focal = "b"
  ),
  list(
    label = "a just-identified bw_cbps ate fit",
    data = quote(sim_binary(200)),
    method = quote(bw_cbps()),
    estimand = "ate",
    focal = NULL
  ),
  list(
    label = "a just-identified bw_cbps categorical ate fit",
    data = quote(sim_categorical(200)),
    method = quote(bw_cbps()),
    estimand = "ate",
    focal = NULL
  ),
  list(
    label = "a just-identified bw_cbps categorical att fit",
    data = quote(sim_categorical(200)),
    method = quote(bw_cbps()),
    estimand = "att",
    focal = "b"
  )
)) {
  for (sampled in c(FALSE, TRUE)) {
    local({
      spec <- spec
      sampled <- sampled
      test_that(
        paste0(
          "the container of ",
          spec$label,
          " re-evaluates its reported weights",
          if (sampled) " with sampling weights" else ""
        ),
        {
          data <- eval(spec$data)
          sampling <- if (sampled) {
            withr::with_seed(505, stats::runif(nrow(data), 0.5, 2))
          } else {
            NULL
          }
          fit <- balance(
            data,
            exposure,
            c(x1, x2),
            method = eval(spec$method),
            estimand = spec$estimand,
            focal_level = spec$focal,
            sampling_weights = sampling
          )
          expect_weights_fn_contract(fit, data)
        }
      )
    })
  }
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

# ---- Outcome response scale -----------------------------------------------

# A binomial outcome model may be fitted on a numeric 0/1 response or on a
# two-level factor, and both fits model the same 0/1 scale: the family's
# initializer maps the factor to zero for its first level and one otherwise.
# The two fits therefore share their coefficients, so ipw() must return the same
# effects and the same standard errors from either one. Reading the response off
# the model frame and coercing it with as.numeric() would put a factor on the
# 1/2 scale, which leaves the coefficient-only point estimates untouched while
# corrupting every sandwich standard error, so the parity assertion is on
# std.err as much as on estimate.

test_that("ipw() gives identical results for numeric and factor binary outcomes", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  data$y_factor <- factor(
    ifelse(data$y == 1, "yes", "no"),
    levels = c("no", "yes")
  )

  numeric_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())
  factor_mod <- fit_outcome(y_factor ~ exposure, data, w, stats::binomial())

  numeric_result <- as.data.frame(propensity::ipw(fit, numeric_mod))
  factor_result <- as.data.frame(propensity::ipw(fit, factor_mod))

  expect_identical(factor_result$effect, numeric_result$effect)
  expect_equal(factor_result$estimate, numeric_result$estimate)
  expect_equal(factor_result$std.err, numeric_result$std.err)

  # The numeric path is the one pinned against the scale-coherent oracle, so
  # anchor the factor path there too rather than only to its numeric twin.
  factor_rd_se <- factor_result$std.err[factor_result$effect == "rd"]
  expect_equal(factor_rd_se, coherent_rd_se(fit, data), tolerance = 1e-8)
})

# A glm stores its modeled response in `$y` by default, which is the reliable
# source for the 0/1 scale. A fit made with `y = FALSE` discards it, and the
# model frame alone leaves the intended scale of a factor response ambiguous, so
# that combination is refused rather than guessed at. Numeric responses and lm
# fits are unaffected: they keep reading the response from the model frame.

test_that("ipw() rejects a factor outcome model fitted without its response", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  data$y_factor <- factor(
    ifelse(data$y == 1, "yes", "no"),
    levels = c("no", "yes")
  )
  data$.wts <- as.numeric(stats::weights(fit))
  factor_mod <- suppressWarnings(stats::glm(
    y_factor ~ exposure,
    data = data,
    family = stats::binomial(),
    weights = .wts,
    y = FALSE
  ))
  expect_null(factor_mod$y)

  expect_error(
    propensity::ipw(fit, factor_mod),
    class = "balancing_ipw_input_error"
  )

  cnd <- rlang::catch_cnd(
    propensity::ipw(fit, factor_mod),
    classes = "balancing_ipw_input_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
})

# The refusal above is narrow: it covers a factor response only. A numeric
# response carries its modeled scale in the model frame whether or not the fit
# stored one, so dropping the stored response must not change anything. Pinning
# this keeps the abort from being widened to every fit made with `y = FALSE`.

test_that("ipw() gives identical results for a numeric response with y = FALSE", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  data$.wts <- as.numeric(stats::weights(fit))

  stored <- suppressWarnings(stats::glm(
    y ~ exposure,
    data = data,
    family = stats::binomial(),
    weights = .wts
  ))
  dropped <- suppressWarnings(stats::glm(
    y ~ exposure,
    data = data,
    family = stats::binomial(),
    weights = .wts,
    y = FALSE
  ))
  expect_false(is.null(stored$y))
  expect_null(dropped$y)

  stored_result <- as.data.frame(propensity::ipw(fit, stored))
  dropped_result <- as.data.frame(propensity::ipw(fit, dropped))

  expect_identical(dropped_result$effect, stored_result$effect)
  expect_equal(dropped_result$estimate, stored_result$estimate)
  expect_equal(dropped_result$std.err, stored_result$std.err)
})

# An lm stores no response by default, so it always reads one from the model
# frame. That response is already on the modeled scale, so a weighted lm and the
# equivalent weighted gaussian glm, which does carry a stored response, must
# agree. This is the other reader of the model-frame path, and it guards the
# path against being dropped once the stored response covers the glm cases.

test_that("ipw() gives identical results for an lm and a gaussian glm", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  data$.wts <- as.numeric(stats::weights(fit))

  glm_mod <- stats::glm(
    y_cont ~ exposure,
    data = data,
    family = stats::gaussian(),
    weights = .wts
  )
  lm_mod <- stats::lm(y_cont ~ exposure, data = data, weights = .wts)
  expect_false(is.null(glm_mod$y))
  expect_null(lm_mod$y)

  glm_result <- as.data.frame(propensity::ipw(fit, glm_mod))
  lm_result <- as.data.frame(propensity::ipw(fit, lm_mod))

  expect_identical(lm_result$effect, "diff")
  expect_identical(lm_result$effect, glm_result$effect)
  expect_equal(lm_result$estimate, glm_result$estimate)
  expect_equal(lm_result$std.err, glm_result$std.err)
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

# ---- The deli-backed stacked sandwich --------------------------------------

# `ipw_deli_sandwich()` replaces the hand-assembled `stacked_sandwich()` with a
# stack built through deli. It evaluates the whole system as one p-by-n
# estimating-function closure and lets `deli::compute_sandwich()` finite
# difference the bread, so the R layer writes no derivative by hand and the
# effect contrasts become parameters of the stack rather than a delta-method
# gradient applied afterwards.
#
# Its interface is
#
#   ipw_deli_sandwich(
#     container,
#     outcome_mod,
#     frame,
#     exposure_name,
#     sampling_weights = NULL
#   )
#
# where `container` is the fit's `balancing_estimating_equations`, `frame` is
# the data frame holding the exposure, and `sampling_weights` is the fit's
# sampling weight vector or `NULL`. The weight parameter count comes from
# `length(container@parameters)`, and the closure reaches the weight path only
# through `container@psi_fn()` and `container@weights_fn()`, so no method math
# is restated here.
#
# It returns `list(theta =, vcov =)`. `theta` is the stacked parameter vector
# and `vcov` its covariance on the standard-error scale, so `sqrt(diag(vcov))`
# is the vector of standard errors. Both carry the same names, in stacked block
# order:
#
#   theta_w1 ... theta_wp   the weight parameters
#   beta_<column>           one per outcome-model design column
#   mu0, mu1                the marginal means
#   rd, log(rr), log(or)    the contrasts for a non-gaussian outcome model
#   diff                    the single contrast for a gaussian outcome model
#
# The contrast names match the `effect` column `ipw()` reports, so an estimates
# table reads its estimate and standard error straight off `theta` and the
# diagonal of `vcov`.

# Call the new sandwich with the pieces read off a fit, the way the `ipw()`
# method will.
call_deli_sandwich <- function(fit, outcome_mod, data) {
  ipw_deli_sandwich(
    container = estimating_equations(fit),
    outcome_mod = outcome_mod,
    frame = data,
    exposure_name = fit@exposure,
    sampling_weights = fit@sampling_weights
  )
}

# The old path, called at the level `ipw()` calls it: the marginal means, the
# marginal-mean covariance block, and the effect standard errors the
# delta-method gradients produce from that block. Going through the internals
# rather than through `ipw()` keeps the comparison on the variance engine alone.
old_stacked_pieces <- function(fit, outcome_mod, data) {
  exposure_name <- fit@exposure
  levels <- sort(unique(data[[exposure_name]]))
  family <- stats::family(outcome_mod)
  fitted <- outcome_model_pieces(outcome_mod, family)
  design0 <- fixed_exposure_pieces(
    outcome_mod,
    data,
    exposure_name,
    levels[[1]]
  )
  design1 <- fixed_exposure_pieces(
    outcome_mod,
    data,
    exposure_name,
    levels[[2]]
  )
  mu0 <- mean(design0$mu)
  mu1 <- mean(design1$mu)
  covariance <- stacked_sandwich(
    container = estimating_equations(fit),
    weights = as.numeric(stats::weights(fit)),
    outcome = resolve_outcome_response(outcome_mod),
    fitted = fitted,
    design0 = design0,
    design1 = design1,
    mu0 = mu0,
    mu1 = mu1
  )
  estimates <- ipw_estimates(
    mu0 = mu0,
    mu1 = mu1,
    covariance = covariance,
    conf_level = 0.95,
    continuous = is_gaussian_outcome(outcome_mod)
  )
  list(
    covariance = covariance,
    mu0 = mu0,
    mu1 = mu1,
    std_err = stats::setNames(estimates$std.err, estimates$effect)
  )
}

# Old-against-new agreement on the two quantities `ipw()` reports from: the
# marginal-mean covariance block and every effect standard error. The tolerance
# is 1e-6 rather than the 1e-8 the oracle comparisons use because the new bread
# is a central finite difference (`deriv_method = "capprox"`) where the old one
# was written analytically. Measured against these fixtures at deli's default
# step the worst case sits near 2e-7, so 1e-6 leaves margin without letting a
# genuinely different bread through.
expect_deli_parity <- function(fit, outcome_mod, data) {
  old <- old_stacked_pieces(fit, outcome_mod, data)
  result <- call_deli_sandwich(fit, outcome_mod, data)

  block <- result$vcov[c("mu0", "mu1"), c("mu0", "mu1")]
  expect_equal(unname(block), unname(old$covariance), tolerance = 1e-6)

  std_err <- sqrt(diag(result$vcov))[names(old$std_err)]
  expect_true(all(is.finite(std_err)))
  expect_true(all(std_err > 0))
  expect_equal(std_err, old$std_err, tolerance = 1e-6)

  invisible(result)
}

test_that("ipw_deli_sandwich() names its blocks for a binary outcome", {
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
  ee <- estimating_equations(fit)
  p <- length(ee@parameters)
  beta <- stats::coef(outcome_mod)

  result <- call_deli_sandwich(fit, outcome_mod, data)
  old <- old_stacked_pieces(fit, outcome_mod, data)

  expected_names <- c(
    paste0("theta_w", seq_len(p)),
    paste0("beta_", colnames(stats::model.matrix(outcome_mod))),
    "mu0",
    "mu1",
    "rd",
    "log(rr)",
    "log(or)"
  )
  expect_named(result$theta, expected_names)
  expect_identical(dimnames(result$vcov), list(expected_names, expected_names))

  expect_equal(unname(result$theta[seq_len(p)]), ee@parameters)
  expect_equal(unname(result$theta[p + seq_along(beta)]), unname(beta))
  expect_equal(result$theta[["mu0"]], old$mu0)
  expect_equal(result$theta[["mu1"]], old$mu1)
  expect_equal(result$theta[["rd"]], old$mu1 - old$mu0)
  expect_equal(result$theta[["log(rr)"]], log(old$mu1) - log(old$mu0))
  expect_equal(
    result$theta[["log(or)"]],
    stats::qlogis(old$mu1) - stats::qlogis(old$mu0)
  )
})

test_that("ipw_deli_sandwich() names its blocks for a continuous outcome", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y_cont ~ exposure, data, w, stats::gaussian())
  ee <- estimating_equations(fit)
  p <- length(ee@parameters)

  result <- call_deli_sandwich(fit, outcome_mod, data)
  old <- old_stacked_pieces(fit, outcome_mod, data)

  expected_names <- c(
    paste0("theta_w", seq_len(p)),
    paste0("beta_", colnames(stats::model.matrix(outcome_mod))),
    "mu0",
    "mu1",
    "diff"
  )
  expect_named(result$theta, expected_names)
  expect_identical(dimnames(result$vcov), list(expected_names, expected_names))
  expect_equal(result$theta[["diff"]], old$mu1 - old$mu0)
})

# The parity grid is the one the hand-assembled sandwich is already pinned on,
# widened to both outcome families. Each case asserts that swapping the variance
# engine changes nothing a caller can see.

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
  for (outcome in list(
    list(
      label = "a binary outcome",
      formula = y ~ exposure,
      family = quote(stats::binomial())
    ),
    list(
      label = "a continuous outcome",
      formula = y_cont ~ exposure,
      family = quote(stats::gaussian())
    )
  )) {
    local({
      spec <- spec
      outcome <- outcome
      test_that(
        paste0(
          "the deli sandwich matches the stacked sandwich for ",
          spec$label,
          " with ",
          outcome$label
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
          outcome_mod <- fit_outcome(
            outcome$formula,
            data,
            w,
            eval(outcome$family)
          )

          expect_deli_parity(fit, outcome_mod, data)
        }
      )
    })
  }
}

test_that("the deli sandwich matches the stacked sandwich with sampling weights", {
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

  expect_deli_parity(fit, outcome_mod, data)
})

# The scale-coherent oracle builds the whole stacked M-estimator at the scale
# the container stores natively, so it shares none of either engine's rescaling
# algebra. The new engine must reach it as closely as the old one does, up to
# the finite-difference bread.

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
        "the deli risk-difference standard error matches the coherent oracle for ",
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

        result <- call_deli_sandwich(fit, outcome_mod, data)
        rd_se <- sqrt(result$vcov[["rd", "rd"]])

        expect_equal(rd_se, oracle_se, tolerance = 1e-6)
      }
    )
  })
}

test_that("the deli risk-difference standard error is coherent with sampling weights", {
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

  result <- call_deli_sandwich(fit, outcome_mod, data)
  rd_se <- sqrt(result$vcov[["rd", "rd"]])

  expect_equal(rd_se, oracle_se, tolerance = 1e-6)
})

# A probit outcome model has a non-canonical link, where the expected (Fisher)
# information the hand-assembled bread uses and the observed information a
# finite-differenced bread produces are different matrices in general. They
# coincide here for a structural reason worth pinning: `ipw()` accepts only the
# marginal outcome model, whose sole predictor is the binary exposure, so the
# model is saturated. The term separating observed from expected information is
# a weighted sum of residuals within each exposure group scaled by a factor that
# is constant within the group, and the score equations of a saturated model set
# exactly those sums to zero. The non-canonical link therefore joins the parity
# grid rather than diverging from it.

test_that("the deli sandwich matches the stacked sandwich for a probit outcome model", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(
    y ~ exposure,
    data,
    w,
    stats::binomial(link = "probit")
  )

  expect_deli_parity(fit, outcome_mod, data)
})

# The bootstrap is the external check that the whole stack, not just its
# agreement with the previous engine, is calibrated. One resampling loop serves
# all three outcome models: the marginal means of a saturated weighted model are
# the weighted group means whatever the link, so the logit and probit fits share
# a risk-difference point estimate and therefore a bootstrap distribution, and
# the continuous outcome rides along on the same replicate fits.

test_that("the deli standard errors track a nonparametric bootstrap", {
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
  logit_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())
  probit_mod <- fit_outcome(
    y ~ exposure,
    data,
    w,
    stats::binomial(link = "probit")
  )
  continuous_mod <- fit_outcome(y_cont ~ exposure, data, w, stats::gaussian())

  logit <- call_deli_sandwich(fit, logit_mod, data)
  probit <- call_deli_sandwich(fit, probit_mod, data)
  continuous <- call_deli_sandwich(fit, continuous_mod, data)
  expect_equal(logit$theta[["rd"]], probit$theta[["rd"]], tolerance = 1e-8)

  n <- nrow(data)
  boot <- withr::with_seed(2024, {
    vapply(
      seq_len(200),
      function(b) {
        idx <- sample.int(n, n, replace = TRUE)
        resampled <- data[idx, , drop = FALSE]
        # A resampled data set can legitimately fail to converge; that replicate
        # drops out through the error handler, and its convergence warning is
        # suppressed so it does not leak into the suite output.
        tryCatch(
          suppressWarnings({
            boot_fit <- balance(
              resampled,
              exposure,
              c(x1, x2),
              method = bw_entropy(),
              estimand = "ate"
            )
            boot_w <- as.numeric(stats::weights(boot_fit))
            binary_mod <- fit_outcome(
              y ~ exposure,
              resampled,
              boot_w,
              stats::binomial()
            )
            gaussian_mod <- fit_outcome(
              y_cont ~ exposure,
              resampled,
              boot_w,
              stats::gaussian()
            )
            binary_means <- marginal_means(binary_mod, resampled)
            gaussian_means <- marginal_means(gaussian_mod, resampled)
            c(
              binary_means$mu1 - binary_means$mu0,
              gaussian_means$mu1 - gaussian_means$mu0
            )
          }),
          error = function(e) c(NA_real_, NA_real_)
        )
      },
      numeric(2)
    )
  })
  boot_rd_se <- stats::sd(boot[1L, ], na.rm = TRUE)
  boot_diff_se <- stats::sd(boot[2L, ], na.rm = TRUE)

  # The bootstrap is noisy at this replicate count, so the agreement is loose.
  expect_equal(sqrt(logit$vcov[["rd", "rd"]]), boot_rd_se, tolerance = 0.15)
  expect_equal(sqrt(probit$vcov[["rd", "rd"]]), boot_rd_se, tolerance = 0.15)
  expect_equal(
    sqrt(continuous$vcov[["diff", "diff"]]),
    boot_diff_se,
    tolerance = 0.15
  )
})

# A container carrying one extra parameter whose estimating function is
# identically zero. Its bread row is an exact zero row, so the stack does not
# identify that parameter. `deli::compute_sandwich()` pseudo-inverts a singular
# bread by default and would return a confident-looking covariance for a system
# that has none, so the call must pass `allow_pinv = FALSE` and let the failure
# surface. Building the deficiency into the container rather than mocking the
# deli call also pins that the weight block is sized from
# `container@parameters` and reached only through `psi_fn()` and `weights_fn()`.
inert_parameter_container <- function(ee) {
  p <- length(ee@parameters)
  jacobian <- matrix(0, p + 1L, p + 1L)
  jacobian[seq_len(p), seq_len(p)] <- ee@jacobian
  balancing_estimating_equations(
    parameters = c(ee@parameters, 0),
    psi = cbind(ee@psi, 0),
    jacobian = jacobian,
    weight_jacobian = cbind(ee@weight_jacobian, 0),
    weights_raw = ee@weights_raw,
    psi_fn = function(theta) cbind(ee@psi_fn(theta[seq_len(p)]), 0),
    weights_fn = function(theta) ee@weights_fn(theta[seq_len(p)])
  )
}

test_that("ipw_deli_sandwich() refuses a rank-deficient stack", {
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
    ipw_deli_sandwich(
      container = inert_parameter_container(estimating_equations(fit)),
      outcome_mod = outcome_mod,
      frame = data,
      exposure_name = fit@exposure,
      sampling_weights = fit@sampling_weights
    ),
    regexp = "singular"
  )
})
