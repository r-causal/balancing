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

# ---- The returned structure -----------------------------------------------

# The result carries the same fields propensity's own `ipw()` returns: which
# standard-error method produced the effect table, and the fitted variance
# system that produced it. The system here is the stacked parameter vector and
# its covariance, which is the whole of what the engine sandwiches, so `fit`
# holds those two rather than a solver object the method never built. The
# covariance is on the standard-error scale, so the reported standard error of
# each contrast is the square root of that contrast's diagonal entry, and the
# reported estimate is that contrast's entry in the parameter vector.

test_that("ipw() reports its standard-error method and the fitted variance system", {
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

  expect_identical(result$se_method, "mestimation")
  expect_named(result$fit, c("theta", "vcov"))

  p <- length(estimating_equations(fit)@parameters)
  expected_names <- c(
    paste0("theta_w", seq_len(p)),
    paste0("beta_", colnames(stats::model.matrix(outcome_mod))),
    "mu0",
    "mu1",
    "rd",
    "log(rr)",
    "log(or)"
  )
  expect_named(result$fit$theta, expected_names)
  expect_identical(
    dimnames(result$fit$vcov),
    list(expected_names, expected_names)
  )

  estimates <- as.data.frame(result)
  expect_equal(
    estimates$estimate,
    unname(result$fit$theta[estimates$effect]),
    tolerance = 1e-12
  )
  expect_equal(
    estimates$std.err,
    unname(sqrt(diag(result$fit$vcov))[estimates$effect]),
    tolerance = 1e-12
  )
})

test_that("the ipw() variance system carries one contrast for a continuous outcome", {
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

  result <- propensity::ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)

  expect_identical(result$se_method, "mestimation")
  expect_identical(
    utils::tail(names(result$fit$theta), 3L),
    c("mu0", "mu1", "diff")
  )
  expect_equal(
    estimates$std.err,
    unname(sqrt(result$fit$vcov[["diff", "diff"]])),
    tolerance = 1e-12
  )
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
    label = "an entropy att fit",
    method = quote(bw_entropy()),
    estimand = "att",
    focal = "1"
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

# ---- Weight consistency between the fit and the outcome model -------------

# The stacked variance differentiates the outcome-model score through the
# weights the fit produced, so it describes the system that was actually solved
# only when the outcome model was fitted at those weights. A model fitted at any
# other weights, or at none, leaves the point estimates and the standard errors
# internally inconsistent, and nothing downstream notices: the result still
# looks like a well-formed effect table. The preflight compares the model's
# prior weights against the fit's composed weights, which already carry the
# sampling weights, at a relative tolerance of 1e-6.

test_that("ipw() rejects an outcome model fitted without weights", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  outcome_mod <- stats::glm(
    y ~ exposure,
    data = data,
    family = stats::binomial()
  )

  # The remedy has to be actionable, so the message names the accessor that
  # produces the weights the fit expects rather than only reporting a mismatch.
  expect_error(
    propensity::ipw(fit, outcome_mod),
    class = "balancing_ipw_input_error",
    regexp = "weights(fit)",
    fixed = TRUE
  )

  cnd <- rlang::catch_cnd(
    propensity::ipw(fit, outcome_mod),
    classes = "balancing_ipw_input_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
})

# An lm fitted without weights records none at all, where a glm records a vector
# of ones. Reading the missing vector as ones is what makes an unweighted lm on
# a weighted fit an error rather than a case the preflight quietly skips.

test_that("ipw() rejects an unweighted lm on a weighted fit", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  outcome_mod <- stats::lm(y_cont ~ exposure, data = data)
  expect_null(stats::weights(outcome_mod))

  expect_error(
    propensity::ipw(fit, outcome_mod),
    class = "balancing_ipw_input_error"
  )
})

test_that("ipw() rejects an outcome model fitted with the wrong weights", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )

  noise <- withr::with_seed(909, stats::runif(nrow(data), 0.5, 2))
  noise_mod <- fit_outcome(y ~ exposure, data, noise, stats::binomial())
  expect_error(
    propensity::ipw(fit, noise_mod),
    class = "balancing_ipw_input_error"
  )

  # The realistic mistake is weights from a neighbouring fit rather than from
  # noise: they are close enough to the right ones that no result looks wrong.
  other_fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  other_mod <- fit_outcome(
    y ~ exposure,
    data,
    as.numeric(stats::weights(other_fit)),
    stats::binomial()
  )
  expect_error(
    propensity::ipw(fit, other_mod),
    class = "balancing_ipw_input_error"
  )
})

# ---- Offsets in the outcome model -----------------------------------------

# An offset is supported: the variance engine carries it through both the
# outcome-model score and the fixed-exposure linear predictors, so the marginal
# means are the g-computation means with each unit's offset held at its observed
# value. The point-estimate check is against a predict()-based g-computation,
# which is the independent statement of what the offset should do, and it is run
# through the public interface so the whole path is covered, including the model
# frame ipw() falls back on when no data are supplied. Both spellings are pinned:
# an offset written into the formula and one passed through the `offset`
# argument reach the fitted model by different routes.

test_that("ipw() honors an offset term in the outcome model", {
  data <- ipw_fixture()
  data$log_time <- withr::with_seed(11, stats::rnorm(nrow(data), 0, 0.3))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  binary_mod <- fit_outcome(
    y ~ exposure + offset(log_time),
    data,
    w,
    stats::binomial()
  )
  continuous_mod <- fit_outcome(
    y_cont ~ exposure + offset(log_time),
    data,
    w,
    stats::gaussian()
  )

  binary <- as.data.frame(propensity::ipw(fit, binary_mod))
  continuous <- as.data.frame(propensity::ipw(fit, continuous_mod))

  binary_means <- marginal_means(binary_mod, data)
  continuous_means <- marginal_means(continuous_mod, data)

  expect_identical(binary$effect, c("rd", "log(rr)", "log(or)"))
  expect_equal(
    binary$estimate[binary$effect == "rd"],
    binary_means$mu1 - binary_means$mu0,
    tolerance = 1e-10
  )
  expect_equal(
    binary$estimate[binary$effect == "log(rr)"],
    log(binary_means$mu1) - log(binary_means$mu0),
    tolerance = 1e-10
  )
  expect_identical(continuous$effect, "diff")
  expect_equal(
    continuous$estimate,
    continuous_means$mu1 - continuous_means$mu0,
    tolerance = 1e-10
  )
})

test_that("ipw() honors an offset argument in the outcome model", {
  data <- ipw_fixture()
  data$log_time <- withr::with_seed(11, stats::rnorm(nrow(data), 0, 0.3))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  data$.wts <- as.numeric(stats::weights(fit))
  outcome_mod <- suppressWarnings(stats::glm(
    y ~ exposure,
    data = data,
    family = stats::binomial(),
    weights = .wts,
    offset = log_time
  ))

  estimates <- as.data.frame(propensity::ipw(fit, outcome_mod))
  means <- marginal_means(outcome_mod, data)

  expect_equal(
    estimates$estimate[estimates$effect == "rd"],
    means$mu1 - means$mu0,
    tolerance = 1e-10
  )
})

# The standard errors an offset model reports are the variance engine's own, so
# they are pinned against the engine rather than re-derived. The engine is
# already bootstrap-pinned for an offset model, which is what makes the identity
# sufficient here.

test_that("ipw() standard errors with an offset come from the variance engine", {
  data <- ipw_fixture()
  data$log_time <- withr::with_seed(11, stats::rnorm(nrow(data), 0, 0.3))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(
    y ~ exposure + offset(log_time),
    data,
    w,
    stats::binomial()
  )

  estimates <- as.data.frame(propensity::ipw(fit, outcome_mod))
  engine <- ipw_deli_sandwich(
    container = estimating_equations(fit),
    outcome_mod = outcome_mod,
    frame = data,
    exposure_name = fit@exposure,
    sampling_weights = fit@sampling_weights
  )

  expect_true(all(is.finite(estimates$std.err)))
  expect_true(all(estimates$std.err > 0))
  expect_equal(
    estimates$std.err,
    unname(sqrt(diag(engine$vcov))[estimates$effect]),
    tolerance = 1e-12
  )
})

# ---- Outcome model family -------------------------------------------------

# The reported effects are contrasts of two marginal means read as
# probabilities: a risk difference, a log risk ratio, and a log odds ratio. A
# count outcome has marginal means that are rates, so the odds ratio is not
# defined for it and the method returns NaN for that row while still labelling
# the first row a risk difference. Rather than report a table whose labels do
# not describe its contents, the accepted families are restricted to the ones
# whose marginal means the effect rows are derived for.

test_that("ipw() rejects a poisson outcome model", {
  data <- ipw_fixture()
  data$y_count <- withr::with_seed(
    13,
    stats::rpois(nrow(data), exp(0.2 + 0.3 * data$exposure))
  )
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y_count ~ exposure, data, w, stats::poisson())

  expect_error(
    propensity::ipw(fit, outcome_mod),
    class = "balancing_ipw_input_error"
  )

  cnd <- rlang::catch_cnd(
    propensity::ipw(fit, outcome_mod),
    classes = "balancing_ipw_input_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
})

# The restriction is an allowlist rather than a list of known-bad families, so
# a quasi family and a family with a positive continuous response are refused
# for the same reason the count family is: their marginal means are not the
# quantities the effect rows contrast.

test_that("ipw() rejects the quasipoisson and inverse gaussian families", {
  data <- ipw_fixture()
  data$y_count <- withr::with_seed(
    13,
    stats::rpois(nrow(data), exp(0.2 + 0.3 * data$exposure))
  )
  data$y_pos <- withr::with_seed(
    17,
    stats::rgamma(nrow(data), shape = 2, rate = 1) + 0.1
  )
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))

  quasipoisson_mod <- fit_outcome(
    y_count ~ exposure,
    data,
    w,
    stats::quasipoisson()
  )
  inverse_mod <- fit_outcome(
    y_pos ~ exposure,
    data,
    w,
    stats::inverse.gaussian(link = "log")
  )

  expect_error(
    propensity::ipw(fit, quasipoisson_mod),
    class = "balancing_ipw_input_error"
  )
  expect_error(
    propensity::ipw(fit, inverse_mod),
    class = "balancing_ipw_input_error"
  )
})

# The restriction must not over-reach. A binomial or quasibinomial glm, a
# gaussian glm, and a plain lm all stay supported, so each is pinned beside the
# rejection. The quasibinomial case is the one most easily lost to a family
# check written against `binomial` alone: it shares the binomial fit's
# coefficients and its variance function, and its dispersion never enters the
# sandwich, so it must return the binomial answer exactly.

test_that("ipw() accepts the binomial, quasibinomial, gaussian, and lm families", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  data$.wts <- w

  binomial_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())
  quasi_mod <- fit_outcome(y ~ exposure, data, w, stats::quasibinomial())
  gaussian_mod <- fit_outcome(y_cont ~ exposure, data, w, stats::gaussian())
  lm_mod <- stats::lm(y_cont ~ exposure, data = data, weights = .wts)

  binomial_result <- as.data.frame(propensity::ipw(fit, binomial_mod))
  quasi_result <- as.data.frame(propensity::ipw(fit, quasi_mod))
  gaussian_result <- as.data.frame(propensity::ipw(fit, gaussian_mod))
  lm_result <- as.data.frame(propensity::ipw(fit, lm_mod))

  expect_identical(binomial_result$effect, c("rd", "log(rr)", "log(or)"))
  expect_identical(quasi_result$effect, c("rd", "log(rr)", "log(or)"))
  expect_identical(gaussian_result$effect, "diff")
  expect_identical(lm_result$effect, "diff")

  expect_equal(quasi_result$estimate, binomial_result$estimate)
  expect_equal(quasi_result$std.err, binomial_result$std.err)
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

# The variance engine reaches the weight path only through the container's two
# re-evaluation hooks, so a container carrying neither cannot be differentiated
# at all. Every fit that reaches ipw() today populates both, and the binary gate
# turns away the one family whose container could arrive without them, so the
# guard is defence in depth. What it buys is a comprehensible refusal in place of
# a failure raised from inside the engine against a NULL it was handed.

test_that("ipw() rejects a fit whose container carries no re-evaluation hooks", {
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
  without_hooks <- function(psi_fn, weights_fn) {
    doctored <- fit
    doctored@estimating_equations <- balancing_estimating_equations(
      parameters = ee@parameters,
      psi = ee@psi,
      jacobian = ee@jacobian,
      weight_jacobian = ee@weight_jacobian,
      weights_raw = ee@weights_raw,
      psi_fn = psi_fn,
      weights_fn = weights_fn
    )
    doctored
  }

  # Either hook alone is not enough, so each absence is pinned separately.
  no_psi <- without_hooks(NULL, ee@weights_fn)
  no_weights <- without_hooks(ee@psi_fn, NULL)
  expect_error(
    propensity::ipw(no_psi, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  )
  expect_error(
    propensity::ipw(no_weights, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  )

  cnd <- rlang::catch_cnd(
    propensity::ipw(without_hooks(NULL, NULL), outcome_mod),
    classes = "balancing_ipw_unsupported_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
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

# `ipw_deli_sandwich()` is the variance engine `ipw()` reads from. It evaluates
# the whole stacked system as one p-by-n estimating-function closure and lets
# `deli::compute_sandwich()` finite difference the bread, so the R layer writes
# no derivative by hand and the effect contrasts are parameters of the stack
# rather than a delta-method gradient applied afterwards.
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

# Call the engine with the pieces read off a fit, the way the `ipw()` method
# does.
call_deli_sandwich <- function(fit, outcome_mod, data) {
  ipw_deli_sandwich(
    container = estimating_equations(fit),
    outcome_mod = outcome_mod,
    frame = data,
    exposure_name = fit@exposure,
    sampling_weights = fit@sampling_weights
  )
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
  means <- marginal_means(outcome_mod, data)

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
  expect_equal(result$theta[["mu0"]], means$mu0)
  expect_equal(result$theta[["mu1"]], means$mu1)
  expect_equal(result$theta[["rd"]], means$mu1 - means$mu0)
  expect_equal(result$theta[["log(rr)"]], log(means$mu1) - log(means$mu0))
  expect_equal(
    result$theta[["log(or)"]],
    stats::qlogis(means$mu1) - stats::qlogis(means$mu0)
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
  means <- marginal_means(outcome_mod, data)

  expected_names <- c(
    paste0("theta_w", seq_len(p)),
    paste0("beta_", colnames(stats::model.matrix(outcome_mod))),
    "mu0",
    "mu1",
    "diff"
  )
  expect_named(result$theta, expected_names)
  expect_identical(dimnames(result$vcov), list(expected_names, expected_names))
  expect_equal(result$theta[["diff"]], means$mu1 - means$mu0)
})

# The scale-coherent oracle builds the whole stacked M-estimator at the scale
# the container stores natively, so it shares none of the engine's rescaling
# algebra. The tolerance is 1e-6 rather than the 1e-8 the public grid uses,
# since the engine is reached here directly and the finite-differenced bread is
# the only source of disagreement.

for (spec in list(
  list(
    label = "an entropy ate fit",
    method = quote(bw_entropy()),
    estimand = "ate",
    focal = NULL
  ),
  list(
    label = "an entropy att fit",
    method = quote(bw_entropy()),
    estimand = "att",
    focal = "1"
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

# The bootstrap is the external check that the whole stack is calibrated. One
# resampling loop serves all three outcome models: the marginal means of a
# saturated weighted model are the weighted group means whatever the link, so
# the logit and probit fits share a risk-difference point estimate and therefore
# a bootstrap distribution, and the continuous outcome rides along on the same
# replicate fits. The probit fit is also the non-canonical link in the suite,
# which the engine handles because it differentiates the score it is given
# rather than assuming an information matrix.

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

# A container whose estimating functions go missing away from the solution
# leaves the finite-differenced bread full of missing values. deli answers that
# with a warning and a `NULL` rather than an error, which would otherwise
# surface much later as a complaint about dimnames applied to a non-array, so
# the engine names the cause where it is still legible. deli's warning is its
# own unclassed one, so it is suppressed rather than pinned: pinning it would
# tie this test to another package's wording.
na_psi_container <- function(ee) {
  balancing_estimating_equations(
    parameters = ee@parameters,
    psi = ee@psi,
    jacobian = ee@jacobian,
    weight_jacobian = ee@weight_jacobian,
    weights_raw = ee@weights_raw,
    psi_fn = function(theta) {
      psi <- ee@psi_fn(theta)
      if (!isTRUE(all.equal(unname(theta), unname(ee@parameters)))) {
        psi[] <- NA_real_
      }
      psi
    },
    weights_fn = ee@weights_fn
  )
}

test_that("ipw_deli_sandwich() refuses a stack whose bread is not finite", {
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
    suppressWarnings(ipw_deli_sandwich(
      container = na_psi_container(estimating_equations(fit)),
      outcome_mod = outcome_mod,
      frame = data,
      exposure_name = fit@exposure,
      sampling_weights = fit@sampling_weights
    )),
    class = "balancing_ipw_unsupported_error"
  )
})

# ---- Hook evaluations through the deli engine ------------------------------

# The stacked closure reaches the weight path only through the container's two
# hooks, and both read the weight block alone. A central difference evaluates
# the closure once at the fitted parameters and twice more per stacked
# coordinate, but a perturbation of a coordinate outside the weight block leaves
# the weight sub-vector at its fitted value. Only the fitted sub-vector and one
# up and one down perturbation per weight coordinate therefore ever reach the
# hooks, so each hook has `2 * p + 1` distinct arguments no matter how long the
# rest of the stack is. Calling them more often than that is repeated work on
# arguments already seen, and it is the container hooks, not the rest of the
# closure, that carry the cost: each one crosses into the method's solver
# entrypoint over the full data.
#
# A container wrapping a real one's hooks with a counter records both how often
# each hook ran and which sub-vectors it saw. The parity assertions against the
# undoctored container are what keep the economy honest: skipping a call is only
# admissible when the answer is bit-for-bit the one the call would have given.
counting_container <- function(ee) {
  calls <- new.env(parent = emptyenv())
  calls$psi <- 0L
  calls$weights <- 0L
  calls$seen <- list()
  container <- balancing_estimating_equations(
    parameters = ee@parameters,
    psi = ee@psi,
    jacobian = ee@jacobian,
    weight_jacobian = ee@weight_jacobian,
    weights_raw = ee@weights_raw,
    psi_fn = function(theta) {
      calls$psi <- calls$psi + 1L
      calls$seen <- c(calls$seen, list(unname(as.numeric(theta))))
      ee@psi_fn(theta)
    },
    weights_fn = function(theta) {
      calls$weights <- calls$weights + 1L
      ee@weights_fn(theta)
    }
  )
  list(container = container, calls = calls)
}

test_that("the deli sandwich evaluates each hook once per distinct weight vector", {
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
  counted <- counting_container(ee)

  plain <- call_deli_sandwich(fit, outcome_mod, data)
  doctored <- ipw_deli_sandwich(
    container = counted$container,
    outcome_mod = outcome_mod,
    frame = data,
    exposure_name = fit@exposure,
    sampling_weights = fit@sampling_weights
  )

  expect_identical(doctored$theta, plain$theta)
  expect_identical(doctored$vcov, plain$vcov)

  expect_identical(length(unique(counted$calls$seen)), 2L * p + 1L)
  expect_identical(counted$calls$psi, 2L * p + 1L)
  expect_identical(counted$calls$weights, 2L * p + 1L)
})

# ---- Outcome families through the deli engine ------------------------------

# deli has no quasibinomial estimating equation and does not need one. The
# quasibinomial variance function is the binomial one, and the dispersion the
# quasi family estimates never enters the sandwich, which is built from the
# score equations alone. Mapping quasibinomial onto deli's binomial therefore
# has to reproduce the binomial result rather than merely approximate it: the
# two fits share their coefficients, so every block of the stack coincides.

test_that("the deli sandwich treats a quasibinomial outcome model as binomial", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  binomial_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())
  quasi_mod <- fit_outcome(y ~ exposure, data, w, stats::quasibinomial())

  binomial_result <- call_deli_sandwich(fit, binomial_mod, data)
  quasi_result <- call_deli_sandwich(fit, quasi_mod, data)

  expect_named(quasi_result$theta, names(binomial_result$theta))
  expect_equal(quasi_result$theta, binomial_result$theta, tolerance = 1e-10)
  expect_equal(quasi_result$vcov, binomial_result$vcov, tolerance = 1e-10)
})

# A gamma fit estimates a dispersion alongside the coefficients in deli's
# estimating equation, which reads the last element of the parameter vector as a
# log dispersion and returns an extra row. Passed the plain coefficient vector
# the stack carries, it would build a wrong-shaped block out of a misread
# parameter, and no downstream check would notice, so the family is refused
# before the closure is ever evaluated.

test_that("the deli sandwich rejects a gamma outcome model", {
  data <- ipw_fixture()
  data$y_pos <- withr::with_seed(
    17,
    stats::rgamma(nrow(data), shape = 2, rate = 1) + 0.1
  )
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(
    y_pos ~ exposure,
    data,
    w,
    stats::Gamma(link = "log")
  )

  expect_error(
    call_deli_sandwich(fit, outcome_mod, data),
    class = "balancing_ipw_unsupported_error"
  )
})

# A negative binomial fit reports its family as "Negative Binomial(theta)", with
# the estimated dispersion spelled into the family name itself. A refusal
# written as an exact string comparison against "negative_binomial" therefore
# never fires, and the call falls through to deli's own unclassed complaint
# about an unknown distribution. The refusal has to match on the family prefix,
# which is what this pins. The dispersion the fit estimates is exactly the
# parameter the stack does not carry, so the family belongs with gamma.

test_that("the deli sandwich rejects a negative binomial outcome model", {
  skip_if_not_installed("MASS")
  data <- ipw_fixture()
  data$y_count <- withr::with_seed(
    13,
    stats::rpois(nrow(data), exp(0.2 + 0.3 * data$exposure))
  )
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  data$.wts <- as.numeric(stats::weights(fit))
  outcome_mod <- MASS::glm.nb(y_count ~ exposure, data = data, weights = .wts)

  expect_match(stats::family(outcome_mod)$family, "^Negative Binomial")
  expect_error(
    call_deli_sandwich(fit, outcome_mod, data),
    class = "balancing_ipw_unsupported_error"
  )
})

# ---- Offsets through the deli engine ---------------------------------------

# The deli engine carries an offset through both the outcome-model score and the
# fixed-exposure linear predictors, which is the capability the user-facing
# rejection of offsets is waiting on. Three pins bracket it. An offset that is
# identically zero must change nothing at all, since it enters only as an
# addition of zero to the linear predictor. A constant offset must move the
# intercept and nothing else, because the model is a reparametrization of the
# same fit. A genuinely unit-varying offset must reach the same marginal means a
# predict()-based g-computation does, which is the check that the offset is
# actually carried into the fixed-exposure predictions rather than dropped.

test_that("the deli sandwich is unchanged by an offset of zero", {
  data <- ipw_fixture()
  data$zero <- 0
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))

  plain <- fit_outcome(y ~ exposure, data, w, stats::binomial())
  offset_mod <- fit_outcome(
    y ~ exposure + offset(zero),
    data,
    w,
    stats::binomial()
  )

  plain_result <- call_deli_sandwich(fit, plain, data)
  offset_result <- call_deli_sandwich(fit, offset_mod, data)

  expect_identical(offset_result$theta, plain_result$theta)
  expect_identical(offset_result$vcov, plain_result$vcov)
})

test_that("the deli sandwich absorbs a constant offset into the intercept", {
  data <- ipw_fixture()
  shift <- 0.4
  data$shift <- shift
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))

  plain <- fit_outcome(y ~ exposure, data, w, stats::binomial())
  offset_mod <- fit_outcome(
    y ~ exposure + offset(shift),
    data,
    w,
    stats::binomial()
  )

  # The reparametrization is exact at the coefficient level, which is what makes
  # every quantity the stack reports comparable between the two fits.
  expect_equal(
    stats::coef(offset_mod)[["(Intercept)"]],
    stats::coef(plain)[["(Intercept)"]] - shift,
    tolerance = 1e-8
  )
  expect_equal(
    stats::coef(offset_mod)[["exposure"]],
    stats::coef(plain)[["exposure"]],
    tolerance = 1e-8
  )

  plain_result <- call_deli_sandwich(fit, plain, data)
  offset_result <- call_deli_sandwich(fit, offset_mod, data)

  # The two fits differ only in where the constant sits, so the means, the
  # contrasts, and their standard errors must agree to well past the
  # finite-difference bread's own accuracy.
  reported <- c("mu0", "mu1", "rd", "log(rr)", "log(or)")
  expect_equal(
    offset_result$theta[reported],
    plain_result$theta[reported],
    tolerance = 1e-8
  )
  expect_equal(
    sqrt(diag(offset_result$vcov))[reported],
    sqrt(diag(plain_result$vcov))[reported],
    tolerance = 1e-8
  )
})

test_that("the deli sandwich carries a unit-varying offset into the means", {
  data <- ipw_fixture()
  data$log_time <- withr::with_seed(11, stats::rnorm(nrow(data), 0, 0.3))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))

  binary_mod <- fit_outcome(
    y ~ exposure + offset(log_time),
    data,
    w,
    stats::binomial()
  )
  continuous_mod <- fit_outcome(
    y_cont ~ exposure + offset(log_time),
    data,
    w,
    stats::gaussian()
  )

  binary_result <- call_deli_sandwich(fit, binary_mod, data)
  continuous_result <- call_deli_sandwich(fit, continuous_mod, data)

  binary_means <- marginal_means(binary_mod, data)
  continuous_means <- marginal_means(continuous_mod, data)

  expect_equal(
    binary_result$theta[["mu0"]],
    binary_means$mu0,
    tolerance = 1e-10
  )
  expect_equal(
    binary_result$theta[["mu1"]],
    binary_means$mu1,
    tolerance = 1e-10
  )
  expect_equal(
    binary_result$theta[["rd"]],
    binary_means$mu1 - binary_means$mu0,
    tolerance = 1e-10
  )
  expect_equal(
    continuous_result$theta[["diff"]],
    continuous_means$mu1 - continuous_means$mu0,
    tolerance = 1e-10
  )

  binary_se <- sqrt(diag(binary_result$vcov))
  continuous_se <- sqrt(diag(continuous_result$vcov))
  expect_true(all(is.finite(binary_se)))
  expect_true(all(binary_se > 0))
  expect_true(all(is.finite(continuous_se)))
  expect_true(all(continuous_se > 0))
})
