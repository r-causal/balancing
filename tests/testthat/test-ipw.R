# Specs for the ipw() integration: the estimating-equations contract that
# variance estimation depends on, and balancing's method on the causalgenerics
# ipw() generic that the package re-exports. ipw() dispatches on the balancing
# fit and computes a stacked-sandwich variance from the fit's
# estimating-equations container.

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

# A three-level categorical fixture with a binary and a continuous outcome, both
# depending on the exposure level and on x1, so every level's marginal mean is
# distinct and an adjusted model has something to adjust for. The exposure comes
# from the shared `sim_categorical()` data-generating process; the outcomes are
# drawn here under their own seed so the fixture is fixed without changing the
# shared helper.
ipw_categorical_fixture <- function(n = 200) {
  data <- sim_categorical(n)
  withr::with_seed(303, {
    is_b <- as.numeric(data$exposure == "b")
    is_c <- as.numeric(data$exposure == "c")
    linear_predictor <- -0.4 + 0.7 * is_b + 1.1 * is_c + 0.5 * data$x1
    data$y <- stats::rbinom(n, 1L, stats::plogis(linear_predictor))
    data$y_cont <- 0.5 +
      0.8 * is_b -
      0.4 * is_c +
      0.6 * data$x1 -
      0.3 * data$x2 +
      stats::rnorm(n)
  })
  data
}

# The K marginal means ipw() reports for a categorical exposure: predict the
# outcome model with the exposure fixed to each level in turn and average over
# the target population. `tilt` is that population's weight, one per unit, which
# is uniform for a pooled estimand and the focal group's indicator for a focal
# one. The result is named by level, in the exposure's own level order, so the
# reference level is its first element.
categorical_marginal_means <- function(
  outcome_mod,
  data,
  tilt = NULL,
  exposure_name = "exposure"
) {
  levels <- levels(data[[exposure_name]])
  weight <- tilt %||% rep(1, nrow(data))
  vapply(
    levels,
    function(level) {
      counterfactual <- data
      counterfactual[[exposure_name]] <- factor(level, levels = levels)
      stats::weighted.mean(
        stats::predict(
          outcome_mod,
          newdata = counterfactual,
          type = "response"
        ),
        weight
      )
    },
    numeric(1)
  )
}

# The categorical counterpart of `coherent_rd_se()`, and coherent in the same
# sense: it builds the whole stacked M-estimator at the scale the container
# stores natively, with an analytic bread, so it shares none of the
# implementation's rescaling algebra. The stack carries one marginal-mean row per
# exposure level, and the risk difference of `level` against the reference level
# is read off the joint covariance through the contrast vector rather than by a
# delta method. It covers the marginal binomial model at a pooled estimand, where
# the model is saturated in the exposure and a per-group rescale of the weights
# therefore leaves the coefficients alone, which is what makes the two scales
# agree exactly.
coherent_categorical_rd_se <- function(fit, data, level, sampling = NULL) {
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
  coefficients <- stats::coef(outcome_mod)
  eta <- as.numeric(design %*% coefficients)
  mu <- family$linkinv(eta)
  mu_eta <- family$mu.eta(eta)
  variance <- family$variance(mu)
  residual <- (data$y - mu) * mu_eta / variance
  working_weight <- composed * mu_eta^2 / variance
  score <- (composed * residual) * design

  levels <- levels(data$exposure)
  k <- length(levels)
  terms <- stats::delete.response(stats::terms(outcome_mod))
  designs <- lapply(levels, function(l) {
    counterfactual <- data
    counterfactual$exposure <- factor(l, levels = levels)
    stats::model.matrix(terms, stats::model.frame(terms, counterfactual))
  })
  predictions <- lapply(designs, function(x) {
    family$linkinv(as.numeric(x %*% coefficients))
  })
  means <- vapply(predictions, mean, numeric(1))

  psi <- ee@psi
  p <- ncol(psi)
  mean_residuals <- do.call(
    cbind,
    lapply(seq_len(k), function(j) predictions[[j]] - means[[j]])
  )
  stacked <- cbind(psi, score, mean_residuals)
  meat <- crossprod(stacked) / n

  size <- p + q + k
  theta <- seq_len(p)
  beta <- p + seq_len(q)
  mu_index <- p + q + seq_len(k)
  bread <- matrix(0, size, size)
  bread[theta, theta] <- ee@jacobian / n
  bread[beta, theta] <- crossprod(design * residual, s * ee@weight_jacobian) / n
  bread[beta, beta] <- -crossprod(design * working_weight, design) / n
  for (j in seq_len(k)) {
    eta_j <- as.numeric(designs[[j]] %*% coefficients)
    bread[mu_index[[j]], beta] <- colSums(family$mu.eta(eta_j) * designs[[j]]) /
      n
    bread[mu_index[[j]], mu_index[[j]]] <- -1
  }

  bread_inv <- solve(bread)
  cov <- bread_inv %*% meat %*% t(bread_inv) / n
  contrast <- numeric(size)
  contrast[mu_index[[match(level, levels)]]] <- 1
  contrast[mu_index[[1]]] <- -1
  sqrt(as.numeric(t(contrast) %*% cov %*% contrast))
}

# A continuous-exposure fixture carrying one outcome per link the continuous
# effect naming distinguishes. The exposure and the balanced covariates come
# from the shared `sim_continuous()` data-generating process; the outcomes are
# drawn here under their own seed so the fixture is fixed without changing the
# shared helper.
#
# `v` is prognostic for the continuous outcome and is not among the balanced
# covariates, which is what makes an adjusted marginal structural model
# distinguishable from a marginal one. A covariate the fit balanced is
# orthogonal to the exposure under the fitted weights, since driving that
# correlation to zero is the constraint, so adjusting for one leaves the
# exposure coefficient and its standard error where the marginal model puts
# them.
#
# `y_rare` is drawn from a log risk model rather than a logistic one because the
# log-link fit that reports a log risk ratio has to keep every fitted risk below
# one, which a rare outcome leaves room for.
ipw_continuous_fixture <- function(n = 300) {
  data <- sim_continuous(n)
  withr::with_seed(404, {
    data$v <- stats::rnorm(n)
    data$y_cont <- 1 +
      0.4 * data$exposure +
      0.5 * data$x1 -
      0.3 * data$x2 +
      0.9 * data$v +
      stats::rnorm(n)
    data$y <- stats::rbinom(
      n,
      1L,
      stats::plogis(-0.4 + 0.5 * data$exposure + 0.3 * data$x1)
    )
    data$y_rare <- stats::rbinom(n, 1L, exp(-1.9 + 0.15 * data$exposure))
  })
  data
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

test_that("categorical fits expose a container carrying both hooks", {
  data <- sim_categorical(200)
  for (method in list(bw_ipt(), bw_cbps(), bw_entropy())) {
    fit <- balance(data, exposure, c(x1, x2), method = method, estimand = "ate")
    ee <- fit@estimating_equations
    expect_false(is.null(ee))
    expect_false(is.null(ee@psi_fn))
    expect_false(is.null(ee@weights_fn))
  }
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

  result <- ipw(fit, outcome_mod)
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

  result <- ipw(fit, outcome_mod)
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

  result <- ipw(fit, outcome_mod)

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

  result <- ipw(fit, outcome_mod)
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

  result <- ipw(fit, outcome_mod)

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

  result <- ipw(fit, outcome_mod)
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

  result <- ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)

  expect_true(all(is.finite(estimates$std.err)))
  expect_true(all(estimates$std.err > 0))
})

test_that("a shift-related covariate leaves the ipw() chain identified", {
  # A covariate that is another covariate plus a constant contributes the same
  # constraint column once the columns are centered, so the expansion drops it
  # and the fit is the fit on the source covariate alone. Carrying both leaves
  # the estimating equations rank deficient and the stacked variance resting on
  # a singular bread.
  data <- ipw_fixture()
  data$x1_shifted <- data$x1 + 5

  fit <- balance(
    data,
    exposure,
    c(x1, x1_shifted),
    method = bw_entropy(),
    estimand = "ate"
  )
  reduced <- balance(
    data,
    exposure,
    x1,
    method = bw_entropy(),
    estimand = "ate"
  )

  expect_length(fit@recipe, 1L)
  jacobian <- estimating_equations(fit)@jacobian
  expect_identical(qr(jacobian)$rank, ncol(jacobian))

  outcome_mod <- fit_outcome(
    y ~ exposure,
    data,
    as.numeric(stats::weights(fit)),
    stats::binomial()
  )
  reduced_mod <- fit_outcome(
    y ~ exposure,
    data,
    as.numeric(stats::weights(reduced)),
    stats::binomial()
  )
  estimates <- as.data.frame(ipw(fit, outcome_mod))
  reduced_estimates <- as.data.frame(ipw(reduced, reduced_mod))

  expect_true(all(is.finite(estimates$std.err)))
  expect_true(all(estimates$std.err > 0))
  expect_equal(estimates$estimate, reduced_estimates$estimate)
  expect_equal(estimates$std.err, reduced_estimates$std.err)
})

test_that("a factor covariate leaves the ipw() sandwich finite", {
  # A factor's level indicators sum to the constant function. Nothing in the raw
  # column rank marks that redundancy, so every level keeps its own constraint
  # column and the entropy estimating equations are rank deficient: the Jacobian
  # block of each solved group has the level-sum direction in its null space. The
  # redundancy is harmless because the weight map is flat along the same
  # direction, so the fit and its influence function are the fit and influence
  # function of the parameterization that drops one level column, and the stacked
  # variance must agree with that reduced fit rather than dissolve into the
  # singularity.
  data <- sim_binary()
  withr::with_seed(11, {
    data$y <- stats::rbinom(
      nrow(data),
      1L,
      stats::plogis(-0.3 + 0.5 * data$exposure + 0.4 * data$x1)
    )
  })
  data$x3_b <- as.numeric(data$x3 == "b")
  data$x3_c <- as.numeric(data$x3 == "c")

  fit <- expect_no_warning(
    balance(
      data,
      exposure,
      c(x1, x2, x3),
      method = bw_entropy(),
      estimand = "ate"
    ),
    class = "balancing_convergence_warning"
  )
  reduced <- balance(
    data,
    exposure,
    c(x1, x2, x3_b, x3_c),
    method = bw_entropy(),
    estimand = "ate"
  )

  jacobian <- estimating_equations(fit)@jacobian
  expect_lt(qr(jacobian)$rank, ncol(jacobian))
  expect_equal(
    as.numeric(stats::weights(fit)),
    as.numeric(stats::weights(reduced)),
    tolerance = 1e-8
  )

  outcome_mod <- fit_outcome(
    y ~ exposure,
    data,
    as.numeric(stats::weights(fit)),
    stats::binomial()
  )
  reduced_mod <- fit_outcome(
    y ~ exposure,
    data,
    as.numeric(stats::weights(reduced)),
    stats::binomial()
  )
  # The deficiency is tolerated rather than refused, and tolerated silently: the
  # weight map is flat along it, so it never reaches the reported effects.
  estimates <- as.data.frame(expect_no_warning(ipw(fit, outcome_mod)))
  reduced_estimates <- as.data.frame(ipw(reduced, reduced_mod))

  expect_true(all(is.finite(estimates$std.err)))
  expect_true(all(estimates$std.err > 0))
  expect_equal(estimates$estimate, reduced_estimates$estimate, tolerance = 1e-6)
  expect_equal(estimates$std.err, reduced_estimates$std.err, tolerance = 1e-6)
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

  result <- ipw(fit, outcome_mod)
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

        result <- ipw(fit, outcome_mod)
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

  result <- ipw(fit, outcome_mod)
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

  result <- ipw(fit, outcome_mod)
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

# ---- Covariate-adjusted outcome models ------------------------------------

# The outcome model may adjust for covariates alongside the exposure. The
# reported effects are still the g-computation contrasts, with one difference
# the marginal form hides. A saturated marginal model predicts a single value
# per exposure level, so every population of units averages those predictions to
# the same pair of means; an adjusted model predicts a value per unit, and the
# population those predictions are averaged over is then part of the estimand.
# The marginal means therefore standardize over the estimand's target
# population: every unit for the average treatment effect, the focal group's
# units for a focal estimand, sampling-weighted in both cases. This is
# propensity's tilted g-computation with the tilt read off the data rather than
# off a model: propensity standardizes with the tilting function h of its fitted
# propensity score, and a balancing fit carries no propensity model, so the
# tilt here is the target population's indicator, which is parameter-free.

# The tilt-standardized g-computation, written with `predict()` and plain
# arithmetic so it shares nothing with the engine's own fixed-exposure design. A
# `NULL` focal level standardizes over every unit, which is the average
# treatment effect; a focal level standardizes over that group alone.
adjusted_means <- function(
  outcome_mod,
  data,
  focal_level = NULL,
  sampling = NULL
) {
  levels <- sort(unique(data$exposure))
  d0 <- d1 <- data
  d0$exposure <- levels[[1]]
  d1$exposure <- levels[[2]]
  indicator <- if (is.null(focal_level)) {
    rep(1, nrow(data))
  } else {
    as.numeric(as.character(data$exposure) == focal_level)
  }
  weight <- indicator * (sampling %||% rep(1, nrow(data)))
  standardize <- function(newdata) {
    prediction <- stats::predict(
      outcome_mod,
      newdata = newdata,
      type = "response"
    )
    sum(weight * prediction) / sum(weight)
  }
  list(mu0 = standardize(d0), mu1 = standardize(d1))
}

# An independent risk-difference standard error for an adjusted outcome model,
# built by hand from the container's stored matrices. The marginal oracle above
# builds its stack at the scale the container stores natively, which it may do
# because a per-group rescale of the weights leaves a marginal model's
# coefficients alone. An adjusted model's coefficients move under that rescale,
# so this oracle has to be built at the scale the outcome model was actually
# fitted at, the reported one, and therefore has to differentiate the reported
# weight map itself.
#
# That map carries each exposure group's weights to a target total that does not
# depend on the parameters, w_i(theta) = raw_i(theta) * target_g / sum_g(theta),
# so its derivative is the raw derivative rescaled, less the term that comes from
# the group sum moving with the parameters:
#
#   d w_i / d theta = c_g * d raw_i / d theta - w_i * (d sum_g / d theta) / sum_g
#
# with c_g the per-group reporting factor. Dropping the second term is the
# tempting shortcut, since it vanishes for a marginal model, and this oracle
# exists to say what dropping it costs: the term is what separates the two
# candidate couplings, and no bootstrap of a feasible size resolves the
# difference between them.
adjusted_rd_se <- function(fit, outcome_mod, data, sampling = NULL) {
  ee <- estimating_equations(fit)
  n <- nrow(data)
  s <- sampling %||% rep(1, n)
  reported <- as.numeric(stats::weights(fit, include_sampling_weights = FALSE))
  composed <- s * reported

  family <- stats::family(outcome_mod)
  design <- stats::model.matrix(outcome_mod)
  q <- ncol(design)
  beta <- stats::coef(outcome_mod)
  eta <- as.numeric(design %*% beta)
  mu <- family$linkinv(eta)
  mu_eta <- family$mu.eta(eta)
  variance <- family$variance(mu)
  residual <- (as.numeric(outcome_mod$y) - mu) * mu_eta / variance
  working_weight <- composed * mu_eta^2 / variance
  score <- (composed * residual) * design

  groups <- split(seq_len(n), as.character(data[[fit@exposure]]))
  raw <- ee@weights_raw
  ratio <- reported / raw
  ratio[!is.finite(ratio)] <- 1
  weight_derivative <- matrix(0, n, ncol(ee@weight_jacobian))
  for (level in names(groups)) {
    idx <- groups[[level]]
    group_sum <- sum(s[idx] * raw[idx])
    group_derivative <- colSums(
      s[idx] * ee@weight_jacobian[idx, , drop = FALSE]
    )
    weight_derivative[idx, ] <- ratio[idx] *
      ee@weight_jacobian[idx, , drop = FALSE] -
      outer(reported[idx], group_derivative) / group_sum
  }

  levels <- sort(unique(data[[fit@exposure]]))
  terms <- stats::delete.response(stats::terms(outcome_mod))
  fixed_design <- function(level) {
    fixed <- data
    fixed[[fit@exposure]] <- level
    stats::model.matrix(terms, stats::model.frame(terms, fixed))
  }
  x0 <- fixed_design(levels[[1]])
  x1 <- fixed_design(levels[[2]])
  eta0 <- as.numeric(x0 %*% beta)
  eta1 <- as.numeric(x1 %*% beta)
  pred0 <- family$linkinv(eta0)
  pred1 <- family$linkinv(eta1)

  indicator <- if (is.null(fit@focal_level)) {
    rep(1, n)
  } else {
    as.numeric(as.character(data[[fit@exposure]]) == fit@focal_level)
  }
  tilt <- indicator * s
  mu0 <- sum(tilt * pred0) / sum(tilt)
  mu1 <- sum(tilt * pred1) / sum(tilt)

  psi <- ee@psi
  p <- ncol(psi)
  stacked <- cbind(psi, score, tilt * (pred0 - mu0), tilt * (pred1 - mu1))
  meat <- crossprod(stacked) / n

  size <- p + q + 2L
  theta <- seq_len(p)
  beta_index <- p + seq_len(q)
  index0 <- p + q + 1L
  index1 <- p + q + 2L
  bread <- matrix(0, size, size)
  bread[theta, theta] <- ee@jacobian / n
  bread[beta_index, theta] <- crossprod(
    design * residual,
    s * weight_derivative
  ) /
    n
  bread[beta_index, beta_index] <- -crossprod(design * working_weight, design) /
    n
  bread[index0, beta_index] <- colSums(tilt * family$mu.eta(eta0) * x0) / n
  bread[index1, beta_index] <- colSums(tilt * family$mu.eta(eta1) * x1) / n
  bread[index0, index0] <- -sum(tilt) / n
  bread[index1, index1] <- -sum(tilt) / n

  bread_inv <- solve(bread)
  cov <- bread_inv %*% meat %*% t(bread_inv) / n
  contrast <- numeric(size)
  contrast[index1] <- 1
  contrast[index0] <- -1
  sqrt(as.numeric(t(contrast) %*% cov %*% contrast))
}

# A nonparametric bootstrap of the adjusted-model risk difference. Each
# replicate refits the balancing weights, refits the adjusted outcome model at
# those weights, and recomputes the tilt-standardized contrast, which is the
# same shape and the same replicate count as the marginal-model bootstrap above
# so the two cost about the same.
adjusted_boot_rd_se <- function(
  data,
  method,
  estimand,
  focal,
  formula,
  replicates = 200
) {
  n <- nrow(data)
  draws <- withr::with_seed(2024, {
    vapply(
      seq_len(replicates),
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
              method = eval(method),
              estimand = estimand,
              focal_level = focal
            )
            boot_w <- as.numeric(stats::weights(boot_fit))
            boot_mod <- fit_outcome(
              formula,
              resampled,
              boot_w,
              stats::binomial()
            )
            means <- adjusted_means(boot_mod, resampled, focal_level = focal)
            means$mu1 - means$mu0
          }),
          error = function(e) NA_real_
        )
      },
      numeric(1)
    )
  })
  stats::sd(draws, na.rm = TRUE)
}

# The whole estimating-equation family accepts an adjusted model, across
# estimands and across both outcome families. Each case pins three things at
# once: that the model is accepted at all, that its marginal means are the
# tilt-standardized g-computation rather than a plain average of the
# predictions, and that every reported standard error is a real number a reader
# could quote. The focal case is the one that separates the two
# standardizations, since its target population is a strict subset of the units.

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
  )
)) {
  local({
    spec <- spec
    test_that(
      paste0(
        "ipw() standardizes a covariate-adjusted outcome model for ",
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
        binary_mod <- fit_outcome(
          y ~ exposure + x1 + x2,
          data,
          w,
          stats::binomial()
        )
        continuous_mod <- fit_outcome(
          y_cont ~ exposure + x1,
          data,
          w,
          stats::gaussian()
        )

        binary <- ipw(fit, binary_mod)
        continuous <- ipw(fit, continuous_mod)
        binary_estimates <- as.data.frame(binary)
        continuous_estimates <- as.data.frame(continuous)

        binary_means <- adjusted_means(
          binary_mod,
          data,
          focal_level = spec$focal
        )
        continuous_means <- adjusted_means(
          continuous_mod,
          data,
          focal_level = spec$focal
        )

        expect_identical(binary_estimates$effect, c("rd", "log(rr)", "log(or)"))
        expect_identical(continuous_estimates$effect, "diff")

        expect_equal(
          binary$fit$theta[["mu0"]],
          binary_means$mu0,
          tolerance = 1e-8
        )
        expect_equal(
          binary$fit$theta[["mu1"]],
          binary_means$mu1,
          tolerance = 1e-8
        )
        expect_equal(
          continuous$fit$theta[["mu0"]],
          continuous_means$mu0,
          tolerance = 1e-8
        )
        expect_equal(
          continuous$fit$theta[["mu1"]],
          continuous_means$mu1,
          tolerance = 1e-8
        )

        expect_equal(
          binary_estimates$estimate,
          c(
            binary_means$mu1 - binary_means$mu0,
            log(binary_means$mu1) - log(binary_means$mu0),
            stats::qlogis(binary_means$mu1) - stats::qlogis(binary_means$mu0)
          ),
          tolerance = 1e-8
        )
        expect_equal(
          continuous_estimates$estimate,
          continuous_means$mu1 - continuous_means$mu0,
          tolerance = 1e-8
        )

        # A focal estimand standardizes over a strict subset of the units, so
        # its means are not the plain average of the predictions. Pinning that
        # the two readings differ here is what makes the assertions above a
        # check on the standardization rather than on the predictions alone.
        if (!is.null(spec$focal)) {
          pooled_means <- adjusted_means(binary_mod, data)
          expect_false(isTRUE(all.equal(pooled_means$mu0, binary_means$mu0)))
          expect_false(isTRUE(all.equal(pooled_means$mu1, binary_means$mu1)))
        }

        expect_true(all(is.finite(binary_estimates$std.err)))
        expect_true(all(binary_estimates$std.err > 0))
        expect_true(all(is.finite(continuous_estimates$std.err)))
        expect_true(all(continuous_estimates$std.err > 0))
      }
    )
  })
}

# The average effect on the untreated is the case where the target population
# could be read off the wrong group. The fit resolves the focal level to the
# group it holds fixed, which is the untreated one here, so the standardization
# population is that group rather than its complement. Both readings are
# available and they differ, so the assertions separate them: the means are
# pinned against the untreated-standardized oracle and against not being the
# treated-standardized one, and the standard error is pinned against the oracle,
# which reads the focal level from the fit and would move with a swap.

test_that("ipw() standardizes an adjusted model over the untreated for an atc fit", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "atc"
  )
  expect_identical(fit@focal_level, "0")

  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure + x1 + x2, data, w, stats::binomial())
  untreated <- adjusted_means(outcome_mod, data, focal_level = "0")
  treated <- adjusted_means(outcome_mod, data, focal_level = "1")
  oracle_se <- adjusted_rd_se(fit, outcome_mod, data)

  result <- ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)

  expect_equal(result$fit$theta[["mu0"]], untreated$mu0, tolerance = 1e-8)
  expect_equal(result$fit$theta[["mu1"]], untreated$mu1, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(treated$mu0, untreated$mu0)))
  expect_false(isTRUE(all.equal(treated$mu1, untreated$mu1)))
  expect_equal(
    estimates$estimate[estimates$effect == "rd"],
    untreated$mu1 - untreated$mu0,
    tolerance = 1e-8
  )
  expect_equal(
    estimates$std.err[estimates$effect == "rd"],
    oracle_se,
    tolerance = 1e-8
  )
})

# Sampling weights enter the standardization as well as the balancing weights: a
# unit standing for more of the target population contributes more of the
# marginal mean. The two readings are far enough apart on this fixture to tell
# apart, so pinning the sampling-weighted one rules out a standardization that
# averages the predictions over the units alone.

test_that("ipw() standardizes an adjusted model over the sampling weights", {
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
  outcome_mod <- fit_outcome(y ~ exposure + x1 + x2, data, w, stats::binomial())

  result <- ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)
  means <- adjusted_means(outcome_mod, data, sampling = data$sw)
  unweighted <- adjusted_means(outcome_mod, data)

  expect_false(isTRUE(all.equal(means$mu0, unweighted$mu0)))
  expect_equal(result$fit$theta[["mu0"]], means$mu0, tolerance = 1e-8)
  expect_equal(result$fit$theta[["mu1"]], means$mu1, tolerance = 1e-8)
  expect_true(all(is.finite(estimates$std.err)))
  expect_true(all(estimates$std.err > 0))
})

# The whole family is compared against the independent oracle, across estimands
# and with sampling weights, at a tolerance the finite-differenced bread clears
# by two orders of magnitude. The tolerance is what gives these cases their
# second job. Entropy reports its weights at the scale its container stores, so
# the reporting factor is one and the group sums do not move with the
# parameters; tilting and the covariate balancing propensity score renormalize
# by a real per-group factor, so for those two the shortcut of holding the
# reporting scale fixed misses the oracle by three to four orders of magnitude
# more than the tolerance allows, while carrying the renormalization through
# meets it.

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
        "the adjusted-model risk-difference standard error matches the oracle for ",
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
        outcome_mod <- fit_outcome(
          y ~ exposure + x1 + x2,
          data,
          w,
          stats::binomial()
        )
        oracle_se <- adjusted_rd_se(fit, outcome_mod, data)

        estimates <- as.data.frame(ipw(fit, outcome_mod))
        rd_se <- estimates$std.err[estimates$effect == "rd"]

        expect_equal(rd_se, oracle_se, tolerance = 1e-8)
      }
    )
  })
}

test_that("the adjusted-model standard error is coherent with sampling weights", {
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
  outcome_mod <- fit_outcome(y ~ exposure + x1 + x2, data, w, stats::binomial())
  oracle_se <- adjusted_rd_se(fit, outcome_mod, data, sampling = data$sw)

  estimates <- as.data.frame(ipw(fit, outcome_mod))
  rd_se <- estimates$std.err[estimates$effect == "rd"]

  expect_equal(rd_se, oracle_se, tolerance = 1e-8)
})

# The bootstrap is the external check that the adjusted stack is calibrated:
# every source of uncertainty the standard error claims to carry is one the
# resampling actually shows. Both a pooled and a focal estimand are pinned,
# since the focal case standardizes over a subset whose composition varies from
# replicate to replicate and the pooled case does not.

test_that("ipw() adjusted-model standard errors track a bootstrap for entropy ate", {
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
  outcome_mod <- fit_outcome(y ~ exposure + x1 + x2, data, w, stats::binomial())

  estimates <- as.data.frame(ipw(fit, outcome_mod))
  rd_se <- estimates$std.err[estimates$effect == "rd"]
  boot_se <- adjusted_boot_rd_se(
    data,
    quote(bw_entropy()),
    "ate",
    NULL,
    y ~ exposure + x1 + x2
  )

  # The bootstrap is noisy at this replicate count, so the agreement is loose.
  expect_equal(rd_se, boot_se, tolerance = 0.15)
})

test_that("ipw() adjusted-model standard errors track a bootstrap for bw_ipt att", {
  skip_on_cran()
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "att",
    focal_level = "1"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure + x1 + x2, data, w, stats::binomial())

  estimates <- as.data.frame(ipw(fit, outcome_mod))
  rd_se <- estimates$std.err[estimates$effect == "rd"]
  boot_se <- adjusted_boot_rd_se(
    data,
    quote(bw_ipt()),
    "att",
    "1",
    y ~ exposure + x1 + x2
  )

  expect_equal(rd_se, boot_se, tolerance = 0.15)
})

# The bootstrap at a larger sample size, where the sandwich and the resampling
# should agree closely and a disagreement is therefore easier to attribute. The
# two candidate weight couplings move this standard error by under a tenth of a
# percent, far less than a bootstrap of any feasible size resolves, so the
# choice between them is settled against the oracle above and what this case
# pins is the calibration of the adjusted stack as a whole.

test_that("ipw() adjusted-model standard errors track a bootstrap for bw_ipt ate", {
  skip_on_cran()
  data <- ipw_fixture(600)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure + x1 + x2, data, w, stats::binomial())

  estimates <- as.data.frame(ipw(fit, outcome_mod))
  rd_se <- estimates$std.err[estimates$effect == "rd"]
  boot_se <- adjusted_boot_rd_se(
    data,
    quote(bw_ipt()),
    "ate",
    NULL,
    y ~ exposure + x1 + x2
  )

  expect_equal(rd_se, boot_se, tolerance = 0.15)
})

# Accepting adjusted models must not disturb the marginal ones. The
# standardization gains a target-population indicator and the weight map gains a
# renormalization, and both reduce to what the marginal path already did: the
# indicator scales a marginal model's mean rows by one constant per row, which
# leaves the sandwich alone because those rows are zero at the solution, and the
# renormalization cancels in the coupling. Pinning the marginal results against
# the frozen coherent oracle rather than against a recorded number is what makes
# this a check on the estimator instead of a check on the last run.

test_that("supporting adjusted outcome models leaves the marginal ones alone", {
  data <- ipw_fixture()
  pooled <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  focal <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "att",
    focal_level = "1"
  )

  for (fit in list(pooled, focal)) {
    w <- as.numeric(stats::weights(fit))
    outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())
    means <- marginal_means(outcome_mod, data)

    estimates <- as.data.frame(ipw(fit, outcome_mod))
    rd <- estimates$estimate[estimates$effect == "rd"]
    rd_se <- estimates$std.err[estimates$effect == "rd"]

    expect_equal(rd, means$mu1 - means$mu0, tolerance = 1e-8)
    expect_equal(rd_se, coherent_rd_se(fit, data), tolerance = 1e-8)
  }
})

# The exposure is what the marginal means are computed by fixing, so a model
# that does not carry it describes no contrast at all: its two fixed-exposure
# designs are the same design, and the effect table would report a risk
# difference of zero as though it were an estimate. Adjusting for covariates is
# what is now allowed, not dropping the exposure, so the two are pinned
# together.

test_that("ipw() requires the exposure among the outcome model's predictors", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  adjusted_mod <- fit_outcome(
    y ~ exposure + x1 + x2,
    data,
    w,
    stats::binomial()
  )
  covariates_only <- fit_outcome(y ~ x1 + x2, data, w, stats::binomial())

  expect_s3_class(ipw(fit, adjusted_mod), "ipw")
  expect_error(
    ipw(fit, covariates_only),
    class = "balancing_ipw_input_error"
  )

  # The refusal has to say which model would be accepted, since the neighbouring
  # mistake is a model that adjusts for the covariates and forgets the exposure.
  cnd <- rlang::catch_cnd(
    ipw(fit, covariates_only),
    classes = "balancing_ipw_input_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
})

# A caller who writes the exposure into an offset has written it into the
# formula, so being told it appears in none of the model's terms reads as a
# contradiction of what they can see. The refusal adds a line naming the offset
# as the reason. Both spellings reach that line, since the terms object records
# only the first, and such a model would be refused by the offset check anyway
# once a real exposure term were added, so the line saves a round trip as well
# as the confusion.

test_that("ipw() names the offset when it carries the only mention of the exposure", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  data$.wts <- as.numeric(stats::weights(fit))
  formula_offset <- suppressWarnings(stats::glm(
    y ~ x1 + x2 + offset(exposure),
    data = data,
    family = stats::binomial(),
    weights = .wts
  ))
  argument_offset <- suppressWarnings(stats::glm(
    y ~ x1 + x2,
    data = data,
    family = stats::binomial(),
    weights = .wts,
    offset = exposure
  ))

  for (outcome_mod in list(formula_offset, argument_offset)) {
    expect_error(
      ipw(fit, outcome_mod),
      "An offset is not a term",
      class = "balancing_ipw_input_error"
    )
  }
  cnd <- rlang::catch_cnd(
    ipw(fit, formula_offset),
    classes = "balancing_ipw_input_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
})

# An interaction between the exposure and a covariate needs no special handling.
# The fixed-exposure design is rebuilt from the model's own terms with the
# exposure column set to one level, so `model.matrix()` recomputes the
# interaction columns from that level, exactly as propensity's stacked system
# builds its counterfactual designs. The effect is then a genuinely
# heterogeneous one averaged over the target population, which is the case where
# the choice of population matters most, so the focal fit is pinned beside the
# pooled one.

test_that("ipw() supports an interaction between the exposure and a covariate", {
  data <- ipw_fixture()
  pooled <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  focal <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "att",
    focal_level = "1"
  )

  for (spec in list(
    list(fit = pooled, focal_level = NULL),
    list(fit = focal, focal_level = "1")
  )) {
    w <- as.numeric(stats::weights(spec$fit))
    outcome_mod <- fit_outcome(y ~ exposure * x1, data, w, stats::binomial())
    means <- adjusted_means(outcome_mod, data, focal_level = spec$focal_level)

    result <- ipw(spec$fit, outcome_mod)
    estimates <- as.data.frame(result)

    expect_equal(result$fit$theta[["mu0"]], means$mu0, tolerance = 1e-8)
    expect_equal(result$fit$theta[["mu1"]], means$mu1, tolerance = 1e-8)
    expect_equal(
      estimates$estimate[estimates$effect == "rd"],
      means$mu1 - means$mu0,
      tolerance = 1e-8
    )
    expect_true(all(is.finite(estimates$std.err)))
    expect_true(all(estimates$std.err > 0))
  }
})

# ---- Categorical exposures -------------------------------------------------

# A categorical exposure with K observed levels reports K marginal means and one
# block of contrasts per non-reference level, each measured against the reference
# level, which is the first level in the fit's own ordering. That is propensity's
# categorical contract, and balancing follows it so that the two packages report
# a categorical effect the same way: the stacked parameter vector names the means
# `mu_<level>` and each contrast `<effect>_<level>`, and the estimates table
# gains a `comparison` column, placed after `effect`, naming the contrast as
# `"<level> vs <reference>"`. The table therefore has one row per effect measure
# per non-reference level: (K - 1) * 3 rows for a binomial outcome and K - 1 for
# a gaussian one.
#
# The outcome model carries the exposure as a factor predictor. Nothing else
# about it changes: it may adjust for covariates, and the marginal means are
# standardized over the estimand's target population exactly as they are for a
# binary exposure, over every unit for a pooled estimand and over the focal
# group for a focal one.

test_that("ipw() computes effects for a categorical bw_ipt ate fit", {
  data <- ipw_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  result <- ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)

  expect_s3_class(result, "ipw")
  expect_identical(result$se_method, "mestimation")
  expect_identical(
    estimates$effect,
    rep(c("rd", "log(rr)", "log(or)"), times = 2)
  )
  expect_identical(
    estimates$comparison,
    rep(c("b vs a", "c vs a"), each = 3)
  )
})

test_that("the categorical estimates table keeps the shared column contract", {
  data <- ipw_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  result <- ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)

  # The comparison column sits between `effect` and `estimate`, and every other
  # column of the binary contract is present and populated as it is there.
  expect_named(
    estimates,
    c(
      "effect",
      "comparison",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
  expect_identical(nrow(estimates), 6L)
  expect_true(all(is.finite(estimates$estimate)))
  expect_true(all(is.finite(estimates$std.err)))
  expect_true(all(estimates$ci.lower < estimates$estimate))
  expect_true(all(estimates$ci.upper > estimates$estimate))
  expect_true(all(estimates$conf.level == 0.95))
  expect_true(all(estimates$p.value >= 0 & estimates$p.value <= 1))
  expect_equal(
    estimates$z,
    estimates$estimate / estimates$std.err,
    tolerance = 1e-12
  )
})

test_that("the categorical variance system names K means and K - 1 contrasts", {
  data <- ipw_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  result <- ipw(fit, outcome_mod)

  p <- length(estimating_equations(fit)@parameters)
  expected_names <- c(
    paste0("theta_w", seq_len(p)),
    paste0("beta_", colnames(stats::model.matrix(outcome_mod))),
    "mu_a",
    "mu_b",
    "mu_c",
    "rd_b",
    "log(rr)_b",
    "log(or)_b",
    "rd_c",
    "log(rr)_c",
    "log(or)_c"
  )
  expect_named(result$fit$theta, expected_names)
  expect_identical(
    dimnames(result$fit$vcov),
    list(expected_names, expected_names)
  )

  # Each reported row is read off the stack at the contrast's own name, so the
  # estimates table and the variance system cannot drift apart.
  estimates <- as.data.frame(result)
  compared_level <- sub(" vs .*$", "", estimates$comparison)
  keys <- paste0(estimates$effect, "_", compared_level)
  expect_equal(
    estimates$estimate,
    unname(result$fit$theta[keys]),
    tolerance = 1e-12
  )
  expect_equal(
    estimates$std.err,
    unname(sqrt(diag(result$fit$vcov))[keys]),
    tolerance = 1e-12
  )
})

test_that("a categorical continuous outcome reports one difference per level", {
  data <- ipw_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y_cont ~ exposure, data, w, stats::gaussian())

  result <- ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)

  expect_identical(estimates$effect, c("diff", "diff"))
  expect_identical(estimates$comparison, c("b vs a", "c vs a"))
  expect_identical(
    utils::tail(names(result$fit$theta), 5L),
    c("mu_a", "mu_b", "mu_c", "diff_b", "diff_c")
  )
})

test_that("a categorical bw_cbps fit works after its over-identified request is ignored", {
  # The over-identified criterion belongs to a binary exposure, so a categorical
  # specification that asks for it is warned and fits the just-identified form.
  # That fit is an ordinary categorical covariate balancing fit, so the stacked
  # variance is available on it exactly as it is for the plain specification.
  data <- ipw_categorical_fixture()
  fit <- suppressWarnings(balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(over_identified = TRUE),
    estimand = "ate"
  ))
  reference <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cbps(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  result <- ipw(fit, outcome_mod)
  expected <- ipw(
    reference,
    fit_outcome(
      y ~ exposure,
      data,
      as.numeric(stats::weights(reference)),
      stats::binomial()
    )
  )

  expect_s3_class(result, "ipw")
  expect_identical(result$se_method, "mestimation")
  expect_equal(result$fit$theta, expected$fit$theta, tolerance = 1e-12)
  expect_equal(result$fit$vcov, expected$fit$vcov, tolerance = 1e-12)
})

test_that("a categorical ipw() result prints its comparisons", {
  data <- ipw_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  result <- ipw(fit, outcome_mod)
  output <- utils::capture.output(print(result))

  expect_true(any(grepl("b vs a", output, fixed = TRUE)))
  expect_true(any(grepl("c vs a", output, fixed = TRUE)))
})

# The point estimates are the weighted g-computation means, one per level, and
# the contrasts are those means combined by the same three formulas a binary
# exposure uses. For a marginal model the means reduce further, to the weighted
# group means, so both readings are pinned.

test_that("categorical marginal means are the weighted g-computation means", {
  data <- ipw_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())
  means <- categorical_marginal_means(outcome_mod, data)

  result <- ipw(fit, outcome_mod)

  expect_equal(
    unname(result$fit$theta[c("mu_a", "mu_b", "mu_c")]),
    unname(means),
    tolerance = 1e-8
  )

  # A model whose only predictor is the exposure is saturated, so each marginal
  # mean is that level's weighted outcome mean.
  group_means <- vapply(
    levels(data$exposure),
    function(level) {
      idx <- data$exposure == level
      stats::weighted.mean(data$y[idx], w[idx])
    },
    numeric(1)
  )
  expect_equal(unname(means), unname(group_means), tolerance = 1e-8)
})

test_that("categorical contrasts follow the formulas against the reference", {
  data <- ipw_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())
  means <- categorical_marginal_means(outcome_mod, data)

  result <- ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)
  estimate_of <- function(effect, comparison) {
    estimates$estimate[
      estimates$effect == effect & estimates$comparison == comparison
    ]
  }

  for (level in c("b", "c")) {
    comparison <- paste0(level, " vs a")
    expect_equal(
      estimate_of("rd", comparison),
      means[[level]] - means[["a"]],
      tolerance = 1e-8
    )
    expect_equal(
      estimate_of("log(rr)", comparison),
      log(means[[level]]) - log(means[["a"]]),
      tolerance = 1e-8
    )
    expect_equal(
      estimate_of("log(or)", comparison),
      stats::qlogis(means[[level]]) - stats::qlogis(means[["a"]]),
      tolerance = 1e-8
    )
  }
})

test_that("the categorical reference level follows the fit's level ordering", {
  data <- ipw_categorical_fixture()
  # The exposure's own level order, not the alphabetical one, decides which
  # level the contrasts are measured against, since that is the order
  # `balance()` groups the data in and the order the fit's weights are reported
  # per group. Declaring the levels in reverse makes the two orders disagree, so
  # an implementation that sorted the levels itself would report the contrasts
  # against the wrong level while every other assertion still held.
  data$exposure <- factor(data$exposure, levels = c("c", "b", "a"))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())
  means <- categorical_marginal_means(outcome_mod, data)

  result <- ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)

  expect_identical(
    estimates$comparison,
    rep(c("b vs c", "a vs c"), each = 3)
  )
  expect_identical(
    utils::tail(names(result$fit$theta), 9L),
    c(
      "mu_c",
      "mu_b",
      "mu_a",
      "rd_b",
      "log(rr)_b",
      "log(or)_b",
      "rd_a",
      "log(rr)_a",
      "log(or)_a"
    )
  )
  expect_equal(
    estimates$estimate[estimates$effect == "rd"],
    c(means[["b"]] - means[["c"]], means[["a"]] - means[["c"]]),
    tolerance = 1e-8
  )
})

# ---- Continuous exposures --------------------------------------------------

# A continuous exposure has no levels to contrast, so there is no pair of
# marginal means to difference and no g-computation step to standardize. What
# ipw() reports instead is the dose-response coefficient of a weighted marginal
# structural model: the balancing weights break the exposure-covariate
# association, the outcome model carries exactly one term in the exposure, and
# that term's coefficient is the effect of a one-unit change in the exposure on
# the model's own link scale. The estimates table therefore holds a single row,
# named for that link, and carries neither the marginal-mean rows nor the
# contrast rows a discrete exposure produces.
#
# Only entropy balancing reaches this path, and only at exact balance and the
# average treatment effect, which is the only estimand a continuous entropy fit
# targets.

# The pieces of the continuous entropy tilt, re-derived here in plain R rather
# than read off the fit. The tilt weights a single group of every unit,
#
#   w_i(theta) = n_eff q_i exp(-c_i' theta) / sum_j s_j q_j exp(-c_j' theta),
#
# with c_i the unit's constraint row, q the base weights, s the sampling
# weights, and n_eff the total the group is normalized to. Its estimating
# functions are psi_ij = s_i w_i(theta) (c_ij - t_j), one per constraint.
#
# The constraint rows are the marginal columns the fit holds at their
# base-measure means, the standardized exposure followed by the covariate
# columns, and then the exposure-covariate product columns it drives to zero.
# The covariate columns come from the fit's own recipe, which records the
# centers and scales they were built with; everything else is assembled here.
# The helpers assume the uniform base weights every fixture in this file uses,
# which is what lets the base measure be the sampling weights alone.
continuous_tilt_pieces <- function(fit, data, sampling = NULL) {
  n <- nrow(data)
  s <- sampling %||% rep(1, n)
  base <- rep(1, n)
  measure <- s * base
  total <- sum(measure)

  exposure <- as.numeric(data[[fit@exposure]])
  covariates <- rebuild_constraint_matrix(fit@recipe, data)

  # The exposure crosses the product columns centered and scaled on the base
  # measure, so that a product column of weighted mean zero is a weighted
  # correlation of zero. The reliability denominator is the one the package's
  # own weighted spread uses, and it reduces to the unweighted standard
  # deviation when the measure is uniform.
  center <- sum(measure * exposure) / total
  spread <- sqrt(
    sum(measure * (exposure - center)^2) / (total - sum(measure^2) / total)
  )

  marginals <- cbind(
    (exposure - mean(exposure)) / stats::sd(exposure),
    covariates
  )
  products <- covariates * ((exposure - center) / spread)
  list(
    covs = unname(cbind(marginals, products)),
    targets = c(
      as.numeric(crossprod(marginals, measure / total)),
      rep(0, ncol(products))
    ),
    base = base,
    sampling = s,
    n_eff = sum(s)
  )
}

continuous_tilt_weights <- function(pieces, theta) {
  tilt <- pieces$base * exp(-as.numeric(pieces$covs %*% theta))
  pieces$n_eff * tilt / sum(pieces$sampling * tilt)
}

continuous_tilt_psi <- function(pieces, theta) {
  (pieces$sampling * continuous_tilt_weights(pieces, theta)) *
    sweep(pieces$covs, 2, pieces$targets, "-")
}

# An independent standard error for the exposure coefficient of a weighted
# marginal structural model, and the continuous counterpart of
# `coherent_rd_se()`. It is the whole stacked M-estimator written out in plain
# R: the tilt's estimating functions above, then the outcome model's weighted
# score, evaluated together at the root each fit already found. The bread is a
# central finite difference of the stacked column sums at that root and the meat
# is their empirical second moment, so the oracle shares no algebra with the
# implementation and reaches nothing the implementation reaches.
#
# The score is written for the general link, `s_i w_i (y_i - mu_i) mu'(eta_i) /
# V(mu_i) x_i`, which is the equation `glm()` solves and which reduces to
# `s_i w_i (y_i - x_i' beta) x_i` for a linear model.
continuous_msm_se <- function(fit, outcome_mod, data, sampling = NULL) {
  pieces <- continuous_tilt_pieces(fit, data, sampling)
  n <- nrow(data)
  s <- pieces$sampling
  design <- stats::model.matrix(outcome_mod)
  family <- stats::family(outcome_mod)
  response <- as.numeric(
    stats::model.response(stats::model.frame(outcome_mod))
  )

  theta <- estimating_equations(fit)@parameters
  p <- length(theta)
  q <- ncol(design)
  root <- c(theta, stats::coef(outcome_mod))

  stacked <- function(parameters) {
    tilt <- parameters[seq_len(p)]
    coefficients <- parameters[p + seq_len(q)]
    tilted <- continuous_tilt_weights(pieces, tilt)
    eta <- as.numeric(design %*% coefficients)
    mu <- family$linkinv(eta)
    score <- (s *
      tilted *
      (response - mu) *
      family$mu.eta(eta) /
      family$variance(mu)) *
      design
    cbind(continuous_tilt_psi(pieces, tilt), score)
  }

  meat <- crossprod(stacked(root)) / n
  size <- p + q
  bread <- vapply(
    seq_len(size),
    function(j) {
      step <- 1e-6 * max(1, abs(root[[j]]))
      up <- root
      down <- root
      up[[j]] <- up[[j]] + step
      down[[j]] <- down[[j]] - step
      (colSums(stacked(up)) - colSums(stacked(down))) / (2 * step * n)
    },
    numeric(size)
  )

  bread_inv <- solve(bread)
  covariance <- bread_inv %*% meat %*% t(bread_inv) / n
  position <- p + match(fit@exposure, colnames(design))
  sqrt(covariance[position, position])
}

# The standard error the same marginal structural model carries when the weights
# are treated as fixed: the empirical sandwich of its score alone, with no
# weight-parameter block above it. It is the continuous counterpart of
# `naive_rd_se()` and is used the same way, to show that accounting for having
# estimated the weights moves the standard error. Written for a linear model,
# whose residuals are the score's own.
naive_msm_se <- function(outcome_mod, wts, exposure_name) {
  design <- stats::model.matrix(outcome_mod)
  n <- nrow(design)
  residual <- as.numeric(stats::residuals(outcome_mod))
  score <- (wts * residual) * design
  bread <- -crossprod(design * wts, design) / n
  bread_inv <- solve(bread)
  covariance <- bread_inv %*% (crossprod(score) / n) %*% t(bread_inv) / n
  position <- match(exposure_name, colnames(design))
  sqrt(covariance[position, position])
}

# A weighted marginal structural model of the shape the continuous path accepts.
# A `NULL` family fits a plain linear model, which is the identity-link case;
# every other family goes through `glm()`, wrapped because balancing weights are
# not counts and a binomial fit says so at every call. The weights ride along as
# a column so the model frame resolves them, as they do for the discrete path.
fit_msm <- function(formula, data, wts, family = NULL, ...) {
  data[[".wts"]] <- wts
  if (is.null(family)) {
    return(stats::lm(formula, data = data, weights = .wts))
  }
  suppressWarnings(
    stats::glm(formula, data = data, family = family, weights = .wts, ...)
  )
}

# The oracle above is only an oracle if its reconstruction of the tilt is the
# tilt the fit solved, so that reconstruction is pinned against the container
# before anything is asked of it. Both the weights and the estimating functions
# are checked, with and without sampling weights, since the sampling weights
# enter the base measure the constraint targets are taken under as well as the
# estimating functions themselves.

test_that("the continuous entropy container matches a hand-rolled tilt", {
  data <- ipw_continuous_fixture()
  data$sw <- withr::with_seed(7, stats::runif(nrow(data), 0.5, 2))

  for (sampled in c(FALSE, TRUE)) {
    fit <- if (sampled) {
      balance(
        data,
        exposure,
        c(x1, x2),
        method = bw_entropy(),
        estimand = "ate",
        sampling_weights = sw
      )
    } else {
      balance(
        data,
        exposure,
        c(x1, x2),
        method = bw_entropy(),
        estimand = "ate"
      )
    }
    pieces <- continuous_tilt_pieces(
      fit,
      data,
      sampling = if (sampled) data$sw else NULL
    )
    theta <- estimating_equations(fit)@parameters

    expect_equal(
      continuous_tilt_weights(pieces, theta),
      as.numeric(stats::weights(fit, include_sampling_weights = FALSE)),
      tolerance = 1e-10
    )
    expect_equal(
      unname(continuous_tilt_psi(pieces, theta)),
      unname(estimating_equations(fit)@psi),
      tolerance = 1e-10
    )
  }
})

# ---- The reported effect ---------------------------------------------------

test_that("ipw() reports the exposure slope for a continuous entropy ate fit", {
  data <- ipw_continuous_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  outcome_mod <- fit_msm(
    y_cont ~ exposure,
    data,
    as.numeric(stats::weights(fit))
  )

  result <- ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)

  expect_s3_class(result, "ipw")
  expect_identical(nrow(estimates), 1L)
  expect_identical(estimates$effect, "slope")
  expect_named(
    estimates,
    c(
      "effect",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
  expect_equal(
    estimates$estimate,
    stats::coef(outcome_mod)[["exposure"]],
    tolerance = 1e-10
  )
  expect_true(is.finite(estimates$std.err))
  expect_gt(estimates$std.err, 0)
})

test_that("ipw() accepts a covariate-adjusted continuous marginal structural model", {
  data <- ipw_continuous_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  marginal_mod <- fit_msm(y_cont ~ exposure, data, w)
  adjusted_mod <- fit_msm(y_cont ~ exposure + v, data, w)

  marginal <- as.data.frame(ipw(fit, marginal_mod))
  adjusted <- as.data.frame(ipw(fit, adjusted_mod))

  expect_identical(nrow(adjusted), 1L)
  expect_identical(adjusted$effect, "slope")
  expect_equal(
    adjusted$estimate,
    stats::coef(adjusted_mod)[["exposure"]],
    tolerance = 1e-10
  )
  expect_true(is.finite(adjusted$std.err))
  expect_gt(adjusted$std.err, 0)

  # `v` is prognostic and unbalanced, so adjusting for it is a real change to
  # the model rather than one the weights have already made irrelevant. Pinning
  # that the two results differ is what makes this a check on the adjusted path.
  expect_false(isTRUE(all.equal(adjusted$estimate, marginal$estimate)))
  expect_false(isTRUE(all.equal(adjusted$std.err, marginal$std.err)))
})

# The single effect row is named for the outcome model's link, since that is the
# scale the exposure coefficient is a one-unit effect on. An identity link
# reports a slope, whether it arrives as a linear model or as a gaussian model
# with the same link; a logit reports a log odds ratio; a log link reports a log
# risk ratio.

test_that("ipw() names the continuous effect for the outcome model's link", {
  data <- ipw_continuous_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))

  models <- list(
    list(effect = "slope", model = fit_msm(y_cont ~ exposure, data, w)),
    list(
      effect = "slope",
      model = fit_msm(y_cont ~ exposure, data, w, stats::gaussian())
    ),
    list(
      effect = "log(or)",
      model = fit_msm(y ~ exposure, data, w, stats::binomial())
    ),
    list(
      effect = "log(rr)",
      model = fit_msm(
        y_rare ~ exposure,
        data,
        w,
        stats::quasibinomial(link = "log"),
        start = c(log(mean(data$y_rare)), 0)
      )
    )
  )

  for (spec in models) {
    outcome_mod <- spec$model
    # A log-link fit has to keep every fitted risk below one, so a run that
    # stopped short would leave the coefficient this case reads meaningless.
    expect_false(isFALSE(outcome_mod$converged))
    estimates <- as.data.frame(ipw(fit, outcome_mod))

    expect_identical(nrow(estimates), 1L)
    expect_identical(estimates$effect, spec$effect)
    expect_equal(
      estimates$estimate,
      stats::coef(outcome_mod)[["exposure"]],
      tolerance = 1e-10
    )
  }
})

# ---- Continuous standard errors --------------------------------------------

# Each case compares the reported standard error against the independent stacked
# oracle at a tolerance the finite-differenced bread clears by several orders of
# magnitude. The three shapes separate what the stack has to carry: the marginal
# linear model, where the design is the exposure alone; the adjusted one, where
# the score couples an unbalanced covariate into the weight block; and the
# binomial one, where the link enters the score and its derivative.

for (spec in list(
  list(
    label = "a marginal linear model",
    formula = quote(y_cont ~ exposure),
    family = NULL
  ),
  list(
    label = "an adjusted linear model",
    formula = quote(y_cont ~ exposure + v),
    family = NULL
  ),
  list(
    label = "a binomial outcome model",
    formula = quote(y ~ exposure),
    family = quote(stats::binomial())
  )
)) {
  local({
    spec <- spec
    test_that(
      paste0(
        "the continuous ipw() standard error matches the stacked oracle for ",
        spec$label
      ),
      {
        data <- ipw_continuous_fixture()
        fit <- balance(
          data,
          exposure,
          c(x1, x2),
          method = bw_entropy(),
          estimand = "ate"
        )
        outcome_mod <- fit_msm(
          eval(spec$formula),
          data,
          as.numeric(stats::weights(fit)),
          eval(spec$family)
        )
        oracle_se <- continuous_msm_se(fit, outcome_mod, data)

        estimates <- as.data.frame(ipw(fit, outcome_mod))

        expect_identical(nrow(estimates), 1L)
        expect_equal(estimates$std.err, oracle_se, tolerance = 1e-8)
      }
    )
  })
}

test_that("the continuous ipw() standard error is coherent with sampling weights", {
  data <- ipw_continuous_fixture()
  data$sw <- withr::with_seed(7, stats::runif(nrow(data), 0.5, 2))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate",
    sampling_weights = sw
  )
  outcome_mod <- fit_msm(
    y_cont ~ exposure,
    data,
    as.numeric(stats::weights(fit))
  )
  oracle_se <- continuous_msm_se(fit, outcome_mod, data, sampling = data$sw)

  estimates <- as.data.frame(ipw(fit, outcome_mod))

  expect_true(is.finite(estimates$std.err))
  expect_gt(estimates$std.err, 0)
  expect_equal(estimates$std.err, oracle_se, tolerance = 1e-8)
})

test_that("the continuous ipw() standard error differs from the weights-fixed sandwich", {
  data <- ipw_continuous_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_msm(y_cont ~ exposure, data, w)
  naive_se <- naive_msm_se(outcome_mod, w, "exposure")

  estimates <- as.data.frame(ipw(fit, outcome_mod))

  # Accounting for weight estimation moves the standard error; it does not equal
  # the sandwich a weighted regression reports when it treats the weights as a
  # design quantity.
  expect_false(isTRUE(all.equal(estimates$std.err, naive_se)))
})

test_that("ipw() reports the mestimation variance system for a continuous fit", {
  data <- ipw_continuous_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  outcome_mod <- fit_msm(
    y_cont ~ exposure,
    data,
    as.numeric(stats::weights(fit))
  )

  result <- ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)

  expect_identical(result$se_method, "mestimation")
  expect_named(result$fit, c("theta", "vcov"))
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

# ---- Continuous validation -------------------------------------------------

# A continuous entropy fit targets the average treatment effect and nothing
# else, so an estimand supplied to ipw() either agrees with the fit or is a
# specification error, on the same terms as for a discrete exposure.

test_that("ipw() takes the estimand from a continuous fit and refuses another", {
  data <- ipw_continuous_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  outcome_mod <- fit_msm(
    y_cont ~ exposure,
    data,
    as.numeric(stats::weights(fit))
  )

  expect_identical(ipw(fit, outcome_mod)$estimand, "ate")
  expect_identical(ipw(fit, outcome_mod, estimand = "ate")$estimand, "ate")
  expect_error(
    ipw(fit, outcome_mod, estimand = "att"),
    class = "balancing_estimand_error"
  )
})

# The reported effect is one coefficient, so the outcome model has to have one
# coefficient to report: exactly one design column that reads the exposure. A
# model carrying a second exposure column describes a dose-response curve rather
# than a slope, and no single coefficient of it is the effect. The three ways to
# arrive at one are pinned together, since they differ in how the second column
# gets in: a second term in the exposure, one term that expands to two columns,
# and an interaction that makes the slope depend on a covariate. Covariates
# alongside the exposure remain fine, so an accepted model is pinned beside
# them.

test_that("ipw() refuses a continuous outcome model with more than one exposure term", {
  data <- ipw_continuous_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  accepted <- fit_msm(y_cont ~ exposure + v, data, w)
  refused <- list(
    quadratic = fit_msm(y_cont ~ exposure + I(exposure^2), data, w),
    polynomial = fit_msm(y_cont ~ poly(exposure, 2), data, w),
    interaction = fit_msm(y_cont ~ exposure * x1, data, w)
  )

  expect_s3_class(ipw(fit, accepted), "ipw")
  for (outcome_mod in refused) {
    expect_error(
      ipw(fit, outcome_mod),
      class = "balancing_ipw_input_error"
    )
  }
})

# What a model missing the exposure is told it may do instead depends on the
# exposure type, and the two readings promise different things. A discrete
# exposure is read through predictions with the exposure fixed to each level, so
# such a model may carry the exposure inside a transformation. A continuous
# exposure reports the exposure's own coefficient, and the one-term contract
# above refuses every transformation, so this refusal has to ask for a term of
# its own rather than offer one that would only be turned away a second time. A
# transformation written into the formula still names the exposure and so reaches
# the one-term contract instead; the model that arrives here is one whose
# transformed exposure was computed into a column of its own beforehand.

test_that("ipw() asks a continuous outcome model for the exposure as its own term", {
  data <- ipw_continuous_fixture()
  data$exposure_scaled <- data$exposure / 10
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  transformed_mod <- fit_msm(y_cont ~ exposure_scaled + x1, data, w)

  cnd <- rlang::catch_cnd(
    ipw(fit, transformed_mod),
    classes = "balancing_ipw_input_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
})

# An offset is a second way for the exposure to reach the linear predictor, and
# it reaches it outside the term labels the exposure-term check reads. An offset
# in the exposure shifts the fitted coefficient by whatever the offset
# contributes, so the reported effect would be a number the model does not
# estimate rather than the dose-response coefficient the table names. Both entry
# points are pinned, since the terms object records only the first: an
# `offset()` term in the formula, and the model's `offset` argument.

test_that("ipw() refuses a continuous outcome model whose offset reads the exposure", {
  data <- ipw_continuous_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  data$.wts <- w
  refused <- list(
    formula_term = fit_msm(y_cont ~ exposure + offset(exposure), data, w),
    formula_transformed = fit_msm(
      y_cont ~ exposure + offset(2 * exposure),
      data,
      w
    ),
    # The `offset` argument is written into the call here rather than passed
    # through `fit_msm()`, because a model function records `..1` for an
    # argument that arrived through another function's dots, and it is the
    # argument's own expression the check reads.
    offset_argument = stats::lm(
      y_cont ~ exposure,
      data = data,
      weights = .wts,
      offset = exposure
    )
  )

  for (outcome_mod in refused) {
    expect_error(
      ipw(fit, outcome_mod),
      class = "balancing_ipw_input_error"
    )
  }
})

# An offset that never reads the exposure is a legitimate part of the model and
# stays supported. A log-link rate model is the case that motivates one: its
# person-time offset makes the exposure coefficient a log rate ratio, which is
# the effect the table already names log(rr). Both entry points are pinned
# again, since the refusal above has to separate the offsets that read the
# exposure from the offsets that do not rather than turn them all away.

test_that("ipw() accepts a continuous outcome model with an exposure-free offset", {
  data <- ipw_continuous_fixture()
  data$followup <- withr::with_seed(91, stats::runif(nrow(data), 0.5, 2))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))

  rate_mod <- fit_msm(
    y_rare ~ exposure + offset(log(followup)),
    data,
    w,
    stats::quasibinomial(link = "log"),
    start = c(log(mean(data$y_rare)), 0)
  )
  # A log-link fit has to keep every fitted risk below one, so a run that
  # stopped short would leave the coefficient this case reads meaningless.
  expect_false(isFALSE(rate_mod$converged))
  data$.wts <- w
  shifted_mod <- stats::lm(
    y_cont ~ exposure,
    data = data,
    weights = .wts,
    offset = followup
  )

  rate <- as.data.frame(ipw(fit, rate_mod))
  shifted <- as.data.frame(ipw(fit, shifted_mod))

  expect_identical(nrow(rate), 1L)
  expect_identical(rate$effect, "log(rr)")
  expect_equal(
    rate$estimate,
    stats::coef(rate_mod)[["exposure"]],
    tolerance = 1e-10
  )
  expect_true(is.finite(rate$std.err))
  expect_gt(rate$std.err, 0)

  expect_identical(nrow(shifted), 1L)
  expect_identical(shifted$effect, "slope")
  expect_equal(
    shifted$estimate,
    stats::coef(shifted_mod)[["exposure"]],
    tolerance = 1e-10
  )
  expect_true(is.finite(shifted$std.err))
  expect_gt(shifted$std.err, 0)
})

# The single effect row is named for the outcome model's link, and only three
# links have a name the coefficient answers to: the identity's slope, the
# logit's log odds ratio, and the log link's log risk ratio. A probit fit runs
# and returns a coefficient like any other, and nothing about that number says
# it is none of the three, so a table naming it would be wrong in a way its
# reader could not see. The refusal is what keeps it out, and the message has to
# name the links that are supported, since the caller's next move is to refit
# under one of them.

test_that("ipw() refuses a continuous outcome model with an unnamed link", {
  data <- ipw_continuous_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  outcome_mod <- fit_msm(
    y ~ exposure,
    data,
    as.numeric(stats::weights(fit)),
    stats::binomial(link = "probit")
  )
  # The refusal is about the link alone, so the model it comes from has to be
  # one every other check accepts.
  expect_identical(
    stats::family(outcome_mod)$link,
    "probit"
  )

  expect_error(
    ipw(fit, outcome_mod),
    class = "balancing_ipw_input_error"
  )

  cnd <- rlang::catch_cnd(
    ipw(fit, outcome_mod),
    classes = "balancing_ipw_input_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
})

# A column name R cannot parse as a symbol is written back-quoted in a formula,
# and the term label and the coefficient name carry those back-quotes while the
# fit records the bare name. The exposure-term check compares the two, so it has
# to compare them on one spelling: a single-term marginal structural model over
# such a column is an ordinary model and is accepted, and only the name it
# reports the coefficient under differs from the syntactic case.

test_that("ipw() accepts a continuous marginal structural model with a non-syntactic exposure name", {
  data <- ipw_continuous_fixture()
  renamed <- data.frame(data, check.names = FALSE)
  names(renamed)[names(renamed) == "exposure"] <- "dose level"

  reference_fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  fit <- balance(
    renamed,
    `dose level`,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  expect_identical(fit@exposure, "dose level")

  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_msm(y_cont ~ `dose level`, renamed, w)
  reference_mod <- fit_msm(
    y_cont ~ exposure,
    data,
    as.numeric(stats::weights(reference_fit))
  )

  estimates <- as.data.frame(ipw(fit, outcome_mod))
  reference <- as.data.frame(ipw(reference_fit, reference_mod))

  expect_identical(nrow(estimates), 1L)
  expect_identical(estimates$effect, "slope")
  expect_equal(
    estimates$estimate,
    stats::coef(outcome_mod)[["`dose level`"]],
    tolerance = 1e-10
  )
  expect_equal(estimates$estimate, reference$estimate, tolerance = 1e-8)
  expect_equal(estimates$std.err, reference$std.err, tolerance = 1e-8)
})

# The continuous path is exact entropy balancing alone. A continuous fit from
# another method solves no smooth estimating equations at all, and an entropy
# fit at a positive tolerance stops solving them the moment the tolerance binds,
# so both raise the shared unsupported condition pointing to the bootstrap
# workflow.

test_that("ipw() rejects continuous fits outside exact entropy balancing", {
  data <- ipw_continuous_fixture()
  energy <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_energy(),
    estimand = "ate"
  )
  # The relaxed fit cannot meet the requested tolerance on this fixture and says
  # so; the warning is the fit's own and is not what this test is about.
  inexact <- suppressWarnings(balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate",
    constraints = balance_terms(tolerance = 0.05)
  ))

  for (fit in list(energy, inexact)) {
    outcome_mod <- fit_msm(
      y_cont ~ exposure,
      data,
      as.numeric(stats::weights(fit))
    )
    expect_error(
      ipw(fit, outcome_mod),
      class = "balancing_ipw_unsupported_error"
    )
  }
})

# ---- Exposure coding -------------------------------------------------------

# A categorical exposure need not be a factor. A character column and an integer
# code describe the same levels in the same order the fit resolves them in, so
# each has to carry the whole chain: the fit, the outcome model, the
# counterfactual designs, and the effects. The factor-coded fit is the reference
# every coding is compared against, since only the labels differ.

test_that("ipw() computes effects for a character categorical exposure", {
  data <- ipw_categorical_fixture()
  data$arm <- as.character(data$exposure)

  reference_fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  fit <- balance(data, arm, c(x1, x2), method = bw_ipt(), estimand = "ate")
  expect_identical(fit@exposure_type, "categorical")
  expect_identical(fit@exposure_levels, c("a", "b", "c"))

  w <- as.numeric(stats::weights(fit))
  expect_equal(w, as.numeric(stats::weights(reference_fit)), tolerance = 1e-8)

  outcome_mod <- fit_outcome(y ~ arm, data, w, stats::binomial())
  reference_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  estimates <- as.data.frame(ipw(fit, outcome_mod))
  reference <- as.data.frame(ipw(reference_fit, reference_mod))

  expect_identical(estimates$comparison, reference$comparison)
  expect_equal(estimates$estimate, reference$estimate, tolerance = 1e-8)
  expect_equal(estimates$std.err, reference$std.err, tolerance = 1e-8)
})

test_that("ipw() computes effects for an integer-coded categorical exposure", {
  # An integer code carries no levels of its own, so the saturated outcome model
  # writes the factor in the formula. The model frame then stores the transformed
  # column rather than the exposure, which is what `.data` supplies.
  data <- ipw_categorical_fixture()
  data$arm <- match(as.character(data$exposure), c("a", "b", "c"))

  reference_fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  fit <- balance(data, arm, c(x1, x2), method = bw_ipt(), estimand = "ate")
  expect_identical(fit@exposure_type, "categorical")
  expect_identical(fit@exposure_levels, c("1", "2", "3"))

  w <- as.numeric(stats::weights(fit))
  expect_equal(w, as.numeric(stats::weights(reference_fit)), tolerance = 1e-8)

  outcome_mod <- fit_outcome(y ~ factor(arm), data, w, stats::binomial())
  reference_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  estimates <- as.data.frame(ipw(fit, outcome_mod, .data = data))
  reference <- as.data.frame(ipw(reference_fit, reference_mod))

  expect_identical(
    estimates$comparison,
    c(rep("2 vs 1", 3L), rep("3 vs 1", 3L))
  )
  expect_equal(estimates$estimate, reference$estimate, tolerance = 1e-8)
  expect_equal(estimates$std.err, reference$std.err, tolerance = 1e-8)
})

test_that("a transformed exposure still needs the exposure column", {
  # The model frame of a transformed exposure stores `factor(arm)`, not `arm`,
  # so the counterfactual designs cannot be built from it alone. The refusal
  # names the argument that carries the untransformed column.
  data <- ipw_categorical_fixture()
  data$arm <- match(as.character(data$exposure), c("a", "b", "c"))
  fit <- balance(data, arm, c(x1, x2), method = bw_ipt(), estimand = "ate")
  outcome_mod <- fit_outcome(
    y ~ factor(arm),
    data,
    as.numeric(stats::weights(fit)),
    stats::binomial()
  )

  expect_error(
    ipw(fit, outcome_mod),
    class = "balancing_ipw_input_error",
    regexp = ".data",
    fixed = TRUE
  )
})

test_that("a categorical adjusted model standardizes over everyone for ate", {
  data <- ipw_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure + x1, data, w, stats::binomial())
  means <- categorical_marginal_means(outcome_mod, data)

  result <- ipw(fit, outcome_mod)

  expect_equal(
    unname(result$fit$theta[c("mu_a", "mu_b", "mu_c")]),
    unname(means),
    tolerance = 1e-8
  )
})

test_that("a categorical att standardizes an adjusted model over the focal group", {
  data <- ipw_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "att",
    focal_level = "b"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure + x1, data, w, stats::binomial())
  focal <- categorical_marginal_means(
    outcome_mod,
    data,
    tilt = as.numeric(data$exposure == "b")
  )
  pooled <- categorical_marginal_means(outcome_mod, data)

  result <- ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)

  expect_equal(
    unname(result$fit$theta[c("mu_a", "mu_b", "mu_c")]),
    unname(focal),
    tolerance = 1e-8
  )

  # The focal group is not the whole sample, so standardizing over it has to
  # move the means. Without that the test would pass on an implementation that
  # ignored the estimand entirely.
  expect_false(isTRUE(all.equal(unname(focal), unname(pooled))))
  expect_true(all(is.finite(estimates$std.err)))
  expect_true(all(estimates$std.err > 0))
})

# The standard errors are checked three ways: they are finite and positive
# throughout, they match an independent analytic oracle exactly, and they track a
# nonparametric bootstrap. The analytic oracle is the categorical extension of
# `coherent_rd_se()`, which generalizes to K levels without new algebra: the
# stack gains one marginal-mean row per level and the contrast is read off the
# joint covariance at the two levels it compares.

test_that("categorical standard errors are finite and positive", {
  data <- ipw_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))

  for (formula in list(y ~ exposure, y ~ exposure + x1, y ~ exposure * x1)) {
    outcome_mod <- fit_outcome(formula, data, w, stats::binomial())
    estimates <- as.data.frame(ipw(fit, outcome_mod))
    expect_true(all(is.finite(estimates$std.err)))
    expect_true(all(estimates$std.err > 0))
  }
})

for (spec in list(
  list(label = "a categorical bw_ipt ate fit", method = quote(bw_ipt())),
  list(
    label = "a categorical just-identified bw_cbps fit",
    method = quote(bw_cbps())
  ),
  list(
    label = "a categorical bw_entropy ate fit",
    method = quote(bw_entropy())
  )
)) {
  local({
    spec <- spec
    test_that(
      paste0(
        "categorical risk-difference standard errors match the coherent oracle for ",
        spec$label
      ),
      {
        data <- ipw_categorical_fixture()
        fit <- balance(
          data,
          exposure,
          c(x1, x2),
          method = eval(spec$method),
          estimand = "ate"
        )
        w <- as.numeric(stats::weights(fit))
        outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

        result <- ipw(fit, outcome_mod)
        estimates <- as.data.frame(result)

        for (level in c("b", "c")) {
          reported <- estimates$std.err[
            estimates$effect == "rd" &
              estimates$comparison == paste0(level, " vs a")
          ]
          expect_equal(
            reported,
            coherent_categorical_rd_se(fit, data, level),
            tolerance = 1e-8
          )
        }
      }
    )
  })
}

test_that("the categorical standard error is coherent with sampling weights", {
  data <- ipw_categorical_fixture()
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

  result <- ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)
  rd_se <- estimates$std.err[
    estimates$effect == "rd" & estimates$comparison == "b vs a"
  ]

  expect_equal(
    rd_se,
    coherent_categorical_rd_se(fit, data, "b", sampling = data$sw),
    tolerance = 1e-8
  )
})

test_that("categorical standard errors track a nonparametric bootstrap", {
  skip_on_cran()
  data <- ipw_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  result <- ipw(fit, outcome_mod)
  estimates <- as.data.frame(result)
  rd_se <- estimates$std.err[
    estimates$effect == "rd" & estimates$comparison == "b vs a"
  ]

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
              method = bw_ipt(),
              estimand = "ate"
            )
            boot_w <- as.numeric(stats::weights(boot_fit))
            boot_mod <- fit_outcome(
              y ~ exposure,
              resampled,
              boot_w,
              stats::binomial()
            )
            means <- categorical_marginal_means(boot_mod, resampled)
            means[["b"]] - means[["a"]]
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

# The lift serves every method whose categorical fit carries a container with
# both re-evaluation hooks, which is the whole estimating-equation family:
# inverse probability tilting, the just-identified covariate balancing propensity
# score, which reuses the tilt's entrypoints for a categorical exposure, and
# entropy balancing at exact balance. Each is pinned on the same shape and point
# estimates so that none of them can drift out of the supported set unnoticed.

for (spec in list(
  list(
    label = "a just-identified bw_cbps categorical ate fit",
    method = quote(bw_cbps()),
    estimand = "ate",
    focal = NULL
  ),
  list(
    label = "a just-identified bw_cbps categorical att fit",
    method = quote(bw_cbps()),
    estimand = "att",
    focal = "b"
  ),
  list(
    label = "a bw_entropy categorical ate fit",
    method = quote(bw_entropy()),
    estimand = "ate",
    focal = NULL
  ),
  list(
    label = "a bw_entropy categorical att fit",
    method = quote(bw_entropy()),
    estimand = "att",
    focal = "b"
  )
)) {
  local({
    spec <- spec
    test_that(paste0("ipw() computes effects for ", spec$label), {
      data <- ipw_categorical_fixture()
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
      tilt <- if (is.null(spec$focal)) {
        NULL
      } else {
        as.numeric(data$exposure == spec$focal)
      }
      means <- categorical_marginal_means(outcome_mod, data, tilt = tilt)

      result <- ipw(fit, outcome_mod)
      estimates <- as.data.frame(result)

      expect_identical(
        estimates$effect,
        rep(c("rd", "log(rr)", "log(or)"), times = 2)
      )
      expect_identical(
        estimates$comparison,
        rep(c("b vs a", "c vs a"), each = 3)
      )
      expect_equal(
        unname(result$fit$theta[c("mu_a", "mu_b", "mu_c")]),
        unname(means),
        tolerance = 1e-8
      )
      expect_true(all(is.finite(estimates$std.err)))
      expect_true(all(estimates$std.err > 0))
    })
  })
}

# Two guards on the level set itself. The contrasts are labeled by the fit's
# levels and the weights were solved per group, so an outcome model fitted on
# data that no longer carries every level describes a different exposure than the
# one the fit weighted: its counterfactual designs would be built from a level
# set the weights never saw, and the result would be an ordinary-looking effect
# table for the wrong contrast. The estimand check is the one a binary fit
# already makes, pinned here because the categorical path reaches it too.

test_that("ipw() rejects an outcome model whose data drop an exposure level", {
  data <- ipw_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  collapsed <- data
  collapsed$exposure <- factor(
    ifelse(
      as.character(data$exposure) == "c",
      "b",
      as.character(data$exposure)
    ),
    levels = c("a", "b")
  )
  outcome_mod <- fit_outcome(y ~ exposure, collapsed, w, stats::binomial())

  expect_error(
    ipw(fit, outcome_mod),
    class = "balancing_ipw_input_error"
  )
})

test_that("ipw() rejects an estimand that contradicts a categorical fit", {
  data <- ipw_categorical_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  expect_error(
    ipw(fit, outcome_mod, estimand = "att"),
    class = "balancing_estimand_error"
  )
})

test_that("a binary fit's estimates table carries no comparison column", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  estimates <- as.data.frame(ipw(fit, outcome_mod))

  # A binary exposure has one comparison, so naming it would add a column that
  # says the same thing on every row. The eight-column contract is what
  # propensity reports there, and lifting the categorical case must not disturb
  # it.
  expect_named(
    estimates,
    c(
      "effect",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
  expect_identical(estimates$effect, c("rd", "log(rr)", "log(or)"))
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

  wide <- as.data.frame(ipw(fit, outcome_mod, conf_level = 0.95))
  narrow <- as.data.frame(ipw(fit, outcome_mod, conf_level = 0.80))

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
    ipw(fit, outcome_mod, estimand = "att"),
    class = "balancing_estimand_error"
  )

  cnd <- rlang::catch_cnd(
    ipw(fit, outcome_mod, estimand = "att"),
    classes = "balancing_estimand_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
})

test_that("ipw() rejects an estimand outside the vocabulary", {
  # A name no estimand carries is a vocabulary error rather than a disagreement
  # with the fit, and the two failures need different messages: a misspelling
  # reported as a mismatch names the fit's estimand and leaves the caller to
  # notice that theirs is not an estimand at all. The vocabulary is the one
  # `balance()` matches against, so a spelling that creates a fit is a spelling
  # that names it here.
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
    ipw(fit, outcome_mod, estimand = "bogus"),
    regexp = "must be one of"
  )
  expect_error(
    ipw(fit, outcome_mod, estimand = "ATE"),
    regexp = "must be one of"
  )

  # A name the vocabulary carries but the fit does not target is still the
  # mismatch it always was.
  expect_error(
    ipw(fit, outcome_mod, estimand = "ato"),
    class = "balancing_estimand_error"
  )

  cnd <- rlang::catch_cnd(ipw(fit, outcome_mod, estimand = "bogus"))
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
})

# The fit stores the untreated target as "atu" however it was spelled, so a
# request spelled "atc" names the same estimand the fit already targets and is
# not a contradiction. The three spellings of that one request, the synonym, the
# canonical name, and letting the fit supply it, therefore have to return the
# same result down to the estimand the result reports.

test_that("ipw() accepts the atc synonym for the estimand the fit stores", {
  data <- ipw_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "atc"
  )
  expect_identical(fit@estimand, "atu")

  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_outcome(y ~ exposure, data, w, stats::binomial())

  synonym <- expect_no_error(ipw(fit, outcome_mod, estimand = "atc"))
  canonical <- ipw(fit, outcome_mod, estimand = "atu")
  inherited <- ipw(fit, outcome_mod)

  expect_identical(synonym, canonical)
  expect_identical(synonym, inherited)
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
    ipw(fit, list(coefficients = 1)),
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
    ipw(fit, outcome_mod, .data = one_level),
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
    ipw(fit, outcome_mod),
    class = "balancing_ipw_input_error",
    regexp = "weights(fit)",
    fixed = TRUE
  )

  cnd <- rlang::catch_cnd(
    ipw(fit, outcome_mod),
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
    ipw(fit, outcome_mod),
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
    ipw(fit, noise_mod),
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
    ipw(fit, other_mod),
    class = "balancing_ipw_input_error"
  )
})

test_that("the weight-mismatch message points at the composed weights", {
  # A fit with sampling weights reports them composed into the balancing weights,
  # and the outcome model has to be fitted with that composed vector. Fitting it
  # with either factor alone is the mistake this message is most often read
  # after, so it names the composition rather than only the accessor.
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
  sampling_mod <- fit_outcome(y ~ exposure, data, data$sw, stats::binomial())

  expect_error(
    ipw(fit, sampling_mod),
    class = "balancing_ipw_input_error",
    regexp = "sampling weights",
    fixed = TRUE
  )

  cnd <- rlang::catch_cnd(
    ipw(fit, sampling_mod),
    classes = "balancing_ipw_input_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
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

  binary <- as.data.frame(ipw(fit, binary_mod))
  continuous <- as.data.frame(ipw(fit, continuous_mod))

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

  estimates <- as.data.frame(ipw(fit, outcome_mod))
  means <- marginal_means(outcome_mod, data)

  expect_equal(
    estimates$estimate[estimates$effect == "rd"],
    means$mu1 - means$mu0,
    tolerance = 1e-10
  )
})

test_that("ipw() honors an offset in a categorical outcome model", {
  # The categorical path assembles one fixed-exposure design per level, so the
  # offset has to survive K designs rather than two. Both spellings are run
  # against the same g-computation oracle, since they reach the fitted model by
  # different routes.
  data <- ipw_categorical_fixture()
  data$log_time <- withr::with_seed(11, stats::rnorm(nrow(data), 0, 0.3))
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  data$.wts <- w
  models <- list(
    formula_term = fit_outcome(
      y ~ exposure + offset(log_time),
      data,
      w,
      stats::binomial()
    ),
    offset_argument = suppressWarnings(stats::glm(
      y ~ exposure,
      data = data,
      family = stats::binomial(),
      weights = .wts,
      offset = log_time
    ))
  )

  # An offset leaves the model unsaturated in the exposure, so the marginal
  # means stop being the weighted group means and start depending on the offset
  # each unit carries. A fixed-exposure design that dropped the offset would
  # predict from the offset-fitted coefficients with no offset at all; the
  # offset-free refit below stands in for that reading to show the gap between it
  # and the offset means is far larger than the tolerance the equalities below
  # are held to.
  offset_free <- fit_outcome(y ~ exposure, data, w, stats::binomial())
  expect_gt(
    max(abs(
      categorical_marginal_means(models$formula_term, data) -
        categorical_marginal_means(offset_free, data)
    )),
    0.005
  )

  for (outcome_mod in models) {
    means <- categorical_marginal_means(outcome_mod, data)
    result <- ipw(fit, outcome_mod)
    estimates <- as.data.frame(result)

    expect_equal(
      unname(result$fit$theta[c("mu_a", "mu_b", "mu_c")]),
      unname(means),
      tolerance = 1e-10
    )
    for (level in c("b", "c")) {
      expect_equal(
        estimates$estimate[
          estimates$effect == "rd" &
            estimates$comparison == paste0(level, " vs a")
        ],
        means[[level]] - means[["a"]],
        tolerance = 1e-10
      )
    }
  }
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

  estimates <- as.data.frame(ipw(fit, outcome_mod))
  engine <- ipw_deli_sandwich(
    container = estimating_equations(fit),
    outcome_mod = outcome_mod,
    frame = data,
    exposure_name = fit@exposure,
    levels = fit@exposure_levels,
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

# An offset that reads the exposure is a different object from the offsets above
# and is refused. The fixed-exposure designs move the exposure to each level
# while the offset is held at its observed value, which is what an offset is: a
# known per-unit quantity the counterfactual does not move. An offset computed
# from the exposure is not known that way, and the predictions such a model
# yields are neither factual nor counterfactual, since each one reads one level
# in the design and another in the offset. The contrast of those means is what
# the table would report, and on this fixture it comes back with the opposite
# sign to the g-computation the same model implies.
#
# Both entry points are pinned on each path, since the terms object records only
# the first, and the categorical path is covered beside the discrete one because
# its per-level designs are assembled the same way.

test_that("ipw() refuses a discrete outcome model whose offset reads the exposure", {
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
  refused <- list(
    formula_term = fit_outcome(
      y ~ exposure + offset(exposure),
      data,
      w,
      stats::quasibinomial()
    ),
    formula_transformed = fit_outcome(
      y ~ exposure + offset(2 * exposure),
      data,
      w,
      stats::quasibinomial()
    ),
    # The `offset` argument is written into the call here rather than passed
    # through `fit_outcome()`, because a model function records `..1` for an
    # argument that arrived through another function's dots, and it is the
    # argument's own expression the check reads.
    offset_argument = suppressWarnings(stats::glm(
      y ~ exposure,
      data = data,
      family = stats::quasibinomial(),
      weights = .wts,
      offset = exposure
    ))
  )

  for (outcome_mod in refused) {
    expect_error(
      ipw(fit, outcome_mod),
      class = "balancing_ipw_input_error"
    )
  }
})

test_that("ipw() refuses a categorical outcome model whose offset reads the exposure", {
  data <- sim_categorical(n = 300)
  data$y <- withr::with_seed(
    5,
    stats::rbinom(nrow(data), 1L, 0.4 + 0.1 * (data$exposure == "b"))
  )
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_ipt(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  data$.wts <- w
  refused <- list(
    formula_term = fit_outcome(
      y ~ exposure + offset(as.numeric(exposure == "b")),
      data,
      w,
      stats::quasibinomial()
    ),
    # As in the discrete case above, the `offset` argument is written into the
    # call here rather than passed through `fit_outcome()`, because a model
    # function records `..1` for an argument that arrived through another
    # function's dots, and it is the argument's own expression the check reads.
    offset_argument = suppressWarnings(stats::glm(
      y ~ exposure,
      data = data,
      family = stats::quasibinomial(),
      weights = .wts,
      offset = as.numeric(exposure == "b")
    ))
  )

  for (outcome_mod in refused) {
    expect_error(
      ipw(fit, outcome_mod),
      class = "balancing_ipw_input_error"
    )
  }
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
    ipw(fit, outcome_mod),
    class = "balancing_ipw_input_error"
  )

  cnd <- rlang::catch_cnd(
    ipw(fit, outcome_mod),
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
    ipw(fit, quasipoisson_mod),
    class = "balancing_ipw_input_error"
  )
  expect_error(
    ipw(fit, inverse_mod),
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

  binomial_result <- as.data.frame(ipw(fit, binomial_mod))
  quasi_result <- as.data.frame(ipw(fit, quasi_mod))
  gaussian_result <- as.data.frame(ipw(fit, gaussian_mod))
  lm_result <- as.data.frame(ipw(fit, lm_mod))

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

  numeric_result <- as.data.frame(ipw(fit, numeric_mod))
  factor_result <- as.data.frame(ipw(fit, factor_mod))

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
    ipw(fit, factor_mod),
    class = "balancing_ipw_input_error"
  )

  cnd <- rlang::catch_cnd(
    ipw(fit, factor_mod),
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

  stored_result <- as.data.frame(ipw(fit, stored))
  dropped_result <- as.data.frame(ipw(fit, dropped))

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

  glm_result <- as.data.frame(ipw(fit, glm_mod))
  lm_result <- as.data.frame(ipw(fit, lm_mod))

  expect_identical(lm_result$effect, "diff")
  expect_identical(lm_result$effect, glm_result$effect)
  expect_equal(lm_result$estimate, glm_result$estimate)
  expect_equal(lm_result$std.err, glm_result$std.err)
})

# A two-column response through `glm()` is the grouped binomial form, one row per
# group of trials rather than one row per Bernoulli draw. It is refused, and the
# refusal has to name the shape: `glm()` records the prior weights as the weights
# it was given times each row's trial count, so the weight preflight sees a
# mismatch and would otherwise tell a caller who supplied exactly the fit's
# weights to supply them. The same check turns away a multivariate `lm()`, whose
# response is a matrix for a different reason, so the message names both.

test_that("ipw() refuses a grouped binomial outcome model", {
  data <- ipw_fixture()
  trials <- 4L
  data$successes <- withr::with_seed(
    505,
    stats::rbinom(
      nrow(data),
      trials,
      stats::plogis(-0.3 + 0.5 * data$exposure + 0.4 * data$x1)
    )
  )
  data$failures <- trials - data$successes
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  grouped_mod <- fit_outcome(
    cbind(successes, failures) ~ exposure,
    data,
    w,
    stats::binomial()
  )

  # The weights the model recorded are the supplied weights times the trial
  # count, which is the mismatch the preflight used to report.
  expect_equal(as.numeric(stats::weights(grouped_mod)), w * trials)

  expect_error(
    ipw(fit, grouped_mod),
    class = "balancing_ipw_unsupported_error"
  )

  cnd <- rlang::catch_cnd(
    ipw(fit, grouped_mod),
    classes = "balancing_ipw_unsupported_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
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
    ipw(fit, outcome_mod),
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
    ipw(fit, outcome_mod),
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
    ipw(fit, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  )
})

# The variance engine reaches the weight path only through the container's two
# re-evaluation hooks, so a container carrying neither cannot be differentiated
# at all. Every fit that reaches ipw() today populates both, so the guard is
# defence in depth. What it buys is a comprehensible refusal in place of a
# failure raised from inside the engine against a NULL it was handed.

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
    ipw(no_psi, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  )
  expect_error(
    ipw(no_weights, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  )

  cnd <- rlang::catch_cnd(
    ipw(without_hooks(NULL, NULL), outcome_mod),
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
    ipw(fit, outcome_mod),
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
#     levels,
#     categorical = FALSE,
#     sampling_weights = NULL
#   )
#
# where `container` is the fit's `balancing_estimating_equations`, `frame` is
# the data frame holding the exposure, `levels` are the exposure levels the fit
# weighted in the fit's own order, `categorical` decides how the mean and
# contrast blocks are named, and `sampling_weights` is the fit's sampling weight
# vector or `NULL`. The weight parameter count comes from
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
#   mu0, mu1                the marginal means of a binary exposure
#   mu_<level>              the marginal means of a categorical exposure
#   rd, log(rr), log(or)    the contrasts for a non-gaussian outcome model
#   diff                    the single contrast for a gaussian outcome model
#
# A categorical exposure suffixes each contrast name with the level it compares
# against the reference level, as `rd_b`. The contrast names key the estimates
# table `ipw()` reports, so it reads each estimate and standard error straight
# off `theta` and the diagonal of `vcov`.

# Call the engine with the pieces read off a fit, the way the `ipw()` method
# does, including the level ordering the fit recorded.
call_deli_sandwich <- function(fit, outcome_mod, data) {
  ipw_deli_sandwich(
    container = estimating_equations(fit),
    outcome_mod = outcome_mod,
    frame = data,
    exposure_name = fit@exposure,
    levels = fit@exposure_levels,
    categorical = identical(fit@exposure_type, "categorical"),
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
      levels = fit@exposure_levels,
      sampling_weights = fit@sampling_weights
    ),
    regexp = "singular"
  )
})

# A container carrying one extra weight parameter whose own estimating equation
# has slope `epsilon`. Analytically the Jacobian is rank deficient at any rank
# tolerance worth the name, while a central difference of the same system returns
# a pivot of that size rather than a zero, and `solve()` accepts it: this is the
# shape that walks past `allow_pinv = FALSE` and answers with finite standard
# errors. A deficiency whose bread row is an exact zero, which the container
# above builds, is the easy case deli already refuses.
#
# `extra` is the extra equation's per-unit value at the fitted parameters. It
# sums to zero there, so the equation is solved, and it gives the meat a row of
# its own, which is what a wrong answer needs: the inverse of an `epsilon` pivot
# multiplied by an exactly zero meat row is still zero. `shift` decides whether
# the extra parameter moves the reported weights, and that is the whole
# difference between a deficiency the reported effects cannot see and one they
# inherit. Measured on this fixture, the moving version returns risk-difference
# standard errors a factor of 1e12 above the truth, finite and confident.
near_singular_container <- function(
  ee,
  epsilon = 1e-12,
  shift = NULL,
  extra = 0
) {
  p <- length(ee@parameters)
  jacobian <- matrix(0, p + 1L, p + 1L)
  jacobian[seq_len(p), seq_len(p)] <- ee@jacobian
  jacobian[p + 1L, p + 1L] <- epsilon
  balancing_estimating_equations(
    parameters = c(ee@parameters, 0),
    psi = cbind(ee@psi, extra),
    jacobian = jacobian,
    weight_jacobian = cbind(
      ee@weight_jacobian,
      ee@weights_raw * (shift %||% 0)
    ),
    weights_raw = ee@weights_raw,
    psi_fn = function(theta) {
      cbind(ee@psi_fn(theta[seq_len(p)]), extra + epsilon * theta[[p + 1L]])
    },
    weights_fn = if (is.null(shift)) {
      function(theta) ee@weights_fn(theta[seq_len(p)])
    } else {
      function(theta) {
        ee@weights_fn(theta[seq_len(p)]) * (1 + theta[[p + 1L]] * shift)
      }
    }
  )
}

test_that("the stacked variance tolerates a deficiency the weights are flat along", {
  # Everything downstream of the weight block reads the weight parameters only
  # through the reported weight map, so a direction the map is flat along cannot
  # reach the score, the means, or the contrasts. It contaminates the weight
  # block of the covariance, which nothing reported is read from, and the effects
  # are the effects of the undoctored stack down to the last bit.
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
  keys <- c("rd", "log(rr)", "log(or)")

  plain <- call_deli_sandwich(fit, outcome_mod, data)
  tolerated <- ipw_deli_sandwich(
    container = near_singular_container(
      estimating_equations(fit),
      extra = data$x1 - mean(data$x1)
    ),
    outcome_mod = outcome_mod,
    frame = data,
    exposure_name = fit@exposure,
    levels = fit@exposure_levels,
    sampling_weights = fit@sampling_weights
  )

  expect_true(all(is.finite(sqrt(diag(tolerated$vcov)[keys]))))
  expect_equal(
    unname(sqrt(diag(tolerated$vcov)[keys])),
    unname(sqrt(diag(plain$vcov)[keys]))
  )
})

test_that("the stacked variance refuses a deficiency that moves the weights", {
  # The same deficiency, with a weight map that is not flat along it. The
  # effects then inherit the inverse of a pivot that is rounding error, so the
  # standard errors come back finite and wrong by whatever factor that pivot
  # happens to take, which is the failure the check exists to refuse.
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
  moving <- near_singular_container(
    estimating_equations(fit),
    shift = 0.1 * as.numeric(scale(data$x1)),
    extra = data$x1 - mean(data$x1)
  )
  refuse <- function() {
    ipw_deli_sandwich(
      container = moving,
      outcome_mod = outcome_mod,
      frame = data,
      exposure_name = fit@exposure,
      levels = fit@exposure_levels,
      sampling_weights = fit@sampling_weights
    )
  }

  expect_error(refuse(), class = "balancing_ipw_unsupported_error")

  cnd <- rlang::catch_cnd(
    refuse(),
    classes = "balancing_ipw_unsupported_error"
  )
  expect_snapshot(error = TRUE, cnd_class = TRUE, stop(cnd))
})

test_that("a movement that is not a number counts as movement", {
  # The rank check reads the weight map through a function it is handed, so a
  # map with no derivative along a deficient direction leaves that direction's
  # movement missing rather than large. Missing is not flat, and a missing
  # comparison must not be left to decide the branch: `if (NA)` stops with base
  # R's own message, which names neither the fit nor the deficiency, and one
  # missing direction would take the classification of every other direction
  # with it.
  jacobian <- diag(c(1, 0.5, 0))
  parameters <- c(0, 0, 0)
  flat <- function(theta) rep(1, 5)
  undefined <- function(theta) {
    if (identical(as.numeric(theta), parameters)) {
      return(flat(theta))
    }
    c(NaN, rep(1, 4))
  }

  expect_error(
    validate_stacked_bread(jacobian, undefined, parameters),
    class = "balancing_ipw_unsupported_error"
  )

  # The same deficiency with a map that is flat along it is still tolerated, so
  # counting the missing movement has not turned every deficiency into a
  # refusal.
  expect_null(validate_stacked_bread(jacobian, flat, parameters))
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
      levels = fit@exposure_levels,
      sampling_weights = fit@sampling_weights
    )),
    class = "balancing_ipw_unsupported_error"
  )
})

# Both refusals above are raised from inside the variance engine, which no caller
# ever writes. Reached the way a caller reaches them, through `ipw()`, they must
# report the call that was made rather than the internal frame the failure
# happened in: `ipw_deli_sandwich()` and `stacked_covariance()` are not functions
# the caller can go and look at, and naming them in the error header sends a
# reader after code that is not theirs. Every other classed refusal `ipw()`
# raises already reports `ipw()`, so these two are the exceptions.
#
# The pin is on the calling frame alone. The rest of the condition is covered by
# the engine-level specs above, and the call is read through `call_name()` so an
# argument the header happens to carry does not decide the test.
test_that("the variance engine's refusals report ipw() as the failing call", {
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

  with_container <- function(container) {
    doctored <- fit
    doctored@estimating_equations <- container
    doctored
  }

  # A rank deficiency the reported weights move along, refused by the bread
  # check the engine makes before it differences anything.
  deficient <- with_container(near_singular_container(
    ee,
    shift = 0.1 * as.numeric(scale(data$x1)),
    extra = data$x1 - mean(data$x1)
  ))
  expect_error(
    ipw(deficient, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  )
  deficient_cnd <- rlang::catch_cnd(
    ipw(deficient, outcome_mod),
    classes = "balancing_ipw_unsupported_error"
  )
  expect_identical(rlang::call_name(conditionCall(deficient_cnd)), "ipw")

  # A stack whose bread is not finite, refused after the difference is taken.
  # deli warns on its own account there, which is suppressed rather than pinned.
  not_finite <- with_container(na_psi_container(ee))
  expect_error(
    suppressWarnings(ipw(not_finite, outcome_mod)),
    class = "balancing_ipw_unsupported_error"
  )
  not_finite_cnd <- rlang::catch_cnd(
    suppressWarnings(ipw(not_finite, outcome_mod)),
    classes = "balancing_ipw_unsupported_error"
  )
  expect_identical(rlang::call_name(conditionCall(not_finite_cnd)), "ipw")
})

# A continuous exposure assembles its stack in a third function, which the two
# refusals above are reached from as well. It carries the same defect, so it is
# pinned on its own rather than left to the discrete path's coverage. One case is
# enough here: the covariance helper is shared with the discrete path and already
# pinned above, so the rank check is the only frame this route adds.
test_that("the continuous variance engine reports ipw() as the failing call", {
  data <- ipw_continuous_fixture()
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_entropy(),
    estimand = "ate"
  )
  w <- as.numeric(stats::weights(fit))
  outcome_mod <- fit_msm(y_cont ~ exposure, data, w)

  deficient <- fit
  deficient@estimating_equations <- near_singular_container(
    estimating_equations(fit),
    shift = 0.1 * as.numeric(scale(data$x1)),
    extra = data$x1 - mean(data$x1)
  )

  expect_error(
    ipw(deficient, outcome_mod),
    class = "balancing_ipw_unsupported_error"
  )
  cnd <- rlang::catch_cnd(
    ipw(deficient, outcome_mod),
    classes = "balancing_ipw_unsupported_error"
  )
  expect_identical(rlang::call_name(conditionCall(cnd)), "ipw")
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
    levels = fit@exposure_levels,
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
