# An independent parity oracle for characteristic function distance balancing.
#
# The energy kernel has an external reference through bw_energy(), but the other
# kernels have no external implementation, so the objective and constraint
# assembly the Rust core solves is checked here against a quadratic program built
# and solved independently in R. The oracle rebuilds the kernel matrix in base R
# (verified against the internal kernel_matrix entry point), assembles the same
# quadratic term, linear term, group-sum equalities, and minimum-weight box that
# the design specifies, and solves it with a different backend: quadprog's
# Goldfarb-Idnani active-set method for the strictly convex kernels, and osqp for
# the indefinite energy kernel. Because the reported weights renormalize each
# group to its estimand total, the comparisons are made at the solver scale: the
# solver objective (fit@objective) against the oracle objective, the recovered
# solver weights against the oracle weights where the ridge makes the solution
# unique, and constraint satisfaction of both solutions so the oracle cannot pass
# on an infeasible point.

# ---- Oracle construction --------------------------------------------------

# Scaled-Euclidean standardization: each covariate column divided by its
# reliability-weighted standard deviation, matching the core transform. The
# denominator 1 - sum(wn^2) with normalized weights reduces to the n - 1 sample
# variance under uniform weights.
scaled_covariates <- function(covs, s) {
  wn <- s / sum(s)
  denom <- 1 - sum(wn^2)
  for (j in seq_len(ncol(covs))) {
    m <- sum(wn * covs[, j])
    v <- (sum(wn * covs[, j]^2) - m^2) / denom
    sdj <- sqrt(v)
    if (sdj > 0) covs[, j] <- covs[, j] / sdj
  }
  covs
}

# The kernel matrix, built independently of the Rust core. The gaussian kernel is
# the squared exponential at the median-heuristic bandwidth; the energy kernel is
# the negative pairwise distance.
oracle_kernel <- function(kernel, covs, s) {
  z <- scaled_covariates(covs, s)
  d <- as.matrix(stats::dist(z))
  dimnames(d) <- NULL
  if (kernel == "energy") {
    return(-d)
  }
  bw <- stats::median(d[lower.tri(d)])
  exp(-0.5 * (d / bw)^2)
}

# Sampling weights normalized to mean one within each exposure level, matching the
# per-group scaling the quadratic-program objectives assume.
group_normalized_weights <- function(s, levels, n_levels) {
  for (g in seq_len(n_levels) - 1L) {
    idx <- which(levels == g)
    m <- mean(s[idx])
    if (m > 0) s[idx] <- s[idx] / m
  }
  s
}

# Assemble the quadratic program the core solves: the group-normalized kernel
# quadratic term with the weight-penalty ridge on its diagonal, the negated
# source-weighted cross-similarity linear term, the group-sum equality rows, and
# the minimum-weight lower bounds. Mirrors methods/cfd.rs and qp_balance.rs.
build_cfd_qp <- function(
  kernel_mat,
  s,
  levels,
  n_levels,
  estimand,
  improved,
  focal,
  lambda,
  min_weight
) {
  n <- length(levels)
  s_norm <- group_normalized_weights(s, levels, n_levels)
  n_t <- tabulate(levels + 1L, nbins = n_levels)
  swnt <- ifelse(n_t[levels + 1L] > 0, s_norm / n_t[levels + 1L], 0)

  if (estimand == "ate") {
    active <- which(levels >= 0)
    group_levels <- seq_len(n_levels) - 1L
  } else {
    active <- which(levels >= 0 & levels != focal)
    group_levels <- setdiff(seq_len(n_levels) - 1L, focal)
  }
  nvar <- length(active)
  active_levels <- levels[active]
  scale <- swnt[active]

  # The group-normalization interaction: the plain average treatment effect and
  # the focal estimand pair same-level units with weight one; the improved variant
  # raises the same-level weight to the number of levels and sets the cross-level
  # weight to minus one.
  same <- outer(active_levels, active_levels, `==`)
  interaction <- if (estimand == "ate" && improved) {
    ifelse(same, n_levels, -1)
  } else {
    ifelse(same, 1, 0)
  }
  quad <- kernel_mat[active, active] * outer(scale, scale) * interaction
  diag(quad) <- diag(quad) + lambda * scale^2 / 2

  # The cross similarity runs over the whole sample for the average treatment
  # effect and over the focal group for a focal estimand.
  if (estimand == "ate") {
    src <- s_norm
    mult <- 2 / n
  } else {
    src <- ifelse(levels == focal, s_norm / max(n_t[focal + 1L], 1), 0)
    mult <- 2
  }
  cross <- as.numeric(-mult * (src %*% kernel_mat[, active]))
  linear <- cross * scale

  eq <- matrix(0, length(group_levels), nvar)
  for (r in seq_along(group_levels)) {
    eq[r, ] <- ifelse(active_levels == group_levels[r], scale, 0)
  }

  list(
    quad = quad,
    linear = linear,
    eq = eq,
    active = active,
    swnt = swnt,
    nvar = nvar,
    group_levels = group_levels,
    min_weight = min_weight
  )
}

# The solver objective x' Q x + c' x, the value the core reports on the fit. The
# spec doubles Q so that its half-quadratic convention leaves this expression.
qp_objective <- function(x, qp) {
  as.numeric(t(x) %*% qp$quad %*% x + qp$linear %*% x)
}

# Recover the core's solver weights from the reported weights. The reported
# weights are a positive per-group multiple of the solver weights, so dividing by
# the group-sum coefficient restores the solver scale on which the objective and
# the equality rows are defined.
recover_solver_weights <- function(reported, qp, levels) {
  x <- numeric(qp$nvar)
  for (t in qp$group_levels) {
    sel <- which(levels[qp$active] == t)
    idx <- qp$active[sel]
    x[sel] <- reported[idx] / sum(qp$swnt[idx] * reported[idx])
  }
  x
}

# The doubled quadratic term is positive definite exactly when the strictly
# convex active-set solver applies; the indefinite energy assembly is not.
qp_is_convex <- function(qp) {
  tryCatch(
    {
      chol(2 * qp$quad)
      TRUE
    },
    error = function(e) FALSE
  )
}

# Solve the quadratic program with a backend independent of the core: quadprog's
# Goldfarb-Idnani method when the program is strictly convex, osqp otherwise.
solve_oracle_qp <- function(qp) {
  meq <- nrow(qp$eq)
  lower <- rep(qp$min_weight, qp$nvar)
  if (qp_is_convex(qp)) {
    a_mat <- t(rbind(qp$eq, diag(qp$nvar)))
    b_vec <- c(rep(1, meq), lower)
    solution <- quadprog::solve.QP(
      2 * qp$quad,
      -qp$linear,
      a_mat,
      b_vec,
      meq = meq
    )
    list(x = solution$solution, solver = "quadprog")
  } else {
    a_mat <- rbind(qp$eq, diag(qp$nvar))
    lo <- c(rep(1, meq), lower)
    up <- c(rep(1, meq), rep(Inf, qp$nvar))
    settings <- osqp::osqpSettings(
      polishing = TRUE,
      verbose = FALSE,
      eps_abs = 1e-9,
      eps_rel = 1e-9,
      max_iter = 40000L
    )
    model <- osqp::osqp(2 * qp$quad, qp$linear, a_mat, lo, up, settings)
    solution <- if (inherits(model, "S7_object")) {
      model@Solve()
    } else {
      model$Solve()
    }
    list(x = solution$x, solver = "osqp")
  }
}

# ---- Comparison -----------------------------------------------------------

# Fit with the core and against the independent oracle, then assert objective
# parity, constraint satisfaction of both solutions, and, when the ridge makes the
# program strictly convex, weight parity as well.
expect_cfd_matches_oracle <- function(
  data,
  kernel,
  estimand,
  oracle_focal,
  n_levels,
  levels_vec,
  focal_level = NULL
) {
  penalty <- 1e-2
  min_weight <- 1e-8
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bw_cfd(kernel = kernel, weight_penalty = penalty),
    estimand = estimand,
    .focal_level = focal_level
  )
  expect_true(fit@converged)

  covs <- cbind(data$x1, data$x2)
  s <- rep(1, nrow(data))

  # The independently built kernel reproduces the internal kernel_matrix entry
  # point, so the oracle rests on the same kernel the core solves.
  kernel_rust <- kernel_matrix(
    as.numeric(covs),
    kernel,
    1,
    1.5,
    numeric(0),
    s,
    logical(0),
    list(threads = 1L)
  )
  kernel_base <- oracle_kernel(kernel, covs, s)
  expect_equal(kernel_base, kernel_rust, tolerance = 1e-9)

  qp <- build_cfd_qp(
    kernel_base,
    s,
    levels_vec,
    n_levels,
    estimand,
    improved = TRUE,
    focal = oracle_focal,
    lambda = penalty,
    min_weight = min_weight
  )
  oracle <- solve_oracle_qp(qp)
  solver_weights <- recover_solver_weights(
    as.numeric(stats::weights(fit)),
    qp,
    levels_vec
  )

  # The solver objective at the core's solution agrees with the oracle objective.
  expect_equal(fit@objective, qp_objective(oracle$x, qp), tolerance = 1e-6)
  # The reconstruction reproduces the reported objective from the recovered
  # solver weights, confirming the assembly rather than the solver alone.
  expect_equal(
    fit@objective,
    qp_objective(solver_weights, qp),
    tolerance = 1e-6
  )

  # Both solutions satisfy the group-sum equalities and the minimum-weight box, so
  # neither objective is read at an infeasible point.
  expect_lt(max(abs(qp$eq %*% oracle$x - 1)), 1e-6)
  expect_lt(max(abs(qp$eq %*% solver_weights - 1)), 1e-6)
  expect_column_all(oracle, "x", function(value) value >= min_weight - 1e-8)
  expect_true(all(solver_weights >= min_weight - 1e-8))

  # A strictly convex program has a unique minimizer, so the weight vectors agree,
  # not only the objective. The indefinite energy assembly is compared on the
  # objective and the constraints alone.
  if (identical(oracle$solver, "quadprog")) {
    expect_equal(oracle$x, solver_weights, tolerance = 1e-5)
  }
}

# ---- Problem instances ----------------------------------------------------

# Small, well-conditioned confounded designs with both groups amply populated so
# the group-sum constraints and the kernel objective are non-degenerate.
oracle_binary <- function(n = 90, seed = 7) {
  withr::with_seed(seed, {
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    exposure <- stats::rbinom(n, 1, stats::plogis(0.7 * x1 - 0.5 * x2))
    data.frame(exposure = exposure, x1 = x1, x2 = x2)
  })
}

oracle_categorical <- function(n = 96, seed = 9) {
  withr::with_seed(seed, {
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    eta_b <- 0.6 * x1 - 0.3 * x2
    eta_c <- -0.4 * x1 + 0.5 * x2
    denom <- 1 + exp(eta_b) + exp(eta_c)
    draw <- stats::runif(n)
    level <- ifelse(
      draw < 1 / denom,
      "a",
      ifelse(draw < (1 + exp(eta_b)) / denom, "b", "c")
    )
    data.frame(exposure = factor(level, c("a", "b", "c")), x1 = x1, x2 = x2)
  })
}

# ---- Binary ---------------------------------------------------------------

test_that("binary ate gaussian matches an independent quadratic program", {
  skip_if_not_installed("quadprog")
  skip_if_not_installed("osqp")
  data <- oracle_binary()
  expect_cfd_matches_oracle(
    data,
    "gaussian",
    "ate",
    oracle_focal = NA,
    n_levels = 2,
    levels_vec = data$exposure
  )
})

test_that("binary att gaussian matches an independent quadratic program", {
  skip_if_not_installed("quadprog")
  skip_if_not_installed("osqp")
  data <- oracle_binary()
  expect_cfd_matches_oracle(
    data,
    "gaussian",
    "att",
    oracle_focal = 1,
    n_levels = 2,
    levels_vec = data$exposure
  )
})

test_that("binary ate energy matches an independent quadratic program", {
  skip_if_not_installed("quadprog")
  skip_if_not_installed("osqp")
  data <- oracle_binary()
  expect_cfd_matches_oracle(
    data,
    "energy",
    "ate",
    oracle_focal = NA,
    n_levels = 2,
    levels_vec = data$exposure
  )
})

test_that("binary att energy matches an independent quadratic program", {
  skip_if_not_installed("quadprog")
  skip_if_not_installed("osqp")
  data <- oracle_binary()
  expect_cfd_matches_oracle(
    data,
    "energy",
    "att",
    oracle_focal = 1,
    n_levels = 2,
    levels_vec = data$exposure
  )
})

# ---- Categorical ----------------------------------------------------------

test_that("categorical ate gaussian matches an independent quadratic program", {
  skip_if_not_installed("quadprog")
  skip_if_not_installed("osqp")
  data <- oracle_categorical()
  levels_vec <- as.integer(data$exposure) - 1L
  expect_cfd_matches_oracle(
    data,
    "gaussian",
    "ate",
    oracle_focal = NA,
    n_levels = 3,
    levels_vec = levels_vec
  )
})

test_that("categorical att gaussian matches an independent quadratic program", {
  skip_if_not_installed("quadprog")
  skip_if_not_installed("osqp")
  data <- oracle_categorical()
  levels_vec <- as.integer(data$exposure) - 1L
  expect_cfd_matches_oracle(
    data,
    "gaussian",
    "att",
    oracle_focal = 1,
    n_levels = 3,
    levels_vec = levels_vec,
    focal_level = "b"
  )
})
