# Characteristic function distance balancing (kernel balancing): weights that
# minimize a kernel measure of the distance between the reweighted exposure
# groups and the target sample. Each kernel is positive semidefinite by
# construction, so the objective is a convex quadratic program with a
# simplex-type constraint set and the method spec carries only tuning parameters.
# fit_method() assembles the covariate columns and, for the t kernel, the
# Monte Carlo frequency projections drawn under R's random number generator,
# calls the Rust solver, renormalizes each group to its estimand target total,
# and reports the solver's dual variables for diagnostics. The energy kernel is
# the negative pairwise distance, so it reproduces energy balancing on the same
# data and constraints. The quadratic-program family has no estimating equations.

#' Characteristic function distance balancing
#'
#' `bal_cfd()` specifies characteristic function distance balancing, also
#' called kernel balancing, for [balance()]. The weights minimize a kernel
#' measure of the distance between the reweighted exposure groups and a target
#' sample, following Wong and Chan. Each supported kernel is positive
#' semidefinite, so the objective is a convex quadratic program, and the
#' reweighting improves multivariate covariate balance without positing a
#' propensity model. Characteristic function distance balancing supports binary
#' and categorical exposures.
#'
#' @details
#' The objective is the sum of each group's kernel distance to the target
#' sample, built from a covariate kernel matrix rather than a pairwise distance
#' matrix. `kernel` selects the kernel. The `"gaussian"`, `"laplace"`, and
#' `"matern"` kernels are distance based: their bandwidth is the median of the
#' pairwise covariate distances, so the fit is invariant to a uniform rescaling
#' of the covariates. `smoothness` sets the Matern order, one of `0.5`, `1.5`,
#' or `2.5`. The `"t"` kernel approximates a heavy-tailed spectral kernel by
#' Monte Carlo: `simulation_draws` frequency vectors are drawn from a
#' multivariate t distribution with `degrees_of_freedom` degrees of freedom, on
#' the R side under R's random number generator, so a fixed seed reproduces the
#' weights. The improved variant for the average treatment effect adds the
#' between-group term of the kernel mean embedding, balancing the groups against
#' one another as well as against the sample.
#'
#' Setting `kernel = "energy"` uses the negative pairwise distance, which
#' reproduces [bal_energy()] with its `"scaled_euclidean"` distance on the
#' same data and constraints. That equivalence is the reference point for the
#' method, since the energy kernel has an external implementation while the other
#' kernels do not.
#'
#' Characteristic function distance balancing belongs to the quadratic-program
#' family, which has no estimating equations, so a fit produces no
#' linearized-inference container and the guarantee is objective-level rather
#' than exact moment balance: without moment constraints the kernel objective
#' drives balance, and with them the constraint rows hold within tolerance. The
#' tolerance in [balance_terms()] relaxes any added moment constraints rather
#' than selecting an inexact solver. A tolerance supplied without moment
#' constraints has nothing to relax, so it is warned and ignored.
#'
#' @param kernel The covariate kernel the objective is built on, one of
#'   `"gaussian"`, `"matern"`, `"laplace"`, `"t"`, or `"energy"`.
#' @param smoothness The Matern smoothness order, one of `0.5`, `1.5`, or `2.5`.
#'   Used only by the Matern kernel.
#' @param degrees_of_freedom The degrees of freedom of the t kernel's frequency
#'   distribution, greater than two. Used only by the `"t"` kernel.
#' @param simulation_draws The number of Monte Carlo frequency projections for
#'   the `"t"` kernel.
#' @param improved Whether to add the between-group term of the improved variant
#'   for the average treatment effect with a discrete exposure.
#' @param weight_penalty The L2 penalty on the weights, which stabilizes the
#'   quadratic program.
#' @param min_weight The smallest permitted weight.
#' @param convergence_tolerance The quadratic-program solver tolerance, or `NULL`
#'   for the core default.
#' @param max_iterations The maximum solver iterations, or `NULL` for the core
#'   default.
#' @param ... Reserved for future extensions; must be empty. Tuning parameters
#'   must be passed by name.
#'
#' @return A `bal_cfd` specification, a [balance_method].
#'
#' @references
#' Wong, R. K. W. and Chan, K. C. G. (2018). Kernel-based covariate functional
#' balancing for observational studies. *Biometrika*, 105(1), 199-213.
#'
#' Huling, J. D. and Mak, S. (2024). Energy balancing of covariate
#' distributions. *Journal of Causal Inference*, 12(1), 20220029.
#'
#' @examples
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' df <- data.frame(
#'   exposure = rbinom(n, 1, plogis(0.5 * x1 - 0.5 * x2)),
#'   x1 = x1,
#'   x2 = x2
#' )
#' fit <- balance(df, exposure, c(x1, x2), method = bal_cfd())
#' fit
#'
#' @export
bal_cfd <- new_class(
  "bal_cfd",
  parent = quadratic_program_method,
  properties = list(
    kernel = class_character,
    smoothness = class_double,
    degrees_of_freedom = class_double,
    simulation_draws = class_integer,
    improved = class_logical
  ),
  constructor = function(
    kernel = c("gaussian", "matern", "laplace", "t", "energy"),
    smoothness = 1.5,
    degrees_of_freedom = 5,
    simulation_draws = 5000,
    improved = TRUE,
    weight_penalty = 1e-4,
    min_weight = 1e-8,
    convergence_tolerance = NULL,
    max_iterations = NULL,
    ...
  ) {
    rlang::check_dots_empty()
    kernel <- rlang::arg_match(kernel)
    smoothness <- vctrs::vec_cast(smoothness, double(), x_arg = "smoothness")
    degrees_of_freedom <- vctrs::vec_cast(
      degrees_of_freedom,
      double(),
      x_arg = "degrees_of_freedom"
    )
    simulation_draws <- vctrs::vec_cast(
      simulation_draws,
      integer(),
      x_arg = "simulation_draws"
    )
    weight_penalty <- vctrs::vec_cast(
      weight_penalty,
      double(),
      x_arg = "weight_penalty"
    )
    min_weight <- vctrs::vec_cast(min_weight, double(), x_arg = "min_weight")
    if (!is.null(convergence_tolerance)) {
      convergence_tolerance <- vctrs::vec_cast(
        convergence_tolerance,
        double(),
        x_arg = "convergence_tolerance"
      )
    }
    if (!is.null(max_iterations)) {
      max_iterations <- vctrs::vec_cast(
        max_iterations,
        integer(),
        x_arg = "max_iterations"
      )
    }
    new_object(
      S7_object(),
      kernel = kernel,
      smoothness = smoothness,
      degrees_of_freedom = degrees_of_freedom,
      simulation_draws = simulation_draws,
      improved = improved,
      weight_penalty = weight_penalty,
      min_weight = min_weight,
      convergence_tolerance = convergence_tolerance,
      max_iterations = max_iterations
    )
  },
  validator = function(self) {
    if (
      length(self@smoothness) != 1 ||
        is.na(self@smoothness) ||
        !any(vapply(
          c(0.5, 1.5, 2.5),
          function(v) isTRUE(all.equal(self@smoothness, v)),
          logical(1)
        ))
    ) {
      return("@smoothness must be one of 0.5, 1.5, or 2.5")
    }
    if (
      length(self@degrees_of_freedom) != 1 ||
        is.na(self@degrees_of_freedom) ||
        self@degrees_of_freedom <= 2
    ) {
      return("@degrees_of_freedom must be a single number greater than two")
    }
    if (
      length(self@simulation_draws) != 1 ||
        is.na(self@simulation_draws) ||
        self@simulation_draws < 1L
    ) {
      return("@simulation_draws must be a single positive whole number")
    }
    if (length(self@improved) != 1 || is.na(self@improved)) {
      return("@improved must be a single logical value")
    }
    if (
      length(self@weight_penalty) != 1 ||
        is.na(self@weight_penalty) ||
        self@weight_penalty < 0
    ) {
      return("@weight_penalty must be a single non-negative number")
    }
    if (
      length(self@min_weight) != 1 ||
        is.na(self@min_weight) ||
        self@min_weight < 0
    ) {
      return("@min_weight must be a single non-negative number")
    }
  }
)

method(method_label, bal_cfd) <- function(method) {
  "Characteristic function distance balancing"
}

method(supported_exposure_types, bal_cfd) <- function(method) {
  c("binary", "categorical")
}

method(supported_estimands, bal_cfd) <- function(method, exposure_type) {
  switch(
    exposure_type,
    # The overlap estimand is legal only for the covariate balancing propensity
    # score, so kernel balancing offers the discrete estimands of the rest of the
    # quadratic-program family and no continuous support.
    binary = c("ate", "att", "atu"),
    categorical = c("ate", "att")
  )
}

# The quadratic-program family has no estimating equations for any exposure type
# or constraint set, so the answer is always FALSE. The context arguments are
# accepted so the generic call shape matches the estimating-equation family.
method(supports_estimating_equations, bal_cfd) <- function(
  method,
  ...,
  exposure_type = NULL,
  constraints = NULL
) {
  rlang::check_dots_empty()
  FALSE
}

# Assemble the Rust option list. The worker-thread count and the quadratic-program
# backend are resolved on the R side and passed on every call; tuning parameters
# left at the core default are dropped so the solver applies its own.
cfd_options <- function(method) {
  options <- list(
    threads = resolve_threads(),
    backend = resolve_qp_backend()
  )
  if (!is.null(method@convergence_tolerance)) {
    options$convergence_tolerance <- method@convergence_tolerance
  }
  if (!is.null(method@max_iterations)) {
    options$max_iterations <- as.integer(method@max_iterations)
  }
  options
}

# The column-major p by n_draws frequency projections for the t kernel, drawn
# under R's random number generator so a fixed seed reproduces them. Each column
# is a multivariate t frequency vector: standard normal coordinates scaled by
# sqrt(df / u) with u a chi-square on df degrees of freedom shared across the
# column, the spectral draw of the t kernel's characteristic function. Returns a
# zero-length vector for the other kernels, which build their kernel internally.
cfd_projection <- function(method, p) {
  if (!identical(method@kernel, "t")) {
    return(numeric(0))
  }
  n_draws <- method@simulation_draws
  df <- method@degrees_of_freedom
  normals <- matrix(stats::rnorm(p * n_draws), nrow = p, ncol = n_draws)
  radial <- sqrt(df / stats::rchisq(n_draws, df))
  as.numeric(sweep(normals, 2, radial, "*"))
}

method(fit_method, bal_cfd) <- function(method, prepared) {
  enforce <- requests_moments(prepared$constraints)
  if (!enforce && has_positive_tolerance(prepared$constraints)) {
    warn_ignored_tolerance()
  }

  n <- prepared$n
  s <- prepared$sampling_weights
  covs <- distance_covariates(prepared$data, prepared$covariates)
  z <- prepared$matrix
  levels <- prepared$exposure_levels
  groups <- prepared$groups
  estimand <- prepared$estimand
  focal <- prepared$focal_level
  key <- prepared$exposure_key

  # The moment-constraint rows hold each reweighted group's weighted mean of the
  # standardized covariate columns within a tolerance band of a target, exactly as
  # the rest of the quadratic-program family assembles them. For the average
  # treatment effect the target is the pooled sample mean; for a focal estimand it
  # is the focal group's mean. Numeric columns cross at unit scale, so their band
  # is the tolerance directly; indicator columns cross raw, so their band scales by
  # the column standard deviation. The kernel objective alone drives balance when
  # no moment constraints are requested.
  if (enforce) {
    if (identical(estimand, "ate")) {
      targets <- target_means(z, seq_len(n), s)
    } else {
      targets <- target_means(z, groups[[focal]], s)
    }
    tols <- solver_box(z, prepared$tolerances)
    moment_covs <- z
  } else {
    targets <- numeric(0)
    tols <- numeric(0)
    moment_covs <- matrix(numeric(0), nrow = n, ncol = 0)
  }

  bw_scale <- 1
  t_proj <- cfd_projection(method, ncol(covs))
  options <- cfd_options(method)

  if (identical(prepared$exposure_type, "binary")) {
    # The focal group is coded one and held fixed; the average treatment effect
    # codes the second level as the reweighting reference. Both focal estimands
    # map to the core's focal solve, which holds the coded-one group.
    if (identical(estimand, "ate")) {
      treat <- as.integer(key == levels[[2]])
      core_estimand <- "ate"
      nvar <- n
      n_group_rows <- 2L
    } else {
      treat <- as.integer(key == focal)
      core_estimand <- "att"
      nvar <- n - length(groups[[focal]])
      n_group_rows <- 1L
    }
    result <- solve_cfd(
      covs,
      treat,
      s,
      method@kernel,
      bw_scale,
      method@smoothness,
      t_proj,
      method@improved,
      core_estimand,
      moment_covs,
      targets,
      tols,
      method@min_weight,
      method@weight_penalty,
      options
    )
  } else {
    treat_idx <- match(key, levels) - 1L
    if (identical(estimand, "ate")) {
      focal_idx <- 0L
      core_estimand <- "ate"
      nvar <- n
      n_group_rows <- length(levels)
    } else {
      focal_idx <- match(focal, levels) - 1L
      core_estimand <- "att"
      nvar <- n - length(groups[[focal]])
      n_group_rows <- length(levels) - 1L
    }
    result <- solve_cfd_multi(
      covs,
      as.integer(treat_idx),
      as.integer(focal_idx),
      s,
      method@kernel,
      bw_scale,
      method@smoothness,
      t_proj,
      method@improved,
      core_estimand,
      moment_covs,
      targets,
      tols,
      method@min_weight,
      method@weight_penalty,
      options
    )
  }

  duals <- energy_duals_frame(result$duals, nvar, n_group_rows)
  assemble_cfd(result, method, prepared, duals, approximate = !enforce)
}

# Renormalize each exposure group to its estimand target total, map the solver
# backend, and pack the fit result. The core returns weights on the per-group
# quadratic-program scale and applies the minimum-weight floor there; the reported
# weights renormalize each group to its estimand target total, a change of
# reporting convention, and the floor is re-applied on that reported scale so the
# documented minimum holds where the specs assert it. For the average treatment
# effect each group is scaled to its own sampling-weighted total; for a focal
# estimand every group is scaled to the focal total, which leaves the focal group
# at its base weight. The automatic backend fallback is recorded as
# "clarabel_fallback" so the behavior difference is never silent, and it is
# announced with an informational alert. The kernel objective drives approximate
# balance when no moment constraints are enforced, so the fit is flagged
# approximate in that case and the balance warning does not fire against a
# tolerance the fit does not target.
assemble_cfd <- function(result, method, prepared, duals, approximate) {
  s <- prepared$sampling_weights
  estimand <- prepared$estimand
  focal <- prepared$focal_level
  groups <- prepared$groups
  levels <- prepared$exposure_levels

  w <- as.numeric(result$weights)
  focal_estimand <- estimand %in% c("att", "atu")
  for (level in levels) {
    idx <- groups[[level]]
    target_sum <- if (focal_estimand) {
      sum(s[groups[[focal]]])
    } else {
      sum(s[idx])
    }
    current <- sum(s[idx] * w[idx])
    if (current > 0) {
      w[idx] <- w[idx] * (target_sum / current)
    }
  }
  w[w < method@min_weight] <- method@min_weight

  fell_back <- isTRUE(result$fell_back)
  backend <- result$solver_status
  solver_status <- if (fell_back) "clarabel_fallback" else backend
  if (fell_back && isTRUE(result$converged)) {
    alert_backend_fallback()
  }

  list(
    weights = w,
    coefficients = NULL,
    duals = duals,
    converged = isTRUE(result$converged),
    interrupted = isTRUE(result$interrupted),
    iterations = as.integer(result$iterations),
    objective = as.numeric(result$objective),
    solver_status = solver_status,
    # The backend's terminal status name, carried so the orchestrator maps an
    # infeasible constraint set and a hard solver failure to their designed
    # conditions rather than the generic convergence warning.
    status = result$status,
    estimating_equations = NULL,
    approximate = approximate,
    groups = groups
  )
}
