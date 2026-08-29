# A reference implementation of the stacked variance system `ipw()` reports,
# built here rather than read out of the package, and the evaluation counter the
# cost of that system is pinned with.
#
# The contrast rows of the stack are deterministic functions of the marginal
# means, so their bread rows are known without differencing anything and their
# meat is zero at the solution. An implementation is free to fill those rows in
# analytically and difference a narrower system, which is cheaper by two closure
# evaluations per contrast. What it is not free to do is change the answer, and
# that is what this file exists to hold it to.
#
# `ipw_reference_stack()` assembles the whole system the way the package
# assembles it today, differences every row of it, and returns the parameter
# vector and the covariance that fall out. It reuses the pieces the arithmetic is
# defined in, the container's own hooks, `deli::ee_glm()` for the outcome score,
# and the contrast formulas, so it restates no method's math; what it does not
# reuse is the assembly and the sandwich call under `ipw()`, which is the code an
# analytic contrast block changes. The blocks are stacked with `rbind()` rather
# than through `stack_psi_blocks()` for that reason, and the two are pinned to
# agree to the bit by `expect_stacked_psi_matches_rbind()` elsewhere in this
# suite.
#
# The deli call is copied argument for argument from `stacked_covariance()`,
# central differences at a step of 1e-6 with no pseudoinverse fallback, and the
# result is divided by the sample size and named the same way, because agreement
# is asserted with `identical()` rather than to a tolerance. Every one of those
# choices moves the last bits of the answer.
#
# The widths of the blocks travel back with the system. The evaluation count an
# analytic contrast block reaches is written from them, so the two readings of
# how many rows are deterministic cannot come apart.
ipw_reference_stack <- function(
  container,
  outcome_mod,
  frame,
  exposure_name,
  levels,
  categorical = FALSE,
  by = NULL,
  joint = NULL,
  sampling_weights = NULL,
  focal_level = NULL
) {
  n <- nrow(frame)
  family <- stats::family(outcome_mod)
  continuous <- is_gaussian_outcome(outcome_mod)
  distribution <- deli_distribution(family)
  outcome <- resolve_outcome_response(outcome_mod)
  design <- stats::model.matrix(outcome_mod)
  offset <- outcome_mod$offset

  pieces <- lapply(
    resolve_level_values(frame[[exposure_name]], levels),
    function(value) {
      fixed_exposure_pieces(
        outcome_mod,
        frame,
        exposure_name,
        value,
        offset = offset
      )
    }
  )

  sampling <- sampling_weights %||% rep(1, n)
  key <- as.character(frame[[exposure_name]])
  groups <- stats::setNames(
    lapply(levels, function(level) which(key == level)),
    levels
  )
  targets <- group_target_sums(sampling, groups, focal_level)
  tilt <- if (is.null(focal_level)) {
    sampling
  } else {
    sampling * (key == focal_level)
  }

  weight_parameters <- container@parameters
  p <- length(weight_parameters)
  coefficients <- stats::coef(outcome_mod)
  q <- length(coefficients)
  means <- vapply(
    pieces,
    function(piece) sum(tilt * piece$mu) / sum(tilt),
    numeric(1)
  )
  m <- length(means)
  contrasts <- ipw_contrast_values(means, continuous)
  effects <- ipw_contrast_names(continuous, if (categorical) levels else NULL)
  if (!is.null(joint)) {
    contrasts <- ipw_joint_values(joint, means, continuous)
    effects <- ipw_joint_names(joint)
  }
  k <- length(effects)

  by_stack <- ipw_by_stack(by, pieces, tilt, continuous, levels, categorical)
  m_by <- length(by_stack$means)
  k_by <- length(by_stack$contrasts)

  theta <- c(
    weight_parameters,
    coefficients,
    means,
    contrasts,
    by_stack$means,
    by_stack$contrasts
  )
  names(theta) <- c(
    paste0("theta_w", seq_len(p)),
    paste0("beta_", colnames(design)),
    ipw_mean_names(levels, categorical),
    effects,
    by_stack$mean_names,
    by_stack$contrast_names
  )

  rescale <- function(weights) {
    renormalize_group_weights(weights, sampling, groups, targets)
  }
  hooks_at <- make_hooks_cache(
    container,
    rescale,
    as.numeric(weight_parameters)
  )

  stacked_equations <- function(theta) {
    beta <- theta[p + seq_len(q)]
    mean_theta <- theta[p + q + seq_len(m)]
    contrast_theta <- theta[p + q + m + seq_len(k)]

    hooks <- hooks_at(as.numeric(theta[seq_len(p)]))
    score <- deli::ee_glm(
      beta,
      X = design,
      y = outcome,
      distribution = distribution,
      link = family$link,
      weights = hooks$weights * sampling,
      offset = offset
    )

    fixed <- lapply(seq_len(m), function(j) {
      eta <- as.numeric(pieces[[j]]$design %*% beta)
      if (!is.null(offset)) {
        eta <- eta + offset
      }
      family$linkinv(eta)
    })
    mean_rows <- do.call(
      rbind,
      lapply(seq_len(m), function(j) tilt * (fixed[[j]] - mean_theta[[j]]))
    )
    contrast_values <- if (is.null(joint)) {
      ipw_contrast_values(mean_theta, continuous) - contrast_theta
    } else {
      ipw_joint_row_values(joint, mean_theta, contrast_theta, continuous)
    }
    contrast_rows <- matrix(contrast_values, nrow = k, ncol = n)
    by_rows <- ipw_by_rows(
      by_stack = by_stack,
      fixed = fixed,
      mean_theta = theta[p + q + m + k + seq_len(m_by)],
      contrast_theta = theta[p + q + m + k + m_by + seq_len(k_by)],
      continuous = continuous,
      n = n
    )

    do.call(
      rbind,
      list(
        hooks$psi,
        score,
        mean_rows,
        contrast_rows,
        by_rows$mean,
        by_rows$contrast
      )
    )
  }

  covariance <- deli::compute_sandwich(
    stacked_equations,
    theta,
    deriv_method = "capprox",
    dx = 1e-6,
    allow_pinv = FALSE
  ) /
    n
  dimnames(covariance) <- list(names(theta), names(theta))

  list(
    theta = theta,
    vcov = covariance,
    width = length(theta),
    # The rows a deterministic function of the means fills: the whole-sample
    # contrast block, and every stratum and stratum-against-stratum block a
    # `.by` request adds.
    deterministic = k + k_by
  )
}

# Assert that `ipw()` reports the system the reference builds, to the bit, and
# that the block of it the reported effects are read from is the block the
# reference produces.
#
# All three surfaces are checked because none of them subsumes the others. The
# fit's copy is the whole system, the reported block is the slice a caller reads
# standard errors and confidence intervals off, and the parameter vector is what
# both are indexed by. `stacked_covariance()` renames the covariance it returns,
# and `attach_effect_covariance()` renames it again, so a block read from the
# wrong rows would still carry the right labels; comparing values under the
# observed dimnames is what makes the assertion about the numbers.
expect_ipw_matches_reference_stack <- function(result, reference, keys) {
  testthat::expect_identical(result$fit$theta, reference$theta)
  testthat::expect_identical(result$fit$vcov, reference$vcov)

  observed <- attr(result$estimates, "ipw_vcov")
  expected <- reference$vcov[keys, keys, drop = FALSE]
  dimnames(expected) <- dimnames(observed)
  testthat::expect_identical(observed, expected)

  invisible(result)
}

# The number of times a call assembles the stacked estimating functions, which
# is the number of times the sandwich evaluates the closure.
#
# deli's central difference evaluates the closure once at the fitted parameters
# and twice more per differenced coordinate, once on each side, so a system of
# width S written out in full costs 2S + 1 evaluations. Rows the implementation
# fills in analytically are neither differenced nor carried, so each of them
# saves two.
#
# The count is read at `stack_psi_blocks()`, the one place every route assembles
# its matrix, and the stand-in delegates to the real helper and returns its
# value, so the call running under it is the call the package performs.
expect_stacked_evaluations <- function(expr, expected) {
  assemble <- stack_psi_blocks
  evaluations <- 0L

  testthat::local_mocked_bindings(
    stack_psi_blocks = function(blocks, n) {
      evaluations <<- evaluations + 1L
      assemble(blocks, n)
    }
  )

  value <- force(expr)

  testthat::expect_identical(evaluations, expected)

  invisible(value)
}
