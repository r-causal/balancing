//! FISTA for composite problems: a smooth convex term plus an L1 penalty.
//!
//! The inexact (tolerance-based) entropy dual is the smooth dual plus a weighted
//! L1 term whose weights are the per-constraint tolerances. Proximal gradient
//! with Nesterov momentum handles it: the prox of the L1 term is coordinatewise
//! soft-thresholding. The Lipschitz constant is discovered by backtracking and
//! carried across iterations, and momentum restarts whenever a step increases
//! the full loss.

/// Maximum Lipschitz doublings per iteration.
const MAX_BACKTRACKS: usize = 10;

/// Result of a FISTA minimization.
#[derive(Debug, Clone, Copy)]
pub struct FistaReport {
    /// Whether the relative-loss stopping rule was met.
    pub converged: bool,
    /// Whether the solve stopped because a user interrupt was pending.
    pub interrupted: bool,
    /// Iterations performed.
    pub iterations: usize,
    /// Full loss (smooth plus L1) at the returned point.
    pub final_loss: f64,
    /// Sup norm of the smooth gradient at the returned point.
    pub grad_norm: f64,
}

fn soft_threshold(value: f64, threshold: f64) -> f64 {
    if value > threshold {
        value - threshold
    } else if value < -threshold {
        value + threshold
    } else {
        0.0
    }
}

fn l1_value(x: &[f64], l1: &[f64]) -> f64 {
    x.iter().zip(l1).map(|(xi, li)| li * xi.abs()).sum()
}

/// Minimize `smooth_value(x) + sum_j l1[j] * |x[j]|` in place from `x`.
///
/// `smooth_grad` writes the gradient of the smooth term. `l1` holds the
/// per-coordinate penalty weights (zero disables the penalty on that
/// coordinate). `interrupt` is polled each iteration.
pub fn minimize<Value, Grad>(
    x: &mut [f64],
    l1: &[f64],
    smooth_value: Value,
    smooth_grad: Grad,
    max_iter: usize,
    rel_tol: f64,
    interrupt: &dyn Fn() -> bool,
) -> FistaReport
where
    Value: Fn(&[f64]) -> f64,
    Grad: Fn(&[f64], &mut [f64]),
{
    let p = x.len();
    let mut y = x.to_vec();
    let mut x_new = vec![0.0; p];
    let mut grad_y = vec![0.0; p];

    let mut t = 1.0_f64;
    let mut lipschitz = 1.0_f64;
    let mut prev_loss = smooth_value(x) + l1_value(x, l1);
    let mut iterations = 0;
    let mut converged = false;
    let mut interrupted = false;

    for iter in 0..max_iter {
        iterations = iter + 1;
        if interrupt() {
            interrupted = true;
            break;
        }

        let f_y = smooth_value(&y);
        smooth_grad(&y, &mut grad_y);

        // Backtrack on the local Lipschitz estimate until the quadratic model
        // upper-bounds the smooth term at the proposed point.
        for _ in 0..MAX_BACKTRACKS {
            for j in 0..p {
                let step = y[j] - grad_y[j] / lipschitz;
                x_new[j] = soft_threshold(step, l1[j] / lipschitz);
            }
            let f_new = smooth_value(&x_new);
            let mut model = f_y;
            let mut quad = 0.0;
            for j in 0..p {
                let diff = x_new[j] - y[j];
                model += grad_y[j] * diff;
                quad += diff * diff;
            }
            model += 0.5 * lipschitz * quad;
            if f_new <= model + 1e-12 * f_y.abs() {
                break;
            }
            lipschitz *= 2.0;
        }

        let loss_new = smooth_value(&x_new) + l1_value(&x_new, l1);
        let t_new = 0.5 * (1.0 + (1.0 + 4.0 * t * t).sqrt());

        if loss_new > prev_loss {
            // Progress reversed: restart momentum from the new point.
            y.copy_from_slice(&x_new);
            t = 1.0;
        } else {
            let momentum = (t - 1.0) / t_new;
            for j in 0..p {
                y[j] = x_new[j] + momentum * (x_new[j] - x[j]);
            }
            t = t_new;
        }
        x.copy_from_slice(&x_new);

        // The stopping rule is relative to the loss, with a floor of one so a loss
        // that passes near zero does not demand an absolute change of zero. The
        // floor carries a calibration assumption: the tolerance was chosen for
        // losses of order one, which the entropy dual satisfies at the weight
        // scales this package normalizes to. It is not scale-free above that. A
        // global rescaling of the sampling weights shifts the dual by the log of
        // the factor, and a factor of a million adds about fourteen to a loss of
        // order one, loosening the rule by roughly the same multiple. The exact
        // path is unaffected, since Newton judges the gradient rather than the
        // loss and the gradient carries its own scale through
        // `EsteqProblem::residual_scale`; only the inexact path's stopping point
        // moves, and it moves in the direction of stopping sooner.
        let denom = prev_loss.abs().max(1.0);
        if (prev_loss - loss_new).abs() <= rel_tol * denom {
            converged = true;
            break;
        }
        prev_loss = loss_new;
    }

    let mut grad = vec![0.0; p];
    smooth_grad(x, &mut grad);
    let grad_norm = grad.iter().fold(0.0_f64, |m, gi| m.max(gi.abs()));

    FistaReport {
        converged,
        interrupted,
        iterations,
        final_loss: smooth_value(x) + l1_value(x, l1),
        grad_norm,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn no_interrupt() -> impl Fn() -> bool {
        || false
    }

    #[test]
    fn recovers_soft_threshold_solution() {
        // min 0.5 (x - 3)^2 + lambda |x| has minimizer 3 - lambda for lambda < 3.
        let lambda = 0.75;
        let mut x = [0.0];
        let report = minimize(
            &mut x,
            &[lambda],
            |v| 0.5 * (v[0] - 3.0).powi(2),
            |v, g| g[0] = v[0] - 3.0,
            5000,
            1e-14,
            &no_interrupt(),
        );
        assert!(report.converged);
        assert!((x[0] - (3.0 - lambda)).abs() < 1e-6, "x = {}", x[0]);
    }

    #[test]
    fn zero_penalty_minimizes_the_smooth_term() {
        // With no penalty FISTA minimizes the quadratic to its unconstrained min.
        let mut x = [0.0, 0.0];
        let report = minimize(
            &mut x,
            &[0.0, 0.0],
            |v| 0.5 * ((v[0] - 1.0).powi(2) + 4.0 * (v[1] + 2.0).powi(2)),
            |v, g| {
                g[0] = v[0] - 1.0;
                g[1] = 4.0 * (v[1] + 2.0);
            },
            5000,
            1e-14,
            &no_interrupt(),
        );
        assert!(report.converged);
        assert!((x[0] - 1.0).abs() < 1e-6, "x0 = {}", x[0]);
        assert!((x[1] + 2.0).abs() < 1e-6, "x1 = {}", x[1]);
    }
}
