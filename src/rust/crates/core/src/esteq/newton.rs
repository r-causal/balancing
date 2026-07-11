//! Damped Newton method with ridge escalation and Armijo backtracking.

use faer::MatMut;

use super::{EsteqProblem, SolveOptions, SolveReport, Solver};
use crate::linalg::solve_symmetric_ridge;

/// Armijo sufficient-decrease constant.
const ARMIJO_C: f64 = 1e-4;
/// Maximum step halvings per iteration.
const MAX_BACKTRACKS: usize = 30;
/// Maximum ridge escalations before an iteration is abandoned.
const MAX_ESCALATIONS: usize = 12;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn sup_norm(a: &[f64]) -> f64 {
    a.iter().fold(0.0_f64, |m, x| m.max(x.abs()))
}

/// Symmetric matrix-vector product for a column-major `p` by `p` slice.
fn sym_matvec(h: &[f64], p: usize, x: &[f64], out: &mut [f64]) {
    for i in 0..p {
        let mut acc = 0.0;
        for j in 0..p {
            acc += h[j * p + i] * x[j];
        }
        out[i] = acc;
    }
}

/// A positive ridge seed scaled to the Hessian magnitude, used the first time a
/// zero ridge fails to yield a usable direction.
fn ridge_seed(h: &[f64], p: usize) -> f64 {
    let mut trace = 0.0;
    for i in 0..p {
        trace += h[i * p + i].abs();
    }
    let scale = (trace / p as f64).max(1.0);
    1e-10 * scale
}

/// Solve a smooth (or root-finding) estimating-equation problem by damped Newton.
pub fn solve<P: EsteqProblem>(
    problem: &P,
    beta: &mut [f64],
    opts: &SolveOptions,
    interrupt: &dyn Fn() -> bool,
) -> SolveReport {
    let p = problem.n_params();
    let smooth = problem.value(beta).is_some();

    let mut g = vec![0.0; p];
    let mut h = vec![0.0; p * p];
    let mut dir = vec![0.0; p];
    let mut beta_trial = vec![0.0; p];
    let mut merit_grad = vec![0.0; p];
    let mut g_trial = vec![0.0; p];

    let mut ridge_max = 0.0_f64;
    let mut iterations = 0;
    let mut converged = false;
    let mut interrupted = false;

    for iter in 0..opts.max_iter {
        iterations = iter;
        if interrupt() {
            interrupted = true;
            break;
        }

        let value = {
            let h_mat = MatMut::from_column_major_slice_mut(&mut h, p, p);
            problem.value_grad_hess(beta, &mut g, h_mat)
        };
        let grad_norm = sup_norm(&g);
        if grad_norm <= opts.grad_tol {
            converged = true;
            break;
        }

        // Newton direction: solve (H + ridge I) dir = -g, growing the ridge from
        // zero until the factorization succeeds and the direction descends.
        let mut ridge = 0.0;
        let mut have_dir = false;
        for escalation in 0..=MAX_ESCALATIONS {
            for k in 0..p {
                dir[k] = -g[k];
            }
            if solve_symmetric_ridge(&h, p, ridge, &mut dir) && dot(&g, &dir) < 0.0 {
                have_dir = true;
                break;
            }
            ridge = if ridge == 0.0 {
                ridge_seed(&h, p)
            } else {
                ridge * 10.0
            };
            ridge_max = ridge_max.max(ridge);
            if escalation == MAX_ESCALATIONS {
                break;
            }
        }
        if !have_dir {
            break;
        }

        // Directional derivative of the merit being reduced. For a smooth
        // objective it is g . dir; for a root-finding merit 0.5||g||^2 it is
        // (H g) . dir.
        let base = value.unwrap_or_else(|| 0.5 * dot(&g, &g));
        let dderiv = if smooth {
            dot(&g, &dir)
        } else {
            sym_matvec(&h, p, &g, &mut merit_grad);
            dot(&merit_grad, &dir)
        };

        let mut step = 1.0;
        let mut accepted = false;
        for _ in 0..MAX_BACKTRACKS {
            for k in 0..p {
                beta_trial[k] = beta[k] + step * dir[k];
            }
            let trial = if smooth {
                problem
                    .value(&beta_trial)
                    .expect("smooth problem returns a value")
            } else {
                problem.gradient(&beta_trial, &mut g_trial);
                0.5 * dot(&g_trial, &g_trial)
            };
            if trial <= base + ARMIJO_C * step * dderiv {
                accepted = true;
                break;
            }
            step *= 0.5;
        }
        // A descent step that no longer reduces the objective means the
        // predicted decrease has fallen below the objective's floating-point
        // resolution: the iterate is at the numerical optimum, so stop rather
        // than crawl in ever-smaller steps.
        if !accepted {
            break;
        }
        for k in 0..p {
            beta[k] += step * dir[k];
        }
    }

    // Report the gradient and objective at the returned parameters.
    problem.gradient(beta, &mut g);
    let grad_norm = sup_norm(&g);
    if grad_norm <= opts.grad_tol {
        converged = true;
    }
    let final_value = problem.value(beta).unwrap_or_else(|| 0.5 * dot(&g, &g));

    SolveReport {
        converged,
        interrupted,
        iterations,
        grad_norm,
        final_value,
        ridge_max,
        solver_used: Solver::Newton,
    }
}
