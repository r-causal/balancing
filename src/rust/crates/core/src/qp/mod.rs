//! Quadratic-program specifications and solver backends.
//!
//! The energy, characteristic-function-distance, and stable-balancing-weights
//! methods all reduce to minimizing `0.5 x' P x + q' x` subject to
//! `l <= A x <= u`. This module owns the specification shared by those methods
//! and the backends that solve it. The OSQP backend is the default and handles
//! every method, including the indefinite quadratic form that energy balancing
//! produces; the Clarabel backend is an interior-point alternative restricted to
//! positive-semidefinite problems, which it enforces by rejecting a spec marked
//! indefinite before any factorization.
//!
//! Following the boundary convention, `P` is stored already doubled: a method
//! whose loss is `x' M x + q' x` places `2 M` in the spec so the solver's
//! `0.5 x' P x` matches the method's `x' M x`.

pub mod osqp;

#[cfg(feature = "qp-clarabel")]
pub mod clarabel;

/// The quadratic term of the objective.
///
/// A dense matrix is stored column-major and only its symmetric content matters;
/// a diagonal matrix stores just the diagonal, which the stable-balancing-weights
/// method uses.
#[derive(Debug, Clone)]
pub enum PMat {
    /// Column-major `n` by `n` dense matrix.
    Dense(Vec<f64>),
    /// The `n` diagonal entries of a diagonal matrix.
    Diagonal(Vec<f64>),
}

/// What is known about the definiteness of the quadratic term.
///
/// The tag is set by the method that assembles the spec, not rediscovered by the
/// backend, because an eigendecomposition of the dense `n` by `n` term is far
/// more expensive than the solve. Energy balancing tags its form `Indefinite`
/// from the mathematics; the diagonal and kernel forms of the other methods are
/// `Psd` by construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Convexity {
    /// Positive semidefinite; both backends may solve it.
    Psd,
    /// Known indefinite; only the ADMM backend is admissible.
    Indefinite,
    /// Definiteness not asserted; an interior-point backend may attempt it and
    /// fail, and the caller falls back.
    Unverified,
}

/// A quadratic program in `n` variables with `m` two-sided linear constraints.
///
/// The constraint matrix `a` is stored in compressed sparse column form as
/// `(indptr, indices, values)` describing an `m` by `n` matrix, the shape both
/// backends consume. `l` and `u` have length `m`.
#[derive(Debug, Clone)]
pub struct QpSpec {
    /// Number of decision variables.
    pub n: usize,
    /// Number of constraints (rows of `a`).
    pub m: usize,
    /// Quadratic term, already doubled.
    pub p: PMat,
    /// Linear term.
    pub q: Vec<f64>,
    /// Column pointers of `a`, length `n + 1`.
    pub a_indptr: Vec<usize>,
    /// Row indices of `a`, ascending within each column.
    pub a_indices: Vec<usize>,
    /// Nonzero values of `a`, aligned with `a_indices`.
    pub a_values: Vec<f64>,
    /// Lower bounds.
    pub l: Vec<f64>,
    /// Upper bounds.
    pub u: Vec<f64>,
    /// Definiteness of `p`.
    pub convexity: Convexity,
}

/// Terminal state a backend reports for a solve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QpStatus {
    /// Optimality tolerances met.
    Solved,
    /// Optimality met only at the looser inaccurate tolerances.
    SolvedInaccurate,
    /// The iteration cap was reached before convergence.
    MaxIter,
    /// The problem was certified primal infeasible.
    PrimalInfeasible,
    /// The problem was certified dual infeasible (unbounded).
    DualInfeasible,
    /// The backend judged the quadratic form non-convex.
    NonConvex,
    /// A user interrupt stopped the solve.
    Interrupted,
    /// A wall-clock limit stopped the solve.
    TimeLimit,
    /// The backend broke down numerically before reaching a verdict.
    NumericalError,
    /// The backend stopped because its iterates stopped making progress.
    InsufficientProgress,
    /// The backend returned without having run a solve at all.
    NotSolved,
}

impl QpStatus {
    /// A lowercase, snake-case name for the status, carried to R for diagnostics.
    pub fn as_str(self) -> &'static str {
        match self {
            QpStatus::Solved => "solved",
            QpStatus::SolvedInaccurate => "solved_inaccurate",
            QpStatus::MaxIter => "max_iter",
            QpStatus::PrimalInfeasible => "primal_infeasible",
            QpStatus::DualInfeasible => "dual_infeasible",
            QpStatus::NonConvex => "non_convex",
            QpStatus::Interrupted => "interrupted",
            QpStatus::TimeLimit => "time_limit",
            QpStatus::NumericalError => "numerical_error",
            QpStatus::InsufficientProgress => "insufficient_progress",
            QpStatus::NotSolved => "not_solved",
        }
    }

    /// Whether the status corresponds to a usable primal solution.
    pub fn is_solved(self) -> bool {
        matches!(self, QpStatus::Solved | QpStatus::SolvedInaccurate)
    }
}

/// Outcome of a quadratic-program solve.
#[derive(Debug, Clone)]
pub struct QpSolution {
    /// The primal solution, length `n`.
    pub x: Vec<f64>,
    /// The constraint dual variables, length `m`.
    pub duals: Vec<f64>,
    /// Terminal status.
    pub status: QpStatus,
    /// Iterations performed.
    pub iterations: usize,
    /// Objective value `0.5 x' P x + q' x` at the solution.
    pub obj: f64,
    /// Primal residual reported by the backend.
    pub pri_res: f64,
    /// Dual residual reported by the backend.
    pub dua_res: f64,
    /// Whether the solve stopped on a user interrupt.
    pub interrupted: bool,
}

/// Reasons a backend declines or fails to produce a solution.
#[derive(Debug, Clone)]
pub enum QpError {
    /// The backend does not accept an indefinite quadratic form.
    Indefinite,
    /// An explicitly requested backend is not compiled into this build; the field
    /// names the missing backend and its Cargo feature.
    BackendUnavailable(&'static str),
    /// Backend setup rejected the data, with a message.
    Setup(String),
}

impl std::fmt::Display for QpError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            QpError::Indefinite => write!(
                f,
                "the quadratic form is indefinite; this backend solves only positive-semidefinite problems"
            ),
            QpError::BackendUnavailable(name) => write!(
                f,
                "the `{name}` backend is not compiled into this build; rebuild with the `qp-clarabel` feature or choose another backend"
            ),
            QpError::Setup(msg) => write!(f, "quadratic-program setup failed: {msg}"),
        }
    }
}

/// Which backend a positive-semidefinite spec routes to.
///
/// The default is `Auto`: osqp solves first, and on an osqp primal-infeasibility
/// certificate the identical spec is re-solved with clarabel, because osqp can
/// falsely certify infeasibility on feasible but large or ill-scaled problems
/// that the interior-point backend handles. An explicit choice pins one backend
/// and disables the fallback.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum QpBackendChoice {
    /// osqp primary with an automatic clarabel fallback on an osqp
    /// primal-infeasibility certificate.
    #[default]
    Auto,
    /// osqp only; a primal-infeasibility certificate is returned as it stands.
    Osqp,
    /// clarabel only.
    Clarabel,
}

/// Tuning shared by the quadratic-program backends.
#[derive(Debug, Clone, Copy)]
pub struct QpOptions {
    /// Absolute convergence tolerance.
    pub eps_abs: f64,
    /// Relative convergence tolerance.
    pub eps_rel: f64,
    /// Maximum solver iterations.
    pub max_iter: usize,
    /// Whether to polish the ADMM solution.
    pub polish: bool,
    /// Whether OSQP adapts its step size.
    pub adaptive_rho: bool,
    /// OSQP iterations solved between interrupt checks. Clarabel takes a
    /// termination callback and polls once per iteration, so this does not
    /// reach it.
    pub chunk_iters: usize,
    /// Backend routing for a positive-semidefinite spec.
    pub backend: QpBackendChoice,
}

impl Default for QpOptions {
    fn default() -> Self {
        Self {
            eps_abs: 1e-8,
            eps_rel: 1e-8,
            max_iter: 200_000,
            polish: true,
            adaptive_rho: true,
            chunk_iters: 2_000,
            backend: QpBackendChoice::Auto,
        }
    }
}

/// A quadratic-program backend.
///
/// The trait lets the benchmark harness drive OSQP and Clarabel through the same
/// interface on the same spec. `interrupt` is polled while the solve runs and a
/// pending interrupt stops it, reporting `Interrupted` and setting the solution's
/// `interrupted` flag. How often it is polled is the backend's own business:
/// Clarabel takes a termination callback and polls every iteration, while OSQP has
/// no mid-solve callback and polls at the boundaries of bounded iteration chunks.
pub trait QpBackend {
    /// A short name carried to R as the solver identity.
    fn name(&self) -> &'static str;

    /// Solve `spec` under `opts`, polling `interrupt` for a pending user
    /// interrupt.
    fn solve(
        &self,
        spec: &QpSpec,
        opts: &QpOptions,
        interrupt: &dyn Fn() -> bool,
    ) -> Result<QpSolution, QpError>;
}

/// A solved positive-semidefinite spec together with the backend that produced
/// the iterate and whether the automatic fallback engaged.
#[derive(Debug, Clone)]
pub struct RoutedSolution {
    /// The backend solution.
    pub solution: QpSolution,
    /// The backend whose iterate is returned, `"osqp"` or `"clarabel"`.
    pub backend: &'static str,
    /// Whether the automatic clarabel fallback engaged after an osqp
    /// primal-infeasibility certificate.
    pub fell_back: bool,
}

/// Solve a positive-semidefinite spec under the routed backend policy.
///
/// `Auto` solves with osqp and, on an osqp primal-infeasibility certificate for a
/// spec that is not tagged indefinite, re-solves the identical spec with clarabel
/// and returns the clarabel result, since osqp can falsely certify infeasibility
/// on feasible but large or ill-scaled instances. An explicit backend disables
/// the fallback: `Osqp` returns the osqp certificate as it stands, and `Clarabel`
/// solves directly. When the `qp-clarabel` feature is disabled the routing
/// degrades to osqp only, preserving the pre-fallback behavior.
///
/// This routing is for the positive-semidefinite family (stable balancing weights
/// and future PSD methods); the indefinite energy form must call osqp directly,
/// since clarabel cannot accept it.
pub fn solve_psd(
    spec: &QpSpec,
    opts: &QpOptions,
    interrupt: &dyn Fn() -> bool,
) -> Result<RoutedSolution, QpError> {
    match opts.backend {
        QpBackendChoice::Osqp => solve_with_osqp(spec, opts, interrupt),
        QpBackendChoice::Clarabel => solve_with_clarabel(spec, opts, interrupt),
        QpBackendChoice::Auto => solve_auto(spec, opts, interrupt),
    }
}

/// Solve with osqp and report it as the producing backend.
fn solve_with_osqp(
    spec: &QpSpec,
    opts: &QpOptions,
    interrupt: &dyn Fn() -> bool,
) -> Result<RoutedSolution, QpError> {
    let solution = osqp::Osqp.solve(spec, opts, interrupt)?;
    Ok(RoutedSolution {
        solution,
        backend: "osqp",
        fell_back: false,
    })
}

/// Solve directly with clarabel when the feature is compiled in; without it, the
/// only backend available is osqp, so the routing degrades to it.
#[cfg(feature = "qp-clarabel")]
fn solve_with_clarabel(
    spec: &QpSpec,
    opts: &QpOptions,
    interrupt: &dyn Fn() -> bool,
) -> Result<RoutedSolution, QpError> {
    let solution = clarabel::Clarabel.solve(spec, opts, interrupt)?;
    Ok(RoutedSolution {
        solution,
        backend: "clarabel",
        fell_back: false,
    })
}

// Without the feature, clarabel is not linked. An explicit clarabel request is a
// deliberate choice, so it errors naming the missing feature rather than silently
// substituting osqp; the automatic path degrades quietly instead, since it did not
// ask for clarabel by name.
#[cfg(not(feature = "qp-clarabel"))]
fn solve_with_clarabel(
    _spec: &QpSpec,
    _opts: &QpOptions,
    _interrupt: &dyn Fn() -> bool,
) -> Result<RoutedSolution, QpError> {
    Err(QpError::BackendUnavailable("clarabel"))
}

/// osqp primary with a clarabel fallback on a primal-infeasibility certificate.
#[cfg(feature = "qp-clarabel")]
fn solve_auto(
    spec: &QpSpec,
    opts: &QpOptions,
    interrupt: &dyn Fn() -> bool,
) -> Result<RoutedSolution, QpError> {
    let osqp_solution = osqp::Osqp.solve(spec, opts, interrupt)?;
    // The fallback fires only on a genuine osqp infeasibility certificate for a
    // convex spec that was not interrupted; an interrupted solve is the user's
    // choice, and an indefinite spec is not clarabel's to accept.
    let certificate = osqp_solution.status == QpStatus::PrimalInfeasible
        && spec.convexity != Convexity::Indefinite
        && !osqp_solution.interrupted;
    if certificate {
        if let Ok(clarabel_solution) = clarabel::Clarabel.solve(spec, opts, interrupt) {
            return Ok(RoutedSolution {
                solution: clarabel_solution,
                backend: "clarabel",
                fell_back: true,
            });
        }
    }
    Ok(RoutedSolution {
        solution: osqp_solution,
        backend: "osqp",
        fell_back: false,
    })
}

#[cfg(not(feature = "qp-clarabel"))]
fn solve_auto(
    spec: &QpSpec,
    opts: &QpOptions,
    interrupt: &dyn Fn() -> bool,
) -> Result<RoutedSolution, QpError> {
    solve_with_osqp(spec, opts, interrupt)
}

/// Evaluate `0.5 x' P x + q' x` for the stored (already doubled) quadratic term.
pub fn objective(spec: &QpSpec, x: &[f64]) -> f64 {
    let n = spec.n;
    let quad = match &spec.p {
        PMat::Diagonal(d) => {
            let mut acc = 0.0;
            for i in 0..n {
                acc += d[i] * x[i] * x[i];
            }
            acc
        }
        PMat::Dense(mat) => {
            let mut acc = 0.0;
            for j in 0..n {
                let xj = x[j];
                if xj == 0.0 {
                    continue;
                }
                for i in 0..n {
                    acc += mat[j * n + i] * x[i] * xj;
                }
            }
            acc
        }
    };
    let lin: f64 = spec
        .q
        .iter()
        .zip(x.iter())
        .take(n)
        .map(|(qi, xi)| qi * xi)
        .sum();
    0.5 * quad + lin
}

/// Compressed sparse column data `(indptr, indices, values)` for the upper
/// triangle of the (already doubled) quadratic term, the form both backends
/// require. Explicit structural zeros are omitted; the diagonal is always
/// emitted so the solvers have a diagonal to regularize.
pub(crate) fn upper_triangular_csc(p: &PMat, n: usize) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
    let mut indptr = Vec::with_capacity(n + 1);
    let mut indices = Vec::new();
    let mut values = Vec::new();
    indptr.push(0);
    match p {
        PMat::Diagonal(d) => {
            for (j, &val) in d.iter().enumerate().take(n) {
                indices.push(j);
                values.push(val);
                indptr.push(indices.len());
            }
        }
        PMat::Dense(mat) => {
            for j in 0..n {
                for i in 0..=j {
                    let val = mat[j * n + i];
                    if i == j || val != 0.0 {
                        indices.push(i);
                        values.push(val);
                    }
                }
                indptr.push(indices.len());
            }
        }
    }
    (indptr, indices, values)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn objective_of_a_diagonal_form() {
        // 0.5 x' diag(2,4) x + q' x at x = (1, 1), q = (-1, -1):
        // 0.5 (2 + 4) - 2 = 1.
        let spec = QpSpec {
            n: 2,
            m: 0,
            p: PMat::Diagonal(vec![2.0, 4.0]),
            q: vec![-1.0, -1.0],
            a_indptr: vec![0, 0, 0],
            a_indices: vec![],
            a_values: vec![],
            l: vec![],
            u: vec![],
            convexity: Convexity::Psd,
        };
        assert!((objective(&spec, &[1.0, 1.0]) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn upper_triangular_csc_keeps_the_upper_triangle() {
        // Symmetric dense [[2, 3], [3, 5]] column-major.
        let p = PMat::Dense(vec![2.0, 3.0, 3.0, 5.0]);
        let (indptr, indices, values) = upper_triangular_csc(&p, 2);
        assert_eq!(indptr, vec![0, 1, 3]);
        assert_eq!(indices, vec![0, 0, 1]);
        assert_eq!(values, vec![2.0, 3.0, 5.0]);
    }

    #[test]
    fn upper_triangular_csc_always_emits_the_diagonal() {
        // A dense matrix with a zero diagonal still lists both diagonal slots so
        // the backends have a diagonal to regularize.
        let p = PMat::Dense(vec![0.0, 1.0, 1.0, 0.0]);
        let (indptr, indices, _values) = upper_triangular_csc(&p, 2);
        assert_eq!(indptr, vec![0, 1, 3]);
        assert_eq!(indices, vec![0, 0, 1]);
    }
}
