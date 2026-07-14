//! Constraint assembly shared by the quadratic-program balancing methods.
//!
//! Energy, characteristic-function-distance, and stable-balancing-weights
//! balancing differ in their quadratic term but share the same constraint
//! structure: a box on each weight with a minimum-weight floor and pinned
//! zero-sampling-weight units, group-normalized sum constraints that fix each
//! reweighted group's total, and optional moment-constraint rows that hold a
//! covariate's weighted mean within a tolerance band. This module builds that
//! constraint matrix and applies the weight penalty to the quadratic diagonal,
//! so each method module supplies only its objective.

use crate::qp::PMat;

/// Threshold below which a sampling weight pins its unit's box to one, matching
/// the reference treatment of zero-sampling-weight units.
pub const ZERO_SW: f64 = 1e-8;

/// The magnitude below which a resolved weight is snapped to zero, applied only
/// when the minimum-weight floor is itself essentially zero.
pub const SHRINK_TO_ZERO: f64 = 1e-10;

/// Group-normalize the sampling weights so each level has mean one, matching the
/// per-group scaling the quadratic-program objectives assume. Units with a
/// negative level (excluded) keep their sampling weight unchanged.
pub fn group_normalized(s: &[f64], levels: &[i32], n_levels: usize) -> Vec<f64> {
    let mut sums = vec![0.0; n_levels];
    let mut counts = vec![0usize; n_levels];
    for (i, &g) in levels.iter().enumerate() {
        if g >= 0 {
            sums[g as usize] += s[i];
            counts[g as usize] += 1;
        }
    }
    let means: Vec<f64> = sums
        .iter()
        .zip(&counts)
        .map(|(&sum, &c)| if c > 0 { sum / c as f64 } else { 1.0 })
        .collect();
    levels
        .iter()
        .enumerate()
        .map(|(i, &g)| {
            if g >= 0 {
                let m = means[g as usize];
                if m > 0.0 { s[i] / m } else { s[i] }
            } else {
                s[i]
            }
        })
        .collect()
}

/// Expand an active-variable solution to the full unit order and apply the
/// reference post-processing: units outside the active set keep unit weight (the
/// focal group a focal estimand holds fixed), each active unit takes its solved
/// value, every weight is floored at `min_weight`, and negligible weights are
/// snapped to zero when the floor is itself essentially zero.
pub fn expand_and_floor(n: usize, active: &[usize], x: &[f64], min_weight: f64) -> Vec<f64> {
    let mut weights = vec![1.0; n];
    for (a, &i) in active.iter().enumerate() {
        weights[i] = x[a];
    }
    for w in weights.iter_mut() {
        if *w < min_weight {
            *w = min_weight;
        }
    }
    if min_weight.abs() < SHRINK_TO_ZERO {
        for w in weights.iter_mut() {
            if w.abs() < SHRINK_TO_ZERO {
                *w = 0.0;
            }
        }
    }
    weights
}

/// The compressed constraint matrix and bounds a builder produces:
/// `(m, indptr, indices, values, l, u)` describing an `m` by `n` matrix in
/// sparse column form together with its lower and upper bounds.
pub type CompressedConstraints = (usize, Vec<usize>, Vec<usize>, Vec<f64>, Vec<f64>, Vec<f64>);

/// A single two-sided constraint row, its coefficients over the decision
/// variables paired with the row bounds.
struct Row {
    coeffs: Vec<(usize, f64)>,
    l: f64,
    u: f64,
}

/// Accumulates constraint rows and compresses them into the sparse column form
/// the quadratic-program spec carries.
pub struct ConstraintBuilder {
    n: usize,
    rows: Vec<Row>,
}

impl ConstraintBuilder {
    /// Start a builder for `n` decision variables.
    pub fn new(n: usize) -> Self {
        Self {
            n,
            rows: Vec::new(),
        }
    }

    /// Add the identity box block: one row per variable bounding it in
    /// `[min_weight, +inf)`, or pinned to exactly one where `pinned` is set (a
    /// unit whose sampling weight is zero contributes nothing and is held fixed).
    pub fn add_box(&mut self, min_weight: f64, pinned: &[bool]) {
        for i in 0..self.n {
            let (l, u) = if pinned.get(i).copied().unwrap_or(false) {
                (1.0, 1.0)
            } else {
                (min_weight, f64::INFINITY)
            };
            self.rows.push(Row {
                coeffs: vec![(i, 1.0)],
                l,
                u,
            });
        }
    }

    /// Add a dense constraint row from a full-length coefficient vector, keeping
    /// only the nonzero entries.
    pub fn add_dense_row(&mut self, coeffs: &[f64], l: f64, u: f64) {
        let sparse: Vec<(usize, f64)> = coeffs
            .iter()
            .enumerate()
            .filter(|&(_, &v)| v != 0.0)
            .map(|(i, &v)| (i, v))
            .collect();
        self.rows.push(Row {
            coeffs: sparse,
            l,
            u,
        });
    }

    /// Compress the accumulated rows into `(m, indptr, indices, values, l, u)`
    /// describing the `m` by `n` constraint matrix in sparse column form.
    ///
    /// A row that is unbounded on both sides constrains nothing and is dropped,
    /// matching the reference assembler; every retained row keeps at least one
    /// finite bound. Row indices are ascending within each column.
    pub fn finish(self) -> CompressedConstraints {
        let kept: Vec<Row> = self
            .rows
            .into_iter()
            .filter(|r| r.l.is_finite() || r.u.is_finite())
            .collect();
        let m = kept.len();

        let mut columns: Vec<Vec<(usize, f64)>> = vec![Vec::new(); self.n];
        let mut l = Vec::with_capacity(m);
        let mut u = Vec::with_capacity(m);
        for (r, row) in kept.iter().enumerate() {
            for &(col, val) in &row.coeffs {
                columns[col].push((r, val));
            }
            l.push(row.l);
            u.push(row.u);
        }

        let mut indptr = Vec::with_capacity(self.n + 1);
        let mut indices = Vec::new();
        let mut values = Vec::new();
        indptr.push(0);
        for col in columns.iter_mut() {
            col.sort_by_key(|&(r, _)| r);
            for &(r, v) in col.iter() {
                indices.push(r);
                values.push(v);
            }
            indptr.push(indices.len());
        }

        (m, indptr, indices, values, l, u)
    }
}

/// Add the weight penalty `lambda * scale_i^2 / 2` to the diagonal of a dense
/// quadratic term stored column-major.
///
/// The penalty grows the effective sample size at the expense of balance by
/// discouraging variable weights; `scale` carries the per-variable group
/// normalization so the penalty is on the same footing as the objective.
pub fn add_diagonal_penalty(p: &mut [f64], n: usize, lambda: f64, scale: &[f64]) {
    if lambda == 0.0 {
        return;
    }
    for i in 0..n {
        p[i * n + i] += lambda * scale[i] * scale[i] / 2.0;
    }
}

/// Build a [`PMat::Dense`] holding twice the loss matrix, the doubled quadratic
/// term the spec convention requires. The method assembles its loss matrix, adds
/// any diagonal penalty, and passes it here to cross into the spec.
pub fn doubled_dense(mut loss: Vec<f64>) -> PMat {
    for v in loss.iter_mut() {
        *v *= 2.0;
    }
    PMat::Dense(loss)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn box_block_bounds_each_variable() {
        let mut b = ConstraintBuilder::new(3);
        b.add_box(1e-8, &[false, true, false]);
        let (m, indptr, indices, values, l, u) = b.finish();
        assert_eq!(m, 3);
        // Identity: column j has one entry in row j.
        assert_eq!(indptr, vec![0, 1, 2, 3]);
        assert_eq!(indices, vec![0, 1, 2]);
        assert_eq!(values, vec![1.0, 1.0, 1.0]);
        assert_eq!(l, vec![1e-8, 1.0, 1e-8]);
        assert_eq!(u[1], 1.0);
        assert!(u[0].is_infinite());
    }

    #[test]
    fn a_fully_unbounded_row_is_dropped() {
        let mut b = ConstraintBuilder::new(2);
        b.add_dense_row(&[1.0, 1.0], f64::NEG_INFINITY, f64::INFINITY);
        b.add_dense_row(&[1.0, 0.0], 0.0, 1.0);
        let (m, _indptr, _indices, _values, l, _u) = b.finish();
        assert_eq!(m, 1);
        assert_eq!(l, vec![0.0]);
    }

    #[test]
    fn dense_rows_compress_by_column_in_row_order() {
        // Box rows 0 and 1, then a dense sum row 2 over both columns.
        let mut b = ConstraintBuilder::new(2);
        b.add_box(0.0, &[false, false]);
        b.add_dense_row(&[2.0, 3.0], 1.0, 1.0);
        let (m, indptr, indices, values, _l, _u) = b.finish();
        assert_eq!(m, 3);
        // Column 0: row 0 (box, 1.0), row 2 (sum, 2.0).
        assert_eq!(indptr, vec![0, 2, 4]);
        assert_eq!(indices, vec![0, 2, 1, 2]);
        assert_eq!(values, vec![1.0, 2.0, 1.0, 3.0]);
    }

    #[test]
    fn penalty_lands_on_the_diagonal() {
        let mut p = vec![0.0; 4];
        add_diagonal_penalty(&mut p, 2, 4.0, &[1.0, 2.0]);
        // diag += lambda * scale^2 / 2 = 4 * 1 / 2 = 2 and 4 * 4 / 2 = 8.
        assert_eq!(p[0], 2.0);
        assert_eq!(p[3], 8.0);
        assert_eq!(p[1], 0.0);
    }
}
