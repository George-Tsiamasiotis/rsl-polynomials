//! Methods for solving linear, quadratic, cubic(todo) and quartic(maybe todo) equations.

use std::cmp::Ordering;

use crate::solve::linear::solve_real_linear;
use crate::{PolyError, Result};

/// Solves a quadratic equation ax²+bx+c = 0 with real coefficients, returning a Vec with the found 0-2
/// real roots. In the case of a=0, solving is passed to the linear equation solver.
pub(crate) fn solve_real_quadratic(a: f64, b: f64, c: f64) -> Result<Vec<f64>> {
    if a == 0.0 {
        return Ok(vec![solve_real_linear(b, c)?]);
    }

    let det = b.powi(2) - 4.0 * a * c;

    let ordering = match det.partial_cmp(&0.0) {
        Some(det) => det,
        None => unreachable!("NaN discriminant"),
    };

    match ordering {
        Ordering::Less => Err(PolyError::NoRealRoots),
        Ordering::Equal => {
            let x = -b / (2.0 * a);
            Ok(vec![x])
        }
        Ordering::Greater => {
            let x1 = (-b + det.sqrt()) / (2.0 * a);
            let x2 = (-b - det.sqrt()) / (2.0 * a);

            Ok(vec![x1, x2])
        }
    }
}
