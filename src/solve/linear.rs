use crate::{PolyError, Result};

/// Solves a linear equation ax+b = 0, with real coefficients, returning a Vec with the
/// found 0-2 real roots.
pub(crate) fn solve_real_linear(a: f64, b: f64) -> Result<f64> {
    match a {
        0.0 => Err(PolyError::ConstantPoly),
        _ => Ok(-b / a),
    }
}
