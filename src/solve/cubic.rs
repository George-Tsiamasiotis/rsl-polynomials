use crate::Result;

/// Solves a **monic** cubic equation x³+bx²+cx+d= 0 with real coefficients, returning a Vec with the found 0-3
/// real roots. In the case of a=0, solving is passed to the quadratic equation solver, and in the
/// case of a=b=0, solving is passed to the linear equation solver.
#[allow(unused_variables)]
pub(crate) fn solve_real_cubic(b: f64, c: f64, d: f64) -> Result<Vec<f64>> {
    todo!()
}
