use std::f64::consts::PI;

use crate::Result;

/// Solves a **depressed** cubic equation  t³+pt+q=0,  where t=x−b/3, awith real coefficients,
/// returning a Vec with the found 0-3 real roots.
///
/// a, b, c correspond to a polynomial x³ + ax² + bx + c.
pub(crate) fn solve_real_cubic(a: f64, b: f64, c: f64) -> Result<Vec<f64>> {
    let q = a.powi(2) - 3.0 * b;
    let r = 2.0 * a.powi(3) - 9.0 * a * b + 27.0 * c;

    let q_cap = q / 9.0;
    let r_cap = r / 54.0;

    let q_cap3 = q_cap.powi(3);
    let r_cap2 = r_cap.powi(2);

    let cq_cap3 = 2916.0 * q.powi(3);
    let cr_cap2 = 729.0 * r.powi(2);

    let mut ans: Vec<f64> = vec![f64::NAN; 3];

    // NOTE: This test is actually `r_cap2==q_cap3`, written in a form suitable for exact
    // computation with integers
    if (r_cap == 0.0) & (q_cap == 0.0) {
        let x = -a / 3.0;
        ans.fill(x);
        return Ok(vec![x, x, x]);
    } else if cr_cap2 == cq_cap3 {
        let sqrtq = q_cap.sqrt();

        if r > 0.0 {
            ans[0] = -2.0 * sqrtq - a / 3.0;
            ans[1] = sqrtq - a / 3.0;
            ans[2] = sqrtq - a / 3.0;
        } else {
            ans[0] = -sqrtq - a / 3.0;
            ans[1] = -sqrtq - a / 3.0;
            ans[2] = 2.0 * sqrtq - a / 3.0;
        }
    } else if r_cap2 < q_cap3 {
        let sgnr = r.signum();
        let ratio = sgnr * (r_cap2 / q_cap3).sqrt();
        let theta = ratio.acos();
        let norm = -2.0 * q_cap.sqrt();

        ans[0] = norm * (theta / 3.0).cos() - a / 3.0;
        ans[1] = norm * ((theta + 2.0 * PI) / 3.0).cos() - a / 3.0;
        ans[2] = norm * ((theta - 2.0 * PI) / 3.0).cos() - a / 3.0;
    } else {
        let sgnr = r.signum();
        let a_cap = -sgnr * (r_cap.abs() + (r_cap2 - q_cap3).sqrt()).powf(1.0 / 3.0);
        let b_cap = q / a_cap;
        let x = a_cap + b_cap - a / 3.0;
        ans.fill(x);
    }

    ans.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Ok(ans)
}
