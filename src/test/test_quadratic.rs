use crate::{PolyError, Polynomial};
use is_close::is_close;
use num::complex::Complex64;
use std::f64::EPSILON;

// GSL's tests use this tolerance
const EPS: f64 = 100.0 * EPSILON;

//  ========================================== GSL Tests ==========================================

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_quadratic1() {
    let p = Polynomial::new(vec![26.0, -20.0, 4.0]);

    assert!(matches!(
        p.solve_real_quadratic().unwrap_err(),
        PolyError::NoRealRoots
    ));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_quadratic2() {
    let p = Polynomial::new(vec![25.0, -20.0, 4.0]);
    let y = p.solve_real_quadratic().unwrap();
    let expected = vec![2.5];

    assert_eq!(y.len(), 1);
    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_quadratic3() {
    let p = Polynomial::new(vec![21.0, -20.0, 4.0]);
    let y = p.solve_real_quadratic().unwrap();
    let expected = vec![3.5, 1.5];

    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
    assert!(is_close!(y[1], expected[1], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_quadratic4() {
    let p = Polynomial::new(vec![0.0, 7.0, 4.0]);
    let y = p.solve_real_quadratic().unwrap();
    let expected = vec![0.0, -1.75];

    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
    assert!(is_close!(y[1], expected[1], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_quadratic5() {
    let p = Polynomial::new(vec![-20.0, 0.0, 5.0]);
    let y = p.solve_real_quadratic().unwrap();
    let expected = vec![2.0, -2.0];

    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
    assert!(is_close!(y[1], expected[1], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_quadratic6() {
    let p = Polynomial::new(vec![-21.0, 3.0, 0.0]);
    let y = p.solve_real_quadratic().unwrap();
    let expected = vec![7.0];

    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_quadratic7() {
    let p = Polynomial::new(vec![1.0, 0.0, 0.0]);

    assert!(matches!(
        p.solve_real_quadratic().unwrap_err(),
        PolyError::ConstantPoly
    ));
}

// ===============================================================================================

#[test]
fn test_solve_real_quadratic_wrong_order() {
    let p = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]);

    assert!(matches!(
        p.solve_real_quadratic().unwrap_err(),
        PolyError::IncorrectOrder(2)
    ));
}

#[test]
fn test_solve_real_quadratic_complex_coefs() {
    let p = Polynomial::new(vec![
        Complex64::new(1.0, 2.0),
        Complex64::new(3.0, 4.0),
        Complex64::new(5.0, 6.0),
    ]);

    assert!(matches!(
        p.solve_real_quadratic().unwrap_err(),
        PolyError::NotRealCoefficients
    ));
}
