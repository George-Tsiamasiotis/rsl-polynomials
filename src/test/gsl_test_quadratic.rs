use crate::{PolyError, Polynomial};
use is_close::is_close;

// GSL's tests use this tolerance
const EPS: f64 = 100.0 * f64::EPSILON;

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_quadratic1() {
    let p = Polynomial::build(&vec![26.0, -20.0, 4.0]).unwrap();

    assert!(matches!(
        p.solve_real_quadratic().unwrap_err(),
        PolyError::NoRealRoots
    ));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_quadratic2() {
    let p = Polynomial::build(&vec![25.0, -20.0, 4.0]).unwrap();
    let y = p.solve_real_quadratic().unwrap();
    let expected = vec![2.5];

    assert_eq!(y.len(), 1);
    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_quadratic3() {
    let p = Polynomial::build(&vec![21.0, -20.0, 4.0]).unwrap();
    let y = p.solve_real_quadratic().unwrap();
    let expected = vec![3.5, 1.5];

    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
    assert!(is_close!(y[1], expected[1], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_quadratic4() {
    let p = Polynomial::build(&vec![0.0, 7.0, 4.0]).unwrap();
    let y = p.solve_real_quadratic().unwrap();
    let expected = vec![0.0, -1.75];

    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
    assert!(is_close!(y[1], expected[1], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_quadratic5() {
    let p = Polynomial::build(&vec![-20.0, 0.0, 5.0]).unwrap();
    let y = p.solve_real_quadratic().unwrap();
    let expected = vec![2.0, -2.0];

    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
    assert!(is_close!(y[1], expected[1], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_quadratic6() {
    let p = Polynomial::build(&vec![-21.0, 3.0, 0.0]).unwrap();
    let y = p.solve_real_quadratic().unwrap();
    let expected = vec![7.0];

    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_quadratic7() {
    let p = Polynomial::build(&vec![1.0, 0.0, 0.0]).unwrap();

    assert!(matches!(
        p.solve_real_quadratic().unwrap_err(),
        PolyError::ConstantPoly
    ));
}
