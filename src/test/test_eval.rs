use crate::Polynomial;
use is_close::is_close;
use std::f64::{EPSILON, INFINITY, NAN};

// GSL's tests use this tolerance
const EPS: f64 = 100.0 * EPSILON;

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_eval1() {
    let p = Polynomial::new(vec![1.0, 0.5, 0.3]);
    let x = 0.5;

    assert_eq!(p.order, 2);
    assert!(is_close!(
        p.eval(x),
        1.0 + 0.5 * x + 0.3 * x * x,
        rel_tol = EPS
    ));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_eval2() {
    let p = Polynomial::new(vec![
        1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0,
    ]);
    let x = 1.0;

    assert_eq!(p.order, 10);
    assert!(is_close!(p.eval(x), 1.0, rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_eval_deriv1() {
    let p = Polynomial::new(vec![1.0, -2.0, 3.0, -4.0, 5.0, -6.0]);
    let x = -0.5;

    let derivs = p.eval_derivs(x, 6);

    assert_eq!(p.order, 5);
    assert!(is_close!(derivs[0], 3.75, rel_tol = EPS));
    assert!(is_close!(derivs[1], -12.375, rel_tol = EPS));
    assert!(is_close!(derivs[2], 48.0, rel_tol = EPS));
    assert!(is_close!(derivs[3], -174.0, rel_tol = EPS));
    assert!(is_close!(derivs[4], 480.0, rel_tol = EPS));
    assert!(is_close!(derivs[5], -720.0, rel_tol = EPS));
}

#[test]
fn test_eval_nan_infinity() {
    let p = Polynomial::new(vec![1.0, 0.5, 0.3]);
    assert!(p.eval(INFINITY).is_infinite());
    assert!(p.eval(NAN).is_nan());
}
