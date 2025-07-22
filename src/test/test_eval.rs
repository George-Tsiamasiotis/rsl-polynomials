use crate::{PolyError, Polynomial};
use is_close::is_close;
use num::complex::Complex64;
use std::f64::EPSILON;

// GSL's tests use this tolerance
const EPS: f64 = 100.0 * EPSILON;

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_eval1() {
    let p = Polynomial::build(&vec![1.0, 0.5, 0.3]).unwrap();
    let x = 0.5;

    assert!(is_close!(
        p.eval(x),
        1.0 + 0.5 * x + 0.3 * x * x,
        rel_tol = EPS
    ));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_eval2() {
    let p = Polynomial::build(&vec![
        1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0,
    ])
    .unwrap();
    let x = 1.0;

    assert!(is_close!(p.eval(x), 1.0, rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_eval_deriv1() {
    let p = Polynomial::build(&vec![1.0, -2.0, 3.0, -4.0, 5.0, -6.0]).unwrap();
    let x = -0.5;

    let derivs = p.eval_derivs(x, 6);

    assert!(is_close!(derivs[0], 3.75, rel_tol = EPS));
    assert!(is_close!(derivs[1], -12.375, rel_tol = EPS));
    assert!(is_close!(derivs[2], 48.0, rel_tol = EPS));
    assert!(is_close!(derivs[3], -174.0, rel_tol = EPS));
    assert!(is_close!(derivs[4], 480.0, rel_tol = EPS));
    assert!(is_close!(derivs[5], -720.0, rel_tol = EPS));
}

#[test]
fn test_build_invalid_coefficients() {
    let poly1 = Polynomial::build(&vec![1.0, 2.0, std::f64::NAN]);
    let poly2 = Polynomial::build(&vec![1.0, 2.0, std::f64::INFINITY]);

    assert!(matches!(poly1.unwrap_err(), PolyError::InvalidCoefficients));
    assert!(matches!(poly2.unwrap_err(), PolyError::InvalidCoefficients));
}

// ===============================================================================================

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_complex_eval1() {
    let coef = Complex64::new(0.3, 0.0);
    let p = Polynomial::build(&vec![coef]).unwrap();
    let x = Complex64::new(0.75, 1.2);
    let y = p.eval(x);

    assert!(is_close!(y.re, 0.3, rel_tol = EPS));
    assert!(is_close!(y.im, 0.0, rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_complex_eval2() {
    let coefs = &vec![
        Complex64::new(2.1, 0.0),
        Complex64::new(-1.34, 0.0),
        Complex64::new(0.76, 0.0),
        Complex64::new(0.45, 0.0),
    ];
    let p = Polynomial::build(coefs).unwrap();
    let x = Complex64::new(0.49, 0.95);
    let y = p.eval(x);

    assert!(is_close!(y.re, 0.3959143, rel_tol = EPS));
    assert!(is_close!(y.im, -0.6433305, rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_complex_eval3() {
    let coef = &vec![Complex64::new(0.674, -1.423)];
    let p = Polynomial::build(coef).unwrap();
    let x = Complex64::new(-1.44, 9.55);
    let y = p.eval(x);

    assert!(is_close!(y.re, 0.674, rel_tol = EPS));
    assert!(is_close!(y.im, -1.423, rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_complex_eval4() {
    let coefs = &vec![
        Complex64::new(-2.31, 0.44),
        Complex64::new(4.21, -3.19),
        Complex64::new(0.93, 1.04),
        Complex64::new(-0.42, 0.68),
    ];
    let p = Polynomial::build(coefs).unwrap();
    let x = Complex64::new(0.49, 0.95);
    let y = p.eval(x);

    assert!(is_close!(y.re, 1.82462012, rel_tol = EPS));
    assert!(is_close!(y.im, 2.30389412, rel_tol = EPS));
}
