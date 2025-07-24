use num::complex::Complex64;

use crate::{PolyError, Polynomial};

#[test]
fn test_build_polynomial() {
    assert!(Polynomial::build(&Vec::<f64>::new()).is_ok());
    assert!(Polynomial::build(&vec![1.0]).is_ok());
    assert!(Polynomial::build(&vec![1.0, 2.0, 3.0]).is_ok());
}

#[test]
fn test_new_polynomial() {
    let float_poly = Polynomial::<f64>::new();
    let complex_poly = Polynomial::<Complex64>::new();

    assert_eq!(float_poly.coef, vec![0.0]);
    assert_eq!(complex_poly.coef, vec![Complex64::new(0.0, 0.0)]);
}

#[test]
fn test_default_polynomial() {
    let default_float_poly = Polynomial::<f64>::default();
    let default_complex_poly = Polynomial::<Complex64>::default();

    assert_eq!(default_float_poly.coef, vec![0.0]);
    assert_eq!(default_complex_poly.coef, vec![Complex64::new(0.0, 0.0)]);
}

#[test]
fn test_build_polynomial_invalid() {
    let poly1 = Polynomial::build(&vec![1.0, 2.0, std::f64::NAN]);
    let poly2 = Polynomial::build(&vec![1.0, 2.0, std::f64::INFINITY]);

    assert!(matches!(poly1.unwrap_err(), PolyError::InvalidCoefficients));
    assert!(matches!(poly2.unwrap_err(), PolyError::InvalidCoefficients));
}

#[test]
#[rustfmt::skip]
fn test_trim_trailing_zeros() {
    let poly0 = Polynomial::build(&vec![0.0]).unwrap().to_trimmed();
    let poly1 = Polynomial::build(&vec![0.0, 1.0, 2.0]).unwrap().to_trimmed();
    let poly2 = Polynomial::build(&vec![0.0, 1.0, 2.0, 0.0, 0.0]).unwrap().to_trimmed();
    let poly3 = Polynomial::build(&vec![1.0, 2.0]).unwrap().to_trimmed();
    let poly4 = Polynomial::build(&vec![1.0, 2.0, 0.0, 0.0]).unwrap().to_trimmed();
    let poly5 = Polynomial::build(&vec![1.0, 0.0, 2.0]).unwrap().to_trimmed();

    assert_eq!(poly0.coef, vec![0.0]);
    assert_eq!(poly1.coef, vec![0.0, 1.0, 2.0]);
    assert_eq!(poly2.coef, vec![0.0, 1.0, 2.0]);
    assert_eq!(poly3.coef, vec![1.0, 2.0]);
    assert_eq!(poly4.coef, vec![1.0, 2.0]);
    assert_eq!(poly5.coef, vec![1.0, 0.0, 2.0]);
}

#[test]
fn test_debug() {
    let poly = Polynomial::build(&vec![0.0]).unwrap();
    let _ = format!("{:?}", poly);
    let _ = format!("{:#?}", poly);
}
