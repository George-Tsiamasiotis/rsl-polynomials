use crate::{PolyError, Polynomial};
use num::complex::Complex64;

#[test]
fn test_solve_real_quadratic_wrong_order() {
    let p = Polynomial::build(&vec![1.0, 2.0, 3.0, 4.0]).unwrap();

    assert!(matches!(
        p.solve_real_quadratic().unwrap_err(),
        PolyError::IncorrectOrder(2)
    ));
}

#[test]
fn test_solve_real_quadratic_complex_coefs() {
    let p = Polynomial::build(&vec![
        Complex64::new(1.0, 2.0),
        Complex64::new(3.0, 4.0),
        Complex64::new(5.0, 6.0),
    ])
    .unwrap();

    assert!(matches!(
        p.solve_real_quadratic().unwrap_err(),
        PolyError::NotRealCoefficients
    ));
}
