use crate::Polynomial;
use is_close::is_close;

// GSL's tests use this tolerance
const EPS: f64 = 100.0 * f64::EPSILON;

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_cubic1() {
    let p = Polynomial::build(&[-27.0, 0.0, 0.0, 1.0]).unwrap();
    let y = p.solve_real_cubic().unwrap();
    let expected = [3.0];

    assert_eq!(y.len(), 3);
    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_cubic2() {
    let p = Polynomial::build(&[-4913.0, 867.0, -51.0, 1.0]).unwrap();
    let y = p.solve_real_cubic().unwrap();
    let expected = [17.0, 17.0, 17.0];

    assert_eq!(y.len(), 3);
    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
    assert!(is_close!(y[1], expected[1], rel_tol = EPS));
    assert!(is_close!(y[2], expected[2], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_cubic3() {
    let p = Polynomial::build(&[-6647.0, 1071.0, -57.0, 1.0]).unwrap();
    let y = p.solve_real_cubic().unwrap();
    let expected = [17.0, 17.0, 23.0];

    assert_eq!(y.len(), 3);
    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
    assert!(is_close!(y[1], expected[1], rel_tol = EPS));
    assert!(is_close!(y[2], expected[2], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_cubic4() {
    let p = Polynomial::build(&[6647.0, -493.0, -11.0, 1.0]).unwrap();
    let y = p.solve_real_cubic().unwrap();
    let expected = [-23.0, 17.0, 17.0];

    assert_eq!(y.len(), 3);
    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
    assert!(is_close!(y[1], expected[1], rel_tol = EPS));
    assert!(is_close!(y[2], expected[2], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_cubic5() {
    let p = Polynomial::build(&[-50065.0, 5087.0, -143.0, 1.0]).unwrap();
    let y = p.solve_real_cubic().unwrap();
    let expected = [17.0, 31.0, 95.0];

    assert_eq!(y.len(), 3);
    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
    assert!(is_close!(y[1], expected[1], rel_tol = EPS));
    assert!(is_close!(y[2], expected[2], rel_tol = EPS));
}

#[test]
/// Source: gsl/poly/test.c
fn test_gsl_cubic6() {
    let p = Polynomial::build(&[50065.0, 803.0, -109.0, 1.0]).unwrap();
    let y = p.solve_real_cubic().unwrap();
    let expected = [-17.0, 31.0, 95.0];

    assert_eq!(y.len(), 3);
    assert!(is_close!(y[0], expected[0], rel_tol = EPS));
    assert!(is_close!(y[1], expected[1], rel_tol = EPS));
    assert!(is_close!(y[2], expected[2], rel_tol = EPS));
}
