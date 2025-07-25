use crate::{PolyError, Result};
use num::{ToPrimitive, Zero};

/// Checks if a polynomial is of the expected order.
pub(crate) fn check_if_correct_order<T>(coef: &[T], expected_order: usize) -> Result<()> {
    if coef.len() != expected_order + 1 {
        return Err(PolyError::IncorrectOrder(expected_order));
    }
    Ok(())
}

/// Checks if all the coefficients of a poly are real, i.e. their imaginary part is 0.
pub(crate) fn check_if_real_coefficients<C: num::complex::ComplexFloat>(coef: &[C]) -> Result<()> {
    for c in coef.iter() {
        let cf = match c.im().to_f64() {
            Some(cf) => cf,
            None => unreachable!("Could not convert imaginary part of ComplexFloat to f64"),
        };
        match cf {
            0.0 => (),
            _ => return Err(PolyError::NotRealCoefficients),
        }
    }
    Ok(())
}

/// Converts a Complex number to f64. Returns an Error if the complex number has an imaginary part.
pub(crate) fn convert_complex_to_real<C>(number: C) -> Result<f64>
where
    C: num::complex::ComplexFloat + std::fmt::Debug,
{
    let err = PolyError::ComplexTof64Conversion(format!("{number:?}").into());

    // Complex64.to_f64() returns the real part, even if the imaginary part is not 0.
    if !number.is_finite() | !number.im().is_zero() {
        return Err(err);
    }

    number.re().to_f64().ok_or(err)
}

#[cfg(test)]
mod test {
    use num::complex::Complex64;

    use crate::Polynomial;

    use super::*;

    #[test]
    fn test_correct_order() {
        let poly1 = Polynomial::build(&[1.0]).unwrap();
        let poly2 = Polynomial::build(&[1.0, 2.0]).unwrap();
        let poly3 = Polynomial::build(&[1.0, 2.0, 3.0]).unwrap();

        assert!(check_if_correct_order(&poly1.coef, 0).is_ok());
        assert!(check_if_correct_order(&poly2.coef, 1).is_ok());
        assert!(check_if_correct_order(&poly3.coef, 2).is_ok());
        assert!(check_if_correct_order(&poly3.coef, 3).is_err());
    }

    #[test]
    fn test_real_coefficients() {
        let poly1 = Polynomial::build(&[1.0]).unwrap();
        let poly2 = Polynomial::build(&[1.0, 2.0]).unwrap();
        let poly3 = Polynomial::build(&[Complex64::new(1.0, 1.0)]).unwrap();
        let poly4 =
            Polynomial::build(&[Complex64::new(1.0, 1.0), Complex64::new(2.0, 2.0)]).unwrap();

        assert!(check_if_real_coefficients(&poly1.coef).is_ok());
        assert!(check_if_real_coefficients(&poly2.coef).is_ok());
        assert!(check_if_real_coefficients(&poly3.coef).is_err());
        assert!(check_if_real_coefficients(&poly4.coef).is_err());
    }

    #[test]
    fn test_complex_to_f64_conversion() {
        let c1 = Complex64::new(1.0, 0.0);
        let c2 = Complex64::new(1.0, 1.0);

        assert_eq!(convert_complex_to_real(c1).unwrap(), 1.0f64);
        assert!(matches!(
            convert_complex_to_real(c2).unwrap_err(),
            PolyError::ComplexTof64Conversion(_)
        ))
    }
}
