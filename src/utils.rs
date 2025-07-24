use crate::{PolyError, Result};
use num::ToPrimitive;

/// Trails all trailing zeros of a Vec.
pub(crate) fn trim_trailing_zeros<T>(values: &Vec<T>) -> Vec<T>
where
    T: num::complex::ComplexFloat,
{
    // Leave [0.0] polynomial as is
    if values.len() == 1 {
        return values.to_owned();
    }

    // Trim leading zeros, then reverse
    let mut keep_trimming = true;
    let mut filtered_values: Vec<T> = values
        .clone()
        .into_iter()
        .rev()
        .filter(|e| {
            if e.is_zero() & keep_trimming {
                return false;
            };
            keep_trimming = false;
            true
        })
        .collect();

    filtered_values.reverse();
    filtered_values
}

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

#[cfg(test)]
mod test {
    use num::complex::Complex64;

    use crate::Polynomial;

    use super::*;

    #[test]
    fn test_correct_order() {
        let poly1 = Polynomial::build(&vec![1.0]).unwrap();
        let poly2 = Polynomial::build(&vec![1.0, 2.0]).unwrap();
        let poly3 = Polynomial::build(&vec![1.0, 2.0, 3.0]).unwrap();

        assert!(check_if_correct_order(&poly1.coef, 0).is_ok());
        assert!(check_if_correct_order(&poly2.coef, 1).is_ok());
        assert!(check_if_correct_order(&poly3.coef, 2).is_ok());
        assert!(check_if_correct_order(&poly3.coef, 3).is_err());
    }

    #[test]
    fn test_real_coefficients() {
        let poly1 = Polynomial::build(&vec![1.0]).unwrap();
        let poly2 = Polynomial::build(&vec![1.0, 2.0]).unwrap();
        let poly3 = Polynomial::build(&vec![Complex64::new(1.0, 1.0)]).unwrap();
        let poly4 =
            Polynomial::build(&vec![Complex64::new(1.0, 1.0), Complex64::new(2.0, 2.0)]).unwrap();

        assert!(check_if_real_coefficients(&poly1.coef).is_ok());
        assert!(check_if_real_coefficients(&poly2.coef).is_ok());
        assert!(check_if_real_coefficients(&poly3.coef).is_err());
        assert!(check_if_real_coefficients(&poly4.coef).is_err());
    }
}
