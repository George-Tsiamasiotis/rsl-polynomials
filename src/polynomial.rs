//! Methods for evaluating a polynomial and its derivatives on a certain point.

use num::{ToPrimitive, Zero};

use crate::{PolyError, Result};

#[allow(rustdoc::broken_intra_doc_links)]
/// Representation of a polynomial.
///
/// A polynomial of degree `n`, represented with a [`Vec`] of length `n+1` containing the coefficients
/// `c[i]`:
///
/// P(x) = c[0] + c[1]x + c[2]x² + ... + c[n−1]xⁿ⁻¹ + c[n]xⁿ
///
/// [`Vec`]: std::vec::Vec
#[derive(Clone)]
pub struct Polynomial<T> {
    /// The polynomial's coefficients.
    pub coef: Vec<T>,
}

impl<T> Polynomial<T>
where
    T: num::complex::ComplexFloat,
{
    /// Creates a new Polynomial with no terms (zero polynomial).
    pub fn new() -> Self {
        Polynomial {
            coef: vec![T::zero()],
        }
    }

    /// Creates a new Polynomial from the given coefficients.
    ///
    /// ## Example
    ///
    /// ```
    /// # use rsl_polynomials::{Polynomial, Result};
    /// # fn main() -> Result<()> {
    /// let poly = Polynomial::build(&vec![1.0, 4.0, 3.0])?; // 1+4x+3x²
    /// # Ok(())
    /// # }
    /// ```
    pub fn build(coef: &[T]) -> Result<Self> {
        if coef.len().is_zero() {
            return Ok(Polynomial::new());
        }

        match coef.iter().any(|x| x.is_nan() | x.is_infinite()) {
            true => Err(PolyError::InvalidCoefficients),
            false => Ok(Polynomial {
                coef: coef.to_vec(),
            }),
        }
    }

    /// Trims the higher order terms with 0 coefficient.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_polynomials::{Polynomial, Result};
    /// # fn main() -> Result<()> {
    /// // 0+x+0+2x³+0 −> x+2x³
    /// let poly1 = Polynomial::build(&vec![0.0, 1.0, 0.0, 2.0, 0.0])?.to_trimmed();
    ///
    /// assert_eq!(poly1.coef, vec![0.0, 1.0, 0.0, 2.0]);
    /// # Ok(())
    /// # }
    /// ```
    pub fn to_trimmed(&self) -> Self {
        // Leave [0.0] polynomial as is
        if self.coef.len() == 1 {
            return self.clone();
        }

        let mut iter = self.coef.iter().rev().copied().peekable();
        let mut new_coeffs = self.coef.clone();

        new_coeffs.reverse();
        while iter.peek().is_some_and(|c| c.is_zero()) {
            new_coeffs.remove(0);
            iter.next();
        }
        new_coeffs.reverse();

        Polynomial { coef: new_coeffs }
    }

    /// Evaluates the polynomial for the value `x`.
    ///
    /// ## Example
    ///
    /// ```
    /// # use rsl_polynomials::{Polynomial, Result};
    /// # fn main() -> Result<()> {
    /// let poly = Polynomial::build(&vec![1.0, 2.0, 3.0])?;
    ///
    /// assert_eq!(poly.eval(1.0), 6.0);
    /// assert_eq!(poly.eval(-1.0), 2.0);
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval(&self, x: T) -> T {
        // NOTE: This evaluates a₀+a₁x+a₂x²+...+aₙx² as if it were in the form
        // a₀+x(a₁+x(a₂+ x(...))), therefore saving a lot of reduntant multiplications.

        // NOTE: `gsl_poly_complex_eval()` and `gsl_complex_poly_complex_eval()` seem to do the
        // same thing as `gsl_poly_eval()`, but perform the complex addition and multiplication
        // manually since its slightly faster.

        self.coef
            .iter()
            .rev()
            .copied()
            .reduce(|res, coef| coef + x * res)
            .unwrap_or(T::zero())
    }

    /// Evaluates the polynomials first `n` derivatives (including the 0-th derivative, i.e. the
    /// polynomial's value) for the value `x`.
    ///
    /// The result is a vector holding the calculated derivatives:
    ///
    /// [d⁰/dx⁰, d¹/dx¹, d²/dx², ..., dⁿ/dxⁿ]
    ///
    /// ## Example
    ///
    /// ```
    /// # use rsl_polynomials::{Polynomial, Result};
    /// # fn main() -> Result<()> {
    /// let poly = Polynomial::build(&vec![1.0, 2.0, 3.0])?;
    ///
    /// assert_eq!(poly.eval_derivs(1.0, 4), vec![6.0, 8.0, 6.0, 0.0]);
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval_derivs(&self, x: T, n: usize) -> Vec<T> {
        let mut res: Vec<T> = vec![T::zero(); n];

        let last_idx = self.coef.len() - 1;
        let nmax = self.coef.len().min(res.len()) - 1;

        // Partially fill res with the dominant term's coefficient
        res.iter_mut()
            .take(nmax + 1)
            .for_each(|e| *e = *self.coef.last().unwrap());

        for i in 0..last_idx {
            let k = last_idx - i;
            res[0] = x * res[0] + self.coef[k - 1];
            let jmax = if nmax < k { nmax } else { k - 1 };
            for j in 1..=jmax {
                res[j] = x * res[j] + res[j - 1];
            }
        }

        // Mutliply each term by the corresponding exponents
        let mut f = T::one();
        for (i, d) in res.iter_mut().enumerate().take(nmax + 1).skip(2) {
            f = f * T::from(i).unwrap();
            *d = *d * f;
        }

        res
    }

    /// Calculates the **real** roots af a quadratic equation `ax²+bx+c`.
    ///
    /// # Error
    ///
    /// Returns an error in 3 cases:
    /// 1. the Polynomial is not of order 2
    /// 2. one of the coefficients is not real
    /// 3. the Polynomial is constant, i.e. a=b=0
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_polynomials::{Polynomial, Result};
    /// #
    /// # fn main() -> Result<()> {
    /// let poly = Polynomial::build(&vec![-20.0, 0.0, 5.0])?; // 5x²-20
    /// let y = poly.solve_real_quadratic()?;
    /// let expected = vec![2.0, -2.0];
    ///
    /// assert_eq!(y, expected);
    /// # Ok(())
    /// # }
    /// ```
    #[doc(alias = "gsl_poly_solve_quadratic")]
    pub fn solve_real_quadratic(&self) -> Result<Vec<f64>> {
        crate::utils::check_if_correct_order(&self.coef, 2)?;
        crate::utils::check_if_real_coefficients(&self.coef)?;

        let a = self.coef[2].re().to_f64().expect("Error converting to f64");
        let b = self.coef[1].re().to_f64().expect("Error converting to f64");
        let c = self.coef[0].re().to_f64().expect("Error converting to f64");

        crate::solve::solve_real_quadratic(a, b, c)
    }
}

impl<T> Default for Polynomial<T>
where
    T: num::complex::ComplexFloat,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T> std::fmt::Debug for Polynomial<T>
where
    T: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Polynomial")
            .field("Coefficients", &self.coef)
            .finish()
    }
}
