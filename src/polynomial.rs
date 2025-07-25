//! Methods for evaluating a polynomial and its derivatives on a certain point.

use num::Zero;

use crate::{
    PolyError, Result, solve,
    utils::{check_if_correct_order, check_if_real_coefficients, convert_complex_to_real},
};

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
pub struct Polynomial<T>
where
    T: std::fmt::Debug,
{
    /// The polynomial's coefficients.
    pub coef: Vec<T>,
}

impl<T> Polynomial<T>
where
    T: num::complex::ComplexFloat + std::fmt::Debug,
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

    /// Converts a general polynomial to a monic polynomial:
    /// ax³ + bx² + cx + d  −>  x³ + ax² + bx + c
    ///
    /// ## Example
    ///
    /// ```
    /// # use rsl_polynomials::{Polynomial, Result};
    /// # fn main() -> Result<()> {
    /// let poly = Polynomial::build(&vec![30.0, 6.0, 3.0])?.to_monic();
    ///
    /// assert_eq!(poly.coef, vec![10.0, 2.0, 1.0]);
    /// # Ok(())
    /// # }
    /// ```
    pub fn to_monic(&self) -> Self {
        // Leave [0.0] polynomial as is
        if self.coef.len() == 1 {
            return self.clone();
        }

        let mut monic = self.to_trimmed();
        let a = *monic.coef.last().unwrap();
        monic.coef.iter_mut().for_each(|e| *e = *e / a);
        monic
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
    #[doc(alias = "gsl_poly_eval")]
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
    #[doc(alias = "gsl_poly_eval_derivs")]
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
        check_if_correct_order(&self.coef, 2)?;
        check_if_real_coefficients(&self.coef)?;

        let mut reals = Vec::<f64>::new();
        for c in self.coef.iter() {
            reals.push(convert_complex_to_real(*c)?);
        }

        solve::solve_real_quadratic(reals[2], reals[1], reals[0])
    }

    /// Calculates the **real** roots af a quadratic equation `ax³+bx²+cx+d`.
    ///
    /// # Error
    ///
    /// Returns an error in 3 cases:
    /// 1. the Polynomial is not of order 3
    /// 2. one of the coefficients is not real
    /// 3. the Polynomial is constant, i.e. a=b=c=0
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_polynomials::{Polynomial, Result};
    /// #
    /// # fn main() -> Result<()> {
    /// let poly = Polynomial::build(&vec![-6.0, 11.0, -6.0, 1.0])?; // x³-6x²+11x-6
    /// let y = poly.solve_real_cubic()?;
    /// let expected = vec![2.0, -2.0]; // TODO:
    ///
    /// assert_eq!(y, expected);
    /// # Ok(())
    /// # }
    /// ```
    #[doc(alias = "gsl_poly_solve_cubic")]
    pub fn solve_real_cubic(&self) -> Result<Vec<f64>> {
        check_if_correct_order(&self.coef, 3)?;
        check_if_real_coefficients(&self.coef)?;

        let monic = self.to_monic();

        let mut reals = Vec::<f64>::new();
        for c in monic.coef.iter().skip(1) {
            reals.push(convert_complex_to_real(*c)?);
        }

        solve::solve_real_cubic(reals[2], reals[1], reals[0])
    }
}

impl<T> Default for Polynomial<T>
where
    T: num::complex::ComplexFloat + std::fmt::Debug,
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
