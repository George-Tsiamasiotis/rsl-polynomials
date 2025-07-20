//! Methods for evaluating a polynomial and its derivatives on a certain point.

#[allow(rustdoc::broken_intra_doc_links)]
/// Representation of a polynomial.
///
/// A polynomial of degree `n`, represented with a [`Vec`] of length `n+1` containing the coefficients
/// `c[i]`:
///
/// P(x) = c[0] + c[1]x + c[2]x² + ... + c[n−1]xⁿ⁻¹ + c[n]xⁿ
///
/// [`Vec`]: std::vec::Vec
pub struct Polynomial<T> {
    /// The polynomial's coefficients.
    pub coef: Vec<T>,
    /// The order of the polynomial.
    pub order: usize,
}

impl<T: num::complex::ComplexFloat> Polynomial<T> {
    /// Creates a new Polynomial from the given coefficients.
    ///
    /// ## Example
    ///
    /// ```
    /// # use rsl_polynomials::Polynomial;
    /// # fn main() {
    /// let poly = Polynomial::new(vec![1.0, 0.4, 3.0]);
    /// # }
    /// ```
    pub fn new(coef: Vec<T>) -> Self {
        let order = coef.len() - 1;
        Polynomial { coef, order }
    }

    /// Evaluates the polynomial for the value `x`.
    ///
    /// ## Example
    ///
    /// ```
    /// # use rsl_polynomials::Polynomial;
    /// # fn main() {
    /// let poly = Polynomial::new(vec![1.0, 2.0, 3.0]);
    ///
    /// assert_eq!(poly.eval(1.0), 6.0);
    /// assert_eq!(poly.eval(-1.0), 2.0);
    /// # }
    /// ```
    pub fn eval(&self, x: T) -> T {
        let last_idx = self.order;
        let mut res = self.coef[last_idx];

        // NOTE: `gsl_poly_complex_eval()` and `gsl_complex_poly_complex_eval()` seem to do the
        // same thing as `gsl_poly_eval()`, but perform the complex addition and multiplication
        // manually since its faster. I don't think this would make any difference here though.
        for i in (1..=last_idx).rev() {
            res = self.coef[i - 1] + x * res
        }
        res
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
    /// # use rsl_polynomials::Polynomial;
    /// # fn main() {
    /// let poly = Polynomial::new(vec![1.0, 2.0, 3.0]);
    ///
    /// assert_eq!(poly.eval_derivs(1.0, 4), vec![6.0, 8.0, 6.0, 0.0]);
    /// # }
    /// ```
    pub fn eval_derivs(&self, x: T, n: usize) -> Vec<T> {
        let mut res: Vec<T> = vec![T::zero(); n];

        let lenres = res.len();
        let lenpoly = self.coef.len(); // lenc
        let last_idx = self.order;
        let nmax = lenpoly.min(lenres) - 1;

        // Partially fill res with the dominant term's coefficient
        for d in res.iter_mut().take(nmax + 1) {
            *d = self.coef[last_idx]
        }

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
}
