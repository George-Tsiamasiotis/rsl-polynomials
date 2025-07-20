//! A re-write of [`GSL's Polynomial Routines`].
//!
//! [`GSL's Polynomial Routines`]: https://www.gnu.org/software/gsl/doc/html/poly.html

mod polynomial;

#[cfg(test)]
mod test;

pub use polynomial::Polynomial;
