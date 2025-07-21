//! A re-write of [`GSL's Polynomial Routines`].
//!
//! [`GSL's Polynomial Routines`]: https://www.gnu.org/software/gsl/doc/html/poly.html

mod error;
mod polynomial;
mod solve;
mod utils;

#[cfg(test)]
mod test;

pub use error::PolyError;
pub use polynomial::Polynomial;

pub type Result<T> = std::result::Result<T, error::PolyError>;
