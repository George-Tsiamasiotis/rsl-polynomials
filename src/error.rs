#[derive(thiserror::Error, std::fmt::Debug)]
/// The error type for Polynomial operations.
pub enum PolyError {
    /// Supplied coefficients are NaN or Infinity
    #[error("Supplied coefficients cannot be NaN or Infinity")]
    InvalidCoefficients,

    /// Supplied Polynomial is trivial.
    #[error("Supplied Polynomial is trivial.")]
    Trivial,

    /// Supplied Polynomial has incorrect order.
    #[error("Supplied Polynomial must be of order {0}")]
    IncorrectOrder(usize),

    /// Supplied Polynomial is constant.
    #[error("Supplied Polynomial is constant.")]
    ConstantPoly,

    /// Supplied Polynomial is not quadratic.
    #[error("Supplied Polynomial is not quadratic: {0}")]
    NotQuadratic(Box<str>),

    /// Supplied Polynomial has no real roots.
    #[error("Supplied Polynomial has no real roots.")]
    NoRealRoots,

    /// Called solve_real() on a poly with complex coefficients.
    #[error("Supplied Polynomial must have real coefficients.")]
    NotRealCoefficients,

    /// Discriminant calculation returned NaN.
    #[error("Discriminant calculation returned NaN.")]
    NanDiscriminant,
}
