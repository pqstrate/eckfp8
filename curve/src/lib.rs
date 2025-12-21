//! Elliptic curve group over the KoalaBear degree-8 extension field.
//!
//! This crate provides affine and projective curve points, a scalar field
//! implementation, and helpers for random sampling. The curve parameters and
//! generators are fixed to the values in the `affine` module.

mod affine;
mod basefield;
mod group;
mod projective;
mod random;
mod scalarfield;

pub use affine::Affine;
pub use basefield::BaseField;
pub use group::{Group, ScalarBits};
pub use p3_koala_bear::KoalaBear;
pub use projective::Projective;
pub use random::RandomField;
pub use scalarfield::ScalarField;
