use rand::distr::{Distribution, StandardUniform};
use rand::Rng;

use crate::{BaseField, ScalarField};

/// Helper trait for sampling random field elements.
pub trait RandomField: Sized {
    fn random<R: Rng + ?Sized>(rng: &mut R) -> Self;
}

impl RandomField for BaseField {
    #[inline]
    fn random<R: Rng + ?Sized>(rng: &mut R) -> Self {
        StandardUniform.sample(rng)
    }
}

impl RandomField for ScalarField {
    #[inline]
    fn random<R: Rng + ?Sized>(rng: &mut R) -> Self {
        StandardUniform.sample(rng)
    }
}
