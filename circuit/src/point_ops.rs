//! Circuit-friendly elliptic curve point operations.
//!
//! This module provides elliptic curve operations optimized for circuit use,
//! where all operations are performed using native KoalaBear field arithmetic.

use curve::{Affine, BaseField, KoalaBear};
use p3_field::PrimeCharacteristicRing;

/// Elliptic curve point in circuit representation.
///
/// Points are represented in affine coordinates (x, y) where each coordinate
/// is a degree-8 extension field element over KoalaBear.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CircuitPoint {
    /// X-coordinate as 8 KoalaBear coefficients
    pub x: [KoalaBear; 8],
    /// Y-coordinate as 8 KoalaBear coefficients
    pub y: [KoalaBear; 8],
    /// Flag indicating if this is the point at infinity
    pub is_infinity: bool,
}

impl CircuitPoint {
    /// Create the point at infinity
    pub fn infinity() -> Self {
        Self {
            x: [KoalaBear::ZERO; 8],
            y: [KoalaBear::ZERO; 8],
            is_infinity: true,
        }
    }

    /// Create a circuit point from an affine point
    pub fn from_affine(point: &Affine) -> Self {
        if point.is_infinity() {
            return Self::infinity();
        }

        let x_coeffs: [KoalaBear; 8] = unsafe { core::mem::transmute(point.x) };
        let y_coeffs: [KoalaBear; 8] = unsafe { core::mem::transmute(point.y) };

        Self {
            x: x_coeffs,
            y: y_coeffs,
            is_infinity: false,
        }
    }

    /// Convert back to an affine point
    pub fn to_affine(&self) -> Affine {
        if self.is_infinity {
            return Affine::INFINITY;
        }

        let x: BaseField = unsafe { core::mem::transmute(self.x) };
        let y: BaseField = unsafe { core::mem::transmute(self.y) };

        Affine::new(x, y)
    }

    /// Get the x-coordinate as a BaseField element
    pub fn x_as_basefield(&self) -> BaseField {
        unsafe { core::mem::transmute(self.x) }
    }

    /// Get the y-coordinate as a BaseField element
    pub fn y_as_basefield(&self) -> BaseField {
        unsafe { core::mem::transmute(self.y) }
    }

    /// Point addition in the circuit
    ///
    /// Uses the underlying Affine point addition for correctness
    pub fn add(&self, other: &Self) -> Self {
        let p1 = self.to_affine();
        let p2 = other.to_affine();
        let result = p1 + p2;
        Self::from_affine(&result)
    }

    /// Point doubling in the circuit
    ///
    /// Uses the underlying Affine point doubling for correctness
    pub fn double(&self) -> Self {
        let p = self.to_affine();
        let result = p.double();
        Self::from_affine(&result)
    }

    /// Negate a point
    pub fn negate(&self) -> Self {
        let p = self.to_affine();
        let result = p.negate();
        Self::from_affine(&result)
    }

    /// Scalar multiplication using double-and-add algorithm
    ///
    /// This is circuit-friendly as it uses only point additions and doublings.
    /// For ZK circuits, you would typically use a windowed or fixed-base version.
    pub fn scalar_mul(&self, scalar_bits: &[bool]) -> Self {
        let mut result = Self::infinity();
        let mut temp = self.clone();

        for &bit in scalar_bits {
            if bit {
                result = result.add(&temp);
            }
            temp = temp.double();
        }

        result
    }

    /// Check if two points are equal
    pub fn equals(&self, other: &Self) -> bool {
        if self.is_infinity && other.is_infinity {
            return true;
        }
        if self.is_infinity || other.is_infinity {
            return false;
        }

        self.x == other.x && self.y == other.y
    }
}

/// Convert a scalar field element to bit representation for scalar multiplication
pub fn scalar_to_bits(scalar: &curve::ScalarField) -> Vec<bool> {
    let limbs = scalar.to_canonical_u64_vec();
    let mut bits = Vec::with_capacity(256);

    for limb in &limbs {
        for i in 0..64 {
            bits.push((*limb >> i) & 1 == 1);
        }
    }

    bits
}

#[cfg(test)]
mod tests {
    use super::*;
    use curve::{Projective, ScalarField};

    #[test]
    fn test_point_conversion() {
        let point = Projective::generator();
        let affine_point = Affine::from(&point);

        let circuit_point = CircuitPoint::from_affine(&affine_point);
        let recovered = circuit_point.to_affine();

        assert_eq!(affine_point, recovered);
    }

    #[test]
    fn test_point_addition() {
        // Use generator and its multiples
        let p1 = Affine::generator();
        let p2_proj = Projective::generator() * ScalarField::from_canonical_u64(2);
        let p2 = Affine::from(&p2_proj);

        let cp1 = CircuitPoint::from_affine(&p1);
        let cp2 = CircuitPoint::from_affine(&p2);

        let circuit_sum = cp1.add(&cp2);
        let expected_sum = Affine::from(&(Projective::from(&p1) + Projective::from(&p2)));

        assert_eq!(circuit_sum.to_affine(), expected_sum);
    }

    #[test]
    fn test_point_doubling() {
        let point = Affine::generator();

        let cp = CircuitPoint::from_affine(&point);
        let doubled = cp.double();
        let expected = Affine::from(&(Projective::from(&point) + Projective::from(&point)));

        assert_eq!(doubled.to_affine(), expected);
    }

    #[test]
    fn test_scalar_multiplication() {
        let point = Affine::generator();
        let scalar = ScalarField::from_canonical_u64(17);

        let cp = CircuitPoint::from_affine(&point);
        let bits = scalar_to_bits(&scalar);
        let result = cp.scalar_mul(&bits);

        let expected = Affine::from(&(Projective::from(&point) * scalar));

        assert_eq!(result.to_affine(), expected);
    }
}
