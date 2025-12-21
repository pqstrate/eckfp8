use crate::affine::Affine;
use crate::basefield::{from_coeffs, BaseField};
use crate::scalarfield::ScalarField;
use core::ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign};
use p3_field::{Field, PrimeCharacteristicRing};
use p3_koala_bear::KoalaBear;
use serde::{Deserialize, Serialize};

/// Projective point on the elliptic curve
/// Represents a point in projective coordinates (X:Y:Z) where (x,y) = (X/Z, Y/Z)
/// The point at infinity is represented as (0:1:0)
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Projective {
    pub x: BaseField,
    pub y: BaseField,
    pub z: BaseField,
}

impl Projective {
    // Curve parameters: y^2 = x^3 + a*x + b
    // a = 3*u where u is the primitive element of the extension
    // b = 42639

    /// Get the 'a' coefficient: 3*u
    #[inline]
    fn curve_a() -> BaseField {
        // Create 3*u: coefficients [0, 3, 0, 0, 0, 0, 0, 0]
        // The second component (index 1) represents u^1
        let zero = KoalaBear::ZERO;
        let three = KoalaBear::new(3);
        from_coeffs([zero, three, zero, zero, zero, zero, zero, zero])
    }

    /// Get the 'b' coefficient: 42639
    #[inline]
    fn curve_b() -> BaseField {
        // Create constant 42639 in the extension field
        let zero = KoalaBear::ZERO;
        let b = KoalaBear::new(42639);
        from_coeffs([b, zero, zero, zero, zero, zero, zero, zero])
    }

    /// The point at infinity (identity element): (0:1:0)
    pub const INFINITY: Self = Projective {
        x: BaseField::ZERO,
        y: BaseField::ONE,
        z: BaseField::ZERO,
    };

    /// Create a new projective point
    pub fn new(x: BaseField, y: BaseField, z: BaseField) -> Self {
        Projective { x, y, z }
    }

    /// Check if this point is the point at infinity
    #[inline]
    pub fn is_infinity(&self) -> bool {
        self.z.is_zero()
    }

    /// Convert to affine coordinates
    pub fn to_affine(&self) -> Affine {
        if self.is_infinity() {
            return Affine::INFINITY;
        }

        let z_inv = self.z.inverse();
        let x = self.x * z_inv;
        let y = self.y * z_inv;

        Affine::new(x, y)
    }

    /// Convert from affine coordinates
    pub fn from_affine(point: &Affine) -> Self {
        if point.is_infinity() {
            return Self::INFINITY;
        }

        Projective::new(point.x, point.y, BaseField::ONE)
    }

    /// Check if a point is on the curve: Y^2*Z = X^3 + a*X*Z^2 + b*Z^3
    pub fn is_on_curve(&self) -> bool {
        if self.is_infinity() {
            return true;
        }

        let y2 = self.y * self.y;
        let x2 = self.x * self.x;
        let x3 = x2 * self.x;
        let z2 = self.z * self.z;
        let z3 = z2 * self.z;

        let lhs = y2 * self.z;
        let rhs = x3 + Self::curve_a() * self.x * z2 + Self::curve_b() * z3;

        lhs == rhs
    }

    /// Generator point from SSWU on 'ZKM2'
    pub fn generator() -> Self {
        Self::from_affine(&Affine::generator())
    }

    /// Alternative generator point from SSWU on 'ZKM2 - Pedersen'
    pub fn generator_pedersen() -> Self {
        Self::from_affine(&Affine::generator_pedersen())
    }

    /// Point doubling: 2*P using projective coordinates
    /// For simplicity and correctness, convert to affine, double, and convert back
    pub fn double(&self) -> Self {
        if self.is_infinity() {
            return *self;
        }

        let affine = self.to_affine();
        let doubled_affine = affine.double();
        Self::from_affine(&doubled_affine)
    }

    /// Negate a point
    pub fn negate(&self) -> Self {
        if self.is_infinity() {
            return *self;
        }
        Projective::new(self.x, -self.y, self.z)
    }

    /// Scalar multiplication using double-and-add algorithm
    pub fn scalar_mul(&self, scalar: &ScalarField) -> Self {
        let scalar_bytes = scalar.to_canonical_u64_vec();
        let mut result = Self::INFINITY;
        let mut temp = *self;

        // Process each limb
        for &limb in scalar_bytes.iter() {
            let mut bits = limb;
            for _ in 0..64 {
                if bits & 1 == 1 {
                    result = result + temp;
                }
                temp = temp.double();
                bits >>= 1;
            }
        }

        result
    }
}

// Implement addition for projective points
impl Add for Projective {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        // Handle infinity cases
        if self.is_infinity() {
            return other;
        }
        if other.is_infinity() {
            return self;
        }

        // For simplicity and correctness, convert to affine, add, and convert back
        let affine1 = self.to_affine();
        let affine2 = other.to_affine();
        let result_affine = affine1 + affine2;
        Self::from_affine(&result_affine)
    }
}

impl AddAssign for Projective {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl Sub for Projective {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self + other.negate()
    }
}

impl SubAssign for Projective {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl Neg for Projective {
    type Output = Self;

    fn neg(self) -> Self {
        self.negate()
    }
}

// Scalar multiplication
impl Mul<ScalarField> for Projective {
    type Output = Self;

    fn mul(self, scalar: ScalarField) -> Self {
        self.scalar_mul(&scalar)
    }
}

impl Mul<&ScalarField> for Projective {
    type Output = Self;

    fn mul(self, scalar: &ScalarField) -> Self {
        self.scalar_mul(scalar)
    }
}

impl Mul<Projective> for ScalarField {
    type Output = Projective;

    fn mul(self, point: Projective) -> Projective {
        point.scalar_mul(&self)
    }
}

impl Mul<&Projective> for ScalarField {
    type Output = Projective;

    fn mul(self, point: &Projective) -> Projective {
        point.scalar_mul(&self)
    }
}

// Conversions
impl From<Affine> for Projective {
    fn from(point: Affine) -> Self {
        Projective::from_affine(&point)
    }
}

impl From<&Affine> for Projective {
    fn from(point: &Affine) -> Self {
        Projective::from_affine(point)
    }
}

impl From<Projective> for Affine {
    fn from(point: Projective) -> Self {
        point.to_affine()
    }
}

impl From<&Projective> for Affine {
    fn from(point: &Projective) -> Self {
        point.to_affine()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_infinity() {
        let inf = Projective::INFINITY;
        assert!(inf.is_infinity());
        assert!(inf.is_on_curve());
    }

    #[test]
    fn test_generator_on_curve() {
        let g = Projective::generator();
        assert!(g.is_on_curve(), "Generator point is not on the curve");
        assert!(!g.is_infinity());
    }

    #[test]
    fn test_generator_pedersen_on_curve() {
        let g = Projective::generator_pedersen();
        assert!(
            g.is_on_curve(),
            "Pedersen generator point is not on the curve"
        );
        assert!(!g.is_infinity());
    }

    #[test]
    fn test_conversion_affine_projective() {
        let affine = Affine::generator();
        let projective = Projective::from_affine(&affine);
        let back_to_affine = projective.to_affine();

        assert_eq!(affine, back_to_affine);
    }

    #[test]
    fn test_point_addition_with_infinity() {
        let g = Projective::generator();
        let inf = Projective::INFINITY;

        assert_eq!(g + inf, g);
        assert_eq!(inf + g, g);
        assert_eq!(inf + inf, inf);
    }

    #[test]
    fn test_point_doubling() {
        let g = Projective::generator();
        let g2 = g.double();

        assert!(g2.is_on_curve(), "Doubled point is not on the curve");
        assert_eq!(g + g, g2);
    }

    #[test]
    fn test_point_negation() {
        let g = Projective::generator();
        let neg_g = g.negate();

        assert!(neg_g.is_on_curve());
        assert_eq!(g + neg_g, Projective::INFINITY);
    }

    #[test]
    fn test_scalar_multiplication() {
        let g = Projective::generator();
        let scalar = ScalarField::from_canonical_u64(5);
        let result = g.scalar_mul(&scalar);

        // 5*G = G + G + G + G + G
        let expected = g + g + g + g + g;
        assert_eq!(result, expected);
        assert!(result.is_on_curve());
    }

    #[test]
    fn test_scalar_mul_zero() {
        let g = Projective::generator();
        let zero = ScalarField::ZERO;
        let result = g.scalar_mul(&zero);

        assert_eq!(result, Projective::INFINITY);
    }

    #[test]
    fn test_scalar_mul_one() {
        let g = Projective::generator();
        let one = ScalarField::ONE;
        let result = g.scalar_mul(&one);

        assert_eq!(result, g);
    }

    #[test]
    fn test_associativity() {
        let g = Projective::generator();
        let a = ScalarField::from_canonical_u64(3);
        let b = ScalarField::from_canonical_u64(5);

        // (a + b) * G = a*G + b*G
        let left = g.scalar_mul(&(a + b));
        let right = g.scalar_mul(&a) + g.scalar_mul(&b);

        assert_eq!(left, right);
    }

    #[test]
    fn test_affine_projective_addition_consistency() {
        let g_affine = Affine::generator();
        let g_projective = Projective::generator();

        // Test that affine and projective addition give the same result
        let affine_sum = g_affine + g_affine;
        let projective_sum = g_projective + g_projective;

        assert_eq!(affine_sum, projective_sum.to_affine());
    }

    #[test]
    fn test_affine_projective_scalar_mul_consistency() {
        let g_affine = Affine::generator();
        let g_projective = Projective::generator();
        let scalar = ScalarField::from_canonical_u64(42);

        // Test that affine and projective scalar multiplication give the same result
        let affine_result = g_affine.scalar_mul(&scalar);
        let projective_result = g_projective.scalar_mul(&scalar);

        assert_eq!(affine_result, projective_result.to_affine());
    }
}
