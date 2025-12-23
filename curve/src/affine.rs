// E(GF((2130706433)^8)) : y^2 = x^3 + 3u*x + 42639
// E generator point (from SSWU on 'ZKM2'): (1195559694*u^7 + 1368232771*u^6 + 438909494*u^5 + 1825476283*u^4 + 1299273209*u^3 + 2115217807*u^2 + 1763905369*u + 1813646457 : 2077084094*u^7 + 434578416*u^6 + 125328769*u^5 + 1286889583*u^4 + 655051022*u^3 + 1365273355*u^2 + 840779000*u + 376996212 : 1)
// E generator point (from SSWU on 'ZKM2 - Pedersen'): (1741425845*u^7 + 1750752810*u^6 + 7156361*u^5 + 1949220725*u^4 + 543192455*u^3 + 358531926*u^2 + 550988532*u + 1709677626 : 71894712*u^7 + 1876016551*u^6 + 1684498755*u^5 + 598910111*u^4 + 156828552*u^3 + 978667041*u^2 + 1399061592*u + 548133034 : 1)
// Curve prime order: 424804331891979973455971894938199991855800421968298112348210302325367590273 (248 bits)
// Curve prime order (hex): 0xf06e44682c2aa440f5f26a5ae1748ff85ccc2efc3068faf2154ff8a2e94d81
// Curve prime order 2-adicity v2(q-1): 7
// Curve cofactor: 1
// Curve security (Pollard-Rho): 123.78
// Curve embedding degree: 13275135371624374170499121716818749745493763186509316010881571947667737196 (>2^242)
// Twist security (Pollard-Rho): 120.86

use crate::basefield::{from_coeffs, BaseField};
use crate::{double_scalar_mul_basepoint_affine, mul_generator_affine, Group, ScalarField};
use core::ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign};
use p3_field::{Field, PrimeCharacteristicRing};
use p3_koala_bear::KoalaBear;
use serde::{Deserialize, Serialize};

/// Affine point on the elliptic curve.
/// Represents a point in affine coordinates (x, y) or the point at infinity.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Affine {
    /// The x-coordinate of the point (Fp8 element)
    pub x: BaseField,
    /// The y-coordinate of the point (Fp8 element)
    pub y: BaseField,
    /// Whether this point is the point at infinity (identity element)
    pub is_infinity: bool,
}

impl Affine {
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

    /// The point at infinity (identity element)
    pub const INFINITY: Self = Affine {
        x: BaseField::ZERO,
        y: BaseField::ZERO,
        is_infinity: true,
    };

    /// Create a new affine point.
    pub fn new(x: BaseField, y: BaseField) -> Self {
        Affine {
            x,
            y,
            is_infinity: false,
        }
    }

    /// Check if this point is the point at infinity.
    #[inline]
    pub fn is_infinity(&self) -> bool {
        self.is_infinity
    }

    /// Check if a point is on the curve: y^2 = x^3 + a*x + b.
    pub fn is_on_curve(&self) -> bool {
        if self.is_infinity {
            return true;
        }

        let y2 = self.y * self.y;
        let x2 = self.x * self.x;
        let x3 = x2 * self.x;
        let ax = Self::curve_a() * self.x;
        let rhs = x3 + ax + Self::curve_b();

        y2 == rhs
    }

    /// Generator point from SSWU on 'ZKM2'.
    pub fn generator() -> Self {
        // (1195559694*u^7 + 1368232771*u^6 + 438909494*u^5 + 1825476283*u^4 +
        //  1299273209*u^3 + 2115217807*u^2 + 1763905369*u + 1813646457 :
        //  2077084094*u^7 + 434578416*u^6 + 125328769*u^5 + 1286889583*u^4 +
        //  655051022*u^3 + 1365273355*u^2 + 840779000*u + 376996212 : 1)

        let x = from_coeffs([
            KoalaBear::new(1813646457),
            KoalaBear::new(1763905369),
            KoalaBear::new(2115217807),
            KoalaBear::new(1299273209),
            KoalaBear::new(1825476283),
            KoalaBear::new(438909494),
            KoalaBear::new(1368232771),
            KoalaBear::new(1195559694),
        ]);

        let y = from_coeffs([
            KoalaBear::new(376996212),
            KoalaBear::new(840779000),
            KoalaBear::new(1365273355),
            KoalaBear::new(655051022),
            KoalaBear::new(1286889583),
            KoalaBear::new(125328769),
            KoalaBear::new(434578416),
            KoalaBear::new(2077084094),
        ]);

        Affine::new(x, y)
    }

    /// Alternative generator point from SSWU on 'ZKM2 - Pedersen'.
    pub fn generator_pedersen() -> Self {
        // (1741425845*u^7 + 1750752810*u^6 + 7156361*u^5 + 1949220725*u^4 +
        //  543192455*u^3 + 358531926*u^2 + 550988532*u + 1709677626 :
        //  71894712*u^7 + 1876016551*u^6 + 1684498755*u^5 + 598910111*u^4 +
        //  156828552*u^3 + 978667041*u^2 + 1399061592*u + 548133034 : 1)

        let x = from_coeffs([
            KoalaBear::new(1709677626),
            KoalaBear::new(550988532),
            KoalaBear::new(358531926),
            KoalaBear::new(543192455),
            KoalaBear::new(1949220725),
            KoalaBear::new(7156361),
            KoalaBear::new(1750752810),
            KoalaBear::new(1741425845),
        ]);

        let y = from_coeffs([
            KoalaBear::new(548133034),
            KoalaBear::new(1399061592),
            KoalaBear::new(978667041),
            KoalaBear::new(156828552),
            KoalaBear::new(598910111),
            KoalaBear::new(1684498755),
            KoalaBear::new(1876016551),
            KoalaBear::new(71894712),
        ]);

        Affine::new(x, y)
    }

    /// Point doubling: 2*P.
    pub fn double(&self) -> Self {
        if self.is_infinity {
            return *self;
        }

        // If y = 0, then 2P = O
        if self.y.is_zero() {
            return Self::INFINITY;
        }

        // Compute slope: λ = (3x^2 + a) / (2y)
        let x2 = self.x * self.x;
        let three_x2 = x2 + x2 + x2;
        let numerator = three_x2 + Self::curve_a();
        let denominator = self.y + self.y;
        let lambda = numerator / denominator;

        // x_r = λ^2 - 2x
        let lambda2 = lambda * lambda;
        let x_r = lambda2 - self.x - self.x;

        // y_r = λ(x - x_r) - y
        let y_r = lambda * (self.x - x_r) - self.y;

        Affine::new(x_r, y_r)
    }

    /// Negate a point.
    pub fn negate(&self) -> Self {
        if self.is_infinity {
            return *self;
        }
        Affine::new(self.x, -self.y)
    }

    /// Multiply the fixed generator using a precomputed table.
    pub fn mul_generator(scalar: &ScalarField) -> Self {
        mul_generator_affine(scalar)
    }

    /// Compute a * G + b * P, where G is the fixed generator.
    pub fn double_scalar_mul_basepoint(a: &ScalarField, b: &ScalarField, point: &Self) -> Self {
        double_scalar_mul_basepoint_affine(a, b, point)
    }
}

impl Group for Affine {
    type Scalar = ScalarField;

    #[inline]
    fn identity() -> Self {
        Self::INFINITY
    }

    #[inline]
    fn is_identity(&self) -> bool {
        self.is_infinity
    }

    #[inline]
    fn generator() -> Self {
        Affine::generator()
    }

    #[inline]
    fn mul_generator(scalar: &ScalarField) -> Self {
        Affine::mul_generator(scalar)
    }

    #[inline]
    fn double(&self) -> Self {
        Self::double(self)
    }

    #[inline]
    fn negate(&self) -> Self {
        Self::negate(self)
    }
}

// Implement addition for affine points
impl Add for Affine {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        // Handle infinity cases
        if self.is_infinity {
            return other;
        }
        if other.is_infinity {
            return self;
        }

        // Check if points are the same
        if self.x == other.x {
            if self.y == other.y {
                // Point doubling
                return self.double();
            } else {
                // Points are inverses, return infinity
                return Self::INFINITY;
            }
        }

        // Regular point addition
        // λ = (y2 - y1) / (x2 - x1)
        let numerator = other.y - self.y;
        let denominator = other.x - self.x;
        let lambda = numerator / denominator;

        // x_r = λ^2 - x1 - x2
        let lambda2 = lambda * lambda;
        let x_r = lambda2 - self.x - other.x;

        // y_r = λ(x1 - x_r) - y1
        let y_r = lambda * (self.x - x_r) - self.y;

        Affine::new(x_r, y_r)
    }
}

impl AddAssign for Affine {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Sub for Affine {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self + other.negate()
    }
}

impl SubAssign for Affine {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl Neg for Affine {
    type Output = Self;

    fn neg(self) -> Self {
        self.negate()
    }
}

// Scalar multiplication
impl Mul<ScalarField> for Affine {
    type Output = Self;

    fn mul(self, scalar: ScalarField) -> Self {
        <Self as Group>::scalar_mul(&self, &scalar)
    }
}

impl Mul<&ScalarField> for Affine {
    type Output = Self;

    fn mul(self, scalar: &ScalarField) -> Self {
        <Self as Group>::scalar_mul(&self, scalar)
    }
}

impl Mul<Affine> for ScalarField {
    type Output = Affine;

    fn mul(self, point: Affine) -> Affine {
        <Affine as Group>::scalar_mul(&point, &self)
    }
}

impl Mul<&Affine> for ScalarField {
    type Output = Affine;

    fn mul(self, point: &Affine) -> Affine {
        <Affine as Group>::scalar_mul(point, &self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Group;

    #[test]
    fn test_infinity() {
        let inf = Affine::INFINITY;
        assert!(inf.is_infinity());
        assert!(inf.is_on_curve());
    }

    #[test]
    fn test_generator_on_curve() {
        let g = Affine::generator();
        assert!(g.is_on_curve(), "Generator point is not on the curve");
        assert!(!g.is_infinity());
    }

    #[test]
    fn test_generator_pedersen_on_curve() {
        let g = Affine::generator_pedersen();
        assert!(
            g.is_on_curve(),
            "Pedersen generator point is not on the curve"
        );
        assert!(!g.is_infinity());
    }

    #[test]
    fn test_point_addition_with_infinity() {
        let g = Affine::generator();
        let inf = Affine::INFINITY;

        assert_eq!(g + inf, g);
        assert_eq!(inf + g, g);
        assert_eq!(inf + inf, inf);
    }

    #[test]
    fn test_point_doubling() {
        let g = Affine::generator();
        let g2 = g.double();

        assert!(g2.is_on_curve(), "Doubled point is not on the curve");
        assert_eq!(g + g, g2);
    }

    #[test]
    fn test_point_negation() {
        let g = Affine::generator();
        let neg_g = g.negate();

        assert!(neg_g.is_on_curve());
        assert_eq!(g + neg_g, Affine::INFINITY);
    }

    #[test]
    fn test_scalar_multiplication() {
        let g = Affine::generator();
        let scalar = ScalarField::from_canonical_u64(5);
        let result = g.scalar_mul(&scalar);

        // 5*G = G + G + G + G + G
        let expected = g + g + g + g + g;
        assert_eq!(result, expected);
        assert!(result.is_on_curve());
    }

    #[test]
    fn test_scalar_mul_zero() {
        let g = Affine::generator();
        let zero = ScalarField::ZERO;
        let result = g.scalar_mul(&zero);

        assert_eq!(result, Affine::INFINITY);
    }

    #[test]
    fn test_scalar_mul_one() {
        let g = Affine::generator();
        let one = ScalarField::ONE;
        let result = g.scalar_mul(&one);

        assert_eq!(result, g);
    }

    #[test]
    fn test_associativity() {
        let g = Affine::generator();
        let a = ScalarField::from_canonical_u64(3);
        let b = ScalarField::from_canonical_u64(5);

        // (a + b) * G = a*G + b*G
        let left = g.scalar_mul(&(a + b));
        let right = g.scalar_mul(&a) + g.scalar_mul(&b);

        assert_eq!(left, right);
    }

    #[test]
    fn test_windowed_scalar_mul() {
        let g = Affine::generator();
        let scalar = ScalarField::from_canonical_u64(123456);

        // Compare windowed and standard scalar multiplication
        let result1 = g.scalar_mul(&scalar);
        let result2 = g.scalar_mul_windowed(&scalar);

        assert_eq!(result1, result2);
        assert!(result1.is_on_curve());
    }

    #[test]
    fn test_mul_generator() {
        let scalar = ScalarField::from_canonical_u64(123456);
        let result = Affine::mul_generator(&scalar);
        let expected = Affine::generator().scalar_mul(&scalar);

        assert_eq!(result, expected);
        assert!(result.is_on_curve());
    }

    #[test]
    fn test_group_mul_generator_matches_scalar_mul() {
        let scalar = ScalarField::from_canonical_u64(424242);
        let result = <Affine as Group>::mul_generator(&scalar);
        let expected = Affine::generator().scalar_mul(&scalar);

        assert_eq!(result, expected);
    }

    #[test]
    fn test_multi_scalar_mul() {
        let g = Affine::generator();
        let h = Affine::generator_pedersen();

        let a = ScalarField::from_canonical_u64(7);
        let b = ScalarField::from_canonical_u64(11);

        let points = vec![g, h];
        let scalars = vec![a, b];

        let result = <Affine as Group>::multi_scalar_mul(&points, &scalars);
        let expected = g.scalar_mul(&a) + h.scalar_mul(&b);

        assert_eq!(result, expected);
        assert!(result.is_on_curve());
    }

    #[test]
    fn test_mul_u64() {
        let g = Affine::generator();
        let n = 42u64;

        let result1 = g.mul_u64(n);
        let result2 = g.scalar_mul(&ScalarField::from_canonical_u64(n));

        assert_eq!(result1, result2);
        assert!(result1.is_on_curve());
    }

    #[test]
    fn test_identity() {
        let id = <Affine as Group>::identity();
        assert!(id.is_identity());
        assert_eq!(id, Affine::INFINITY);

        let g = Affine::generator();
        assert_eq!(g + id, g);
        assert_eq!(id + g, g);
    }

    #[test]
    fn test_group_properties() {
        let g = Affine::generator();

        // Test that doubling is the same as adding to itself
        assert_eq!(g.double(), g + g);

        // Test that triple is correct
        let triple1 = g + g + g;
        let triple2 = g.mul_u64(3);
        assert_eq!(triple1, triple2);

        // Test inverse property
        let h = g.mul_u64(5);
        let neg_h = -h;
        assert_eq!(h + neg_h, Affine::INFINITY);
    }
}
