//! Scalar field of the curve. p = 0xf06e44682c2aa440f5f26a5ae1748ff85ccc2efc3068faf2154ff8a2e94d81
//!
//! This implementation uses Montgomery form for efficient modular arithmetic.
//! The field element is represented as [u64; 4] in little-endian order.

extern crate alloc;

use alloc::vec::Vec;
use core::fmt::{self, Debug, Display, Formatter};
use core::hash::{Hash, Hasher};
use core::iter::{Product, Sum};
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use num_bigint::BigUint;
use p3_field::integers::QuotientMap;
use p3_field::{Field, Packable, PrimeCharacteristicRing, PrimeField, RawDataSerializable};
use rand::distr::{Distribution, StandardUniform};
use rand::Rng;
use serde::{Deserialize, Serialize};

/// Scalar field element for the curve
/// Represented in Montgomery form with [u64; 4]
#[derive(Copy, Clone, Default, Eq, PartialEq, Ord, PartialOrd, Serialize, Deserialize)]
pub struct ScalarField {
    /// Montgomery form: value * R mod p, where R = 2^256
    limbs: [u64; 4],
}

// Field modulus: p = 0xf06e44682c2aa440f5f26a5ae1748ff85ccc2efc3068faf2154ff8a2e94d81
const MODULUS: [u64; 4] = [
    0xf2154ff8a2e94d81,
    0xf85ccc2efc3068fa,
    0x40f5f26a5ae1748f,
    0x00f06e44682c2aa4,
];

// R = 2^256 mod p (Montgomery parameter)
const R: [u64; 4] = [
    0xc95b07d2e81da6f0,
    0x1d670e140c90755e,
    0xfaae6eff70742708,
    0x008ad7515112b17a,
];

// R^2 = 2^512 mod p (for Montgomery conversion)
const R2: [u64; 4] = [
    0x23eabb3eaf3c12e3,
    0xefbc3b2088f7b0f7,
    0x0943bc9a31f37148,
    0x004497b874228e49,
];

// R^3 = 2^768 mod p (for efficient conversion)
#[allow(dead_code)]
const R3: [u64; 4] = [
    0xcfda2e5499aadbfe,
    0x3b9752a61aed3bc4,
    0xddf904d7687cb0cb,
    0x001e3196ff5393fa,
];

// -p^{-1} mod 2^64 (Montgomery parameter mu)
const MU: u64 = 0x921d21f874d30d7f;

impl ScalarField {
    /// Zero element (in Montgomery form)
    pub const ZERO: Self = ScalarField {
        limbs: [0, 0, 0, 0],
    };

    /// One element (in Montgomery form: R mod p)
    pub const ONE: Self = ScalarField { limbs: R };

    /// Create a new scalar field element from a u64 value
    #[inline]
    pub fn from_canonical_u64(val: u64) -> Self {
        // Convert to Montgomery form: val * R^2 * R^{-1} = val * R
        let result = ScalarField {
            limbs: [val, 0, 0, 0],
        };
        montgomery_mul(result, ScalarField { limbs: R2 })
    }

    /// Convert from Montgomery form to canonical form
    #[inline]
    pub fn to_canonical_u64_vec(&self) -> [u64; 4] {
        // Multiply by 1 to get out of Montgomery form
        let one = ScalarField {
            limbs: [1, 0, 0, 0],
        };
        montgomery_mul(*self, one).limbs
    }

    #[inline]
    fn from_canonical_limbs(limbs: [u64; 4]) -> Self {
        montgomery_mul(ScalarField { limbs }, ScalarField { limbs: R2 })
    }
}

/// Helper: Add two 256-bit numbers mod p
#[inline]
const fn add_mod(a: [u64; 4], b: [u64; 4]) -> [u64; 4] {
    let (r0, carry) = a[0].overflowing_add(b[0]);
    let (r1, carry) = carrying_add(a[1], b[1], carry);
    let (r2, carry) = carrying_add(a[2], b[2], carry);
    let (r3, carry) = carrying_add(a[3], b[3], carry);

    // Subtract modulus if we overflowed or result >= p
    let (s0, borrow) = r0.overflowing_sub(MODULUS[0]);
    let (s1, borrow) = borrowing_sub(r1, MODULUS[1], borrow);
    let (s2, borrow) = borrowing_sub(r2, MODULUS[2], borrow);
    let (s3, borrow) = borrowing_sub(r3, MODULUS[3], borrow);

    // If no borrow and no carry, or if carry, we need to subtract
    if carry || !borrow {
        [s0, s1, s2, s3]
    } else {
        [r0, r1, r2, r3]
    }
}

/// Helper: Subtract two 256-bit numbers mod p
#[inline]
const fn sub_mod(a: [u64; 4], b: [u64; 4]) -> [u64; 4] {
    let (r0, borrow) = a[0].overflowing_sub(b[0]);
    let (r1, borrow) = borrowing_sub(a[1], b[1], borrow);
    let (r2, borrow) = borrowing_sub(a[2], b[2], borrow);
    let (r3, borrow) = borrowing_sub(a[3], b[3], borrow);

    // Add modulus if we underflowed
    if borrow {
        let (r0, carry) = r0.overflowing_add(MODULUS[0]);
        let (r1, carry) = carrying_add(r1, MODULUS[1], carry);
        let (r2, carry) = carrying_add(r2, MODULUS[2], carry);
        let (r3, _) = carrying_add(r3, MODULUS[3], carry);
        [r0, r1, r2, r3]
    } else {
        [r0, r1, r2, r3]
    }
}

/// Helper: Negate a 256-bit number mod p
#[inline]
const fn neg_mod(a: [u64; 4]) -> [u64; 4] {
    if a[0] == 0 && a[1] == 0 && a[2] == 0 && a[3] == 0 {
        return [0, 0, 0, 0];
    }
    sub_mod(MODULUS, a)
}

#[inline]
const fn is_canonical(limbs: [u64; 4]) -> bool {
    let (_, borrow) = limbs[0].overflowing_sub(MODULUS[0]);
    let (_, borrow) = borrowing_sub(limbs[1], MODULUS[1], borrow);
    let (_, borrow) = borrowing_sub(limbs[2], MODULUS[2], borrow);
    let (_, borrow) = borrowing_sub(limbs[3], MODULUS[3], borrow);
    borrow
}

/// Helper: Carrying addition
#[inline]
const fn carrying_add(a: u64, b: u64, carry: bool) -> (u64, bool) {
    let (sum, overflow1) = a.overflowing_add(b);
    let (sum, overflow2) = sum.overflowing_add(carry as u64);
    (sum, overflow1 || overflow2)
}

/// Helper: Borrowing subtraction
#[inline]
const fn borrowing_sub(a: u64, b: u64, borrow: bool) -> (u64, bool) {
    let (diff, overflow1) = a.overflowing_sub(b);
    let (diff, overflow2) = diff.overflowing_sub(borrow as u64);
    (diff, overflow1 || overflow2)
}

/// Montgomery multiplication: (a * b * R^{-1}) mod p
#[inline]
fn montgomery_mul(a: ScalarField, b: ScalarField) -> ScalarField {
    // Compute a * b
    let mut t = [0u64; 8];

    for i in 0..4 {
        let mut carry = 0u128;
        for j in 0..4 {
            let product = (a.limbs[i] as u128) * (b.limbs[j] as u128) + (t[i + j] as u128) + carry;
            t[i + j] = product as u64;
            carry = product >> 64;
        }
        t[i + 4] = carry as u64;
    }

    // Montgomery reduction
    for i in 0..4 {
        let k = t[i].wrapping_mul(MU);
        let mut carry = 0u128;

        for j in 0..4 {
            let product = (k as u128) * (MODULUS[j] as u128) + (t[i + j] as u128) + carry;
            t[i + j] = product as u64;
            carry = product >> 64;
        }

        for j in 4..8 - i {
            let sum = (t[i + j] as u128) + carry;
            t[i + j] = sum as u64;
            carry = sum >> 64;
        }
    }

    // Extract high half and conditionally subtract p
    let result = [t[4], t[5], t[6], t[7]];

    // Check if result >= p
    let (_, borrow) = result[0].overflowing_sub(MODULUS[0]);
    let (_, borrow) = borrowing_sub(result[1], MODULUS[1], borrow);
    let (_, borrow) = borrowing_sub(result[2], MODULUS[2], borrow);
    let (_, borrow) = borrowing_sub(result[3], MODULUS[3], borrow);

    if !borrow {
        ScalarField {
            limbs: sub_mod(result, MODULUS),
        }
    } else {
        ScalarField { limbs: result }
    }
}

// Implement PrimeCharacteristicRing trait
impl PrimeCharacteristicRing for ScalarField {
    type PrimeSubfield = Self;

    const ZERO: Self = Self::ZERO;
    const ONE: Self = Self::ONE;
    const TWO: Self = ScalarField {
        limbs: add_mod(R, R),
    };
    const NEG_ONE: Self = ScalarField {
        limbs: sub_mod(MODULUS, R),
    };

    #[inline]
    fn from_prime_subfield(elem: Self::PrimeSubfield) -> Self {
        elem
    }

    #[inline]
    fn halve(&self) -> Self {
        // Compute (self + p) / 2 if odd, else self / 2
        let is_odd = self.limbs[0] & 1 == 1;

        if is_odd {
            // Add p then shift right
            let sum = add_mod(self.limbs, MODULUS);
            let mut result = [0u64; 4];
            result[0] = (sum[0] >> 1) | (sum[1] << 63);
            result[1] = (sum[1] >> 1) | (sum[2] << 63);
            result[2] = (sum[2] >> 1) | (sum[3] << 63);
            result[3] = sum[3] >> 1;
            ScalarField { limbs: result }
        } else {
            // Just shift right
            let mut result = [0u64; 4];
            result[0] = (self.limbs[0] >> 1) | (self.limbs[1] << 63);
            result[1] = (self.limbs[1] >> 1) | (self.limbs[2] << 63);
            result[2] = (self.limbs[2] >> 1) | (self.limbs[3] << 63);
            result[3] = self.limbs[3] >> 1;
            ScalarField { limbs: result }
        }
    }
}

// Implement Packable trait - required for Field
impl Packable for ScalarField {}

// Implement RawDataSerializable trait
impl RawDataSerializable for ScalarField {
    const NUM_BYTES: usize = 32;

    fn into_bytes(self) -> impl IntoIterator<Item = u8> {
        let canonical = self.to_canonical_u64_vec();
        let mut bytes = Vec::with_capacity(32);
        for &limb in &canonical {
            bytes.extend_from_slice(&limb.to_le_bytes());
        }
        bytes
    }
}

impl Distribution<ScalarField> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> ScalarField {
        loop {
            let mut bytes: [u8; 32] = rng.random();
            bytes[31] = 0;

            let limbs = [
                u64::from_le_bytes(bytes[0..8].try_into().unwrap()),
                u64::from_le_bytes(bytes[8..16].try_into().unwrap()),
                u64::from_le_bytes(bytes[16..24].try_into().unwrap()),
                u64::from_le_bytes(bytes[24..32].try_into().unwrap()),
            ];

            if is_canonical(limbs) {
                return ScalarField::from_canonical_limbs(limbs);
            }
        }
    }
}

// Implement Field trait from p3_field
impl Field for ScalarField {
    type Packing = Self;

    // Generator g = 5 in Montgomery form
    const GENERATOR: Self = ScalarField {
        limbs: [
            0x0a9c872d42c1a7ae,
            0xa249ae06467178e4,
            0x637c46287c81da08,
            0x00d5580dc505221e,
        ],
    };

    fn try_inverse(&self) -> Option<Self> {
        if *self == Self::ZERO {
            None
        } else {
            Some(self.inverse())
        }
    }

    fn order() -> BigUint {
        // Return the field modulus as a BigUint
        let mut bytes = alloc::vec::Vec::with_capacity(32);
        for &limb in &MODULUS {
            bytes.extend_from_slice(&limb.to_le_bytes());
        }
        BigUint::from_bytes_le(&bytes)
    }
}

// Implement PrimeField trait
impl PrimeField for ScalarField {
    fn as_canonical_biguint(&self) -> BigUint {
        let canonical = self.to_canonical_u64_vec();
        let mut bytes = Vec::with_capacity(32);
        for &limb in &canonical {
            bytes.extend_from_slice(&limb.to_le_bytes());
        }
        BigUint::from_bytes_le(&bytes)
    }
}

// Implement QuotientMap for common integer types
impl QuotientMap<u8> for ScalarField {
    fn from_int(int: u8) -> Self {
        Self::from_canonical_u64(int as u64)
    }

    fn from_canonical_checked(int: u8) -> Option<Self> {
        Some(Self::from_canonical_u64(int as u64))
    }

    unsafe fn from_canonical_unchecked(int: u8) -> Self {
        Self::from_canonical_u64(int as u64)
    }
}

impl QuotientMap<u16> for ScalarField {
    fn from_int(int: u16) -> Self {
        Self::from_canonical_u64(int as u64)
    }

    fn from_canonical_checked(int: u16) -> Option<Self> {
        Some(Self::from_canonical_u64(int as u64))
    }

    unsafe fn from_canonical_unchecked(int: u16) -> Self {
        Self::from_canonical_u64(int as u64)
    }
}

impl QuotientMap<u32> for ScalarField {
    fn from_int(int: u32) -> Self {
        Self::from_canonical_u64(int as u64)
    }

    fn from_canonical_checked(int: u32) -> Option<Self> {
        Some(Self::from_canonical_u64(int as u64))
    }

    unsafe fn from_canonical_unchecked(int: u32) -> Self {
        Self::from_canonical_u64(int as u64)
    }
}

impl QuotientMap<u64> for ScalarField {
    fn from_int(int: u64) -> Self {
        Self::from_canonical_u64(int)
    }

    fn from_canonical_checked(int: u64) -> Option<Self> {
        // For simplicity, accept all u64 values (they're much smaller than our modulus)
        Some(Self::from_canonical_u64(int))
    }

    unsafe fn from_canonical_unchecked(int: u64) -> Self {
        Self::from_canonical_u64(int)
    }
}

impl QuotientMap<u128> for ScalarField {
    fn from_int(int: u128) -> Self {
        // Take the lower 64 bits for u128 (mod p)
        Self::from_canonical_u64(int as u64)
    }

    fn from_canonical_checked(int: u128) -> Option<Self> {
        if int <= u64::MAX as u128 {
            Some(Self::from_canonical_u64(int as u64))
        } else {
            None
        }
    }

    unsafe fn from_canonical_unchecked(int: u128) -> Self {
        Self::from_canonical_u64(int as u64)
    }
}

impl QuotientMap<i8> for ScalarField {
    fn from_int(int: i8) -> Self {
        if int >= 0 {
            Self::from_canonical_u64(int as u64)
        } else {
            -Self::from_canonical_u64((-int) as u64)
        }
    }

    fn from_canonical_checked(int: i8) -> Option<Self> {
        Some(Self::from_int(int))
    }

    unsafe fn from_canonical_unchecked(int: i8) -> Self {
        Self::from_int(int)
    }
}

impl QuotientMap<i16> for ScalarField {
    fn from_int(int: i16) -> Self {
        if int >= 0 {
            Self::from_canonical_u64(int as u64)
        } else {
            -Self::from_canonical_u64((-int) as u64)
        }
    }

    fn from_canonical_checked(int: i16) -> Option<Self> {
        Some(Self::from_int(int))
    }

    unsafe fn from_canonical_unchecked(int: i16) -> Self {
        Self::from_int(int)
    }
}

impl QuotientMap<i32> for ScalarField {
    fn from_int(int: i32) -> Self {
        if int >= 0 {
            Self::from_canonical_u64(int as u64)
        } else {
            -Self::from_canonical_u64((-int) as u64)
        }
    }

    fn from_canonical_checked(int: i32) -> Option<Self> {
        Some(Self::from_int(int))
    }

    unsafe fn from_canonical_unchecked(int: i32) -> Self {
        Self::from_int(int)
    }
}

impl QuotientMap<i64> for ScalarField {
    fn from_int(int: i64) -> Self {
        if int >= 0 {
            Self::from_canonical_u64(int as u64)
        } else {
            -Self::from_canonical_u64((-int) as u64)
        }
    }

    fn from_canonical_checked(int: i64) -> Option<Self> {
        Some(Self::from_int(int))
    }

    unsafe fn from_canonical_unchecked(int: i64) -> Self {
        Self::from_int(int)
    }
}

impl QuotientMap<i128> for ScalarField {
    fn from_int(int: i128) -> Self {
        if int >= 0 {
            Self::from_canonical_u64(int as u64)
        } else {
            -Self::from_canonical_u64((-int) as u64)
        }
    }

    fn from_canonical_checked(int: i128) -> Option<Self> {
        Some(Self::from_int(int))
    }

    unsafe fn from_canonical_unchecked(int: i128) -> Self {
        Self::from_int(int)
    }
}

// Arithmetic operations
impl Add for ScalarField {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        ScalarField {
            limbs: add_mod(self.limbs, rhs.limbs),
        }
    }
}

impl AddAssign for ScalarField {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for ScalarField {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        ScalarField {
            limbs: sub_mod(self.limbs, rhs.limbs),
        }
    }
}

impl SubAssign for ScalarField {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Neg for ScalarField {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        ScalarField {
            limbs: neg_mod(self.limbs),
        }
    }
}

impl Mul for ScalarField {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        montgomery_mul(self, rhs)
    }
}

impl MulAssign for ScalarField {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Div for ScalarField {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self {
        self * rhs.inverse()
    }
}

impl DivAssign for ScalarField {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl Sum for ScalarField {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::ZERO, |acc, x| acc + x)
    }
}

impl Product for ScalarField {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::ONE, |acc, x| acc * x)
    }
}

// Display and Debug
impl Display for ScalarField {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let canonical = self.to_canonical_u64_vec();
        write!(
            f,
            "0x{:016x}{:016x}{:016x}{:016x}",
            canonical[3], canonical[2], canonical[1], canonical[0]
        )
    }
}

impl Debug for ScalarField {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "ScalarField({})", self)
    }
}

impl Hash for ScalarField {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.limbs.hash(state);
    }
}

// Inverse using Fermat's little theorem: a^{-1} = a^{p-2}
impl ScalarField {
    /// Compute multiplicative inverse using binary exponentiation
    pub fn inverse(&self) -> Self {
        // p - 2 for Fermat's little theorem
        let exp = sub_mod(MODULUS, [2, 0, 0, 0]);
        self.pow_vartime(exp)
    }

    /// Variable-time exponentiation
    fn pow_vartime(&self, exp: [u64; 4]) -> Self {
        if self.is_zero() {
            return Self::ZERO;
        }

        let mut result = Self::ONE;
        let mut base = *self;

        // Process bits from least significant to most significant
        for &limb in exp.iter() {
            let mut remaining = limb;
            for _ in 0..64 {
                if remaining & 1 == 1 {
                    result = result * base;
                }
                base = base * base;
                remaining >>= 1;
            }
        }

        result
    }

    /// Check if this field element is zero
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.limbs == [0, 0, 0, 0]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_one() {
        assert_eq!(ScalarField::ZERO + ScalarField::ZERO, ScalarField::ZERO);
        assert_eq!(ScalarField::ONE * ScalarField::ONE, ScalarField::ONE);
        assert_eq!(ScalarField::ZERO * ScalarField::ONE, ScalarField::ZERO);
        assert_eq!(ScalarField::ONE + ScalarField::ZERO, ScalarField::ONE);
    }

    #[test]
    fn test_addition() {
        let a = ScalarField::from_canonical_u64(5);
        let b = ScalarField::from_canonical_u64(7);
        let c = a + b;
        assert_eq!(c, ScalarField::from_canonical_u64(12));
    }

    #[test]
    fn test_subtraction() {
        let a = ScalarField::from_canonical_u64(10);
        let b = ScalarField::from_canonical_u64(3);
        let c = a - b;
        assert_eq!(c, ScalarField::from_canonical_u64(7));
    }

    #[test]
    fn test_multiplication() {
        let a = ScalarField::from_canonical_u64(6);
        let b = ScalarField::from_canonical_u64(7);
        let c = a * b;
        assert_eq!(c, ScalarField::from_canonical_u64(42));
    }

    #[test]
    fn test_negation() {
        let a = ScalarField::from_canonical_u64(5);
        let b = -a;
        assert_eq!(a + b, ScalarField::ZERO);
    }

    #[test]
    fn test_inverse() {
        let a = ScalarField::from_canonical_u64(5);
        let a_inv = a.inverse();
        assert_eq!(a * a_inv, ScalarField::ONE);
    }
}
