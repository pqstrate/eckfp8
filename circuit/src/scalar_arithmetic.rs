//! Scalar field arithmetic for circuit implementation.
//!
//! This module provides a scalar representation in KoalaBear limbs for circuit use.
//! Since the scalar field is larger than the KoalaBear field, we represent scalar
//! elements as multiple KoalaBear field elements.

use curve::{KoalaBear, ScalarField};
use p3_field::{PrimeCharacteristicRing, PrimeField32};

/// Number of KoalaBear limbs needed to represent a scalar field element.
/// The scalar field is ~252 bits, KoalaBear is 31 bits, so we need 9 limbs.
pub const SCALAR_LIMBS: usize = 9;

/// Number of bits per limb (using 28 bits for easier arithmetic)
pub const LIMB_BITS: u32 = 28;

/// Scalar field element represented as KoalaBear limbs for circuit operations.
///
/// This representation breaks down a 252-bit scalar into 9 limbs of 28 bits each.
/// This allows us to perform non-native arithmetic using native KoalaBear operations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CircuitScalar {
    /// Limbs in little-endian order, each limb should be < 2^28
    pub limbs: [KoalaBear; SCALAR_LIMBS],
}

impl CircuitScalar {
    /// Create a scalar from a ScalarField element
    pub fn from_scalar_field(scalar: ScalarField) -> Self {
        let canonical = scalar.to_canonical_u64_vec();
        let mut limbs = [KoalaBear::ZERO; SCALAR_LIMBS];

        // Convert 4 u64 limbs to 9 28-bit limbs
        // Process bit by bit to avoid overflow
        let mut bit_position = 0u32;

        for limb_idx in 0..SCALAR_LIMBS {
            if bit_position >= 256 {
                break;
            }

            let word_idx = (bit_position / 64) as usize;
            let bit_offset = bit_position % 64;

            if word_idx < 4 {
                // Extract 28 bits starting at bit_position
                let mut limb_value =
                    (canonical[word_idx] >> bit_offset) & ((1u64 << LIMB_BITS) - 1);

                // If we need more bits from the next word
                if bit_offset + LIMB_BITS > 64 && word_idx + 1 < 4 {
                    let bits_from_next = bit_offset + LIMB_BITS - 64;
                    let next_bits = canonical[word_idx + 1] & ((1u64 << bits_from_next) - 1);
                    limb_value |= next_bits << (LIMB_BITS - bits_from_next);
                }

                limbs[limb_idx] = KoalaBear::from_u32(limb_value as u32);
            }

            bit_position += LIMB_BITS;
        }

        Self { limbs }
    }

    /// Convert back to ScalarField
    pub fn to_scalar_field(&self) -> ScalarField {
        let mut result = [0u64; 4];

        // Accumulate bits limb by limb
        let mut bit_position = 0u32;

        for limb in &self.limbs {
            let limb_value = limb.as_canonical_u32() as u64;
            let word_idx = (bit_position / 64) as usize;
            let bit_offset = bit_position % 64;

            if word_idx < 4 {
                result[word_idx] |= limb_value << bit_offset;

                // Handle overflow to next word
                if bit_offset + LIMB_BITS > 64 && word_idx + 1 < 4 {
                    let overflow_bits = bit_offset + LIMB_BITS - 64;
                    result[word_idx + 1] |= limb_value >> (LIMB_BITS - overflow_bits);
                }
            }

            bit_position += LIMB_BITS;
        }

        ScalarField::from_canonical_limbs(result)
    }
}
