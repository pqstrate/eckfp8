//! Signature types and challenge hashing for the Schnorr signature scheme.

use curve::{Affine, KoalaBear, ScalarField};
use p3_baby_bear::{BabyBear, Poseidon2BabyBear, default_babybear_poseidon2_16};
use p3_field::{PrimeCharacteristicRing, PrimeField32};
use p3_symmetric::{CryptographicHasher, PaddingFreeSponge};
use serde::{Deserialize, Serialize};

use crate::constants::{POSEIDON2_OUT, POSEIDON2_RATE, POSEIDON2_WIDTH};
use crate::errors::SchnorrError;

/// A Schnorr signature consisting of a curve point and a scalar.
///
/// The signature is a pair `(R, s)` where:
/// - `R` is a point on the KoalaBear elliptic curve (the commitment)
/// - `s` is a scalar in the scalar field (the response)
///
/// # Structure
///
/// The signature satisfies the verification equation: `G * s == R + pk * e`
/// where `e = H(R || pk || msg)` is the Fiat-Shamir challenge.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Signature {
    /// The commitment point R = G * k, where k is the signing nonce
    pub r: Affine,
    /// The response scalar s = k + e * sk, where sk is the signing key
    pub s: ScalarField,
}

/// Computes the Fiat-Shamir challenge for the Schnorr signature scheme.
///
/// The challenge is computed as `e = H(R || pk || msg)` using the Poseidon2
/// hash function over the KoalaBear field.
///
/// # Arguments
///
/// * `r` - The commitment point R from the signature
/// * `pk` - The public verifying key
/// * `msg` - The message being signed/verified, encoded as KoalaBear field elements
///
/// # Returns
///
/// Returns a `Result` containing:
/// - `Ok(e)` where `e` is the challenge scalar
/// - `Err(SchnorrError::InvalidPoint)` if either `r` or `pk` is the point at infinity
///
/// # Implementation Details
///
/// 1. Points are encoded as 16 KoalaBear field elements (8 for x-coordinate, 8 for y-coordinate)
/// 2. The input is `R || pk || msg` concatenated
/// 3. Poseidon2 with width 16, rate 8, and output 8 is used for hashing
/// 4. The first 5 digest elements are packed into a scalar field element
pub(crate) fn hash_challenge(
    r: &Affine,
    pk: &Affine,
    msg: &[BabyBear],
) -> Result<ScalarField, SchnorrError> {
    if r.is_infinity() || pk.is_infinity() {
        return Err(SchnorrError::InvalidPoint);
    }

    let mut input = Vec::with_capacity(msg.len() + 32);
    input.extend_from_slice(&encode_point(r));
    input.extend_from_slice(&encode_point(pk));
    input.extend_from_slice(msg);

    let perm = default_babybear_poseidon2_16();
    let sponge = PaddingFreeSponge::<
        Poseidon2BabyBear<POSEIDON2_WIDTH>,
        POSEIDON2_WIDTH,
        POSEIDON2_RATE,
        POSEIDON2_OUT,
    >::new(perm);

    let digest = sponge.hash_iter(input);
    let d0 = digest[0].as_canonical_u32() as u64;
    let d1 = digest[1].as_canonical_u32() as u64;
    let d2 = digest[2].as_canonical_u32() as u64;
    let d3 = digest[3].as_canonical_u32() as u64;
    let d4 = digest[4].as_canonical_u32() as u64;

    let limb0 = d0 | (d1 << 31);
    let limb1 = d2 | (d3 << 31);
    let limb2 = d4;

    Ok(ScalarField::from_canonical_limbs([limb0, limb1, limb2, 0]))
}

/// Encodes an elliptic curve point as an array of KoalaBear field elements.
///
/// The KoalaBear curve is defined over an Fp8 extension field, where each
/// coordinate consists of 8 base field elements. This function extracts
/// those elements and converts them to KoalaBear field elements for use in Poseidon2.
///
/// # Arguments
///
/// * `point` - An affine point on the KoalaBear curve
///
/// # Returns
///
/// An array of 16 KoalaBear field elements:
/// - Elements 0-7: x-coordinate coefficients
/// - Elements 8-15: y-coordinate coefficients
fn encode_point(point: &Affine) -> [BabyBear; 16] {
    let x_coeffs: [KoalaBear; 8] = unsafe { core::mem::transmute(point.x) };
    let y_coeffs: [KoalaBear; 8] = unsafe { core::mem::transmute(point.y) };

    let mut out = [BabyBear::ZERO; 16];
    for (i, coeff) in x_coeffs.into_iter().enumerate() {
        out[i] = BabyBear::from_u32(coeff.as_canonical_u32());
    }
    for (i, coeff) in y_coeffs.into_iter().enumerate() {
        out[i + 8] = BabyBear::from_u32(coeff.as_canonical_u32());
    }

    out
}
