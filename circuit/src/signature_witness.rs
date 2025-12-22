//! Signature witness data for Schnorr verification.

use crate::point_ops::CircuitPoint;
use crate::scalar_arithmetic::CircuitScalar;
use p3_baby_bear::BabyBear;
use schnorr::{hash_challenge, Signature, VerifyingKey};

/// Witness data for the signature verification circuit.
///
/// This contains all the data needed to verify a signature within a circuit.
/// In a ZK-SNARK, some of these fields would be public inputs (like the message
/// and public key) while others would be private witnesses (like the signature).
#[derive(Clone, Debug)]
pub struct SignatureWitness {
    /// The signature commitment point R
    pub r: CircuitPoint,
    /// The signature scalar s
    pub s: CircuitScalar,
    /// The public key (verifying key)
    pub public_key: CircuitPoint,
    /// The message (as BabyBear field elements)
    pub message: Vec<BabyBear>,
    /// The challenge e = H(R || pk || msg)
    pub challenge: CircuitScalar,
}

impl SignatureWitness {
    /// Create a witness from a signature, public key, and message
    pub fn new(
        signature: &Signature,
        public_key: &VerifyingKey,
        message: &[BabyBear],
    ) -> Result<Self, String> {
        let pk_point = public_key.as_affine();

        // Compute the challenge
        let challenge = hash_challenge(&signature.r, &pk_point, message)
            .map_err(|e| format!("Failed to compute challenge: {:?}", e))?;

        Ok(Self {
            r: CircuitPoint::from_affine(&signature.r),
            s: CircuitScalar::from_scalar_field(signature.s),
            public_key: CircuitPoint::from_affine(&pk_point),
            message: message.to_vec(),
            challenge: CircuitScalar::from_scalar_field(challenge),
        })
    }
}
