//! Schnorr signature scheme over the KoalaBear elliptic curve.
//!
//! This library implements a Schnorr signature scheme using:
//! - The KoalaBear elliptic curve (Fp8 curve)
//! - Poseidon2-BabyBear hash function for the Fiat-Shamir challenge
//! - BabyBear field elements for message encoding
//!
//! # Overview
//!
//! The Schnorr signature scheme is a digital signature scheme that provides:
//! - Unforgeability: Only the holder of the secret key can produce valid signatures
//! - Non-repudiation: The signer cannot deny having signed a message
//! - Verification efficiency: Signatures can be verified efficiently
//!
//! # Example
//!
//! ```
//! use schnorr::{SigningKey, VerifyingKey, Signature};
//! use p3_baby_bear::BabyBear;
//! use p3_field::PrimeCharacteristicRing;
//! use rand::thread_rng;
//!
//! // Generate a random signing key
//! let mut rng = thread_rng();
//! let signing_key = SigningKey::random(&mut rng);
//!
//! // Derive the corresponding verifying key
//! let verifying_key = signing_key.verifying_key();
//!
//! // Create a message (using BabyBear field elements)
//! let message = [
//!     BabyBear::from_u32(1),
//!     BabyBear::from_u32(2),
//!     BabyBear::from_u32(3),
//! ];
//!
//! // Sign the message
//! let signature = signing_key.sign(&mut rng, &message).expect("signing failed");
//!
//! // Verify the signature
//! let is_valid = verifying_key.verify(&message, &signature).expect("verification failed");
//! assert!(is_valid);
//! ```
//!
//! # Security Considerations
//!
//! - Always use a cryptographically secure random number generator (CSRNG)
//! - Each signature must use a fresh random nonce
//! - Protect the signing key from unauthorized access
//! - Messages are encoded as BabyBear field elements

mod constants;
mod errors;
mod keys;
mod signatures;

#[cfg(test)]
mod tests;

pub use errors::SchnorrError;
pub use keys::{SigningKey, VerifyingKey};
pub use signatures::Signature;
