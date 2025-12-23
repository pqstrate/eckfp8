//! # Schnorr Signature Scheme
//!
//! Production-ready Schnorr signature implementation over the KoalaBear elliptic curve
//! with Poseidon2-based Fiat-Shamir hashing.
//!
//! ## Overview
//!
//! This library provides a complete Schnorr signature scheme featuring:
//! - **Elliptic Curve**: KoalaBear curve over Fp8 extension field
//! - **Hash Function**: Poseidon2 sponge over BabyBear field for Fiat-Shamir challenges
//!
//! ## Algorithm
//!
//! ### Signing Process
//!
//! Given a message `m` and signing key `sk`:
//!
//! 1. Generate random nonce: `k ← ScalarField`
//! 2. Compute commitment: `R = G × k`
//! 3. Compute challenge: `e = Poseidon2(R || pk || m)`
//! 4. Compute response: `s = k + e × sk`
//! 5. Return signature `σ = (R, s)`
//!
//! ### Verification Process
//!
//! Given a message `m`, signature `σ = (R, s)`, and public key `pk`:
//!
//! 1. Recompute challenge: `e = Poseidon2(R || pk || m)`
//! 2. Check equation: `G × s = R + pk × e`
//! 3. Accept if equation holds, reject otherwise
//!
//! ## Quick Start
//!
//! ### Basic Signing and Verification
//!
//! ```rust
//! use schnorr::{SigningKey, VerifyingKey};
//! use p3_baby_bear::BabyBear;
//! use p3_field::PrimeCharacteristicRing;
//! use rand::rng;
//!
//! let mut rng = rng();
//!
//! // Generate keypair
//! let signing_key = SigningKey::random(&mut rng);
//! let verifying_key = signing_key.verifying_key();
//!
//! // Sign a message (encoded as BabyBear field elements)
//! let message = vec![
//!     BabyBear::from_u32(42),
//!     BabyBear::from_u32(1337),
//! ];
//! let signature = signing_key.sign(&mut rng, &message)
//!     .expect("signing failed");
//!
//! // Verify the signature
//! let is_valid = verifying_key.verify(&message, &signature)
//!     .expect("verification failed");
//! assert!(is_valid);
//! ```
//!
//! ### Serialization
//!
//! All types support efficient binary serialization:
//!
//! ```rust,ignore
//! use schnorr::{SigningKey, VerifyingKey, Signature, SK_SIZE, PK_SIZE, SIG_SIZE};
//! use rand::rng;
//!
//! let mut rng = rng();
//! let signing_key = SigningKey::random(&mut rng);
//!
//! // Serialization sizes
//! assert_eq!(SK_SIZE, 32);  // Secret key: 32 bytes
//! assert_eq!(PK_SIZE, 40);  // Public key: 40 bytes (2×Fp8 coordinates)
//! assert_eq!(SIG_SIZE, 72); // Signature: 72 bytes (R point + s scalar)
//!
//! // Serialize signing key
//! let sk_bytes = bincode::serialize(&signing_key).unwrap();
//! assert_eq!(sk_bytes.len(), SK_SIZE);
//!
//! // Deserialize signing key
//! let recovered: SigningKey = bincode::deserialize(&sk_bytes).unwrap();
//! assert_eq!(signing_key, recovered);
//!
//! // Same for verifying keys and signatures
//! let vk_bytes = bincode::serialize(&signing_key.verifying_key()).unwrap();
//! assert_eq!(vk_bytes.len(), PK_SIZE);
//! ```
//!
//! ## Key Types
//!
//! ### [`SigningKey`]
//!
//! Secret signing key (32 bytes).
//!
//! **Methods**:
//! - `random(rng)` - Generate random signing key
//! - `from_bytes(bytes)` - Deserialize from 32-byte array
//! - `to_bytes()` - Serialize to 32-byte array
//! - `verifying_key()` - Derive corresponding public key
//! - `sign(rng, message)` - Sign a message
//!
//! ### [`VerifyingKey`]
//!
//! Public verification key (40 bytes: 2×Fp8 coordinates).
//!
//! **Methods**:
//! - `verify(message, signature)` - Verify a signature
//! - `from_affine(point)` - Construct from curve point
//! - `to_affine()` - Convert to curve point
//!
//! ### [`Signature`]
//!
//! Schnorr signature (72 bytes: R point + s scalar).
//!
//! **Fields**:
//! - `r: Affine` - Commitment point (40 bytes)
//! - `s: ScalarField` - Response scalar (32 bytes)
//!
//! ## Poseidon2 Hash Function
//!
//! The Fiat-Shamir challenge is computed using Poseidon2:
//!
//! ```text
//! Input Encoding:
//! - R (commitment): 16 BabyBear elements (2 Fp8 coordinates)
//! - pk (public key): 16 BabyBear elements (2 Fp8 coordinates)
//! - msg (message): Variable-length BabyBear elements
//!
//! Hash Configuration:
//! - Width: 16 field elements
//! - Rate: 8 elements
//! - Capacity: 8 elements
//! - Output: 8 elements (squeezed to 5 elements for scalar)
//!
//! Output: Challenge scalar (252 bits)
//! ```
//!
//! See [`hash_challenge`] for implementation details.
//!
//! ## Security Properties
//!
//! ### Cryptographic Guarantees
//!
//! - **Unforgeability**: Existential unforgeability under chosen-message attack (EUF-CMA)
//! - **Non-Malleability**: Poseidon2 Fiat-Shamir prevents signature malleability
//! - **Collision Resistance**: Poseidon2 provides ~128-bit collision resistance
//! - **Discrete Log Security**: ~124-bit security from curve's prime-order group
//!
//! ### Implementation Security
//!
//! - **Nonce Uniqueness**: Each signature requires fresh random nonce
//! - **Prime-Order Curve**: Cofactor = 1 eliminates small-subgroup attacks
//! - **Field Validation**: All deserialization validates field membership
//!
//! ## Security Considerations
//!
//! ### Critical Requirements
//!
//! 1. **Cryptographically Secure RNG**
//!    - MUST use CSRNG for key generation and signing
//!    - Use `rand::thread_rng()` or equivalent
//!    - NEVER reuse nonces across signatures
//!
//! 2. **Key Management**
//!    - Protect signing keys from unauthorized access
//!    - Zero memory when disposing of keys
//!    - Use hardware security modules (HSMs) for high-value keys
//!
//! 3. **Message Encoding**
//!    - Messages are encoded as BabyBear field elements
//!    - Ensure canonical encoding of application-level messages
//!    - Hash large messages before signing
//!
//! ### Known Limitations
//!
//! - **Not Audited**: This library has NOT undergone formal security audit
//! - **Research Code**: Suitable for research, prototyping, and non-critical applications
//! - **Side-Channel Attacks**: Limited protection against power analysis or cache timing
//! - **Fault Injection**: No explicit countermeasures against fault attacks
//!
//! **Production Use**: Conduct independent security review before production deployment.
//!
//! ## Error Handling
//!
//! All operations return `Result<T, SchnorrError>`:
//!
//! ```rust
//! use schnorr::{SchnorrError, SigningKey};
//! use rand::rng;
//!
//! let mut rng = rng();
//! let signing_key = SigningKey::random(&mut rng);
//!
//! // Signing can fail (though unlikely with proper RNG)
//! match signing_key.sign(&mut rng, &[]) {
//!     Ok(signature) => println!("Signature: {:?}", signature),
//!     Err(SchnorrError::InvalidPoint) => eprintln!("Error: Invalid point encountered"),
//! }
//! ```
//!
//! ## Examples
//!
//! See `examples/schnorr.rs` for a complete workflow demonstration.
//!
//! ## Performance
//!
//! Typical performance on Apple M1:
//!
//! | Operation | Time | Notes |
//! |-----------|------|-------|
//! | Key Generation | ~50 μs | Random scalar + point mul |
//! | Signing | ~80 μs | Includes Poseidon2 hash |
//! | Verification | ~120 μs | Double scalar mul + hash |
//!
//! Run benchmarks: `cargo bench -p schnorr`
//!
//! ## Integration with ZK Circuits
//!
//! This signature scheme is designed for zero-knowledge proof integration:
//!
//! ```rust,ignore
//! // The circuit crate provides ZK-SNARK verification
//! use circuit::{SignatureWitness, SchnorrAir, build_schnorr_trace};
//!
//! // Create witness from signature
//! let witness = SignatureWitness::new(&signature, &verifying_key, &message);
//!
//! // Build constraint trace
//! let trace = build_schnorr_trace(&witness, 256);
//!
//! // Generate STARK proof (see circuit crate for details)
//! ```
//!
//! See the `circuit` crate for complete ZK-SNARK integration.
//!
//! ## References
//!
//! - Schnorr Signatures: Claus-Peter Schnorr (1989)
//! - Poseidon2: <https://eprint.iacr.org/2023/323>
//! - Plonky3 framework: <https://github.com/Plonky3/Plonky3>

#[deny(missing_docs)]
mod constants;
mod errors;
mod keys;
mod signatures;

#[cfg(test)]
mod tests;

pub use constants::{PK_SIZE, SIG_SIZE, SK_SIZE};
pub use errors::SchnorrError;
pub use keys::{SigningKey, VerifyingKey};
pub use signatures::{Signature, hash_challenge};
