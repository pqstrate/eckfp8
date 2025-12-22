//! Circuit implementation for Schnorr signature verification.
//!
//! This crate provides a ZK-SNARK friendly circuit for verifying Schnorr signatures
//! over the KoalaBear elliptic curve. The circuit uses native field operations
//! (KoalaBear field) for all elliptic curve operations, with special handling for
//! scalar field arithmetic using non-native field techniques.
//!
//! # Overview
//!
//! The signature verification circuit implements the Schnorr verification equation:
//! ```text
//! G * s - pk * e == R
//! ```
//! where:
//! - `G` is the generator point
//! - `s` is the signature response scalar
//! - `R` is the signature commitment point
//! - `pk` is the public key
//! - `e` is the Fiat-Shamir challenge: `e = H(R || pk || msg)`
//!
//! # Example
//!
//! ```rust
//! use circuit::{SignatureWitness, SchnorrAir, build_schnorr_trace};
//! use schnorr::SigningKey;
//! use p3_baby_bear::BabyBear;
//! use p3_field::PrimeCharacteristicRing;
//! use p3_matrix::Matrix;
//! use rand::rng;
//!
//! // Generate a key pair
//! let mut rng = rng();
//! let signing_key = SigningKey::random(&mut rng);
//! let verifying_key = signing_key.verifying_key();
//!
//! // Create and sign a message
//! let message = vec![BabyBear::from_u32(1), BabyBear::from_u32(2)];
//! let signature = signing_key.sign(&mut rng, &message).unwrap();
//!
//! // Create witness and verify in circuit
//! let witness = SignatureWitness::new(&signature, &verifying_key, &message).unwrap();
//! let trace = build_schnorr_trace(&witness);
//! let _air = SchnorrAir::new(trace.trace.height());
//! assert!(trace.trace.height().is_power_of_two());
//! ```
//!
//! # Circuit Design
//!
//! ## Native Field Operations
//!
//! All elliptic curve point operations (addition, doubling) are performed using
//! native KoalaBear field arithmetic. Points are represented in affine coordinates
//! where each coordinate is a degree-8 extension field element.
//!
//! ## Non-Native Scalar Arithmetic
//!
//! Since the scalar field is different from the base field, scalar field operations
//! require non-native arithmetic. We represent scalar field elements as multiple
//! KoalaBear limbs (9 limbs of 28 bits each) and implement limb-based arithmetic.
//!
//! ## Optimization Opportunities
//!
//! - **Fixed-base scalar multiplication**: Since G is fixed, we can use precomputed tables
//! - **Windowed methods**: For variable-base scalar multiplication, windowing reduces doublings
//! - **Shamir's trick**: Double scalar multiplication (G*s + pk*(-e)) can be optimized
//! - **Batch verification**: Multiple signatures can be verified more efficiently together

mod point_ops;
pub mod poseidon2_hash_air;
mod scalar_arithmetic;
pub mod scalar_mul_air;
pub mod schnorr_air;
mod signature_witness;

pub use point_ops::{scalar_to_bits, CircuitPoint};
pub use poseidon2_hash_air::{
    build_poseidon2_hash_trace, Poseidon2HashAir, Poseidon2HashTrace, POSEIDON2_INPUT_LEN,
    POSEIDON2_NUM_PERMS, POSEIDON2_OUT, POSEIDON2_PACKED_LIMBS, POSEIDON2_RATE, POSEIDON2_WIDTH,
};
pub use scalar_arithmetic::{CircuitScalar, LIMB_BITS, SCALAR_LIMBS};
pub use signature_witness::SignatureWitness;

// Re-export commonly used types
pub use curve::{Affine, BaseField, KoalaBear, ScalarField};
pub use scalar_mul_air::{
    build_generator_mul_trace, build_scalar_mul_trace, ScalarMulAir, ScalarMulTrace, ACC_X_START,
    ACC_Y_START, BASE_X_START, BASE_Y_START, NUM_COLUMNS as SCALAR_MUL_NUM_COLUMNS,
    PUBLIC_BASE_LIMBS, PUBLIC_OUT_LIMBS,
};
pub use schnorr::{Signature, SigningKey, VerifyingKey};
pub use schnorr_air::{build_schnorr_trace, SchnorrAir, SchnorrTrace, SCHNORR_COLUMNS};
