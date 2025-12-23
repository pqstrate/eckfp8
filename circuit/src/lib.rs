//! # Zero-Knowledge Schnorr Signature Verification Circuit
//!
//! Production-grade ZK-SNARK circuits for verifying Schnorr signatures over the KoalaBear
//! elliptic curve using the Plonky3 STARK framework.
//!
//! ## Overview
//!
//! This crate provides a complete ZK-SNARK implementation for proving knowledge of
//! valid Schnorr signatures without revealing the signature itself. Built on Plonky3's
//! Algebraic Intermediate Representation (AIR) and STARK proof system.
//!
//! **Key Features**:
//! - **Native Field Operations**: KoalaBear field arithmetic for elliptic curve operations
//! - **Non-Native Scalar Arithmetic**: Efficient 9×28-bit limb representation for scalars
//! - **Optimized Constraints**: 49-column trace with degree-4 constraints
//! - **Preprocessing Support**: Generator multiplication uses preprocessed trace
//! - **Modular Design**: Composable sub-circuits (scalar mul, point ops, hashing)
//!
//! ## Quick Start
//!
//! ### Basic Circuit Usage
//!
//! ```rust,ignore
//! use circuit::{SignatureWitness, SchnorrAir, build_schnorr_trace};
//! use schnorr::SigningKey;
//! use p3_baby_bear::BabyBear;
//! use p3_field::PrimeCharacteristicRing;
//! use rand::rng;
//!
//! let mut rng = rng();
//!
//! // Generate signature
//! let signing_key = SigningKey::random(&mut rng);
//! let verifying_key = signing_key.verifying_key();
//! let message = vec![BabyBear::from_u32(42)];
//! let signature = signing_key.sign(&mut rng, &message).unwrap();
//!
//! // Create witness
//! let witness = SignatureWitness::new(&signature, &verifying_key, &message);
//!
//! // Build constraint trace (height must be power of 2)
//! let trace = build_schnorr_trace(&witness, 256);
//!
//! // Create AIR
//! let air = SchnorrAir::new(256);
//!
//! // Generate STARK proof (see examples/prove_verify_signature.rs)
//! ```
//!
//! ### Complete STARK Proof Example
//!
//! See [`examples/prove_verify_signature.rs`](../examples/prove_verify_signature.rs)
//! for a full proof generation and verification workflow including:
//! - Trace generation
//! - AIR constraint evaluation
//! - FRI commitment and proof generation
//! - Proof verification
//! - Performance metrics
//!
//! ## Verification Equation
//!
//! The circuit implements the Schnorr verification equation:
//!
//! ```text
//! G × s - pk × e = R
//! ```
//!
//! Where:
//! - `G`: Generator point (fixed, preprocessed)
//! - `s`: Signature response scalar
//! - `R`: Signature commitment point
//! - `pk`: Public key (verifying key)
//! - `e`: Fiat-Shamir challenge `e = Poseidon2(R || pk || msg)`
//!
//! This is computed as:
//! 1. `P₁ = G × s` (fixed-base scalar multiplication, preprocessed)
//! 2. `P₂ = pk × e` (variable-base scalar multiplication)
//! 3. `P₃ = P₁ - P₂` (point subtraction)
//! 4. Assert `P₃ = R` (public input check)
//!
//! ## Key Types
//!
//! ### [`SignatureWitness`]
//!
//! Container for all witness data needed to build the circuit trace.
//!
//! **Fields**:
//! - `r: CircuitPoint` - Signature commitment R
//! - `s: CircuitScalar` - Signature response scalar
//! - `public_key: CircuitPoint` - Verifying key
//! - `message: Vec<BabyBear>` - Message field elements
//! - `challenge: CircuitScalar` - Fiat-Shamir challenge e
//!
//! **Methods**:
//! - `new(signature, verifying_key, message)` - Construct from signature components
//!
//! ### [`CircuitPoint`]
//!
//! Elliptic curve point representation for circuit constraints.
//!
//! **Representation**: Each coordinate (x, y) is an Fp8 element = 8 KoalaBear field elements
//!
//! **Methods**:
//! - `from_affine(point)` - Convert from curve::Affine
//! - Total size: 16 KoalaBear elements (2 coordinates × 8 elements each)
//!
//! ### [`CircuitScalar`]
//!
//! Scalar field element representation for non-native arithmetic.
//!
//! **Representation**: 9 limbs × 28 bits = 252 bits (scalar field size)
//!
//! **Constants**:
//! - `SCALAR_LIMBS = 9` - Number of limbs
//! - `LIMB_BITS = 28` - Bits per limb
//!
//! **Methods**:
//! - `from_scalar_field(scalar)` - Convert from curve::ScalarField
//!
//! ## Circuit Components
//!
//! ### [`SchnorrAir`] - Main Circuit
//!
//! Complete Schnorr verification constraint system.
//!
//! **Trace Structure**:
//! - Main trace: 49 columns (point coordinates, scalars, intermediate values)
//! - Preprocessed trace: Generator multiplication table
//! - Public inputs: 32 KoalaBear elements (pk: 16 + R: 16)
//!
//! **Usage**:
//! ```rust,ignore
//! use circuit::{SchnorrAir, build_schnorr_trace, SignatureWitness};
//! # use schnorr::SigningKey;
//! # use p3_baby_bear::BabyBear;
//! # use rand::rng;
//! # let mut rng = rng();
//! # let signing_key = SigningKey::random(&mut rng);
//! # let message = vec![BabyBear::from_u32(42)];
//! # let signature = signing_key.sign(&mut rng, &message).unwrap();
//! # let verifying_key = signing_key.verifying_key();
//! # let witness = SignatureWitness::new(&signature, &verifying_key, &message);
//!
//! let air = SchnorrAir::new(256); // trace height (power of 2)
//! let trace = build_schnorr_trace(&witness, 256);
//!
//! // Integrate with Plonky3 prover
//! ```
//!
//! ### [`ScalarMulAir`] - Scalar Multiplication Sub-Circuit
//!
//! Implements scalar multiplication `k × P` using binary double-and-add.
//!
//! **Columns**: `NUM_COLUMNS` (defined in module)
//! - Accumulator point (ACC_X_START, ACC_Y_START): 16 limbs
//! - Base point (BASE_X_START, BASE_Y_START): 16 limbs
//! - Intermediate computation columns
//!
//! **Functions**:
//! - `build_scalar_mul_trace(scalar, base, height)` - Variable-base multiplication
//! - `build_generator_mul_trace(scalar, height)` - Fixed-base multiplication (preprocessed)
//!
//! ### [`Poseidon2HashAir`] - Hash Function Sub-Circuit
//!
//! Verifies Poseidon2 sponge hash computation.
//!
//! **Configuration**:
//! - `POSEIDON2_WIDTH = 16` - Sponge width
//! - `POSEIDON2_RATE = 8` - Absorption rate
//! - `POSEIDON2_OUT = 8` - Output size
//!
//! **Function**:
//! - `build_poseidon2_hash_trace(input, output)` - Generate hash constraint trace
//!
//! ## Circuit Architecture
//!
//! ### Native Field Operations
//!
//! All elliptic curve arithmetic uses **native KoalaBear field operations**:
//! - Point addition: Complete formulas handling all edge cases
//! - Point doubling: Efficient tangent-based formulas
//! - Fp8 arithmetic: Binomial extension field operations
//!
//! Benefits:
//! - Minimal constraint overhead
//! - Direct representation (no limb decomposition for coordinates)
//! - Efficient Plonky3 AIR constraints
//!
//! ### Non-Native Scalar Arithmetic
//!
//! Scalar field operations use **9×28-bit limb decomposition**:
//!
//! ```text
//! ScalarField (252 bits) → [l₀, l₁, ..., l₈] where each lᵢ ∈ [0, 2²⁸)
//! ```
//!
//! Operations constrained:
//! - Bit decomposition: Scalar → 252 bits
//! - Range checks: Each limb < 2²⁸
//! - Carry propagation: Multi-limb arithmetic
//!
//! Benefits:
//! - Efficient representation in KoalaBear constraints
//! - 28 bits per limb fits comfortably in 31-bit field
//! - Minimal overhead for range checking
//!
//! ## Trace Structure
//!
//! ### Main Trace (49 columns × height rows)
//!
//! **Layout**:
//! - Columns 0-15: Generator multiplication result (P₁ = G × s)
//! - Columns 16-31: Public key multiplication result (P₂ = pk × e)
//! - Columns 32-47: Point subtraction result (P₃ = P₁ - P₂)
//! - Column 48: Control/selector flags
//!
//! **Height**: Must be power of 2 (typically 256-4096)
//! - Larger heights: More computation, better amortization
//! - Smaller heights: Less memory, faster proving for simple cases
//!
//! ### Preprocessed Trace
//!
//! Generator multiplication table precomputed at circuit setup:
//! - Fixed generator G (known at compile time)
//! - Precomputed multiples: [G, 2G, 4G, 8G, ...]
//! - Reduces online computation during proving
//!
//! ### Public Inputs (32 elements)
//!
//! **Layout**:
//! - Elements 0-15: Public key (pk.x: 8, pk.y: 8)
//! - Elements 16-31: Commitment point (R.x: 8, R.y: 8)
//!
//! Verifier checks these match the claimed values.
//!
//! Run benchmarks: `cargo bench -p circuit`
//!
//! ## Constraint Analysis
//!
//! ### Constraint Degree
//!
//! - **Maximum degree**: 4 (quadratic constraints over extension field)
//! - Most constraints: Degree 2-3
//! - Multiplication constraints: Degree 2
//! - Point addition formulas: Degree 3-4
//!
//! ### Constraint Count
//!
//! Per row:
//! - Scalar multiplication: ~20-30 constraints
//! - Point operations: ~15-20 constraints
//! - Limb range checks: ~10-15 constraints
//!
//! Total: ~50-60 constraints per trace row
//!
//! ## Optimization Techniques
//!
//! ### Implemented
//!
//! 1. **Fixed-Base Preprocessing**: Generator multiplication table precomputed
//! 2. **Efficient Limb Width**: 28-bit limbs minimize range check overhead
//! 3. **Native Field Ops**: Elliptic curve ops use native KoalaBear arithmetic
//! 4. **Complete Formulas**: Unified addition avoiding edge-case branching
//!
//! ### Future Optimizations
//!
//! 1. **Windowed Multiplication**: 4-bit windows reduce doublings
//! 2. **Shamir's Trick**: Simultaneous double scalar multiplication
//! 3. **Batch Verification**: Amortize constraints across multiple signatures
//! 4. **Lookup Tables**: Replace range checks with table lookups
//!
//! **Production Use**: Conduct independent security review before production deployment.
//!
//! ## Examples
//!
//! ### Generating a STARK Proof
//!
//! Complete example in [`examples/prove_verify_signature.rs`](../examples/prove_verify_signature.rs):
//!
//! ```rust,ignore
//! use circuit::{SignatureWitness, SchnorrAir, build_schnorr_trace};
//! use p3_uni_stark::{prove, verify, StarkConfig};
//! // ... (full imports in example file)
//!
//! // 1. Create signature witness
//! let witness = SignatureWitness::new(&signature, &verifying_key, &message);
//!
//! // 2. Build trace
//! let trace = build_schnorr_trace(&witness, 256);
//!
//! // 3. Create AIR
//! let air = SchnorrAir::new(256);
//!
//! // 4. Configure STARK
//! let config = StarkConfig::default();
//!
//! // 5. Generate proof
//! let proof = prove(&air, &trace, &config);
//!
//! // 6. Verify proof
//! let is_valid = verify(&air, &proof, &config);
//! assert!(is_valid);
//! ```
//!
//! ## Integration with Schnorr Crate
//!
//! This crate tightly integrates with the `schnorr` crate:
//!
//! ```rust,ignore
//! use schnorr::{SigningKey, Signature, VerifyingKey};
//! use circuit::SignatureWitness;
//! use rand::rng;
//! use p3_baby_bear::BabyBear;
//!
//! // Generate signature using schnorr crate
//! let mut rng = rng();
//! let signing_key = SigningKey::random(&mut rng);
//! let message = vec![BabyBear::from_u32(42)];
//! let signature = signing_key.sign(&mut rng, &message).unwrap();
//!
//! // Create circuit witness
//! let witness = SignatureWitness::new(
//!     &signature,
//!     &signing_key.verifying_key(),
//!     &message
//! );
//!
//! // Now use witness in circuit...
//! ```
//!
//! ## References
//!
//! - Plonky3 framework: <https://github.com/Plonky3/Plonky3>
//! - STARK proofs: <https://eprint.iacr.org/2018/046>
//! - FRI protocol: <https://eccc.weizmann.ac.il/report/2017/134/>
//! - Poseidon2: <https://eprint.iacr.org/2023/323>
//! - Non-native field arithmetic in circuits: <https://eprint.iacr.org/2019/458>

#[deny(missing_docs)]
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
