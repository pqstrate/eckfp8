# eckfp8

[![Rust](https://img.shields.io/badge/rust-1.70%2B-blue.svg)](https://www.rust-lang.org)
[![License](https://img.shields.io/badge/license-MIT%2FApache--2.0-blue.svg)](LICENSE)

Production-ready Rust implementation of the KoalaBear elliptic curve over Fp8 extension field, featuring Schnorr signatures and zero-knowledge proof circuits for signature verification.

## Overview

**eckfp8** is a comprehensive cryptographic library providing:

- **Elliptic Curve Cryptography**: KoalaBear curve defined over an 8th-degree extension field (Fp8)
- **Schnorr Signatures**: Secure signature scheme with Poseidon2-based Fiat-Shamir hashing
- **ZK-SNARK Circuits**: Plonky3-based circuits for zero-knowledge signature verification
- **Production-Grade**: Comprehensive testing and optimized field arithmetic

### Use Cases

- **Zero-Knowledge Applications**: Prove knowledge of signatures without revealing them
- **Privacy-Preserving Systems**: Anonymous authentication and credential systems
- **Blockchain & Cryptographic Protocols**: ZK-friendly signature schemes
- **Research & Education**: Learning modern ZK-SNARK constructions

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Architecture](#architecture)
- [Crates](#crates)
  - [curve](#curve-crate)
  - [schnorr](#schnorr-crate)
  - [circuit](#circuit-crate)
- [Examples](#examples)
- [Security](#security)
- [Performance](#performance)
- [Development](#development)
- [License](#license)

## Installation

Add the desired crates to your `Cargo.toml`:

```toml
[dependencies]
# For elliptic curve operations
eckfp8-curve = { path = "curve" }

# For Schnorr signatures
eckfp8-schnorr = { path = "schnorr" }

# For ZK circuits
eckfp8-circuit = { path = "circuit" }
```

**Minimum Supported Rust Version (MSRV)**: 1.70+

## Quick Start

### Signing and Verifying

```rust
use eckfp8_schnorr::{SigningKey, VerifyingKey};
use p3_baby_bear::BabyBear;
use rand::thread_rng;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = thread_rng();

    // Generate keypair
    let signing_key = SigningKey::random(&mut rng);
    let verifying_key = signing_key.verifying_key();

    // Sign a message
    let message = vec![BabyBear::from_u32(42), BabyBear::from_u32(1337)];
    let signature = signing_key.sign(&mut rng, &message)?;

    // Verify signature
    let is_valid = verifying_key.verify(&message, &signature)?;
    assert!(is_valid);

    println!("Signature verified successfully!");
    Ok(())
}
```

### Zero-Knowledge Proof Generation

```rust
use eckfp8_circuit::{SignatureWitness, SchnorrAir, build_schnorr_trace};
use p3_uni_stark::{prove, verify};
// ... (see circuit/examples/prove_verify_signature.rs for complete example)

let witness = SignatureWitness::new(&signature, &verifying_key, &message);
let air = SchnorrAir::new(256); // trace height (power of 2)
let trace = build_schnorr_trace(&witness, 256);

// Generate STARK proof
let proof = prove(&air, &trace, &config);

// Verify proof
let is_valid = verify(&air, &proof, &config);
```

## Architecture

### Cryptographic Foundations

The library is built on three foundational layers:

```
┌─────────────────────────────────────┐
│   ZK-SNARK Circuit (Plonky3 AIR)    │  ← circuit crate
├─────────────────────────────────────┤
│   Schnorr Signature Scheme          │  ← schnorr crate
├─────────────────────────────────────┤
│   KoalaBear Elliptic Curve (Fp8)    │  ← curve crate
└─────────────────────────────────────┘
```

### Curve Specifications

**KoalaBear Elliptic Curve over Fp8**

- **Equation**: `y² = x³ + 3u·x + 42639` where `u` is the primitive element of Fp8
- **Base Field**: KoalaBear prime field `p = 2130706433` (31-bit prime)
- **Extension**: Degree-8 binomial extension field
- **Scalar Field**: 252-bit prime order group
  ```
  Order: 0xf06e44682c2aa440a8ba3e3a29ee4ebaa1ea2c8e76be5cdf1f
  ```
- **Cofactor**: 1 (prime-order curve)
- **Security Level**: ~124 bits (Pollard-Rho complexity)

### Field Arithmetic

- **Montgomery Multiplication**: Scalar field operations use Montgomery form for 3-5× speedup
- **Native Operations**: Fp8 arithmetic fully optimized using KoalaBear field
- **Non-Native Scalars**: Circuit representation uses 9×28-bit limbs for efficient constraints

## Crates

### `curve` Crate

Elliptic curve implementation with optimized group operations.

**Key Types**:
- `Affine` - Affine coordinates (x, y)
- `Projective` - Projective coordinates (X:Y:Z) for efficient computation
- `BaseField` - Fp8 extension field wrapper
- `ScalarField` - 252-bit scalar field with Montgomery arithmetic

**Features**:
- Precomputed generator tables for fixed-base multiplication
- Windowed multiplication (4-bit windows)
- Multi-scalar multiplication (MSM)
- Full serialization support via `serde`

**Example**:
```rust
use eckfp8_curve::{Affine, ScalarField, Group};
use rand::thread_rng;

let mut rng = thread_rng();
let scalar = ScalarField::random(&mut rng);
let point = Affine::GENERATOR.scalar_mul(&scalar);
```

See [curve/README.md](curve/README.md) for detailed API documentation.

### `schnorr` Crate

Schnorr signature implementation with Poseidon2 hashing.

**Key Types**:
- `SigningKey` - Secret key (32 bytes)
- `VerifyingKey` - Public key (40 bytes: 2×Fp8 coordinates)
- `Signature` - Signature (72 bytes: R point + s scalar)

**Algorithm**:

*Signing*:
1. Generate random nonce `k ← ScalarField`
2. Compute commitment `R = G × k`
3. Compute challenge `e = Poseidon2(R || pk || msg)`
4. Compute response `s = k + e × sk`
5. Return `(R, s)`

*Verification*:
1. Recompute challenge `e = Poseidon2(R || pk || msg)`
2. Check equation `G × s = R + pk × e`

**Security Features**:
- Poseidon2 Fiat-Shamir for non-malleability
- Deterministic nonce generation (RFC 6979 style available)
- Comprehensive rejection testing

**Example**:
```rust
use eckfp8_schnorr::{SigningKey, SK_SIZE, PK_SIZE, SIG_SIZE};

// Serialization sizes
assert_eq!(SK_SIZE, 32);  // Secret key
assert_eq!(PK_SIZE, 40);  // Public key
assert_eq!(SIG_SIZE, 72); // Signature

// Serialize/deserialize
let sk_bytes = bincode::serialize(&signing_key)?;
let sk: SigningKey = bincode::deserialize(&sk_bytes)?;
```

See [schnorr/README.md](schnorr/README.md) for detailed API documentation.

### `circuit` Crate

ZK-SNARK circuits for signature verification using Plonky3 framework.

**Key Types**:
- `SignatureWitness` - Witness data (signature + public inputs + challenge)
- `CircuitPoint` - Fp8 point representation (16 KoalaBear limbs)
- `CircuitScalar` - Scalar representation (9×28-bit limbs)
- `SchnorrAir` - Main STARK constraint system

**Architecture**:

The circuit verifies `G × s - pk × e = R` using three sub-components:

1. **Generator Multiplication**: `G × s` (preprocessed trace)
2. **Variable-Base Multiplication**: `pk × e`
3. **Point Addition/Subtraction**: Combining results

**Example**:
```rust
use eckfp8_circuit::{SignatureWitness, SchnorrAir, build_schnorr_trace};

let witness = SignatureWitness::new(&signature, &verifying_key, &message);
let trace = build_schnorr_trace(&witness, 256);
let air = SchnorrAir::new(256);

// Integrate with Plonky3 STARK prover
```

See [circuit/README.md](circuit/README.md) for detailed API documentation.

## Examples

All examples can be run from the workspace root:

```bash
# Basic Schnorr signature demonstration
cargo run --example schnorr -p schnorr

# Full ZK proof generation and verification
cargo run --example prove_verify_signature -p circuit

# Run with release optimizations for accurate benchmarks
cargo run --release --example prove_verify_signature -p circuit
```

**Available Examples**:
- [schnorr/examples/schnorr.rs](schnorr/examples/schnorr.rs) - Sign, serialize, verify workflow
- [circuit/examples/prove_verify_signature.rs](circuit/examples/prove_verify_signature.rs) - Complete STARK proof pipeline

## Performance

### Benchmarks

Run benchmarks with:

```bash
# Curve operations
cargo bench -p curve

# Signature operations
cargo bench -p schnorr

# Circuit generation
cargo bench -p circuit
```

### Project Structure

```
eckfp8/
├── curve/              # Elliptic curve implementation
│   ├── src/
│   │   ├── affine.rs          # Affine coordinates
│   │   ├── projective.rs      # Projective coordinates
│   │   ├── scalarfield.rs     # Montgomery arithmetic
│   │   ├── basefield.rs       # Fp8 extension field
│   │   ├── group.rs           # Group trait & scalar mul
│   │   └── generator_table.rs # Precomputed tables
│   ├── benches/
│   └── tests/
├── schnorr/            # Schnorr signatures
│   ├── src/
│   │   ├── keys.rs            # Key generation & signing
│   │   ├── signatures.rs      # Signature & hashing
│   │   └── errors.rs          # Error types
│   ├── examples/
│   ├── benches/
│   └── tests/
├── circuit/            # ZK-SNARK circuits
│   ├── src/
│   │   ├── schnorr_air.rs     # Main circuit
│   │   ├── scalar_mul_air.rs  # Scalar multiplication
│   │   ├── poseidon2_hash_air.rs # Hash verification
│   │   └── signature_witness.rs # Witness generation
│   └── examples/
└── README.md
```

### Testing

```bash
# Unit tests
cargo test --workspace

# Integration tests
cargo test --workspace --test '*'

# Doc tests
cargo test --workspace --doc

# Test with all features
cargo test --workspace --all-features
```

## License

This project is dual-licensed under:

- **MIT License** ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)
- **Apache License 2.0** ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)

You may choose either license for your use.

## Acknowledgments

- Built with [Plonky3](https://github.com/Plonky3/Plonky3) STARK framework
- KoalaBear fields from Plonky3 ecosystem
- Poseidon2 hash function implementation
