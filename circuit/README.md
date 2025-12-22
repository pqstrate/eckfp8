# Circuit: Signature Verification Circuit for KoalaBear Schnorr Signatures

A ZK-SNARK-friendly circuit implementation for verifying Schnorr signatures over the KoalaBear elliptic curve. This circuit uses native field operations (KoalaBear field) for all elliptic curve computations, with special handling for scalar field arithmetic using non-native techniques.

## Overview

This crate provides a circuit-based implementation of the Schnorr signature verification equation:

```
G * s == R + pk * e
```

where:
- `G` is the generator point
- `s` is the signature response scalar
- `R` is the signature commitment point
- `pk` is the public key
- `e` is the Fiat-Shamir challenge: `e = H(R || pk || msg)`

## Architecture

### Native Field Operations

All elliptic curve point operations (addition, doubling) are performed using native KoalaBear field arithmetic. Points are represented in affine coordinates where each coordinate is a degree-8 extension field element (BaseField = BinomialExtensionField<KoalaBear, 8>).

### Non-Native Scalar Arithmetic

Since the scalar field is different from the base field, scalar field operations require non-native arithmetic. We represent scalar field elements as multiple KoalaBear limbs (9 limbs of 28 bits each) and implement limb-based arithmetic.

## Modules

- **`scalar_arithmetic.rs`**: Scalar representation in KoalaBear limbs
- **`point_ops.rs`**: Circuit-friendly elliptic curve point operations
- **`signature_witness.rs`**: Signature witness construction (inputs for AIR proofs)

## Usage

```rust
use circuit::{build_schnorr_trace, SchnorrAir, SignatureWitness};
use schnorr::{SigningKey, VerifyingKey};
use p3_baby_bear::BabyBear;
use p3_field::PrimeCharacteristicRing;
use rand::rng;

// Generate a key pair
let mut rng = rng();
let signing_key = SigningKey::random(&mut rng);
let verifying_key = signing_key.verifying_key();

// Create and sign a message
let message = vec![BabyBear::from_u32(1), BabyBear::from_u32(2)];
let signature = signing_key.sign(&mut rng, &message).unwrap();

// Create witness and verify in circuit
let witness = SignatureWitness::new(&signature, &verifying_key, &message).unwrap();
let trace = build_schnorr_trace(&witness);
let air = SchnorrAir::new(trace.trace.height());

assert!(is_valid);
```

## Example

Run the example to see the circuit in action:

```bash
cargo run --example verify_signature
```

This demonstrates:
1. Key generation
2. Message signing
3. Standard verification
4. Circuit-based verification
5. Circuit complexity analysis
6. Invalid signature detection

## Circuit Complexity

The circuit implements signature verification with the following approximate complexity:

- **Scalar multiplications**: 2 (G*s and pk*e)
- **Point additions**: 1 (R + pk*e)
- **Point doublings**: ~512 (256 bits × 2 scalar muls)
- **Field multiplications**: ~5,120 (estimated)

## Optimization Opportunities

Several optimizations can be applied to reduce circuit size:

1. **Fixed-base scalar multiplication**: Precomputed tables for G*s
2. **Windowed methods**: Reduce doublings in variable-base scalar multiplication
3. **Shamir's trick**: Optimize double scalar multiplication (G*s - pk*e)
4. **Batch verification**: Verify multiple signatures more efficiently together
5. **Lookup tables**: For common sub-operations

## Design Choices

### Why Affine Coordinates?

The circuit uses affine coordinates for points because:
- Simpler formulas with fewer field operations
- No need for division-free projective formulas
- Division is relatively cheap in SNARK circuits

### Why Limb-based Scalars?

Scalar field elements are represented as limbs because:
- The scalar field is different from the circuit field (non-native)
- Limb-based arithmetic allows range checks and proper overflow handling
- Compatible with SNARK constraint systems

### Circuit vs. Arithmetic

This implementation prioritizes clarity and correctness over optimization. For production use, consider:
- Using optimized lookup tables for scalar multiplication
- Implementing Plonky3-specific optimizations
- Batching multiple verifications

## Integration with ZK Systems

This circuit is designed to be compatible with:
- **Plonky3**: Native support for KoalaBear and BabyBear fields
- **STARKs**: Efficient with the native field arithmetic
- **SNARKs**: Can be adapted for R1CS or other constraint systems

To integrate with a proof system:
1. Convert circuit operations to constraints
2. Add range checks for limb-based arithmetic
3. Implement proper witness generation
4. Add public input handling for messages/signatures

## Testing

Run the test suite:

```bash
cargo test
```

Tests cover:
- Scalar field arithmetic (conversion, addition, multiplication)
- Point operations (conversion, addition, doubling, scalar multiplication)
- Signature verification (valid and invalid signatures)
- Circuit complexity analysis

## Status

✅ **All tests passing!** The circuit correctly verifies Schnorr signatures with all test cases passing.

## Future Work

- [ ] Implement windowed scalar multiplication
- [ ] Add precomputed tables for fixed-base mul
- [ ] Optimize for specific proof systems (Plonky3, etc.)
- [ ] Add batch verification support
- [ ] Implement constraint generation for R1CS
- [ ] Add benchmarks comparing circuit vs. native verification

## License

This project is part of the eckfp8 workspace.
