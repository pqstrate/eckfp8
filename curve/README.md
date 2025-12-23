# KoalaBear Elliptic Curve (Fp8)

A Rust implementation of an elliptic curve defined over the KoalaBear Fp8 extension field, optimized for zero-knowledge proof systems and cryptographic applications.

## Table of Contents

- [Overview](#overview)
- [Curve Parameters](#curve-parameters)
- [Curve Selection Rationale](#curve-selection-rationale)
- [Security Properties](#security-properties)
- [Field Structure](#field-structure)
- [Features](#features)
- [Usage](#usage)
- [Performance](#performance)
- [References](#references)

## Overview

This library provides an elliptic curve over the degree-8 extension field of the KoalaBear base field. The curve is designed for efficient cryptographic operations and zero-knowledge proof systems, particularly those using the Plonky3 framework.

### Curve Equation

```
E(Fp8): y² = x³ + a·x + b
```

where:
- **a = 3u** (where u is the primitive element of the extension field)
- **b = 42639**

The base field is GF(p) where **p = 2130706433 = 2³¹ - 2²⁴ + 1** (the KoalaBear prime).

## Curve Parameters

### Base Field
- **Prime modulus**: p = 2,130,706,433 (2³¹ - 2²⁴ + 1)
- **Extension degree**: 8
- **Field size**: p⁸ ≈ 2²⁴⁸

### Curve Coefficients
- **a**: 3u (coefficient of x term)
  - In extension field representation: [0, 3, 0, 0, 0, 0, 0, 0]
- **b**: 42639 (constant term)
  - In extension field representation: [42639, 0, 0, 0, 0, 0, 0, 0]

### Group Parameters
- **Curve order**: 424,804,331,891,979,973,455,971,894,938,199,991,855,800,421,968,298,112,348,210,302,325,367,590,273
- **Order (hex)**: 0xf06e44682c2aa440f5f26a5ae1748ff85ccc2efc3068faf2154ff8a2e94d81
- **Order (bits)**: 248 bits
- **Cofactor**: 1 (prime order curve)
- **2-adicity**: 7 (v₂(q-1) = 7)

### Generator Points

The library provides two generator points derived using Simplified SWU (SSWU) hash-to-curve:

#### Primary Generator (from 'ZKM2')
```
x = 1195559694·u⁷ + 1368232771·u⁶ + 438909494·u⁵ + 1825476283·u⁴
    + 1299273209·u³ + 2115217807·u² + 1763905369·u + 1813646457

y = 2077084094·u⁷ + 434578416·u⁶ + 125328769·u⁵ + 1286889583·u⁴
    + 655051022·u³ + 1365273355·u² + 840779000·u + 376996212
```

#### Pedersen Generator (from 'ZKM2 - Pedersen')
```
x = 1741425845·u⁷ + 1750752810·u⁶ + 7156361·u⁵ + 1949220725·u⁴
    + 543192455·u³ + 358531926·u² + 550988532·u + 1709677626

y = 71894712·u⁷ + 1876016551·u⁶ + 1684498755·u⁵ + 598910111·u⁴
    + 156828552·u³ + 978667041·u² + 1399061592·u + 548133034
```

## Curve Selection Rationale

### KoalaBear field

1. **STARK-Friendly Construction**: Extension fields work well with STARK proof systems and algebraic hash functions
2. **Integration with Plonky3**: Seamless compatibility with the Plonky3 ZK-STARK framework which uses KoalaBear fields

### Efficient Coefficients

- **a = 3u**:
  - Simplifies point addition formulas

- **b = 42639**:
  - Chosen to ensure the curve has prime order

### Hash-to-Curve Method

Generators are derived using **Simplified SWU (SSWU)** mapping

## Security Properties

### Discrete Logarithm Security

- **Pollard-Rho Security**: ~123.78 bits
  - Based on group order of ~2²⁴⁸
  - Exceeds 128-bit security target

- **Twist Security**: ~120.86 bits (Pollard-Rho)
  - Protects against invalid point attacks
  - Ensures robust security even with twist points

### Embedding Degree

```
k > 2²⁴²
```

### Prime Order

- **Cofactor = 1**: Every non-identity point generates the full group

### 2-Adicity

- **v₂(q-1) = 7**: Enables 128-point FFTs

## Field Structure

### Base Field (KoalaBear)

```
GF(2130706433) where p = 2^31 - 2^24 + 1
```

### Extension Field (Fp8)

The degree-8 extension is constructed as:

```
Fp8 = Fp[u] / (irreducible polynomial)
```

Elements are represented as:
```
a₀ + a₁u + a₂u² + a₃u³ + a₄u⁴ + a₅u⁵ + a₆u⁶ + a₇u⁷
```

where each aᵢ ∈ Fp (KoalaBear base field).

### Scalar Field

The scalar field for the curve is a 248-bit prime field:

```
Fr: 0xf06e44682c2aa440f5f26a5ae1748ff85ccc2efc3068faf2154ff8a2e94d81
```

Scalars are represented in Montgomery form using four 64-bit limbs for efficient modular arithmetic.

## Usage

### Basic Example

```rust
use curve::{Affine, Group, ScalarField, RandomField};
use rand::thread_rng;

fn main() {
    let mut rng = thread_rng();

    // Generate random scalar
    let scalar = ScalarField::random(&mut rng);

    // Scalar multiplication with generator
    let point = <Affine as Group>::mul_generator(&scalar);

    // Point addition
    let point2 = <Affine as Group>::mul_generator(&ScalarField::ONE);
    let sum = point + point2;

    // Verify point is on curve
    assert!(sum.is_on_curve());
}
```

### Scalar Multiplication

```rust
use curve::{Affine, Group, ScalarField};

// Fixed-base multiplication (using precomputed table)
let result = <Affine as Group>::mul_generator(&scalar);

// Variable-base multiplication
let result = point.scalar_mul(&scalar);

// Multi-scalar multiplication
let scalars = vec![scalar1, scalar2, scalar3];
let points = vec![point1, point2, point3];
let result = <Affine as Group>::multi_scalar_mul(&scalars, &points);
```

### Coordinate Conversions

```rust
use curve::{Affine, Projective};

// Affine to Projective
let proj = Projective::from_affine(&affine_point);

// Projective to Affine
let affine = proj.to_affine();

// Batch normalization (more efficient)
let projective_points = vec![p1, p2, p3];
let affine_points = Projective::batch_normalize(&projective_points);
```

### Point Arithmetic

```rust
use curve::Affine;

// Point addition
let sum = point1 + point2;

// Point doubling
let doubled = point.double();

// Point negation
let negated = point.negate();

// Point subtraction
let diff = point1 - point2;
```

## Performance

### Benchmarks

Run benchmarks with:

```bash
cargo bench
```

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
curve = { path = "path/to/curve" }
p3-koala-bear = { git = "https://github.com/Plonky3/Plonky3.git" }
rand = "0.9"
```

## Building

```bash
# Build the library
cargo build --release

# Run tests
cargo test

# Run benchmarks
cargo bench

# Generate documentation
cargo doc --open
```

## API Documentation

Generate and view the full API documentation:

```bash
cargo doc --open
```

## Dependencies

- **p3-koala-bear**: KoalaBear field implementation from Plonky3
- **p3-field**: Field trait definitions
- **rand**: Random number generation
- **serde**: Serialization support
- **num-bigint**: Big integer arithmetic for field operations

## References

1. **Plonky3**: https://github.com/Plonky3/Plonky3
2. **KoalaBear Field**: "Mersenne31 Field" in STARK literature
3. **Hash-to-Curve**: [RFC 9380](https://datatracker.ietf.org/doc/rfc9380/)
4. **ECC Handbook**: "Guide to Elliptic Curve Cryptography" by Hankerson, Menezes, and Vanstone
5. **Montgomery Arithmetic**: "Modular Multiplication Without Trial Division" by Montgomery (1985)

## License

This project is part of the eckfp8 library suite.
