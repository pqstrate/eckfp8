# Schnorr Signatures over KoalaBear Curve

A Rust implementation of the Schnorr signature scheme using the KoalaBear elliptic curve (defined over an Fp8 extension field) and Poseidon2 hash function over the KoalaBear field.

## Features

- **Schnorr Signatures**: Classic Schnorr signature scheme with provable security
- **KoalaBear Curve**: Elliptic curve over an Fp8 extension field optimized for zero-knowledge proofs
- **Poseidon2 Hashing**: Uses Poseidon2 over the KoalaBear field for the Fiat-Shamir transform
- **Serialization**: Built-in support for `serde` serialization/deserialization
- **Type Safety**: Separate types for signing and verifying keys

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
schnorr = { path = "path/to/schnorr" }
p3-baby-bear = { git = "https://github.com/Plonky3/Plonky3.git" }
p3-field = { git = "https://github.com/Plonky3/Plonky3.git" }
rand = "0.9"
```

## Usage

### Basic Example

```rust
use schnorr::{SigningKey, VerifyingKey, Signature};
use p3_baby_bear::BabyBear;
use p3_field::PrimeCharacteristicRing;
use rand::thread_rng;

fn main() {
    // Generate a random signing key
    let mut rng = thread_rng();
    let signing_key = SigningKey::random(&mut rng);

    // Derive the corresponding verifying key
    let verifying_key = signing_key.verifying_key();

    // Create a message (using KoalaBear field elements)
    let message = [
        BabyBear::from_u32(1),
        BabyBear::from_u32(2),
        BabyBear::from_u32(3),
    ];

    // Sign the message
    let signature = signing_key.sign(&mut rng, &message).expect("signing failed");

    // Verify the signature
    let is_valid = verifying_key.verify(&message, &signature).expect("verification failed");
    assert!(is_valid);
}
```

### Key Generation

```rust
use schnorr::SigningKey;
use rand::thread_rng;

let mut rng = thread_rng();
let signing_key = SigningKey::random(&mut rng);
let verifying_key = signing_key.verifying_key();
```

### Signing

```rust
use p3_baby_bear::BabyBear;
use p3_field::PrimeCharacteristicRing;

let message = [BabyBear::from_u32(42)];
let signature = signing_key.sign(&mut rng, &message)?;
```

### Verification

```rust
let is_valid = verifying_key.verify(&message, &signature)?;
if is_valid {
    println!("Signature is valid!");
} else {
    println!("Signature is invalid!");
}
```

## How It Works

### Schnorr Signature Scheme

The implementation follows the standard Schnorr signature protocol:

**Key Generation:**
- Secret key: random scalar `sk` in the scalar field
- Public key: point `pk = G * sk` where `G` is the curve generator

**Signing:**
1. Generate random nonce `k`
2. Compute commitment `R = G * k`
3. Compute challenge `e = H(R || pk || msg)` using Poseidon2
4. Compute response `s = k + e * sk`
5. Output signature `(R, s)`

**Verification:**
1. Compute challenge `e = H(R || pk || msg)`
2. Check that `G * s == R + pk * e`

### Hash Function

The implementation uses Poseidon2 over the KoalaBear field for the Fiat-Shamir challenge:
- **Width**: 16 field elements
- **Rate**: 8 field elements
- **Output**: 8 field elements

Points are encoded as 16 KoalaBear field elements (8 for x-coordinate, 8 for y-coordinate) before hashing.

## Security Considerations

- **Random Number Generator**: Always use a cryptographically secure RNG (CSRNG) like `rand::thread_rng()` or `rand::rngs::OsRng`
- **Nonce Uniqueness**: Each signature must use a fresh random nonce. Reusing nonces can leak the secret key
- **Key Protection**: Keep signing keys secret and secure. Only share verifying keys
- **Message Encoding**: Messages must be properly encoded as KoalaBear field elements

## API Documentation

### Types

- **`SigningKey`**: Secret key for creating signatures
- **`VerifyingKey`**: Public key for verifying signatures
- **`Signature`**: A signature consisting of a point `R` and scalar `s`
- **`SchnorrError`**: Error type for signature operations

### Constants

- **`PK_SIZE`**: Size of serialized public key (40 bytes)
- **`SK_SIZE`**: Size of serialized secret key (32 bytes)
- **`SIG_SIZE`**: Size of serialized signature (72 bytes)

## Testing

Run the test suite:

```bash
cargo test
```

## Benchmarking

Run benchmarks:

```bash
cargo bench
```

## License

This project is part of the eckfp8 library suite.

## References

- [Schnorr Signatures](https://en.wikipedia.org/wiki/Schnorr_signature)
- [Poseidon2 Hash Function](https://eprint.iacr.org/2023/323)
- [Plonky3 Library](https://github.com/Plonky3/Plonky3)
