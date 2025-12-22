# eckfp8

Rust workspace for the KoalaBear Fp8 elliptic curve, Schnorr signatures, and a ZK-friendly verification circuit.

## Workspace crates

- `curve`: KoalaBear elliptic curve over an Fp8 extension field, plus scalar field arithmetic.
- `schnorr`: Schnorr signature scheme over the curve with Poseidon2 Fiat-Shamir.
- `circuit`: Signature verification circuit using native field ops and non-native scalar arithmetic.

## Quick start

```bash
cargo build
cargo test
```

## Examples

```bash
# Run the circuit verification example
cargo run --example verify_signature -p circuit
```

## Layout

```
.
├── circuit
├── curve
└── schnorr
```

## License

See `LICENSE`.
