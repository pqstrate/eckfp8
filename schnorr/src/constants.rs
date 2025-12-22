//! Constants used in the Schnorr signature scheme implementation.

/// Size of a serialized public verifying key in bytes.
///
/// A verifying key is a point on the KoalaBear curve (Fp8 extension),
/// requiring 40 bytes when serialized.
pub const PK_SIZE: usize = 40;

/// Size of a serialized secret signing key in bytes.
///
/// A signing key is a scalar in the scalar field, requiring 32 bytes
/// when serialized.
pub const SK_SIZE: usize = 32;

/// Size of a serialized signature in bytes.
///
/// A signature consists of:
/// - A point R (40 bytes)
/// - A scalar s (32 bytes)
/// Total: 72 bytes
pub const SIG_SIZE: usize = 72;

/// Width parameter for the Poseidon2 permutation.
///
/// This is the total state size of the Poseidon2 sponge construction
/// used for hashing in the Fiat-Shamir transform.
pub(crate) const POSEIDON2_WIDTH: usize = 16;

/// Rate parameter for the Poseidon2 sponge.
///
/// This is the number of field elements absorbed per permutation call.
pub(crate) const POSEIDON2_RATE: usize = 8;

/// Output size for the Poseidon2 hash.
///
/// This is the number of field elements in the hash digest.
pub(crate) const POSEIDON2_OUT: usize = 8;
