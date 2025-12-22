//! Error types for the Schnorr signature scheme.

/// Errors that can occur during signing and verification operations.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum SchnorrError {
    /// A point at infinity was encountered where a valid curve point was expected.
    ///
    /// This error occurs when:
    /// - The verifying key is the point at infinity
    /// - The signature commitment point R is the point at infinity
    /// - A nonce generates the point at infinity (extremely unlikely)
    ///
    /// In practice, this error should be extremely rare for randomly generated keys
    /// and nonces, as the probability of generating the point at infinity is negligible.
    InvalidPoint,
}
