//! Signing and verifying keys for the Schnorr signature scheme.

use curve::{Affine, Group, RandomField, ScalarField};
use p3_baby_bear::BabyBear;
use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::errors::SchnorrError;
use crate::signatures::{Signature, hash_challenge};

/// A secret signing key for creating Schnorr signatures.
///
/// The signing key is a random scalar in the scalar field of the KoalaBear curve.
/// It must be kept secret and protected from unauthorized access.
///
/// # Example
///
/// ```
/// use schnorr::SigningKey;
/// use rand::thread_rng;
///
/// let mut rng = thread_rng();
/// let signing_key = SigningKey::random(&mut rng);
/// ```
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SigningKey {
    scalar: ScalarField,
}

/// A public verifying key for verifying Schnorr signatures.
///
/// The verifying key is a point on the KoalaBear elliptic curve, derived from
/// the signing key by multiplying the curve generator by the secret scalar.
///
/// # Example
///
/// ```
/// use schnorr::SigningKey;
/// use rand::thread_rng;
///
/// let mut rng = thread_rng();
/// let signing_key = SigningKey::random(&mut rng);
/// let verifying_key = signing_key.verifying_key();
/// ```
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct VerifyingKey {
    point: Affine,
}

impl SigningKey {
    /// Generates a random signing key using the provided random number generator.
    ///
    /// # Arguments
    ///
    /// * `rng` - A cryptographically secure random number generator
    ///
    /// # Example
    ///
    /// ```
    /// use schnorr::SigningKey;
    /// use rand::thread_rng;
    ///
    /// let mut rng = thread_rng();
    /// let signing_key = SigningKey::random(&mut rng);
    /// ```
    pub fn random<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self {
            scalar: ScalarField::random(rng),
        }
    }

    /// Derives the public verifying key from this signing key.
    ///
    /// The verifying key is computed as `G * sk` where `G` is the generator
    /// of the KoalaBear curve and `sk` is the secret scalar.
    ///
    /// # Example
    ///
    /// ```
    /// use schnorr::SigningKey;
    /// use rand::thread_rng;
    ///
    /// let mut rng = thread_rng();
    /// let signing_key = SigningKey::random(&mut rng);
    /// let verifying_key = signing_key.verifying_key();
    /// ```
    pub fn verifying_key(&self) -> VerifyingKey {
        VerifyingKey {
            point: <Affine as Group>::mul_generator(&self.scalar),
        }
    }

    /// Signs a message using this signing key.
    ///
    /// The signature is computed using the Schnorr signature algorithm:
    /// 1. Generate a random nonce `k`
    /// 2. Compute `R = G * k`
    /// 3. Compute challenge `e = H(R || pk || msg)` using Poseidon2
    /// 4. Compute `s = k + e * sk`
    /// 5. Return signature `(R, s)`
    ///
    /// # Arguments
    ///
    /// * `rng` - A cryptographically secure random number generator for the nonce
    /// * `msg` - The message to sign, encoded as BabyBear field elements
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing the signature on success, or a `SchnorrError` if:
    /// - The nonce point is at infinity (extremely unlikely)
    /// - The public key point is at infinity (extremely unlikely)
    ///
    /// # Example
    ///
    /// ```
    /// use schnorr::SigningKey;
    /// use p3_baby_bear::BabyBear;
    /// use p3_field::PrimeCharacteristicRing;
    /// use rand::thread_rng;
    ///
    /// let mut rng = thread_rng();
    /// let signing_key = SigningKey::random(&mut rng);
    /// let message = [BabyBear::from_u32(1), BabyBear::from_u32(2)];
    /// let signature = signing_key.sign(&mut rng, &message).expect("signing failed");
    /// ```
    pub fn sign<R: Rng + ?Sized>(
        &self,
        rng: &mut R,
        msg: &[BabyBear],
    ) -> Result<Signature, SchnorrError> {
        let nonce = ScalarField::random(rng);
        let r = <Affine as Group>::mul_generator(&nonce);
        let pk = self.verifying_key();

        let e = hash_challenge(&r, &pk.point, msg)?;
        let s = nonce + e * self.scalar;

        Ok(Signature { r, s })
    }
}

impl VerifyingKey {
    /// Verifies a signature on a message using this verifying key.
    ///
    /// The verification checks whether the signature equation holds:
    /// `G * s == R + pk * e`, where:
    /// - `G` is the curve generator
    /// - `s` is the signature scalar
    /// - `R` is the signature point
    /// - `pk` is this verifying key
    /// - `e = H(R || pk || msg)` is the challenge hash
    ///
    /// # Arguments
    ///
    /// * `msg` - The message that was signed, encoded as BabyBear field elements
    /// * `sig` - The signature to verify
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing:
    /// - `Ok(true)` if the signature is valid
    /// - `Ok(false)` if the signature is invalid
    /// - `Err(SchnorrError::InvalidPoint)` if the verifying key or signature contains an invalid point
    ///
    /// # Example
    ///
    /// ```
    /// use schnorr::SigningKey;
    /// use p3_baby_bear::BabyBear;
    /// use p3_field::PrimeCharacteristicRing;
    /// use rand::thread_rng;
    ///
    /// let mut rng = thread_rng();
    /// let signing_key = SigningKey::random(&mut rng);
    /// let verifying_key = signing_key.verifying_key();
    /// let message = [BabyBear::from_u32(1), BabyBear::from_u32(2)];
    ///
    /// let signature = signing_key.sign(&mut rng, &message).expect("signing failed");
    /// let is_valid = verifying_key.verify(&message, &signature).expect("verification failed");
    /// assert!(is_valid);
    /// ```
    pub fn verify(&self, msg: &[BabyBear], sig: &Signature) -> Result<bool, SchnorrError> {
        if self.point.is_infinity() || sig.r.is_infinity() {
            return Err(SchnorrError::InvalidPoint);
        }

        let e = hash_challenge(&sig.r, &self.point, msg)?;
        let lhs = Affine::double_scalar_mul_basepoint(&sig.s, &-e, &self.point);

        Ok(lhs == sig.r)
    }
}

impl From<&SigningKey> for VerifyingKey {
    /// Converts a reference to a signing key into a verifying key.
    ///
    /// This is equivalent to calling `signing_key.verifying_key()`.
    fn from(sk: &SigningKey) -> Self {
        sk.verifying_key()
    }
}
