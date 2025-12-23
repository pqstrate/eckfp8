//! # KoalaBear Elliptic Curve (Fp8)
//!
//! Production-grade elliptic curve implementation over the KoalaBear degree-8 extension field.
//!
//! ## Overview
//!
//! This crate provides a complete elliptic curve implementation featuring:
//! - **Affine and Projective Coordinates**: Optimized point representations
//! - **Montgomery Scalar Field**: 252-bit scalar arithmetic with ~3-5× speedup
//! - **Precomputed Tables**: Fixed-base multiplication optimization
//! - **Serialization**: Full `serde` support for all types
//!
//! ## Curve Specification
//!
//! **KoalaBear Elliptic Curve over Fp8**
//!
//! - **Equation**: `y² = x³ + 3u·x + 42639` where `u` is the primitive element of Fp8
//! - **Base Field**: KoalaBear prime field `p = 2130706433` (31-bit prime)
//! - **Extension**: Degree-8 binomial extension field
//! - **Scalar Field**: 252-bit prime order
//!   ```text
//!   Order: 0xf06e44682c2aa440a8ba3e3a29ee4ebaa1ea2c8e76be5cdf1f
//!   ```
//! - **Cofactor**: 1 (prime-order curve)
//! - **Security**: ~124 bits (Pollard-Rho complexity)
//!
//! ## Quick Start
//!
//! ```rust
//! use curve::{Affine, ScalarField, Group, RandomField};
//! use rand::rng;
//!
//! let mut rng = rng();
//!
//! // Generate random scalar
//! let scalar = ScalarField::random(&mut rng);
//!
//! // Scalar multiplication: P = G × scalar
//! let point = Affine::generator().scalar_mul(&scalar);
//!
//! // Point addition
//! let point2 = point + Affine::generator();
//!
//! // Verify point is on curve
//! assert!(point.is_on_curve());
//! assert!(point2.is_on_curve());
//! ```
//!
//! ## Key Types
//!
//! - [`Affine`]: Affine coordinate representation (x, y)
//! - [`Projective`]: Projective coordinate representation (X:Y:Z) for efficient computation
//! - [`ScalarField`]: 252-bit scalar field with Montgomery arithmetic
//! - [`BaseField`]: Fp8 extension field (8 KoalaBear elements)
//! - [`Group`]: Trait providing scalar multiplication algorithms
//! - [`ScalarBits`]: Trait exposing canonical bit representation
//!
//! ## Performance Optimizations
//!
//! ### Fixed-Base Multiplication
//!
//! For scalar multiplication by the generator, use precomputed tables:
//!
//! ```rust
//! use curve::{mul_generator_affine, ScalarField, RandomField};
//! use rand::rng;
//!
//! let mut rng = rng();
//! let scalar = ScalarField::random(&mut rng);
//!
//! // Fast: uses precomputed 4-bit windowed table
//! let point = mul_generator_affine(&scalar);
//! ```
//!
//! ### Double Scalar Multiplication
//!
//! For computing `a·G + b·P` (common in signature verification):
//!
//! ```rust
//! use curve::{double_scalar_mul_basepoint_affine, ScalarField, Affine, RandomField, Group};
//! use rand::rng;
//!
//! let mut rng = rng();
//! let a = ScalarField::random(&mut rng);
//! let b = ScalarField::random(&mut rng);
//! let public_key = Affine::generator().scalar_mul(&b);
//!
//! // Optimized: computes a·G + b·P efficiently
//! let result = double_scalar_mul_basepoint_affine(&a, &b, &public_key);
//! assert!(result.is_on_curve());
//! ```
//!
//! ### Windowed Multiplication
//!
//! For variable-base scalar multiplication with improved performance:
//!
//! ```rust
//! use curve::{Affine, ScalarField, Group, RandomField};
//! use rand::rng;
//!
//! let mut rng = rng();
//! let scalar = ScalarField::random(&mut rng);
//! let base = Affine::generator();
//!
//! // Uses 4-bit windows to reduce point doublings
//! let result = base.scalar_mul_windowed(&scalar);
//! ```
//!
//! ## Field Arithmetic
//!
//! ### Scalar Field (Montgomery Form)
//!
//! The scalar field uses Montgomery multiplication for efficiency:
//!
//! ```rust
//! use curve::{ScalarField, RandomField};
//! use rand::rng;
//! use p3_field::Field;
//!
//! let mut rng = rng();
//! let a = ScalarField::random(&mut rng);
//! let b = ScalarField::random(&mut rng);
//!
//! let sum = a + b;
//! let product = a * b;
//! let inverse = a.inverse(); // Extended Euclidean algorithm
//! ```
//!
//! ### Base Field (Fp8 Extension)
//!
//! Fp8 is represented as 8 KoalaBear field elements:
//!
//! ```rust,ignore
//! use curve::{BaseField, KoalaBear};
//! use p3_field::AbstractField;
//!
//! let fp8_element = BaseField([
//!     KoalaBear::one(),
//!     KoalaBear::two(),
//!     // ... 6 more elements
//!     KoalaBear::zero(),
//! ]);
//! ```
//!
//! ## Examples
//!
//! See the `examples/` directory for complete usage examples.
//!
//! ## References
//!
//! - Plonky3 framework: <https://github.com/Plonky3/Plonky3>
//! - KoalaBear field specification
//! - Montgomery arithmetic: Peter Montgomery (1985)

#[deny(missing_docs)]
mod affine;
mod basefield;
mod generator_table;
mod group;
mod msm;
mod projective;
mod random;
mod scalarfield;

pub use affine::Affine;
pub use basefield::BaseField;
pub use basefield::{to_bytes, to_u32s};
pub use generator_table::mul_generator_affine;
pub use group::{Group, ScalarBits};
pub use msm::double_scalar_mul_basepoint_affine;
pub use p3_koala_bear::KoalaBear;
pub use projective::Projective;
pub use random::RandomField;
pub use scalarfield::ScalarField;
