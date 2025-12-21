use core::ops::{Add, AddAssign, Neg, Sub, SubAssign};

/// A scalar type that can expose its canonical 256-bit representation.
pub trait ScalarBits {
    fn to_u64_limbs(&self) -> [u64; 4];
}

/// Basic additive group behavior for curve points.
///
/// This trait centralizes scalar multiplication and related utilities so
/// point types can share one correct implementation.
pub trait Group:
    Sized + Copy + Add<Output = Self> + AddAssign + Sub<Output = Self> + SubAssign + Neg<Output = Self>
{
    type Scalar: ScalarBits;

    /// Return the identity element.
    fn identity() -> Self;
    /// Return true if this element is the identity.
    fn is_identity(&self) -> bool;
    /// A fixed generator for the curve group.
    fn generator() -> Self;
    /// Return 2 * self.
    fn double(&self) -> Self;
    /// Return -self.
    fn negate(&self) -> Self;

    /// Double-and-add scalar multiplication.
    #[inline]
    fn scalar_mul(&self, scalar: &Self::Scalar) -> Self {
        let scalar_limbs = scalar.to_u64_limbs();
        let mut result = Self::identity();
        let mut temp = *self;

        for &limb in scalar_limbs.iter() {
            let mut bits = limb;
            for _ in 0..64 {
                if bits & 1 == 1 {
                    result = result + temp;
                }
                temp = temp.double();
                bits >>= 1;
            }
        }

        result
    }

    /// Fixed-window (width = 4) scalar multiplication.
    fn scalar_mul_windowed(&self, scalar: &Self::Scalar) -> Self {
        if self.is_identity() {
            return Self::identity();
        }

        let mut table = [Self::identity(); 16];
        table[1] = *self;

        for i in 2..16 {
            table[i] = if i % 2 == 0 {
                table[i / 2].double()
            } else {
                table[i - 1] + table[1]
            };
        }

        let scalar_limbs = scalar.to_u64_limbs();
        let mut result = Self::identity();

        for &limb in scalar_limbs.iter().rev() {
            for shift in (0..64).step_by(4).rev() {
                result = result.double();
                result = result.double();
                result = result.double();
                result = result.double();

                let window = ((limb >> shift) & 0xF) as usize;
                if window != 0 {
                    result = result + table[window];
                }
            }
        }

        result
    }

    /// Multiply by a small `u64` scalar.
    fn mul_u64(&self, n: u64) -> Self {
        if n == 0 {
            return Self::identity();
        }
        if n == 1 {
            return *self;
        }

        let mut result = Self::identity();
        let mut temp = *self;
        let mut bits = n;

        while bits > 0 {
            if bits & 1 == 1 {
                result = result + temp;
            }
            temp = temp.double();
            bits >>= 1;
        }

        result
    }

    /// Naive multi-scalar multiplication.
    fn multi_scalar_mul(points: &[Self], scalars: &[Self::Scalar]) -> Self {
        assert_eq!(
            points.len(),
            scalars.len(),
            "Points and scalars must have same length"
        );

        let mut result = Self::identity();
        for (point, scalar) in points.iter().zip(scalars.iter()) {
            result = result + point.scalar_mul(scalar);
        }
        result
    }
}
