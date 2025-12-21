use core::ops::{Add, AddAssign, Neg, Sub, SubAssign};

pub trait ScalarBits {
    fn to_u64_limbs(&self) -> [u64; 4];
}

pub trait Group:
    Sized + Copy + Add<Output = Self> + AddAssign + Sub<Output = Self> + SubAssign + Neg<Output = Self>
{
    type Scalar: ScalarBits;

    fn identity() -> Self;
    fn is_identity(&self) -> bool;
    fn generator() -> Self;
    fn double(&self) -> Self;
    fn negate(&self) -> Self;

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
