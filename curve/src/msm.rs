use crate::generator_table::affine_table;
use crate::group::ScalarBits;
use crate::{Affine, ScalarField};

/// Compute a * G + b * P using precomputed generator table and a point table.
pub fn double_scalar_mul_basepoint_affine(
    a: &ScalarField,
    b: &ScalarField,
    point: &Affine,
) -> Affine {
    let base_table = affine_table();
    let mut point_table = [Affine::INFINITY; 256];
    point_table[1] = *point;
    for i in 2..256 {
        point_table[i] = point_table[i - 1] + point_table[1];
    }

    let a_limbs = a.to_u64_limbs();
    let b_limbs = b.to_u64_limbs();
    let mut result = Affine::INFINITY;

    for limb_idx in (0..4).rev() {
        let a_limb = a_limbs[limb_idx];
        let b_limb = b_limbs[limb_idx];
        for shift in (0..64).step_by(8).rev() {
            for _ in 0..8 {
                result = result.double();
            }

            let a_window = ((a_limb >> shift) & 0xFF) as usize;
            if a_window != 0 {
                result += base_table[a_window];
            }

            let b_window = ((b_limb >> shift) & 0xFF) as usize;
            if b_window != 0 {
                result += point_table[b_window];
            }
        }
    }

    result
}
