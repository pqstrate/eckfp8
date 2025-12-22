//! Poseidon2 hash AIR for proving hash_challenge inputs in BabyBear.
//!
//! This module proves a Poseidon2 permutation trace and constrains the final
//! digest (first OUT elements of the last permutation output) as public values.

use core::borrow::Borrow;

use p3_air::{Air, AirBuilder, AirBuilderWithPublicValues, BaseAir, BaseAirWithPublicValues};
use p3_baby_bear::{
    default_babybear_poseidon2_16, BabyBear, GenericPoseidon2LinearLayersBabyBear,
    BABYBEAR_RC16_EXTERNAL_FINAL, BABYBEAR_RC16_EXTERNAL_INITIAL, BABYBEAR_RC16_INTERNAL,
};
use p3_field::PrimeCharacteristicRing;
use p3_matrix::{dense::RowMajorMatrix, Matrix};
use p3_poseidon2_air::{num_cols, Poseidon2Cols, RoundConstants, VectorizedPoseidon2Air};
use p3_symmetric::Permutation;

pub const POSEIDON2_WIDTH: usize = 16;
pub const POSEIDON2_RATE: usize = 8;
pub const POSEIDON2_OUT: usize = 8;
pub const POSEIDON2_SBOX_DEGREE: u64 = 7;
pub const POSEIDON2_SBOX_REGISTERS: usize = 1;
pub const POSEIDON2_HALF_FULL_ROUNDS: usize = 4;
pub const POSEIDON2_PARTIAL_ROUNDS: usize = 13;
pub const POSEIDON2_INPUT_LEN: usize = 35;
pub const POSEIDON2_NUM_PERMS: usize = (POSEIDON2_INPUT_LEN + POSEIDON2_RATE - 1) / POSEIDON2_RATE;
pub const POSEIDON2_PACKED_LIMBS: usize = 3;

#[derive(Clone, Debug)]
pub struct Poseidon2HashTrace {
    pub trace: RowMajorMatrix<BabyBear>,
    pub digest: [BabyBear; POSEIDON2_OUT],
    pub num_permutations: usize,
}

pub struct Poseidon2HashAir {
    inner: VectorizedPoseidon2Air<
        BabyBear,
        GenericPoseidon2LinearLayersBabyBear,
        POSEIDON2_WIDTH,
        POSEIDON2_SBOX_DEGREE,
        POSEIDON2_SBOX_REGISTERS,
        POSEIDON2_HALF_FULL_ROUNDS,
        POSEIDON2_PARTIAL_ROUNDS,
        POSEIDON2_NUM_PERMS,
    >,
}

impl Poseidon2HashAir {
    pub fn new() -> Self {
        let constants = RoundConstants::new(
            BABYBEAR_RC16_EXTERNAL_INITIAL,
            BABYBEAR_RC16_INTERNAL,
            BABYBEAR_RC16_EXTERNAL_FINAL,
        );
        Self {
            inner: VectorizedPoseidon2Air::new(constants),
        }
    }
}

impl Default for Poseidon2HashAir {
    fn default() -> Self {
        Self::new()
    }
}

impl BaseAir<BabyBear> for Poseidon2HashAir {
    fn width(&self) -> usize {
        self.inner.width()
    }
}

impl BaseAirWithPublicValues<BabyBear> for Poseidon2HashAir {
    fn num_public_values(&self) -> usize {
        POSEIDON2_INPUT_LEN + POSEIDON2_OUT + POSEIDON2_PACKED_LIMBS
    }
}

impl<AB> Air<AB> for Poseidon2HashAir
where
    AB: AirBuilder<F = BabyBear> + AirBuilderWithPublicValues,
{
    fn eval(&self, builder: &mut AB) {
        self.inner.eval(builder);

        let main = builder.main();
        let local = main.row_slice(0).expect("Poseidon2 hash trace is empty");
        // Constrain the digest (first OUT elements of final state) to public values on last row.
        let public = builder.public_values().to_vec();
        let (public_inputs, rest) = public.split_at(POSEIDON2_INPUT_LEN);
        let (public_digest, public_packed) = rest.split_at(POSEIDON2_OUT);
        let mut builder = builder.when_last_row();
        let row = (*local).as_ref();
        let perm_width = num_cols::<
            POSEIDON2_WIDTH,
            POSEIDON2_SBOX_DEGREE,
            POSEIDON2_SBOX_REGISTERS,
            POSEIDON2_HALF_FULL_ROUNDS,
            POSEIDON2_PARTIAL_ROUNDS,
        >();

        for perm_idx in 0..POSEIDON2_NUM_PERMS {
            let start = perm_idx * perm_width;
            let perm_slice = &row[start..start + perm_width];
            let perm_cols: &Poseidon2Cols<
                AB::Var,
                POSEIDON2_WIDTH,
                POSEIDON2_SBOX_DEGREE,
                POSEIDON2_SBOX_REGISTERS,
                POSEIDON2_HALF_FULL_ROUNDS,
                POSEIDON2_PARTIAL_ROUNDS,
            > = perm_slice.borrow();
            let inputs = &perm_cols.inputs;
            if perm_idx == 0 {
                for i in 0..POSEIDON2_RATE {
                    builder.assert_eq(inputs[i].clone(), public_inputs[i]);
                }
                for i in POSEIDON2_RATE..POSEIDON2_WIDTH {
                    builder.assert_eq(inputs[i].clone(), BabyBear::ZERO);
                }
                continue;
            }

            let prev_start = (perm_idx - 1) * perm_width;
            let prev_slice = &row[prev_start..prev_start + perm_width];
            let prev_cols: &Poseidon2Cols<
                AB::Var,
                POSEIDON2_WIDTH,
                POSEIDON2_SBOX_DEGREE,
                POSEIDON2_SBOX_REGISTERS,
                POSEIDON2_HALF_FULL_ROUNDS,
                POSEIDON2_PARTIAL_ROUNDS,
            > = prev_slice.borrow();
            let prev_output = &prev_cols.ending_full_rounds[POSEIDON2_HALF_FULL_ROUNDS - 1].post;
            for i in 0..POSEIDON2_WIDTH {
                let msg_index = perm_idx * POSEIDON2_RATE + i;
                if i < POSEIDON2_RATE && msg_index < POSEIDON2_INPUT_LEN {
                    builder.assert_eq(inputs[i].clone(), public_inputs[msg_index]);
                } else {
                    builder.assert_eq(inputs[i].clone(), prev_output[i].clone());
                }
            }
        }

        let output_start = (POSEIDON2_NUM_PERMS - 1) * perm_width;
        let output_slice = &row[output_start..output_start + perm_width];
        let output_cols: &Poseidon2Cols<
            AB::Var,
            POSEIDON2_WIDTH,
            POSEIDON2_SBOX_DEGREE,
            POSEIDON2_SBOX_REGISTERS,
            POSEIDON2_HALF_FULL_ROUNDS,
            POSEIDON2_PARTIAL_ROUNDS,
        > = output_slice.borrow();
        let output = &output_cols.ending_full_rounds[POSEIDON2_HALF_FULL_ROUNDS - 1].post;
        for i in 0..POSEIDON2_OUT {
            builder.assert_eq(output[i].clone(), public_digest[i]);
        }

        let two_pow_31 = BabyBear::from_u32(1 << 31);
        builder.assert_eq(
            output[0].clone() + output[1].clone() * two_pow_31,
            public_packed[0],
        );
        builder.assert_eq(
            output[2].clone() + output[3].clone() * two_pow_31,
            public_packed[1],
        );
        builder.assert_eq(output[4].clone(), public_packed[2]);
    }
}

pub fn build_poseidon2_hash_trace(input: &[BabyBear]) -> Result<Poseidon2HashTrace, String> {
    if input.is_empty() {
        return Err("Poseidon2 hash input must be non-empty".to_string());
    }
    if input.len() != POSEIDON2_INPUT_LEN {
        return Err(format!(
            "Poseidon2 hash input length {} does not match expected {}",
            input.len(),
            POSEIDON2_INPUT_LEN
        ));
    }

    let perm = default_babybear_poseidon2_16();
    let mut state = [BabyBear::ZERO; POSEIDON2_WIDTH];
    let mut inputs = Vec::new();

    for chunk in input.chunks(POSEIDON2_RATE) {
        for (i, &value) in chunk.iter().enumerate() {
            state[i] = value;
        }
        inputs.push(state);
        perm.permute_mut(&mut state);
    }

    let digest: [BabyBear; POSEIDON2_OUT] = state[..POSEIDON2_OUT]
        .try_into()
        .map_err(|_| "Failed to extract Poseidon2 digest".to_string())?;

    let constants = RoundConstants::new(
        BABYBEAR_RC16_EXTERNAL_INITIAL,
        BABYBEAR_RC16_INTERNAL,
        BABYBEAR_RC16_EXTERNAL_FINAL,
    );
    let num_permutations = inputs.len();
    let trace = p3_poseidon2_air::generate_vectorized_trace_rows::<
        BabyBear,
        GenericPoseidon2LinearLayersBabyBear,
        POSEIDON2_WIDTH,
        POSEIDON2_SBOX_DEGREE,
        POSEIDON2_SBOX_REGISTERS,
        POSEIDON2_HALF_FULL_ROUNDS,
        POSEIDON2_PARTIAL_ROUNDS,
        POSEIDON2_NUM_PERMS,
    >(inputs, &constants, 0);

    Ok(Poseidon2HashTrace {
        trace,
        digest,
        num_permutations,
    })
}
