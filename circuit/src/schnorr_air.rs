//! Schnorr verification AIR using two scalar-mul traces and a final add.

use p3_air::{Air, AirBuilder, AirBuilderWithPublicValues, BaseAir, BaseAirWithPublicValues};
use p3_field::{Field, PrimeCharacteristicRing};
use p3_matrix::{dense::RowMajorMatrix, Matrix};

use crate::scalar_mul_air::{
    assert_fp8_eq, enforce_on_curve, eval_scalar_mul_core, fp8_mul, fp8_one, fp8_sub, read_fp8,
    ACC_X_START, ACC_Y_START, BASE_X_START, BASE_Y_START, COORD_LIMBS, NUM_COLUMNS,
};
use crate::{scalar_to_bits, CircuitPoint, SignatureWitness};
use curve::{BaseField, KoalaBear};

pub const SCHNORR_BASE_PUBLIC: usize = COORD_LIMBS * 2; // pk
pub const SCHNORR_R_PUBLIC: usize = COORD_LIMBS * 2; // R
pub const SCHNORR_PUBLIC_VALUES: usize = SCHNORR_BASE_PUBLIC + SCHNORR_R_PUBLIC;

pub const GS_OFFSET: usize = 0;
pub const PKE_OFFSET: usize = NUM_COLUMNS;
pub const ADD_OFFSET: usize = NUM_COLUMNS * 2;

pub const ADD_ACC_X_START: usize = 0;
pub const ADD_ACC_Y_START: usize = ADD_ACC_X_START + COORD_LIMBS;
pub const ADD_BASE_X_START: usize = ADD_ACC_Y_START + COORD_LIMBS;
pub const ADD_BASE_Y_START: usize = ADD_BASE_X_START + COORD_LIMBS;
pub const ADD_OUT_X_START: usize = ADD_BASE_Y_START + COORD_LIMBS;
pub const ADD_OUT_Y_START: usize = ADD_OUT_X_START + COORD_LIMBS;
pub const ADD_NUM_START: usize = ADD_OUT_Y_START + COORD_LIMBS;
pub const ADD_DEN_START: usize = ADD_NUM_START + COORD_LIMBS;
pub const ADD_INV_START: usize = ADD_DEN_START + COORD_LIMBS;
pub const ADD_SLOPE_START: usize = ADD_INV_START + COORD_LIMBS;
pub const ADD_BIT_COL: usize = ADD_SLOPE_START + COORD_LIMBS;
pub const ADD_ACC_INF_COL: usize = ADD_BIT_COL + 1;
pub const ADD_BLOCK_COLS: usize = ADD_ACC_INF_COL + 1;

pub const SCHNORR_COLUMNS: usize = NUM_COLUMNS * 2 + ADD_BLOCK_COLS;

#[derive(Clone, Debug)]
pub struct SchnorrTrace {
    pub trace: RowMajorMatrix<KoalaBear>,
}

#[derive(Clone, Debug)]
pub struct SchnorrAir {
    pub num_rows: usize,
}

impl SchnorrAir {
    pub fn new(num_rows: usize) -> Self {
        assert!(num_rows.is_power_of_two(), "num_rows must be power of 2");
        Self { num_rows }
    }
}

impl BaseAir<KoalaBear> for SchnorrAir {
    fn width(&self) -> usize {
        SCHNORR_COLUMNS
    }
}

impl BaseAirWithPublicValues<KoalaBear> for SchnorrAir {
    fn num_public_values(&self) -> usize {
        SCHNORR_PUBLIC_VALUES
    }
}

impl<AB> Air<AB> for SchnorrAir
where
    AB: AirBuilder<F = KoalaBear> + AirBuilderWithPublicValues,
{
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local = main.row_slice(0).expect("schnorr trace is empty");
        let row = (*local).as_ref();
        let next_row = main.row_slice(1).expect("schnorr next row missing");
        let next_row = (*next_row).as_ref();

        let public = builder.public_values().to_vec();
        let (pk_public, r_public) = public.split_at(SCHNORR_BASE_PUBLIC);

        eval_scalar_mul_core(builder, row, next_row, GS_OFFSET);
        eval_scalar_mul_core(builder, row, next_row, PKE_OFFSET);

        let mut first = builder.when_first_row();
        for i in 0..COORD_LIMBS {
            first.assert_eq(row[PKE_OFFSET + BASE_X_START + i].clone(), pk_public[i]);
            first.assert_eq(
                row[PKE_OFFSET + BASE_Y_START + i].clone(),
                pk_public[i + COORD_LIMBS],
            );
        }

        // Enforce generator base for GS block.
        let generator = CircuitPoint::from_affine(&curve::Projective::generator().to_affine());
        for i in 0..COORD_LIMBS {
            first.assert_eq(row[GS_OFFSET + BASE_X_START + i].clone(), generator.x[i]);
            first.assert_eq(row[GS_OFFSET + BASE_Y_START + i].clone(), generator.y[i]);
        }

        let mut last = builder.when_last_row();
        // Bind R (acc of add block) to public R.
        for i in 0..COORD_LIMBS {
            last.assert_eq(row[ADD_OFFSET + ADD_ACC_X_START + i].clone(), r_public[i]);
            last.assert_eq(
                row[ADD_OFFSET + ADD_ACC_Y_START + i].clone(),
                r_public[i + COORD_LIMBS],
            );
        }

        // Bind add base to pk*e accumulator and add result to g*s accumulator.
        for i in 0..COORD_LIMBS {
            last.assert_eq(
                row[ADD_OFFSET + ADD_BASE_X_START + i].clone(),
                row[PKE_OFFSET + ACC_X_START + i].clone(),
            );
            last.assert_eq(
                row[ADD_OFFSET + ADD_BASE_Y_START + i].clone(),
                row[PKE_OFFSET + ACC_Y_START + i].clone(),
            );
            last.assert_eq(
                row[ADD_OFFSET + ADD_OUT_X_START + i].clone(),
                row[GS_OFFSET + ACC_X_START + i].clone(),
            );
            last.assert_eq(
                row[ADD_OFFSET + ADD_OUT_Y_START + i].clone(),
                row[GS_OFFSET + ACC_Y_START + i].clone(),
            );
        }

        // Force the add block bit to 1 and acc infinity to 0 on last row.
        last.assert_eq(row[ADD_OFFSET + ADD_BIT_COL].clone(), KoalaBear::ONE);
        last.assert_eq(row[ADD_OFFSET + ADD_ACC_INF_COL].clone(), KoalaBear::ZERO);

        // Enforce add constraints on last row for R + pk*e.
        let mut add_builder = builder.when_last_row();
        enforce_add_constraints_compact(&mut add_builder, row, ADD_OFFSET);
        enforce_on_curve(
            &mut add_builder,
            row,
            ADD_OFFSET + ADD_ACC_X_START,
            ADD_OFFSET + ADD_ACC_Y_START,
        );
        enforce_on_curve(
            &mut add_builder,
            row,
            ADD_OFFSET + ADD_BASE_X_START,
            ADD_OFFSET + ADD_BASE_Y_START,
        );
        enforce_on_curve(
            &mut add_builder,
            row,
            ADD_OFFSET + ADD_OUT_X_START,
            ADD_OFFSET + ADD_OUT_Y_START,
        );
    }
}

pub fn build_schnorr_trace(witness: &SignatureWitness) -> SchnorrTrace {
    let s_bits = scalar_to_bits(&witness.s.to_scalar_field());
    let e_bits = scalar_to_bits(&witness.challenge.to_scalar_field());

    let generator = CircuitPoint::from_affine(&curve::Projective::generator().to_affine());
    let gs_trace = crate::build_scalar_mul_trace(&generator, &s_bits);
    let pke_trace = crate::build_scalar_mul_trace(&witness.public_key, &e_bits);

    let height = gs_trace.trace.height().max(pke_trace.trace.height());
    let mut trace = Vec::with_capacity(height * SCHNORR_COLUMNS);

    for row_idx in 0..height {
        let mut row = vec![KoalaBear::ZERO; SCHNORR_COLUMNS];
        let gs_row = gs_trace.trace.row_slice(row_idx).unwrap();
        let pke_row = pke_trace.trace.row_slice(row_idx).unwrap();

        row[GS_OFFSET..GS_OFFSET + NUM_COLUMNS].copy_from_slice(&gs_row);
        row[PKE_OFFSET..PKE_OFFSET + NUM_COLUMNS].copy_from_slice(&pke_row);

        if row_idx == height - 1 {
            fill_add_block(&mut row, witness, &gs_row, &pke_row);
        }

        trace.extend_from_slice(&row);
    }

    SchnorrTrace {
        trace: RowMajorMatrix::new(trace, SCHNORR_COLUMNS),
    }
}

fn fill_add_block(
    row: &mut [KoalaBear],
    witness: &SignatureWitness,
    gs_row: &[KoalaBear],
    pke_row: &[KoalaBear],
) {
    let r = &witness.r;
    for i in 0..COORD_LIMBS {
        row[ADD_OFFSET + ADD_ACC_X_START + i] = r.x[i];
        row[ADD_OFFSET + ADD_ACC_Y_START + i] = r.y[i];
        row[ADD_OFFSET + ADD_BASE_X_START + i] = pke_row[ACC_X_START + i];
        row[ADD_OFFSET + ADD_BASE_Y_START + i] = pke_row[ACC_Y_START + i];
        row[ADD_OFFSET + ADD_OUT_X_START + i] = gs_row[ACC_X_START + i];
        row[ADD_OFFSET + ADD_OUT_Y_START + i] = gs_row[ACC_Y_START + i];
    }

    // Compute add intermediates using base field arithmetic.
    let acc_x = coeffs_to_base(r.x);
    let acc_y = coeffs_to_base(r.y);
    let base_x = coeffs_from_row(pke_row, ACC_X_START);
    let base_y = coeffs_from_row(pke_row, ACC_Y_START);

    let add_num = base_y - acc_y;
    let add_den = base_x - acc_x;
    let add_inv = add_den.inverse();
    let add_slope = add_num * add_inv;

    write_base(row, ADD_OFFSET + ADD_NUM_START, add_num);
    write_base(row, ADD_OFFSET + ADD_DEN_START, add_den);
    write_base(row, ADD_OFFSET + ADD_INV_START, add_inv);
    write_base(row, ADD_OFFSET + ADD_SLOPE_START, add_slope);
    row[ADD_OFFSET + ADD_BIT_COL] = KoalaBear::ONE;
    row[ADD_OFFSET + ADD_ACC_INF_COL] = KoalaBear::ZERO;
}

fn enforce_add_constraints_compact<AB: AirBuilder<F = KoalaBear>>(
    builder: &mut AB,
    row: &[AB::Var],
    offset: usize,
) {
    let acc_x = read_fp8::<AB>(row, offset + ADD_ACC_X_START);
    let acc_y = read_fp8::<AB>(row, offset + ADD_ACC_Y_START);
    let base_x = read_fp8::<AB>(row, offset + ADD_BASE_X_START);
    let base_y = read_fp8::<AB>(row, offset + ADD_BASE_Y_START);

    let add_num = read_fp8::<AB>(row, offset + ADD_NUM_START);
    let add_den = read_fp8::<AB>(row, offset + ADD_DEN_START);
    let add_inv = read_fp8::<AB>(row, offset + ADD_INV_START);
    let add_slope = read_fp8::<AB>(row, offset + ADD_SLOPE_START);
    let add_x = read_fp8::<AB>(row, offset + ADD_OUT_X_START);
    let add_y = read_fp8::<AB>(row, offset + ADD_OUT_Y_START);

    let num_expected = fp8_sub::<AB>(&base_y, &acc_y);
    let den_expected = fp8_sub::<AB>(&base_x, &acc_x);
    assert_fp8_eq(builder, &add_num, &num_expected);
    assert_fp8_eq(builder, &add_den, &den_expected);
    assert_fp8_eq(
        builder,
        &fp8_mul::<AB>(&add_den, &add_inv),
        &fp8_one::<AB>(),
    );
    assert_fp8_eq(builder, &add_slope, &fp8_mul::<AB>(&add_num, &add_inv));

    let slope2 = fp8_mul::<AB>(&add_slope, &add_slope);
    let x3 = fp8_sub::<AB>(&fp8_sub::<AB>(&slope2, &acc_x), &base_x);
    let y3 = fp8_sub::<AB>(
        &fp8_mul::<AB>(&add_slope, &fp8_sub::<AB>(&acc_x, &x3)),
        &acc_y,
    );
    assert_fp8_eq(builder, &add_x, &x3);
    assert_fp8_eq(builder, &add_y, &y3);
}

fn coeffs_to_base(coeffs: [KoalaBear; COORD_LIMBS]) -> BaseField {
    unsafe { core::mem::transmute(coeffs) }
}

fn coeffs_from_row(row: &[KoalaBear], start: usize) -> BaseField {
    let mut coeffs = [KoalaBear::ZERO; COORD_LIMBS];
    coeffs.copy_from_slice(&row[start..start + COORD_LIMBS]);
    unsafe { core::mem::transmute(coeffs) }
}

fn write_base(row: &mut [KoalaBear], start: usize, value: BaseField) {
    let coeffs: [KoalaBear; COORD_LIMBS] = unsafe { core::mem::transmute(value) };
    row[start..start + COORD_LIMBS].copy_from_slice(&coeffs[..COORD_LIMBS]);
}
