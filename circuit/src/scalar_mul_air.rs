//! AIR and trace builder for a scalar multiplication trace (scaffold).

use p3_air::{Air, AirBuilder, AirBuilderWithPublicValues, BaseAir, BaseAirWithPublicValues};
use p3_field::{Field, PrimeCharacteristicRing};
use p3_matrix::{dense::RowMajorMatrix, Matrix};

use crate::{scalar_to_bits, CircuitPoint, SignatureWitness};
use curve::{BaseField, KoalaBear};

pub const COORD_LIMBS: usize = 8;
pub const ACC_X_START: usize = 0;
pub const ACC_Y_START: usize = ACC_X_START + COORD_LIMBS;
pub const BASE_X_START: usize = ACC_Y_START + COORD_LIMBS;
pub const BASE_Y_START: usize = BASE_X_START + COORD_LIMBS;
pub const ADD_X_START: usize = BASE_Y_START + COORD_LIMBS;
pub const ADD_Y_START: usize = ADD_X_START + COORD_LIMBS;
pub const DOUBLE_X_START: usize = ADD_Y_START + COORD_LIMBS;
pub const DOUBLE_Y_START: usize = DOUBLE_X_START + COORD_LIMBS;
pub const ADD_NUM_START: usize = DOUBLE_Y_START + COORD_LIMBS;
pub const ADD_DEN_START: usize = ADD_NUM_START + COORD_LIMBS;
pub const ADD_INV_START: usize = ADD_DEN_START + COORD_LIMBS;
pub const ADD_SLOPE_START: usize = ADD_INV_START + COORD_LIMBS;
pub const DOUBLE_NUM_START: usize = ADD_SLOPE_START + COORD_LIMBS;
pub const DOUBLE_DEN_START: usize = DOUBLE_NUM_START + COORD_LIMBS;
pub const DOUBLE_INV_START: usize = DOUBLE_DEN_START + COORD_LIMBS;
pub const DOUBLE_SLOPE_START: usize = DOUBLE_INV_START + COORD_LIMBS;
pub const BIT_COL: usize = DOUBLE_SLOPE_START + COORD_LIMBS;
pub const ACC_INF_COL: usize = BIT_COL + 1;
pub const NUM_COLUMNS: usize = ACC_INF_COL + 1;

pub const PUBLIC_BASE_LIMBS: usize = COORD_LIMBS * 2;
pub const PUBLIC_OUT_LIMBS: usize = COORD_LIMBS * 2;

#[derive(Clone, Debug)]
pub struct ScalarMulTrace {
    pub trace: RowMajorMatrix<KoalaBear>,
}

#[derive(Clone, Debug)]
pub struct ScalarMulAir {
    pub num_rows: usize,
}

impl ScalarMulAir {
    pub fn new(num_rows: usize) -> Self {
        assert!(num_rows.is_power_of_two(), "num_rows must be power of 2");
        Self { num_rows }
    }
}

impl BaseAir<KoalaBear> for ScalarMulAir {
    fn width(&self) -> usize {
        NUM_COLUMNS
    }
}

impl BaseAirWithPublicValues<KoalaBear> for ScalarMulAir {
    fn num_public_values(&self) -> usize {
        PUBLIC_BASE_LIMBS + PUBLIC_OUT_LIMBS
    }
}

impl<AB> Air<AB> for ScalarMulAir
where
    AB: AirBuilder<F = KoalaBear> + AirBuilderWithPublicValues,
{
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local = main.row_slice(0).expect("scalar mul trace is empty");
        let row = (*local).as_ref();
        let public = builder.public_values().to_vec();
        let (base_public, out_public) = public.split_at(PUBLIC_BASE_LIMBS);

        let mut first = builder.when_first_row();
        for i in 0..COORD_LIMBS {
            first.assert_eq(row[BASE_X_START + i].clone(), base_public[i]);
            first.assert_eq(row[BASE_Y_START + i].clone(), base_public[i + COORD_LIMBS]);
            first.assert_eq(row[ACC_X_START + i].clone(), KoalaBear::ZERO);
            first.assert_eq(row[ACC_Y_START + i].clone(), KoalaBear::ZERO);
        }
        first.assert_bool(row[ACC_INF_COL].clone());

        let mut last = builder.when_last_row();
        for i in 0..COORD_LIMBS {
            last.assert_eq(row[ACC_X_START + i].clone(), out_public[i]);
            last.assert_eq(row[ACC_Y_START + i].clone(), out_public[i + COORD_LIMBS]);
        }

        let next_row = main.row_slice(1).expect("scalar mul trace row missing");
        let next_row = (*next_row).as_ref();
        eval_scalar_mul_core(builder, row, next_row, 0);
    }
}

pub(crate) fn eval_scalar_mul_core<AB: AirBuilder<F = KoalaBear>>(
    builder: &mut AB,
    row: &[AB::Var],
    next_row: &[AB::Var],
    offset: usize,
) {
    let bit = row[offset + BIT_COL].clone();
    builder.assert_bool(bit.clone());
    builder.assert_bool(row[offset + ACC_INF_COL].clone());

    let mut next = builder.when_transition();
    for i in 0..COORD_LIMBS {
        let acc = row[offset + ACC_X_START + i].clone();
        let acc_add = row[offset + ADD_X_START + i].clone();
        let acc_next = next_row[offset + ACC_X_START + i].clone();
        next.assert_eq(acc_next, acc.clone() + bit.clone() * (acc_add - acc));

        let accy = row[offset + ACC_Y_START + i].clone();
        let accy_add = row[offset + ADD_Y_START + i].clone();
        let accy_next = next_row[offset + ACC_Y_START + i].clone();
        next.assert_eq(accy_next, accy.clone() + bit.clone() * (accy_add - accy));

        let base_next = next_row[offset + BASE_X_START + i].clone();
        let base_double = row[offset + DOUBLE_X_START + i].clone();
        next.assert_eq(base_next, base_double);

        let basey_next = next_row[offset + BASE_Y_START + i].clone();
        let basey_double = row[offset + DOUBLE_Y_START + i].clone();
        next.assert_eq(basey_next, basey_double);
    }

    let acc_inf = row[offset + ACC_INF_COL].clone();
    let next_acc_inf = next_row[offset + ACC_INF_COL].clone();
    next.assert_bool(next_acc_inf.clone());
    let bit_expr: AB::Expr = bit.clone().into();
    let acc_inf_expr: AB::Expr = acc_inf.clone().into();
    next.assert_eq(
        next_acc_inf.clone(),
        acc_inf_expr.clone() * (AB::Expr::ONE - bit_expr),
    );

    let mut add_builder = builder.when(AB::Expr::ONE - acc_inf_expr);
    enforce_add_constraints(&mut add_builder, row, offset);
    enforce_on_curve(
        &mut add_builder,
        row,
        offset + ACC_X_START,
        offset + ACC_Y_START,
    );
    enforce_on_curve(
        &mut add_builder,
        row,
        offset + ADD_X_START,
        offset + ADD_Y_START,
    );

    let mut base_builder = builder.when(KoalaBear::ONE);
    enforce_double_constraints(&mut base_builder, row, offset);
    enforce_on_curve(
        &mut base_builder,
        row,
        offset + BASE_X_START,
        offset + BASE_Y_START,
    );
    enforce_on_curve(
        &mut base_builder,
        row,
        offset + DOUBLE_X_START,
        offset + DOUBLE_Y_START,
    );
}

pub fn build_scalar_mul_trace(base: &CircuitPoint, scalar_bits: &[bool]) -> ScalarMulTrace {
    let mut acc = CircuitPoint::infinity();
    let mut current = base.clone();
    let num_rows = scalar_bits.len().next_power_of_two();
    let mut trace = Vec::with_capacity(num_rows * NUM_COLUMNS);

    for row_idx in 0..num_rows {
        let mut row = vec![KoalaBear::ZERO; NUM_COLUMNS];
        let bit = if row_idx < scalar_bits.len() {
            scalar_bits[row_idx]
        } else {
            false
        };

        write_point(&mut row, ACC_X_START, &acc);
        write_point(&mut row, BASE_X_START, &current);

        let acc_plus = if acc.is_infinity {
            current.clone()
        } else {
            acc.add(&current)
        };
        let base_double = current.double();
        write_point(&mut row, ADD_X_START, &acc_plus);
        write_point(&mut row, DOUBLE_X_START, &base_double);

        let acc_inf = acc.is_infinity;
        row[ACC_INF_COL] = KoalaBear::from_u32(acc_inf as u32);

        if !acc_inf {
            let acc_x = coeffs_to_base(acc.x);
            let acc_y = coeffs_to_base(acc.y);
            let base_x = coeffs_to_base(current.x);
            let base_y = coeffs_to_base(current.y);

            let add_num = base_y - acc_y;
            let add_den = base_x - acc_x;
            let add_inv = add_den.inverse();
            let add_slope = add_num * add_inv;
            write_base(&mut row, ADD_NUM_START, add_num);
            write_base(&mut row, ADD_DEN_START, add_den);
            write_base(&mut row, ADD_INV_START, add_inv);
            write_base(&mut row, ADD_SLOPE_START, add_slope);
        }

        let base_x = coeffs_to_base(current.x);
        let base_y = coeffs_to_base(current.y);
        let double_num = {
            let three = KoalaBear::new(3);
            let a = coeffs_to_base([
                KoalaBear::ZERO,
                three,
                KoalaBear::ZERO,
                KoalaBear::ZERO,
                KoalaBear::ZERO,
                KoalaBear::ZERO,
                KoalaBear::ZERO,
                KoalaBear::ZERO,
            ]);
            let base_x2 = base_x * base_x;
            base_x2
                * coeffs_to_base([
                    three,
                    KoalaBear::ZERO,
                    KoalaBear::ZERO,
                    KoalaBear::ZERO,
                    KoalaBear::ZERO,
                    KoalaBear::ZERO,
                    KoalaBear::ZERO,
                    KoalaBear::ZERO,
                ])
                + a
        };
        let double_den = base_y + base_y;
        let double_inv = double_den.inverse();
        let double_slope = double_num * double_inv;
        write_base(&mut row, DOUBLE_NUM_START, double_num);
        write_base(&mut row, DOUBLE_DEN_START, double_den);
        write_base(&mut row, DOUBLE_INV_START, double_inv);
        write_base(&mut row, DOUBLE_SLOPE_START, double_slope);

        row[BIT_COL] = KoalaBear::from_u32(bit as u32);

        trace.extend_from_slice(&row);

        if bit {
            acc = acc_plus;
        }
        current = base_double;
    }

    ScalarMulTrace {
        trace: RowMajorMatrix::new(trace, NUM_COLUMNS),
    }
}

pub fn build_generator_mul_trace(witness: &SignatureWitness) -> ScalarMulTrace {
    let s_bits = scalar_to_bits(&witness.s.to_scalar_field());
    let generator = CircuitPoint::from_affine(&curve::Projective::generator().to_affine());
    build_scalar_mul_trace(&generator, &s_bits)
}

fn write_point(row: &mut [KoalaBear], start: usize, point: &CircuitPoint) {
    for i in 0..COORD_LIMBS {
        row[start + i] = point.x[i];
        row[start + i + COORD_LIMBS] = point.y[i];
    }
}

fn write_base(row: &mut [KoalaBear], start: usize, value: BaseField) {
    let coeffs: [KoalaBear; COORD_LIMBS] = unsafe { core::mem::transmute(value) };
    row[start..start + COORD_LIMBS].copy_from_slice(&coeffs[..COORD_LIMBS]);
}

pub(crate) fn enforce_add_constraints<AB: AirBuilder<F = KoalaBear>>(
    builder: &mut AB,
    row: &[AB::Var],
    offset: usize,
) {
    let acc_x = read_fp8::<AB>(row, offset + ACC_X_START);
    let acc_y = read_fp8::<AB>(row, offset + ACC_Y_START);
    let base_x = read_fp8::<AB>(row, offset + BASE_X_START);
    let base_y = read_fp8::<AB>(row, offset + BASE_Y_START);

    let add_num = read_fp8::<AB>(row, offset + ADD_NUM_START);
    let add_den = read_fp8::<AB>(row, offset + ADD_DEN_START);
    let add_inv = read_fp8::<AB>(row, offset + ADD_INV_START);
    let add_slope = read_fp8::<AB>(row, offset + ADD_SLOPE_START);
    let add_x = read_fp8::<AB>(row, offset + ADD_X_START);
    let add_y = read_fp8::<AB>(row, offset + ADD_Y_START);

    let num_expected = fp8_sub::<AB>(&base_y, &acc_y);
    let den_expected = fp8_sub::<AB>(&base_x, &acc_x);
    assert_fp8_eq(builder, &add_num, &num_expected);
    assert_fp8_eq(builder, &add_den, &den_expected);
    let one = fp8_one::<AB>();
    assert_fp8_eq(builder, &fp8_mul::<AB>(&add_den, &add_inv), &one);
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

pub(crate) fn enforce_double_constraints<AB: AirBuilder<F = KoalaBear>>(
    builder: &mut AB,
    row: &[AB::Var],
    offset: usize,
) {
    let base_x = read_fp8::<AB>(row, offset + BASE_X_START);
    let base_y = read_fp8::<AB>(row, offset + BASE_Y_START);

    let double_num = read_fp8::<AB>(row, offset + DOUBLE_NUM_START);
    let double_den = read_fp8::<AB>(row, offset + DOUBLE_DEN_START);
    let double_inv = read_fp8::<AB>(row, offset + DOUBLE_INV_START);
    let double_slope = read_fp8::<AB>(row, offset + DOUBLE_SLOPE_START);
    let double_x = read_fp8::<AB>(row, offset + DOUBLE_X_START);
    let double_y = read_fp8::<AB>(row, offset + DOUBLE_Y_START);

    let three = KoalaBear::from_u32(3);
    let a = fp8_a::<AB>();
    let x2 = fp8_mul::<AB>(&base_x, &base_x);
    let num_expected = fp8_add::<AB>(&fp8_mul_scalar::<AB>(&x2, three), &a);
    let den_expected = fp8_mul_scalar::<AB>(&base_y, KoalaBear::from_u32(2));

    assert_fp8_eq(builder, &double_num, &num_expected);
    assert_fp8_eq(builder, &double_den, &den_expected);
    assert_fp8_eq(
        builder,
        &fp8_mul::<AB>(&double_den, &double_inv),
        &fp8_one::<AB>(),
    );
    assert_fp8_eq(
        builder,
        &double_slope,
        &fp8_mul::<AB>(&double_num, &double_inv),
    );

    let slope2 = fp8_mul::<AB>(&double_slope, &double_slope);
    let x3 = fp8_sub::<AB>(
        &slope2,
        &fp8_mul_scalar::<AB>(&base_x, KoalaBear::from_u32(2)),
    );
    let y3 = fp8_sub::<AB>(
        &fp8_mul::<AB>(&double_slope, &fp8_sub::<AB>(&base_x, &x3)),
        &base_y,
    );
    assert_fp8_eq(builder, &double_x, &x3);
    assert_fp8_eq(builder, &double_y, &y3);
}

pub(crate) fn read_fp8<AB: AirBuilder<F = KoalaBear>>(
    row: &[AB::Var],
    start: usize,
) -> [AB::Expr; COORD_LIMBS] {
    core::array::from_fn(|i| row[start + i].clone().into())
}

pub(crate) fn assert_fp8_eq<AB: AirBuilder<F = KoalaBear>>(
    builder: &mut AB,
    lhs: &[AB::Expr; COORD_LIMBS],
    rhs: &[AB::Expr; COORD_LIMBS],
) {
    for i in 0..COORD_LIMBS {
        builder.assert_eq(lhs[i].clone(), rhs[i].clone());
    }
}

pub(crate) fn fp8_one<AB: AirBuilder<F = KoalaBear>>() -> [AB::Expr; COORD_LIMBS] {
    let mut out = [AB::Expr::ZERO; COORD_LIMBS];
    out[0] = AB::Expr::ONE;
    out
}

pub(crate) fn fp8_a<AB: AirBuilder<F = KoalaBear>>() -> [AB::Expr; COORD_LIMBS] {
    let mut out = [AB::Expr::ZERO; COORD_LIMBS];
    out[1] = AB::Expr::from(KoalaBear::from_u32(3));
    out
}

pub(crate) fn fp8_add<AB: AirBuilder<F = KoalaBear>>(
    a: &[AB::Expr; COORD_LIMBS],
    b: &[AB::Expr; COORD_LIMBS],
) -> [AB::Expr; COORD_LIMBS] {
    core::array::from_fn(|i| a[i].clone() + b[i].clone())
}

pub(crate) fn fp8_sub<AB: AirBuilder<F = KoalaBear>>(
    a: &[AB::Expr; COORD_LIMBS],
    b: &[AB::Expr; COORD_LIMBS],
) -> [AB::Expr; COORD_LIMBS] {
    core::array::from_fn(|i| a[i].clone() - b[i].clone())
}

pub(crate) fn fp8_mul_scalar<AB: AirBuilder<F = KoalaBear>>(
    a: &[AB::Expr; COORD_LIMBS],
    scalar: KoalaBear,
) -> [AB::Expr; COORD_LIMBS] {
    core::array::from_fn(|i| a[i].clone() * scalar)
}

pub(crate) fn fp8_mul<AB: AirBuilder<F = KoalaBear>>(
    a: &[AB::Expr; COORD_LIMBS],
    b: &[AB::Expr; COORD_LIMBS],
) -> [AB::Expr; COORD_LIMBS] {
    let mut t = vec![AB::Expr::ZERO; 2 * COORD_LIMBS - 1];
    for (i, a_i) in a.iter().enumerate().take(COORD_LIMBS) {
        for (j, b_j) in b.iter().enumerate().take(COORD_LIMBS) {
            let idx = i + j;
            t[idx] = t[idx].clone() + a_i.clone() * b_j.clone();
        }
    }

    let w = KoalaBear::from_u32(3);
    let mut out = [AB::Expr::ZERO; COORD_LIMBS];
    for k in 0..COORD_LIMBS {
        let mut acc = t[k].clone();
        if k + COORD_LIMBS < t.len() {
            acc += t[k + COORD_LIMBS].clone() * w;
        }
        out[k] = acc;
    }
    out
}

pub(crate) fn enforce_on_curve<AB: AirBuilder<F = KoalaBear>>(
    builder: &mut AB,
    row: &[AB::Var],
    x_start: usize,
    y_start: usize,
) {
    let x = read_fp8::<AB>(row, x_start);
    let y = read_fp8::<AB>(row, y_start);
    let y2 = fp8_mul::<AB>(&y, &y);
    let x2 = fp8_mul::<AB>(&x, &x);
    let x3 = fp8_mul::<AB>(&x2, &x);
    let ax = fp8_mul::<AB>(&x, &fp8_a::<AB>());
    let rhs = fp8_add::<AB>(&fp8_add::<AB>(&x3, &ax), &fp8_b::<AB>());
    assert_fp8_eq(builder, &y2, &rhs);
}

pub(crate) fn fp8_b<AB: AirBuilder<F = KoalaBear>>() -> [AB::Expr; COORD_LIMBS] {
    let mut out = [AB::Expr::ZERO; COORD_LIMBS];
    out[0] = AB::Expr::from(KoalaBear::from_u32(42639));
    out
}

fn coeffs_to_base(coeffs: [KoalaBear; COORD_LIMBS]) -> BaseField {
    unsafe { core::mem::transmute(coeffs) }
}
