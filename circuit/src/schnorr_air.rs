//! Schnorr verification AIR using two scalar-mul traces and a final add.

use p3_air::{
    Air, AirBuilder, AirBuilderWithPublicValues, BaseAir, BaseAirWithPublicValues, PairBuilder,
};
use p3_field::{Field, PrimeCharacteristicRing};
use p3_matrix::{dense::RowMajorMatrix, Matrix};

use crate::scalar_mul_air::{
    assert_fp8_eq, enforce_on_curve, fp8_a, fp8_add, fp8_mul, fp8_mul_scalar, fp8_one, fp8_sub,
    read_fp8, COORD_LIMBS,
};
use crate::{scalar_to_bits, CircuitPoint, SignatureWitness};
use curve::{BaseField, KoalaBear};

pub const SCHNORR_BASE_PUBLIC: usize = COORD_LIMBS * 2; // pk
pub const SCHNORR_R_PUBLIC: usize = COORD_LIMBS * 2; // R
pub const SCHNORR_PUBLIC_VALUES: usize = SCHNORR_BASE_PUBLIC + SCHNORR_R_PUBLIC;

pub const GS_OFFSET: usize = 0;
pub const DS_ACC_X_START: usize = 0;
pub const DS_ACC_Y_START: usize = DS_ACC_X_START + COORD_LIMBS;
pub const DS_PK_X_START: usize = DS_ACC_Y_START + COORD_LIMBS;
pub const DS_PK_Y_START: usize = DS_PK_X_START + COORD_LIMBS;
pub const DS_PK_DOUBLE_X_START: usize = DS_PK_Y_START + COORD_LIMBS;
pub const DS_PK_DOUBLE_Y_START: usize = DS_PK_DOUBLE_X_START + COORD_LIMBS;
pub const DS_PK_DOUBLE_NUM_START: usize = DS_PK_DOUBLE_Y_START + COORD_LIMBS;
pub const DS_PK_DOUBLE_DEN_START: usize = DS_PK_DOUBLE_NUM_START + COORD_LIMBS;
pub const DS_PK_DOUBLE_INV_START: usize = DS_PK_DOUBLE_DEN_START + COORD_LIMBS;
pub const DS_PK_DOUBLE_SLOPE_START: usize = DS_PK_DOUBLE_INV_START + COORD_LIMBS;
pub const DS_SUM_X_START: usize = DS_PK_DOUBLE_SLOPE_START + COORD_LIMBS;
pub const DS_SUM_Y_START: usize = DS_SUM_X_START + COORD_LIMBS;
pub const DS_SUM_NUM_START: usize = DS_SUM_Y_START + COORD_LIMBS;
pub const DS_SUM_DEN_START: usize = DS_SUM_NUM_START + COORD_LIMBS;
pub const DS_SUM_INV_START: usize = DS_SUM_DEN_START + COORD_LIMBS;
pub const DS_SUM_SLOPE_START: usize = DS_SUM_INV_START + COORD_LIMBS;
pub const DS_ADDEND_X_START: usize = DS_SUM_SLOPE_START + COORD_LIMBS;
pub const DS_ADDEND_Y_START: usize = DS_ADDEND_X_START + COORD_LIMBS;
pub const DS_ADD_X_START: usize = DS_ADDEND_Y_START + COORD_LIMBS;
pub const DS_ADD_Y_START: usize = DS_ADD_X_START + COORD_LIMBS;
pub const DS_ADD_NUM_START: usize = DS_ADD_Y_START + COORD_LIMBS;
pub const DS_ADD_DEN_START: usize = DS_ADD_NUM_START + COORD_LIMBS;
pub const DS_ADD_INV_START: usize = DS_ADD_DEN_START + COORD_LIMBS;
pub const DS_ADD_SLOPE_START: usize = DS_ADD_INV_START + COORD_LIMBS;
pub const DS_S_BIT_COL: usize = DS_ADD_SLOPE_START + COORD_LIMBS;
pub const DS_E_BIT_COL: usize = DS_S_BIT_COL + 1;
pub const DS_ACC_INF_COL: usize = DS_E_BIT_COL + 1;
pub const SCHNORR_COLUMNS: usize = DS_ACC_INF_COL + 1;
pub const GS_PREP_BASE_X_START: usize = 0;
pub const GS_PREP_BASE_Y_START: usize = GS_PREP_BASE_X_START + COORD_LIMBS;
pub const GS_PREP_COLS: usize = GS_PREP_BASE_Y_START + COORD_LIMBS;

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

    fn preprocessed_trace(&self) -> Option<RowMajorMatrix<KoalaBear>> {
        Some(build_gs_preprocessed_trace(self.num_rows))
    }
}

impl BaseAirWithPublicValues<KoalaBear> for SchnorrAir {
    fn num_public_values(&self) -> usize {
        SCHNORR_PUBLIC_VALUES
    }
}

impl<AB> Air<AB> for SchnorrAir
where
    AB: AirBuilder<F = KoalaBear> + AirBuilderWithPublicValues + PairBuilder,
{
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local = main.row_slice(0).expect("schnorr trace is empty");
        let row = (*local).as_ref();
        let next_row = main.row_slice(1).expect("schnorr next row missing");
        let next_row = (*next_row).as_ref();
        let preprocessed = builder.preprocessed();
        let preprocessed_row = preprocessed
            .row_slice(0)
            .expect("schnorr preprocessed is empty");
        let preprocessed_row = (*preprocessed_row).as_ref();

        let public = builder.public_values().to_vec();
        let (pk_public, r_public) = public.split_at(SCHNORR_BASE_PUBLIC);

        eval_double_scalar_core(builder, row, next_row, preprocessed_row, GS_OFFSET);

        let mut first = builder.when_first_row();
        for i in 0..COORD_LIMBS {
            first.assert_eq(row[DS_PK_X_START + i].clone(), pk_public[i]);
            first.assert_eq(row[DS_PK_Y_START + i].clone(), pk_public[i + COORD_LIMBS]);
        }

        // Bind preprocessed GS base to the generator on the first row.
        let generator = CircuitPoint::from_affine(&curve::Projective::generator().to_affine());
        for i in 0..COORD_LIMBS {
            first.assert_eq(
                preprocessed_row[GS_PREP_BASE_X_START + i].clone(),
                generator.x[i],
            );
            first.assert_eq(
                preprocessed_row[GS_PREP_BASE_Y_START + i].clone(),
                generator.y[i],
            );
        }

        let mut last = builder.when_last_row();
        // Bind accumulator to public R on the last row.
        for i in 0..COORD_LIMBS {
            last.assert_eq(row[DS_ACC_X_START + i].clone(), r_public[i]);
            last.assert_eq(row[DS_ACC_Y_START + i].clone(), r_public[i + COORD_LIMBS]);
        }
    }
}

pub fn build_schnorr_trace(witness: &SignatureWitness) -> SchnorrTrace {
    let s_bits = scalar_to_bits(&witness.s.to_scalar_field());
    let neg_e = -witness.challenge.to_scalar_field();
    let neg_e_bits = scalar_to_bits(&neg_e);
    let trace = build_double_scalar_trace(&s_bits, &neg_e_bits, &witness.public_key);

    SchnorrTrace {
        trace: RowMajorMatrix::new(trace, SCHNORR_COLUMNS),
    }
}

fn build_gs_preprocessed_trace(num_rows: usize) -> RowMajorMatrix<KoalaBear> {
    let mut trace = Vec::with_capacity(num_rows * GS_PREP_COLS);
    let mut current = CircuitPoint::from_affine(&curve::Projective::generator().to_affine());

    for _ in 0..num_rows {
        let mut row = vec![KoalaBear::ZERO; GS_PREP_COLS];
        write_preprocessed_point(&mut row, GS_PREP_BASE_X_START, &current);
        trace.extend_from_slice(&row);
        current = current.double();
    }

    RowMajorMatrix::new(trace, GS_PREP_COLS)
}

fn write_preprocessed_point(row: &mut [KoalaBear], start: usize, point: &CircuitPoint) {
    for i in 0..COORD_LIMBS {
        row[start + i] = point.x[i];
        row[start + i + COORD_LIMBS] = point.y[i];
    }
}

fn build_double_scalar_trace(
    s_bits: &[bool],
    e_bits: &[bool],
    pk: &CircuitPoint,
) -> Vec<KoalaBear> {
    let mut acc = CircuitPoint::infinity();
    let mut acc_inf = true;
    let mut pk_current = pk.clone();
    let mut g_current = CircuitPoint::from_affine(&curve::Projective::generator().to_affine());
    let num_rows = s_bits.len().max(e_bits.len()).next_power_of_two();
    let mut trace = Vec::with_capacity(num_rows * SCHNORR_COLUMNS);

    for row_idx in 0..num_rows {
        let mut row = vec![KoalaBear::ZERO; SCHNORR_COLUMNS];
        let s_bit = if row_idx < s_bits.len() {
            s_bits[row_idx]
        } else {
            false
        };
        let e_bit = if row_idx < e_bits.len() {
            e_bits[row_idx]
        } else {
            false
        };

        write_point(&mut row, DS_ACC_X_START, &acc);
        write_point(&mut row, DS_PK_X_START, &pk_current);

        let pk_double = pk_current.double();
        write_point(&mut row, DS_PK_DOUBLE_X_START, &pk_double);
        fill_double_intermediates(&mut row, DS_PK_DOUBLE_NUM_START, &pk_current);

        let sum = pk_current.add(&g_current);
        write_point(&mut row, DS_SUM_X_START, &sum);
        fill_add_intermediates(&mut row, DS_SUM_NUM_START, &pk_current, &g_current);

        let addend = match (s_bit, e_bit) {
            (false, false) => CircuitPoint::infinity(),
            (true, false) => g_current.clone(),
            (false, true) => pk_current.clone(),
            (true, true) => sum.clone(),
        };
        write_point(&mut row, DS_ADDEND_X_START, &addend);

        let addend_inf = addend.is_infinity;

        if !acc_inf && !addend_inf {
            let acc_add = acc.add(&addend);
            write_point(&mut row, DS_ADD_X_START, &acc_add);
            fill_add_intermediates(&mut row, DS_ADD_NUM_START, &acc, &addend);
        } else if acc_inf && !addend_inf {
            write_point(&mut row, DS_ADD_X_START, &addend);
        } else {
            write_point(&mut row, DS_ADD_X_START, &acc);
        }

        row[DS_S_BIT_COL] = KoalaBear::from_u32(s_bit as u32);
        row[DS_E_BIT_COL] = KoalaBear::from_u32(e_bit as u32);
        row[DS_ACC_INF_COL] = KoalaBear::from_u32(acc_inf as u32);

        trace.extend_from_slice(&row);

        if acc_inf {
            if !addend_inf {
                acc = addend;
                acc_inf = false;
            }
        } else if !addend_inf {
            acc = acc.add(&addend);
        }

        pk_current = pk_double;
        g_current = g_current.double();
    }

    trace
}

fn write_point(row: &mut [KoalaBear], start: usize, point: &CircuitPoint) {
    for i in 0..COORD_LIMBS {
        row[start + i] = point.x[i];
        row[start + i + COORD_LIMBS] = point.y[i];
    }
}

fn fill_add_intermediates(
    row: &mut [KoalaBear],
    start: usize,
    acc: &CircuitPoint,
    base: &CircuitPoint,
) {
    let acc_x = coeffs_to_base(acc.x);
    let acc_y = coeffs_to_base(acc.y);
    let base_x = coeffs_to_base(base.x);
    let base_y = coeffs_to_base(base.y);

    let add_num = base_y - acc_y;
    let add_den = base_x - acc_x;
    let add_inv = add_den.inverse();
    let add_slope = add_num * add_inv;

    write_base(row, start, add_num);
    write_base(row, start + COORD_LIMBS, add_den);
    write_base(row, start + 2 * COORD_LIMBS, add_inv);
    write_base(row, start + 3 * COORD_LIMBS, add_slope);
}

fn fill_double_intermediates(row: &mut [KoalaBear], start: usize, base: &CircuitPoint) {
    let base_x = coeffs_to_base(base.x);
    let base_y = coeffs_to_base(base.y);
    let three = KoalaBear::from_u32(3);
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
    let double_num = base_x2
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
        + a;
    let double_den = base_y + base_y;
    let double_inv = double_den.inverse();
    let double_slope = double_num * double_inv;

    write_base(row, start, double_num);
    write_base(row, start + COORD_LIMBS, double_den);
    write_base(row, start + 2 * COORD_LIMBS, double_inv);
    write_base(row, start + 3 * COORD_LIMBS, double_slope);
}

fn enforce_add_constraints_with_base<AB: AirBuilder<F = KoalaBear>>(
    builder: &mut AB,
    acc_x: &[AB::Expr; COORD_LIMBS],
    acc_y: &[AB::Expr; COORD_LIMBS],
    base_x: &[AB::Expr; COORD_LIMBS],
    base_y: &[AB::Expr; COORD_LIMBS],
    add_num: &[AB::Expr; COORD_LIMBS],
    add_den: &[AB::Expr; COORD_LIMBS],
    add_inv: &[AB::Expr; COORD_LIMBS],
    add_slope: &[AB::Expr; COORD_LIMBS],
    add_x: &[AB::Expr; COORD_LIMBS],
    add_y: &[AB::Expr; COORD_LIMBS],
) {
    let num_expected = fp8_sub::<AB>(base_y, acc_y);
    let den_expected = fp8_sub::<AB>(base_x, acc_x);
    assert_fp8_eq(builder, add_num, &num_expected);
    assert_fp8_eq(builder, add_den, &den_expected);
    assert_fp8_eq(builder, &fp8_mul::<AB>(add_den, add_inv), &fp8_one::<AB>());
    assert_fp8_eq(builder, add_slope, &fp8_mul::<AB>(add_num, add_inv));

    let slope2 = fp8_mul::<AB>(add_slope, add_slope);
    let x3 = fp8_sub::<AB>(&fp8_sub::<AB>(&slope2, acc_x), base_x);
    let y3 = fp8_sub::<AB>(&fp8_mul::<AB>(add_slope, &fp8_sub::<AB>(acc_x, &x3)), acc_y);
    assert_fp8_eq(builder, add_x, &x3);
    assert_fp8_eq(builder, add_y, &y3);
}

fn enforce_double_constraints_with_base<AB: AirBuilder<F = KoalaBear>>(
    builder: &mut AB,
    base_x: &[AB::Expr; COORD_LIMBS],
    base_y: &[AB::Expr; COORD_LIMBS],
    double_num: &[AB::Expr; COORD_LIMBS],
    double_den: &[AB::Expr; COORD_LIMBS],
    double_inv: &[AB::Expr; COORD_LIMBS],
    double_slope: &[AB::Expr; COORD_LIMBS],
    double_x: &[AB::Expr; COORD_LIMBS],
    double_y: &[AB::Expr; COORD_LIMBS],
) {
    let three = KoalaBear::from_u32(3);
    let a = fp8_a::<AB>();
    let x2 = fp8_mul::<AB>(base_x, base_x);
    let num_expected = fp8_add::<AB>(&fp8_mul_scalar::<AB>(&x2, three), &a);
    let den_expected = fp8_mul_scalar::<AB>(base_y, KoalaBear::from_u32(2));

    assert_fp8_eq(builder, double_num, &num_expected);
    assert_fp8_eq(builder, double_den, &den_expected);
    assert_fp8_eq(
        builder,
        &fp8_mul::<AB>(double_den, double_inv),
        &fp8_one::<AB>(),
    );
    assert_fp8_eq(
        builder,
        double_slope,
        &fp8_mul::<AB>(double_num, double_inv),
    );

    let slope2 = fp8_mul::<AB>(double_slope, double_slope);
    let x3 = fp8_sub::<AB>(
        &slope2,
        &fp8_mul_scalar::<AB>(base_x, KoalaBear::from_u32(2)),
    );
    let y3 = fp8_sub::<AB>(
        &fp8_mul::<AB>(double_slope, &fp8_sub::<AB>(base_x, &x3)),
        base_y,
    );
    assert_fp8_eq(builder, double_x, &x3);
    assert_fp8_eq(builder, double_y, &y3);
}

pub(crate) fn eval_double_scalar_core<AB: AirBuilder<F = KoalaBear>>(
    builder: &mut AB,
    row: &[AB::Var],
    next_row: &[AB::Var],
    preprocessed_row: &[AB::Var],
    offset: usize,
) {
    let s_bit = row[offset + DS_S_BIT_COL].clone();
    let e_bit = row[offset + DS_E_BIT_COL].clone();
    builder.assert_bool(s_bit.clone());
    builder.assert_bool(e_bit.clone());
    builder.assert_bool(row[offset + DS_ACC_INF_COL].clone());

    let addend_inf_expr: AB::Expr =
        (AB::Expr::ONE - s_bit.clone().into()) * (AB::Expr::ONE - e_bit.clone().into());
    let add_sel = AB::Expr::ONE - addend_inf_expr.clone();

    let mut next = builder.when_transition();
    for i in 0..COORD_LIMBS {
        let acc = row[offset + DS_ACC_X_START + i].clone();
        let acc_add = row[offset + DS_ADD_X_START + i].clone();
        let acc_next = next_row[offset + DS_ACC_X_START + i].clone();
        next.assert_eq(acc_next, acc.clone() + add_sel.clone() * (acc_add - acc));

        let accy = row[offset + DS_ACC_Y_START + i].clone();
        let accy_add = row[offset + DS_ADD_Y_START + i].clone();
        let accy_next = next_row[offset + DS_ACC_Y_START + i].clone();
        next.assert_eq(
            accy_next,
            accy.clone() + add_sel.clone() * (accy_add - accy),
        );

        let pk_next = next_row[offset + DS_PK_X_START + i].clone();
        let pk_double = row[offset + DS_PK_DOUBLE_X_START + i].clone();
        next.assert_eq(pk_next, pk_double);

        let pky_next = next_row[offset + DS_PK_Y_START + i].clone();
        let pky_double = row[offset + DS_PK_DOUBLE_Y_START + i].clone();
        next.assert_eq(pky_next, pky_double);
    }

    let acc_inf = row[offset + DS_ACC_INF_COL].clone();
    let next_acc_inf = next_row[offset + DS_ACC_INF_COL].clone();
    next.assert_bool(next_acc_inf.clone());
    let acc_inf_expr: AB::Expr = acc_inf.clone().into();
    next.assert_eq(next_acc_inf.clone(), acc_inf_expr.clone() * addend_inf_expr);

    let g_x = read_fp8::<AB>(preprocessed_row, GS_PREP_BASE_X_START);
    let g_y = read_fp8::<AB>(preprocessed_row, GS_PREP_BASE_Y_START);
    let pk_x = read_fp8::<AB>(row, offset + DS_PK_X_START);
    let pk_y = read_fp8::<AB>(row, offset + DS_PK_Y_START);
    let pk_double_num = read_fp8::<AB>(row, offset + DS_PK_DOUBLE_NUM_START);
    let pk_double_den = read_fp8::<AB>(row, offset + DS_PK_DOUBLE_DEN_START);
    let pk_double_inv = read_fp8::<AB>(row, offset + DS_PK_DOUBLE_INV_START);
    let pk_double_slope = read_fp8::<AB>(row, offset + DS_PK_DOUBLE_SLOPE_START);
    let pk_double_x = read_fp8::<AB>(row, offset + DS_PK_DOUBLE_X_START);
    let pk_double_y = read_fp8::<AB>(row, offset + DS_PK_DOUBLE_Y_START);
    enforce_double_constraints_with_base(
        builder,
        &pk_x,
        &pk_y,
        &pk_double_num,
        &pk_double_den,
        &pk_double_inv,
        &pk_double_slope,
        &pk_double_x,
        &pk_double_y,
    );

    let sum_x = read_fp8::<AB>(row, offset + DS_SUM_X_START);
    let sum_y = read_fp8::<AB>(row, offset + DS_SUM_Y_START);
    let sum_num = read_fp8::<AB>(row, offset + DS_SUM_NUM_START);
    let sum_den = read_fp8::<AB>(row, offset + DS_SUM_DEN_START);
    let sum_inv = read_fp8::<AB>(row, offset + DS_SUM_INV_START);
    let sum_slope = read_fp8::<AB>(row, offset + DS_SUM_SLOPE_START);
    enforce_add_constraints_with_base(
        builder, &pk_x, &pk_y, &g_x, &g_y, &sum_num, &sum_den, &sum_inv, &sum_slope, &sum_x, &sum_y,
    );

    let addend_x = read_fp8::<AB>(row, offset + DS_ADDEND_X_START);
    let addend_y = read_fp8::<AB>(row, offset + DS_ADDEND_Y_START);

    let s_expr: AB::Expr = s_bit.clone().into();
    let e_expr: AB::Expr = e_bit.clone().into();
    let sel00 = (AB::Expr::ONE - s_expr.clone()) * (AB::Expr::ONE - e_expr.clone());
    let sel10 = s_expr.clone() * (AB::Expr::ONE - e_expr.clone());
    let sel01 = (AB::Expr::ONE - s_expr.clone()) * e_expr.clone();
    let sel11 = s_expr * e_expr;
    for i in 0..COORD_LIMBS {
        builder.assert_eq(
            addend_x[i].clone(),
            sel00.clone() * AB::Expr::ZERO
                + sel10.clone() * g_x[i].clone()
                + sel01.clone() * pk_x[i].clone()
                + sel11.clone() * sum_x[i].clone(),
        );
        builder.assert_eq(
            addend_y[i].clone(),
            sel00.clone() * AB::Expr::ZERO
                + sel10.clone() * g_y[i].clone()
                + sel01.clone() * pk_y[i].clone()
                + sel11.clone() * sum_y[i].clone(),
        );
    }

    let mut add_builder = builder.when((AB::Expr::ONE - acc_inf_expr.clone()) * add_sel.clone());
    let acc_x = read_fp8::<AB>(row, offset + DS_ACC_X_START);
    let acc_y = read_fp8::<AB>(row, offset + DS_ACC_Y_START);
    let add_num = read_fp8::<AB>(row, offset + DS_ADD_NUM_START);
    let add_den = read_fp8::<AB>(row, offset + DS_ADD_DEN_START);
    let add_inv = read_fp8::<AB>(row, offset + DS_ADD_INV_START);
    let add_slope = read_fp8::<AB>(row, offset + DS_ADD_SLOPE_START);
    let add_x = read_fp8::<AB>(row, offset + DS_ADD_X_START);
    let add_y = read_fp8::<AB>(row, offset + DS_ADD_Y_START);
    enforce_add_constraints_with_base(
        &mut add_builder,
        &acc_x,
        &acc_y,
        &addend_x,
        &addend_y,
        &add_num,
        &add_den,
        &add_inv,
        &add_slope,
        &add_x,
        &add_y,
    );

    let mut init_builder = builder.when(acc_inf_expr.clone() * add_sel.clone());
    for i in 0..COORD_LIMBS {
        init_builder.assert_eq(
            row[offset + DS_ADD_X_START + i].clone(),
            row[offset + DS_ADDEND_X_START + i].clone(),
        );
        init_builder.assert_eq(
            row[offset + DS_ADD_Y_START + i].clone(),
            row[offset + DS_ADDEND_Y_START + i].clone(),
        );
    }

    let mut acc_curve_builder = builder.when(AB::Expr::ONE - acc_inf_expr);
    enforce_on_curve(
        &mut acc_curve_builder,
        row,
        offset + DS_ACC_X_START,
        offset + DS_ACC_Y_START,
    );
    let mut addend_curve_builder = builder.when(add_sel);
    enforce_on_curve(
        &mut addend_curve_builder,
        row,
        offset + DS_ADDEND_X_START,
        offset + DS_ADDEND_Y_START,
    );
    let mut base_curve_builder = builder.when(KoalaBear::ONE);
    enforce_on_curve(
        &mut base_curve_builder,
        row,
        offset + DS_SUM_X_START,
        offset + DS_SUM_Y_START,
    );
    enforce_on_curve(
        &mut base_curve_builder,
        row,
        offset + DS_PK_X_START,
        offset + DS_PK_Y_START,
    );
    enforce_on_curve(
        &mut base_curve_builder,
        preprocessed_row,
        GS_PREP_BASE_X_START,
        GS_PREP_BASE_Y_START,
    );
}

fn coeffs_to_base(coeffs: [KoalaBear; COORD_LIMBS]) -> BaseField {
    unsafe { core::mem::transmute(coeffs) }
}

fn write_base(row: &mut [KoalaBear], start: usize, value: BaseField) {
    let coeffs: [KoalaBear; COORD_LIMBS] = unsafe { core::mem::transmute(value) };
    row[start..start + COORD_LIMBS].copy_from_slice(&coeffs[..COORD_LIMBS]);
}
