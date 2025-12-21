use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use curve::{Affine, Group, Projective, RandomField, ScalarField};
use rand::rngs::StdRng;
use rand::SeedableRng;

fn random_scalar(rng: &mut StdRng) -> ScalarField {
    ScalarField::random(rng)
}

fn bench_affine_double(c: &mut Criterion) {
    let g = Affine::generator();
    c.bench_function("affine_double", |bencher| {
        bencher.iter(|| black_box(black_box(g).double()))
    });
}

fn bench_projective_double(c: &mut Criterion) {
    let g = Projective::generator();
    c.bench_function("projective_double", |bencher| {
        bencher.iter(|| black_box(black_box(g).double()))
    });
}

fn bench_affine_add(c: &mut Criterion) {
    let g = Affine::generator();
    let h = Affine::generator_pedersen();
    c.bench_function("affine_add", |bencher| {
        bencher.iter(|| black_box(black_box(g) + black_box(h)))
    });
}

fn bench_projective_add(c: &mut Criterion) {
    let g = Projective::generator();
    let h = Projective::generator_pedersen();
    c.bench_function("projective_add", |bencher| {
        bencher.iter(|| black_box(black_box(g) + black_box(h)))
    });
}

fn bench_affine_scalar_mul(c: &mut Criterion) {
    let g = Affine::generator();
    let mut rng = StdRng::seed_from_u64(42);
    let scalar = random_scalar(&mut rng);

    c.bench_function("affine_scalar_mul", |bencher| {
        bencher.iter(|| black_box(black_box(g).scalar_mul(black_box(&scalar))))
    });
}

fn bench_projective_scalar_mul(c: &mut Criterion) {
    let g = Projective::generator();
    let mut rng = StdRng::seed_from_u64(42);
    let scalar = random_scalar(&mut rng);

    c.bench_function("projective_scalar_mul", |bencher| {
        bencher.iter(|| black_box(black_box(g).scalar_mul(black_box(&scalar))))
    });
}

fn bench_affine_scalar_mul_windowed(c: &mut Criterion) {
    let g = Affine::generator();
    let mut rng = StdRng::seed_from_u64(42);
    let scalar = random_scalar(&mut rng);

    c.bench_function("affine_scalar_mul_windowed", |bencher| {
        bencher.iter(|| black_box(black_box(g).scalar_mul_windowed(black_box(&scalar))))
    });
}

fn bench_projective_scalar_mul_windowed(c: &mut Criterion) {
    let g = Projective::generator();
    let mut rng = StdRng::seed_from_u64(42);
    let scalar = random_scalar(&mut rng);

    c.bench_function("projective_scalar_mul_windowed", |bencher| {
        bencher.iter(|| black_box(black_box(g).scalar_mul_windowed(black_box(&scalar))))
    });
}

fn bench_affine_msm(c: &mut Criterion) {
    let mut group = c.benchmark_group("affine_msm");

    for size in [2, 4, 8, 16, 32, 64, 128].iter() {
        let mut rng = StdRng::seed_from_u64(12345);
        let g = Affine::generator();
        let h = Affine::generator_pedersen();

        // Create points by multiplying the generators
        let points: Vec<Affine> = (0..*size)
            .map(|i| {
                if i % 2 == 0 {
                    g.mul_u64(i as u64 + 1)
                } else {
                    h.mul_u64(i as u64 + 1)
                }
            })
            .collect();

        // Generate random scalars
        let scalars: Vec<ScalarField> = (0..*size).map(|_| random_scalar(&mut rng)).collect();

        group.bench_with_input(BenchmarkId::from_parameter(size), size, |bencher, _| {
            bencher.iter(|| {
                black_box(<Affine as Group>::multi_scalar_mul(
                    black_box(&points),
                    black_box(&scalars),
                ))
            })
        });
    }
    group.finish();
}

fn bench_projective_msm(c: &mut Criterion) {
    let mut group = c.benchmark_group("projective_msm");

    for size in [2, 4, 8, 16, 32, 64, 128].iter() {
        let mut rng = StdRng::seed_from_u64(12345);
        let g = Projective::generator();
        let h = Projective::generator_pedersen();

        // Create points by multiplying the generators
        let points: Vec<Projective> = (0..*size)
            .map(|i| {
                if i % 2 == 0 {
                    g.mul_u64(i as u64 + 1)
                } else {
                    h.mul_u64(i as u64 + 1)
                }
            })
            .collect();

        // Generate random scalars
        let scalars: Vec<ScalarField> = (0..*size).map(|_| random_scalar(&mut rng)).collect();

        group.bench_with_input(BenchmarkId::from_parameter(size), size, |bencher, _| {
            bencher.iter(|| {
                black_box(<Projective as Group>::multi_scalar_mul(
                    black_box(&points),
                    black_box(&scalars),
                ))
            })
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_affine_double,
    bench_projective_double,
    bench_affine_add,
    bench_projective_add,
    bench_affine_scalar_mul,
    bench_projective_scalar_mul,
    bench_affine_scalar_mul_windowed,
    bench_projective_scalar_mul_windowed,
    bench_affine_msm,
    bench_projective_msm
);
criterion_main!(benches);
