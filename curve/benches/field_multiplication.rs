use criterion::{black_box, criterion_group, criterion_main, Criterion};
use curve::{BaseField, KoalaBear};

fn bench_koalabear_mul(c: &mut Criterion) {
    c.bench_function("koalabear_mul", |bencher| {
        let a = KoalaBear::new(123456789);
        let b = KoalaBear::new(987654321);
        bencher.iter(|| black_box(black_box(a) * black_box(b)))
    });
}

fn bench_koalabear_ext8_mul(c: &mut Criterion) {
    c.bench_function("koalabear_ext8_mul", |bencher| {
        // Create two extension field elements with non-trivial coefficients
        let a = BaseField::from([
            KoalaBear::new(1),
            KoalaBear::new(2),
            KoalaBear::new(3),
            KoalaBear::new(4),
            KoalaBear::new(5),
            KoalaBear::new(6),
            KoalaBear::new(7),
            KoalaBear::new(8),
        ]);
        let b = BaseField::from([
            KoalaBear::new(8),
            KoalaBear::new(7),
            KoalaBear::new(6),
            KoalaBear::new(5),
            KoalaBear::new(4),
            KoalaBear::new(3),
            KoalaBear::new(2),
            KoalaBear::new(1),
        ]);
        bencher.iter(|| black_box(black_box(a) * black_box(b)))
    });
}

criterion_group!(benches, bench_koalabear_mul, bench_koalabear_ext8_mul);
criterion_main!(benches);
