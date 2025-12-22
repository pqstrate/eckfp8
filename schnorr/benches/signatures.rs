use criterion::{Criterion, black_box, criterion_group, criterion_main};
use p3_baby_bear::BabyBear;
use p3_field::PrimeCharacteristicRing;
use rand::SeedableRng;
use rand::rngs::StdRng;
use schnorr::{SigningKey, VerifyingKey};

fn bench_sign(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(42);
    let sk = SigningKey::random(&mut rng);
    let msg = [
        BabyBear::from_u32(1),
        BabyBear::from_u32(2),
        BabyBear::from_u32(3),
    ];

    c.bench_function("schnorr_sign", |bencher| {
        bencher.iter(|| {
            let sig = sk.sign(&mut rng, black_box(&msg)).expect("sign");
            black_box(sig);
        })
    });
}

fn bench_verify(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(42);
    let sk = SigningKey::random(&mut rng);
    let vk = VerifyingKey::from(&sk);
    let msg = [
        BabyBear::from_u32(1),
        BabyBear::from_u32(2),
        BabyBear::from_u32(3),
    ];
    let sig = sk.sign(&mut rng, &msg).expect("sign");

    c.bench_function("schnorr_verify", |bencher| {
        bencher.iter(|| {
            let ok = vk.verify(black_box(&msg), black_box(&sig)).expect("verify");
            black_box(ok);
        })
    });
}

criterion_group!(benches, bench_sign, bench_verify);
criterion_main!(benches);
