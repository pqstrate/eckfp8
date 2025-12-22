use super::*;
use p3_baby_bear::BabyBear;
use p3_field::PrimeCharacteristicRing;
use rand::SeedableRng;
use rand::rngs::StdRng;

#[test]
fn test_sign_verify() {
    let mut rng = StdRng::seed_from_u64(42);
    let sk = SigningKey::random(&mut rng);
    let pk = sk.verifying_key();
    let msg = [
        BabyBear::from_u32(1),
        BabyBear::from_u32(2),
        BabyBear::from_u32(3),
    ];

    let sig = sk.sign(&mut rng, &msg).expect("sign");
    let ok = pk.verify(&msg, &sig).expect("verify");
    assert!(ok);
}

#[test]
fn test_verify_rejects_wrong_message() {
    let mut rng = StdRng::seed_from_u64(42);
    let sk = SigningKey::random(&mut rng);
    let pk = sk.verifying_key();
    let msg = [
        BabyBear::from_u32(10),
        BabyBear::from_u32(11),
        BabyBear::from_u32(12),
    ];
    let sig = sk.sign(&mut rng, &msg).expect("sign");

    let wrong_msg = [
        BabyBear::from_u32(10),
        BabyBear::from_u32(11),
        BabyBear::from_u32(13),
    ];

    let ok = pk.verify(&wrong_msg, &sig).expect("verify");
    assert!(!ok);
}

#[test]
fn test_verify_rejects_wrong_key() {
    let mut rng = StdRng::seed_from_u64(42);
    let sk = SigningKey::random(&mut rng);
    let _pk = sk.verifying_key();
    let msg = [
        BabyBear::from_u32(21),
        BabyBear::from_u32(22),
        BabyBear::from_u32(23),
    ];
    let sig = sk.sign(&mut rng, &msg).expect("sign");

    let wrong_sk = SigningKey::random(&mut rng);
    let wrong_pk = wrong_sk.verifying_key();

    let ok = wrong_pk.verify(&msg, &sig).expect("verify");
    assert!(!ok);
}
