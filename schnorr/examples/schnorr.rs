use p3_baby_bear::BabyBear;
use p3_field::PrimeCharacteristicRing;
use rand::SeedableRng;
use rand::rngs::StdRng;
use schnorr::{SigningKey, VerifyingKey};

fn main() {
    let mut rng = StdRng::seed_from_u64(42);
    let sk = SigningKey::random(&mut rng);
    let vk = VerifyingKey::from(&sk);

    let sk_bytes = bincode::serialize(&sk).expect("serialize sk");
    let vk_bytes = bincode::serialize(&vk).expect("serialize vk");

    let msg_bytes = b"hello schnorr";
    let msg_field: Vec<BabyBear> = msg_bytes
        .iter()
        .copied()
        .map(|b| BabyBear::from_u32(b as u32))
        .collect();

    let sig = sk.sign(&mut rng, &msg_field).expect("sign");
    let sig_bytes = bincode::serialize(&sig).expect("serialize sig");

    let sk2: SigningKey = bincode::deserialize(&sk_bytes).expect("deserialize sk");
    let vk2: VerifyingKey = bincode::deserialize(&vk_bytes).expect("deserialize vk");
    let sig2 = bincode::deserialize(&sig_bytes).expect("deserialize sig");

    let ok = vk2.verify(&msg_field, &sig2).expect("verify");
    assert!(ok);

    let _ = sk2;
}
