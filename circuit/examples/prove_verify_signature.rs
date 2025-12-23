//! Example proving and verifying the Schnorr AIR over a full trace.

use circuit::{build_schnorr_trace, CircuitPoint, KoalaBear, SchnorrAir, SignatureWitness};
use p3_baby_bear::BabyBear;
use p3_challenger::{HashChallenger, SerializingChallenger32};
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::PrimeCharacteristicRing;
use p3_fri::{FriParameters, TwoAdicFriPcs};
use p3_keccak::{Keccak256Hash, KeccakF};
use p3_matrix::Matrix;
use p3_merkle_tree::MerkleTreeMmcs;
use p3_symmetric::{CompressionFunctionFromHasher, PaddingFreeSponge, SerializingHasher};
use p3_uni_stark::{
    prove_with_preprocessed, setup_preprocessed, verify_with_preprocessed, StarkConfig,
};
use p3_util::log2_strict_usize;
use rand::rng;
use schnorr::SigningKey;
use std::time::Instant;

fn encode_point(point: &CircuitPoint) -> [KoalaBear; 16] {
    let mut out = [KoalaBear::ZERO; 16];
    out[..8].copy_from_slice(&point.x);
    out[8..16].copy_from_slice(&point.y);
    out
}

fn main() {
    let mut rng = rng();
    let signing_key = SigningKey::random(&mut rng);
    let verifying_key = signing_key.verifying_key();

    let message = vec![
        BabyBear::from_u32(0xBEEF),
        BabyBear::from_u32(0x1234),
        BabyBear::from_u32(0xCAFE),
    ];

    let signature = signing_key.sign(&mut rng, &message).expect("sign");
    let witness = SignatureWitness::new(&signature, &verifying_key, &message).expect("witness");

    let trace = build_schnorr_trace(&witness);
    let height = trace.trace.height();
    let air = SchnorrAir::new(height);
    let width = trace.trace.width();
    let gates = height * width;

    type Val = KoalaBear;
    type Challenge = BinomialExtensionField<Val, 4>;

    type ByteHash = Keccak256Hash;
    let byte_hash = ByteHash {};

    type U64Hash = PaddingFreeSponge<KeccakF, 25, 17, 4>;
    let u64_hash = U64Hash::new(KeccakF {});
    type FieldHash = SerializingHasher<U64Hash>;
    let field_hash = FieldHash::new(u64_hash);

    type MyCompress = CompressionFunctionFromHasher<U64Hash, 2, 4>;
    let compress = MyCompress::new(u64_hash);

    type ValMmcs = MerkleTreeMmcs<
        [Val; p3_keccak::VECTOR_LEN],
        [u64; p3_keccak::VECTOR_LEN],
        FieldHash,
        MyCompress,
        4,
    >;
    let val_mmcs = ValMmcs::new(field_hash, compress);

    type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());

    type Challenger = SerializingChallenger32<Val, HashChallenger<u8, ByteHash, 32>>;
    let challenger = Challenger::from_hasher(vec![], byte_hash);

    let fri_params = FriParameters {
        log_blowup: 3,
        log_final_poly_len: 0,
        num_queries: 2,
        commit_proof_of_work_bits: 1,
        query_proof_of_work_bits: 1,
        mmcs: challenge_mmcs,
    };
    type Dft = Radix2DitParallel<Val>;
    type Pcs = TwoAdicFriPcs<Val, Dft, ValMmcs, ChallengeMmcs>;
    let pcs = Pcs::new(Dft::default(), val_mmcs, fri_params);

    type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;
    let config = MyConfig::new(pcs, challenger);

    let mut public_values = Vec::with_capacity(32);
    public_values.extend_from_slice(&encode_point(&witness.public_key));
    public_values.extend_from_slice(&encode_point(&witness.r));

    println!("Trace rows: {}", height);
    println!("Trace columns: {}", width);
    println!("Gates (rows * columns): {}", gates);

    let prove_start = Instant::now();
    let degree_bits = log2_strict_usize(height);
    let (preprocessed_prover_data, preprocessed_vk) =
        setup_preprocessed::<MyConfig, _>(&config, &air, degree_bits).expect("preprocessed");
    let proof = prove_with_preprocessed(
        &config,
        &air,
        trace.trace,
        &public_values,
        Some(&preprocessed_prover_data),
    );
    let proof_bytes = postcard::to_allocvec(&proof).expect("serialize proof");
    let prove_time = prove_start.elapsed();

    let verify_start = Instant::now();
    verify_with_preprocessed(
        &config,
        &air,
        &proof,
        &public_values,
        Some(&preprocessed_vk),
    )
    .expect("verify");
    let verify_time = verify_start.elapsed();

    println!("Proving time: {:?}", prove_time);
    println!("Proof size (bytes): {}", proof_bytes.len());
    println!("Verification time: {:?}", verify_time);
    println!("Schnorr proof verified.");
}
