use ark_bls12_381::Bls12_381;
use dory_pcs::backends::arkworks::{
    ArkworksPolynomial, Blake2bTranscript, G1Routines, G2Routines, BN254,
};
use dory_pcs::backends::{ArkFr, ArkG1, ArkG2, ArkGT};
use dory_pcs::primitives::arithmetic::DoryRoutines;
use dory_pcs::{prove, setup, verify};
use std::time::Instant;
use ark_bls12_381::Fr;
use ark_std::{test_rng, UniformRand};
use dory_pcs::primitives::arithmetic::Field;
use dory_pcs::Polynomial;
use quarks_zk::{DoryPCS, PolynomialCommitmentScheme};
use ark_serialize::CanonicalSerialize;

fn run_dory_BN254_once(log_n: usize) -> (u128, u128, usize) {
        let mut rng = rand::thread_rng();
        let n = 1 << log_n;
        let nu = log_n / 2; // Rows log2
        let sigma = log_n - nu; // Cols log2

        println!("Benchmarking Dory PCS for n = 2^{}", log_n);

        // --- SETUP ---
        let (prover_setup, verifier_setup) = setup::<BN254, _>(&mut rng, log_n);

        // Create a random polynomial and evaluation point
        let coefficients: Vec<_> = (0..n).map(|_| ArkFr::random(&mut rng)).collect();
        let polynomial = ArkworksPolynomial::new(coefficients);
        let point: Vec<_> = (0..log_n).map(|_| ArkFr::random(&mut rng)).collect();

        // --- 1. COMMITMENT TIME ---
        let start_commit = Instant::now();
        let (tier_2_comm, row_commitments) = polynomial
            .commit::<BN254, G1Routines>(nu, sigma, &prover_setup)
            .unwrap();
        let duration_commit_ms = start_commit.elapsed().as_millis();
        println!("Commitment Time:      {:?}", duration_commit_ms);

        // --- 2. EVALUATION PROOF TIME ---
        let mut prover_transcript = Blake2bTranscript::new(b"dory-benchmark");
        let start_prove = Instant::now();
        let proof = prove::<_, BN254, G1Routines, G2Routines, _, _>(
            &polynomial,
            &point,
            row_commitments,
            nu,
            sigma,
            &prover_setup,
            &mut prover_transcript,
        )
        .unwrap();
        let duration_prove_ms = start_prove.elapsed().as_millis();
        println!("Evaluation Proof Time: {:?}", duration_prove_ms);

        // --- PROOF SIZE ---
        let group_elems = 2 * log_n;   // G1 + G2 elements
        let field_elems = 2 * log_n;   // scalar challenges
        let proof_elems = group_elems + field_elems;
        let proof_size = 32 * proof_elems;

        // --- 3. VERIFIER TIME ---
        let evaluation = polynomial.evaluate(&point);
        let mut verifier_transcript = Blake2bTranscript::new(b"dory-benchmark");

        let start_verify = Instant::now();
        verify::<_, BN254, G1Routines, G2Routines, _>(
            tier_2_comm,
            evaluation,
            &point,
            &proof,
            verifier_setup,
            &mut verifier_transcript,
        )
        .unwrap();
        let duration_verify_ms = start_verify.elapsed().as_millis();
        println!("Verifier Time:        {:?}", duration_verify_ms);
        

        let total_prover_ms = duration_commit_ms + duration_prove_ms;

        (total_prover_ms, duration_verify_ms, proof_size)

    }

    fn run_dory_BLS_once(log_n: usize) -> (u128, u128, u128, usize) {
        let mut rng = test_rng();
        
        let n: usize = 1 << log_n;
        let params = DoryPCS::setup(log_n, &mut rng);

        let evals: Vec<Fr> = (0..n).map(|i| Fr::rand(&mut rng)).collect();
        let start_commit = Instant::now();
        let commitment = DoryPCS::commit(&params, &evals);
        let duration_commit_ms = start_commit.elapsed().as_millis();
        println!("Commitment Time:      {:?}", duration_commit_ms);

        // Random evaluation point
        let point: Vec<Fr> = (0..log_n).map(|_| Fr::rand(&mut rng)).collect();

        // Generate proof
        let start_prove = Instant::now();
        let (value, proof) = DoryPCS::prove_eval(&params, &evals, &point, &mut rng);
        let duration_prove_ms = start_prove.elapsed().as_millis() - duration_commit_ms; //prove_eval recomputes the commitment, so taking the difference
        println!("Evaluation Proof Time: {:?}", duration_prove_ms);

        // proof size (BLS12-381)
        // let gt_elems = 6 * log_n;
        // let g1_elems = 3 * log_n;
        // let g2_elems = 3 * log_n;

        // // group sizes
        // let gt_bytes = 128;
        // let g1_bytes = 48;
        // let g2_bytes = 96;

        // let proof_size =
        //     gt_elems * gt_bytes +
        //     g1_elems * g1_bytes +
        //     g2_elems * g2_bytes;

        let mut proof_bytes = Vec::new();
        proof.serialize_compressed(&mut proof_bytes).unwrap();
        let proof_size = proof_bytes.len();

        // Verify
        let start_verify = Instant::now();
        let valid = DoryPCS::verify_eval(&params, &commitment, &point, value, &proof);
        let duration_verify_ms = start_verify.elapsed().as_millis();
        println!("Verifier Time:        {:?}", duration_verify_ms);
        
        //let total_prover_ms = duration_commit_ms + duration_prove_ms;

        (duration_commit_ms, duration_prove_ms, duration_verify_ms, proof_size)
    }

mod tests {
    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_bls12_381::Fr;
    use ark_std::{test_rng, UniformRand};
    use dory_pcs::backends::ArkFr;
    use dory_pcs::primitives::arithmetic::Field;
    use dory_pcs::Polynomial;
    use quarks_zk::{DoryPCS, PolynomialCommitmentScheme};
    #[test]
    fn benchmark_dory() {
        let mut rng = rand::thread_rng();

        // Configuration: 2^10 coefficients
        let log_n = 20;
        let n = 1 << log_n;
        let nu = log_n / 2; // Rows log2
        let sigma = log_n - nu; // Cols log2

        println!("Benchmarking Dory PCS for n = 2^{}", log_n);

        // --- SETUP ---
        let (prover_setup, verifier_setup) = setup::<BN254, _>(&mut rng, log_n);

        // Create a random polynomial and evaluation point
        let coefficients: Vec<_> = (0..n).map(|_| ArkFr::random(&mut rng)).collect();
        let polynomial = ArkworksPolynomial::new(coefficients);
        let point: Vec<_> = (0..log_n).map(|_| ArkFr::random(&mut rng)).collect();

        // --- 1. COMMITMENT TIME ---
        let start_commit = Instant::now();
        let (tier_2_comm, row_commitments) = polynomial
            .commit::<BN254, G1Routines>(nu, sigma, &prover_setup)
            .unwrap();
        let duration_commit = start_commit.elapsed();
        println!("Commitment Time:      {:?}", duration_commit);

        // --- 2. EVALUATION PROOF TIME ---
        let mut prover_transcript = Blake2bTranscript::new(b"dory-benchmark");
        let start_prove = Instant::now();
        let proof = prove::<_, BN254, G1Routines, G2Routines, _, _>(
            &polynomial,
            &point,
            row_commitments,
            nu,
            sigma,
            &prover_setup,
            &mut prover_transcript,
        )
        .unwrap();
        let duration_prove = start_prove.elapsed();
        println!("Evaluation Proof Time: {:?}", duration_prove);

        // --- 3. VERIFIER TIME ---
        let evaluation = polynomial.evaluate(&point);
        let mut verifier_transcript = Blake2bTranscript::new(b"dory-benchmark");

        let start_verify = Instant::now();
        verify::<_, BN254, G1Routines, G2Routines, _>(
            tier_2_comm,
            evaluation,
            &point,
            &proof,
            verifier_setup,
            &mut verifier_transcript,
        )
        .unwrap();
        let duration_verify = start_verify.elapsed();
        println!("Verifier Time:        {:?}", duration_verify);
    }

    #[test]
    fn test_dory_pcs_prove_verify() {
        let mut rng = test_rng();
        let log_n = 20usize;
        let n: usize = 1 << log_n;
        let params = DoryPCS::setup(log_n, &mut rng);

        let evals: Vec<Fr> = (0..n).map(|i| Fr::rand(&mut rng)).collect();
        let start = Instant::now();
        let commitment = DoryPCS::commit(&params, &evals);
        println!("Time to commit polynomial = {} msec", start.elapsed().as_millis());

        // Random evaluation point
        let point: Vec<Fr> = (0..log_n).map(|_| Fr::rand(&mut rng)).collect();

        // Generate proof
        let mut start = Instant::now();
        let (value, proof) = DoryPCS::prove_eval(&params, &evals, &point, &mut rng);
        println!("Time to generate proof: {:?}", start.elapsed().as_millis());

        // Verify
        start = Instant::now();
        let valid = DoryPCS::verify_eval(&params, &commitment, &point, value, &proof);
        println!("Time to verify proof: {:?}", start.elapsed().as_millis());
        assert!(valid, "Dory-PCS proof should verify");
    }

    #[test]
    #[ignore]
    fn dory_BN254_scaling_experiment() {
        println!("m,commit_ms,prove_ms,verifier_ms,proof_bytes");

        for log_n in (18..=26).step_by(2) {
            let runs = 10;
            let mut p_sum = 0;
            let mut v_sum = 0;
            let mut proof_size = 0;

            for _ in 0..runs {
                let (p, v, size) = run_dory_BN254_once(log_n);
                p_sum += p;
                v_sum += v;
                proof_size = size;
            }

            println!(
                "{},{},{},{}",
                log_n,
                p_sum / runs,
                v_sum / runs,
                proof_size
            );
        }
    }

    #[test]
    #[ignore]
    fn dory_BLS_scaling_experiment() {
        println!("m,commit_ms,prove_ms,verifier_ms,proof_bytes");

        for log_n in (18..=26).step_by(2) {
            let runs = 1;
            let mut c_sum = 0;
            let mut p_sum = 0;
            let mut v_sum = 0;
            let mut proof_size = 0;

            for _ in 0..runs {
                let (c, p, v, size) = run_dory_BLS_once(log_n);
                c_sum += c;
                p_sum += p;
                v_sum += v;
                proof_size += size;
            }

            println!(
                "{},{},{},{},{}",
                log_n,
                c_sum / runs,
                p_sum / runs,
                v_sum / runs,
                proof_size / runs as usize,
            );
        }
    }
}
