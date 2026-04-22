use ark_bls12_381::Fr;
use ark_serialize::CanonicalSerialize;
use ark_std::{test_rng, UniformRand};
use quarks_zk::{KopisPCS, PolynomialCommitmentScheme};
use std::time::Instant;

fn run_kopis_once(log_n: usize) -> (u128, u128, u128, usize) {
    let mut rng = test_rng();
    let n = 1 << log_n;

    println!("Benchmarking Kopis PCS for n = 2^{}", log_n);

    // setup
    let params = KopisPCS::setup(log_n, &mut rng);

    // Polynomial evaluations
    let evals: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
    let point: Vec<Fr> = (0..log_n).map(|_| Fr::rand(&mut rng)).collect();

    // commit
    let start_commit = Instant::now();
    let commitment = KopisPCS::commit(&params, &evals);
    let commit_ms = start_commit.elapsed().as_millis();

    // prove
    let start_prove = Instant::now();
    let (value, proof) =KopisPCS::prove_eval(&params, &evals, &point, &mut rng);
    let prove_ms = start_prove.elapsed().as_millis() - commit_ms;

    // verify
    let start_verify = Instant::now();
    let ok = KopisPCS::verify_eval(
        &params,
        &commitment,
        &point,
        value,
        &proof,
    );
    let verify_ms = start_verify.elapsed().as_millis();
    assert!(ok, "Kopis PCS proof should verify");

    // proof size
    // let gt_elems = log_n;
    // let gt_bytes = 128; // BLS12-381 GT
    // let proof_size = gt_elems * gt_bytes;
    // //let total_prover_ms = commit_ms + prove_ms;

    let mut proof_bytes = Vec::new();
    proof.serialize_compressed(&mut proof_bytes).unwrap();
    let proof_size = proof_bytes.len();

    println!(
        "commit = {} ms, prove = {} ms, verifier = {} ms, proof = {} bytes",
        commit_ms, prove_ms, verify_ms, proof_size
    );

    (commit_ms, prove_ms, verify_ms, proof_size)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore] 
    fn kopis_scaling_experiment() {
        println!("log_n,commit_ms,prove_ms,verifier_ms,proof_bytes");

        for log_n in (18..=26).step_by(2) {
            let runs = 5;
            let mut c_sum = 0;
            let mut p_sum = 0;
            let mut v_sum = 0;
            let mut proof_size = 0;

            for _ in 0..runs {
                let (c, p, v, size) = run_kopis_once(log_n);
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