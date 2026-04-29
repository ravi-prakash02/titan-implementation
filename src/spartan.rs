#[cfg(test)]
mod tests {
    use libspartan::{Instance, InputsAssignment, SNARKGens, SNARK, NIZKGens, NIZK, VarsAssignment};
    use merlin3::Transcript;
    use std::time::Instant;
    use serde::Serialize;

    #[test]
    fn test_snark_synthetic() {
        let num_vars = 1usize << 20;
        let num_cons = 1usize << 20;
        let num_inputs = 10;
        let num_non_zero_entries = 1usize << 20;

        let gens = SNARKGens::new(num_cons, num_vars, num_inputs, num_non_zero_entries);

        let (inst, vars, inputs) = Instance::produce_synthetic_r1cs(num_cons, num_vars, num_inputs);

        let start = Instant::now();
        let (comm, decomm) = SNARK::encode(&inst, &gens);
        println!("SNARK encode: {} ms", start.elapsed().as_millis());

        let mut prover_transcript = Transcript::new(b"snark_test");
        let start = Instant::now();
        let proof = SNARK::prove(&inst, &comm, &decomm, vars, &inputs, &gens, &mut prover_transcript);
        println!("SNARK prove: {} ms", start.elapsed().as_millis());
        let proof_bytes = bincode::serialize(&proof).unwrap();
        println!("SNARK proof size: {} bytes", proof_bytes.len());
        let mut verifier_transcript = Transcript::new(b"snark_test");
        let start = Instant::now();
        let res = proof.verify(&comm, &inputs, &mut verifier_transcript, &gens);
        println!("SNARK verify: {} ms", start.elapsed().as_millis());
        assert!(res.is_ok(), "SNARK proof verification failed");
    }

    #[test]
    fn test_nizk_synthetic() {
        let num_vars = 1024;
        let num_cons = 1024;
        let num_inputs = 10;

        let gens = NIZKGens::new(num_cons, num_vars, num_inputs);

        let (inst, vars, inputs) = Instance::produce_synthetic_r1cs(num_cons, num_vars, num_inputs);

        let mut prover_transcript = Transcript::new(b"nizk_test");
        let proof = NIZK::prove(&inst, vars, &inputs, &gens, &mut prover_transcript);

        let mut verifier_transcript = Transcript::new(b"nizk_test");
        assert!(
            proof.verify(&inst, &inputs, &mut verifier_transcript, &gens).is_ok(),
            "NIZK proof verification failed"
        );
    }
}
