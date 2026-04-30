use crate::group_sumcheck::{
    compute_Sl_poly, compute_all_S_tables, compute_gi_values, eval_triple_at_alpha,
};
use crate::group_whir_committer::{
    BatchOpeningProof, GroupWhirCommitment, GroupWhirProverMerkleState,
    GroupWhirVerifierMerkleState, MerkleProofStrategy, OpeningProof,
};
use crate::merkle_tree::{Sha256Compress, Sha256LeafHash, Sha256MerkleTreeParams};
use crate::multilinear::MultilinearPoly;
use crate::titantranscript::{
    append_to_transcript, get_challenge, get_challenge_u64, transcript_append_point,
    transcript_append_scalar, transcript_challenge_scalar,
};
use crate::traits::{ByteSerializable, InnerProduct, Linear};
use ark_crypto_primitives::crh::pedersen::CRH;
use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::sponge::Absorb;
use ark_std::log2;
use merlin::Transcript;
use rand::prelude::StdRng;
use rand::{thread_rng, SeedableRng};
use std::fmt::Debug;
use std::marker::PhantomData;
use std::time::Instant;
use crate::utils::generate_power_vec;
//type Sha256Config = Sha256MerkleTreeParams;

pub trait SumCheckGroup {
    type Instance;
    type Proof;
    type Witness;
    type VerifierState;

    // prove and verify take reference to transcript to allow composition of protocols
    fn prove(
        instance: &Self::Instance,
        witness: &Self::Witness,
        transcript: &mut Transcript,
    ) -> Self::Proof;
    fn verify(
        proof: &Self::Proof,
        instance: &Self::Instance,
        transcript: &mut Transcript,
        state: &Self::VerifierState,
    ) -> bool;
}

// this encapsulates the protocol functions
pub struct ProtoSumCheckGroup<C, G, F> {
    phantom_data: PhantomData<(C, G, F)>,
}

pub struct SumCheckGroupInstance<C, G, F>
where
    C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>>,
{
    pub m: usize,                   // size of original polynomial
    pub l: usize,                   // folding parameter
    pub domain_g: Vec<F>,           // domain for folded function
    pub commitment: C::InnerDigest, // commitment to original function
    pub alpha: Vec<F>,              // evaluation point
    pub sigma: G,                   // claimed evaluation
    pub num_queries: usize,         // number of queries to oracles
}

pub struct SumCheckGroupWitness<C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>>, G, F>
{
    pub whir_commitment: GroupWhirCommitment<C, G>,
    pub phantom_data: PhantomData<(G, F)>,
}

pub struct SumCheckGroupProof<C, G, F>
where
    C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>>,
{
    pub round_messages: Vec<[G; 3]>,             // sum-check round messages
    pub g_poly: MultilinearPoly<G>,              // reduced (folded) polynomial
    pub opening_proofs: BatchOpeningProof<C, G>, // merkle opening proofs
    pub evaluations: Vec<G>,
    phantom_data: PhantomData<F>,
}

impl<C, G, F> SumCheckGroupProof<C, G, F>
where
    C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>>,
{
    /// Number of internal Merkle tree nodes transmitted in the proof.
    pub fn merkle_transmitted_nodes(&self) -> usize {
        self.opening_proofs.multi_proof.transmitted_node_count()
    }

    pub fn get_proof_size(&self) -> usize {
        let mut proof_size = 0usize;
        proof_size += self.merkle_transmitted_nodes();
        proof_size += 3usize * self.round_messages.len();
        let total: usize = self.opening_proofs.leaves.iter().map(Vec::len).sum();
        proof_size += total;
        proof_size += self.g_poly.coeffs.len();
        proof_size
    }
}

// C captures crypto config including leaf hasher, node hasher etc.
// G denotes group which is linear over F
// F is a field.
impl<C, G, F> SumCheckGroup for ProtoSumCheckGroup<C, G, F>
where
    C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>> + Send + Sync,
    G: Linear<F> + InnerProduct<F, Output = G> + ByteSerializable + Send + Debug + PartialEq,
    G::Affine: InnerProduct<F, Output = G>,
    F: ff::Field + ByteSerializable,
{
    type Instance = SumCheckGroupInstance<C, G, F>;
    type Proof = SumCheckGroupProof<C, G, F>;
    type Witness = SumCheckGroupWitness<C, G, F>;
    type VerifierState = GroupWhirVerifierMerkleState<C, G>;

    // This will generate a non-interactive proof for proving correct folding of the function
    // instance.commitment denotes the input function oracle.
    fn prove(
        instance: &Self::Instance,
        witness: &Self::Witness,
        transcript: &mut Transcript,
    ) -> Self::Proof {
        assert_eq!(
            instance.m, witness.whir_commitment.poly.num_vars,
            "Number of variables in polynomial and instance must match"
        );
        assert!(
            instance.l < instance.m,
            "Folding must be less than size of polynomial"
        );
        println!("Folding factor = {}", instance.l);

        // append statement to transcript
        transcript.append_message(b"protocol", b"sumcheck");
        transcript.append_u64(b"instance.m", instance.m as u64);
        transcript.append_u64(b"instance.l", instance.l as u64);
        transcript.append_message(b"commitment", &instance.commitment.to_sponge_bytes_as_vec());
        for a in instance.alpha.iter() {
            append_to_transcript::<F>(transcript, b"alpha", a);
        }
        append_to_transcript::<G>(transcript, b"sigma", &instance.sigma);

        // Precompute
        let precompute_time = Instant::now();
        let S_l = compute_Sl_poly(
            instance.m,
            instance.l,
            &witness.whir_commitment.poly,
            &instance.alpha,
        );
        let S_tables = compute_all_S_tables(S_l.coeffs, instance.l);
        println!(
            "Precompute time: {:?}",
            precompute_time.elapsed().as_millis()
        );

        // Send messages as the first l messages of group-sumcheck
        let mut rho_prefix: Vec<F> = Vec::new();
        let mut r_vec: Vec<F> = Vec::with_capacity(instance.l);
        let mut gi_triples: Vec<[G; 3]> = Vec::with_capacity(instance.m);
        let u_values = [F::ZERO, F::ONE, F::ONE + F::ONE];

        // Sumcheck rounds
        for i in 1..=instance.l {
            let start_gi_time = Instant::now();
            let gi_01 = compute_gi_values(
                instance.m,
                i,
                instance.l,
                &instance.alpha,
                &rho_prefix,
                &S_tables,
                &u_values,
            );
            append_to_transcript(transcript, b"g_i(0)", &gi_01[0]);
            append_to_transcript(transcript, b"g_i(1)", &gi_01[1]);
            append_to_transcript(transcript, b"g_i(2)", &gi_01[2]);

            // Verifier challenge r_i
            let r_i = get_challenge::<F>(transcript, b"r_i");
            rho_prefix.push(r_i);
            r_vec.push(r_i);
            println!("prover r_{}: {:?}", i, r_i);
            gi_triples.push([gi_01[0], gi_01[1], gi_01[2]]);

            // Check correctness
            if i == 1 {
                assert_eq!(gi_triples[0][0] + gi_triples[0][1], instance.sigma);
            } else {
                assert_eq!(
                    gi_triples[i - 1][0] + gi_triples[i - 1][1],
                    eval_triple_at_alpha(&gi_triples[i - 2], r_vec[i - 2])
                );
            }
            println!("prover gi_{}: {:?}", i, start_gi_time.elapsed().as_millis());
        }

        // At this stage prover has sent messages g_1,\ldots,g_l to verifier, while challenges r_1,\ldots,r_l
        // have been sent by the verifier. The reduced sum-check claim is now:
        // sum(h(r_1,\ldots,r_l, b)) = g_l(r_l)
        // The prover sends the folded poly h_poly
        let restrict_poly_time = Instant::now();
        let h_poly = witness.whir_commitment.poly.restrict(&r_vec);
        println!(
            "Restrict time: {:?}",
            restrict_poly_time.elapsed().as_millis()
        );
        for i in 0..h_poly.coeffs.len() {
            append_to_transcript(transcript, b"h_poly", &h_poly.coeffs[i]);
        }

        // sample num_queries points from domain
        let y_vec: Vec<usize> = (0..instance.num_queries)
            .into_iter()
            .map(|_| {
                (get_challenge_u64(transcript, b"challenge_y") as usize) % instance.domain_g.len()
            })
            .collect();

        // generate a multi-proof for above queries.
        let prover_state = GroupWhirProverMerkleState::new(MerkleProofStrategy::Compressed);

        // Open at multiple points
        let mut openings = Vec::new();

        for leaf_idx in 0..instance.num_queries {
            openings.push((y_vec[leaf_idx], r_vec.clone()));
        }

        // Now returns (BatchOpeningProof, Vec<G>)
        let (opening_proof, evaluations) = witness
            .whir_commitment
            .open_at_field_points_batch(&openings, &prover_state)
            .expect("Failed to open");

        Self::Proof {
            round_messages: gi_triples,
            g_poly: h_poly,
            opening_proofs: opening_proof,
            evaluations,
            phantom_data: PhantomData,
        }
    }

    // fn verify(
    //     proof: &Self::Proof,
    //     instance: &Self::Instance,
    //     transcript: &mut Transcript,
    //     state: &Self::VerifierState,
    // ) -> bool {
    //     // Add instance to transcript
    //     // append statement to transcript
    //     transcript.append_message(b"protocol", b"sumcheck");
    //     transcript.append_u64(b"instance.m", instance.m as u64);
    //     transcript.append_u64(b"instance.l", instance.l as u64);
    //     transcript.append_message(b"commitment", &instance.commitment.to_sponge_bytes_as_vec());
    //     for a in instance.alpha.iter() {
    //         append_to_transcript::<F>(transcript, b"alpha", a);
    //     }
    //     append_to_transcript::<G>(transcript, b"sigma", &instance.sigma);

    //     // verify sum-check rounds
    //     let mut r_vec: Vec<F> = Vec::new();
    //     for i in 1..=instance.l {
    //         let gi_triples = &proof.round_messages[i - 1].clone();
    //         append_to_transcript(transcript, b"g_i(0)", &gi_triples[0]);
    //         append_to_transcript(transcript, b"g_i(1)", &gi_triples[1]);
    //         append_to_transcript(transcript, b"g_i(2)", &gi_triples[2]);

    //         // Verifier challenge r_i
    //         let r_i = get_challenge::<F>(transcript, b"r_i");
    //         r_vec.push(r_i);
    //         println!("verifier r_{}: {:?}", i, r_i);

    //         // Check correctness
    //         if i == 1 {
    //             assert_eq!(gi_triples[0] + gi_triples[1], instance.sigma);
    //         } else {
    //             assert_eq!(
    //                 gi_triples[0] + gi_triples[1],
    //                 eval_triple_at_alpha(&proof.round_messages[i - 2], r_vec[i - 2])
    //             );
    //         }
    //     }

    //     // Append final folded polynomial
    //     for i in 0..proof.g_poly.coeffs.len() {
    //         append_to_transcript(transcript, b"h_poly", &proof.g_poly.coeffs[i]);
    //     }

    //     // get opening challenges
    //     let y_vec: Vec<usize> = (0..instance.num_queries)
    //         .into_iter()
    //         .map(|_| {
    //             (get_challenge_u64(transcript, b"challenge_y") as usize) % instance.domain_g.len()
    //         })
    //         .collect();

    //     // verify merkle opening proofs
    //     let mut openings = Vec::new();
    //     for leaf_idx in 0..instance.num_queries {
    //         openings.push((y_vec[leaf_idx], r_vec.clone()));
    //     }

    //     let proof_refs: Vec<_> = proof.opening_proofs.iter().map(|(p, _)| p).collect();
    //     let alphas: Vec<_> = openings.iter().map(|(_, a)| a.clone()).collect();
    //     let evals: Vec<_> = proof.opening_proofs.iter().map(|(_, e)| *e).collect();

    //     let start = Instant::now();
    //     let is_valid = GroupWhirCommitment::<C, G>::verify_openings_batch(
    //         proof_refs.as_slice(),
    //         &alphas,
    //         &evals,
    //         &instance.commitment,
    //         state,
    //     )
    //     .expect("Batch verification failed");
    //     println!("Batch verification time {}", start.elapsed().as_millis());

    //     //todo - check that g_folded is consistent with the evaluations

    //     assert!(is_valid);
    //     is_valid
    // }
    fn verify(
        proof: &Self::Proof,
        instance: &Self::Instance,
        transcript: &mut Transcript,
        state: &Self::VerifierState,
    ) -> bool {
        // Append statement to transcript
        transcript.append_message(b"protocol", b"sumcheck");
        transcript.append_u64(b"instance.m", instance.m as u64);
        transcript.append_u64(b"instance.l", instance.l as u64);
        transcript.append_message(b"commitment", &instance.commitment.to_sponge_bytes_as_vec());
        for a in instance.alpha.iter() {
            append_to_transcript::<F>(transcript, b"alpha", a);
        }
        append_to_transcript::<G>(transcript, b"sigma", &instance.sigma);

        // Verify sum-check rounds
        let mut r_vec: Vec<F> = Vec::new();
        for i in 1..=instance.l {
            let gi_triples = &proof.round_messages[i - 1].clone();
            append_to_transcript(transcript, b"g_i(0)", &gi_triples[0]);
            append_to_transcript(transcript, b"g_i(1)", &gi_triples[1]);
            append_to_transcript(transcript, b"g_i(2)", &gi_triples[2]);

            // Verifier challenge r_i
            let r_i = get_challenge::<F>(transcript, b"r_i");
            r_vec.push(r_i);
            println!("verifier r_{}: {:?}", i, r_i);

            // Check correctness
            if i == 1 {
                assert_eq!(gi_triples[0] + gi_triples[1], instance.sigma);
            } else {
                assert_eq!(
                    gi_triples[0] + gi_triples[1],
                    eval_triple_at_alpha(&proof.round_messages[i - 2], r_vec[i - 2])
                );
            }
        }

        // Append final folded polynomial
        for i in 0..proof.g_poly.coeffs.len() {
            append_to_transcript(transcript, b"h_poly", &proof.g_poly.coeffs[i]);
        }

        // Get opening challenges
        let y_vec: Vec<usize> = (0..instance.num_queries)
            .into_iter()
            .map(|_| {
                (get_challenge_u64(transcript, b"challenge_y") as usize) % instance.domain_g.len()
            })
            .collect();

        // Recompute evaluations from leaves
        let start = Instant::now();
        let eq_vec = MultilinearPoly::init_with_eq(&r_vec).coeffs;
        // Verify merkle opening proofs
        let start = Instant::now();
        let is_valid = GroupWhirCommitment::<C, G>::verify_openings_batch(
            &proof.opening_proofs,
            &r_vec,
            &proof.evaluations,
            &instance.commitment,
            state,
        )
        .expect("Batch verification failed");
        println!("Batch verification time {}", start.elapsed().as_millis());

        // Batched check: g_poly consistency with evaluations.
        // Instead of num_queries individual MSMs, combine with random γ powers
        // and perform a single MSM: <g_poly, Σ_i γ^i · eq(pow_y_i)> == Σ_i γ^i · eval_i
        let leaf_indices = proof.opening_proofs.leaf_indices();
        let num_tail_vars = instance.m - instance.l;
        let poly_size = 1usize << num_tail_vars;
        let gamma = F::random(&mut thread_rng()); // @todo make gamma part of transcript

        let mut combined_scalars = vec![F::ZERO; poly_size];
        let mut combined_eval = G::zero();
        let mut gamma_power = F::ONE;

        //@todo must check leaf_indices occur in y_vec

        for i in 0..leaf_indices.len() {
            let pow_y_vec = generate_power_vec(instance.domain_g[leaf_indices[i]], num_tail_vars, false);
            let eq_i = MultilinearPoly::init_with_eq(&pow_y_vec).coeffs;
            for j in 0..poly_size {
                combined_scalars[j] = combined_scalars[j] + gamma_power * eq_i[j];
            }
            combined_eval = combined_eval + proof.evaluations[i] * gamma_power;
            gamma_power = gamma_power * gamma;
        }

        let lhs = G::inner_product_msm(&proof.g_poly.coeffs, &combined_scalars);
        assert_eq!(lhs, combined_eval, "Batched g_poly evaluation check failed");
        assert!(is_valid);
        is_valid
    }
}

// ============================================================================
// Batch sumcheck: k polynomials with independent commitments, shared eval point
// ============================================================================

pub struct SumCheckGroupBatchInstance<C, G, F>
where
    C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>>,
{
    pub m: usize,                         // size of each polynomial
    pub l: usize,                         // folding parameter
    pub domain_g: Vec<F>,                 // domain for folded function
    pub commitments: Vec<C::InnerDigest>, // k Merkle roots (one per polynomial)
    pub alpha: Vec<F>,                    // shared evaluation point
    pub sigmas: Vec<G>,                   // k claimed evaluations
    pub num_queries: usize,               // number of queries to oracles
}

pub struct SumCheckGroupBatchWitness<C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>>, G, F>
{
    pub whir_commitments: Vec<GroupWhirCommitment<C, G>>,
    pub phantom_data: PhantomData<(G, F)>,
}

pub struct SumCheckGroupBatchProof<C, G, F>
where
    C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>>,
{
    pub round_messages: Vec<[G; 3]>,                  // sumcheck messages for BATCHED polynomial
    pub g_poly: MultilinearPoly<G>,                   // single folded batched polynomial
    pub opening_proofs: Vec<BatchOpeningProof<C, G>>,  // one per oracle (k total)
    pub evaluations: Vec<Vec<G>>,                      // evaluations[i][j] = oracle i at query j
    pub phantom_data: PhantomData<F>,
}

impl<C, G, F> SumCheckGroupBatchProof<C, G, F>
where
    C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>>,
{
    /// Total number of internal Merkle tree nodes transmitted across all k oracles.
    pub fn merkle_transmitted_nodes(&self) -> usize {
        self.opening_proofs
            .iter()
            .map(|p| p.multi_proof.transmitted_node_count())
            .sum()
    }

    pub fn get_proof_size(&self) -> usize {
        let mut proof_size = self.merkle_transmitted_nodes();
        proof_size += 3*self.round_messages.len();
        proof_size += self.g_poly.coeffs.len();
        proof_size
    }
}

impl<C, G, F> ProtoSumCheckGroup<C, G, F>
where
    C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>> + Send + Sync,
    G: Linear<F> + InnerProduct<F, Output = G> + ByteSerializable + Send + Debug + PartialEq,
    G::Affine: InnerProduct<F, Output = G>,
    F: ff::Field + ByteSerializable,
{
    /// Prove batch: k polynomials each with their own commitment, shared evaluation point.
    /// Returns (proof, gamma) where gamma is the random linear combination challenge.
    pub fn prove_batch(
        instance: &SumCheckGroupBatchInstance<C, G, F>,
        witness: &SumCheckGroupBatchWitness<C, G, F>,
        transcript: &mut Transcript,
    ) -> (SumCheckGroupBatchProof<C, G, F>, F) {
        let k = instance.commitments.len();
        assert_eq!(k, instance.sigmas.len());
        assert_eq!(k, witness.whir_commitments.len());
        for w in &witness.whir_commitments {
            assert_eq!(instance.m, w.poly.num_vars);
        }
        assert!(instance.l < instance.m);

        // Transcript setup
        transcript.append_message(b"protocol", b"batch_sumcheck");
        transcript.append_u64(b"k", k as u64);
        transcript.append_u64(b"instance.m", instance.m as u64);
        transcript.append_u64(b"instance.l", instance.l as u64);
        for root in &instance.commitments {
            transcript.append_message(b"commitment", &root.to_sponge_bytes_as_vec());
        }
        for a in instance.alpha.iter() {
            append_to_transcript::<F>(transcript, b"alpha", a);
        }
        for sigma in &instance.sigmas {
            append_to_transcript::<G>(transcript, b"sigma", sigma);
        }

        // Derive gamma
        let gamma = get_challenge::<F>(transcript, b"gamma");

        // Compute batched sigma
        let mut sigma_batched = G::zero();
        let mut gamma_power = F::ONE;
        for sigma in &instance.sigmas {
            sigma_batched = sigma_batched + *sigma * gamma_power;
            gamma_power = gamma_power * gamma;
        }

        // Compute S_l tables per-polynomial, then combine with gamma
        let mut S_l_batched = vec![G::zero(); 1usize << instance.l];
        gamma_power = F::ONE;
        for i in 0..k {
            let S_l_i = compute_Sl_poly(
                instance.m,
                instance.l,
                &witness.whir_commitments[i].poly,
                &instance.alpha,
            );
            for j in 0..S_l_batched.len() {
                S_l_batched[j] = S_l_batched[j] + S_l_i.coeffs[j] * gamma_power;
            }
            gamma_power = gamma_power * gamma;
        }
        let S_tables = compute_all_S_tables(S_l_batched, instance.l);

        // Sumcheck rounds
        let mut rho_prefix: Vec<F> = Vec::new();
        let mut r_vec: Vec<F> = Vec::with_capacity(instance.l);
        let mut gi_triples: Vec<[G; 3]> = Vec::with_capacity(instance.l);
        let u_values = [F::ZERO, F::ONE, F::ONE + F::ONE];

        for i in 1..=instance.l {
            let gi_01 = compute_gi_values(
                instance.m,
                i,
                instance.l,
                &instance.alpha,
                &rho_prefix,
                &S_tables,
                &u_values,
            );
            append_to_transcript(transcript, b"g_i(0)", &gi_01[0]);
            append_to_transcript(transcript, b"g_i(1)", &gi_01[1]);
            append_to_transcript(transcript, b"g_i(2)", &gi_01[2]);

            let r_i = get_challenge::<F>(transcript, b"r_i");
            rho_prefix.push(r_i);
            r_vec.push(r_i);
            gi_triples.push([gi_01[0], gi_01[1], gi_01[2]]);

            // Check correctness
            if i == 1 {
                assert_eq!(gi_triples[0][0] + gi_triples[0][1], sigma_batched);
            } else {
                assert_eq!(
                    gi_triples[i - 1][0] + gi_triples[i - 1][1],
                    eval_triple_at_alpha(&gi_triples[i - 2], r_vec[i - 2])
                );
            }
        }

        // Compute folded h_poly: per-polynomial restrict then combine with gamma
        let num_tail_vars = instance.m - instance.l;
        let poly_size = 1usize << num_tail_vars;
        let mut h_coeffs = vec![G::zero(); poly_size];
        gamma_power = F::ONE;
        for i in 0..k {
            let h_i = witness.whir_commitments[i].poly.restrict(&r_vec);
            for j in 0..poly_size {
                h_coeffs[j] = h_coeffs[j] + h_i.coeffs[j] * gamma_power;
            }
            gamma_power = gamma_power * gamma;
        }
        let h_poly = MultilinearPoly::new(h_coeffs);

        // Append h_poly to transcript
        for coeff in &h_poly.coeffs {
            append_to_transcript(transcript, b"h_poly", coeff);
        }

        // Sample query points from transcript
        let y_vec: Vec<usize> = (0..instance.num_queries)
            .map(|_| {
                (get_challenge_u64(transcript, b"challenge_y") as usize) % instance.domain_g.len()
            })
            .collect();

        // Open all k oracles at same query points
        let prover_state = GroupWhirProverMerkleState::new(MerkleProofStrategy::Compressed);
        let mut openings = Vec::new();
        for leaf_idx in 0..instance.num_queries {
            openings.push((y_vec[leaf_idx], r_vec.clone()));
        }

        let mut all_opening_proofs = Vec::with_capacity(k);
        let mut all_evaluations = Vec::with_capacity(k);
        for i in 0..k {
            let (opening_proof, evaluations) = witness
                .whir_commitments[i]
                .open_at_field_points_batch(&openings, &prover_state)
                .expect("Failed to open");
            all_opening_proofs.push(opening_proof);
            all_evaluations.push(evaluations);
        }

        let proof = SumCheckGroupBatchProof {
            round_messages: gi_triples,
            g_poly: h_poly,
            opening_proofs: all_opening_proofs,
            evaluations: all_evaluations,
            phantom_data: PhantomData,
        };

        (proof, gamma)
    }

    /// Verify batch sumcheck proof.
    /// Returns (is_valid, gamma) where gamma is the random linear combination challenge.
    pub fn verify_batch(
        proof: &SumCheckGroupBatchProof<C, G, F>,
        instance: &SumCheckGroupBatchInstance<C, G, F>,
        transcript: &mut Transcript,
        state: &GroupWhirVerifierMerkleState<C, G>,
    ) -> (bool, F)
    where
        G: ByteSerializable,
    {
        let k = instance.commitments.len();

        // Replay transcript identically
        transcript.append_message(b"protocol", b"batch_sumcheck");
        transcript.append_u64(b"k", k as u64);
        transcript.append_u64(b"instance.m", instance.m as u64);
        transcript.append_u64(b"instance.l", instance.l as u64);
        for root in &instance.commitments {
            transcript.append_message(b"commitment", &root.to_sponge_bytes_as_vec());
        }
        for a in instance.alpha.iter() {
            append_to_transcript::<F>(transcript, b"alpha", a);
        }
        for sigma in &instance.sigmas {
            append_to_transcript::<G>(transcript, b"sigma", sigma);
        }

        // Derive gamma
        let gamma = get_challenge::<F>(transcript, b"gamma");

        // Compute batched sigma
        let mut sigma_batched = G::zero();
        let mut gamma_power = F::ONE;
        for sigma in &instance.sigmas {
            sigma_batched = sigma_batched + *sigma * gamma_power;
            gamma_power = gamma_power * gamma;
        }

        // Verify sumcheck rounds
        let mut r_vec: Vec<F> = Vec::new();
        for i in 1..=instance.l {
            let gi_triples = &proof.round_messages[i - 1];
            append_to_transcript(transcript, b"g_i(0)", &gi_triples[0]);
            append_to_transcript(transcript, b"g_i(1)", &gi_triples[1]);
            append_to_transcript(transcript, b"g_i(2)", &gi_triples[2]);

            let r_i = get_challenge::<F>(transcript, b"r_i");
            r_vec.push(r_i);

            if i == 1 {
                assert_eq!(gi_triples[0] + gi_triples[1], sigma_batched);
            } else {
                assert_eq!(
                    gi_triples[0] + gi_triples[1],
                    eval_triple_at_alpha(&proof.round_messages[i - 2], r_vec[i - 2])
                );
            }
        }

        // Append h_poly to transcript
        for coeff in &proof.g_poly.coeffs {
            append_to_transcript(transcript, b"h_poly", coeff);
        }

        // Get opening challenges
        let y_vec: Vec<usize> = (0..instance.num_queries)
            .map(|_| {
                (get_challenge_u64(transcript, b"challenge_y") as usize) % instance.domain_g.len()
            })
            .collect();

        // Verify k Merkle proofs
        for i in 0..k {
            let is_valid = GroupWhirCommitment::<C, G>::verify_openings_batch(
                &proof.opening_proofs[i],
                &r_vec,
                &proof.evaluations[i],
                &instance.commitments[i],
                state,
            )
            .expect("Batch verification failed");
            assert!(is_valid, "Merkle verification failed for oracle {}", i);
        }

        // Batched oracle consistency check
        let leaf_indices = proof.opening_proofs[0].leaf_indices();
        let num_tail_vars = instance.m - instance.l;
        let poly_size = 1usize << num_tail_vars;
        let delta = F::random(&mut thread_rng());

        let mut combined_scalars = vec![F::ZERO; poly_size];
        let mut combined_eval = G::zero();
        let mut delta_power = F::ONE;

        for j in 0..leaf_indices.len() {
            let pow_y_vec = generate_power_vec(instance.domain_g[leaf_indices[j]], num_tail_vars, false);
            let eq_j = MultilinearPoly::init_with_eq(&pow_y_vec).coeffs;

            // combined_eval_j = sum_i gamma^i * evaluations[i][j]
            let mut combined_eval_j = G::zero();
            gamma_power = F::ONE;
            for i in 0..k {
                combined_eval_j = combined_eval_j + proof.evaluations[i][j] * gamma_power;
                gamma_power = gamma_power * gamma;
            }

            for idx in 0..poly_size {
                combined_scalars[idx] = combined_scalars[idx] + delta_power * eq_j[idx];
            }
            combined_eval = combined_eval + combined_eval_j * delta_power;
            delta_power = delta_power * delta;
        }

        let lhs = G::inner_product_msm(&proof.g_poly.coeffs, &combined_scalars);
        assert_eq!(lhs, combined_eval, "Batched g_poly evaluation check failed");

        (true, gamma)
    }
}

mod tests {
    use crate::merkle_tree::{Sha256Compress, Sha256LeafHash};
    use crate::pastatypes::{Point as G, Scalar as F};
    use ff::Field;
    use pasta_curves::group::Group;
    use rand::prelude::StdRng;
    use rand::SeedableRng;

    use super::*;
    use crate::utils::create_smooth_domain;
    #[test]
    fn test_sumcheck_prover() {
        let mut rng = StdRng::from_entropy();
        // Very large polynomial: 2^10 = 1024 coefficients
        let m = 8;
        let coeffs: Vec<G> = (0..1 << m).map(|_| G::random(&mut rng)).collect();

        let poly = MultilinearPoly::new(coeffs);
        let l = 2;

        let domain_points = create_smooth_domain(9);
        let (leaf_hash_params, two_to_one_params) =
            crate::merkle_tree::default_config::<G, Sha256LeafHash<G>, Sha256Compress>(&mut rng);
        let verifier_state = GroupWhirVerifierMerkleState::<Sha256MerkleTreeParams<G>, G>::new(
            MerkleProofStrategy::Compressed,
            leaf_hash_params,
            two_to_one_params,
        );

        let start = Instant::now();
        let commitment = GroupWhirCommitment::<Sha256MerkleTreeParams<G>, G>::new(
            &poly,
            l,
            domain_points.1.clone(),
            &leaf_hash_params,
            &two_to_one_params,
        )
        .expect("Failed to create commitment");
        let root = commitment.root();
        let duration = start.elapsed();
        println!("time to compute commitment is {:?}", duration);

        let alpha: Vec<F> = (0..m).map(|_| F::random(&mut rng)).collect();

        let start = Instant::now();
        let sigma = poly.evaluate_msm(&alpha);
        println!(
            "time to evaluate polynomial {} msec",
            start.elapsed().as_millis()
        );

        let instance = SumCheckGroupInstance::<Sha256MerkleTreeParams<G>, G, F> {
            m: m,
            l: l,
            domain_g: domain_points.1,
            commitment: root.clone(),
            alpha: alpha,
            sigma: sigma,
            num_queries: 80,
        };

        let witness = SumCheckGroupWitness {
            whir_commitment: commitment,
            phantom_data: Default::default(),
        };

        let mut transcript = Transcript::new(b"sumcheck");
        let start = Instant::now();
        let proof = ProtoSumCheckGroup::prove(&instance, &witness, &mut transcript);
        println!(
            "Time to generate group sumcheck proof = {} msec",
            start.elapsed().as_millis()
        );

        let mut transcript = Transcript::new(b"sumcheck");

        let _res = ProtoSumCheckGroup::verify(&proof, &instance, &mut transcript, &verifier_state);
    }

    #[test]
    fn benchmark_compressed_vs_uncompressed() {
        let mut rng = StdRng::from_entropy();

        let m = 8;
        let coeffs: Vec<G> = (0..1 << m).map(|_| G::random(&mut rng)).collect();
        let poly = MultilinearPoly::new(coeffs);

        let k = 2;
        let domain_points = create_smooth_domain(9);

        let (leaf_hash_params, two_to_one_params) =
            crate::merkle_tree::default_config::<G, Sha256LeafHash<G>, Sha256Compress>(&mut rng);

        let commitment = GroupWhirCommitment::<Sha256MerkleTreeParams<G>, G>::new(
            &poly,
            k,
            domain_points.1.clone(),
            &leaf_hash_params,
            &two_to_one_params,
        )
        .expect("Failed to create commitment");

        let root = commitment.root();
        let alpha: Vec<F> = (0..k).map(|_| F::random(&mut rng)).collect();

        // Same openings for both
        let openings: Vec<_> = (0..80)
            .map(|i| (i % domain_points.1.len(), alpha.clone()))
            .collect();

        println!("\n========== COMPRESSED ==========");
        let prover_compressed = GroupWhirProverMerkleState::new(MerkleProofStrategy::Compressed);
        let verifier_compressed = GroupWhirVerifierMerkleState::<Sha256MerkleTreeParams<G>, G>::new(
            MerkleProofStrategy::Compressed,
            leaf_hash_params.clone(),
            two_to_one_params.clone(),
        );

        let start = Instant::now();
        let (batch_proof_compressed, evals_compressed) = commitment
            .open_at_field_points_batch(&openings, &prover_compressed)
            .expect("Failed");
        println!("Proof generation time: {:?}", start.elapsed());

        // Measure proof size
        let compressed_size = match &batch_proof_compressed.multi_proof {
            crate::group_whir_committer::MerkleMultiProof::Compressed(mp) => {
                println!("Number of leaf_indexes: {}", mp.leaf_indexes.len());
                println!(
                    "Auth path nodes: {}",
                    mp.auth_paths_suffixes
                        .iter()
                        .map(|p| p.len())
                        .sum::<usize>()
                );
                mp.leaf_indexes.len()
                    + mp.auth_paths_suffixes
                        .iter()
                        .map(|p| p.len())
                        .sum::<usize>()
            }
            _ => 0,
        };
        println!("Total proof elements: {}", compressed_size);

        let start = Instant::now();
        let is_valid = GroupWhirCommitment::<Sha256MerkleTreeParams<G>, G>::verify_openings_batch(
            &batch_proof_compressed,
            &alpha,
            &evals_compressed,
            &root,
            &verifier_compressed,
        )
        .expect("Verification failed");
        println!("Verification time: {:?}", start.elapsed());
        assert!(is_valid);

        println!("\n========== UNCOMPRESSED ==========");
        let prover_uncompressed =
            GroupWhirProverMerkleState::new(MerkleProofStrategy::Uncompressed);
        let verifier_uncompressed =
            GroupWhirVerifierMerkleState::<Sha256MerkleTreeParams<G>, G>::new(
                MerkleProofStrategy::Uncompressed,
                leaf_hash_params,
                two_to_one_params,
            );

        let start = Instant::now();
        let (batch_proof_uncompressed, evals_uncompressed) = commitment
            .open_at_field_points_batch(&openings, &prover_uncompressed)
            .expect("Failed");
        println!("Proof generation time: {:?}", start.elapsed());

        // Measure proof size
        let uncompressed_size = match &batch_proof_uncompressed.multi_proof {
            crate::group_whir_committer::MerkleMultiProof::Uncompressed(paths) => {
                let total_nodes: usize = paths.iter().map(|p| p.auth_path.len()).sum();
                println!("Number of paths: {}", paths.len());
                println!("Total auth path nodes: {}", total_nodes);
                paths.len() + total_nodes
            }
            _ => 0,
        };
        println!("Total proof elements: {}", uncompressed_size);

        let start = Instant::now();
        let is_valid = GroupWhirCommitment::<Sha256MerkleTreeParams<G>, G>::verify_openings_batch(
            &batch_proof_uncompressed,
            &alpha,
            &evals_uncompressed,
            &root,
            &verifier_uncompressed,
        )
        .expect("Verification failed");
        println!("Verification time: {:?}", start.elapsed());
        assert!(is_valid);

        println!("\n========== COMPARISON ==========");
        println!("Compressed proof elements: {}", compressed_size);
        println!("Uncompressed proof elements: {}", uncompressed_size);
        println!(
            "Size reduction: {:.1}%",
            (1.0 - compressed_size as f64 / uncompressed_size as f64) * 100.0
        );
    }
}
