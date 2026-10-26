use crate::group_sumcheck::{
    compute_Sl_poly, compute_all_S_tables, compute_gi_values, eval_triple_at_alpha,
};
use crate::group_whir_committer::{
    ood_point, BatchOpeningProof, GroupWhirCommitment, GroupWhirProverMerkleState,
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
    pub ood_eval: G,                // out-of-domain evaluation y_0 bound in the commitment (Sec. 5.4)
    pub num_queries: usize,         // number of queries to oracles
}

/// Appends the statement, including the out-of-domain claim, to the transcript and returns the
/// challenge xi batching the claim at alpha with the claim at w_0 = ood_point(root) (Sec. 5.4).
fn absorb_statement<C, G, F>(instance: &SumCheckGroupInstance<C, G, F>, transcript: &mut Transcript) -> F
where
    C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>>,
    G: ByteSerializable,
    F: ByteSerializable,
{
    transcript.append_message(b"protocol", b"sumcheck");
    transcript.append_u64(b"instance.m", instance.m as u64);
    transcript.append_u64(b"instance.l", instance.l as u64);
    transcript.append_message(b"commitment", &instance.commitment.to_sponge_bytes_as_vec());
    for a in instance.alpha.iter() {
        append_to_transcript::<F>(transcript, b"alpha", a);
    }
    append_to_transcript::<G>(transcript, b"sigma", &instance.sigma);
    append_to_transcript::<G>(transcript, b"ood_eval", &instance.ood_eval);
    get_challenge::<F>(transcript, b"xi")
}

/// Round message for the batched claim sigma + xi * y_0: g_i = g_i^{alpha} + xi * g_i^{w_0}.
pub(crate) fn compute_batched_gi_values<C, F>(
    m: usize,
    i: usize,
    l: usize,
    alpha: &[F],
    S_tables_alpha: &Vec<Vec<C>>,
    w0: &[F],
    S_tables_w0: &Vec<Vec<C>>,
    xi: F,
    rho_prefix: &[F],
    u_values: &[F],
) -> Vec<C>
where
    C: Linear<F> + InnerProduct<F, Output = C>,
    F: ff::Field,
{
    let g_alpha = compute_gi_values(m, i, l, alpha, rho_prefix, S_tables_alpha, u_values);
    let g_w0 = compute_gi_values(m, i, l, w0, rho_prefix, S_tables_w0, u_values);
    g_alpha.into_iter().zip(g_w0).map(|(a, b)| a + b * xi).collect()
}

/// Scalars e with <h, e> = sum_b h(b) * (eq(alpha, (r, b)) + xi * eq(w_0, (r, b))), i.e. the
/// value g_l(r_l) must take for the reduced claim (Figure 2, Step 6, with the OOD claim batched in).
fn terminal_scalars<F: ff::Field>(alpha: &[F], w0: &[F], xi: F, r_vec: &[F]) -> Vec<F> {
    let l = r_vec.len();
    let eq_prefix = |p: &[F]| {
        p.iter().zip(r_vec).fold(F::ONE, |acc, (a, r)| acc * (*a * *r + (F::ONE - *a) * (F::ONE - *r)))
    };
    let (c_alpha, c_w0) = (eq_prefix(&alpha[..l]), eq_prefix(&w0[..l]) * xi);
    let eq_alpha = MultilinearPoly::init_with_eq(&alpha[l..]).coeffs;
    let eq_w0 = MultilinearPoly::init_with_eq(&w0[l..]).coeffs;
    eq_alpha.iter().zip(eq_w0.iter()).map(|(a, b)| *a * c_alpha + *b * c_w0).collect()
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

        // append statement to transcript; xi batches the claim at alpha with the OOD claim
        let xi = absorb_statement(instance, transcript);
        let w0 = ood_point::<C, F>(&instance.commitment, instance.m);
        let sigma_batched = instance.sigma + instance.ood_eval * xi;

        // Precompute
        let precompute_time = Instant::now();
        let S_l = compute_Sl_poly(
            instance.m,
            instance.l,
            &witness.whir_commitment.poly,
            &instance.alpha,
        );
        let S_tables = compute_all_S_tables(S_l.coeffs, instance.l);
        let S_l_w0 = compute_Sl_poly(instance.m, instance.l, &witness.whir_commitment.poly, &w0);
        let S_tables_w0 = compute_all_S_tables(S_l_w0.coeffs, instance.l);
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
            let gi_01 = compute_batched_gi_values(
                instance.m,
                i,
                instance.l,
                &instance.alpha,
                &S_tables,
                &w0,
                &S_tables_w0,
                xi,
                &rho_prefix,
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
                assert_eq!(gi_triples[0][0] + gi_triples[0][1], sigma_batched);
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
        let num_tail_vars = instance.m - instance.l;
        let poly_size = 1usize << num_tail_vars;

        // Reject malformed instances/proofs before touching the transcript
        if instance.l == 0 || instance.l >= instance.m || instance.alpha.len() != instance.m {
            println!("Reject: malformed sum-check instance");
            return false;
        }
        if proof.round_messages.len() != instance.l {
            println!("Reject: expected {} round messages, got {}", instance.l, proof.round_messages.len());
            return false;
        }
        if proof.g_poly.coeffs.len() != poly_size {
            println!("Reject: folded polynomial has wrong size");
            return false;
        }

        // Append statement to transcript; batch the claim at alpha with the OOD claim (Sec. 5.4)
        let xi = absorb_statement(instance, transcript);
        let w0 = ood_point::<C, F>(&instance.commitment, instance.m);

        // Verify sum-check rounds
        let mut r_vec: Vec<F> = Vec::new();
        for i in 1..=instance.l {
            let gi_triples = &proof.round_messages[i - 1];
            append_to_transcript(transcript, b"g_i(0)", &gi_triples[0]);
            append_to_transcript(transcript, b"g_i(1)", &gi_triples[1]);
            append_to_transcript(transcript, b"g_i(2)", &gi_triples[2]);

            // Verifier challenge r_i
            let r_i = get_challenge::<F>(transcript, b"r_i");
            r_vec.push(r_i);
            println!("verifier r_{}: {:?}", i, r_i);

            // Check correctness
            let expected = if i == 1 {
                instance.sigma + instance.ood_eval * xi
            } else {
                eval_triple_at_alpha(&proof.round_messages[i - 2], r_vec[i - 2])
            };
            if gi_triples[0] + gi_triples[1] != expected {
                println!("Reject: sum-check round {} check failed", i);
                return false;
            }
        }

        // Reduced sum-check terminal check (Figure 2, Step 6), with the OOD claim batched in:
        // g_l(r_l) = eq(alpha_{..l}, r) * h(alpha_{l..}) + xi * eq(w0_{..l}, r) * h(w0_{l..})
        let g_l_at_r_l = eval_triple_at_alpha(&proof.round_messages[instance.l - 1], r_vec[instance.l - 1]);
        let scalars = terminal_scalars(&instance.alpha, &w0, xi, &r_vec);
        if g_l_at_r_l != G::inner_product_msm(&proof.g_poly.coeffs, &scalars) {
            println!("Reject: reduced sum-check terminal check failed");
            return false;
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

        // Bind the opened leaves to the transcript queries y_vec
        if !proof.opening_proofs.opens_exactly(&y_vec) {
            println!("Reject: opened leaves do not match queried indices");
            return false;
        }
        let leaf_indices = proof.opening_proofs.leaf_indices();
        let coset_size = 1usize << instance.l;
        if proof.evaluations.len() != leaf_indices.len()
            || proof.opening_proofs.leaves.iter().any(|leaf| leaf.len() != coset_size)
        {
            println!("Reject: malformed opening proof");
            return false;
        }

        // Verify merkle opening proofs
        let start = Instant::now();
        let is_valid = GroupWhirCommitment::<C, G>::verify_openings_batch(
            &proof.opening_proofs,
            &r_vec,
            &proof.evaluations,
            &instance.commitment,
            state,
        )
        .unwrap_or(false);
        println!("Batch verification time {}", start.elapsed().as_millis());
        if !is_valid {
            println!("Reject: Merkle opening verification failed");
            return false;
        }

        // Batched check: g_poly consistency with evaluations.
        // Instead of num_queries individual MSMs, combine with random γ powers
        // and perform a single MSM: <g_poly, Σ_i γ^i · eq(pow_y_i)> == Σ_i γ^i · eval_i
        let gamma = F::random(&mut thread_rng()); // private verifier coin

        let mut combined_scalars = vec![F::ZERO; poly_size];
        let mut combined_eval = G::zero();
        let mut gamma_power = F::ONE;

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
        if lhs != combined_eval {
            println!("Reject: batched g_poly evaluation check failed");
            return false;
        }
        true
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
        let num_tail_vars = instance.m - instance.l;
        let poly_size = 1usize << num_tail_vars;

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

        // Reject malformed proofs
        if k == 0
            || instance.l == 0
            || instance.l >= instance.m
            || instance.alpha.len() != instance.m
            || instance.sigmas.len() != k
            || proof.round_messages.len() != instance.l
            || proof.g_poly.coeffs.len() != poly_size
            || proof.opening_proofs.len() != k
            || proof.evaluations.len() != k
        {
            println!("Reject: malformed batch sum-check proof");
            return (false, gamma);
        }

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

            let expected = if i == 1 {
                sigma_batched
            } else {
                eval_triple_at_alpha(&proof.round_messages[i - 2], r_vec[i - 2])
            };
            if gi_triples[0] + gi_triples[1] != expected {
                println!("Reject: batch sum-check round {} check failed", i);
                return (false, gamma);
            }
        }

        // Reduced sum-check terminal check (Figure 2, Step 6) on the batched folded polynomial
        let mut eq_alpha_r = F::ONE;
        for j in 0..instance.l {
            let (a, r) = (instance.alpha[j], r_vec[j]);
            eq_alpha_r = eq_alpha_r * (a * r + (F::ONE - a) * (F::ONE - r));
        }
        let g_l_at_r_l = eval_triple_at_alpha(&proof.round_messages[instance.l - 1], r_vec[instance.l - 1]);
        if g_l_at_r_l != proof.g_poly.evaluate_msm(&instance.alpha[instance.l..]) * eq_alpha_r {
            println!("Reject: batch reduced sum-check terminal check failed");
            return (false, gamma);
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

        // Every oracle must open exactly the queried leaves y_vec, in the same order, so that
        // evaluations[i][j] for different i refer to the same domain point.
        let leaf_indices = proof.opening_proofs[0].leaf_indices();
        let coset_size = 1usize << instance.l;
        for i in 0..k {
            let op = &proof.opening_proofs[i];
            if !op.opens_exactly(&y_vec)
                || op.leaf_indices() != leaf_indices
                || proof.evaluations[i].len() != leaf_indices.len()
                || op.leaves.iter().any(|leaf| leaf.len() != coset_size)
            {
                println!("Reject: oracle {} openings do not match queried indices", i);
                return (false, gamma);
            }
        }

        // Verify k Merkle proofs
        for i in 0..k {
            let is_valid = GroupWhirCommitment::<C, G>::verify_openings_batch(
                &proof.opening_proofs[i],
                &r_vec,
                &proof.evaluations[i],
                &instance.commitments[i],
                state,
            )
            .unwrap_or(false);
            if !is_valid {
                println!("Reject: Merkle verification failed for oracle {}", i);
                return (false, gamma);
            }
        }

        // Batched oracle consistency check
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
        if lhs != combined_eval {
            println!("Reject: batched g_poly evaluation check failed");
            return (false, gamma);
        }

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
    use std::collections::BTreeSet;
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
            ood_eval: commitment.ood_eval,
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

    // ------------------------------------------------------------------------
    // Soundness / rejection tests
    // ------------------------------------------------------------------------

    type Cfg = Sha256MerkleTreeParams<G>;

    struct Fixture {
        polys: Vec<MultilinearPoly<G>>,
        commitments: Vec<GroupWhirCommitment<Cfg, G>>,
        domain: Vec<F>,
        alpha: Vec<F>,
        verifier_state: GroupWhirVerifierMerkleState<Cfg, G>,
    }

    const M: usize = 6;
    const L: usize = 2;
    const NUM_QUERIES: usize = 16;

    fn fixture(k: usize) -> Fixture {
        let mut rng = StdRng::seed_from_u64(7);
        let domain = create_smooth_domain(M - L + 3).1;
        let (leaf_hash_params, two_to_one_params) =
            crate::merkle_tree::default_config::<G, Sha256LeafHash<G>, Sha256Compress>(&mut rng);
        let polys: Vec<_> = (0..k)
            .map(|_| MultilinearPoly::new((0..1 << M).map(|_| G::random(&mut rng)).collect()))
            .collect();
        let commitments = polys
            .iter()
            .map(|p| {
                GroupWhirCommitment::<Cfg, G>::new(p, L, domain.clone(), &leaf_hash_params, &two_to_one_params)
                    .unwrap()
            })
            .collect();
        let alpha = (0..M).map(|_| F::random(&mut rng)).collect();
        let verifier_state = GroupWhirVerifierMerkleState::new(
            MerkleProofStrategy::Compressed,
            leaf_hash_params,
            two_to_one_params,
        );
        Fixture { polys, commitments, domain, alpha, verifier_state }
    }

    fn instance(fx: &Fixture, sigma: G) -> SumCheckGroupInstance<Cfg, G, F> {
        SumCheckGroupInstance {
            m: M,
            l: L,
            domain_g: fx.domain.clone(),
            commitment: fx.commitments[0].root(),
            alpha: fx.alpha.clone(),
            sigma,
            ood_eval: fx.commitments[0].ood_eval,
            num_queries: NUM_QUERIES,
        }
    }

    fn honest_proof(fx: &Fixture) -> (SumCheckGroupInstance<Cfg, G, F>, SumCheckGroupProof<Cfg, G, F>) {
        let inst = instance(fx, fx.polys[0].evaluate_msm(&fx.alpha));
        let witness = SumCheckGroupWitness { whir_commitment: fx.commitments[0].clone(), phantom_data: PhantomData };
        let proof = ProtoSumCheckGroup::prove(&inst, &witness, &mut Transcript::new(b"test"));
        (inst, proof)
    }

    fn verify(fx: &Fixture, inst: &SumCheckGroupInstance<Cfg, G, F>, proof: &SumCheckGroupProof<Cfg, G, F>) -> bool {
        ProtoSumCheckGroup::verify(proof, inst, &mut Transcript::new(b"test"), &fx.verifier_state)
    }

    /// Replays the verifier transcript up to the query challenges, returning (r_vec, y_vec).
    fn replay_challenges(inst: &SumCheckGroupInstance<Cfg, G, F>, proof: &SumCheckGroupProof<Cfg, G, F>) -> (Vec<F>, Vec<usize>) {
        let mut t = Transcript::new(b"test");
        t.append_message(b"protocol", b"sumcheck");
        t.append_u64(b"instance.m", inst.m as u64);
        t.append_u64(b"instance.l", inst.l as u64);
        t.append_message(b"commitment", &inst.commitment.to_sponge_bytes_as_vec());
        for a in &inst.alpha {
            append_to_transcript(&mut t, b"alpha", a);
        }
        append_to_transcript(&mut t, b"sigma", &inst.sigma);
        append_to_transcript(&mut t, b"ood_eval", &inst.ood_eval);
        let _xi = get_challenge::<F>(&mut t, b"xi");
        let mut r_vec = Vec::new();
        for g in &proof.round_messages {
            append_to_transcript(&mut t, b"g_i(0)", &g[0]);
            append_to_transcript(&mut t, b"g_i(1)", &g[1]);
            append_to_transcript(&mut t, b"g_i(2)", &g[2]);
            r_vec.push(get_challenge::<F>(&mut t, b"r_i"));
        }
        for c in &proof.g_poly.coeffs {
            append_to_transcript(&mut t, b"h_poly", c);
        }
        let y_vec = (0..inst.num_queries)
            .map(|_| (get_challenge_u64(&mut t, b"challenge_y") as usize) % inst.domain_g.len())
            .collect();
        (r_vec, y_vec)
    }

    /// Cheating prover for the claim sigma' = 2 * f(alpha): sends the sum-check round messages
    /// of 2f (which are consistent with sigma'), but the honest folded polynomial f(r, .) and
    /// honest openings of the committed f. Every round check, every Merkle check and the
    /// g_poly/opening consistency check pass; only the terminal check (Figure 2, Step 6) fails.
    fn forge_doubled_sigma(fx: &Fixture) -> (SumCheckGroupInstance<Cfg, G, F>, SumCheckGroupProof<Cfg, G, F>) {
        let two = F::ONE + F::ONE;
        let f = &fx.polys[0];
        let f2 = MultilinearPoly::new(f.coeffs.iter().map(|c| *c * two).collect());
        let inst = instance(fx, f.evaluate_msm(&fx.alpha) * two);

        let mut t = Transcript::new(b"test");
        t.append_message(b"protocol", b"sumcheck");
        t.append_u64(b"instance.m", inst.m as u64);
        t.append_u64(b"instance.l", inst.l as u64);
        t.append_message(b"commitment", &inst.commitment.to_sponge_bytes_as_vec());
        for a in &inst.alpha {
            append_to_transcript(&mut t, b"alpha", a);
        }
        append_to_transcript(&mut t, b"sigma", &inst.sigma);
        append_to_transcript(&mut t, b"ood_eval", &inst.ood_eval);
        let xi = get_challenge::<F>(&mut t, b"xi");

        // alpha-part from 2f (consistent with sigma'), OOD part from the committed f (consistent with y_0)
        let w0 = ood_point::<Cfg, F>(&inst.commitment, M);
        let s_tables = compute_all_S_tables(compute_Sl_poly(M, L, &f2, &fx.alpha).coeffs, L);
        let s_tables_w0 = compute_all_S_tables(compute_Sl_poly(M, L, f, &w0).coeffs, L);
        let u_values = [F::ZERO, F::ONE, two];
        let mut r_vec = Vec::new();
        let mut round_messages = Vec::new();
        for i in 1..=L {
            let g = compute_batched_gi_values(M, i, L, &fx.alpha, &s_tables, &w0, &s_tables_w0, xi, &r_vec, &u_values);
            append_to_transcript(&mut t, b"g_i(0)", &g[0]);
            append_to_transcript(&mut t, b"g_i(1)", &g[1]);
            append_to_transcript(&mut t, b"g_i(2)", &g[2]);
            r_vec.push(get_challenge::<F>(&mut t, b"r_i"));
            round_messages.push([g[0], g[1], g[2]]);
        }
        let g_poly = f.restrict(&r_vec);
        for c in &g_poly.coeffs {
            append_to_transcript(&mut t, b"h_poly", c);
        }
        let openings: Vec<_> = (0..NUM_QUERIES)
            .map(|_| ((get_challenge_u64(&mut t, b"challenge_y") as usize) % fx.domain.len(), r_vec.clone()))
            .collect();
        let (opening_proofs, evaluations) = fx.commitments[0]
            .open_at_field_points_batch(&openings, &GroupWhirProverMerkleState::new(MerkleProofStrategy::Compressed))
            .unwrap();
        let proof = SumCheckGroupProof { round_messages, g_poly, opening_proofs, evaluations, phantom_data: PhantomData };
        (inst, proof)
    }

    /// Re-opens the commitment at `indices` (instead of the transcript-derived queries).
    fn reopen(fx: &Fixture, idx: usize, indices: &[usize], r_vec: &[F], strategy: MerkleProofStrategy) -> (BatchOpeningProof<Cfg, G>, Vec<G>) {
        let openings: Vec<_> = indices.iter().map(|&y| (y, r_vec.to_vec())).collect();
        fx.commitments[idx]
            .open_at_field_points_batch(&openings, &GroupWhirProverMerkleState::new(strategy))
            .unwrap()
    }

    #[test]
    fn test_sumcheck_accepts_honest_proof() {
        let fx = fixture(1);
        let (inst, proof) = honest_proof(&fx);
        assert!(verify(&fx, &inst, &proof));
    }

    #[test]
    fn test_sumcheck_rejects_wrong_ood_evaluation() {
        let fx = fixture(1);
        let (inst, proof) = honest_proof(&fx);
        // The commitment's y_0 really is G(w_0) for w_0 derived from the root
        let w0 = ood_point::<Cfg, F>(&inst.commitment, M);
        assert_eq!(inst.ood_eval, fx.polys[0].evaluate_msm(&w0));
        // A commitment carrying a different y_0 must not verify against the same proof
        let mut wrong = instance(&fx, inst.sigma);
        wrong.ood_eval = inst.ood_eval + G::random(&mut StdRng::seed_from_u64(2));
        assert!(!verify(&fx, &wrong, &proof));
    }

    #[test]
    fn test_sumcheck_rejects_false_evaluation_claim() {
        let fx = fixture(1);
        let (inst, proof) = forge_doubled_sigma(&fx);
        // Sanity: the forged openings really are the transcript queries and are consistent
        // with g_poly, so this proof is rejected by the terminal check alone.
        let (_, y_vec) = replay_challenges(&inst, &proof);
        assert!(proof.opening_proofs.opens_exactly(&y_vec));
        assert!(!verify(&fx, &inst, &proof), "false claim 2*sigma must be rejected");
    }

    #[test]
    fn test_sumcheck_rejects_single_opening_for_many_queries() {
        let fx = fixture(1);
        let (inst, mut proof) = honest_proof(&fx);
        let (r_vec, y_vec) = replay_challenges(&inst, &proof);
        assert!(y_vec.iter().collect::<BTreeSet<_>>().len() > 1);
        let (op, ev) = reopen(&fx, 0, &y_vec[..1], &r_vec, MerkleProofStrategy::Compressed);
        proof.opening_proofs = op;
        proof.evaluations = ev;
        assert!(!verify(&fx, &inst, &proof), "1 opening for {} queries must be rejected", NUM_QUERIES);
    }

    #[test]
    fn test_sumcheck_rejects_openings_at_unqueried_indices() {
        let fx = fixture(1);
        let (inst, mut proof) = honest_proof(&fx);
        let (r_vec, y_vec) = replay_challenges(&inst, &proof);
        let queried: BTreeSet<usize> = y_vec.iter().copied().collect();
        let chosen: Vec<usize> = (0..fx.domain.len()).filter(|y| !queried.contains(y)).take(queried.len()).collect();
        let (op, ev) = reopen(&fx, 0, &chosen, &r_vec, MerkleProofStrategy::Compressed);
        proof.opening_proofs = op;
        proof.evaluations = ev;
        assert!(!verify(&fx, &inst, &proof));
    }

    #[test]
    fn test_sumcheck_uncompressed_openings_bound_to_queries() {
        let fx = fixture(1);
        let (inst, mut proof) = honest_proof(&fx);
        let (r_vec, y_vec) = replay_challenges(&inst, &proof);
        let (op, ev) = reopen(&fx, 0, &y_vec, &r_vec, MerkleProofStrategy::Uncompressed);
        proof.opening_proofs = op;
        proof.evaluations = ev;
        let state = GroupWhirVerifierMerkleState::<Cfg, G>::new(
            MerkleProofStrategy::Uncompressed,
            fx.verifier_state.leaf_hash_params.clone(),
            fx.verifier_state.two_to_one_params.clone(),
        );
        assert!(ProtoSumCheckGroup::verify(&proof, &inst, &mut Transcript::new(b"test"), &state));

        let (op, ev) = reopen(&fx, 0, &y_vec[..NUM_QUERIES - 1], &r_vec, MerkleProofStrategy::Uncompressed);
        proof.opening_proofs = op;
        proof.evaluations = ev;
        assert!(!ProtoSumCheckGroup::verify(&proof, &inst, &mut Transcript::new(b"test"), &state));
    }

    #[test]
    fn test_sumcheck_rejects_malformed_proofs() {
        let fx = fixture(1);
        let (inst, honest) = honest_proof(&fx);
        let bump = G::random(&mut StdRng::seed_from_u64(1));
        let tamper = |f: &dyn Fn(&mut SumCheckGroupProof<Cfg, G, F>)| {
            let mut p = SumCheckGroupProof {
                round_messages: honest.round_messages.clone(),
                g_poly: honest.g_poly.clone(),
                opening_proofs: honest.opening_proofs.clone(),
                evaluations: honest.evaluations.clone(),
                phantom_data: PhantomData,
            };
            f(&mut p);
            verify(&fx, &inst, &p)
        };

        assert!(!tamper(&|p| { p.round_messages.pop(); }), "truncated round messages");
        assert!(!tamper(&|p| p.round_messages.push(p.round_messages[0])), "extra round message");
        assert!(!tamper(&|p| p.round_messages[0][2] = p.round_messages[0][2] + bump), "tampered round message");
        assert!(!tamper(&|p| p.round_messages[L - 1][2] = p.round_messages[L - 1][2] + bump), "tampered last round");
        assert!(!tamper(&|p| p.g_poly.coeffs[0] = p.g_poly.coeffs[0] + bump), "tampered folded polynomial");
        assert!(!tamper(&|p| p.g_poly = MultilinearPoly::new(p.g_poly.coeffs[..2].to_vec())), "short folded polynomial");
        assert!(!tamper(&|p| p.evaluations[0] = p.evaluations[0] + bump), "tampered evaluation");
        assert!(!tamper(&|p| { p.evaluations.pop(); }), "missing evaluation");
        assert!(!tamper(&|p| p.opening_proofs.leaves[0][0] = p.opening_proofs.leaves[0][0] + bump), "tampered leaf");
        assert!(!tamper(&|p| { p.opening_proofs.leaves.pop(); }), "missing leaf");
        assert!(!tamper(&|p| p.opening_proofs.leaves[0].push(bump)), "oversized leaf");

        let wrong_sigma = instance(&fx, inst.sigma + bump);
        assert!(!verify(&fx, &wrong_sigma, &honest), "wrong claimed sigma");
        let mut more_queries = instance(&fx, inst.sigma);
        more_queries.num_queries = NUM_QUERIES + 8;
        assert!(!verify(&fx, &more_queries, &honest), "proof generated for fewer queries");
    }

    // ---------------------------- batch path --------------------------------

    fn batch_instance(fx: &Fixture, sigmas: Vec<G>) -> SumCheckGroupBatchInstance<Cfg, G, F> {
        SumCheckGroupBatchInstance {
            m: M,
            l: L,
            domain_g: fx.domain.clone(),
            commitments: fx.commitments.iter().map(|c| c.root()).collect(),
            alpha: fx.alpha.clone(),
            sigmas,
            num_queries: NUM_QUERIES,
        }
    }

    fn honest_batch(fx: &Fixture) -> (SumCheckGroupBatchInstance<Cfg, G, F>, SumCheckGroupBatchProof<Cfg, G, F>) {
        let inst = batch_instance(fx, fx.polys.iter().map(|p| p.evaluate_msm(&fx.alpha)).collect());
        let witness = SumCheckGroupBatchWitness { whir_commitments: fx.commitments.clone(), phantom_data: PhantomData };
        let (proof, _) = ProtoSumCheckGroup::prove_batch(&inst, &witness, &mut Transcript::new(b"test"));
        (inst, proof)
    }

    fn verify_b(fx: &Fixture, inst: &SumCheckGroupBatchInstance<Cfg, G, F>, proof: &SumCheckGroupBatchProof<Cfg, G, F>) -> bool {
        ProtoSumCheckGroup::verify_batch(proof, inst, &mut Transcript::new(b"test"), &fx.verifier_state).0
    }

    /// Replays the batch verifier transcript, returning (r_vec, y_vec).
    fn replay_batch(inst: &SumCheckGroupBatchInstance<Cfg, G, F>, proof: &SumCheckGroupBatchProof<Cfg, G, F>) -> (Vec<F>, Vec<usize>) {
        let mut t = Transcript::new(b"test");
        t.append_message(b"protocol", b"batch_sumcheck");
        t.append_u64(b"k", inst.commitments.len() as u64);
        t.append_u64(b"instance.m", inst.m as u64);
        t.append_u64(b"instance.l", inst.l as u64);
        for root in &inst.commitments {
            t.append_message(b"commitment", &root.to_sponge_bytes_as_vec());
        }
        for a in &inst.alpha {
            append_to_transcript(&mut t, b"alpha", a);
        }
        for s in &inst.sigmas {
            append_to_transcript(&mut t, b"sigma", s);
        }
        let _gamma = get_challenge::<F>(&mut t, b"gamma");
        let mut r_vec = Vec::new();
        for g in &proof.round_messages {
            append_to_transcript(&mut t, b"g_i(0)", &g[0]);
            append_to_transcript(&mut t, b"g_i(1)", &g[1]);
            append_to_transcript(&mut t, b"g_i(2)", &g[2]);
            r_vec.push(get_challenge::<F>(&mut t, b"r_i"));
        }
        for c in &proof.g_poly.coeffs {
            append_to_transcript(&mut t, b"h_poly", c);
        }
        let y_vec = (0..inst.num_queries)
            .map(|_| (get_challenge_u64(&mut t, b"challenge_y") as usize) % inst.domain_g.len())
            .collect();
        (r_vec, y_vec)
    }

    /// Batch analogue of `forge_doubled_sigma`: claims 2 * f_0(alpha) for the first oracle.
    fn forge_batch_doubled_sigma(fx: &Fixture) -> (SumCheckGroupBatchInstance<Cfg, G, F>, SumCheckGroupBatchProof<Cfg, G, F>) {
        let two = F::ONE + F::ONE;
        let k = fx.polys.len();
        let mut sigmas: Vec<G> = fx.polys.iter().map(|p| p.evaluate_msm(&fx.alpha)).collect();
        sigmas[0] = sigmas[0] * two;
        let inst = batch_instance(fx, sigmas);

        let mut t = Transcript::new(b"test");
        t.append_message(b"protocol", b"batch_sumcheck");
        t.append_u64(b"k", k as u64);
        t.append_u64(b"instance.m", inst.m as u64);
        t.append_u64(b"instance.l", inst.l as u64);
        for root in &inst.commitments {
            t.append_message(b"commitment", &root.to_sponge_bytes_as_vec());
        }
        for a in &inst.alpha {
            append_to_transcript(&mut t, b"alpha", a);
        }
        for s in &inst.sigmas {
            append_to_transcript(&mut t, b"sigma", s);
        }
        let gamma = get_challenge::<F>(&mut t, b"gamma");

        // Round messages of sum_i gamma^i f_i with f_0 replaced by 2 f_0
        let mut cheat = vec![G::zero(); 1 << M];
        let mut honest = vec![G::zero(); 1 << M];
        let mut gp = F::ONE;
        for (i, p) in fx.polys.iter().enumerate() {
            let scale = if i == 0 { gp * two } else { gp };
            for j in 0..1 << M {
                cheat[j] = cheat[j] + p.coeffs[j] * scale;
                honest[j] = honest[j] + p.coeffs[j] * gp;
            }
            gp = gp * gamma;
        }
        let cheat = MultilinearPoly::new(cheat);
        let honest = MultilinearPoly::new(honest);
        let s_tables = compute_all_S_tables(compute_Sl_poly(M, L, &cheat, &fx.alpha).coeffs, L);
        let u_values = [F::ZERO, F::ONE, two];
        let mut r_vec = Vec::new();
        let mut round_messages = Vec::new();
        for i in 1..=L {
            let g = compute_gi_values(M, i, L, &fx.alpha, &r_vec, &s_tables, &u_values);
            append_to_transcript(&mut t, b"g_i(0)", &g[0]);
            append_to_transcript(&mut t, b"g_i(1)", &g[1]);
            append_to_transcript(&mut t, b"g_i(2)", &g[2]);
            r_vec.push(get_challenge::<F>(&mut t, b"r_i"));
            round_messages.push([g[0], g[1], g[2]]);
        }
        let g_poly = honest.restrict(&r_vec);
        for c in &g_poly.coeffs {
            append_to_transcript(&mut t, b"h_poly", c);
        }
        let y_vec: Vec<usize> = (0..NUM_QUERIES)
            .map(|_| (get_challenge_u64(&mut t, b"challenge_y") as usize) % fx.domain.len())
            .collect();
        let (opening_proofs, evaluations) = (0..k)
            .map(|i| reopen(fx, i, &y_vec, &r_vec, MerkleProofStrategy::Compressed))
            .unzip();
        let proof = SumCheckGroupBatchProof { round_messages, g_poly, opening_proofs, evaluations, phantom_data: PhantomData };
        (inst, proof)
    }

    #[test]
    fn test_batch_sumcheck_accepts_honest_proof() {
        let fx = fixture(3);
        let (inst, proof) = honest_batch(&fx);
        assert!(verify_b(&fx, &inst, &proof));
    }

    #[test]
    fn test_batch_sumcheck_rejects_false_evaluation_claim() {
        let fx = fixture(3);
        let (inst, proof) = forge_batch_doubled_sigma(&fx);
        let (_, y_vec) = replay_batch(&inst, &proof);
        assert!(proof.opening_proofs.iter().all(|op| op.opens_exactly(&y_vec)));
        assert!(!verify_b(&fx, &inst, &proof));
    }

    #[test]
    fn test_batch_sumcheck_rejects_unbound_openings() {
        let fx = fixture(3);
        let (inst, honest) = honest_batch(&fx);
        let (r_vec, y_vec) = replay_batch(&inst, &honest);
        let clone = || SumCheckGroupBatchProof {
            round_messages: honest.round_messages.clone(),
            g_poly: honest.g_poly.clone(),
            opening_proofs: honest.opening_proofs.clone(),
            evaluations: honest.evaluations.clone(),
            phantom_data: PhantomData,
        };

        // Every oracle opened at only the first query
        let mut p = clone();
        for i in 0..fx.polys.len() {
            let (op, ev) = reopen(&fx, i, &y_vec[..1], &r_vec, MerkleProofStrategy::Compressed);
            p.opening_proofs[i] = op;
            p.evaluations[i] = ev;
        }
        assert!(!verify_b(&fx, &inst, &p), "1 opening per oracle for {} queries", NUM_QUERIES);

        // Only a non-first oracle is under-opened (previously only oracle 0's indices were read)
        let mut p = clone();
        let (op, ev) = reopen(&fx, 2, &y_vec[..1], &r_vec, MerkleProofStrategy::Compressed);
        p.opening_proofs[2] = op;
        p.evaluations[2] = ev;
        assert!(!verify_b(&fx, &inst, &p), "under-opened oracle 2");

        // Missing oracle opening
        let mut p = clone();
        p.opening_proofs.pop();
        p.evaluations.pop();
        assert!(!verify_b(&fx, &inst, &p), "missing oracle opening");

        // Truncated round messages
        let mut p = clone();
        p.round_messages.pop();
        assert!(!verify_b(&fx, &inst, &p), "truncated round messages");
    }
}
