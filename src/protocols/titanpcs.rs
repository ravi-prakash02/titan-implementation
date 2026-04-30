use crate::group_whir_committer::{
    GroupWhirCommitment, GroupWhirProverMerkleState, GroupWhirVerifierMerkleState,
    MerkleProofStrategy,
};
use crate::merkle_tree::{Sha256Compress, Sha256LeafHash, Sha256MerkleTreeParams};
use crate::multilinear::MultilinearPoly;
use crate::pastatypes;
use crate::protocols::groupbulletproof::{
    BulletProof, BulletProofGroup, BulletProofInstance, BulletProofParams, BulletProofWitness,
    ProtoBulletProofGroup,
};
use crate::protocols::sumcheck::{
    ProtoSumCheckGroup, SumCheckGroup, SumCheckGroupBatchInstance, SumCheckGroupBatchProof,
    SumCheckGroupBatchWitness, SumCheckGroupInstance, SumCheckGroupProof, SumCheckGroupWitness,
};
use crate::protocols::field_sumcheck;
use crate::titantranscript::{append_to_transcript, get_challenge};
use crate::traits::{ByteSerializable, InnerProduct, Linear};
use crate::utils::create_smooth_domain;
use merlin::Transcript;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;
use std::fmt::Debug;
use std::hash::Hash;
use std::marker::PhantomData;
use std::time::Instant;

type PallasG = crate::pastatypes::Point;
type PallasF = crate::pastatypes::Scalar;

type ArkPallasG = crate::arkpallastypes::ArkPallasPoint;

pub trait PolyCommitment<G, F> {
    type CommitmentWithOpen;
    type Commitment;
    type ProverParams;
    type VerifierParams;
    type EvalProof;

    //fn setup(config: Self::SetupConfig) -> (Self::ProverParams, Self::VerifierParams);

    fn commit(pp: &Self::ProverParams, poly: &MultilinearPoly<F>) -> Self::CommitmentWithOpen;

    fn prove(
        pp: &Self::ProverParams,
        poly: &MultilinearPoly<F>,
        comm_with_open: &Self::CommitmentWithOpen,
        alpha: &[F],
    ) -> Self::EvalProof;

    fn verify(
        pp_v: &Self::VerifierParams,
        comm: &Self::Commitment,
        alpha: &[F],
        sigma: F,
        proof: &Self::EvalProof,
    ) -> bool;
}

pub struct TitanPolyCommitment<G, F> {
    phantom: PhantomData<(G, F)>,
}

#[derive(Clone)]
pub struct TitanSetupConfig {
    pub m: usize,  // size of polynomial
    pub m1: usize, // size of group polynomial
    pub l1: usize, // folding factor for group poly
    // m2 = m - m1
    pub domain_g1_size: usize, // eval domain for folded polynomial
    pub num_queries: usize,    // number of queries to group oracle
    pub l2: usize,             // folding factor for generator oracle
    pub domain_g2_size: usize, // eval domain for folded generator oracle
    pub num_merkle_nodes: usize, // number of merkle nodes transmitted
}

impl TitanSetupConfig {
    // returns only proof_size for now
    pub fn estimate_metrics(&self) -> usize {
        // proof size
        let mut proof_elems: usize = 0;
        let m2 = self.m - self.m1;
        let coset_size = 1usize << self.l1;
        let folded_size = 1usize << (self.m1 - self.l1);
        proof_elems += 3 * self.l1; // sum check messages
        proof_elems += self.num_queries * coset_size; // openings for Group oracle
        proof_elems += folded_size; // size of folded poly
        // contribution from bullet proofs
        proof_elems += 1usize << (m2 - self.l2); // size of folded poly
        proof_elems += 4 * self.l2; // bullet proof messages

        // contribution of merkle authentication
        proof_elems += self.num_merkle_nodes;

        println!("Proof size = {}", 32 * proof_elems);
        println!(
            "Merkle tree size, queries = {}, {}",
            self.domain_g1_size, self.num_queries
        );
        32 * proof_elems // returns proof size in bytes
    }
}

#[derive(Clone)]
pub struct TitanSetup<G, F> {
    pub m: usize,  // size of field polynomial
    pub m1: usize, // number of variables of group polynomial
    pub l1: usize, // folding factor for group polynomial
    pub num_queries: usize,
    pub domain_g1: Vec<F>,                 // domain for folded group oracle
    pub bp_setup: BulletProofParams<G, F>, // bulletproof setup
}

pub struct TitanProverParams<G, F>
where
    G: Linear<F> + InnerProduct<F, Output = G> + ByteSerializable + PartialEq + Send + Sync,
    F: ff::Field,
{
    pub pp: TitanSetup<G, F>,
    pub state: GroupWhirProverMerkleState,
    pub verifier_state: GroupWhirVerifierMerkleState<Sha256MerkleTreeParams<G>, G>,
}

pub struct TitanVerifierParams<G, F>
where
    G: Linear<F> + InnerProduct<F, Output = G> + ByteSerializable + PartialEq + Send + Sync,
    F: ff::Field,
{
    pub pp: TitanSetup<G, F>,
    pub state: GroupWhirVerifierMerkleState<Sha256MerkleTreeParams<G>, G>,
}

pub struct TitanEvalProof<G, F>
where
    G: Linear<F> + InnerProduct<F, Output = G> + ByteSerializable + PartialEq + Send + Sync,
    F: ff::Field,
{
    pub comm_partial: G,
    pub partial_eval_proof: SumCheckGroupProof<Sha256MerkleTreeParams<G>, G, F>,
    pub bullet_proof: BulletProof<G, F>,
}

impl<G, F> TitanEvalProof<G, F>
where
    G: Linear<F> + InnerProduct<F, Output=G> + ByteSerializable + PartialEq + Send + Sync,
    F: ff::Field,
{
    pub fn get_proof_size(&self) -> usize {
        let mut proof_size = 0usize;
        proof_size += 32; //for G
        proof_size += self.partial_eval_proof.get_proof_size()*32; //*32 to convert it inot bytes
        proof_size += self.bullet_proof.get_proof_size()*32; //*32 to convert it inot bytes
        proof_size
    }
}

pub struct AggregateEvalProof<G, F>
where
    G: Linear<F> + InnerProduct<F, Output = G> + ByteSerializable + PartialEq + Send + Sync,
    F: ff::Field,
{
    pub evaluations: Vec<F>,
    pub sumcheck_proof: field_sumcheck::SumCheckProof<F>,
    pub p_at_r: F,
    pub eval_proof: TitanEvalProof<G, F>,
}

impl<G, F> PolyCommitment<G, F> for TitanPolyCommitment<G, F>
where
    G: Linear<F> + InnerProduct<F, Output = G> + ByteSerializable + Debug + PartialEq + Send + Sync,
    G::Affine: InnerProduct<F, Output = G> + Send + Sync,
    F: ff::Field + ByteSerializable + Linear<F> + InnerProduct<F, Output = F>,
{
    type CommitmentWithOpen = GroupWhirCommitment<Sha256MerkleTreeParams<G>, G>;
    type Commitment =
        <Sha256MerkleTreeParams<G> as ark_crypto_primitives::merkle_tree::Config>::InnerDigest;
    type ProverParams = TitanProverParams<G, F>;
    type VerifierParams = TitanVerifierParams<G, F>;
    type EvalProof = TitanEvalProof<G, F>;

    fn commit(
        prover_params: &Self::ProverParams,
        poly: &MultilinearPoly<F>,
    ) -> Self::CommitmentWithOpen {
        // use pedersen commitment to commit to q=(1 << m1) sized chunks.
        assert_eq!(prover_params.pp.m, poly.num_vars);
        let chunk_size = 1usize << (prover_params.pp.m - prover_params.pp.m1);
        let chunks: Vec<&[F]> = poly.coeffs.chunks_exact(chunk_size).collect();
        let start = Instant::now();
        let c_poly_affine = G::to_affine(&prover_params.pp.bp_setup.c_poly.coeffs);
        let partial_comms: Vec<G> = chunks
            .into_iter()
            .map(|chunk| G::Affine::inner_product_msm(&c_poly_affine, chunk))
            .collect();
        //let partial_comms:Vec<G> = chunks.into_par_iter().map(|chunk| G::inner_product_orbit(&prover_params.pp.bp_setup.c_poly.coeffs, chunk)).collect();

        println!(
            "Time to compute G polynomial: {}",
            start.elapsed().as_millis()
        );
        let G_poly = MultilinearPoly::new(partial_comms);
        // construct whir commitment
        let start = Instant::now();
        let whir_commitment = GroupWhirCommitment::<Sha256MerkleTreeParams<G>, G>::new(
            &G_poly,
            prover_params.pp.l1,
            prover_params.pp.domain_g1.clone(),
            &prover_params.verifier_state.leaf_hash_params,
            &prover_params.verifier_state.two_to_one_params,
        );
        println!(
            "Time to compute whir commitment: {}",
            start.elapsed().as_millis()
        );

        whir_commitment.unwrap()
    }

    fn prove(
        prover_params: &Self::ProverParams,
        poly: &MultilinearPoly<F>,
        comm_with_open: &Self::CommitmentWithOpen,
        alpha: &[F],
    ) -> Self::EvalProof {
        // Slice alpha as (alpha_x, alpha_y)
        // Compute evaluation of poly sigma
        // Compute commitment sigma_partial to partial evaluated poly g(x) = \sum_{i} f(X,i)eq(i,\alpha_y)
        // Generate evaluation proof for sigma_partial
        // Generate evaluation proof for sigma
        let m1 = prover_params.pp.m1;
        let m2 = prover_params.pp.m - m1;
        let l1 = prover_params.pp.l1;

        let start = Instant::now();
        let (alpha_x, alpha_y) = (alpha[..m2].to_vec(), alpha[m2..].to_vec());
        let a_poly = poly.fold_msb(&alpha_y);
        let sigma = a_poly.evaluate(&alpha_x);
        println!("Time to Evaluation {} msec", start.elapsed().as_millis());
        let start = Instant::now();
        let sigma_partial = comm_with_open.poly.evaluate_msm(&alpha_y);
        println!(
            "Time to compute partial commitment = {} msec",
            start.elapsed().as_millis()
        );

        // Create a Group sumcheck instance for proving partial evaluation
        let instance = SumCheckGroupInstance::<Sha256MerkleTreeParams<G>, G, F> {
            m: m1,
            l: l1,
            domain_g: prover_params.pp.domain_g1.clone(),
            commitment: comm_with_open.root().clone(),
            alpha: alpha_y.clone(),
            sigma: sigma_partial,
            num_queries: prover_params.pp.num_queries,
        };

        let witness = SumCheckGroupWitness {
            whir_commitment: comm_with_open.clone(),
            phantom_data: Default::default(),
        };

        let mut transcript = Transcript::new(b"sumcheck");
        let start = Instant::now();
        let proof = ProtoSumCheckGroup::prove(&instance, &witness, &mut transcript);
        println!(
            "Time to generate group sumcheck proof = {} msec",
            start.elapsed().as_millis()
        );

        // Next generate the bulletproof proof
        let bp_instance = BulletProofInstance::<G, F> {
            pp: prover_params.pp.bp_setup.clone(),
            commitment: sigma_partial,
            alpha: alpha_x,
            sigma: sigma,
        };

        let bp_witness = BulletProofWitness::<G, F> {
            pp: prover_params.pp.bp_setup.clone(),
            a_poly,
        };

        let start = Instant::now();
        let proof_bp =
            ProtoBulletProofGroup::<G, F>::prove(&bp_instance, &bp_witness, &mut transcript);
        println!(
            "Time to do bulletproofs {} msec",
            start.elapsed().as_millis()
        );

        TitanEvalProof {
            comm_partial: sigma_partial,
            partial_eval_proof: proof,
            bullet_proof: proof_bp,
        }
    }

    fn verify(
        pp_v: &Self::VerifierParams,
        comm: &Self::Commitment,
        alpha: &[F],
        sigma: F,
        proof: &Self::EvalProof,
    ) -> bool {
        let mut transcript = Transcript::new(b"sumcheck");
        // first verify the group sum-check
        let m1 = pp_v.pp.m1;
        let l1 = pp_v.pp.l1;
        let m2 = pp_v.pp.m - m1;
        let (alpha_x, alpha_y) = (alpha[..m2].to_vec(), alpha[m2..].to_vec());
        let start = Instant::now();
        let instance = SumCheckGroupInstance::<Sha256MerkleTreeParams<G>, G, F> {
            m: m1,
            l: l1,
            domain_g: pp_v.pp.domain_g1.clone(),
            commitment: comm.clone(),
            alpha: alpha_y.clone(),
            sigma: proof.comm_partial,
            num_queries: pp_v.pp.num_queries,
        };
        println!("Time to create instance {}", start.elapsed().as_millis());

        let start = Instant::now();
        let is_valid = ProtoSumCheckGroup::verify(
            &proof.partial_eval_proof,
            &instance,
            &mut transcript,
            &pp_v.state,
        );
        assert!(is_valid);
        let merkle_nodes = proof.partial_eval_proof.merkle_transmitted_nodes();
        println!(
            "Group sumcheck verification took {} msec, merkle proof transmitted nodes = {}",
            start.elapsed().as_millis(),
            merkle_nodes,
        );

        // now verify bulletproof
        let start = Instant::now();
        let bp_instance = BulletProofInstance::<G, F> {
            pp: pp_v.pp.bp_setup.clone(),
            commitment: proof.comm_partial,
            alpha: alpha_x,
            sigma: sigma,
        };
        println!("Time to create bp instance {}", start.elapsed().as_millis());

        let start = Instant::now();
        let is_valid =
            ProtoBulletProofGroup::verify(&proof.bullet_proof, &bp_instance, &mut transcript);
        assert!(is_valid);
        println!("Bulletproof verification took {} msec", start.elapsed().as_millis());

        is_valid
    }
}

impl<G, F> TitanPolyCommitment<G, F>
where
    G: Linear<F> + InnerProduct<F, Output = G> + ByteSerializable + Debug + PartialEq + Send + Sync,
    G::Affine: InnerProduct<F, Output = G> + Send + Sync,
    F: ff::Field + ByteSerializable + Linear<F> + InnerProduct<F, Output = F>,
{
    /// This function commits to a list of polynomials p_1,...,p_k as aggregate oracle.
    /// First we flatten p_0,...,p_{k-1} into a multilinear polynomial p
    /// We generate a WHIR commitment of p given the setup parameter
    /// We have p(x,0000)=p_0, p(x,1000)=p_1,...,p(x,1111) = p_{k-1} assuming k is power of 2.
    pub fn aggregate_commit(
        prover_params: &<TitanPolyCommitment<G, F> as PolyCommitment<G, F>>::ProverParams,
        polys: &Vec<MultilinearPoly<F>>,
    ) -> <TitanPolyCommitment<G, F> as PolyCommitment<G, F>>::CommitmentWithOpen {
        let all_coeffs: Vec<Vec<F>> = polys.iter().map(|p| p.coeffs.clone()).collect();
        let mut p_coeffs: Vec<F> = all_coeffs.into_iter().flatten().collect();
        p_coeffs.resize(p_coeffs.len().next_power_of_two(), F::ZERO);
        let p_agg = MultilinearPoly::new(p_coeffs);
        <TitanPolyCommitment<G,F> as PolyCommitment<G,F>>::commit(prover_params, &p_agg)
    }

    /// Generate evaluation proof given aggregate polynomial p for polynomials p_1,...,p_k
    /// Output values v_0,...,v_{max_idx-1}
    /// Proof for p_i(point) = v_i for i=0...max_idx-1
    ///
    /// The aggregate polynomial has n = m + s variables, where m = point.len() and
    /// s = n - m (selector bits). Layout: p(x, bin(i)) = p_i(x).
    ///
    /// We aggregate k claims into a single degree-2 field sum-check:
    ///   Σ_x p(x) · q(x) = v_agg
    /// where q(x,y) = eq(x,point) · h(y), h(bin(i)) = γ^i, v_agg = Σ_i γ^i v_i.
    /// After the sum-check yields challenge point r, we provide a TitanPCS eval proof for p(r).
    pub fn aggregate_prove(
        prover_params: &<TitanPolyCommitment<G, F> as PolyCommitment<G, F>>::ProverParams,
        poly: &MultilinearPoly<F>,
        comm_with_open: &<TitanPolyCommitment<G,F> as PolyCommitment<G,F>>::CommitmentWithOpen,
        point: &[F],
        max_idx: usize,
    ) -> AggregateEvalProof<G,F>
    {
        let n = poly.num_vars;
        let m = point.len();
        let s = n - m; // number of selector bits

        // Step 1: Compute evaluations v_i = p(point, bin(i)) for i = 0..max_idx
        // Fold by point (first m variables) to get a poly in s variables whose
        // coefficients are exactly the evaluations at the Boolean hypercube.
        let mut folded = poly.clone();
        for &r in point.iter() {
            folded.fold_first(r);
        }
        let evaluations: Vec<F> = folded.coeffs[..max_idx].to_vec();

        // Step 2: Fiat-Shamir — append evaluations to transcript, derive γ
        let mut transcript = Transcript::new(b"aggregate_sumcheck");
        for v in &evaluations {
            append_to_transcript(&mut transcript, b"eval", v);
        }
        let gamma: F = get_challenge::<F>(&mut transcript, b"gamma");

        // Step 3: Compute aggregated claim v_agg = Σ_i γ^i · v_i
        let mut v_agg = F::ZERO;
        let mut gamma_power = F::ONE;
        for v in &evaluations {
            v_agg += gamma_power * *v;
            gamma_power *= gamma;
        }

        // Step 4: Construct q polynomial
        // q(x, y) = eq(x, point) · h(y) where h(bin(i)) = γ^i for i < max_idx, else 0
        let eq_point = MultilinearPoly::<F>::init_with_eq(point);

        let h_size = 1usize << s;
        let mut h_coeffs = vec![F::ZERO; h_size];
        gamma_power = F::ONE;
        for i in 0..max_idx {
            h_coeffs[i] = gamma_power;
            gamma_power *= gamma;
        }

        // Build q via tensor product: q.coeffs[j * 2^m + i] = eq_point[i] · h[j]
        let m_size = 1usize << m;
        let mut q_coeffs = vec![F::ZERO; 1 << n];
        for j in 0..h_size {
            for i in 0..m_size {
                q_coeffs[j * m_size + i] = eq_point.coeffs[i] * h_coeffs[j];
            }
        }
        let q_poly = MultilinearPoly::new(q_coeffs);

        // Step 5: Field sum-check for Σ_x p(x) · q(x) = v_agg (degree 2)
        let combine_fn = |vals: &[F]| -> F { vals[0] * vals[1] };
        let (sumcheck_proof, challenges) = field_sumcheck::prove(
            vec![poly.clone(), q_poly],
            &combine_fn,
            2,
            &mut transcript,
        );

        // Step 6: Evaluate p at the challenge point and generate TitanPCS eval proof
        let p_at_r = poly.evaluate(&challenges);
        let eval_proof = <TitanPolyCommitment<G, F> as PolyCommitment<G, F>>::prove(
            prover_params,
            poly,
            comm_with_open,
            &challenges,
        );

        AggregateEvalProof {
            evaluations,
            sumcheck_proof,
            p_at_r,
            eval_proof,
        }
    }

    /// Verify an aggregate evaluation proof.
    ///
    /// Checks that for i = 0..max_idx-1: p_i(point) = proof.evaluations[i],
    /// where the p_i are the sub-polynomials committed inside the aggregate oracle.
    ///
    /// Parameters:
    ///   - pp_v: verifier params (set up for the aggregate polynomial size)
    ///   - comm: Merkle root of the aggregate polynomial commitment
    ///   - point: evaluation point for each sub-polynomial (m variables)
    ///   - num_vars: total number of variables in the aggregate polynomial (m + s)
    ///   - max_idx: number of sub-polynomials evaluated
    ///   - proof: the aggregate evaluation proof
    pub fn aggregate_verify(
        pp_v: &<TitanPolyCommitment<G, F> as PolyCommitment<G, F>>::VerifierParams,
        comm: &<TitanPolyCommitment<G, F> as PolyCommitment<G, F>>::Commitment,
        point: &[F],
        num_vars: usize,
        max_idx: usize,
        proof: &AggregateEvalProof<G, F>,
    ) -> bool
    {
        let m = point.len();
        let s = num_vars - m;

        // Step 1: Replay Fiat-Shamir — append evaluations to transcript, derive γ
        let mut transcript = Transcript::new(b"aggregate_sumcheck");
        for v in &proof.evaluations {
            append_to_transcript(&mut transcript, b"eval", v);
        }
        let gamma: F = get_challenge::<F>(&mut transcript, b"gamma");

        // Step 2: Compute aggregated claim v_agg = Σ_i γ^i · v_i
        let mut v_agg = F::ZERO;
        let mut gamma_power = F::ONE;
        for v in &proof.evaluations {
            v_agg += gamma_power * *v;
            gamma_power *= gamma;
        }

        // Step 3: Verify the field sum-check
        let subclaim = field_sumcheck::verify(
            &proof.sumcheck_proof,
            v_agg,
            num_vars,
            2, // degree
            &mut transcript,
        );

        let subclaim = match subclaim {
            Ok(sc) => sc,
            Err(e) => {
                println!("Field sum-check verification failed: {}", e);
                return false;
            }
        };

        let r = &subclaim.point;

        // Step 4: Compute q(r) directly
        // q(r) = eq(r[..m], point) · h(r[m..])
        // where h(r_y) = Σ_{i=0}^{max_idx-1} γ^i · eq(r_y, bin(i))
        let r_x = &r[..m];
        let r_y = &r[m..];

        // eq(r_x, point) = Π_j (r_x[j] · point[j] + (1 - r_x[j]) · (1 - point[j]))
        let mut eq_rx_point = F::ONE;
        for j in 0..m {
            eq_rx_point *= r_x[j] * point[j] + (F::ONE - r_x[j]) * (F::ONE - point[j]);
        }

        // h(r_y) = Σ_{i=0}^{max_idx-1} γ^i · eq(r_y, bin(i))
        let mut h_at_ry = F::ZERO;
        gamma_power = F::ONE;
        for i in 0..max_idx {
            // eq(r_y, bin(i)) = Π_{bit=0}^{s-1} (r_y[bit] · bit_val + (1 - r_y[bit]) · (1 - bit_val))
            let mut eq_ry_bin_i = F::ONE;
            for bit in 0..s {
                let bit_val = if (i >> bit) & 1 == 1 { F::ONE } else { F::ZERO };
                eq_ry_bin_i *= r_y[bit] * bit_val + (F::ONE - r_y[bit]) * (F::ONE - bit_val);
            }
            h_at_ry += gamma_power * eq_ry_bin_i;
            gamma_power *= gamma;
        }

        let q_at_r = eq_rx_point * h_at_ry;

        // Step 5: Check sum-check oracle consistency
        if proof.p_at_r * q_at_r != subclaim.expected_evaluation {
            println!("Oracle check failed: p(r) * q(r) != expected_evaluation");
            return false;
        }

        // Step 6: Verify TitanPCS eval proof for p(r) = p_at_r
        let is_valid = <TitanPolyCommitment<G, F> as PolyCommitment<G, F>>::verify(
            pp_v,
            comm,
            &subclaim.point,
            proof.p_at_r,
            &proof.eval_proof,
        );

        if !is_valid {
            println!("TitanPCS eval proof verification failed");
            return false;
        }

        true
    }
}

pub struct TitanBatchEvalProof<G, F>
where
    G: Linear<F> + InnerProduct<F, Output = G> + ByteSerializable + PartialEq + Send + Sync,
    F: ff::Field,
{
    pub comm_partials: Vec<G>,
    pub partial_eval_proof: SumCheckGroupBatchProof<Sha256MerkleTreeParams<G>, G, F>,
    pub bullet_proof: BulletProof<G, F>,
}

impl<G, F> TitanBatchEvalProof<G, F>
where
    G: Linear<F> + InnerProduct<F, Output=G> + ByteSerializable + PartialEq + Send + Sync,
    F: ff::Field,
{
    pub fn get_prooof_size(&self) -> usize {
        self.partial_eval_proof.get_proof_size() + self.bullet_proof.get_proof_size()
    }

    pub fn get_honest_oracle_proof_size(&self) -> usize {
        self.partial_eval_proof.g_poly.coeffs.len()
    }
}

impl<G, F> TitanPolyCommitment<G, F>
where
    G: Linear<F> + InnerProduct<F, Output = G> + ByteSerializable + Debug + PartialEq + Send + Sync,
    G::Affine: InnerProduct<F, Output = G> + Send + Sync,
    F: ff::Field + ByteSerializable + Linear<F> + InnerProduct<F, Output = F>,
{
    /// Batch prove: k polynomials each with their own independent commitment,
    /// sharing a single evaluation point alpha.
    pub fn batch_prove(
        prover_params: &<TitanPolyCommitment<G, F> as PolyCommitment<G, F>>::ProverParams,
        polys: &[MultilinearPoly<F>],
        comm_with_opens: &[<TitanPolyCommitment<G, F> as PolyCommitment<G, F>>::CommitmentWithOpen],
        alpha: &[F],
    ) -> TitanBatchEvalProof<G, F> {
        let k = polys.len();
        assert_eq!(k, comm_with_opens.len());
        let m1 = prover_params.pp.m1;
        let m2 = prover_params.pp.m - m1;
        let l1 = prover_params.pp.l1;

        let (alpha_x, alpha_y) = (alpha[..m2].to_vec(), alpha[m2..].to_vec());

        // For each polynomial: compute a_poly_i, sigma_i, sigma_partial_i
        let mut a_polys: Vec<MultilinearPoly<F>> = Vec::with_capacity(k);
        let mut sigmas: Vec<F> = Vec::with_capacity(k);
        let mut comm_partials: Vec<G> = Vec::with_capacity(k);
        let mut sigma_partials: Vec<G> = Vec::with_capacity(k);

        for i in 0..k {
            let a_poly_i = polys[i].fold_msb(&alpha_y);
            let sigma_i = a_poly_i.evaluate(&alpha_x);
            let sigma_partial_i = comm_with_opens[i].poly.evaluate_msm(&alpha_y);
            a_polys.push(a_poly_i);
            sigmas.push(sigma_i);
            comm_partials.push(sigma_partial_i);
            sigma_partials.push(sigma_partial_i);
        }

        // Build batch instance + witness
        let commitments: Vec<_> = comm_with_opens.iter().map(|c| c.root()).collect();
        let instance = SumCheckGroupBatchInstance::<Sha256MerkleTreeParams<G>, G, F> {
            m: m1,
            l: l1,
            domain_g: prover_params.pp.domain_g1.clone(),
            commitments,
            alpha: alpha_y.clone(),
            sigmas: sigma_partials.clone(),
            num_queries: prover_params.pp.num_queries,
        };

        let witness = SumCheckGroupBatchWitness {
            whir_commitments: comm_with_opens.to_vec(),
            phantom_data: Default::default(),
        };

        let mut transcript = Transcript::new(b"batch_sumcheck");
        let (batch_proof, gamma) =
            ProtoSumCheckGroup::prove_batch(&instance, &witness, &mut transcript);

        // One bulletproof on batched claim
        let mut sigma_partial_batched = G::zero();
        let mut a_poly_batched_coeffs = vec![F::ZERO; a_polys[0].coeffs.len()];
        let mut sigma_batched = F::ZERO;
        let mut gamma_power = F::ONE;

        for i in 0..k {
            sigma_partial_batched = sigma_partial_batched + sigma_partials[i] * gamma_power;
            sigma_batched = sigma_batched + sigmas[i] * gamma_power;
            for j in 0..a_poly_batched_coeffs.len() {
                a_poly_batched_coeffs[j] = a_poly_batched_coeffs[j] + a_polys[i].coeffs[j] * gamma_power;
            }
            gamma_power = gamma_power * gamma;
        }
        let a_poly_batched = MultilinearPoly::new(a_poly_batched_coeffs);

        let bp_instance = BulletProofInstance::<G, F> {
            pp: prover_params.pp.bp_setup.clone(),
            commitment: sigma_partial_batched,
            alpha: alpha_x,
            sigma: sigma_batched,
        };

        let bp_witness = BulletProofWitness::<G, F> {
            pp: prover_params.pp.bp_setup.clone(),
            a_poly: a_poly_batched,
        };

        let proof_bp =
            ProtoBulletProofGroup::<G, F>::prove(&bp_instance, &bp_witness, &mut transcript);

        TitanBatchEvalProof {
            comm_partials,
            partial_eval_proof: batch_proof,
            bullet_proof: proof_bp,
        }
    }

    /// Batch verify: verify that k polynomials (each with their own commitment) all
    /// evaluate correctly at a shared point alpha.
    pub fn batch_verify(
        pp_v: &<TitanPolyCommitment<G, F> as PolyCommitment<G, F>>::VerifierParams,
        comms: &[<TitanPolyCommitment<G, F> as PolyCommitment<G, F>>::Commitment],
        alpha: &[F],
        sigmas: &[F],
        proof: &TitanBatchEvalProof<G, F>,
    ) -> bool {
        let k = comms.len();
        assert_eq!(k, sigmas.len());
        assert_eq!(k, proof.comm_partials.len());
        let m1 = pp_v.pp.m1;
        let l1 = pp_v.pp.l1;
        let m2 = pp_v.pp.m - m1;
        let (alpha_x, alpha_y) = (alpha[..m2].to_vec(), alpha[m2..].to_vec());

        // Build batch instance
        let instance = SumCheckGroupBatchInstance::<Sha256MerkleTreeParams<G>, G, F> {
            m: m1,
            l: l1,
            domain_g: pp_v.pp.domain_g1.clone(),
            commitments: comms.to_vec(),
            alpha: alpha_y.clone(),
            sigmas: proof.comm_partials.clone(),
            num_queries: pp_v.pp.num_queries,
        };

        let mut transcript = Transcript::new(b"batch_sumcheck");
        let (is_valid, gamma) = ProtoSumCheckGroup::verify_batch(
            &proof.partial_eval_proof,
            &instance,
            &mut transcript,
            &pp_v.state,
        );
        assert!(is_valid);

        // One bulletproof verification on batched claim
        let mut sigma_partial_batched = G::zero();
        let mut sigma_batched = F::ZERO;
        let mut gamma_power = F::ONE;
        for i in 0..k {
            sigma_partial_batched = sigma_partial_batched + proof.comm_partials[i] * gamma_power;
            sigma_batched = sigma_batched + sigmas[i] * gamma_power;
            gamma_power = gamma_power * gamma;
        }

        let bp_instance = BulletProofInstance::<G, F> {
            pp: pp_v.pp.bp_setup.clone(),
            commitment: sigma_partial_batched,
            alpha: alpha_x,
            sigma: sigma_batched,
        };

        let is_valid =
            ProtoBulletProofGroup::verify(&proof.bullet_proof, &bp_instance, &mut transcript);
        assert!(is_valid);

        is_valid
    }
}

use crate::rand_scalar;
use ark_std::UniformRand;
use ff::PrimeField;

fn titan_config_for_m(m: usize) -> TitanSetupConfig {
    match m {
        18 => TitanSetupConfig {
            m,
            m1: 8,
            l1: 1,
            domain_g1_size: 11,
            num_queries: 70,
            l2: 5,
            domain_g2_size: 9,
            num_merkle_nodes: 0,
        },
        20 => TitanSetupConfig {
            m,
            m1: 8,
            l1: 1,
            domain_g1_size: 11,
            num_queries: 70,
            l2: 6,
            domain_g2_size: 11,
            num_merkle_nodes: 0,
        },
        22 => TitanSetupConfig {
            m,
            m1: 9,
            l1: 2,
            domain_g1_size: 11,
            num_queries: 70,
            l2: 6,
            domain_g2_size: 12,
            num_merkle_nodes: 0,
        },
        24 => TitanSetupConfig {
            m,
            m1: 10,
            l1: 2,
            domain_g1_size: 12,
            num_queries: 70,
            l2: 7,
            domain_g2_size: 12,
            num_merkle_nodes: 0,
        },
        26 => TitanSetupConfig {
            m,
            m1: 11,
            l1: 3,
            domain_g1_size: 12,
            num_queries: 70,
            l2: 8,
            domain_g2_size: 12,
            num_merkle_nodes: 0,
        },
        _ => panic!("Unsupported Titan parameter m = {}", m),
    }
}

fn run_titan_once(m: usize) -> (u128, u128, u128, usize) {
    let mut rng = StdRng::from_entropy();
    // Load config
    let mut config = titan_config_for_m(m);

    let titan_pp = TitanPolyCommitment::<PallasG, PallasF>::setup(&config);

    // generate random polynomial
    //let f_coeffs = (0..(1usize << config.m)).into_iter().map(|_| rand_scalar(&mut rng)).collect::<Vec<_>>();
    //let f_coeffs = (0..(1usize << config.m))
    //     .into_iter()
    //     .map(|_| PallasF::from_u128(u32::rand(&mut rng) as u128))
    //     .collect::<Vec<_>>();
    let f_coeffs = (0..(1usize << config.m))
            .into_iter()
            .map(|_| rand_scalar(&mut rng))
            .collect::<Vec<_>>();

    let f_poly = MultilinearPoly::new(f_coeffs);
    let alpha = (0..config.m)
        .map(|_| rand_scalar(&mut rng))
        .collect::<Vec<_>>();
    let sigma = f_poly.evaluate(&alpha);
    // Compute commitment to the polynomial
    let start = Instant::now();
    let comm = TitanPolyCommitment::commit(&titan_pp.0, &f_poly);
    let commit_ms = start.elapsed().as_millis();
    //println!("Time to commit to polynomial {} msec", start.elapsed().as_millis());

    // Generate proof
    //let mut transcript = Transcript::new(b"TitanProver");
    let start = Instant::now();
    let proof = TitanPolyCommitment::prove(&titan_pp.0, &f_poly, &comm, &alpha);
    let prove_ms = start.elapsed().as_millis();
    //println!("Time to generate proof = {} msec", start.elapsed().as_millis());

    let start = Instant::now();
    let _res = TitanPolyCommitment::verify(&titan_pp.1, &comm.root, &alpha, sigma, &proof);
    let verify_ms = start.elapsed().as_millis();
    //println!("Time to verify proof = {} msec", start.elapsed().as_millis());

    //let prover_ms = commit_ms + prove_ms;

    config.num_merkle_nodes = proof.partial_eval_proof.merkle_transmitted_nodes();

    // Proof size
    //let proof_bytes = config.estimate_metrics();
    let proof_size = proof.get_proof_size();

    (commit_ms, prove_ms, verify_ms, proof_size)
}

mod tests {
    use super::*;
    use crate::rand_scalar;
    use ark_std::UniformRand;
    use ff::{Field, PrimeField};

    #[test]
    fn test_titan_prover() {
        let mut rng = StdRng::from_entropy();
        let m = 20usize;
        let m1 = 8usize;
        let l1 = 1usize;
        let domain_g1_size = m1 - l1 + 4;
        let m2 = m - m1;
        let l2 = 5usize;
        let domain_g2_size = m2 - l2 + 5;
        let orbits = 6usize;
        let mut config: TitanSetupConfig = TitanSetupConfig {
            m: m,
            m1: m1,
            l1: l1,
            domain_g1_size, // this should ideally be m1 - l1 + 3
            num_queries: 32,
            l2: l2,
            domain_g2_size, // this should ideally be m2 - l2 + 4, where m2 = m - m1
            num_merkle_nodes: 0,    // we'll set it before calling estimate metrics
        };

        let titan_pp = TitanPolyCommitment::<PallasG, PallasF>::setup(&config);
        let table: Vec<PallasF> = (0..(1usize << orbits))
            .into_iter()
            .map(|_| PallasF::random(&mut rng))
            .collect();

        // generate random polynomial
        let f_coeffs = (0..(1usize << config.m))
            .into_iter()
            .map(|_| rand_scalar(&mut rng))
            .collect::<Vec<_>>();
        //let f_coeffs = (0..(1usize << config.m)).into_iter().map(|_| table[usize::rand(&mut rng) % table.len()]).collect::<Vec<_>>();

        let f_poly = MultilinearPoly::new(f_coeffs);
        let alpha = (0..config.m)
            .map(|_| rand_scalar(&mut rng))
            .collect::<Vec<_>>();
        let sigma = f_poly.evaluate(&alpha);
        // Compute commitment to the polynomial
        let start = Instant::now();
        let comm = TitanPolyCommitment::commit(&titan_pp.0, &f_poly);
        println!(
            "Time to commit to polynomial {} msec",
            start.elapsed().as_millis()
        );

        // Generate proof
        //let mut transcript = Transcript::new(b"TitanProver");
        let start = Instant::now();
        let proof = TitanPolyCommitment::prove(&titan_pp.0, &f_poly, &comm, &alpha);
        println!(
            "Time to generate proof = {} msec",
            start.elapsed().as_millis()
        );

        let start = Instant::now();
        let _res = TitanPolyCommitment::verify(&titan_pp.1, &comm.root, &alpha, sigma, &proof);
        println!(
            "Time to verify proof = {} msec",
            start.elapsed().as_millis()
        );

        config.num_merkle_nodes = proof.partial_eval_proof.merkle_transmitted_nodes();

        config.estimate_metrics();
    }

    #[test]
    fn test_titan_aggregate_proof() {
        let mut rng = StdRng::from_entropy();
        let m = 22usize;
        let m1 = 9usize;
        let l1 = 2usize;
        let domain_g1_size = m1 - l1 + 4;
        let m2 = m - m1 +2;
        let l2 = 6usize;
        let domain_g2_size = m2 - l2 + 5;
        let orbits = 6usize;
        let mut config: TitanSetupConfig = TitanSetupConfig {
            m: m + 2,
            m1: m1,
            l1: l1,
            domain_g1_size, // this should ideally be m1 - l1 + 3
            num_queries: 64,
            l2: l2,
            domain_g2_size, // this should ideally be m2 - l2 + 4, where m2 = m - m1
            num_merkle_nodes: 0,    // we'll set it later before calling estimate metrics
        };

        let titan_pp = TitanPolyCommitment::<PallasG, PallasF>::setup(&config);
        // Generate 4 polynomial, 2 with short coefficients and 2 with random coefficients.
        let coeffs_1 = (0..(1usize << m)).into_iter().map(|_| PallasF::random(&mut rng)).collect::<Vec<_>>();
        let coeffs_2 = (0..(1usize << m)).into_iter().map(|_| PallasF::random(&mut rng)).collect::<Vec<_>>();
        let coeffs_3 = (0..(1usize << m)).into_iter().map(|_| PallasF::from_u128(u32::rand(&mut rng) as u128)).collect::<Vec<_>>();
        let coeffs_4 = (0..(1usize << m)).into_iter().map(|_| PallasF::from_u128(u32::rand(&mut rng) as u128)).collect::<Vec<_>>();

        let poly_1 = MultilinearPoly::new(coeffs_1.clone());
        let poly_2 = MultilinearPoly::new(coeffs_2.clone());
        let poly_3 = MultilinearPoly::new(coeffs_3.clone());
        let poly_4 = MultilinearPoly::new(coeffs_4.clone());

        let coeffs: Vec<PallasF> = [coeffs_1.as_slice(), coeffs_2.as_slice(), coeffs_3.as_slice(), coeffs_4.as_slice()].concat();
        let agg_poly = MultilinearPoly::new(coeffs);

        let start = Instant::now();
        let agg_comm = TitanPolyCommitment::<PallasG, PallasF>::aggregate_commit(
            &titan_pp.0,
            &vec![poly_1, poly_2, poly_3, poly_4],
        );
        println!("Time to genenerate aggregate commitment = {}", start.elapsed().as_millis());

        let point = (0usize..m).into_iter().map(|_| PallasF::random(&mut rng)).collect::<Vec<_>>();
        let start = Instant::now();
        let agg_proof = TitanPolyCommitment::aggregate_prove(
            &titan_pp.0,
            &agg_poly,
            &agg_comm,
            &point,
            4
        );
        println!("Time to prove evaluations = {}", start.elapsed().as_millis());

        let start = Instant::now();
        let res = TitanPolyCommitment::<PallasG, PallasF>::aggregate_verify(
            &titan_pp.1,
            &agg_comm.root,
            &point,
            m+2,
            4,
            &agg_proof
        );
        println!("Time to verify agg proof = {}", start.elapsed().as_millis());
        config.num_merkle_nodes = agg_proof.eval_proof.partial_eval_proof.merkle_transmitted_nodes();
        config.estimate_metrics();

        assert_eq!(res, true, "failed to verify proof");
    }


    #[test]
    fn test_titan_batch_proof() {
        let mut rng = StdRng::from_entropy();
        let m = 22usize;
        let m1 = 9usize;
        let l1 = 2usize;
        let domain_g1_size = m1 - l1 + 4;
        let m2 = m - m1;
        let l2 = 6usize;
        let domain_g2_size = m2 - l2 + 5;
        let mut config = TitanSetupConfig {
            m,
            m1,
            l1,
            domain_g1_size,
            num_queries: 64,
            l2,
            domain_g2_size,
            num_merkle_nodes: 0,
        };

        let titan_pp = TitanPolyCommitment::<PallasG, PallasF>::setup(&config);
        let k = 2usize;

        // Generate k polynomials and commit each independently
        let mut polys = Vec::with_capacity(k);
        let mut comms = Vec::with_capacity(k);
        for _ in 0..k {
            let coeffs: Vec<PallasF> = (0..(1usize << m))
                .map(|_| PallasF::ONE)
                .collect();
            let poly = MultilinearPoly::new(coeffs);
            let comm = TitanPolyCommitment::commit(&titan_pp.0, &poly);
            polys.push(poly);
            comms.push(comm);
        }

        // Shared random evaluation point
        let alpha: Vec<PallasF> = (0..m)
            .map(|_| rand_scalar(&mut rng))
            .collect();

        // Compute expected evaluations
        let sigmas: Vec<PallasF> = polys.iter().map(|p| p.evaluate(&alpha)).collect();

        // Batch prove
        let start = Instant::now();
        let batch_proof = TitanPolyCommitment::batch_prove(
            &titan_pp.0,
            &polys,
            &comms,
            &alpha,
        );
        println!(
            "Time to generate batch proof = {} msec",
            start.elapsed().as_millis()
        );

        // Batch verify
        let comm_roots: Vec<_> = comms.iter().map(|c| c.root.clone()).collect();
        let start = Instant::now();
        let res = TitanPolyCommitment::<PallasG, PallasF>::batch_verify(
            &titan_pp.1,
            &comm_roots,
            &alpha,
            &sigmas,
            &batch_proof,
        );
        println!(
            "Time to verify batch proof = {} msec",
            start.elapsed().as_millis()
        );

        //for eval_proof in batch_proof.partial_eval_proof.merkle_transmitted_nodes()
        config.num_merkle_nodes = 1100;
        config.estimate_metrics();
        assert!(res, "Batch proof verification failed");
    }

    #[test]
    #[ignore]
    fn titan_scaling_experiment() {
        let mut results = Vec::new();

        for m in (18..=26).step_by(2) {
            let runs = 4;
            let mut c_sum: u128 = 0;
            let mut p_sum: u128 = 0;
            let mut v_sum: u128 = 0;
            let mut proof_size: usize = 0;

            for _ in 0..runs {
                let (c, p, v, size) = run_titan_once(m);
                c_sum += c;
                p_sum += p;
                v_sum += v;
                proof_size += size; 
            }

            results.push((m, c_sum / runs as u128, p_sum / runs as u128, v_sum / runs as u128, proof_size / runs as usize));
        }

        println!("m,commit_ms,prove_ms,verifier_ms,proof_bytes");
        for (m, c, p, v, size) in &results {
            println!("{},{},{},{},{}", m, c, p, v, size);
        }
    }
}
