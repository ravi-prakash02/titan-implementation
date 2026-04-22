use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::merkle_tree::{MerkleTree as ArkMerkleTree, MultiPath, Path};
use pasta_curves::group::ff::Field;
use pasta_curves::group::{Curve, Group};
use rayon::iter::ParallelIterator;
/// GroupWhir Coset-wise Polynomial Commitment Scheme
///
/// Main implementation of the GroupWhir commitment scheme with support for:
/// - Multiple hash functions (SHA256, Keccak, Blake3)
/// - Flexible proof strategies (Compressed, Uncompressed)
/// - Batch opening operations
/// - Group element multilinear polynomials
use std::borrow::Borrow;
use std::time::Instant;

use crate::merkle_tree::{MerkleTreeParams, Sha256Compress, Sha256Digest, Sha256LeafHash};
use crate::multilinear::MultilinearPoly;
use crate::traits::{ByteSerializable, InnerProduct, Linear};
use crate::utils::{eq_all, eq_from_index, multilinear_fft};
use ark_serialize::CanonicalSerialize;
use ark_std::log2;
use rand::thread_rng;
use rayon::prelude::{IntoParallelIterator, IntoParallelRefIterator};
use std::collections::HashMap;
// ============================================================================
// Merkle Proof Strategy
// ============================================================================

#[derive(Debug, Clone, Copy)]
pub enum MerkleProofStrategy {
    Uncompressed,
    Compressed,
}

// ============================================================================
// Merkle State Managers
// ============================================================================

/// Prover-side state for generating Merkle proofs with a fixed strategy
pub struct GroupWhirProverMerkleState {
    strategy: MerkleProofStrategy,
}
#[derive(Clone)]
pub enum MerkleMultiProof<C: ark_crypto_primitives::merkle_tree::Config> {
    Compressed(MultiPath<C>),
    Uncompressed(Vec<Path<C>>),
}

impl<C: ark_crypto_primitives::merkle_tree::Config> MerkleMultiProof<C> {
    /// Count the number of internal tree nodes transmitted in the proof.
    /// For Compressed (MultiPath): leaf_siblings + suffix nodes (prefix-encoded).
    /// For Uncompressed (Vec<Path>): leaf_sibling + auth_path nodes per path.
    pub fn transmitted_node_count(&self) -> usize {
        match self {
            MerkleMultiProof::Compressed(mp) => {
                let suffix_nodes: usize = mp.auth_paths_suffixes.iter().map(|s| s.len()).sum();
                let leaf_siblings = mp.leaf_siblings_hashes.len();
                suffix_nodes + leaf_siblings
            }
            MerkleMultiProof::Uncompressed(paths) => {
                paths.iter().map(|p| 1 + p.auth_path.len()).sum()
            }
        }
    }
}

impl GroupWhirProverMerkleState {
    pub fn new(strategy: MerkleProofStrategy) -> Self {
        Self { strategy }
    }

    pub fn generate_proof<C: ark_crypto_primitives::merkle_tree::Config>(
        &self,
        tree: &ArkMerkleTree<C>,
        leaf_idx: usize,
    ) -> Result<ark_crypto_primitives::merkle_tree::Path<C>, String> {
        tree.generate_proof(leaf_idx)
            .map_err(|e| format!("Merkle proof generation failed: {:?}", e))
    }
    pub fn generate_multi_proof<C: ark_crypto_primitives::merkle_tree::Config>(
        &self,
        tree: &ArkMerkleTree<C>,
        indices: &[usize],
    ) -> Result<MerkleMultiProof<C>, String> {
        match self.strategy {
            MerkleProofStrategy::Compressed => tree
                .generate_multi_proof(indices.iter().copied())
                .map(MerkleMultiProof::Compressed)
                .map_err(|e| format!("Compressed multi-proof generation failed: {:?}", e)),
            MerkleProofStrategy::Uncompressed => indices
                .iter()
                .map(|&idx| self.generate_proof(tree, idx))
                .collect::<Result<Vec<_>, _>>()
                .map(MerkleMultiProof::Uncompressed),
        }
    }
}

/// Verifier-side state for verifying Merkle proofs with a fixed strategy
pub struct GroupWhirVerifierMerkleState<C, G>
where
    C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>>,
{
    pub strategy: MerkleProofStrategy,
    pub leaf_hash_params: <C::LeafHash as ark_crypto_primitives::crh::CRHScheme>::Parameters,
    pub two_to_one_params:
        <C::TwoToOneHash as ark_crypto_primitives::crh::TwoToOneCRHScheme>::Parameters,
    _phantom: std::marker::PhantomData<(C, G)>,
}

impl<C, G> GroupWhirVerifierMerkleState<C, G>
where
    C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>>,
    G: ByteSerializable,
{
    pub fn new(
        strategy: MerkleProofStrategy,
        leaf_hash_params: <C::LeafHash as ark_crypto_primitives::crh::CRHScheme>::Parameters,
        two_to_one_params: <C::TwoToOneHash as ark_crypto_primitives::crh::TwoToOneCRHScheme>::Parameters,
    ) -> Self {
        Self {
            strategy,
            leaf_hash_params,
            two_to_one_params,
            _phantom: std::marker::PhantomData,
        }
    }

    pub fn verify_proof(
        &self,
        leaf: &Vec<G>,
        auth_path: &ark_crypto_primitives::merkle_tree::Path<C>,
        root: &C::InnerDigest,
    ) -> Result<bool, String> {
        auth_path
            .verify(&self.leaf_hash_params, &self.two_to_one_params, root, leaf)
            .map_err(|e| format!("Path verification failed: {:?}", e))
    }

    pub fn verify_multi_proof(
        &self,
        leaves: &[Vec<G>],
        multi_proof: &MerkleMultiProof<C>,
        root: &C::InnerDigest,
    ) -> Result<bool, String> {
        match multi_proof {
            MerkleMultiProof::Compressed(multi_path) => {
                /*
                println!("DEBUG: Compressed verification");
                println!("DEBUG: Number of leaves passed: {}", leaves.len());
                println!(
                    "DEBUG: Number of leaf_indexes in MultiPath: {}",
                    multi_path.leaf_indexes.len()
                );
                println!("DEBUG: leaf_indexes: {:?}", multi_path.leaf_indexes);
                
                 */
                multi_path
                    .verify(
                        &self.leaf_hash_params,
                        &self.two_to_one_params,
                        root,
                        leaves.iter(),
                    )
                    .map_err(|e| format!("Compressed multi-proof verification failed: {:?}", e))
            }
            MerkleMultiProof::Uncompressed(auth_paths) => {
                println!("DEBUG: Uncompressed verification");
                println!("DEBUG: Number of leaves: {}", leaves.len());
                println!("DEBUG: Number of paths: {}", auth_paths.len());

                if leaves.len() != auth_paths.len() {
                    return Err("Leaf and path count mismatch".to_string());
                }

                for (leaf, path) in leaves.iter().zip(auth_paths.iter()) {
                    if !self.verify_proof(leaf, path, root)? {
                        return Ok(false);
                    }
                }

                Ok(true)
            }
        }
    }
}

// ============================================================================
// Commitment Proof and Main Commitment Structure
// ============================================================================

/// Opening proof for a single commitment
#[derive(Debug, Clone)]
pub struct OpeningProof<C: ark_crypto_primitives::merkle_tree::Config, G> {
    pub leaf: Vec<G>,
    pub auth_path: ark_crypto_primitives::merkle_tree::Path<C>,
    pub leaf_idx: usize,
}

#[derive(Clone)]
pub struct BatchOpeningProof<C: ark_crypto_primitives::merkle_tree::Config, G> {
    pub leaves: Vec<Vec<G>>,
    pub multi_proof: MerkleMultiProof<C>,
}

impl<C: ark_crypto_primitives::merkle_tree::Config, G> BatchOpeningProof<C, G> {
    /// Returns the leaf indices in the order used by the Merkle proof.
    /// For Compressed proofs, these are deduplicated and sorted by the BTreeSet
    /// inside arkworks' `generate_multi_proof`. For Uncompressed proofs, they
    /// preserve the original query order.
    pub fn leaf_indices(&self) -> Vec<usize> {
        match &self.multi_proof {
            MerkleMultiProof::Compressed(mp) => mp.leaf_indexes.clone(),
            MerkleMultiProof::Uncompressed(paths) => {
                paths.iter().map(|p| p.leaf_index).collect()
            }
        }
    }
}

/// The GroupWhir commitment scheme for multilinear polynomials
///
/// Generic over any Merkle tree configuration supporting group element leaves.
#[derive(Clone)]
pub struct GroupWhirCommitment<C, G>
where
    C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>>,
{
    pub poly: MultilinearPoly<G>, // multilinear polynomial
    pub k: usize,                 // domain folding dimension
    pub tree: ArkMerkleTree<C>,   // Merkle tree over coset-wise leaves
    pub root: C::InnerDigest,     // Merkle root
    pub leaves: Vec<Vec<G>>,      // leaves for each coset
}

impl<C, G> GroupWhirCommitment<C, G>
where
    C: ark_crypto_primitives::merkle_tree::Config<Leaf = Vec<G>> + Send,
    G: ByteSerializable,
{
    /// Create a commitment to a multilinear polynomial with coset-wise structure
    /// we require additional trait bound on G, as restrict function requires certain properties
    /// with respect to underlying field.
    pub fn new<F: ff::Field>(
        poly: &MultilinearPoly<G>,
        k: usize,
        domain_points: Vec<F>,
        leaf_hash_params: &<C::LeafHash as CRHScheme>::Parameters,
        two_to_one_params: &<C::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
    ) -> Result<Self, String>
    where
        G: Linear<F> + InnerProduct<F, Output = G> + Sync + Send,
        G::Affine: InnerProduct<F, Output = G>,
    {
        assert!(k < poly.num_vars, "k must be less than number of variables");

        let mut leaves = vec![vec![G::zero(); 1usize << k]; domain_points.len()];
        // We re-model the polynomial so that restrictions f(x,y) are easy to compute
        let slice_size = 1usize << (poly.num_vars - k);
        let num_slices = 1usize << k;
        let mut remodel_coeffs = vec![G::zero(); 1usize << poly.num_vars];
        for i in 0..num_slices {
            for j in 0..slice_size {
                remodel_coeffs[i * slice_size + j] = poly.coeffs[j * num_slices + i];
            }
        }

        let all_leaves: Vec<Vec<G>> = (0..num_slices)
            .into_par_iter()
            .map(|i| {
                let leaf = multilinear_fft(
                    &remodel_coeffs[i * slice_size..(i + 1) * slice_size],
                    &domain_points,
                    poly.num_vars - k,
                    log2(domain_points.len()) as usize,
                );
                leaf
            })
            .collect();

        for i in 0..domain_points.len() {
            for j in 0..(1usize << k) {
                leaves[i][j] = all_leaves[j][i];
            }
        }

        // Build Merkle tree from leaves
        let tree = ArkMerkleTree::new(leaf_hash_params, two_to_one_params, &leaves)
            .map_err(|e| format!("Failed to create Merkle tree: {:?}", e))?;

        let root = tree.root();

        Ok(GroupWhirCommitment {
            poly: MultilinearPoly::new(poly.coeffs.clone()),
            k,
            tree,
            root,
            leaves,
        })
    }

    pub fn root(&self) -> C::InnerDigest {
        self.root.clone()
    }

    /// Open the commitment at an arbitrary field point
    pub fn open_at_field_point<F>(
        &self,
        y_idx: usize,
        alpha: &[F],
        prover_state: &GroupWhirProverMerkleState,
    ) -> Result<(OpeningProof<C, G>, G), String>
    where
        G: Linear<F> + InnerProduct<F, Output = G>,
        F: ff::Field,
    {
        if alpha.len() != self.k {
            return Err(format!(
                "Alpha length {} does not match k={}",
                alpha.len(),
                self.k
            ));
        }

        if y_idx >= self.leaves.len() {
            return Err(format!("Invalid y index: {}", y_idx));
        }

        let leaf = self.leaves[y_idx].clone();
        let auth_path = prover_state.generate_proof(&self.tree, y_idx)?;

        let proof = OpeningProof {
            leaf: leaf.clone(),
            auth_path,
            leaf_idx: y_idx,
        };

        // Compute: f(α, y, y^2, ...) = Σ_{b∈{0,1}^k} eq̃(α, b) · leaf[b]
        let eq_values = MultilinearPoly::init_with_eq(alpha).coeffs;
        let mut evaluation = G::zero();

        for b_idx in 0..(1 << self.k) {
            evaluation = evaluation + leaf[b_idx] * eq_values[b_idx];
        }

        Ok((proof, evaluation))
    }

    /// Open at multiple field points (batch)
    pub fn open_at_field_points_batch<F: ff::Field>(
        &self,
        openings: &[(usize, Vec<F>)],
        prover_state: &GroupWhirProverMerkleState,
    ) -> Result<(BatchOpeningProof<C, G>, Vec<G>), String>
    where
        G: Linear<F> + InnerProduct<F, Output = G>,
    {
        // All indices (including duplicates)
        let indices: Vec<usize> = openings.iter().map(|(idx, _)| *idx).collect();

        // Generate multi-proof for all indices
        let multi_proof = prover_state.generate_multi_proof(&self.tree, &indices)?;

        let eq_values = MultilinearPoly::init_with_eq(&openings[0].1).coeffs;

        let leaf_indices: Vec<usize> = match &multi_proof {
            MerkleMultiProof::Compressed(mp) => mp.leaf_indexes.clone(),
            MerkleMultiProof::Uncompressed(paths) => paths.iter().map(|p| p.leaf_index).collect(),
        };

        let mut leaves = Vec::new();
        let mut evaluations = Vec::new();

        for idx in &leaf_indices {
            let leaf = self.leaves[*idx].clone();
            let evaluation = G::inner_product_msm(&leaf, &eq_values);
            leaves.push(leaf);
            evaluations.push(evaluation);
        }

        Ok((
            BatchOpeningProof {
                leaves,
                multi_proof,
            },
            evaluations,
        ))
    }

    /// Verify a single opening proof
    pub fn verify_opening<F: ff::Field>(
        proof: &OpeningProof<C, G>,
        alpha: &[F],
        claimed_evaluation: &G,
        root: &C::InnerDigest,
        verifier_state: &GroupWhirVerifierMerkleState<C, G>,
    ) -> Result<bool, String> {
        // Step 1: Verify Merkle path
        if !verifier_state.verify_proof(&proof.leaf, &proof.auth_path, root)? {
            return Ok(false);
        }

        Ok(true)
    }

    /// Verify multiple opening proofs (batch)
    pub fn verify_openings_batch<F: ff::Field>(
        batch_proof: &BatchOpeningProof<C, G>,
        alpha: &[F],
        claimed_evaluations: &[G],
        root: &C::InnerDigest,
        verifier_state: &GroupWhirVerifierMerkleState<C, G>,
    ) -> Result<bool, String>
    where
        G: Linear<F> + InnerProduct<F, Output = G> + PartialEq,
    {
        let num_openings = batch_proof.leaves.len();

        if num_openings != claimed_evaluations.len() {
            return Err("Proof, alpha, and evaluation counts must match".to_string());
        }
        println!("DEBUG verify_openings_batch:");
        println!("  num_openings (leaves): {}", num_openings);
        println!("  claimed_evaluations.len(): {}", claimed_evaluations.len());
        println!("  alpha.len(): {}", alpha.len());

        //  Verify Merkle multi-proof
        if !verifier_state.verify_multi_proof(
            &batch_proof.leaves,
            &batch_proof.multi_proof,
            root,
        )? {
            return Ok(false);
        }

        // Batch verify claimed evaluations
        let eq_vec = MultilinearPoly::init_with_eq(alpha).coeffs;
        let r = F::random(&mut thread_rng()); //this needs to come from the prover transcript
        //println!("DEBUG: eq_vec.len(): {}", eq_vec.len());
        //println!("DEBUG: leaf[0].len(): {}", batch_proof.leaves[0].len());

        let leaf_size = eq_vec.len();
        let mut combined_leaves: Vec<G> = Vec::with_capacity(num_openings * leaf_size);
        let mut combined_eq: Vec<F> = Vec::with_capacity(num_openings * leaf_size);
        let mut r_powers: Vec<F> = Vec::with_capacity(num_openings);

        let mut r_power = F::ONE;
        for i in 0..num_openings {
            combined_leaves.extend(batch_proof.leaves[i].iter().cloned());
            let scaled_eq: Vec<F> = eq_vec.iter().map(|x| *x * r_power).collect();
            combined_eq.extend(scaled_eq);
            r_powers.push(r_power);
            r_power = r_power * r;
        }

        let lhs = G::inner_product_msm(&combined_leaves, &combined_eq);
        let rhs = G::inner_product_msm(claimed_evaluations, &r_powers);

        //println!("DEBUG: Evaluation check lhs == rhs: {}", lhs == rhs);

        Ok(lhs == rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::merkle_tree::Sha256MerkleTreeParams;
    use crate::utils::create_smooth_domain;
    use pasta_curves::pallas::{Point as G, Scalar as F};
    use rand::prelude::StdRng;
    use rand::SeedableRng;

    type Sha256Config = Sha256MerkleTreeParams<G>;

    #[test]
    fn test_GroupWhir_single_opening() {
        let mut rng = StdRng::from_entropy();

        let m = 3;
        let coeffs: Vec<G> = (0..1 << m).map(|_| G::random(&mut rng)).collect();

        let poly = MultilinearPoly::new(coeffs);
        let k = 1;

        let domain_points: Vec<F> = (0..4).map(|i| F::from(i as u64 + 1)).collect();

        //@todo this should be internal to GroupWhirCommitment
        let (leaf_hash_params, two_to_one_params) =
            crate::merkle_tree::default_config::<G, Sha256LeafHash<G>, Sha256Compress>(&mut rng);

        // Pass parameters to new()
        let commitment = GroupWhirCommitment::<Sha256Config, G>::new(
            &poly,
            k,
            domain_points,
            &leaf_hash_params,
            &two_to_one_params,
        )
        .expect("Failed to create commitment");
        let root = commitment.root();

        let alpha: Vec<F> = (0..k).map(|_| F::random(&mut rng)).collect();

        //@todo prover and verifier state seem like internal
        let prover_state = GroupWhirProverMerkleState::new(MerkleProofStrategy::Compressed);
        let verifier_state = GroupWhirVerifierMerkleState::<Sha256Config, G>::new(
            MerkleProofStrategy::Compressed,
            leaf_hash_params,
            two_to_one_params,
        );

        let (proof, evaluation) = commitment
            .open_at_field_point(0, &alpha, &prover_state)
            .expect("Failed to open");

        let is_valid = GroupWhirCommitment::<Sha256Config, G>::verify_opening(
            &proof,
            &alpha,
            &evaluation,
            &root,
            &verifier_state,
        )
        .expect("Verification failed");

        assert!(is_valid);
    }

    #[test]
    fn test_GroupWhir_batch_openings() {
        let mut rng = StdRng::from_entropy();

        let m = 3;
        let coeffs: Vec<G> = (0..1 << m).map(|_| G::random(&mut rng)).collect();
        let poly = MultilinearPoly::new(coeffs);

        let k = 2;
        let domain_points: Vec<F> = (0..4).map(|i| F::from(i as u64 + 1)).collect();

        let (leaf_hash_params, two_to_one_params) =
            crate::merkle_tree::default_config::<G, Sha256LeafHash<G>, Sha256Compress>(&mut rng);

        let commitment = GroupWhirCommitment::<Sha256Config, G>::new(
            &poly,
            k,
            domain_points,
            &leaf_hash_params,
            &two_to_one_params,
        )
        .expect("Failed to create commitment");

        let root = commitment.root();

        let prover_state = GroupWhirProverMerkleState::new(MerkleProofStrategy::Compressed);
        let verifier_state = GroupWhirVerifierMerkleState::<Sha256Config, G>::new(
            MerkleProofStrategy::Compressed,
            leaf_hash_params,
            two_to_one_params,
        );

        // Use same alpha for all openings
        let alpha: Vec<F> = (0..k).map(|_| F::random(&mut rng)).collect();
        let openings = vec![(0, alpha.clone()), (1, alpha.clone()), (2, alpha.clone())];

        let (batch_proof, evaluations) = commitment
            .open_at_field_points_batch(&openings, &prover_state)
            .expect("Failed to batch open");

        let is_valid = GroupWhirCommitment::<Sha256Config, G>::verify_openings_batch(
            &batch_proof,
            &alpha,
            &evaluations,
            &root,
            &verifier_state,
        )
        .expect("Batch verification failed");

        assert!(is_valid);
    }
    #[test]
    fn test_GroupWhir_larger_polynomial_m6() {
        let mut rng = StdRng::from_entropy();

        // Larger polynomial: 2^6 = 64 coefficients
        let m = 6;
        let coeffs: Vec<G> = (0..1 << m).map(|_| G::random(&mut rng)).collect();

        let poly = MultilinearPoly::new(coeffs);
        let k = 3;

        let domain_points: Vec<F> = (0..8).map(|i| F::from(i as u64 + 1)).collect();

        let (leaf_hash_params, two_to_one_params) =
            crate::merkle_tree::default_config::<G, Sha256LeafHash<G>, Sha256Compress>(&mut rng);

        let start = std::time::Instant::now();
        let commitment = GroupWhirCommitment::<Sha256Config, G>::new(
            &poly,
            k,
            domain_points,
            &leaf_hash_params,
            &two_to_one_params,
        )
        .expect("Failed to create commitment");
        let root = commitment.root();
        let duration = start.elapsed();
        println!("time to compute commitment is {:?}", duration);

        let alpha: Vec<F> = (0..k).map(|_| F::random(&mut rng)).collect();

        let start = std::time::Instant::now();
        let prover_state = GroupWhirProverMerkleState::new(MerkleProofStrategy::Uncompressed);
        let verifier_state = GroupWhirVerifierMerkleState::<Sha256Config, G>::new(
            MerkleProofStrategy::Uncompressed,
            leaf_hash_params,
            two_to_one_params,
        );

        let (proof, evaluation) = commitment
            .open_at_field_point(0, &alpha, &prover_state)
            .expect("Failed to open");
        let duration = start.elapsed();
        println!("time to compute opening proof is {:?}", duration);

        let is_valid = GroupWhirCommitment::<Sha256Config, G>::verify_opening(
            &proof,
            &alpha,
            &evaluation,
            &root,
            &verifier_state,
        )
        .expect("Verification failed");

        assert!(is_valid);
    }
    #[test]
    fn test_GroupWhir_large_batch_m8() {
        let mut rng = StdRng::from_entropy();

        // Very large polynomial: 2^8 = 256 coefficients
        let m = 8;
        let coeffs: Vec<G> = (0..1 << m).map(|_| G::random(&mut rng)).collect();
        let poly = MultilinearPoly::new(coeffs);

        let k = 4;
        let domain_points: Vec<F> = (0..16).map(|i| F::from(i as u64 + 1)).collect();

        let (leaf_hash_params, two_to_one_params) =
            crate::merkle_tree::default_config::<G, Sha256LeafHash<G>, Sha256Compress>(&mut rng);

        let start = std::time::Instant::now();
        let commitment = GroupWhirCommitment::<Sha256Config, G>::new(
            &poly,
            k,
            domain_points,
            &leaf_hash_params,
            &two_to_one_params,
        )
        .expect("Failed to create commitment");

        let root = commitment.root();
        let duration = start.elapsed();
        println!("time to generate commitment {:?}", duration);

        let prover_state = GroupWhirProverMerkleState::new(MerkleProofStrategy::Compressed);
        let verifier_state = GroupWhirVerifierMerkleState::<Sha256Config, G>::new(
            MerkleProofStrategy::Compressed,
            leaf_hash_params,
            two_to_one_params,
        );

        // Use same alpha for all openings
        let alpha: Vec<F> = (0..k).map(|_| F::random(&mut rng)).collect();
        let mut openings = Vec::new();
        for leaf_idx in 0..8 {
            openings.push((leaf_idx, alpha.clone()));
        }

        let (batch_proof, evaluations) = commitment
            .open_at_field_points_batch(&openings, &prover_state)
            .expect("Failed to batch open");

        let is_valid = GroupWhirCommitment::<Sha256Config, G>::verify_openings_batch(
            &batch_proof,
            &alpha,
            &evaluations,
            &root,
            &verifier_state,
        )
        .expect("Batch verification failed");

        assert!(is_valid);
    }

    #[test]
    fn test_GroupWhir_m10_single_opening() {
        let mut rng = StdRng::from_entropy();

        // Very large polynomial: 2^10 = 1024 coefficients
        let m = 8;
        let coeffs: Vec<G> = (0..1 << m).map(|_| G::random(&mut rng)).collect();
        let poly = MultilinearPoly::new(coeffs);

        let k = 2;
        let domain_points = create_smooth_domain(9);

        let (leaf_hash_params, two_to_one_params) =
            crate::merkle_tree::default_config::<G, Sha256LeafHash<G>, Sha256Compress>(&mut rng);

        let start = Instant::now();
        let commitment = GroupWhirCommitment::<Sha256Config, G>::new(
            &poly,
            k,
            domain_points.1,
            &leaf_hash_params,
            &two_to_one_params,
        )
        .expect("Failed to create commitment");

        let root = commitment.root();
        let duration = start.elapsed();
        println!("time to compute commitment is {:?}", duration);

        let prover_state = GroupWhirProverMerkleState::new(MerkleProofStrategy::Uncompressed);
        let verifier_state = GroupWhirVerifierMerkleState::<Sha256Config, G>::new(
            MerkleProofStrategy::Uncompressed,
            leaf_hash_params,
            two_to_one_params,
        );

        // Open at multiple points with same alpha
        let alpha: Vec<F> = (0..k).map(|_| F::random(&mut rng)).collect();
        let mut openings = Vec::new();
        for leaf_idx in 0..80 {
            openings.push((leaf_idx, alpha.clone()));
        }

        let start = Instant::now();
        let (batch_proof, evaluations) = commitment
            .open_at_field_points_batch(&openings, &prover_state)
            .expect("Failed to open");
        let duration = start.elapsed();
        println!("time to compute opening proof is {:?}", duration);

        let start = Instant::now();
        let is_valid = GroupWhirCommitment::<Sha256Config, G>::verify_openings_batch(
            &batch_proof,
            &alpha,
            &evaluations,
            &root,
            &verifier_state,
        )
        .expect("Batch verification failed");
        println!("Batch verification time {}", start.elapsed().as_millis());

        assert!(is_valid);
    }
}
