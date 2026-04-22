use std::borrow::Borrow;

use ark_crypto_primitives::{
    merkle_tree::{Config, LeafParam, Path, TwoToOneParam},
    Error,
};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

/// Wrapper for multiple Merkle proofs (paths) for batch verification.
///
/// This struct stores a collection of individual `Path` proofs and provides
/// batch verification functionality. Useful for verifying multiple leaf openings
/// from the same Merkle tree in a single operation.
#[derive(Clone, Default, CanonicalSerialize, CanonicalDeserialize)]
pub struct FullMultiPath<P: Config> {
    pub proofs: Vec<Path<P>>,
}

impl<P: Config> FullMultiPath<P> {
    /// Extract the leaf indices from all proofs
    pub fn indices(&self) -> Vec<usize> {
        self.proofs.iter().map(|p| p.leaf_index).collect()
    }

    /// Verify all proofs against the given root
    ///
    /// # Arguments
    /// * `leaf_hash_params` - Parameters for leaf hashing
    /// * `two_to_one_params` - Parameters for internal node hashing
    /// * `root_hash` - The expected Merkle tree root
    /// * `leaves` - Iterator of leaf values to verify
    ///
    /// # Returns
    /// `Ok(true)` if all proofs are valid, `Ok(false)` if any proof is invalid,
    /// or `Err` if verification encounters an error.
    pub fn verify<L: Borrow<P::Leaf> + Clone>(
        &self,
        leaf_hash_params: &LeafParam<P>,
        two_to_one_params: &TwoToOneParam<P>,
        root_hash: &P::InnerDigest,
        leaves: impl IntoIterator<Item = L>,
    ) -> Result<bool, Error> {
        let leaves_vec: Vec<L> = leaves.into_iter().collect();
        self.proofs
            .iter()
            .enumerate()
            .map(|(i, proof)| {
                proof.verify(
                    leaf_hash_params,
                    two_to_one_params,
                    root_hash,
                    leaves_vec[i].borrow(),
                )
            })
            .try_fold(true, |acc, res| res.map(|b| acc && b))
    }
}

impl<P: Config> From<Vec<Path<P>>> for FullMultiPath<P> {
    fn from(proofs: Vec<Path<P>>) -> Self {
        Self { proofs }
    }
}