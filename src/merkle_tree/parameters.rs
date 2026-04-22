// ============================================================================
// parameters.rs
// ============================================================================

use std::{hash::Hash, marker::PhantomData};

use ark_crypto_primitives::{
    crh::{CRHScheme, TwoToOneCRHScheme},
    merkle_tree::{Config, DigestConverter},
    Error, sponge::Absorb
};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

/// A trivial converter where digest of previous layer's hash is the same as next layer's input.
pub struct IdentityDigestConverter<T> {
    _prev_layer_digest: T,
}

impl<T> DigestConverter<T, T> for IdentityDigestConverter<T> {
    type TargetType = T;
    fn convert(item: T) -> Result<T, Error> {
        Ok(item)
    }
}

/// A generic Merkle tree config usable across hash types (e.g., SHA256, Blake3, Keccak).
///
/// # Type Parameters:
/// - `G`: Group element used in the leaves (TITAN specific)
/// - `LeafH`: Leaf hash function
/// - `CompressH`: Internal node hasher
/// - `Digest`: Digest type
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]

pub struct MerkleTreeParams<G, LeafH, CompressH, Digest> {
    _marker: PhantomData<(G, LeafH, CompressH, Digest)>,
}

impl<G, LeafH, CompressH, Digest> Config for MerkleTreeParams<G, LeafH, CompressH, Digest>
where
    G:  Send,
    LeafH: CRHScheme<Input = Vec<G>, Output = Digest>,
    CompressH: TwoToOneCRHScheme<Input = Digest, Output = Digest>,
    Digest: Clone
        + std::fmt::Debug
        + Default
        + CanonicalSerialize
        + CanonicalDeserialize
        + Eq
        + PartialEq
        + Hash
        + Send
        + Absorb,
{
    type Leaf = Vec<G>;

    type LeafDigest = Digest;
    type LeafInnerDigestConverter = IdentityDigestConverter<Digest>;
    type InnerDigest = Digest;

    type LeafHash = LeafH;
    type TwoToOneHash = CompressH;
}

/// Returns the `(leaf_hash_params, two_to_one_hash_params)` for any compatible Merkle tree.
///
/// # Type Parameters
/// - `G`: The leaf group element type
/// - `LeafH`: The leaf hash function
/// - `CompressH`: The two-to-one internal hash function
///
/// # Panics
/// Panics if `setup()` fails (which should not happen for deterministic hashers).
pub fn default_config<G, LeafH, CompressH>(
    rng: &mut impl rand::RngCore,
) -> (
    <LeafH as CRHScheme>::Parameters,
    <CompressH as TwoToOneCRHScheme>::Parameters,
)
where
    G:  Send,
    LeafH: CRHScheme<Input = Vec<G>> + Send,
    CompressH: TwoToOneCRHScheme + Send,
{
    (
        LeafH::setup(rng).expect("Failed to setup Leaf hash"),
        CompressH::setup(rng).expect("Failed to setup Compress hash"),
    )
}