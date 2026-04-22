use std::sync::atomic::{AtomicUsize, Ordering};

pub mod blake3;
pub mod digest;
pub mod keccak;
pub mod parameters;
pub mod sha256;
pub mod proof;

// Re-export commonly used types
pub use blake3::{Blake3Compress, Blake3Digest, Blake3LeafHash, Blake3MerkleTreeParams};
pub use digest::GenericDigest;
pub use keccak::{KeccakCompress, KeccakDigest, KeccakLeafHash, KeccakMerkleTreeParams};
pub use parameters::{default_config, IdentityDigestConverter, MerkleTreeParams};
pub use sha256::{Sha256Compress, Sha256Digest, Sha256LeafHash, Sha256MerkleTreeParams};

/// Hash counter utility for benchmarking and profiling hash operations.
#[derive(Debug, Default)]
pub struct HashCounter;

static HASH_COUNTER: AtomicUsize = AtomicUsize::new(0);

impl HashCounter {
    /// Increment the hash counter
    pub(crate) fn add() -> usize {
        HASH_COUNTER.fetch_add(1, Ordering::SeqCst)
    }

    /// Reset the hash counter to zero
    pub fn reset() {
        HASH_COUNTER.store(0, Ordering::SeqCst);
    }

    /// Get the current hash counter value
    pub fn get() -> usize {
        HASH_COUNTER.load(Ordering::SeqCst)
    }
}