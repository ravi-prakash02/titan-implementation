use std::{borrow::Borrow, marker::PhantomData};
use ark_crypto_primitives::{
    crh::{CRHScheme, TwoToOneCRHScheme},
    Error,
};

use pasta_curves::pallas::{Scalar as F, Point as G } ;
use pasta_curves::group::{ff::Field, Curve, Group, GroupEncoding} ;
use pasta_curves::pallas::Affine as GAffine;
use ark_serialize::CanonicalSerialize;
use rand::RngCore;
use serde::{Deserialize, Serialize};
use sha2::{Sha256, Digest as Sha2Digest};
use crate::traits::ByteSerializable;
use super::{digest::GenericDigest, parameters::MerkleTreeParams, HashCounter};
use crate::utils::{ark_to_pasta, pasta_to_ark, point_to_bytes};

/// Digest type used in SHA256-based TITAN Merkle trees.
///
/// Alias for a 32-byte generic digest.
pub type Sha256Digest = GenericDigest<32>;

/// Merkle tree configuration using SHA256 as both leaf and node hasher for TITAN.
pub type Sha256MerkleTreeParams<C> =
    MerkleTreeParams<C, Sha256LeafHash<C>, Sha256Compress, Sha256Digest>;

/// Leaf hash function using SHA256 over compressed `[G]` input (group elements).
///
/// This struct implements `CRHScheme` where the input is a slice of
/// canonical-serializable group elements `[G]`, and the output is a
/// 32-byte SHA256 digest.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct Sha256LeafHash<C>(#[serde(skip)] PhantomData<C>);

impl<C: Send + ByteSerializable> CRHScheme for Sha256LeafHash<C> {
    type Input = Vec<C>;
    type Output = Sha256Digest;
    type Parameters = ();

    fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        (): &Self::Parameters,
        input: T,
    ) -> Result<Self::Output, Error> {

        let mut buf = vec![] ;


        //this is an inefficient way to get each points and convert them to affine points
        //and then bytes

        //.............................TODO.....................//
        //better efficiency
        for pt in input.borrow() {
               buf.extend(pt.serialize_to_bytes());
        }
        

        //input.borrow()(&mut buf)?;
        
        let mut hasher = Sha256::new();
        hasher.update(&buf);
        let output: [_; 32] = hasher.finalize().into();
        
        HashCounter::add();
        Ok(output.into())
    }
}

/// Node compression function using SHA256 over two 32-byte digests.
///
/// This struct implements `TwoToOneCRHScheme`, combining two digests
/// by concatenation and hashing with SHA256.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Sha256Compress;

impl TwoToOneCRHScheme for Sha256Compress {
    type Input = Sha256Digest;
    type Output = Sha256Digest;
    type Parameters = ();

    fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        (): &Self::Parameters,
        left_input: T,
        right_input: T,
    ) -> Result<Self::Output, Error> {
        let mut hasher = Sha256::new();
        hasher.update(&left_input.borrow().0);
        hasher.update(&right_input.borrow().0);
        let output: [_; 32] = hasher.finalize().into();
        
        HashCounter::add();
        Ok(output.into())
    }

    fn compress<T: Borrow<Self::Output>>(
        parameters: &Self::Parameters,
        left_input: T,
        right_input: T,
    ) -> Result<Self::Output, Error> {
        Self::evaluate(parameters, left_input, right_input)
    }
}