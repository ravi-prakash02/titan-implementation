use crate::traits::*;
use ark_ff::{BigInteger, PrimeField as ArkPrimeField};
use ark_pallas::{Fr as ArkFr};
use ff::PrimeField as PastaPrimeField;
use pasta_curves::group::ff::Field;
use pasta_curves::group::prime::PrimeCurveAffine;
use pasta_curves::group::{Curve, Group, GroupEncoding};

use pasta_msm::pallas;
use std::convert::TryInto;
use std::hash::Hash;
use std::ops::{Add, AddAssign};
use std::time::Instant;
use crate::utils::{ark_to_pasta_batch, as_array_32, group_positions_vec, pasta_to_ark_batch};
// Implement traits for pasta curve types
pub type Point =  pasta_curves::pallas::Point;
pub type GAffine = pasta_curves::pallas::Affine;
pub type Scalar = pasta_curves::pallas::Scalar;

use Point as G;
use Scalar as F;

impl Linear<F> for F {
    fn zero() -> Self {
        F::zero()
    }
}

impl Linear<F> for G {
    fn zero() -> Self {
        G::identity()
    }
}

impl InnerProduct<F> for F {
    type Output = F;
    type Affine = F;

    fn inner_product(bases: &[Self], scalars: &[F]) -> Self::Output {
        bases.iter().zip(scalars.iter()).map(|(x, y)| *x * *y).sum()
    }

    fn inner_product_msm(bases: &[Self], scalars: &[F]) -> Self::Output {
        Self::inner_product(bases, scalars)
    }

    fn inner_product_orbit(bases: &[Self], scalars: &[Scalar]) -> Self::Output {
        Self::inner_product_msm(bases, scalars)
    }

    fn to_affine(bases: &[Self]) -> Vec<Self::Affine> {
        Vec::from(bases)
    }
}

impl InnerProduct<F> for G {
    type Output = G;
    type Affine = GAffine;

    fn inner_product(bases: &[Self], scalars: &[F]) -> Self::Output {
        bases.iter().zip(scalars.iter()).map(|(x, y)| *x * *y).sum()
    }

    fn inner_product_msm(bases: &[Self], scalars: &[F]) -> Self::Output {
        let mut bases_affine = vec![GAffine::identity(); bases.len()];
        G::batch_normalize(bases, &mut bases_affine);
        pallas(&bases_affine, scalars)
    }

    fn inner_product_orbit(bases: &[Self], scalars: &[Scalar]) -> Self::Output {
        let ark_scalars = pasta_to_ark_batch(scalars);
        let orbits = group_positions_vec(&ark_scalars);
        // compute generator sum for each orbit
        let start = Instant::now();
        let g_sums = (0..orbits.len()).into_iter().map(|i| {
            let gsum = orbits[i].1.iter().fold(Self::identity(), |acc, x| (acc + bases[*x]));
            (orbits[i].0,gsum)
        }).collect::<Vec<_>>();
        //println!("Time to compute gsums = {}", start.elapsed().as_micros());

        let (orbit_scalars, orbit_points): (Vec<ArkFr>, Vec<Self>) = g_sums.into_iter().unzip();
        //println!("orbit scalars len = {}", orbit_scalars.len());
        let orbit_scalars = ark_to_pasta_batch(&orbit_scalars);
        Self::inner_product_msm(&orbit_points, &orbit_scalars)
    }


    fn to_affine(bases: &[Self]) -> Vec<Self::Affine> {
        let mut bases_affine: Vec<Self::Affine> = vec![GAffine::identity(); bases.len()];
        G::batch_normalize(bases, &mut bases_affine);
        bases_affine
    }
}

impl InnerProduct<F> for GAffine {
    type Output = G;
    type Affine = GAffine;

    fn inner_product(bases: &[Self], scalars: &[F]) -> Self::Output {
        bases.iter().zip(scalars.iter()).map(|(x, y)| *x * *y).sum()
    }

    fn inner_product_msm(bases: &[Self], scalars: &[F]) -> Self::Output {
        pallas(bases, scalars)
    }

    fn inner_product_orbit(bases: &[Self], scalars: &[Scalar]) -> Self::Output {
        let ark_scalars = pasta_to_ark_batch(scalars);
        let orbits = group_positions_vec(&ark_scalars);
        // compute generator sum for each orbit
        let start = Instant::now();
        let g_sums = (0..orbits.len()).into_iter().map(|i| {
            let gsum = orbits[i].1.iter().fold(Self::identity(), |acc, x| (acc + bases[*x]).into());
            (orbits[i].0,gsum)
        }).collect::<Vec<_>>();
        println!("Time to compute gsums = {}", start.elapsed().as_micros());

        let (orbit_scalars, orbit_points): (Vec<ArkFr>, Vec<Self>) = g_sums.into_iter().unzip();
        println!("orbit scalars len = {}", orbit_scalars.len());
        let orbit_scalars = ark_to_pasta_batch(&orbit_scalars);
        Self::inner_product_msm(&orbit_points, &orbit_scalars)
    }

    fn to_affine(bases: &[Self]) -> Vec<Self::Affine> {
        Vec::from(bases)
    }
}

impl ByteSerializable for F {
    fn serialize_to_bytes(&self) -> Vec<u8> {
        self.to_repr().to_vec()
    }

    fn deserialize_from_bytes(bytes: &[u8]) -> Option<Self> {

        let mut limbs = [0u64; 4];
        for i in 0..4 {
            limbs[i] = u64::from_le_bytes(bytes[i * 8..(i + 1) * 8].try_into().unwrap());
        }
        Some(F::from_raw(limbs))

        /*
        let repr = <F as PastaPrimeField>::Repr::try_from(bytes).ok()?;
        let fq = Scalar::from_repr(repr);
        return if fq.is_none().into() {
            None
        } else {
            Some(fq.unwrap())
        };

         */
    }
}

impl ByteSerializable for GAffine {
    fn serialize_to_bytes(&self) -> Vec<u8> {
        self.to_bytes().to_vec()
    }

    fn deserialize_from_bytes(bytes: &[u8]) -> Option<Self> {
        let repr = <GAffine as GroupEncoding>::Repr::try_from(bytes).ok()?;
        let g_affine = GAffine::from_bytes(&repr);
        if g_affine.is_none().into() {
            return None;
        }
        Some(g_affine.unwrap())
    }
}

impl ByteSerializable for G {
    fn serialize_to_bytes(&self) -> Vec<u8> {
        self.to_affine().to_bytes().to_vec()
    }

    fn deserialize_from_bytes(bytes: &[u8]) -> Option<Self> {
        let repr = <GAffine as GroupEncoding>::Repr::try_from(bytes).ok()?;
        let g_affine = GAffine::from_bytes(&repr);
        if g_affine.is_none().into() {
            return None;
        }
        Some(G::from(g_affine.unwrap()))
    }
}

