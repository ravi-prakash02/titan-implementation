use ark_ec::{CurveGroup, VariableBaseMSM};
use ark_ff::{BigInteger, PrimeField, Zero};
use ark_pallas;
use ark_pallas::Fr as ArkFr;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use pasta_curves::group::GroupEncoding;
use pasta_msm::pallas;
use std::ops::{Add, AddAssign, Mul, Sub};
use std::time::Instant;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ArkPallasPoint(pub ark_pallas::Projective);
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ArkAffine(pub ark_pallas::Affine);
use crate::pastatypes::{Scalar as F, Scalar};
use crate::traits::{ByteSerializable, InnerProduct, Linear};
use crate::utils::{ark_to_pasta_batch, group_positions_vec, pasta_to_ark_batch};
use ArkAffine as GAffine;
use ArkPallasPoint as G;
// Implement traits for ark pallas curve types
impl Add for ArkPallasPoint {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl AddAssign for ArkPallasPoint {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
    }
}

impl Mul<F> for ArkPallasPoint {
    type Output = Self;

    fn mul(self, rhs: F) -> Self::Output {
        let arkF = pasta_to_ark_batch(&[rhs]);
        Self(self.0 * arkF[0])
    }
}

impl Sub for ArkPallasPoint {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        ArkPallasPoint(self.0 - rhs.0)
    }
}

impl ArkPallasPoint {
    pub fn normalize_batch(points: &[Self]) -> Vec<ark_pallas::Affine> {
        let inner: Vec<ark_pallas::Projective> = points.iter().map(|p| p.0).collect();
        ark_pallas::Projective::normalize_batch(&inner)
    }

    pub fn msm(bases: &[Self], scalars: &[F]) -> Option<Self> {
        let scalars = pasta_to_ark_batch(scalars);
        let bases_affine = Self::normalize_batch(bases);
        Some(ArkPallasPoint(
            ark_pallas::Projective::msm(&bases_affine, &scalars).unwrap(),
        ))
    }
}

impl Mul<F> for ArkAffine {
    type Output = ArkPallasPoint;

    fn mul(self, rhs: F) -> Self::Output {
        let arkF = pasta_to_ark_batch(&[rhs]);
        ArkPallasPoint(self.0 * arkF[0])
    }
}

impl ArkAffine {
    pub fn msm(bases: &[Self], scalars: &[F]) -> Option<ArkPallasPoint> {
        let inner: Vec<ark_pallas::Affine> = bases.iter().map(|p| p.0).collect();
        let scalars = pasta_to_ark_batch(scalars);
        Some(ArkPallasPoint(
            ark_pallas::Projective::msm(&inner, &scalars).unwrap(),
        ))
    }
}

impl Linear<F> for G {
    fn zero() -> Self {
        ArkPallasPoint(ark_pallas::Projective::zero())
    }
}

impl InnerProduct<F> for G {
    type Output = G;
    type Affine = GAffine;

    fn inner_product(bases: &[Self], scalars: &[F]) -> Self::Output {
        bases
            .iter()
            .zip(scalars.iter())
            .fold(G::zero(), |acc, (x, y)| *x * *y)
    }

    fn inner_product_msm(bases: &[Self], scalars: &[F]) -> Self::Output {
        G::msm(&bases, scalars).unwrap()
    }

    fn inner_product_orbit(bases: &[Self], scalars: &[F]) -> Self::Output {
        let ark_scalars = pasta_to_ark_batch(scalars);
        let orbits = group_positions_vec(&ark_scalars);
        // compute generator sum for each orbit
        let start = Instant::now();
        let g_sums = (0..orbits.len())
            .into_iter()
            .map(|i| {
                let gsum = orbits[i]
                    .1
                    .iter()
                    .fold(Self::zero(), |acc, x| (acc + bases[*x]));
                (orbits[i].0, gsum)
            })
            .collect::<Vec<_>>();
        //println!("Time to compute gsums = {}", start.elapsed().as_micros());

        let (orbit_scalars, orbit_points): (Vec<ArkFr>, Vec<Self>) = g_sums.into_iter().unzip();
        //println!("orbit scalars len = {}", orbit_scalars.len());
        let orbit_scalars = ark_to_pasta_batch(&orbit_scalars);
        Self::inner_product_msm(&orbit_points, &orbit_scalars)
    }

    fn to_affine(bases: &[Self]) -> Vec<Self::Affine> {
        G::normalize_batch(bases)
            .iter()
            .map(|p| ArkAffine(*p))
            .collect()
    }
}

impl InnerProduct<F> for GAffine {
    type Output = G;
    type Affine = GAffine;

    fn inner_product(bases: &[Self], scalars: &[F]) -> Self::Output {
        bases
            .iter()
            .zip(scalars.iter())
            .fold(G::zero(), |acc, (x, y)| acc + *x * *y)
    }

    fn inner_product_msm(bases: &[Self], scalars: &[F]) -> Self::Output {
        GAffine::msm(bases, scalars).unwrap()
    }

    fn inner_product_orbit(bases: &[Self], scalars: &[F]) -> Self::Output {
        Self::inner_product_msm(bases, scalars)
    }

    fn to_affine(bases: &[Self]) -> Vec<Self::Affine> {
        Vec::from(bases)
    }
}

impl ByteSerializable for GAffine {
    fn serialize_to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(self.0.compressed_size());
        self.0.serialize_compressed(&mut bytes).unwrap();
        bytes
    }

    fn deserialize_from_bytes(bytes: &[u8]) -> Option<Self> {
        Some(ArkAffine(
            ark_pallas::Affine::deserialize_compressed(bytes).unwrap(),
        ))
    }
}

impl ByteSerializable for G {
    fn serialize_to_bytes(&self) -> Vec<u8> {
        ArkAffine(self.0.into_affine()).serialize_to_bytes()
    }
    fn deserialize_from_bytes(bytes: &[u8]) -> Option<Self> {
        let g_affine = GAffine::deserialize_from_bytes(bytes);
        if g_affine.is_none().into() {
            return None;
        }
        let g_proj = ark_pallas::Projective::from(g_affine.unwrap().0);
        Some(ArkPallasPoint(g_proj))
    }
}
