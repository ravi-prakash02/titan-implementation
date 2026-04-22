use std::hash::Hash;
// Define key traits
use std::ops::{Add, AddAssign, Mul, Sub};

// trait to capture F modules
pub trait Linear<F>: Copy + Clone + Add<Output = Self> + AddAssign + Sub<Output=Self> + Mul<F, Output = Self> {
    fn zero() -> Self;
}

// trait to capture inner product with scalars
pub trait InnerProduct<F, Affine=F>: Sized {
    type Output;
    type Affine;
    fn inner_product(bases: &[Self], scalars: &[F]) -> Self::Output;
    fn inner_product_msm(bases: &[Self], scalars: &[F]) -> Self::Output;
    fn inner_product_orbit(bases:&[Self], scalars: &[F]) -> Self::Output;
    fn to_affine(bases: &[Self]) -> Vec<Self::Affine>;
}

// trait to capture serialization and deserialization to bytes
pub trait ByteSerializable: Sized {
    fn serialize_to_bytes(&self) -> Vec<u8>;
    fn deserialize_from_bytes(bytes: &[u8]) -> Option<Self>;
}
