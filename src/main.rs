#![allow(non_snake_case)]
#![allow(unused_imports)]
#![allow(dead_code)]
mod sumcheck;
mod multilinear;
mod group_sumcheck;
mod traits;
mod pastatypes;
mod titantranscript;
mod group_whir_committer;
mod utils;
mod merkle_tree;
mod dory;
mod protocols;
mod arkpallastypes;
mod r1cs;
mod spartan;

use anyhow::Result;
use pasta_curves::group::{Curve, Group};
use pasta_curves::group::ff::FromUniformBytes;
use rand::{rngs::StdRng, RngCore, SeedableRng};
use ark_poly::{univariate::DensePolynomial, GeneralEvaluationDomain};
use crate::pastatypes::{Point, GAffine, Scalar};
// Re-exported types from pasta_curves

//use pasta_msm::{pallas as pasta_msm_api, pallas};

pub fn rand_scalar(rng: &mut StdRng) -> Scalar {
    // Construct a scalar from 32 random bytes (reduce modulo group order)
    let mut b = [0u8; 64];
    rng.fill_bytes(&mut b);
    Scalar::from_uniform_bytes(&b)
}

pub fn rand_point(rng: &mut StdRng) -> Point {
    // Use the generator multiplied by a random scalar to get a random point.
    let s = rand_scalar(rng);
    Point::generator() * &s
}

pub fn naive_msm(points: &[Point], scalars: &[Scalar]) -> Point {
    assert_eq!(points.len(), scalars.len());
    let mut acc = Point::identity();
    for (P, s) in points.iter().zip(scalars.iter()) {
        acc += P * s;
    }
    acc
}


fn main() {
    println!("Hello, world!");
    // Example: compute MSM of random points/scalars and compare with pasta-msm result
    let n = 64*64usize; // works for small tests; pasta-msm shines on large n
    let mut rng = StdRng::from_entropy();


    // generate inputs
    let mut points = Vec::with_capacity(n);
    let mut scalars = Vec::with_capacity(n);
    for _ in 0..n {
        points.push(rand_point(&mut rng));
        scalars.push(rand_scalar(&mut rng));
    }

    // naive reference (single-threaded)
    let  start = std::time::Instant::now();
    let ref_res = naive_msm(&points, &scalars);
    println!("Naive MSM took {:?}", start.elapsed());

    // convert points/scalars to the representation expected by pasta-msm
    // NOTE: the pasta-msm API may require Affine points or a specific point type; adapt if necessary.
    let  affine_points: Vec<GAffine> = points.iter().map(|p| p.to_affine()).collect();
    //Curve::batch_normalize(&points, &mut affine_points);
    // Call into the pasta-msm crate. This is a best-effort example; check the actual
    // function name in the crate (e.g. `pasta_msm::msm::msm`, `pasta_msm::multi_msm`, etc.).
    let  start = std::time::Instant::now();
    //let msm_res = pallas(&affine_points, &scalars);
    println!("MSM took {:?}", start.elapsed());

    // compare
    //assert_eq!(ref_res.to_affine(), msm_res.to_affine());
    println!("MSM result matches reference for n = {}", n);


}
