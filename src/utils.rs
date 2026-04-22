use crate::pastatypes::{GAffine, Point as G, Scalar as F, Scalar};
/// Helper functions and utilities for TITAN commitment scheme
use ark_pallas::Fr as ArkFr;
use ff::PrimeField as PastaPrimeField;
use ff::{FromUniformBytes, PrimeFieldBits};
use pasta_curves::group::{ff::Field, Curve, Group, GroupEncoding};
//use pasta_curves::group::Curve ;

use crate::multilinear::MultilinearPoly;
use crate::traits::{ByteSerializable, Linear};
use ark_ff::{BigInteger, FftField, PrimeField};
use ark_pallas::Fr as ArkPallasScalar; // arkworks Pallas scalar
use ark_poly::domain::radix2::Radix2EvaluationDomain;
use ark_poly::EvaluationDomain;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use std::collections::HashSet;
use std::convert::TryInto;
use std::fmt::Debug;
use std::hash::Hash;

/// Utility functions for interfacing with Arkworks field types.
pub fn as_array_32(bytes: &[u8]) -> Option<&[u8; 32]> {
    bytes.try_into().ok()
}

/// Convert from pasta_curves Pallas to arkworks Pallas
pub fn pasta_to_ark(scalar: &F) -> ArkPallasScalar {
    let bits = scalar.to_le_bits(); //this returns a bit array
    let mut bytes = vec![0u8; 64];

    // Convert bits to bytes
    for (i, bit) in bits.iter().enumerate() {
        if *bit {
            bytes[i / 8] |= 1 << (i % 8);
        }
    }

    ArkPallasScalar::from_le_bytes_mod_order(&bytes)
}

/// Convert from arkworks Pallas back to pasta_curves
pub fn ark_to_pasta(scalar: &ArkPallasScalar) -> F {
    let mut bytes = [0u8; 64];

    scalar.serialize_compressed(&mut bytes.as_mut()).unwrap();
    F::from_uniform_bytes(&bytes)
}

pub fn point_to_bytes(pt: &G) -> Vec<u8> {
    pt.serialize_to_bytes()
    //let affine_pt = pt.to_affine();
    //let buf = affine_pt.to_bytes();
    //buf.to_vec()
}

/// Batch conversion variants
pub fn ark_to_pasta_batch(input: &[ArkFr]) -> Vec<F> {
    input
        .iter()
        .map(|f| {
            // Ark -> canonical bytes
            let bigint = f.into_bigint();
            let bytes = bigint.to_bytes_le();
            assert_eq!(bytes.len(), 32);
            let bytes32 = as_array_32(&bytes).unwrap();

            // Reinterpret as Pasta repr
            let pasta_repr = <F as PastaPrimeField>::Repr::from(*bytes32);
            // Canonical parse (cannot fail)
            Scalar::from_repr(pasta_repr).unwrap()
        })
        .collect()
}

pub fn pasta_to_ark_batch(input: &[F]) -> Vec<ArkFr> {
    input
        .iter()
        .map(|s| {
            let repr = s.to_repr(); // Pasta repr (32 bytes)
            ArkFr::from_le_bytes_mod_order(repr.as_ref())
        })
        .collect()
}

/// Montgomery's batch inversion: 1 field inversion + O(n) multiplications.
pub fn batch_invert<F: ff::Field>(input: &[F]) -> Vec<F> {
    let n = input.len();
    if n == 0 { return vec![]; }

    // Prefix products
    let mut prefix = Vec::with_capacity(n);
    prefix.push(input[0]);
    for i in 1..n {
        prefix.push(prefix[i - 1] * input[i]);
    }

    // Single inversion
    let mut acc = prefix[n - 1].invert().unwrap();

    // Back-propagate
    let mut result = vec![F::ZERO; n];
    for i in (1..n).rev() {
        result[i] = acc * prefix[i - 1];
        acc *= input[i];
    }
    result[0] = acc;

    result
}

/// Create a smooth evaluation domain using arkworks
pub fn create_smooth_domain(k: usize) -> (Radix2EvaluationDomain<ArkPallasScalar>, Vec<F>) {
    let domain_size = 1 << k;
    let ark_domain = Radix2EvaluationDomain::<ArkPallasScalar>::new(domain_size)
        .expect("domain size must be power of 2");

    // Convert domain elements back to pasta_curves if needed
    //let pasta_elements: Vec<F> = ark_domain
    //    .elements()
    //    .map(|e| ark_to_pasta(&e))
    //    .collect();
    let ark_elements = ark_domain.elements().collect::<Vec<_>>();
    let pasta_elements = ark_to_pasta_batch(&ark_elements);

    (ark_domain, pasta_elements)
}

// returns the vector y, y^2,\ldots, y^{2^{l-1}}.
// reverses the vector if the flag is set.
pub fn generate_power_vec<F: ff::Field>(y: F, l: usize, reverse: bool) -> Vec<F> {
    let mut y_pow: Vec<F> = Vec::new();
    let mut yp = y;
    for _i in 0..l {
        y_pow.push(yp);
        yp = yp * yp;
    }
    if !reverse {
        y_pow
    } else {
        y_pow.reverse();
        y_pow
    }
}

/// For a domain L with generator g, the folded domain L^(2^k) has generator g^(2^k).
/// This function generates all points in the folded domain
/// The folded domain has size (ρ^{-1} · 2^m) / 2^k = ρ^{-1} · 2^(m-k) points.
/// # Arguments
/// * domain_gen - Generator of the original domain L
/// * inv_rate - Inverse rate
/// * num_variables - Number of variables m
/// * k - Folding parameter (generates L^(2^k))
/// # Returns
/// Vector of folded domain points: {1, g^(2^k), (g^(2^k))^2, ..., (g^(2^k))^(|L^(2^k)|-1)}
pub fn generate_folded_domain_points(
    domain_gen: F,
    inv_rate: usize,
    num_variables: usize,
    folding_param: usize,
) -> Vec<F> {
    // |L| = ρ^{-1} · 2^m
    let domain_size = inv_rate << num_variables;
    // |L^(2^k)| = |L| / 2^k = ρ^{-1} · 2^(m-k)
    let folded_size = domain_size >> folding_param;
    // Compute the generator of the folded domain L^(2^k)
    let folded_gen = domain_gen.pow(&[(1u64 << folding_param as u64)]);
    let mut points = Vec::with_capacity(folded_size);
    points.push(F::one());
    let mut current = folded_gen;
    for _ in 1..folded_size {
        points.push(current);
        current *= folded_gen;
    }
    points
}

/// Multilinear equality polynomial structure
/// Stores the monomial basis coefficients of eq̃(x, α)
#[derive(Debug, Clone)]
pub struct EqualityPolynomial {
    pub coeffs: Vec<F>,
    pub num_vars: usize,
}

impl EqualityPolynomial {
    pub fn new(alpha: &[F]) -> Self {
        let num_vars = alpha.len();
        let mut coeffs: Vec<F> = vec![F::one(); 1 << num_vars]; // 2^n coefficients initialized to 1

        // Compute product ∏(1 - α_i)
        let prod = alpha.iter().fold(F::one(), |acc, x| acc * (F::one() - *x));

        // Compute val[i] = α_i / (1 - α_i) for each i
        let val: Vec<F> = alpha
            .iter()
            .map(|a| *a * (F::one() - *a).invert().unwrap())
            .collect();

        for i in 1..=alpha.len() {
            let n = 1usize << i; // n = 2^i
            for j in 0..(n / 2) {
                coeffs[j + n / 2] = coeffs[j] * val[i - 1];
            }
        }

        // Multiply all coefficients by the product term
        coeffs = coeffs.iter().map(|c| *c * prod).collect();

        Self { coeffs, num_vars }
    }

    pub fn evaluate(&self, point: &[F]) -> F {
        assert_eq!(
            point.len(),
            self.num_vars,
            "Point dimension must match number of variables"
        );
        let mut result = F::zero();
        for (i, coeff) in self.coeffs.iter().enumerate() {
            let mut term = *coeff;
            for j in 0..self.num_vars {
                if (i >> j) & 1 == 1 {
                    term *= point[j];
                } else {
                    term *= F::one() - point[j];
                }
            }
            result += term;
        }
        result
    }

    /// Evaluate at all points in {0,1}^k
    ///
    /// Evaluates the polynomial at all 2^k boolean points.
    ///
    /// # Returns
    /// Vector of all 2^k equality values indexed by binary representation
    /// results[i] = eq̃(b, α) where b is the binary representation of i
    pub fn evaluate_all(&self) -> Vec<F> {
        let size = 1 << self.num_vars;
        (0..size)
            .map(|i| {
                let mut point = vec![F::zero(); self.num_vars];
                for j in 0..self.num_vars {
                    if (i >> j) & 1 == 1 {
                        point[j] = F::one();
                    }
                }
                self.evaluate(&point)
            })
            .collect()
    }
}

/// Batch compute eq̃(α, b) for all b ∈ {0,1}^k (convenience function)
///
/// Creates an equality polynomial and evaluates it at all boolean points.
///
/// # Arguments
/// * `alpha` - Challenge vector of scalars
///
/// # Returns
/// Vector of all 2^k equality values
pub fn eq_all(alpha: &[F]) -> Vec<F> {
    let poly = MultilinearPoly::init_with_eq(alpha);
    poly.coeffs
}

/// Compute eq̃(α, b) from a binary index
///
/// Convenience function that creates the equality polynomial and evaluates at a single point.
///
/// # Arguments
/// * `alpha` - Challenge vector of scalars
/// * `b_idx` - Binary index (will be converted to binary vector of length alpha.len())
///
/// # Returns
/// The value of eq̃(α, b) where b is the binary representation of b_idx
pub fn eq_from_index(alpha: &[F], b_idx: usize) -> F {
    let k = alpha.len();
    let mut b = vec![F::zero(); k];

    for i in 0..k {
        if (b_idx >> i) & 1 == 1 {
            b[i] = F::one();
        }
    }

    // compute eq(alpha,b) in alpha.len() multiplications.
    let mut K = F::one();
    for i in 0..k {
        K = K * (b[i] * alpha[i] + (F::one() - alpha[i]) * (F::one() - b[i]));
    }
    K
}

/// Reverse the least significant `num_bits` bits of `x`.
fn reverse_bits(x: usize, num_bits: usize) -> usize {
    let mut result = 0;
    let mut x = x;
    for _ in 0..num_bits {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    result
}

/// Bit-reverse permutation: reorder so that element at index `i` moves to
/// index `reverse_bits(i, m)`. Applied to multilinear polynomial coefficients
/// indexed as `coeffs[b_1 + 2·b_2 + … + 2^{m-1}·b_m]`, this maps the
/// polynomial f(x_1, …, x_m) to f(x_m, …, x_1).
pub fn bit_reverse_permutation<T>(data: &mut [T], m: usize) {
    let n = 1usize << m;
    assert_eq!(data.len(), n);
    for i in 0..n {
        let j = reverse_bits(i, m);
        if i < j {
            data.swap(i, j);
        }
    }
}

/// Compute fft evaluations f(x,x^2,..x^{2^{m-1}}) over domain D
/// Given the m variate multilinear polynomial f in evaluation basis
/// The algorithm works in O(m|D|) operations
pub fn multilinear_fft<G, F>(f_coeffs: &[G], domain: &[F], m: usize, d: usize) -> Vec<G>
where
    F: ff::Field,
    G: Linear<F>,
{
    assert_eq!(
        f_coeffs.len(),
        1usize << m,
        "Polynomial size does not match m"
    );
    assert_eq!(domain.len(), 1usize << d, "Domain size does not match d");
    assert!(d >= m, "Domain must be larger than polynomial");

    // Bit-reverse the coefficients before the butterfly. The butterfly processes
    // variables from bit m-1 down to bit 0, but the coefficient index encodes
    // x_1 at bit 0 (LSB). After bit-reversal, bit 0 holds x_m, so the butterfly
    // processes x_m first with the highest power y^{2^{m-1}}, yielding the
    // standard evaluation f(y, y^2, …, y^{2^{m-1}}).
    let mut permuted = f_coeffs.to_vec();
    bit_reverse_permutation(&mut permuted, m);

    let mut init_vec: Vec<G> = Vec::new();
    for i in 0..permuted.len() {
        for _ in 0..(1usize << (d - m)) {
            init_vec.push(permuted[i]);
        }
    }
    let D = 1usize << d;
    let mut new_vec: Vec<G> = vec![G::zero(); D];

    for p in 0..m {
        let chunk_size = 1usize << (d - (m - 1 - p));
        let step_size = chunk_size >> 1;
        for chunk_start in (0..D).step_by(chunk_size) {
            for i in 0..step_size {
                let rootidx: usize = ((1usize << (m-1-p))*i) % D;
                let root = domain[rootidx];
                let factor =
                    (init_vec[chunk_start + step_size + i] - init_vec[chunk_start + i]) * root;
                new_vec[chunk_start + i] = init_vec[chunk_start + i] + factor;
                new_vec[chunk_start + step_size + i] = init_vec[chunk_start + i] - factor;
            }
        }
        std::mem::swap(&mut init_vec, &mut new_vec);

    }

    init_vec
}


pub fn group_positions_vec<T: Eq + Hash + Clone>(
    v: &[T],
) -> Vec<(T,  Vec<usize>)> {
    let mut map = std::collections::HashMap::new();

    for (i, x) in v.iter().enumerate() {
        map.entry(x.clone())
            .or_insert_with(Vec::new)
            .push(i);
    }

    map.into_iter().collect()
}




#[cfg(test)]
mod tests {
    use std::time::Instant;
    use ark_std::iterable::Iterable;
    use ark_std::UniformRand;
    use super::*;
    use rand::prelude::StdRng;
    use rand::SeedableRng;
    use crate::traits::InnerProduct;

    #[test]
    fn test_eq_function() {
        let mut rng = StdRng::from_entropy();
        let k = 3;
        let alpha: Vec<F> = (0..k).map(|_| F::random(&mut rng)).collect();

        let poly = EqualityPolynomial::new(&alpha);

        // Evaluate at a specific binary point
        let b = vec![F::one(), F::zero(), F::one()];
        let eq_val = poly.evaluate(&b);

        // Manually compute to verify (using the definition)
        let mut expected = F::one();
        for i in 0..k {
            expected *= b[i] * alpha[i] + (F::one() - b[i]) * (F::one() - alpha[i]);
        }

        assert_eq!(eq_val, expected);
    }

    #[test]
    fn test_eq_from_index_consistency() {
        let mut rng = StdRng::from_entropy();
        let k = 4;
        let alpha: Vec<F> = (0..k).map(|_| F::random(&mut rng)).collect();

        let poly = EqualityPolynomial::new(&alpha);

        // Test that eq_from_index and evaluate produce same results
        for b_idx in 0..(1 << k) {
            let mut b = vec![F::zero(); k];
            for i in 0..k {
                if (b_idx >> i) & 1 == 1 {
                    b[i] = F::one();
                }
            }

            let eq_direct = poly.evaluate(&b);
            let eq_from_idx = eq_from_index(&alpha, b_idx);

            assert_eq!(eq_direct, eq_from_idx);
        }
    }

    #[test]
    fn test_eq_all_consistency() {
        let mut rng = StdRng::from_entropy();
        let k = 3;
        let alpha: Vec<F> = (0..k).map(|_| F::random(&mut rng)).collect();

        let poly = EqualityPolynomial::new(&alpha);
        let eq_all_vals = poly.evaluate_all();

        for b_idx in 0..(1 << k) {
            let mut b = vec![F::zero(); k];
            for i in 0..k {
                if (b_idx >> i) & 1 == 1 {
                    b[i] = F::one();
                }
            }

            let eq_individual = poly.evaluate(&b);
            assert_eq!(eq_all_vals[b_idx], eq_individual);
        }
    }

    #[test]
    fn test_folded_domain_points_generation() {
        let inv_rate = 4;
        let num_variables = 4;
        let k = 2;

        let domain_size = inv_rate << num_variables;
        let folded_size = domain_size >> k;

        assert_eq!(domain_size, 64);
        assert_eq!(folded_size, 16);
    }

    #[test]
    fn test_ark_pasta_generation() {
        let k: usize = 3;
        let domain = create_smooth_domain(k);
        // domain should be of the type: 1, a, a^2,...,a^{size-1}
        for i in 0..(domain.1.len() - 1) {
            assert_eq!(domain.1[1] * domain.1[i], domain.1[i + 1]);
        }

        // check consistency between batch conversion and single conversion.
        assert_eq!(ark_to_pasta(&domain.0.element(1)), domain.1[1]);
    }

    #[test]
    fn test_multilinear_fft() {
        let k: usize = 10;
        let m: usize = 7;
        let mut rng = StdRng::from_entropy();

        let domain = create_smooth_domain(k);
        let evals:Vec<G> = (0..(1usize << m)).into_iter().map(|_| G::random(&mut rng)).collect();
        let poly = MultilinearPoly::new(evals.clone());
        let start = Instant::now();
        let fft_evals = multilinear_fft(&evals, &domain.1, m, k);
        println!("Time for multilinear fft = {} millis", start.elapsed().as_millis());
    }

    #[test]
    fn test_orbit_msm() {
        let n: usize = 12;
        let m: usize = 5;
        let N: usize = 1usize << n;
        let M: usize = 1usize << m;
        let mut rng = StdRng::from_entropy();
        // Sample random generators
        let gens: Vec<G> = (0..N).into_iter().map(|_| G::random(&mut rng)).collect();
        let table: Vec<F> = (0..M).into_iter().map(|_| F::random(&mut rng)).collect();
        let scalars: Vec<F> = (0..N).into_iter().map(|_| table[usize::rand(&mut rng) % M]).collect();

        let start = Instant::now();
        let res1 = G::inner_product_msm(&gens, &scalars);
        println!("Time for generic msm {}", start.elapsed().as_micros());
        let start = Instant::now();
        let res2 = G::inner_product_orbit(&gens, &scalars);
        println!("Time for orbit msm {}", start.elapsed().as_micros());

        assert_eq!(res1, res2, "Results don't match");
    }
}
