use crate::multilinear::MultilinearPoly;
use crate::traits::{ByteSerializable, InnerProduct, Linear};
use ff::{Field, PrimeField};
use merlin::Transcript;
use pasta_curves::group::prime::PrimeCurveAffine;
use pasta_curves::group::{Curve, Group, GroupEncoding};
use pasta_curves::pallas::{Affine as PallasAffine, Point as PallasPoint, Scalar as PallasScalar};
use pasta_msm::pallas;
use std::ops::AddAssign;
use std::time::Instant;
use ark_std::{end_timer, start_timer};
use ff::BatchInvert;
use crate::titantranscript;
use crate::titantranscript::{transcript_append_point, transcript_append_scalar, transcript_challenge_scalar};
use crate::utils::batch_invert;

pub fn compute_Sl_poly<C, F>(
    m: usize,
    l: usize,
    f_poly: &MultilinearPoly<C>,
    alpha: &[F],
) -> MultilinearPoly<C>
where
    C: Linear<F> + InnerProduct<F, Output = C>,
    F: ff::Field,
{
    assert!(l <= m);
    assert_eq!(f_poly.coeffs.len(), 1usize << m);
    assert_eq!(alpha.len(), m);

    let sl_poly_size = 1usize << l;
    let slice_size = 1usize << (m - l);

    let eq_poly = MultilinearPoly::init_with_eq(alpha);
    let mut eq_poly_transpose: Vec<F> = vec![F::ZERO; 1 << m];
    let mut f_poly_transpose: Vec<C> = vec![C::zero(); 1 << m];

    // To compute sl_poly[b] for each b, we need MSM over slices of size slice_size.
    // It will be convenient to first compute transpose of f_poly and eq_poly so that slices occur
    // as continuous indices in one dimensional representation.
    for b in 0..sl_poly_size {
        for x in 0..slice_size {
            eq_poly_transpose[b * slice_size + x] = eq_poly.coeffs[x * sl_poly_size + b];
            f_poly_transpose[b * slice_size + x] = f_poly.coeffs[x * sl_poly_size + b];
        }
    }

    let mut sl_poly_coeffs: Vec<C> = vec![C::zero(); sl_poly_size];
    for b in 0..sl_poly_size {
        sl_poly_coeffs[b] = C::inner_product_msm(
            &f_poly_transpose[b * slice_size..(b + 1) * slice_size],
            &eq_poly_transpose[b * slice_size..(b + 1) * slice_size],
        );
    }
    MultilinearPoly::new(sl_poly_coeffs)
}

/// From S_l compute all S tables: [S_l, S_{l-1}, ..., S_1]
pub fn compute_all_S_tables<C, F>(S_l: Vec<C>, l: usize) -> Vec<Vec<C>>
where C: Linear<F> + InnerProduct<F, Output=C>, F: ff::Field,
{
    let mut tables = Vec::with_capacity(l);
    tables.push(S_l);

    let mut current = tables[0].clone();

    for k in (1..=l - 1).rev() {
        let stride = 1usize << k; // size of lower half
        let next: Vec<C> = current[..stride]
            .iter()
            .zip(current[stride..].iter())
            .map(|(x, y)| *x + *y)
            .collect();
        tables.push(next.clone());
        current = next;
    }

    tables
}

/// Serialize a PallasPoint into bytes deterministically.
/// `to_affine().to_compressed()` or `to_affine().to_bytes()`.
fn point_to_bytes(p: &PallasPoint) -> Vec<u8> {
    p.serialize_to_bytes()
}

/// Serialize a scalar deterministically.
fn scalar_to_bytes(s: &PallasScalar) -> Vec<u8> {
    s.serialize_to_bytes()
}

pub fn eval_triple_at_alpha<C, F>(s: &[C; 3], alpha: F) -> C
where C: Linear<F>, F: ff::Field,
{
    let two = F::ONE + F::ONE;
    let two_inv = two.invert().unwrap();
    let one = F::ONE;
    let f1 = (alpha - one)*(alpha - two)*two_inv;
    let f2 = (two - alpha)*alpha;
    let f3 = (alpha - one)*alpha*two_inv;
    s[0] * f1 + s[1] * f2 + s[2] * f3
}


/// Compute g_i(u) values for given i and u_values for rounds 1\leq i\leq \ell
/// S_i from S_tables (S_tables[0]=S_l, S_tables[l-i] = S_i)
/// Folklore sum-check is used for rounds i > \ell.
pub fn compute_gi_values<C, F>(
    m: usize,
    i: usize,
    l: usize,
    alpha: &[F],
    rho_prefix: &[F],
    S_tables: &Vec<Vec<C>>,
    u_values: &[F],
) -> Vec<C> where
C: Linear<F> + InnerProduct<F, Output = C>, F: ff::Field,
 {
    assert!(i >= 1 && i <= l);
    assert_eq!(alpha.len(), m);
    assert_eq!(rho_prefix.len(), i - 1);

    let Si = &S_tables[l - i]; //because S_i is at l-i th position
    assert_eq!(Si.len(), 1usize << i);
    let beta_i = &alpha[..(i-1)];   // this corresponds to vector \alpha_{i-1} in write-up.
    let chunk_size = 1usize << (i-1);

    /*
     * This is an optimized computation for calculation in Equation 5 of group sum-check
     * write-up. We isolate u dependent parts. This allows us to save on MSMs.
     * We compute two MSMs of size 2^{i-1} which are as:
     * 1. H0 = \sum_{b\in \{0,1\}^{i-1}} eq(r,b)/eq(alpha_{i-1},b) S_i(b,0)
     * 2. H1 = \sum_{b\in \{0,1\}^{i-1}} eq(r,b)/eq(alpha_{i-1},b) S_i(b,1)
     * 3. Global factor K = eq(r,u,alpha_i)
     * 4. Then g_i(u) = K * (1-u)/(1-\alpha_i) * H0 + K * u/alpha_i * H1
     */

    // Parts not dependent on u
    let eq_r_b_poly = MultilinearPoly::init_with_eq(rho_prefix);
    let eq_alpha_b_poly = MultilinearPoly::init_with_eq(beta_i);
    let numerators = eq_r_b_poly.coeffs.clone();
    let denominators = batch_invert::<F>(&eq_alpha_b_poly.coeffs);
    let scalars: Vec<F> = numerators.iter().zip(denominators.iter()).map(|(n, d)| *n * *d).collect();
    let H0 = C::inner_product_msm(&Si[..chunk_size], &scalars);
    let H1 = C::inner_product_msm(&Si[chunk_size..], &scalars);

    let mut g_i_u = vec![C::zero(); u_values.len()];

    for j in 0..u_values.len() {
        let u = u_values[j];
        // compute global multiplier K for u
        let mut K = u * alpha[i-1] + (F::ONE - u) * (F::ONE - alpha[i-1]);
        for k in 0..i-1 {
            let z = rho_prefix[k];
            K = K * (z * alpha[k] + (F::ONE - z) * (F::ONE - alpha[k]));
        }
        let mut left_factor = (F::ONE - u)*((F::ONE - alpha[i-1]).invert().unwrap());
        let mut right_factor = u * alpha[i-1].invert().unwrap();
        left_factor = left_factor * K;
        right_factor = right_factor * K;
        g_i_u[j] = H0 * left_factor + H1 * right_factor;
    }
    g_i_u
}

pub fn run_prover_noninteractive(
    m: usize,
    l: usize,
    f_table: &[PallasPoint],
    alpha: &[PallasScalar],
    sigma: PallasPoint,
) -> (Vec<PallasScalar>, Vec<[PallasPoint; 3]>)
{
    // Precompute S tables
    //let S_l = compute_S_l(m, l, f_table, alpha, &fast_msm);
    let f_poly = MultilinearPoly::new(Vec::from(f_table));

    let precompute_time = Instant::now();
    let S_l = compute_Sl_poly(m, l, &f_poly, alpha);
    let S_tables = compute_all_S_tables(S_l.coeffs, l);
    println!("Precompute time: {:?}", precompute_time.elapsed().as_millis());


    // Initialize transcript
    let mut transcript = Transcript::new(b"sumcheck-pallas");

    // bind public inputs
    transcript.append_u64(b"m", m as u64);
    transcript.append_u64(b"l", l as u64);

    for a in alpha.iter() {
        transcript_append_scalar(&mut transcript, b"alpha", a);
    }

    transcript_append_point(&mut transcript, b"sigma", &sigma);

    // prover state
    let mut rho_prefix: Vec<PallasScalar> = Vec::new();
    let mut r_vec: Vec<PallasScalar> = Vec::with_capacity(m);
    let mut gi_triples: Vec<[PallasPoint; 3]> = Vec::with_capacity(m);
    let u_values = [PallasScalar::from(0u64), PallasScalar::from(1u64), PallasScalar::from(2u64)];

    // Sumcheck rounds
    for i in 1..=l {
        let start_gi_time = Instant::now();
        let gi_01 = compute_gi_values(
            m,
            i,
            l,
            alpha,
            &rho_prefix,
            &S_tables,
            &u_values,
        );
        transcript_append_point(&mut transcript, b"g_i(0)", &gi_01[0]);
        transcript_append_point(&mut transcript, b"g_i(1)", &gi_01[1]);
        transcript_append_point(&mut transcript, b"g_i(2)", &gi_01[2]);

        // Verifier challenge r_i
        let r_i = transcript_challenge_scalar(&mut transcript, b"r_i");
        rho_prefix.push(r_i);
        r_vec.push(r_i);
        println!("prover r_{}: {:?}", i, r_i);
        gi_triples.push([gi_01[0], gi_01[1], gi_01[2]]);

        // Check correctness
        if i == 1 {
            assert_eq!(gi_triples[0][0] + gi_triples[0][1], sigma);
        } else {
            assert_eq!(gi_triples[i - 1][0] + gi_triples[i - 1][1], eval_triple_at_alpha(&gi_triples[i-2], r_vec[i-2]));
        }
        println!("prover gi_{}: {:?}", i, start_gi_time.elapsed().as_millis());

    }

    // At this stage prover has sent messages g_1,\ldots,g_l to verifier, while challenges r_1,\ldots,r_l
    // have been sent by the verifier. The reduced sum-check claim is now:
    // sum(h(r_1,\ldots,r_l, b)) = g_l(r_l)
    let restrict_poly_time = Instant::now();
    let mut h_poly = f_poly.restrict(&r_vec);
    let eq_poly = MultilinearPoly::init_with_eq(alpha);
    let mut eq_poly_reduced = eq_poly.restrict(&r_vec);
    println!("Restrict time: {:?}", restrict_poly_time.elapsed().as_millis());

    for i in (l+1)..=m {
        let start_gi_time = Instant::now();
        let stride = h_poly.coeffs.len() / 2;
        let points_even = (0..stride).map(|i| h_poly.coeffs[2*i]).collect::<Vec<_>>();
        let points_odd = (0..stride).map(|i| h_poly.coeffs[2*i+1]).collect::<Vec<_>>();
        let scalars_even = (0..stride).map(|i| eq_poly_reduced.coeffs[2*i]).collect::<Vec<_>>();
        let scalars_odd = (0..stride).map(|i| eq_poly_reduced.coeffs[2*i+1]).collect::<Vec<_>>();
        let sum1 = <PallasPoint as InnerProduct<PallasScalar>>::inner_product_msm(&points_even, &scalars_even);
        let sum2 = <PallasPoint as InnerProduct<PallasScalar>>::inner_product_msm(&points_even, &scalars_odd);
        let sum3 = <PallasPoint as InnerProduct<PallasScalar>>::inner_product_msm(&points_odd, &scalars_even);
        let sum4 = <PallasPoint as InnerProduct<PallasScalar>>::inner_product_msm(&points_odd, &scalars_odd);
        let g0 = sum1;
        let g1 = sum4;
        let g2 = sum1 - (sum2 + sum2 + sum3 + sum3) + (sum4 + sum4 + sum4 + sum4);

        transcript_append_point(&mut transcript, b"g_i(0)", &g0);
        transcript_append_point(&mut transcript, b"g_i(1)", &g1);
        transcript_append_point(&mut transcript, b"g_i(2)", &g2);

        // Verifier challenge r_i
        let r_i = transcript_challenge_scalar(&mut transcript, b"r_i");
        rho_prefix.push(r_i);
        r_vec.push(r_i);
        println!("prover r_{} = {:?}", i, r_i);


        gi_triples.push([g0, g1, g2]);
        assert_eq!(gi_triples[i - 1][0] + gi_triples[i - 1][1], eval_triple_at_alpha(&gi_triples[i-2], r_vec[i-2]));
        h_poly.fold_first(r_i);
        eq_poly_reduced.fold_first(r_i);
        println!("prover round {} = {:?}", i, start_gi_time.elapsed().as_millis());
    }


    (r_vec, gi_triples)
}


/// Verifier for non-interactive sumcheck
pub fn run_verifier_noninteractive(
    m: usize,
    l: usize,
    f_table: &[PallasPoint],
    alpha: &[PallasScalar],
    sigma: &PallasPoint,
    gi_triples: &Vec<[PallasPoint; 3]>,
) -> Result<(), &'static str> {
    assert_eq!(gi_triples.len(), m);

    let mut transcript = Transcript::new(b"sumcheck-pallas");

    // bind public inputs
    transcript.append_u64(b"m", m as u64);
    transcript.append_u64(b"l", l as u64);

    for a in alpha.iter() {
        transcript_append_scalar(&mut transcript, b"alpha", a);
    }

    transcript_append_point(&mut transcript, b"sigma", sigma);

    let mut prev_claim = *sigma;
    let mut r_vec: Vec<PallasScalar> = Vec::with_capacity(m);

    let one = PallasScalar::from(1u64);

    //sumcheck rounds
    for i in 1..m {
        println!("i={}", i);
        let [g0, g1, gu] = gi_triples[i - 1];

        // Check sum g_i(0) + g_i(1) = prev_claim = g_{i-1}(r_{i-1})
        let mut sum = g0;
        sum.add_assign(&g1);
        println!("sum check: g0 + g1 == prev_claim ? {}", sum == prev_claim);
        if sum != prev_claim {
            return Err("Sumcheck consistency failed");
        }
        // Bind g_i(0), g_i(1)
        transcript_append_point(&mut transcript, b"g_i(0)", &g0);
        transcript_append_point(&mut transcript, b"g_i(1)", &g1);

        // Derive random u (Fiat Shamir)
        let u = transcript_challenge_scalar(&mut transcript, b"u");

        // Degree-1 check
        let mut rhs = g0 * (one - u);
        rhs.add_assign(&(g1 * u));
        println!("degree-1: gu == rhs ? {}", gu == rhs);
        if gu != rhs {
            return Err("Degree-1 check failed");
        }
        // Bind g_i(u)
        transcript_append_point(&mut transcript, b"g_i(u)", &gu);

        // Derive random r_i (Fiat Shamir)
        let r_i = transcript_challenge_scalar(&mut transcript, b"r_i");
        transcript_append_scalar(&mut transcript, b"r_i_bind", &r_i);
        let mut gi_ri = g0 * (one - r_i);
        gi_ri.add_assign(&(g1 * r_i)); //gi_ri = g_i(r_i) = g_i(0)*(1-r_i) + g_i(1)*r_i

        r_vec.push(r_i);
        println!("r_i: {:?}", r_i);
        prev_claim = gi_ri;
    }

    println!("performing final mult check");
    // Final multilinear check
    let final_eval = PallasPoint::identity(); //evaluate_f_at_r(f_table, &r_vec);
    println!("prev_claim={:?}", prev_claim);
    println!("final_eval={:?}", final_eval);

    if prev_claim != final_eval {
        return Err("Final polynomial evaluation mismatch");
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{rngs::StdRng, Rng, SeedableRng};

    fn random_scalar(rng: &mut impl Rng) -> PallasScalar {
        PallasScalar::random(rng)
    }

    #[test]
    fn test_group_sumcheck() {
        let mut rng = StdRng::seed_from_u64(42); //safe for testing
        //let mut rng = StdRng::from_entropy();

        let m = 16; //number of variables in f
        let l = 8; //split for MSM and folklore method

        let table_size = 1usize << m;

        // Random f-table: evaluations of f
        let mut f_table = Vec::with_capacity(table_size);
        for _ in 0..table_size {
            let s = random_scalar(&mut rng);
            f_table.push(PallasPoint::generator() * s);
        }

        let f_poly = MultilinearPoly::new(f_table.clone());

        // Random alpha: f(alpha) = sigma
        let mut alpha = Vec::with_capacity(m);
        for _ in 0..m {
            alpha.push(random_scalar(&mut rng));
        }

        let sigma = f_poly.evaluate_msm(&alpha);

        // Run prover
        let start = std::time::Instant::now();
        let (r_vec, gi_triples) =
            run_prover_noninteractive(m, l, &f_table, &alpha, sigma);
        println!("prover took {:?}", start.elapsed());

        // Run verifier
        //let start = std::time::Instant::now();
        //let res = run_verifier_noninteractive(m, l, &f_table, &alpha, &sigma, &gi_triples);
        //println!("verifier took {:?}", start.elapsed());
        //assert!(res.is_ok(), "Sumcheck verification failed");
    }
}
