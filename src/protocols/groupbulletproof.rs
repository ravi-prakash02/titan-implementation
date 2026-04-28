use std::fmt::Debug;
use std::marker::PhantomData;
use std::time::Instant;
use ark_std::iterable::Iterable;
use ark_std::UniformRand;
use merlin::Transcript;
use pasta_msm::pallas;
use rand::rngs::StdRng;
use rand::{thread_rng, SeedableRng};
use crate::multilinear::MultilinearPoly;
use crate::{arkpallastypes, pastatypes};
use crate::arkpallastypes::ArkPallasPoint;
use crate::rand_point;
use crate::titantranscript::{append_to_transcript, get_challenge, get_challenge_u64};
use crate::traits::{ByteSerializable, InnerProduct, Linear};
use crate::utils::{create_smooth_domain, generate_power_vec};

pub trait BulletProofGroup {
    type Instance;
    type Proof;
    type Witness;

    // prove and verify take reference to transcript to allow composition of protocols
    fn prove(
        instance: &Self::Instance,
        witness: &Self::Witness,
        transcript: &mut Transcript,
    ) -> Self::Proof;
    fn verify(
        proof: &Self::Proof,
        instance: &Self::Instance,
        transcript: &mut Transcript,
    ) -> bool;
}

// We reinterpret bulletproof inner product argument as sum-check.
// View initial vectors as polynomials in evaluation basis.
// Each challenge reduces a claim to claim over polynomials folded by the challenge.
// Classical bulletproof corresponds to folding by the "most" significant variable.
// WHIR folded oracles correspond to fixing least significant variable.
// So, the WHIR oracle to generators given by g(x_1,\ldots,x_m) is actually constructed for t = g(x_m,\ldots,x_1).
// restricted polynomial sent by the prover after l challenges is g(
// To enable efficient querying, we commit to g_poly(x_1,...,x_m) as t_poly(x_1,...,x_m)=g_poly(x_m,...,x_1).
// Then g_poly(x,rev(ch_vec)) is simply t_poly(ch_vec,rev(x)). The prover sends g_fold(x) = t(ch_vec, rev(x))
// Thus for y used to check consistency of g_fold with t oracle, we check g_fold(rev(power_vec(y))) = fold(t,ch_vec)(power_vec(y))
#[derive(Clone, Debug)]
pub struct BulletProofParams<G, F> {
    pub m: usize,                   // Number of generators = 1 << m
    pub l: usize,                   // Size of each coset in coset-wise commitment of generators.
    pub g_poly: MultilinearPoly<G>, // polynomial interpolating the generators
    pub c_poly: MultilinearPoly<G>, // canonical polynomial in non-permuted order
    pub g_by_cosets: Vec<Vec<G>>,   // vector of cosets of g
    pub domain_g: Vec<F>,           // domain consisting of the cosets.
}


impl BulletProofParams<pastatypes::Point, pastatypes::Scalar>
{
    pub fn new(m: usize, l: usize) -> Self {
        let mut rng = StdRng::from_entropy();
        let g_vec: Vec<pastatypes::Point> = (0..(1 << m)).into_iter().map(|_| rand_point(&mut rng)).collect();
        // compute canonical polynomial for generators
        let c_poly = MultilinearPoly::new(g_vec.clone());
        // Next we compute permuted polynomial which will be used to commit generators
        let n = 1usize << m;
        let mut perm: Vec<usize> = vec![0; n];
        for j in 0..n {
            perm[j] = (0..m).into_iter().fold(0usize, |acc, i| acc + (((j >> i) % 2) << (m - i - 1)));
            assert!(perm[j] < n);
        }
        let t_vec: Vec<pastatypes::Point> = (0..(1 << m)).into_iter().map(|i| g_vec[perm[i]]).collect();
        let g_poly = MultilinearPoly::<pastatypes::Point>::new(t_vec);
        let mut g_by_cosets: Vec<Vec<pastatypes::Point>> = Vec::new();
        let domain = create_smooth_domain(l);
        let domain_g: Vec<pastatypes::Scalar> = domain.1;

        for idx in 0..domain_g.len() {
            // build powers of y vector
            let y = domain_g[idx];
            let pow_y = generate_power_vec(y, m-l, false);
            g_by_cosets.push(g_poly.restrict_msb(&pow_y).coeffs);
        }

        Self {
            m,
            l,
            g_poly: g_poly,
            c_poly: c_poly,
            g_by_cosets: g_by_cosets,
            domain_g: domain_g,
        }
    }
}

impl BulletProofParams<arkpallastypes::ArkPallasPoint, pastatypes::Scalar>
{
    pub fn new(m: usize, l: usize) -> Self {
        let mut rng = StdRng::from_entropy();
        let g_vec: Vec<ArkPallasPoint> = (0..(1 << m)).into_iter().map(|_| ArkPallasPoint(ark_pallas::Projective::rand(&mut rng))).collect();
        // compute canonical polynomial for generators
        let c_poly = MultilinearPoly::new(g_vec.clone());
        // Next we compute permuted polynomial which will be used to commit generators
        let n = 1usize << m;
        let mut perm: Vec<usize> = vec![0; n];
        for j in 0..n {
            perm[j] = (0..m).into_iter().fold(0usize, |acc, i| acc + (((j >> i) % 2) << (m - i - 1)));
            assert!(perm[j] < n);
        }
        let t_vec = (0..(1 << m)).into_iter().map(|i| g_vec[perm[i]]).collect();
        let g_poly = MultilinearPoly::<ArkPallasPoint>::new(t_vec);
        let mut g_by_cosets: Vec<Vec<ArkPallasPoint>> = Vec::new();
        let domain = create_smooth_domain(l);
        let domain_g: Vec<pastatypes::Scalar> = domain.1;

        for idx in 0..domain_g.len() {
            // build powers of y vector
            let y = domain_g[idx];
            let pow_y = generate_power_vec(y, m-l, false);
            g_by_cosets.push(g_poly.restrict_msb(&pow_y).coeffs);
        }

        Self {
            m,
            l,
            g_poly: g_poly,
            c_poly: c_poly,
            g_by_cosets: g_by_cosets,
            domain_g: domain_g,
        }
    }
}


// Instance is the statement Comm(a; g) = commitment and a(alpha) = sigma.
pub struct BulletProofInstance<G, F> {
    pub pp: BulletProofParams<G, F>,
    pub commitment: G,
    pub alpha: Vec<F>,
    pub sigma: F,
}

pub struct BulletProof<G, F> {
    // vectors containing round messages
    pub cm_L: Vec<G>,
    pub cm_R: Vec<G>,
    pub lf_L: Vec<F>,
    pub lf_R: Vec<F>,
    // folded polynomials for final message.
    g_folded_poly: MultilinearPoly<G>,
    a_folded_poly: MultilinearPoly<F>,
    phantom_data: PhantomData<(G, F)>,
}

impl<G, F> BulletProof<G, F> {
    pub fn get_proof_size(&self) -> usize {
        let mut proof_size = 0usize;
        proof_size += (4*self.cm_L.len() + self.g_folded_poly.coeffs.len() + self.a_folded_poly.coeffs.len());
        proof_size
    }
}

pub struct BulletProofWitness<G, F> {
    pub pp: BulletProofParams<G, F>,
    pub a_poly: MultilinearPoly<F>,
}

pub struct BulletProofVerifier<F> {
    phantom_data: PhantomData<F>,
}

pub struct ProtoBulletProofGroup<G, F> {
    phantom_data: PhantomData<(G, F)>,
}


impl<G, F> BulletProofGroup for ProtoBulletProofGroup<G, F>
where
    G: Linear<F> + InnerProduct<F, Output=G> + ByteSerializable + PartialEq + Debug,
    G::Affine: InnerProduct<F, Output=G>,
    F: ff::Field + Linear<F> + InnerProduct<F, Output=F> + ByteSerializable,
{
    type Instance = BulletProofInstance<G, F>;
    type Proof = BulletProof<G, F>;
    type Witness = BulletProofWitness<G, F>;

    fn prove(instance: &Self::Instance, witness: &Self::Witness, transcript: &mut Transcript) -> Self::Proof {
        let mut cm_L: Vec<G> = vec![];
        let mut cm_R: Vec<G> = vec![];
        let mut lf_L: Vec<F> = vec![];
        let mut lf_R: Vec<F> = vec![];
        let mut a_poly = MultilinearPoly::new(witness.a_poly.coeffs.clone());
        let mut eq_poly = MultilinearPoly::init_with_eq(&instance.alpha);

        transcript.append_message(b"protocol", b"bulletproof");
        let mut n: usize = 1usize << instance.pp.m; //number of coefficients
        let mut ch_vec: Vec<F> = Vec::new();
        let mut P: G = instance.commitment.clone();
        let mut V: F = instance.sigma;
        // canonical representation is used in loop below
        let ck_affine = G::to_affine(&instance.pp.c_poly.coeffs.clone());

        // add instance to transcript
        transcript.append_u64(b"m", instance.pp.m as u64);
        transcript.append_u64(b"l", instance.pp.l as u64);
        append_to_transcript(transcript, b"commitment", &instance.commitment);
        for a in instance.alpha.clone() {
            append_to_transcript(transcript, b"alpha", &a);
        }
        append_to_transcript(transcript, b"sigma", &instance.sigma);


        for i in 0..instance.pp.l {
            // construct round commitments and values
            let ch_rev = ch_vec.iter().rev().cloned().collect::<Vec<F>>();
            let eq_ch: MultilinearPoly<F> = MultilinearPoly::init_with_eq(&ch_rev);

            let mut l_factors: Vec<G> = vec![G::zero(); 1usize << i];
            let mut r_factors: Vec<G> = vec![G::zero(); 1usize << i];
            for b in 0..(1usize << i) {
                //splitting the commitment key into left and right chunks
                let start_l = b * (1 << (instance.pp.m - i));
                let end_l = start_l + (1 << (instance.pp.m - i - 1));
                let start_r = end_l;
                let end_r = start_r + (1 << (instance.pp.m - i - 1));
                //print!("i = {}, start_l = {}, end_l = {}, start_r = {}, end_r = {}\n", i, start_l, end_l, start_r, end_r);
                l_factors[b] = G::Affine::inner_product_msm(&ck_affine[start_l..end_l], &a_poly.coeffs[(n / 2)..n]);
                r_factors[b] = G::Affine::inner_product_msm(&ck_affine[start_r..end_r], &a_poly.coeffs[0..(n / 2)]);
            }

            let cL = G::inner_product_msm(&l_factors, &eq_ch.coeffs);
            let cR = G::inner_product_msm(&r_factors, &eq_ch.coeffs);
            let aL = F::inner_product(&a_poly.coeffs[..(n / 2)], &eq_poly.coeffs[(n / 2)..]);
            let aR = F::inner_product(&a_poly.coeffs[(n / 2)..], &eq_poly.coeffs[..(n / 2)]);
            // aL=<aL,eq_R>, aR=<aR,eq_L>
            cm_L.push(cL);
            cm_R.push(cR);
            lf_L.push(aL);
            lf_R.push(aR);

            append_to_transcript(transcript, b"cmL", &cL);
            append_to_transcript(transcript, b"cmR", &cR);
            append_to_transcript(transcript, b"lfL", &aL);
            append_to_transcript(transcript, b"lfR", &aR);

            //get new challenge
            let x = get_challenge::<F>(transcript, b"x");
            ch_vec.push(x);
            // update the witness and linear form polynomials by folding them
            a_poly.fold(F::ONE - x);
            eq_poly.fold(x);

            // compute the new commitment P and value V
            // P = x(1-x)P + (1-x)^2 cL + x^2 cR
            // V = x(1-x)V + (1-x)^2 aL + x^2 aR
            V = V * x * (F::ONE - x) + aL * (F::ONE - x) * (F::ONE - x) + aR * x * x;
            P = P * x * (F::ONE - x) + cL * (F::ONE - x) * (F::ONE - x) + cR * x * x;
            n = n / 2;
        }

        // At this stage, the prover sends t_poly(c_1,\ldots,c_l,...), a_poly
        let ch_rev = ch_vec.iter().rev().cloned().collect::<Vec<F>>();

        // send the folding of canonical polynomial. During consistency check with the original
        // oracle (which was shuffled before committing), we will permute the query vector (reverse it) to
        // simulate the correct check.
        let g_poly_reduced = instance.pp.c_poly.restrict_msb(&ch_rev);
        // append polynomials to transcript
        assert_eq!(g_poly_reduced.coeffs.len(), a_poly.coeffs.len(), "Polynomials don't match in size");
        for i in 0..g_poly_reduced.coeffs.len() {
            append_to_transcript(transcript, b"g_poly", &g_poly_reduced.coeffs[i]);
        }
        for i in 0..a_poly.coeffs.len() {
            append_to_transcript(transcript, b"a_poly", &a_poly.coeffs[i]);
        }

        BulletProof {
            cm_L,
            cm_R,
            lf_L,
            lf_R,
            g_folded_poly: g_poly_reduced,
            a_folded_poly: a_poly,
            phantom_data: Default::default(),
        }
    }

    fn verify(proof: &Self::Proof, instance: &Self::Instance, transcript: &mut Transcript) -> bool {
        transcript.append_message(b"protocol", b"bulletproof");
        let mut n: usize = 1usize << instance.pp.m; //number of coefficients
        let mut ch_vec: Vec<F> = Vec::new();
        let mut P: G = instance.commitment.clone();
        let mut V: F = instance.sigma;

        // add instance to transcript
        transcript.append_u64(b"m", instance.pp.m as u64);
        transcript.append_u64(b"l", instance.pp.l as u64);
        append_to_transcript(transcript, b"commitment", &instance.commitment);
        for a in instance.alpha.clone() {
            append_to_transcript(transcript, b"alpha", &a);
        }
        append_to_transcript(transcript, b"sigma", &instance.sigma);

        // first verify the round messages
        for i in 0usize..instance.pp.l {
            let cL = proof.cm_L[i];
            let cR = proof.cm_R[i];
            let aL = proof.lf_L[i];
            let aR = proof.lf_R[i];
            append_to_transcript(transcript, b"cmL", &cL);
            append_to_transcript(transcript, b"cmR", &cR);
            append_to_transcript(transcript, b"lfL", &aL);
            append_to_transcript(transcript, b"lfR", &aR);

            //get new challenge
            let x = get_challenge::<F>(transcript, b"x");
            println!("Round: {} ch: {:?}", i, x);
            ch_vec.push(x);
            // compute the new commitment P and value V
            // P = x(1-x)P + (1-x)^2 cL + x^2 cR
            // V = x(1-x)V + (1-x)^2 aL + x^2 aR
            V = V * x * (F::ONE - x) + aL * (F::ONE - x) * (F::ONE - x) + aR * x * x;
            P = P * x * (F::ONE - x) + cL * (F::ONE - x) * (F::ONE - x) + cR * x * x;
        }

        // Compute query locations for consistency check
        for i in 0..proof.g_folded_poly.coeffs.len() {
            append_to_transcript(transcript, b"g_poly", &proof.g_folded_poly.coeffs[i]);
        }
        for i in 0..proof.a_folded_poly.coeffs.len() {
            append_to_transcript(transcript, b"a_poly", &proof.a_folded_poly.coeffs[i]);
        }

        // get opening challenges
        let num_queries = 35usize;
        let y_vec: Vec<usize> = (0usize..num_queries)
            .into_iter()
            .map(|_| {
                (get_challenge_u64(transcript, b"challenge_y") as usize) % instance.pp.domain_g.len()
            })
            .collect();

        // Batch verification of consistency check relations:
        // We need to verify the following:
        // For each y in query set: inner_prod(g_fold, eq(\cdot,rev(y))) = inner_prod(leaf(y), eq(ch_vec,\cdot))
        // We aggregate the above inner products using random challenge r. One the left, we aggregate vectors rev(y)
        // On the right, we concatenate leaf(y) for each y, and r powers scaled eq(ch_vec,\cdot).

        let eq_vec = MultilinearPoly::init_with_eq(&ch_vec).coeffs;
        let r = F::random(&mut thread_rng()); // @todo make r part of proof transcript.
        // init_leaf builds the concatenation l_1||l_2||..||l_t,
        // init_eq_vec builds the concatenation eq_vec|| r.eq_vec || r^2.eq_vec || r^{t-1} eq_vec
        let mut init_leaf = instance.pp.g_by_cosets[y_vec[0]].clone();
        let mut init_eq_vec = eq_vec.clone();

        // init_ypow_vec aggregates the query dependent vectors for the left hand side
        let y_val = instance.pp.domain_g[y_vec[0]];
        let ypow_vec: Vec<F> = generate_power_vec(y_val, instance.pp.m - instance.pp.l, true);
        let mut init_ypow_vec = MultilinearPoly::init_with_eq(&ypow_vec).coeffs;

        let mut r_powers: Vec<F> = Vec::new();
        let mut r_power = F::ONE;
        r_powers.push(r_power);
        let start = Instant::now();
        for i in 1..num_queries {
            r_power = r_power * r;
            init_leaf.extend(instance.pp.g_by_cosets[y_vec[i]].iter().cloned());
            let scaled_eq_vec: Vec<F> = eq_vec.iter().map(|x| *x * r_power).collect();
            init_eq_vec.extend(scaled_eq_vec.iter().cloned());

            let y_val = instance.pp.domain_g[y_vec[i]];
            let ypow_vec: Vec<F> = generate_power_vec(y_val, instance.pp.m - instance.pp.l, true);
            let scaled_ypow_vec: Vec<F> = MultilinearPoly::init_with_eq(&ypow_vec).coeffs.iter().map(|x| *x * r_power).collect();
            // update aggregated power vector
            for (x, y) in init_ypow_vec.iter_mut().zip(scaled_ypow_vec.iter()) {
                *x += y;
            }
            r_powers.push(r_power);
        }
        //
        //assert_eq!(G::inner_product_msm(&proof.g_folded_poly.coeffs, &init_ypow_vec),
        //    G::inner_product_msm(&init_leaf, &init_eq_vec));
        // check that reduced statement holds
        // assert_eq!(G::inner_product_msm(&proof.g_folded_poly.coeffs, &proof.a_folded_poly.coeffs), P);
        // Batch the above MSM checks
        let mut generators = proof.g_folded_poly.coeffs.clone();
        generators.extend(init_leaf.as_slice());
        let mut batched_scalars =
            init_ypow_vec.iter().zip(proof.a_folded_poly.coeffs.iter()).map(|x| *x.0 + *x.1 * r).collect::<Vec<_>>();
        let _ = &mut batched_scalars.extend(init_eq_vec.iter().map(|x| F::ZERO - *x));

        assert_eq!(G::inner_product_msm(&generators, &batched_scalars), P * r);


        true
    }
}

mod tests {
    use std::time::Instant;
    use crate::rand_scalar;
    use super::*;

    #[test]
    fn test_bulletproof_prover()
    {
        let m = 9usize;
        let l = 3usize;
        let mut rng = StdRng::from_entropy();

        let bp_params = BulletProofParams::<pastatypes::Point, pastatypes::Scalar>::new(m, l);
        let a_coeffs: Vec<pastatypes::Scalar> = (0..(1usize << m)).into_iter().map(|_| rand_scalar(&mut rng)).collect();
        let a_poly = MultilinearPoly::new(a_coeffs);
        // commitment according to canonical polynomial
        let commitment = pastatypes::Point::inner_product_msm(&bp_params.c_poly.coeffs, &a_poly.coeffs);
        let alpha: Vec<pastatypes::Scalar> = (0..m).into_iter().map(|_| rand_scalar(&mut rng)).collect();
        let sigma = a_poly.evaluate(&alpha);

        let bp_instance = BulletProofInstance::<pastatypes::Point, pastatypes::Scalar> {
            pp: bp_params.clone(),
            commitment,
            alpha: alpha,
            sigma: sigma,
        };

        let bp_witness = BulletProofWitness::<pastatypes::Point, pastatypes::Scalar> {
            pp: bp_params,
            a_poly,
        };

        let mut transcript = Transcript::new(b"Bulletproof");
        let start = Instant::now();
        let proof = ProtoBulletProofGroup::<pastatypes::Point, pastatypes::Scalar>::prove(
            &bp_instance,
            &bp_witness,
            &mut transcript,
        );
        println!("Time to do bulletproofs {} msec", start.elapsed().as_millis());

        let mut transcript = Transcript::new(b"Bulletproof");
        let start = Instant::now();
        let res = ProtoBulletProofGroup::<pastatypes::Point, pastatypes::Scalar>::verify(
            &proof,
            &bp_instance,
            &mut transcript,
        );
        println!("Time to verify bulletproofs {} msec", start.elapsed().as_millis());
    }
}



