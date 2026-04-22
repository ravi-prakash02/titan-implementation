use crate::multilinear::MultilinearPoly;
use crate::titantranscript::{append_to_transcript, get_challenge};
use crate::traits::{ByteSerializable, Linear};
use merlin::Transcript;

/// Proof for the field sum-check protocol.
/// Each round j produces d+1 evaluations g_j(0), g_j(1), ..., g_j(degree).
pub struct SumCheckProof<F> {
    pub round_evals: Vec<Vec<F>>,
}

impl<F> SumCheckProof<F> {
    pub fn get_proof_size(&self) -> usize {
        self.round_evals.iter().map(|x| x.len()).sum()
    }
}

/// Verifier output: the evaluation point and expected value for the oracle check.
#[derive(Debug)]
pub struct SubClaim<F> {
    pub point: Vec<F>,
    pub expected_evaluation: F,
}

/// Lagrange interpolation: given evaluations at {0, 1, ..., d}, evaluate the
/// unique degree-d polynomial at `r`.
fn interpolate_and_eval<F: ff::Field>(evals: &[F], r: F) -> F {
    let d = evals.len(); // number of points = degree + 1
    // Build domain points 0, 1, 2, ... as field elements
    let mut domain = Vec::with_capacity(d);
    let mut acc = F::ZERO;
    for _ in 0..d {
        domain.push(acc);
        acc += F::ONE;
    }

    let mut result = F::ZERO;
    for i in 0..d {
        // Lagrange basis polynomial L_i(r)
        let mut num = F::ONE;
        let mut den = F::ONE;
        for j in 0..d {
            if j != i {
                num *= r - domain[j];
                den *= domain[i] - domain[j];
            }
        }
        result += evals[i] * num * den.invert().unwrap();
    }
    result
}

/// Prove the sum-check relation:
///   Σ_{x ∈ {0,1}^n} G(p1(x), ..., pk(x)) = claimed_sum
///
/// `degree` is the total degree of `combine_fn` in each variable (e.g., 1 for
/// identity, 2 for product of two polynomials).
///
/// Returns `(proof, challenges)` where challenges = [r_0, ..., r_{n-1}].
pub fn prove<F>(
    mut polys: Vec<MultilinearPoly<F>>,
    combine_fn: &dyn Fn(&[F]) -> F,
    degree: usize,
    transcript: &mut Transcript,
) -> (SumCheckProof<F>, Vec<F>)
where
    F: ff::Field + Linear<F> + ByteSerializable,
{
    let num_vars = polys[0].num_vars;
    let k = polys.len();

    // Compute claimed sum
    let n = polys[0].coeffs.len();
    let mut claimed_sum = F::ZERO;
    let mut evals_buf = vec![F::ZERO; k];
    for idx in 0..n {
        for p in 0..k {
            evals_buf[p] = polys[p].coeffs[idx];
        }
        claimed_sum += combine_fn(&evals_buf);
    }

    // Append claimed sum to transcript
    append_to_transcript(transcript, b"claimed_sum", &claimed_sum);

    let mut round_evals = Vec::with_capacity(num_vars);
    let mut challenges = Vec::with_capacity(num_vars);

    // Build domain points 0, 1, ..., degree as field elements
    let mut t_values = Vec::with_capacity(degree + 1);
    let mut acc = F::ZERO;
    for _ in 0..=degree {
        t_values.push(acc);
        acc += F::ONE;
    }

    for _round in 0..num_vars {
        let half = polys[0].coeffs.len() / 2;

        // Compute g_j(t) for t = 0, 1, ..., degree
        let mut g_evals = vec![F::ZERO; degree + 1];
        for t in 0..=degree {
            let t_f = t_values[t];
            let one_minus_t = F::ONE - t_f;
            let mut sum = F::ZERO;
            for b in 0..half {
                // For each polynomial, interpolate at t between coeffs[2b] and coeffs[2b+1]
                for p in 0..k {
                    evals_buf[p] =
                        polys[p].coeffs[2 * b] * one_minus_t + polys[p].coeffs[2 * b + 1] * t_f;
                }
                sum += combine_fn(&evals_buf);
            }
            g_evals[t] = sum;
        }

        // Append round evaluations to transcript
        for val in &g_evals {
            append_to_transcript(transcript, b"g_eval", val);
        }

        // Derive challenge
        let r_j: F = get_challenge::<F>(transcript, b"r_j");
        challenges.push(r_j);

        // Fold all polynomials
        for poly in polys.iter_mut() {
            poly.fold_first(r_j);
        }

        round_evals.push(g_evals);
    }

    (SumCheckProof { round_evals }, challenges)
}

/// Verify a sum-check proof. Returns a `SubClaim` containing the evaluation
/// point and the expected value that the oracle must confirm.
pub fn verify<F>(
    proof: &SumCheckProof<F>,
    claimed_sum: F,
    num_vars: usize,
    degree: usize,
    transcript: &mut Transcript,
) -> Result<SubClaim<F>, &'static str>
where
    F: ff::Field + Linear<F> + ByteSerializable,
{
    if proof.round_evals.len() != num_vars {
        return Err("wrong number of rounds in proof");
    }

    // Append claimed sum to transcript (must match prover)
    append_to_transcript(transcript, b"claimed_sum", &claimed_sum);

    let mut prev_claim = claimed_sum;
    let mut challenges = Vec::with_capacity(num_vars);

    for round in 0..num_vars {
        let evals = &proof.round_evals[round];
        if evals.len() != degree + 1 {
            return Err("wrong number of evaluations in round");
        }

        // Check: g_j(0) + g_j(1) == prev_claim
        if evals[0] + evals[1] != prev_claim {
            return Err("round check failed: g(0) + g(1) != prev_claim");
        }

        // Append round evaluations to transcript
        for val in evals {
            append_to_transcript(transcript, b"g_eval", val);
        }

        // Derive challenge (must match prover)
        let r_j: F = get_challenge::<F>(transcript, b"r_j");
        challenges.push(r_j);

        // Update claim via interpolation
        prev_claim = interpolate_and_eval(evals, r_j);
    }

    Ok(SubClaim {
        point: challenges,
        expected_evaluation: prev_claim,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rand_scalar;
    use ff::Field;
    use merlin::Transcript;
    use pasta_curves::pallas::Scalar as F;
    use rand::prelude::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_identity_sumcheck() {
        // G(a) = a, degree 1, single polynomial
        let mut rng = StdRng::seed_from_u64(42);
        let num_vars = 4;
        let n = 1 << num_vars;
        let coeffs: Vec<F> = (0..n).map(|_| rand_scalar(&mut rng)).collect();
        let poly = MultilinearPoly::new(coeffs.clone());

        let combine_fn = |vals: &[F]| -> F { vals[0] };

        // Compute expected sum
        let expected_sum: F = coeffs.iter().copied().sum();

        let mut prover_transcript = Transcript::new(b"field_sumcheck_test");
        let (proof, challenges) =
            prove(vec![poly.clone()], &combine_fn, 1, &mut prover_transcript);

        // Verify
        let mut verifier_transcript = Transcript::new(b"field_sumcheck_test");
        let subclaim =
            verify(&proof, expected_sum, num_vars, 1, &mut verifier_transcript).unwrap();

        assert_eq!(subclaim.point, challenges);

        // Oracle check: evaluate poly at the challenge point
        let oracle_eval = poly.evaluate(&subclaim.point);
        assert_eq!(oracle_eval, subclaim.expected_evaluation);
    }

    #[test]
    fn test_product_sumcheck() {
        // G(a, b) = a * b, degree 2, two polynomials
        let mut rng = StdRng::seed_from_u64(42);
        let num_vars = 5;
        let n = 1 << num_vars;
        let coeffs_p: Vec<F> = (0..n).map(|_| rand_scalar(&mut rng)).collect();
        let coeffs_q: Vec<F> = (0..n).map(|_| rand_scalar(&mut rng)).collect();
        let poly_p = MultilinearPoly::new(coeffs_p.clone());
        let poly_q = MultilinearPoly::new(coeffs_q.clone());

        let combine_fn = |vals: &[F]| -> F { vals[0] * vals[1] };

        // Compute expected sum
        let expected_sum: F = coeffs_p
            .iter()
            .zip(coeffs_q.iter())
            .map(|(a, b)| *a * *b)
            .sum();

        let mut prover_transcript = Transcript::new(b"field_sumcheck_product");
        let (proof, challenges) = prove(
            vec![poly_p.clone(), poly_q.clone()],
            &combine_fn,
            2,
            &mut prover_transcript,
        );

        // Verify
        let mut verifier_transcript = Transcript::new(b"field_sumcheck_product");
        let subclaim =
            verify(&proof, expected_sum, num_vars, 2, &mut verifier_transcript).unwrap();

        assert_eq!(subclaim.point, challenges);

        // Oracle check
        let eval_p = poly_p.evaluate(&subclaim.point);
        let eval_q = poly_q.evaluate(&subclaim.point);
        assert_eq!(eval_p * eval_q, subclaim.expected_evaluation);
    }

    #[test]
    fn test_soundness_rejection() {
        // Tamper with proof and check that verifier rejects
        let mut rng = StdRng::seed_from_u64(42);
        let num_vars = 3;
        let n = 1 << num_vars;
        let coeffs: Vec<F> = (0..n).map(|_| rand_scalar(&mut rng)).collect();
        let poly = MultilinearPoly::new(coeffs.clone());

        let combine_fn = |vals: &[F]| -> F { vals[0] };
        let expected_sum: F = coeffs.iter().copied().sum();

        let mut prover_transcript = Transcript::new(b"field_sumcheck_soundness");
        let (mut proof, _challenges) =
            prove(vec![poly.clone()], &combine_fn, 1, &mut prover_transcript);

        // Tamper: add ONE to the first evaluation of round 0
        proof.round_evals[0][0] += F::ONE;

        let mut verifier_transcript = Transcript::new(b"field_sumcheck_soundness");
        let result = verify(&proof, expected_sum, num_vars, 1, &mut verifier_transcript);

        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "round check failed: g(0) + g(1) != prev_claim"
        );
    }
}
