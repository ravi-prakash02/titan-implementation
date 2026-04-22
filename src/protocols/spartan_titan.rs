use std::time::Instant;
use ark_crypto_primitives::sponge::Absorb;
use ff::{Field, PrimeField};
use merlin::Transcript;
use crate::group_whir_committer::GroupWhirCommitment;
use crate::merkle_tree::{Sha256Digest, Sha256MerkleTreeParams};
use crate::multilinear::MultilinearPoly;
use crate::protocols::field_sumcheck;
use crate::protocols::titanpcs::{AggregateEvalProof, PolyCommitment, TitanBatchEvalProof, TitanPolyCommitment, TitanProverParams, TitanSetupConfig, TitanVerifierParams};
use crate::r1cs::{SparseMatrix, R1CS};
use crate::titantranscript::{append_to_transcript, get_challenge};
use crate::traits::{ByteSerializable, InnerProduct, Linear};

type G = crate::pastatypes::Point;
type F = crate::pastatypes::Scalar;

pub struct TitanSnarkConfig {
    pub single_oracle_config: TitanSetupConfig,
    pub trusted_oracle_config: TitanSetupConfig,
    pub aggregate_oracle_config: TitanSetupConfig,
}

pub struct TitanSnarkInstance {
    pub a_matrix: SparseMatrix<F>,
    pub b_matrix: SparseMatrix<F>,
    pub c_matrix: SparseMatrix<F>,
}

pub struct TitanIndexPolys
{
    pub row_poly: MultilinearPoly<F>,
    pub col_poly: MultilinearPoly<F>,
    pub val_A_poly: MultilinearPoly<F>,
    pub val_B_poly: MultilinearPoly<F>,
    pub val_C_poly: MultilinearPoly<F>,

    // Commit with opens (we assume verifier knows these)
    pub comm_row_poly: GroupWhirCommitment<Sha256MerkleTreeParams<G>, G>,
    pub comm_col_poly: GroupWhirCommitment<Sha256MerkleTreeParams<G>, G>,
    pub comm_val_A_poly: GroupWhirCommitment<Sha256MerkleTreeParams<G>, G>,
    pub comm_val_B_poly: GroupWhirCommitment<Sha256MerkleTreeParams<G>, G>,
    pub comm_val_C_poly: GroupWhirCommitment<Sha256MerkleTreeParams<G>, G>,
}

/// Spartan proof containing three sum-check sub-proofs.
///
/// First sum-check:  Σ_x eq(x,τ)·(Az(x)·Bz(x) − Cz(x)) = 0
/// Second sum-check: Σ_y (eA(y) + β·eB(y) + β²·eC(y))·z(y) = a_at_rx + β·b_at_rx + β²·c_at_rx
/// Third sum-check:  Batched memory-checking + matrix decomposition (Eq 11, 12, 13 of main.pdf)
pub struct SpartanProof {
    // --- First sum-check ---
    pub sumcheck_proof_1: field_sumcheck::SumCheckProof<F>,
    pub a_at_rx: F,
    pub b_at_rx: F,
    pub c_at_rx: F,

    // --- Second sum-check ---
    pub sumcheck_proof_2: field_sumcheck::SumCheckProof<F>,
    pub eA_at_ry: F,
    pub eB_at_ry: F,
    pub eC_at_ry: F,
    pub z_at_ry: F,

    // --- Commitment roots ---
    pub comm_pqm_root: Sha256Digest,
    pub comm_g_root: Sha256Digest,
    pub comm_z_root: Sha256Digest,

    // --- Third batched sum-check (degree 6) ---
    pub sumcheck_proof_3: field_sumcheck::SumCheckProof<F>,
    /// Oracle evaluations at the third sum-check challenge point r_3
    pub val_A_at_r3: F,
    pub val_B_at_r3: F,
    pub val_C_at_r3: F,
    pub p_at_r3: F,
    pub q_at_r3: F,
    pub g_at_r3: F,
    pub row_at_r3: F,
    pub col_at_r3: F,
    pub m_r_at_r3: F,
    pub m_c_at_r3: F,

    // --- Fourth sum-check (reducing z(r_y) and g(r_3) to common point r_4) ---
    pub sumcheck_proof_4: field_sumcheck::SumCheckProof<F>,
    pub z_at_r4: F,
    pub g_at_r4: F,

    // --- PCS evaluation proofs ---
    pub agg_pqm_proof: AggregateEvalProof<G, F>,
    pub batch_index_proof: TitanBatchEvalProof<G, F>,
    pub batch_zg_proof: TitanBatchEvalProof<G, F>,
}

impl SpartanProof {
    pub fn get_proof_size(&self) -> usize {
        let mut proof_size = 0usize;
        proof_size += self.sumcheck_proof_1.get_proof_size();
        proof_size += 3;
        proof_size += self.sumcheck_proof_2.get_proof_size();
        proof_size += 7;
        proof_size += self.sumcheck_proof_3.get_proof_size();
        proof_size += 10;
        proof_size += self.sumcheck_proof_4.get_proof_size();
        proof_size += 2;
        proof_size += self.agg_pqm_proof.eval_proof.get_proof_size();
        proof_size += self.batch_zg_proof.get_prooof_size();
        proof_size += self.batch_index_proof.get_honest_oracle_proof_size();
        proof_size
    }
}

pub struct TitanSnark {}

impl TitanSnark {
    // setup calls the setup of underlying PCS
    pub fn setup(config: &TitanSnarkConfig) -> Vec<(TitanProverParams<G, F>, TitanVerifierParams<G, F>)> {
        vec![
            TitanPolyCommitment::<G, F>::setup(&config.single_oracle_config),
            TitanPolyCommitment::<G,F>::setup(&config.trusted_oracle_config),
            TitanPolyCommitment::<G,F>::setup(&config.aggregate_oracle_config),
        ]

    }

    pub fn index(
        instance: &R1CS<F>,
        prover_params: &TitanProverParams<G, F>,
        max_nonzeros: usize,
    ) -> TitanIndexPolys {
        let nA = instance.a.entries.len();
        let nB = instance.b.entries.len();
        let nC = instance.c.entries.len();
        let N = nA + nB + nC;
        assert!(N <= max_nonzeros, "Total entries {} exceeds max_nonzeros {}", N, max_nonzeros);

        // Use max_nonzeros as the padded size (must be power of 2)
        let N_padded = max_nonzeros.next_power_of_two();

        let mut row = vec![F::ZERO; N_padded];
        let mut col = vec![F::ZERO; N_padded];
        let mut val_A = vec![F::ZERO; N_padded];
        let mut val_B = vec![F::ZERO; N_padded];
        let mut val_C = vec![F::ZERO; N_padded];

        // Fill from A entries: indices [0, nA)
        for (i, &(r, c, v)) in instance.a.entries.iter().enumerate() {
            row[i] = F::from(r as u64);
            col[i] = F::from(c as u64);
            val_A[i] = v;
        }

        // Fill from B entries: indices [nA, nA + nB)
        // val_B must be at the SAME positions as row/col for Eq 10 decomposition
        for (i, &(r, c, v)) in instance.b.entries.iter().enumerate() {
            row[nA + i] = F::from(r as u64);
            col[nA + i] = F::from(c as u64);
            val_B[nA + i] = v;
        }

        // Fill from C entries: indices [nA + nB, N)
        for (i, &(r, c, v)) in instance.c.entries.iter().enumerate() {
            row[nA + nB + i] = F::from(r as u64);
            col[nA + nB + i] = F::from(c as u64);
            val_C[nA + nB + i] = v;
        }

        // Construct multilinear polynomials as MLE extensions
        let row_poly = MultilinearPoly::new(row);
        let col_poly = MultilinearPoly::new(col);
        let val_A_poly = MultilinearPoly::new(val_A);
        let val_B_poly = MultilinearPoly::new(val_B);
        let val_C_poly = MultilinearPoly::new(val_C);

        // Commit each polynomial using TitanPCS
        let comm_row_poly = TitanPolyCommitment::<G, F>::commit(prover_params, &row_poly);
        let comm_col_poly = TitanPolyCommitment::<G, F>::commit(prover_params, &col_poly);
        let comm_val_A_poly = TitanPolyCommitment::<G, F>::commit(prover_params, &val_A_poly);
        let comm_val_B_poly = TitanPolyCommitment::<G, F>::commit(prover_params, &val_B_poly);
        let comm_val_C_poly = TitanPolyCommitment::<G, F>::commit(prover_params, &val_C_poly);

        TitanIndexPolys {
            row_poly,
            col_poly,
            val_A_poly,
            val_B_poly,
            val_C_poly,
            comm_row_poly,
            comm_col_poly,
            comm_val_A_poly,
            comm_val_B_poly,
            comm_val_C_poly,
        }
    }

    /// We assume that witness is 0-extended to match max non-zero entries in instance.
    pub fn prove(
        instance: &R1CS<F>,
        index: &TitanIndexPolys,
        prover_params_single: &TitanProverParams<G, F>,
        prover_params_trusted: &TitanProverParams<G,F>,
        prover_params_agg: &TitanProverParams<G,F>,
        witness: &[F],
    ) -> SpartanProof {
        let start = Instant::now();
        // Compute Az, Bz, Cz vectors
        let az = instance.a.mult_vec(witness);
        let bz = instance.b.mult_vec(witness);
        let cz = instance.c.mult_vec(witness);

        for i in 0..az.len() {
            assert_eq!(az[i] * bz[i], cz[i], "R1CS not satisfied at constraint {}", i);
        }

        // Pad constraints to next power of 2
        let num_constraints = instance.m;
        let m_padded = num_constraints.next_power_of_two();
        let num_vars_x = (m_padded as f64).log2() as usize;

        // Pad witness to next power of 2
        let n_padded = witness.len().next_power_of_two();
        let num_vars_y = (n_padded as f64).log2() as usize;

        let mut az_padded = az;
        let mut bz_padded = bz;
        let mut cz_padded = cz;
        az_padded.resize(m_padded, F::ZERO);
        bz_padded.resize(m_padded, F::ZERO);
        cz_padded.resize(m_padded, F::ZERO);

        let a_poly = MultilinearPoly::new(az_padded);
        let b_poly = MultilinearPoly::new(bz_padded);
        let c_poly = MultilinearPoly::new(cz_padded);

        // =====================================================================
        // First sum-check: Σ_x eq(x, τ) · (Az(x) · Bz(x) − Cz(x)) = 0
        // =====================================================================
        let mut transcript = Transcript::new(b"spartan");
        transcript.append_u64(b"num_vars_x", num_vars_x as u64);
        transcript.append_u64(b"num_vars_y", num_vars_y as u64);

        let tau: Vec<F> = (0..num_vars_x)
            .map(|_| get_challenge::<F>(&mut transcript, b"tau"))
            .collect();

        let eq_tau = MultilinearPoly::<F>::init_with_eq(&tau);

        let combine_fn_1 = |vals: &[F]| -> F {
            vals[0] * (vals[1] * vals[2] - vals[3])
        };

        let (sumcheck_proof_1, challenges_1) = field_sumcheck::prove(
            vec![eq_tau, a_poly.clone(), b_poly.clone(), c_poly.clone()],
            &combine_fn_1,
            3,
            &mut transcript,
        );
        println!("Time for Spartan Sum-check 1: {}", start.elapsed().as_millis());

        let start = Instant::now();
        let r_x = challenges_1;
        let a_at_rx = a_poly.evaluate(&r_x);
        let b_at_rx = b_poly.evaluate(&r_x);
        let c_at_rx = c_poly.evaluate(&r_x);

        // =====================================================================
        // Second sum-check:
        //   Σ_y (eA(y) + β·eB(y) + β²·eC(y)) · z(y) = a_at_rx + β·b_at_rx + β²·c_at_rx
        //
        // where eA[y] = Σ_x eq(r_x, x) · A(x, y)  (sparse accumulation)
        // =====================================================================

        // Append first sum-check oracle evals to transcript, derive beta
        append_to_transcript(&mut transcript, b"a_at_rx", &a_at_rx);
        append_to_transcript(&mut transcript, b"b_at_rx", &b_at_rx);
        append_to_transcript(&mut transcript, b"c_at_rx", &c_at_rx);
        let beta: F = get_challenge::<F>(&mut transcript, b"beta");
        let beta2 = beta * beta;

        // Precompute eq(r_x, ·) on the boolean hypercube {0,1}^{num_vars_x}
        let eq_rx = MultilinearPoly::<F>::init_with_eq(&r_x);

        // eM = M^T · eq(r_x, ·)  via sparse transpose multiplication
        let eA = instance.a.mult_vec_transpose(&eq_rx.coeffs);
        let eB = instance.b.mult_vec_transpose(&eq_rx.coeffs);
        let eC = instance.c.mult_vec_transpose(&eq_rx.coeffs);

        // Combine: combined(y) = eA(y) + β·eB(y) + β²·eC(y)
        let combined: Vec<F> = (0..n_padded)
            .map(|i| eA[i] + beta * eB[i] + beta2 * eC[i])
            .collect();
        let combined_poly = MultilinearPoly::new(combined);

        // z polynomial (witness padded to power of 2)
        let mut z_padded = witness.to_vec();
        z_padded.resize(n_padded, F::ZERO);
        let z_poly = MultilinearPoly::new(z_padded);

        // Claimed sum: a_at_rx + β·b_at_rx + β²·c_at_rx
        let claimed_sum_2 = a_at_rx + beta * b_at_rx + beta2 * c_at_rx;

        // Degree-2 sum-check: Σ_y combined(y) · z(y)
        let combine_fn_2 = |vals: &[F]| -> F { vals[0] * vals[1] };

        let (sumcheck_proof_2, challenges_2) = field_sumcheck::prove(
            vec![combined_poly.clone(), z_poly.clone()],
            &combine_fn_2,
            2,
            &mut transcript,
        );
        println!("Time for Spartan Sumcheck 2: {}", start.elapsed().as_millis());


        let start = Instant::now();
        let r_y = &challenges_2;
        let eA_poly = MultilinearPoly::new(eA);
        let eB_poly = MultilinearPoly::new(eB);
        let eC_poly = MultilinearPoly::new(eC);
        let eA_at_ry = eA_poly.evaluate(r_y);
        let eB_at_ry = eB_poly.evaluate(r_y);
        let eC_at_ry = eC_poly.evaluate(r_y);
        let z_at_ry = z_poly.evaluate(r_y);

        // =====================================================================
        // Memory checking: compute p, q, m_r, m_c
        // =====================================================================

        // Index polynomial dimension
        let N_padded = 1usize << index.row_poly.num_vars;

        // Tx = eq(r_x, ·) already computed as eq_rx (size m_padded)
        // Ty = eq(r_y, ·)
        let eq_ry = MultilinearPoly::<F>::init_with_eq(r_y);

        // Convert index polynomial field elements to usize for table lookup
        let to_usize = |f: &F| -> usize {
            let repr = f.to_repr();
            u64::from_le_bytes(repr.as_ref()[..8].try_into().unwrap()) as usize
        };

        // p[i] = Tx[row[i]] = eq(r_x, row[i]), q[i] = Ty[col[i]] = eq(r_y, col[i])
        let p: Vec<F> = (0..N_padded).map(|i| eq_rx.coeffs[to_usize(&index.row_poly.coeffs[i])]).collect();
        let q: Vec<F> = (0..N_padded).map(|i| eq_ry.coeffs[to_usize(&index.col_poly.coeffs[i])]).collect();
        let p_poly = MultilinearPoly::new(p);
        let q_poly = MultilinearPoly::new(q);

        // m_r[x] = number of indices i where row(i) == x (multiplicity of row value x)
        let mut m_r = vec![F::ZERO; N_padded];
        for i in 0..N_padded {
            m_r[to_usize(&index.row_poly.coeffs[i])] += F::ONE;
        }
        let m_r_poly = MultilinearPoly::new(m_r);

        // m_c[x] = number of indices i where col(i) == x (multiplicity of col value x)
        let mut m_c = vec![F::ZERO; N_padded];
        for i in 0..N_padded {
            m_c[to_usize(&index.col_poly.coeffs[i])] += F::ONE;
        }
        let m_c_poly = MultilinearPoly::new(m_c);
        println!("Time to compute auxiliary polynomials in round 2: {}", start.elapsed().as_millis());

        let start = Instant::now();
        // Commit to p, q, m_r, m_c and add commitment roots to transcript
        let comm_pqm = TitanPolyCommitment::<G,F>::aggregate_commit(
            prover_params_agg,
            &vec![p_poly.clone(), q_poly.clone(), m_r_poly.clone(), m_c_poly.clone()],
        );
        println!("Time to commit aggregate polynomial: {}", start.elapsed().as_millis());

        let start = Instant::now();
        transcript.append_message(b"comm_p_q_m", &comm_pqm.root.to_sponge_bytes_as_vec());

        // Obtain alpha, beta, epsilon as challenges
        let alpha_3: F = get_challenge::<F>(&mut transcript, b"alpha");
        let beta_3: F = get_challenge::<F>(&mut transcript, b"beta_3");
        let epsilon: F = get_challenge::<F>(&mut transcript, b"epsilon");

        // =====================================================================
        // Third batched sum-check: memory checking + matrix decomposition
        // (Eq 11, 12, 13 from main.pdf, page 26-27)
        //
        // Batches five sub-instances into one degree-6 sum-check:
        //   SC1: Σ_x val_A(x)·p(x)·q(x) = eA_at_ry
        //   SC2: Σ_x val_B(x)·p(x)·q(x) = eB_at_ry
        //   SC3: Σ_x val_C(x)·p(x)·q(x) = eC_at_ry
        //   SC4a: Σ_x eq(x,η)·(t(x) − g(x)·d1·d2·d3·d4) = 0
        //   SC4b: Σ_x g(x) = 0
        // =====================================================================

        let num_vars_idx = index.row_poly.num_vars;

        // Compute denominators d1..d4 pointwise (setting X=α, Y=β from the PDF)
        // d1(x) = α + β·id(x) + eq(r_x, x)
        // d2(x) = α + β·row(x) + p(x)
        // d3(x) = α + β·id(x) + eq(r_y, x)
        // d4(x) = α + β·col(x) + q(x)
        let d1: Vec<F> = (0..N_padded).map(|i| {
            alpha_3 + beta_3 * F::from(i as u64) + eq_rx.coeffs[i]
        }).collect();
        let d2: Vec<F> = (0..N_padded).map(|i| {
            alpha_3 + beta_3 * index.row_poly.coeffs[i] + p_poly.coeffs[i]
        }).collect();
        let d3: Vec<F> = (0..N_padded).map(|i| {
            alpha_3 + beta_3 * F::from(i as u64) + eq_ry.coeffs[i]
        }).collect();
        let d4: Vec<F> = (0..N_padded).map(|i| {
            alpha_3 + beta_3 * index.col_poly.coeffs[i] + q_poly.coeffs[i]
        }).collect();

        // Batch-invert all four denominator vectors
        let inv_d1 = crate::utils::batch_invert::<F>(&d1);
        let inv_d2 = crate::utils::batch_invert::<F>(&d2);
        let inv_d3 = crate::utils::batch_invert::<F>(&d3);
        let inv_d4 = crate::utils::batch_invert::<F>(&d4);

        // g(x) = m_r(x)·d1⁻¹(x) − d2⁻¹(x) + ε·m_c(x)·d3⁻¹(x) − ε·d4⁻¹(x)
        let g: Vec<F> = (0..N_padded).map(|i| {
            m_r_poly.coeffs[i] * inv_d1[i]
            - inv_d2[i]
            + epsilon * m_c_poly.coeffs[i] * inv_d3[i]
            - epsilon * inv_d4[i]
        }).collect();
        let g_poly = MultilinearPoly::new(g);

        // Commit g̃ and add root to transcript
        let comm_g = TitanPolyCommitment::<G, F>::commit(prover_params_single, &g_poly);
        transcript.append_message(b"comm_g", &comm_g.root.to_sponge_bytes_as_vec());

        // Derive η (vector, for Eq 12) and δ (scalar, for batching)
        let eta: Vec<F> = (0..num_vars_idx)
            .map(|_| get_challenge::<F>(&mut transcript, b"eta"))
            .collect();
        let eq_eta = MultilinearPoly::<F>::init_with_eq(&eta);

        let delta: F = get_challenge::<F>(&mut transcript, b"delta");
        let delta2 = delta * delta;
        let delta3 = delta2 * delta;
        let delta4 = delta3 * delta;

        let d1_poly = MultilinearPoly::new(d1);
        let d2_poly = MultilinearPoly::new(d2);
        let d3_poly = MultilinearPoly::new(d3);
        let d4_poly = MultilinearPoly::new(d4);

        // Batched combine function (degree 6)
        // polys: [0:val_A, 1:val_B, 2:val_C, 3:p, 4:q, 5:eq_eta, 6:g,
        //         7:d1, 8:d2, 9:d3, 10:d4, 11:m_r, 12:m_c]
        let eps = epsilon;
        let combine_fn_3 = move |vals: &[F]| -> F {
            let pq = vals[3] * vals[4];

            // SC1-3: val_K · p · q
            let sc1 = vals[0] * pq;
            let sc2 = vals[1] * pq;
            let sc3 = vals[2] * pq;

            // SC4a: eq_η · (t − g·d1·d2·d3·d4)
            // where t = m_r·d2·d3·d4 − d1·d3·d4 + ε·m_c·d1·d2·d4 − ε·d1·d2·d3
            let (d1, d2, d3, d4) = (vals[7], vals[8], vals[9], vals[10]);
            let d1d2 = d1 * d2;
            let d3d4 = d3 * d4;
            let t = vals[11] * d2 * d3d4
                  - d1 * d3d4
                  + eps * vals[12] * d1d2 * d4
                  - eps * d1d2 * d3;
            let sc4a = vals[5] * (t - vals[6] * d1d2 * d3d4);

            // SC4b: g
            let sc4b = vals[6];

            sc1 + delta * sc2 + delta2 * sc3 + delta3 * sc4a + delta4 * sc4b
        };

        let (sumcheck_proof_3, challenges_3) = field_sumcheck::prove(
            vec![
                index.val_A_poly.clone(), index.val_B_poly.clone(), index.val_C_poly.clone(),
                p_poly.clone(), q_poly.clone(),
                eq_eta, g_poly.clone(),
                d1_poly, d2_poly, d3_poly, d4_poly,
                m_r_poly.clone(), m_c_poly.clone(),
            ],
            &combine_fn_3,
            6,
            &mut transcript,
        );

        // Oracle evaluations at challenge point r_3
        let r_3 = &challenges_3;
        let val_A_at_r3 = index.val_A_poly.evaluate(r_3);
        let val_B_at_r3 = index.val_B_poly.evaluate(r_3);
        let val_C_at_r3 = index.val_C_poly.evaluate(r_3);
        let p_at_r3 = p_poly.evaluate(r_3);
        let q_at_r3 = q_poly.evaluate(r_3);
        let g_at_r3 = g_poly.evaluate(r_3);
        let row_at_r3 = index.row_poly.evaluate(r_3);
        let col_at_r3 = index.col_poly.evaluate(r_3);
        let m_r_at_r3 = m_r_poly.evaluate(r_3);
        let m_c_at_r3 = m_c_poly.evaluate(r_3);
        println!("Time taken for Sumcheck 3: {}", start.elapsed().as_millis());

        // =====================================================================
        // Fourth sum-check: reduce z(r_y) and g(r_3) to common point r_4
        //   Σ_x [z(x)·eq(x,r_y) + ρ·g(x)·eq(x,r_3)] = z_at_ry + ρ·g_at_r3
        // =====================================================================

        let start = Instant::now();
        // Append third sum-check oracle evaluations to transcript
        append_to_transcript(&mut transcript, b"val_A_at_r3", &val_A_at_r3);
        append_to_transcript(&mut transcript, b"val_B_at_r3", &val_B_at_r3);
        append_to_transcript(&mut transcript, b"val_C_at_r3", &val_C_at_r3);
        append_to_transcript(&mut transcript, b"p_at_r3", &p_at_r3);
        append_to_transcript(&mut transcript, b"q_at_r3", &q_at_r3);
        append_to_transcript(&mut transcript, b"g_at_r3", &g_at_r3);
        append_to_transcript(&mut transcript, b"row_at_r3", &row_at_r3);
        append_to_transcript(&mut transcript, b"col_at_r3", &col_at_r3);
        append_to_transcript(&mut transcript, b"m_r_at_r3", &m_r_at_r3);
        append_to_transcript(&mut transcript, b"m_c_at_r3", &m_c_at_r3);

        // Commit witness polynomial z
        let comm_z = TitanPolyCommitment::<G, F>::commit(prover_params_single, &z_poly);
        transcript.append_message(b"comm_z", &comm_z.root.to_sponge_bytes_as_vec());

        let rho: F = get_challenge::<F>(&mut transcript, b"rho");

        let eq_ry_4 = MultilinearPoly::<F>::init_with_eq(r_y);
        let eq_r3_4 = MultilinearPoly::<F>::init_with_eq(r_3);

        let combine_fn_4 = move |vals: &[F]| -> F {
            vals[0] * vals[1] + rho * vals[2] * vals[3]
        };

        let (sumcheck_proof_4, challenges_4) = field_sumcheck::prove(
            vec![z_poly.clone(), eq_ry_4, g_poly.clone(), eq_r3_4],
            &combine_fn_4,
            2,
            &mut transcript,
        );

        let r_4 = &challenges_4;
        let z_at_r4 = z_poly.evaluate(r_4);
        let g_at_r4 = g_poly.evaluate(r_4);
        println!("Time taken for Sumcheck 4: {}", start.elapsed().as_millis());

        // =====================================================================
        // PCS evaluation proofs
        // =====================================================================

        // 1. Aggregate prove for p, q, m_r, m_c at r_3
        let start = Instant::now();
        let mut agg_pqm_coeffs: Vec<F> = [
            p_poly.coeffs.as_slice(),
            q_poly.coeffs.as_slice(),
            m_r_poly.coeffs.as_slice(),
            m_c_poly.coeffs.as_slice(),
        ].concat();
        agg_pqm_coeffs.resize(agg_pqm_coeffs.len().next_power_of_two(), F::ZERO);
        let agg_pqm_poly = MultilinearPoly::new(agg_pqm_coeffs);

        let agg_pqm_proof = TitanPolyCommitment::<G, F>::aggregate_prove(
            prover_params_agg,
            &agg_pqm_poly,
            &comm_pqm,
            r_3,
            4,
        );
        println!("Time for aggregate prove (p,q,m_r,m_c): {}", start.elapsed().as_millis());

        // 2. Batch prove for public index polynomials at r_3
        let start = Instant::now();
        let index_polys = vec![
            index.row_poly.clone(),
            index.col_poly.clone(),
            index.val_A_poly.clone(),
            index.val_B_poly.clone(),
            index.val_C_poly.clone(),
        ];
        let index_comms = vec![
            index.comm_row_poly.clone(),
            index.comm_col_poly.clone(),
            index.comm_val_A_poly.clone(),
            index.comm_val_B_poly.clone(),
            index.comm_val_C_poly.clone(),
        ];
        let batch_index_proof = TitanPolyCommitment::<G, F>::batch_prove(
            prover_params_trusted,
            &index_polys,
            &index_comms,
            r_3,
        );
        println!("Time for batch prove (index polys): {}", start.elapsed().as_millis());

        // 3. Batch prove for z, g at r_4
        let start = Instant::now();
        let batch_zg_proof = TitanPolyCommitment::<G, F>::batch_prove(
            prover_params_single,
            &[z_poly, g_poly],
            &[comm_z.clone(), comm_g.clone()],
            r_4,
        );
        println!("Time for batch prove (z, g): {}", start.elapsed().as_millis());

        SpartanProof {
            sumcheck_proof_1,
            a_at_rx,
            b_at_rx,
            c_at_rx,
            sumcheck_proof_2,
            eA_at_ry,
            eB_at_ry,
            eC_at_ry,
            z_at_ry,
            comm_pqm_root: comm_pqm.root.clone(),
            comm_g_root: comm_g.root.clone(),
            comm_z_root: comm_z.root.clone(),
            sumcheck_proof_3,
            val_A_at_r3,
            val_B_at_r3,
            val_C_at_r3,
            p_at_r3,
            q_at_r3,
            g_at_r3,
            row_at_r3,
            col_at_r3,
            m_r_at_r3,
            m_c_at_r3,
            sumcheck_proof_4,
            z_at_r4,
            g_at_r4,
            agg_pqm_proof,
            batch_index_proof,
            batch_zg_proof,
        }
    }

    pub fn verify(
        instance: &R1CS<F>,
        index: &TitanIndexPolys,
        verifier_params_single: &TitanVerifierParams<G, F>,
        verifier_params_trusted: &TitanVerifierParams<G, F>,
        verifier_params_agg: &TitanVerifierParams<G, F>,
        proof: &SpartanProof,
    ) -> bool {
        let num_constraints = instance.m;
        let m_padded = num_constraints.next_power_of_two();
        let num_vars_x = (m_padded as f64).log2() as usize;

        let n_padded = instance.n.next_power_of_two();
        let num_vars_y = (n_padded as f64).log2() as usize;

        // =====================================================================
        // Replay transcript identically to prover
        // =====================================================================
        let mut transcript = Transcript::new(b"spartan");
        transcript.append_u64(b"num_vars_x", num_vars_x as u64);
        transcript.append_u64(b"num_vars_y", num_vars_y as u64);

        let tau: Vec<F> = (0..num_vars_x)
            .map(|_| get_challenge::<F>(&mut transcript, b"tau"))
            .collect();

        // =====================================================================
        // Verify first sum-check: Σ_x eq(x,τ)·(Az(x)·Bz(x) − Cz(x)) = 0
        // =====================================================================
        let subclaim_1 = field_sumcheck::verify(
            &proof.sumcheck_proof_1,
            F::ZERO,
            num_vars_x,
            3,
            &mut transcript,
        );

        let subclaim_1 = match subclaim_1 {
            Ok(sc) => sc,
            Err(e) => {
                println!("First sum-check verification failed: {}", e);
                return false;
            }
        };

        // Oracle consistency check 1:
        // eq(r_x, τ) · (a_at_rx · b_at_rx − c_at_rx) == expected_evaluation
        let r_x = &subclaim_1.point;
        let mut eq_rx_tau = F::ONE;
        for j in 0..num_vars_x {
            eq_rx_tau *= r_x[j] * tau[j] + (F::ONE - r_x[j]) * (F::ONE - tau[j]);
        }

        let lhs_1 = eq_rx_tau * (proof.a_at_rx * proof.b_at_rx - proof.c_at_rx);
        if lhs_1 != subclaim_1.expected_evaluation {
            println!("Oracle check 1 failed: eq(r_x,τ)·(a·b−c) != expected");
            return false;
        }

        // =====================================================================
        // Verify second sum-check:
        //   Σ_y (eA(y) + β·eB(y) + β²·eC(y)) · z(y) = a_at_rx + β·b_at_rx + β²·c_at_rx
        // =====================================================================

        // Append first sum-check oracle evals to transcript, derive beta
        append_to_transcript(&mut transcript, b"a_at_rx", &proof.a_at_rx);
        append_to_transcript(&mut transcript, b"b_at_rx", &proof.b_at_rx);
        append_to_transcript(&mut transcript, b"c_at_rx", &proof.c_at_rx);
        let beta: F = get_challenge::<F>(&mut transcript, b"beta");
        let beta2 = beta * beta;

        let claimed_sum_2 = proof.a_at_rx + beta * proof.b_at_rx + beta2 * proof.c_at_rx;

        let subclaim_2 = field_sumcheck::verify(
            &proof.sumcheck_proof_2,
            claimed_sum_2,
            num_vars_y,
            2,
            &mut transcript,
        );

        let subclaim_2 = match subclaim_2 {
            Ok(sc) => sc,
            Err(e) => {
                println!("Second sum-check verification failed: {}", e);
                return false;
            }
        };

        // Oracle consistency check 2:
        // (eA_at_ry + β·eB_at_ry + β²·eC_at_ry) · z_at_ry == expected_evaluation
        let combined_at_ry = proof.eA_at_ry + beta * proof.eB_at_ry + beta2 * proof.eC_at_ry;
        let lhs_2 = combined_at_ry * proof.z_at_ry;
        if lhs_2 != subclaim_2.expected_evaluation {
            println!("Oracle check 2 failed: combined(r_y)·z(r_y) != expected");
            return false;
        }

        let r_y = &subclaim_2.point;

        // =====================================================================
        // Verify third batched sum-check (memory checking + matrix decomposition)
        // =====================================================================

        // Replay commitment roots to transcript
        transcript.append_message(b"comm_p_q_m", &proof.comm_pqm_root.to_sponge_bytes_as_vec());

        let alpha_3: F = get_challenge::<F>(&mut transcript, b"alpha");
        let beta_3: F = get_challenge::<F>(&mut transcript, b"beta_3");
        let epsilon: F = get_challenge::<F>(&mut transcript, b"epsilon");

        transcript.append_message(b"comm_g", &proof.comm_g_root.to_sponge_bytes_as_vec());

        let num_vars_idx = num_vars_x; // N_padded = m_padded in our setup
        let eta: Vec<F> = (0..num_vars_idx)
            .map(|_| get_challenge::<F>(&mut transcript, b"eta"))
            .collect();
        let delta: F = get_challenge::<F>(&mut transcript, b"delta");
        let delta2 = delta * delta;
        let delta3 = delta2 * delta;
        let delta4 = delta3 * delta;

        // Claimed sum = eA_at_ry + δ·eB_at_ry + δ²·eC_at_ry (SC4a, SC4b have sum 0)
        let claimed_sum_3 = proof.eA_at_ry + delta * proof.eB_at_ry + delta2 * proof.eC_at_ry;

        let subclaim_3 = field_sumcheck::verify(
            &proof.sumcheck_proof_3,
            claimed_sum_3,
            num_vars_idx,
            6,
            &mut transcript,
        );

        let subclaim_3 = match subclaim_3 {
            Ok(sc) => sc,
            Err(e) => {
                println!("Third sum-check verification failed: {}", e);
                return false;
            }
        };

        // Oracle consistency check 3
        let r_3 = &subclaim_3.point;

        // Helper: eq(a, b) = ∏_j (a_j·b_j + (1−a_j)·(1−b_j))
        let eq_eval = |a: &[F], b: &[F]| -> F {
            a.iter().zip(b.iter())
                .fold(F::ONE, |acc, (&ai, &bi)| {
                    acc * (ai * bi + (F::ONE - ai) * (F::ONE - bi))
                })
        };

        let eq_eta_at_r3 = eq_eval(r_3, &eta);
        let eq_rx_at_r3 = eq_eval(r_3, r_x);
        let eq_ry_at_r3 = eq_eval(r_3, r_y);

        // id(r_3) = Σ_j 2^j · r_3[j]
        let two = F::ONE + F::ONE;
        let mut id_at_r3 = F::ZERO;
        let mut pow2 = F::ONE;
        for j in 0..num_vars_idx {
            id_at_r3 += pow2 * r_3[j];
            pow2 *= two;
        }

        // Reconstruct d1..d4 at r_3
        let d1_r3 = alpha_3 + beta_3 * id_at_r3 + eq_rx_at_r3;
        let d2_r3 = alpha_3 + beta_3 * proof.row_at_r3 + proof.p_at_r3;
        let d3_r3 = alpha_3 + beta_3 * id_at_r3 + eq_ry_at_r3;
        let d4_r3 = alpha_3 + beta_3 * proof.col_at_r3 + proof.q_at_r3;

        // Evaluate the combine function from oracle values
        let pq = proof.p_at_r3 * proof.q_at_r3;
        let sc1 = proof.val_A_at_r3 * pq;
        let sc2 = proof.val_B_at_r3 * pq;
        let sc3 = proof.val_C_at_r3 * pq;

        let d1d2 = d1_r3 * d2_r3;
        let d3d4 = d3_r3 * d4_r3;
        let t_r3 = proof.m_r_at_r3 * d2_r3 * d3d4
                  - d1_r3 * d3d4
                  + epsilon * proof.m_c_at_r3 * d1d2 * d4_r3
                  - epsilon * d1d2 * d3_r3;
        let sc4a = eq_eta_at_r3 * (t_r3 - proof.g_at_r3 * d1d2 * d3d4);
        let sc4b = proof.g_at_r3;

        let lhs_3 = sc1 + delta * sc2 + delta2 * sc3 + delta3 * sc4a + delta4 * sc4b;
        if lhs_3 != subclaim_3.expected_evaluation {
            println!("Oracle check 3 failed");
            return false;
        }

        // =====================================================================
        // Verify fourth sum-check: reduce z(r_y) and g(r_3) to common point r_4
        //   Σ_x [z(x)·eq(x,r_y) + ρ·g(x)·eq(x,r_3)] = z_at_ry + ρ·g_at_r3
        // =====================================================================

        // Append third sum-check oracle evaluations to transcript
        append_to_transcript(&mut transcript, b"val_A_at_r3", &proof.val_A_at_r3);
        append_to_transcript(&mut transcript, b"val_B_at_r3", &proof.val_B_at_r3);
        append_to_transcript(&mut transcript, b"val_C_at_r3", &proof.val_C_at_r3);
        append_to_transcript(&mut transcript, b"p_at_r3", &proof.p_at_r3);
        append_to_transcript(&mut transcript, b"q_at_r3", &proof.q_at_r3);
        append_to_transcript(&mut transcript, b"g_at_r3", &proof.g_at_r3);
        append_to_transcript(&mut transcript, b"row_at_r3", &proof.row_at_r3);
        append_to_transcript(&mut transcript, b"col_at_r3", &proof.col_at_r3);
        append_to_transcript(&mut transcript, b"m_r_at_r3", &proof.m_r_at_r3);
        append_to_transcript(&mut transcript, b"m_c_at_r3", &proof.m_c_at_r3);

        // Replay witness commitment
        transcript.append_message(b"comm_z", &proof.comm_z_root.to_sponge_bytes_as_vec());

        let rho: F = get_challenge::<F>(&mut transcript, b"rho");

        let claimed_sum_4 = proof.z_at_ry + rho * proof.g_at_r3;
        let subclaim_4 = field_sumcheck::verify(
            &proof.sumcheck_proof_4,
            claimed_sum_4,
            num_vars_idx,
            2,
            &mut transcript,
        );

        let subclaim_4 = match subclaim_4 {
            Ok(sc) => sc,
            Err(e) => {
                println!("Fourth sum-check verification failed: {}", e);
                return false;
            }
        };

        let r_4 = &subclaim_4.point;

        // Oracle consistency check 4:
        // z(r_4)·eq(r_4,r_y) + ρ·g(r_4)·eq(r_4,r_3) == expected
        let eq_r4_ry = eq_eval(r_4, r_y);
        let eq_r4_r3 = eq_eval(r_4, r_3);
        let lhs_4 = proof.z_at_r4 * eq_r4_ry + rho * proof.g_at_r4 * eq_r4_r3;
        if lhs_4 != subclaim_4.expected_evaluation {
            println!("Oracle check 4 failed");
            return false;
        }

        // =====================================================================
        // PCS verification
        // =====================================================================

        // 1. Verify aggregate proof for p, q, m_r, m_c at r_3
        let res = TitanPolyCommitment::<G, F>::aggregate_verify(
            verifier_params_agg,
            &proof.comm_pqm_root,
            r_3,
            num_vars_idx + 2,
            4,
            &proof.agg_pqm_proof,
        );
        if !res {
            println!("Aggregate PCS verification for p,q,m_r,m_c failed");
            return false;
        }

        // Cross-check evaluations from aggregate proof against claimed oracle values
        if proof.agg_pqm_proof.evaluations[0] != proof.p_at_r3
            || proof.agg_pqm_proof.evaluations[1] != proof.q_at_r3
            || proof.agg_pqm_proof.evaluations[2] != proof.m_r_at_r3
            || proof.agg_pqm_proof.evaluations[3] != proof.m_c_at_r3
        {
            println!("Aggregate proof evaluations don't match claimed oracle values");
            return false;
        }

        // 2. Verify batch proof for index polys (row, col, val_A, val_B, val_C) at r_3
        let index_roots = vec![
            index.comm_row_poly.root.clone(),
            index.comm_col_poly.root.clone(),
            index.comm_val_A_poly.root.clone(),
            index.comm_val_B_poly.root.clone(),
            index.comm_val_C_poly.root.clone(),
        ];
        let index_sigmas = vec![
            proof.row_at_r3,
            proof.col_at_r3,
            proof.val_A_at_r3,
            proof.val_B_at_r3,
            proof.val_C_at_r3,
        ];
        let res = TitanPolyCommitment::<G, F>::batch_verify(
            verifier_params_trusted,
            &index_roots,
            r_3,
            &index_sigmas,
            &proof.batch_index_proof,
        );
        if !res {
            println!("Batch PCS verification for index polys failed");
            return false;
        }

        // 3. Verify batch proof for z, g at r_4
        let zg_roots = vec![proof.comm_z_root.clone(), proof.comm_g_root.clone()];
        let zg_sigmas = vec![proof.z_at_r4, proof.g_at_r4];
        let res = TitanPolyCommitment::<G, F>::batch_verify(
            verifier_params_single,
            &zg_roots,
            r_4,
            &zg_sigmas,
            &proof.batch_zg_proof,
        );
        if !res {
            println!("Batch PCS verification for z, g failed");
            return false;
        }

        true
    }
}

mod tests {
    use std::cmp::max;
    use super::*;
    use crate::r1cs::{generate_random_r1cs, verify_r1cs};
    use std::time::Instant;

    #[test]
    fn test_spartan_titan_index() {
        let r1cs_m = 1usize << 20;
        let r1cs_n = 1usize << 20;
        let max_nonzeros = 1usize << 22;

        // Generate a random satisfiable R1CS instance
        let mut r1cs = generate_random_r1cs::<F>(r1cs_m, r1cs_n, max_nonzeros);
        r1cs.witness.resize(max_nonzeros, F::ZERO);
        assert!(verify_r1cs(&r1cs), "R1CS instance should be satisfiable");
        println!(
            "R1CS: {} constraints, {} variables, A={} B={} C={} entries",
            r1cs.m, r1cs.n,
            r1cs.a.entries.len(), r1cs.b.entries.len(), r1cs.c.entries.len()
        );

        // Index polynomials are extended to size max_nonzeros (padded to power of 2)
        let N_padded = max_nonzeros.next_power_of_two();
        let poly_vars = (N_padded as f64).log2() as usize;
        println!("Index poly size: max_nonzeros={}, N_padded={}, num_vars={}", max_nonzeros, N_padded, poly_vars);

        // Setup PCS with m = poly_vars (size of index polynomials)
        let poly_vars_extended = poly_vars + 2;
        let m1 = poly_vars/2 - 2;
        let m1_ext = (poly_vars_extended / 2) - 2;
        let l1 = 1usize;
        let l1_ext = 2usize;
        let m2 = poly_vars - m1;
        let m2_ext = poly_vars_extended - m1_ext;
        let l2 = m2/2 - 1;
        let l2_ext = m2_ext/2 - 1;
        let config = TitanSetupConfig {
            m: poly_vars,
            m1,
            l1,
            domain_g1_size: m1 - l1 + 4,
            num_queries: 64,
            l2,
            domain_g2_size: m2 - l2 + 4,
            num_merkle_nodes: 0,
        };

        let mut config_agg = config.clone();
        let mut config_trusted = config.clone();
        config_trusted.num_queries = 32;

        config_agg.m += 2;
        config_agg.m1 = m1_ext;
        config_agg.l1 = l1_ext;
        config_agg.domain_g1_size = m1_ext - l1_ext + 4;
        config_agg.l2 = l2_ext;
        config_agg.domain_g2_size = m2_ext - l2_ext + 4;

        let snark_config = TitanSnarkConfig {single_oracle_config: config, trusted_oracle_config: config_trusted, aggregate_oracle_config: config_agg};

        let start = Instant::now();
        let pp_vec = TitanSnark::setup(&snark_config);
        let prover_params_single = &pp_vec[0].0;
        let prover_params_trusted = &pp_vec[1].0;
        let prover_params_agg = &pp_vec[2].0;
        let verifier_params_single = &pp_vec[0].1;
        let verifier_params_trusted = &pp_vec[1].1;
        let verifier_params_agg = &pp_vec[2].1;

        println!("Setup time: {} msec", start.elapsed().as_millis());

        // Index
        let start = Instant::now();
        let index = TitanSnark::index(&r1cs, prover_params_trusted, max_nonzeros);
        println!("Index time: {} msec", start.elapsed().as_millis());

        // Verify index polynomials have correct size
        assert_eq!(index.row_poly.num_vars, poly_vars);
        assert_eq!(index.col_poly.num_vars, poly_vars);
        assert_eq!(index.val_A_poly.num_vars, poly_vars);
        assert_eq!(index.val_B_poly.num_vars, poly_vars);
        assert_eq!(index.val_C_poly.num_vars, poly_vars);

        // Prove
        let start = Instant::now();
        let proof = TitanSnark::prove(
            &r1cs, &index,
            prover_params_single,
            prover_params_trusted,
            prover_params_agg,
            &r1cs.witness);
        println!("Prove time: {} msec", start.elapsed().as_millis());

        // Verify
        let start = Instant::now();
        let res = TitanSnark::verify(
            &r1cs, &index,
            verifier_params_single,
            verifier_params_trusted,
            verifier_params_agg,
            &proof,
        );
        println!("Verify time: {} msec", start.elapsed().as_millis());
        assert!(res, "Spartan proof verification failed");

        println!("Proof Size = {} KB", proof.get_proof_size()/32);
    }
}

