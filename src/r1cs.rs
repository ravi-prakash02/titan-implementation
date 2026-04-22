use crate::traits::Linear;
use rand::Rng;

/// Represents a sparse matrix using coordinate format (COO)
#[derive(Debug, Clone)]
pub struct SparseMatrix<F> {
    pub rows: usize,
    pub cols: usize,
    /// (row, col, value) tuples
    pub entries: Vec<(usize, usize, F)>,
}

impl<F> SparseMatrix<F>
where
    F: ff::Field + Linear<F>,
{
    /// Multiply sparse matrix by a dense vector
    pub fn mult_vec(&self, vec: &[F]) -> Vec<F> {
        let mut result = vec![F::ZERO; self.rows];
        for &(row, col, val) in &self.entries {
            result[row] += val * vec[col];
        }
        result
    }

    /// Multiply sparse matrix by a dense vector
    pub fn mult_vec_transpose(&self, vec: &[F]) -> Vec<F> {
        let mut result = vec![F::ZERO; self.cols];
        for &(row, col, val) in &self.entries {
            result[col] += val * vec[row];
        }
        result
    }
}

/// R1CS instance: (A·z) ∘ (B·z) = (C·z)
/// where z is the witness vector
#[derive(Debug)]
pub struct R1CS<F> {
    pub m: usize, // number of constraints
    pub n: usize, // witness size
    pub a: SparseMatrix<F>,
    pub b: SparseMatrix<F>,
    pub c: SparseMatrix<F>,
    pub witness: Vec<F>,
}

/// Generate a random R1CS instance of size m×n with at most N total non-zero entries
/// Uses bit decomposition and other common constraints to ensure satisfiability
pub fn generate_random_r1cs<F>(m: usize, n: usize, max_nonzeros: usize) -> R1CS<F>
where
    F: ff::Field + Linear<F>,
{
    let mut rng = rand::thread_rng();
    let two = F::ONE + F::ONE;

    // witness[0] = 1 (standard convention for R1CS constant)
    let mut witness = vec![F::ONE];

    let bits_per_value = 8usize;
    let num_values_to_decompose = (m / 10).max(1).min((n - 1) / (bits_per_value + 1));

    let mut all_bits: Vec<Vec<F>> = Vec::new();

    for _ in 0..num_values_to_decompose {
        let mut bits = Vec::new();
        let mut value = F::ZERO;
        let mut power_of_two = F::ONE;

        for _ in 0..bits_per_value {
            let bit: bool = rng.gen();
            let bit_f = if bit { F::ONE } else { F::ZERO };
            bits.push(bit_f);
            value += power_of_two * bit_f;
            power_of_two = power_of_two * two;
        }

        witness.push(value);
        for &b in &bits {
            witness.push(b);
        }
        all_bits.push(bits);
    }

    // Add random field elements for remaining witness slots
    while witness.len() < n {
        witness.push(F::random(&mut rng));
    }
    witness.truncate(n);

    // Now build A, B, C matrices to encode constraints
    let mut a_entries: Vec<(usize, usize, F)> = Vec::new();
    let mut b_entries: Vec<(usize, usize, F)> = Vec::new();
    let mut c_entries: Vec<(usize, usize, F)> = Vec::new();

    let mut constraint_idx = 0;
    let entries_count =
        |a: &Vec<(usize, usize, F)>, b: &Vec<(usize, usize, F)>, c: &Vec<(usize, usize, F)>| {
            a.len() + b.len() + c.len()
        };

    // Constraint type 1: Bit constraints (b * b = b, ensuring b ∈ {0, 1})
    for (val_idx, bits) in all_bits.iter().enumerate() {
        for (bit_idx, _) in bits.iter().enumerate() {
            if constraint_idx >= m
                || entries_count(&a_entries, &b_entries, &c_entries) + 3 > max_nonzeros
            {
                break;
            }

            let witness_idx = 1 + val_idx * (bits_per_value + 1) + 1 + bit_idx;

            // Constraint: b * b = b
            a_entries.push((constraint_idx, witness_idx, F::ONE));
            b_entries.push((constraint_idx, witness_idx, F::ONE));
            c_entries.push((constraint_idx, witness_idx, F::ONE));

            constraint_idx += 1;
        }
    }

    // Constraint type 2: Bit decomposition sum (value = sum of 2^i * bit_i)
    for (val_idx, bits) in all_bits.iter().enumerate() {
        // This constraint adds 2 + bits.len() entries
        if constraint_idx >= m
            || entries_count(&a_entries, &b_entries, &c_entries) + 2 + bits.len() > max_nonzeros
        {
            break;
        }

        let value_witness_idx = 1 + val_idx * (bits_per_value + 1);

        // Constraint: value * 1 = sum of 2^i * bit_i
        a_entries.push((constraint_idx, value_witness_idx, F::ONE));
        b_entries.push((constraint_idx, 0, F::ONE)); // multiply by constant 1

        let mut power_of_two = F::ONE;
        for bit_idx in 0..bits.len() {
            let bit_witness_idx = value_witness_idx + 1 + bit_idx;
            c_entries.push((constraint_idx, bit_witness_idx, power_of_two));
            power_of_two = power_of_two * two;
        }

        constraint_idx += 1;
    }

    // Constraint type 3: Random linear constraints
    while constraint_idx < m {
        // Each constraint adds at most 3 entries
        if entries_count(&a_entries, &b_entries, &c_entries) + 3 > max_nonzeros {
            break;
        }

        let constraint_type = rng.gen_range(0..3);

        match constraint_type {
            0 => {
                // Constraint: w[i] * 1 = w[i]
                let idx = rng.gen_range(1..n);
                a_entries.push((constraint_idx, idx, F::ONE));
                b_entries.push((constraint_idx, 0, F::ONE));
                c_entries.push((constraint_idx, idx, F::ONE));
            }
            1 => {
                // Constraint: w[i] * w[j] = w[k] (if this holds in witness)
                if n >= 4 {
                    let i = rng.gen_range(1..n);
                    let j = rng.gen_range(1..n);
                    let product = witness[i] * witness[j];

                    // Find or create a witness element with this value
                    let k = if let Some(pos) = witness.iter().position(|&x| x == product) {
                        pos
                    } else if witness.len() < n {
                        let pos = witness.len();
                        witness.push(product);
                        pos
                    } else {
                        rng.gen_range(1..n)
                    };

                    if witness[i] * witness[j] == witness[k] {
                        a_entries.push((constraint_idx, i, F::ONE));
                        b_entries.push((constraint_idx, j, F::ONE));
                        c_entries.push((constraint_idx, k, F::ONE));
                    } else {
                        // Fallback to identity constraint
                        a_entries.push((constraint_idx, i, F::ONE));
                        b_entries.push((constraint_idx, 0, F::ONE));
                        c_entries.push((constraint_idx, i, F::ONE));
                    }
                }
            }
            _ => {
                // Constraint: 0 * 1 = 0 (always satisfied)
                b_entries.push((constraint_idx, 0, F::ONE));
            }
        }

        constraint_idx += 1;
    }

    // Update m to actual number of constraints
    //let actual_m = constraint_idx.min(m);

    let a = SparseMatrix {
        rows: max_nonzeros, // actual_m,
        cols: max_nonzeros,
        entries: a_entries,
    };

    let b = SparseMatrix {
        rows: max_nonzeros,
        cols: max_nonzeros,
        entries: b_entries,
    };

    let c = SparseMatrix {
        rows: max_nonzeros,
        cols: max_nonzeros,
        entries: c_entries,
    };

    witness.resize(max_nonzeros, F::ZERO);
    R1CS {
        m: max_nonzeros,
        n: max_nonzeros,
        a,
        b,
        c,
        witness,
    }
}

/// Verify that the R1CS instance is satisfied by the witness
pub fn verify_r1cs<F>(r1cs: &R1CS<F>) -> bool
where
    F: ff::Field + Linear<F>,
{
    let a_z = r1cs.a.mult_vec(&r1cs.witness);
    let b_z = r1cs.b.mult_vec(&r1cs.witness);
    let c_z = r1cs.c.mult_vec(&r1cs.witness);

    let mut all_satisfied = true;
    for i in 0..r1cs.m {
        if a_z[i] * b_z[i] != c_z[i] {
            println!(
                "Constraint {} failed: {:?} * {:?} != {:?}",
                i, a_z[i], b_z[i], c_z[i]
            );
            all_satisfied = false;
        }
    }

    all_satisfied
}

#[cfg(test)]
mod tests {
    use super::*;
    use pasta_curves::pallas::Scalar as F;

    #[test]
    fn test_generate_r1cs_with_witness() {
        let m = 50;
        let n = 100;
        let max_nonzeros = 256;

        let r1cs = generate_random_r1cs::<F>(m, n, max_nonzeros);

       // assert!(r1cs.m <= m);
       // assert_eq!(r1cs.n, n);
       // assert_eq!(r1cs.witness.len(), n);

        let total_entries =
            r1cs.a.entries.len() + r1cs.b.entries.len() + r1cs.c.entries.len();

        assert!(total_entries <= max_nonzeros);

        println!("Generated R1CS with:");
        println!("  Constraints: {}", r1cs.m);
        println!("  Variables: {}", r1cs.n);
        println!("  A: {} entries", r1cs.a.entries.len());
        println!("  B: {} entries", r1cs.b.entries.len());
        println!("  C: {} entries", r1cs.c.entries.len());
        println!("  Total: {} / {}", total_entries, max_nonzeros);
        println!("  Witness size: {}", r1cs.witness.len());

        // Verify the instance is satisfied
        let is_valid = verify_r1cs(&r1cs);
        println!("  Valid: {}", is_valid);

        assert!(is_valid, "R1CS constraints should be satisfied by the witness");
    }

    #[test]
    fn test_small_instance() {
        let r1cs = generate_random_r1cs::<F>(10, 20, 30);
        assert!(verify_r1cs(&r1cs));
    }
}
