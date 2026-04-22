use crate::traits::*;
use ff::BitViewSized;

#[derive(Clone, Debug)]
pub struct MultilinearPoly<C> {
    pub coeffs: Vec<C>,  // Coefficients of the multilinear polynomial
    pub num_vars: usize, // Number of variables (n)
}

impl<F> MultilinearPoly<F> {
    // Create a new multilinear polynomial with given coefficients
    pub fn new(coeffs: Vec<F>) -> Self {
        let num_vars = (coeffs.len() as f64).log2() as usize;
        assert_eq!(
            coeffs.len(),
            1 << num_vars,
            "Coefficient length must be a power of 2"
        );
        Self { coeffs, num_vars }
    }
}

// This block specifies operations on multilinear polynomials over F-vector space C
// Some functions require an additional trait bound InnerProduct<F, Output = C> to be supported.
// which specifies that inner product can be defined for vector over C and vector over F.
impl<C> MultilinearPoly<C> {
    pub fn init_with_eq(alpha: &[C]) -> Self
    where
        C: ff::Field,
    {
        let num_vars = alpha.len();
        let mut coeffs: Vec<C> = vec![C::ONE; 1 << num_vars]; // 2^n coefficients initialized to 1
        let mut filled:usize = 1;
        let mut var_idx: usize = 0;
        while filled < coeffs.len() {
            for i in 0..filled {
                coeffs[i+filled] = alpha[var_idx] * coeffs[i];
            }
            for i in 0..filled {
                coeffs[i] = coeffs[i] * (C::ONE - alpha[var_idx]);
            }
            filled += filled;
            var_idx += 1;
        }

        Self { coeffs, num_vars }
    }

    // folds the multilinear polynomial by substituting the last variable with r
    pub fn fold<F>(&mut self, r: F)
    where
        C: Linear<F>,
        F: ff::Field,
    {
        let n = self.coeffs.len();
        assert_eq!(
            n,
            1 << self.num_vars,
            "Coefficient length must be a power of 2"
        );
        let mut new_coeffs = vec![C::zero(); n / 2];
        for i in 0..(n / 2) {
            new_coeffs[i] = self.coeffs[i] * (F::ONE - r) + self.coeffs[i + n / 2] * r;
        }
        self.coeffs = new_coeffs;
        self.num_vars -= 1;
    }

    // folds the multilinear polynomial by substituting the first variable with r
    pub fn fold_first<F>(&mut self, r: F)
    where
        C: Linear<F>,
        F: ff::Field,
    {
        let n = self.coeffs.len();
        assert_eq!(n, 1 << self.num_vars, "Size not a power of 2");

        if self.num_vars == 0 {
            return;
        }

        let mut new_coeffs = vec![C::zero(); n / 2];
        for i in 0..(n / 2) {
            new_coeffs[i] = self.coeffs[2 * i] * (F::ONE - r) + self.coeffs[1 + 2 * i] * r;
        }
        self.coeffs = new_coeffs;
        self.num_vars -= 1;
    }

    // This function returns restriction of a multilinear polynomial
    // rho: fixes the values of first rho.len() variables to vector rho
    // output: f(rho,x) where f is the multilinear polynomial
    pub fn restrict<F>(&self, rho: &[F]) -> Self
    where
        C: Linear<F> + InnerProduct<F, Output = C>,
        C::Affine: InnerProduct<F, Output = C>,
        F: ff::Field,
    {
        assert_eq!(
            rho.len() < self.num_vars,
            true,
            "rho length must be less than number of variables"
        );
        let num_vars_restricted = self.num_vars - rho.len();
        let mut coeffs = vec![C::zero(); 1 << num_vars_restricted];
        let affine_coeffs = C::to_affine(&self.coeffs);
        let restricted_poly_size = 1 << num_vars_restricted;
        let eq_poly = MultilinearPoly::init_with_eq(rho);
        let slice_size = 1 << (rho.len());
        for i in 0..restricted_poly_size {
            coeffs[i] = C::Affine::inner_product_msm(
                &affine_coeffs[slice_size * i..slice_size * (i + 1)],
                &eq_poly.coeffs,
            );
            /*
            coeffs[i] = eq_poly
                .coeffs
                .iter()
                .zip(self.coeffs.iter().skip((1 << rho.len()) * i))
                .map(|(a, b)| *b * *a)
                .fold(C::zero(), |acc, x| acc + x);

             */
        }
        Self {
            coeffs,
            num_vars: num_vars_restricted,
        }
    }

    pub fn fold_msb<F>(&self, rho: &[F]) -> Self
    where C: Linear<F>, F: ff::Field,
    {
        let mut poly = MultilinearPoly::new(self.coeffs.clone());
        for i in 0..rho.len() {
            poly.fold(rho[rho.len() - 1 - i]);
        }
        poly
    }

    // This function returns restriction of a multilinear polynomial
    // rho: fixes the values of first rho.len() variables to vector rho
    // output: f(x,rho) where f is the multilinear polynomial
    pub fn restrict_msb<F>(&self, rho: &[F]) -> Self
    where
        C: Linear<F> + InnerProduct<F, Output = C>,
        //C::Affine: InnerProduct<F, Output=C>,
        F: ff::Field,
    {
        assert_eq!(
            rho.len() <= self.num_vars,
            true,
            "rho length must be less than number of variables"
        );
        let num_vars_restricted = self.num_vars - rho.len();
        let mut coeffs = vec![C::zero(); 1 << num_vars_restricted];
        let restricted_poly_size = 1 << num_vars_restricted;
        let eq_poly = MultilinearPoly::init_with_eq(rho);
        let slice_size = 1 << (rho.len());
        // we first build a flattened version of transposed polynomial for convenience.
        let mut flattened_coeffs = vec![C::zero(); restricted_poly_size * slice_size];
        for x in 0..restricted_poly_size {
            for b in 0..slice_size {
                flattened_coeffs[x * slice_size + b] = self.coeffs[b * restricted_poly_size + x];
            }
        }

        for i in 0..restricted_poly_size {
            coeffs[i] = C::inner_product_msm(
                &flattened_coeffs[slice_size * i..slice_size * (i + 1)],
                &eq_poly.coeffs,
            );
            /*
            coeffs[i] = eq_poly
                .coeffs
                .iter()
                .zip(self.coeffs.iter().skip((1 << rho.len()) * i))
                .map(|(a, b)| *b * *a)
                .fold(C::zero(), |acc, x| acc + x);

             */
        }
        Self {
            coeffs,
            num_vars: num_vars_restricted,
        }
    }



    // Evaluate the multilinear polynomial at a point
    pub fn evaluate<F>(&self, point: &[F]) -> C
    where
        C: Linear<F>,
        F: ff::Field,
    {
        assert_eq!(
            point.len(),
            self.num_vars,
            "Point dimension must match number of variables"
        );
       let mut poly = MultilinearPoly::new(self.coeffs.clone());
       for p in point {
        poly.fold_first(*p);
       };
       poly.coeffs[0]

    }

    // Evaluate using msm. Requires InnerProduct trait to be supported
    pub fn evaluate_msm<F>(&self, point: &[F]) -> C
    where
        C: InnerProduct<F, Output = C>,
        F: ff::Field,
    {
        assert_eq!(
            point.len(),
            self.num_vars,
            "Point dimension must match number of variables"
        );
        let eq_poly = MultilinearPoly::init_with_eq(point);
        C::inner_product_msm(&self.coeffs, &eq_poly.coeffs)
    }
}

mod tests {
    use crate::rand_scalar;
    use pasta_curves::group::ff::Field;
    use pasta_curves::group::Group;
    use pasta_curves::pallas::Point as G;
    use pasta_curves::pallas::Scalar as F;
    use rand::prelude::StdRng;
    use rand::SeedableRng;

    use super::*;
    #[test]
    fn test_init_with_eq() {
        let mut rng = StdRng::from_entropy();

        let p_alpha = MultilinearPoly::init_with_eq(&vec![F::from(3u64), F::from(2u64)]);
        // eq((a1,a2),(x1,x2)) = ((1-a1)*(1-x1) + a1*x1)) * ((1-a2)*(1-x2) + a2*x2))
        // p_alpha(x1,x2) = (5x1 - 2)*(3x2 - 1)
        let x1 = rand_scalar(&mut rng);
        let x2 = rand_scalar(&mut rng);
        let eval = p_alpha.evaluate(&vec![x1, x2]);
        let eval_ref = (F::from(5u64) * x1 - F::from(2u64)) * (F::from(3u64) * x2 - F::one());
        assert_eq!(eval, eval_ref);
    }

    #[test]
    fn test_restrict() {
        let mut rng = StdRng::from_entropy();
        let num_vars: usize = 4;
        let coeffs: Vec<F> = (0..(1 << num_vars))
            .map(|_| rand_scalar(&mut rng))
            .collect();
        let poly_p = MultilinearPoly::new(coeffs);
        let r1 = rand_scalar(&mut rng);
        let r2 = rand_scalar(&mut rng);
        let poly_q = poly_p.restrict(&vec![r1, r2]);
        assert_eq!(poly_q.num_vars, num_vars - 2);

        let x1 = rand_scalar(&mut rng);
        let x2 = rand_scalar(&mut rng);
        let eval_p = poly_p.evaluate(&vec![r1, r2, x1, x2]);
        let eval_q = poly_q.evaluate(&vec![x1, x2]);
        assert_eq!(eval_p, eval_q);

        let poly_r = poly_p.restrict_msb(&vec![r2, x1, x2]);
        assert_eq!(poly_r.num_vars, num_vars - 3);
        let eval_r = poly_r.evaluate(&vec![r1]);
        assert_eq!(eval_p, eval_r);
    }

    #[test]
    fn test_restrict_group_poly() {
        let mut rng = StdRng::from_entropy();
        let num_vars: usize = 4;
        let coeffs: Vec<F> = (0..(1 << num_vars))
            .map(|_| rand_scalar(&mut rng))
            .collect();
        // reuse coeffs as group elements
        let coeffs: Vec<G> = coeffs.iter().map(|x| G::generator() * x).collect();
        let poly_p = MultilinearPoly::new(coeffs);
        let r1 = rand_scalar(&mut rng);
        let r2 = rand_scalar(&mut rng);
        let x1 = rand_scalar(&mut rng);
        let x2 = rand_scalar(&mut rng);
        let poly_q = poly_p.restrict(&vec![r1, r2]);
        assert_eq!(poly_q.num_vars, num_vars - 2);
        let eval_p = poly_p.evaluate(&vec![r1, r2, x1, x2]);
        let eval_q = poly_q.evaluate(&vec![x1, x2]);
        assert_eq!(eval_p, eval_q);
    }

    #[test]
    fn test_folding() {
        let mut rng = StdRng::from_entropy();
        let num_vars: usize = 4;
        let coeffs: Vec<F> = (0..(1 << num_vars))
            .map(|_| rand_scalar(&mut rng))
            .collect();
        let mut poly_p = MultilinearPoly::new(coeffs);
        let r1 = rand_scalar(&mut rng);
        let r2 = rand_scalar(&mut rng);
        let poly_q = poly_p.restrict(&vec![r1, r2]);
        poly_p.fold_first(r1);
        poly_p.fold_first(r2);

        assert_eq!(poly_q.coeffs, poly_p.coeffs);
    }
}
