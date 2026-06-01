//! Cup product: protocol composition as H¹ × H¹ → H².
//!
//! The cup product in sheaf cohomology models how protocol compositions
//! interact: combining two protocols P and Q yields a new protocol whose
//! obstruction theory lives in H².

use serde::{Deserialize, Serialize};

use crate::cohomology::Cohomology;

/// Cup product computation for sheaf cohomology.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CupProduct {
    /// The bilinear form defining the cup product.
    product_matrix: Vec<Vec<f64>>,
    h1_dim_a: usize,
    h1_dim_b: usize,
    h2_dim: usize,
}

impl CupProduct {
    /// Create a cup product from two H¹ spaces into H².
    pub fn new(h1_dim_a: usize, h1_dim_b: usize, h2_dim: usize) -> Self {
        // The cup product is a bilinear map H¹(A) × H¹(B) → H²
        // Represented as a matrix of shape (h1_a × h1_b) → h2
        let total_pairs = h1_dim_a * h1_dim_b;
        let product_matrix = if total_pairs > 0 && h2_dim > 0 {
            vec![vec![0.0; total_pairs]; h2_dim]
        } else {
            vec![]
        };
        CupProduct { product_matrix, h1_dim_a, h1_dim_b, h2_dim }
    }

    /// Set the cup product coefficients.
    pub fn set_coefficients(&mut self, coeffs: Vec<Vec<f64>>) {
        self.product_matrix = coeffs;
    }

    /// Compute the cup product of two obstruction classes.
    pub fn compute(&self, a: &[f64], b: &[f64]) -> Vec<f64> {
        if self.h2_dim == 0 || self.product_matrix.is_empty() {
            return vec![];
        }
        // Tensor product a ⊗ b
        let tensor: Vec<f64> = (0..self.h1_dim_a)
            .flat_map(|i| (0..self.h1_dim_b).map(move |j| a.get(i).copied().unwrap_or(0.0) * b.get(j).copied().unwrap_or(0.0)))
            .collect();
        // Apply product matrix
        (0..self.h2_dim)
            .map(|i| {
                self.product_matrix[i].iter().zip(tensor.iter())
                    .map(|(m, t)| m * t)
                    .sum()
            })
            .collect()
    }

    /// Compose two protocols: H¹(P) × H¹(Q) → H²(P ∪ Q).
    pub fn compose_protocols(
        cohomology_a: &Cohomology,
        cohomology_b: &Cohomology,
    ) -> (CupProduct, Cohomology) {
        let h1a = cohomology_a.h1_dimension();
        let h1b = cohomology_b.h1_dimension();
        // For composition, H² of the combined system is bounded by
        // the cup product structure
        let h2_combined = h1a * h1b; // Upper bound
        let mut cup = CupProduct::new(h1a, h1b, h2_combined);

        // Default: anti-commutative cup product ⌣(α, β) = -⌣(β, α)
        for i in 0..h1a {
            for j in 0..h1b {
                let idx = i * h1b + j;
                let target = i * h1b + j; // Map to H² component
                if target < h2_combined {
                    cup.product_matrix[target][idx] = 1.0;
                }
            }
        }

        // Combined cohomology
        let h0_combined = cohomology_a.h0_dimension() + cohomology_b.h0_dimension() - 1;
        let h1_combined = cohomology_a.h1_dimension() + cohomology_b.h1_dimension();
        let combined = Cohomology::with_dimensions(
            h0_combined.max(1),
            h1_combined,
            h2_combined,
        );

        (cup, combined)
    }

    /// Check if the cup product of two obstruction classes vanishes.
    pub fn vanishes(&self, a: &[f64], b: &[f64]) -> bool {
        let result = self.compute(a, b);
        result.iter().all(|&x| x.abs() < 1e-10)
    }

    /// Compute all cup products of basis elements.
    pub fn product_table(&self) -> Vec<Vec<Vec<f64>>> {
        (0..self.h1_dim_a)
            .map(|i| {
                let a = {
                    let mut v = vec![0.0; self.h1_dim_a];
                    v[i] = 1.0;
                    v
                };
                (0..self.h1_dim_b)
                    .map(|j| {
                        let b = {
                            let mut v = vec![0.0; self.h1_dim_b];
                            v[j] = 1.0;
                            v
                        };
                        self.compute(&a, &b)
                    })
                    .collect()
            })
            .collect()
    }

    /// The anti-commutativity check: ⌣(α, β) = -⌣(β, α) in cohomology.
    pub fn check_anticommutativity(&self) -> bool {
        for i in 0..self.h1_dim_a {
            for j in 0..self.h1_dim_b {
                let mut a = vec![0.0; self.h1_dim_a];
                a[i] = 1.0;
                let mut b = vec![0.0; self.h1_dim_b];
                b[j] = 1.0;

                let ab = self.compute(&a, &b);
                let ba = self.compute(&b, &a);

                // ⌣(α, β) + ⌣(β, α) should be 0 (up to sign)
                for k in 0..ab.len().max(ba.len()) {
                    let v1 = ab.get(k).copied().unwrap_or(0.0);
                    let v2 = ba.get(k).copied().unwrap_or(0.0);
                    if (v1 + v2).abs() > 1e-10 {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Dimensions of the cup product.
    pub fn dimensions(&self) -> (usize, usize, usize) {
        (self.h1_dim_a, self.h1_dim_b, self.h2_dim)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cup_product_creation() {
        let cp = CupProduct::new(2, 2, 4);
        assert_eq!(cp.dimensions(), (2, 2, 4));
    }

    #[test]
    fn test_cup_product_compute() {
        let mut cp = CupProduct::new(2, 2, 1);
        cp.set_coefficients(vec![vec![1.0, 0.0, 0.0, 1.0]]);
        let a = vec![1.0, 0.0];
        let b = vec![0.0, 1.0];
        let result = cp.compute(&a, &b);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_cup_product_vanishes() {
        let cp = CupProduct::new(2, 2, 0);
        assert!(cp.vanishes(&[1.0, 0.0], &[0.0, 1.0]));
    }

    #[test]
    fn test_compose_protocols() {
        let ca = Cohomology::with_dimensions(1, 1, 0);
        let cb = Cohomology::with_dimensions(1, 1, 0);
        let (cup, combined) = CupProduct::compose_protocols(&ca, &cb);
        assert_eq!(combined.h1_dimension(), 2);
        assert_eq!(combined.h2_dimension(), 1);
    }

    #[test]
    fn test_compose_deadlock_free() {
        let ca = Cohomology::trivial();
        let cb = Cohomology::trivial();
        let (cup, combined) = CupProduct::compose_protocols(&ca, &cb);
        assert!(combined.is_deadlock_free());
        assert_eq!(cup.h2_dim, 0);
    }

    #[test]
    fn test_product_table() {
        let mut cp = CupProduct::new(1, 1, 1);
        cp.set_coefficients(vec![vec![1.0]]);
        let table = cp.product_table();
        assert_eq!(table.len(), 1);
        assert_eq!(table[0].len(), 1);
    }

    #[test]
    fn test_bilinear() {
        let mut cp = CupProduct::new(2, 2, 1);
        cp.set_coefficients(vec![vec![1.0, 1.0, 1.0, 1.0]]);
        // ⌣(α₁ + α₂, β) should equal ⌣(α₁, β) + ⌣(α₂, β)
        let a_sum = vec![1.0, 1.0];
        let a1 = vec![1.0, 0.0];
        let a2 = vec![0.0, 1.0];
        let b = vec![1.0, 0.0];
        let sum_result = cp.compute(&a_sum, &b);
        let a1_result = cp.compute(&a1, &b);
        let a2_result = cp.compute(&a2, &b);
        assert_eq!(sum_result.len(), 1);
        assert!((sum_result[0] - (a1_result[0] + a2_result[0])).abs() < 1e-10);
    }
}
