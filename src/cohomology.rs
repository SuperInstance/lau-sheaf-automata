//! Sheaf cohomology: H⁰, H¹ computation and obstruction theory.
//!
//! Kimi's Theorem 2: Protocol P is deadlock-free ⟺ H¹(Sh(States); L(P)) = 0.

use serde::{Deserialize, Serialize};

/// Sheaf cohomology groups for a protocol space.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Cohomology {
    /// H⁰ dimension (global sections).
    h0_dim: usize,
    /// H¹ dimension (obstruction classes = deadlocks).
    h1_dim: usize,
    /// H² dimension (for cup product target).
    h2_dim: usize,
    /// The coboundary map d0: C⁰ → C¹.
    d0: Vec<Vec<f64>>,
    /// The coboundary map d1: C¹ → C².
    d1: Vec<Vec<f64>>,
    /// Representatives of H⁰ (global sections).
    h0_representatives: Vec<Vec<f64>>,
    /// Representatives of H¹ (obstruction classes).
    h1_representatives: Vec<Vec<f64>>,
}

impl Cohomology {
    /// Compute cohomology from coboundary maps.
    pub fn from_coboundaries(
        d0: Vec<Vec<f64>>,
        d1: Vec<Vec<f64>>,
        c0_dim: usize,
        c1_dim: usize,
    ) -> Self {
        // H⁰ = ker(d0)
        let h0_representatives = Self::kernel(&d0, c0_dim);
        let h0_dim = h0_representatives.len();

        // H¹ = ker(d1) / im(d0)
        let ker_d1 = Self::kernel(&d1, c1_dim);
        let im_d0 = Self::image(&d0, c0_dim);

        // Compute H¹ = ker(d1) / im(d0) by finding vectors in ker(d1) not in im(d0)
        let h1_representatives = Self::quotient_space(&ker_d1, &im_d0);
        let h1_dim = h1_representatives.len();

        // H² = cokernel of d1 (simplified)
        let h2_dim = if d1.is_empty() { 0 } else {
            // Dimension of target - rank(d1)
            let d1_target_dim: usize = d1.first().map(|_| 0).unwrap_or(0);
            let rank_d1 = Self::rank(&d1);
            d1_target_dim.saturating_sub(rank_d1)
        };

        Cohomology {
            h0_dim,
            h1_dim,
            h2_dim,
            d0,
            d1,
            h0_representatives,
            h1_representatives,
        }
    }

    /// Trivial cohomology (H⁰ = 1, everything else = 0).
    pub fn trivial() -> Self {
        Cohomology {
            h0_dim: 1,
            h1_dim: 0,
            h2_dim: 0,
            d0: vec![],
            d1: vec![],
            h0_representatives: vec![vec![1.0]],
            h1_representatives: vec![],
        }
    }

    /// Create cohomology with explicit dimensions (for testing/construction).
    pub fn with_dimensions(h0_dim: usize, h1_dim: usize, h2_dim: usize) -> Self {
        let h0_reps: Vec<Vec<f64>> = (0..h0_dim)
            .map(|i| {
                let mut v = vec![0.0; h0_dim.max(1)];
                if i < v.len() { v[i] = 1.0; }
                v
            })
            .collect();
        let h1_reps: Vec<Vec<f64>> = (0..h1_dim)
            .map(|i| {
                let mut v = vec![0.0; h1_dim.max(1)];
                if i < v.len() { v[i] = 1.0; }
                v
            })
            .collect();
        Cohomology {
            h0_dim,
            h1_dim,
            h2_dim,
            d0: vec![],
            d1: vec![],
            h0_representatives: h0_reps,
            h1_representatives: h1_reps,
        }
    }

    /// Dimension of H⁰ (global sections).
    pub fn h0_dimension(&self) -> usize {
        self.h0_dim
    }

    /// Dimension of H¹ (obstruction classes / deadlocks).
    pub fn h1_dimension(&self) -> usize {
        self.h1_dim
    }

    /// Dimension of H².
    pub fn h2_dimension(&self) -> usize {
        self.h2_dim
    }

    /// Get H⁰ representatives (global sections).
    pub fn global_sections(&self) -> &[Vec<f64>] {
        &self.h0_representatives
    }

    /// Get H¹ representatives (obstruction classes).
    pub fn obstruction_classes(&self) -> &[Vec<f64>] {
        &self.h1_representatives
    }

    /// Check if the protocol is deadlock-free (H¹ = 0).
    pub fn is_deadlock_free(&self) -> bool {
        self.h1_dim == 0
    }

    /// Euler characteristic: χ = dim(H⁰) - dim(H¹) + dim(H²).
    pub fn euler_characteristic(&self) -> isize {
        self.h0_dim as isize - self.h1_dim as isize + self.h2_dim as isize
    }

    /// Compute kernel of a matrix (basis for null space).
    fn kernel(matrix: &[Vec<f64>], input_dim: usize) -> Vec<Vec<f64>> {
        if matrix.is_empty() {
            // Everything is in the kernel
            return (0..input_dim)
                .map(|i| {
                    let mut v = vec![0.0; input_dim];
                    v[i] = 1.0;
                    v
                })
                .collect();
        }

        let n = matrix.len();       // rows
        let m = input_dim;           // columns

        // Augmented matrix for RREF
        let mut aug: Vec<Vec<f64>> = matrix.to_vec();

        let mut pivot_cols: Vec<usize> = vec![];
        let mut row = 0;
        for col in 0..m {
            // Find pivot
            let mut pivot_row = None;
            for r in row..n {
                if aug[r].len() > col && aug[r][col].abs() > 1e-10 {
                    pivot_row = Some(r);
                    break;
                }
            }
            if let Some(pr) = pivot_row {
                aug.swap(row, pr);
                let scale = aug[row][col];
                for j in 0..aug[row].len() {
                    aug[row][j] /= scale;
                }
                for r in 0..n {
                    if r != row && aug[r].len() > col && aug[r][col].abs() > 1e-10 {
                        let factor = aug[r][col];
                        for j in 0..aug[r].len() {
                            aug[r][j] -= factor * aug[row][j];
                        }
                    }
                }
                pivot_cols.push(col);
                row += 1;
            }
        }

        // Free variables are columns not in pivot_cols
        let free_cols: Vec<usize> = (0..m).filter(|c| !pivot_cols.contains(c)).collect();

        free_cols.iter().map(|&fc| {
            let mut v = vec![0.0; m];
            v[fc] = 1.0;
            // Set pivot variables
            for (i, &pc) in pivot_cols.iter().enumerate() {
                if i < aug.len() && aug[i].len() > fc {
                    v[pc] = -aug[i][fc];
                }
            }
            v
        }).collect()
    }

    /// Compute image of a matrix (column space basis).
    fn image(matrix: &[Vec<f64>], input_dim: usize) -> Vec<Vec<f64>> {
        if matrix.is_empty() || input_dim == 0 {
            return vec![];
        }

        let output_dim = matrix.len();
        // Transpose to work with columns
        let n_cols = input_dim;
        let mut basis: Vec<Vec<f64>> = vec![];

        for col in 0..n_cols {
            let mut col_vec: Vec<f64> = vec![0.0; output_dim];
            for row in 0..output_dim {
                if matrix[row].len() > col {
                    col_vec[row] = matrix[row][col];
                }
            }
            // Check if linearly independent from current basis
            if Self::is_independent(&col_vec, &basis) {
                basis.push(col_vec);
            }
        }
        basis
    }

    /// Check if a vector is linearly independent from a set.
    fn is_independent(v: &[f64], basis: &[Vec<f64>]) -> bool {
        if basis.is_empty() {
            return v.iter().any(|&x| x.abs() > 1e-10);
        }
        // Simple Gram-Schmidt check
        let mut w = v.to_vec();
        for b in basis {
            let dot: f64 = w.iter().zip(b.iter()).map(|(a, b)| a * b).sum();
            let norm_sq: f64 = b.iter().map(|x| x * x).sum();
            if norm_sq > 1e-20 {
                for (wi, bi) in w.iter_mut().zip(b.iter()) {
                    *wi -= (dot / norm_sq) * bi;
                }
            }
        }
        w.iter().any(|&x| x.abs() > 1e-10)
    }

    /// Compute rank of a matrix.
    fn rank(matrix: &[Vec<f64>]) -> usize {
        if matrix.is_empty() { return 0; }
        let mut m = matrix.to_vec();
        let rows = m.len();
        let cols = m.first().map(|r| r.len()).unwrap_or(0);
        let mut rank = 0;
        let mut row = 0;
        for col in 0..cols {
            let mut pivot = None;
            for r in row..rows {
                if m[r].len() > col && m[r][col].abs() > 1e-10 {
                    pivot = Some(r);
                    break;
                }
            }
            if let Some(pr) = pivot {
                m.swap(row, pr);
                let scale = m[row][col];
                for j in 0..m[row].len() {
                    m[row][j] /= scale;
                }
                for r in 0..rows {
                    if r != row && m[r].len() > col && m[r][col].abs() > 1e-10 {
                        let factor = m[r][col];
                        for j in 0..m[r].len() {
                            m[r][j] -= factor * m[row][j];
                        }
                    }
                }
                rank += 1;
                row += 1;
            }
        }
        rank
    }

    /// Compute quotient space ker/im (representatives).
    fn quotient_space(ker: &[Vec<f64>], im: &[Vec<f64>]) -> Vec<Vec<f64>> {
        if ker.is_empty() {
            return vec![];
        }
        if im.is_empty() {
            return ker.to_vec();
        }

        // Project each kernel vector away from image span
        let mut result: Vec<Vec<f64>> = vec![];
        for k in ker {
            let mut residual = k.clone();
            // Gram-Schmidt: remove image components
            for i_vec in im {
                let dot: f64 = residual.iter().zip(i_vec.iter()).map(|(a, b)| a * b).sum();
                let norm_sq: f64 = i_vec.iter().map(|x| x * x).sum();
                if norm_sq > 1e-20 {
                    for (r, iv) in residual.iter_mut().zip(i_vec.iter()) {
                        *r -= (dot / norm_sq) * iv;
                    }
                }
            }
            // Also remove previously found quotient components
            for prev in &result {
                let dot: f64 = residual.iter().zip(prev.iter()).map(|(a, b)| a * b).sum();
                let norm_sq: f64 = prev.iter().map(|x| x * x).sum();
                if norm_sq > 1e-20 {
                    for (r, pv) in residual.iter_mut().zip(prev.iter()) {
                        *r -= (dot / norm_sq) * pv;
                    }
                }
            }
            if residual.iter().any(|&x| x.abs() > 1e-10) {
                result.push(residual);
            }
        }
        result
    }

    /// Betti numbers: (h0, h1, h2).
    pub fn betti_numbers(&self) -> (usize, usize, usize) {
        (self.h0_dim, self.h1_dim, self.h2_dim)
    }
}

/// A cochain complex: C⁰ → C¹ → C² with coboundary maps.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CochainComplex {
    pub c0_dim: usize,
    pub c1_dim: usize,
    pub c2_dim: usize,
    pub d0: Vec<Vec<f64>>,
    pub d1: Vec<Vec<f64>>,
}

impl CochainComplex {
    pub fn new(c0_dim: usize, c1_dim: usize, c2_dim: usize) -> Self {
        CochainComplex {
            c0_dim,
            c1_dim,
            c2_dim,
            d0: vec![vec![0.0; c0_dim]; c1_dim],
            d1: vec![vec![0.0; c1_dim]; c2_dim],
        }
    }

    /// Set the d0 coboundary map.
    pub fn set_d0(&mut self, d0: Vec<Vec<f64>>) {
        self.d0 = d0;
    }

    /// Set the d1 coboundary map.
    pub fn set_d1(&mut self, d1: Vec<Vec<f64>>) {
        self.d1 = d1;
    }

    /// Compute cohomology from this cochain complex.
    pub fn compute_cohomology(&self) -> Cohomology {
        Cohomology::from_coboundaries(
            self.d0.clone(),
            self.d1.clone(),
            self.c0_dim,
            self.c1_dim,
        )
    }

    /// Verify the cochain complex property: d1 ∘ d0 = 0.
    pub fn verify_complex(&self) -> bool {
        // d1 * d0 should be zero
        let n = self.d1.len();
        let m = self.d0.first().map(|r| r.len()).unwrap_or(0);
        for i in 0..n {
            for j in 0..m {
                let val: f64 = self.d1[i].iter()
                    .zip(self.d0.iter().map(|row| row.get(j).copied().unwrap_or(0.0)))
                    .map(|(a, b)| a * b)
                    .sum();
                if val.abs() > 1e-10 {
                    return false;
                }
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trivial_cohomology() {
        let c = Cohomology::trivial();
        assert_eq!(c.h0_dimension(), 1);
        assert_eq!(c.h1_dimension(), 0);
        assert!(c.is_deadlock_free());
    }

    #[test]
    fn test_with_dimensions() {
        let c = Cohomology::with_dimensions(2, 1, 0);
        assert_eq!(c.h0_dimension(), 2);
        assert_eq!(c.h1_dimension(), 1);
        assert!(!c.is_deadlock_free());
    }

    #[test]
    fn test_euler_characteristic() {
        let c = Cohomology::with_dimensions(3, 1, 0);
        assert_eq!(c.euler_characteristic(), 2);
    }

    #[test]
    fn test_betti_numbers() {
        let c = Cohomology::with_dimensions(2, 3, 1);
        assert_eq!(c.betti_numbers(), (2, 3, 1));
    }

    #[test]
    fn test_kernel_identity() {
        let matrix = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let ker = Cohomology::kernel(&matrix, 2);
        assert_eq!(ker.len(), 0); // Only trivial kernel
    }

    #[test]
    fn test_kernel_zero_matrix() {
        let matrix = vec![vec![0.0, 0.0], vec![0.0, 0.0]];
        let ker = Cohomology::kernel(&matrix, 2);
        assert_eq!(ker.len(), 2); // Full kernel
    }

    #[test]
    fn test_kernel_rank_1() {
        let matrix = vec![vec![1.0, 1.0]];
        let ker = Cohomology::kernel(&matrix, 2);
        assert_eq!(ker.len(), 1);
        // Should be proportional to (1, -1)
        assert!((ker[0][0] - 1.0).abs() < 1e-10 || (ker[0][0] + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_image_full_rank() {
        let matrix = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let im = Cohomology::image(&matrix, 2);
        assert_eq!(im.len(), 2);
    }

    #[test]
    fn test_image_rank_1() {
        let matrix = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        let im = Cohomology::image(&matrix, 2);
        assert_eq!(im.len(), 1);
    }

    #[test]
    fn test_cochain_complex() {
        let mut cc = CochainComplex::new(3, 2, 1);
        assert!(cc.verify_complex()); // Zero maps compose to zero
    }

    #[test]
    fn test_cochain_compute_cohomology() {
        let mut cc = CochainComplex::new(2, 2, 0);
        cc.set_d0(vec![vec![1.0, 0.0], vec![0.0, 1.0]]);
        let c = cc.compute_cohomology();
        assert_eq!(c.h0_dimension(), 0); // d0 is injective
    }

    #[test]
    fn test_rank() {
        let m = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        assert_eq!(Cohomology::rank(&m), 2);

        let m2 = vec![vec![1.0, 2.0], vec![3.0, 6.0]];
        assert_eq!(Cohomology::rank(&m2), 1);
    }

    #[test]
    fn test_obstruction_classes_empty() {
        let c = Cohomology::trivial();
        assert!(c.obstruction_classes().is_empty());
    }

    #[test]
    fn test_obstruction_classes_nontrivial() {
        let c = Cohomology::with_dimensions(1, 2, 0);
        assert_eq!(c.obstruction_classes().len(), 2);
    }
}
