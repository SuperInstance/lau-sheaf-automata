//! Sheaf Laplacian for protocol verification.
//!
//! The sheaf Laplacian Δ_L has kernel = H⁰ (global sections).
//! This provides an efficient way to verify protocol consistency.

use serde::{Deserialize, Serialize};

/// Sheaf Laplacian Δ_L for computing global sections.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SheafLaplacian {
    /// The Laplacian matrix.
    pub matrix: Vec<Vec<f64>>,
    /// Dimension of the matrix.
    pub dimension: usize,
    /// Stalk dimensions at each vertex.
    pub stalk_dims: Vec<usize>,
}

impl SheafLaplacian {
    /// Build the sheaf Laplacian from restriction maps and adjacency.
    ///
    /// For vertices i, j with edge (i,j):
    ///   Δ_L = D - A_sheaf
    /// where D is the degree matrix and A_sheaf = ρ^T ρ for each edge.
    pub fn from_sheaf_data(
        _n_vertices: usize,
        stalk_dims: Vec<usize>,
        _edges: &[(usize, usize)],
        restriction_norms: &[(usize, usize, f64)],
    ) -> Self {
        let total_dim: usize = stalk_dims.iter().sum();
        let mut matrix = vec![vec![0.0; total_dim]; total_dim];

        // Build diagonal degree blocks
        for &(i, j, norm) in restriction_norms {
            let dim_i = stalk_dims.get(i).copied().unwrap_or(0);
            let dim_j = stalk_dims.get(j).copied().unwrap_or(0);

            let offset_i: usize = stalk_dims.iter().take(i).sum();
            let offset_j: usize = stalk_dims.iter().take(j).sum();

            // D[i,i] += norm² on diagonal blocks
            for k in 0..dim_i {
                matrix[offset_i + k][offset_i + k] += norm * norm;
            }
            // D[j,j] += norm²
            for k in 0..dim_j {
                matrix[offset_j + k][offset_j + k] += norm * norm;
            }

            // A_sheaf[i,j] = -ρ^T ρ
            let min_dim = dim_i.min(dim_j);
            for k in 0..min_dim {
                matrix[offset_i + k][offset_j + k] -= norm * norm;
                matrix[offset_j + k][offset_i + k] -= norm * norm;
            }
        }

        SheafLaplacian { matrix, dimension: total_dim, stalk_dims }
    }

    /// Identity Laplacian for trivial sheaf on n vertices.
    pub fn identity(n_vertices: usize, stalk_dim: usize) -> Self {
        let total_dim = n_vertices * stalk_dim;
        let stalk_dims = vec![stalk_dim; n_vertices];
        let matrix = (0..total_dim)
            .map(|i| (0..total_dim).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        SheafLaplacian { matrix, dimension: total_dim, stalk_dims }
    }

    /// Zero Laplacian.
    pub fn zero(n_vertices: usize, stalk_dim: usize) -> Self {
        let total_dim = n_vertices * stalk_dim;
        let stalk_dims = vec![stalk_dim; n_vertices];
        SheafLaplacian {
            matrix: vec![vec![0.0; total_dim]; total_dim],
            dimension: total_dim,
            stalk_dims,
        }
    }

    /// Compute the kernel of the Laplacian (= H⁰ = global sections).
    pub fn kernel(&self) -> Vec<Vec<f64>> {
        // RREF on the Laplacian matrix
        let n = self.dimension;
        if n == 0 { return vec![]; }

        let mut aug = self.matrix.clone();
        let mut pivot_cols: Vec<usize> = vec![];
        let mut row = 0;

        for col in 0..n {
            let mut pivot_row = None;
            for r in row..n {
                if aug[r][col].abs() > 1e-10 {
                    pivot_row = Some(r);
                    break;
                }
            }
            if let Some(pr) = pivot_row {
                aug.swap(row, pr);
                let scale = aug[row][col];
                for j in 0..n {
                    aug[row][j] /= scale;
                }
                for r in 0..n {
                    if r != row && aug[r][col].abs() > 1e-10 {
                        let factor = aug[r][col];
                        for j in 0..n {
                            aug[r][j] -= factor * aug[row][j];
                        }
                    }
                }
                pivot_cols.push(col);
                row += 1;
            }
        }

        let free_cols: Vec<usize> = (0..n).filter(|c| !pivot_cols.contains(c)).collect();
        free_cols.iter().map(|&fc| {
            let mut v = vec![0.0; n];
            v[fc] = 1.0;
            for (i, &pc) in pivot_cols.iter().enumerate() {
                if i < aug.len() {
                    v[pc] = -aug[i][fc];
                }
            }
            v
        }).collect()
    }

    /// Dimension of kernel (= dim H⁰).
    pub fn kernel_dimension(&self) -> usize {
        self.kernel().len()
    }

    /// Apply the Laplacian to a vector.
    pub fn apply(&self, v: &[f64]) -> Vec<f64> {
        if self.matrix.is_empty() { return vec![]; }
        self.matrix.iter()
            .map(|row| {
                row.iter().zip(v.iter())
                    .map(|(m, vi)| m * vi)
                    .sum()
            })
            .collect()
    }

    /// Check if a section is in the kernel (global section).
    pub fn is_global_section(&self, v: &[f64]) -> bool {
        let result = self.apply(v);
        result.iter().all(|&x| x.abs() < 1e-10)
    }

    /// Compute the eigenvalues (approximate, sorted ascending).
    pub fn eigenvalues(&self) -> Vec<f64> {
        let n = self.dimension;
        if n == 0 { return vec![]; }
        if n == 1 { return vec![self.matrix[0][0]]; }

        // QR iteration: repeatedly decompose M = QR, then set M = RQ
        let mut m = self.matrix.clone();
        for _ in 0..200 {
            let (q, r) = Self::qr_decompose(&m);
            // M_new = R * Q
            m = Self::mat_mul(&r, &q);
        }

        let mut eigs: Vec<f64> = (0..n).map(|i| m[i][i]).collect();
        eigs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        eigs
    }

    /// QR decomposition using Householder reflections.
    fn qr_decompose(m: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
        let n = m.len();
        let mut q = Self::identity_matrix(n);
        let mut r = m.to_vec();

        for col in 0..n {
            // Build Householder reflector for column `col`
            let mut v = vec![0.0; n];
            let norm: f64 = (col..n).map(|i| r[i][col] * r[i][col]).sum::<f64>().sqrt();
            if norm < 1e-15 { continue; }
            let sign = if r[col][col] >= 0.0 { 1.0 } else { -1.0 };
            v[col] = r[col][col] + sign * norm;
            for i in (col + 1)..n {
                v[i] = r[i][col];
            }
            let v_norm_sq: f64 = v.iter().map(|x| x * x).sum();
            if v_norm_sq < 1e-30 { continue; }

            // Apply H = I - 2vv^T/v^Tv to r and q
            // r = H * r
            for j in col..n {
                let dot: f64 = v.iter().zip(r.iter().map(|row| row.get(j).copied().unwrap_or(0.0))).map(|(a, b)| a * b).sum();
                for i in col..n {
                    r[i][j] -= 2.0 * dot * v[i] / v_norm_sq;
                }
            }
            // q = q * H
            for i in 0..n {
                let dot: f64 = q[i].iter().zip(v.iter()).map(|(a, b)| a * b).sum();
                for j in col..n {
                    q[i][j] -= 2.0 * dot * v[j] / v_norm_sq;
                }
            }
        }

        (q, r)
    }

    fn identity_matrix(n: usize) -> Vec<Vec<f64>> {
        (0..n).map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect()).collect()
    }

    /// Spectral gap (smallest non-zero eigenvalue).
    pub fn spectral_gap(&self) -> f64 {
        let eigs = self.eigenvalues();
        eigs.iter()
            .find(|&&e| e > 1e-10)
            .copied()
            .unwrap_or(0.0)
    }

    fn transpose(m: &[Vec<f64>]) -> Vec<Vec<f64>> {
        if m.is_empty() { return vec![]; }
        let rows = m.len();
        let cols = m[0].len();
        (0..cols)
            .map(|j| (0..rows).map(|i| m[i][j]).collect())
            .collect()
    }

    fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
        if a.is_empty() || b.is_empty() { return vec![]; }
        let n = a.len();
        let m = b[0].len();
        let k = b.len();
        (0..n)
            .map(|i| (0..m)
                .map(|j| (0..k).map(|l| a[i][l] * b[l][j]).sum())
                .collect())
            .collect()
    }

}

/// Verify a protocol using the sheaf Laplacian.
pub fn verify_protocol(laplacian: &SheafLaplacian) -> ProtocolVerification {
    let kernel = laplacian.kernel();
    let h0_dim = kernel.len();
    let _eigs = laplacian.eigenvalues();
    let spectral_gap = laplacian.spectral_gap();

    ProtocolVerification {
        h0_dimension: h0_dim,
        global_sections: kernel,
        spectral_gap,
        is_consistent: h0_dim > 0,
    }
}

/// Result of protocol verification.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProtocolVerification {
    pub h0_dimension: usize,
    pub global_sections: Vec<Vec<f64>>,
    pub spectral_gap: f64,
    pub is_consistent: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity_laplacian() {
        let lap = SheafLaplacian::identity(3, 2);
        assert_eq!(lap.dimension, 6);
        assert_eq!(lap.kernel_dimension(), 0); // Full rank
    }

    #[test]
    fn test_zero_laplacian() {
        let lap = SheafLaplacian::zero(2, 2);
        assert_eq!(lap.kernel_dimension(), 4); // Everything in kernel
    }

    #[test]
    fn test_laplacian_from_sheaf_data() {
        let lap = SheafLaplacian::from_sheaf_data(
            2,
            vec![1, 1],
            &[(0, 1)],
            &[(0, 1, 1.0)],
        );
        assert_eq!(lap.dimension, 2);
        // Graph Laplacian on 2 vertices: [[1, -1], [-1, 1]]
        assert_eq!(lap.matrix[0][0], 1.0);
        assert_eq!(lap.matrix[0][1], -1.0);
        assert_eq!(lap.matrix[1][0], -1.0);
        assert_eq!(lap.matrix[1][1], 1.0);
    }

    #[test]
    fn test_laplacian_kernel_graph() {
        let lap = SheafLaplacian::from_sheaf_data(
            2,
            vec![1, 1],
            &[(0, 1)],
            &[(0, 1, 1.0)],
        );
        let ker = lap.kernel();
        assert_eq!(ker.len(), 1); // Connected graph has 1-dimensional kernel (constants)
        // Should be proportional to (1, 1)
        assert!((ker[0][0] - ker[0][1]).abs() < 1e-10);
    }

    #[test]
    fn test_is_global_section() {
        let lap = SheafLaplacian::from_sheaf_data(
            2,
            vec![1, 1],
            &[(0, 1)],
            &[(0, 1, 1.0)],
        );
        assert!(lap.is_global_section(&[1.0, 1.0]));
        assert!(!lap.is_global_section(&[1.0, 0.0]));
    }

    #[test]
    fn test_apply_laplacian() {
        let lap = SheafLaplacian::identity(2, 1);
        let result = lap.apply(&[3.0, 5.0]);
        assert_eq!(result, vec![3.0, 5.0]);
    }

    #[test]
    fn test_laplacian_triangle() {
        let lap = SheafLaplacian::from_sheaf_data(
            3,
            vec![1, 1, 1],
            &[(0, 1), (1, 2), (0, 2)],
            &[(0, 1, 1.0), (1, 2, 1.0), (0, 2, 1.0)],
        );
        let ker = lap.kernel();
        assert_eq!(ker.len(), 1); // One global section
    }

    #[test]
    fn test_laplacian_disconnected() {
        let lap = SheafLaplacian::from_sheaf_data(
            4,
            vec![1, 1, 1, 1],
            &[(0, 1), (2, 3)], // Two disconnected components
            &[(0, 1, 1.0), (2, 3, 1.0)],
        );
        let ker = lap.kernel();
        assert_eq!(ker.len(), 2); // Two global sections
    }

    #[test]
    fn test_spectral_gap() {
        let lap = SheafLaplacian::from_sheaf_data(
            2,
            vec![1, 1],
            &[(0, 1)],
            &[(0, 1, 1.0)],
        );
        let gap = lap.spectral_gap();
        assert!(gap > 0.0);
    }

    #[test]
    fn test_laplacian_higher_dim_stalks() {
        let lap = SheafLaplacian::from_sheaf_data(
            2,
            vec![2, 2],
            &[(0, 1)],
            &[(0, 1, 1.0)],
        );
        assert_eq!(lap.dimension, 4);
        let ker = lap.kernel();
        assert!(ker.len() >= 1);
    }
}
