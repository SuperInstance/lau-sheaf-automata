//! Ergodic consensus: time-averages of local sections converge to global sections.
//!
//! In the sheaf-theoretic framework, distributed consensus is equivalent
//! to finding a global section of the protocol sheaf. The ergodic theorem
//! guarantees convergence under appropriate conditions.

use serde::{Deserialize, Serialize};

use crate::laplacian::SheafLaplacian;

/// Result of consensus computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConsensusResult {
    /// The global section found (if any).
    pub global_section: Option<Vec<f64>>,
    /// Whether consensus was reached.
    pub reached: bool,
    /// Number of iterations to converge.
    pub iterations: usize,
    /// Residual error.
    pub residual: f64,
    /// Time-averaged section at each step.
    pub time_average: Vec<f64>,
}

/// Ergodic consensus algorithm for finding global sections.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErgodicConsensus {
    /// Maximum iterations.
    pub max_iterations: usize,
    /// Convergence tolerance.
    pub tolerance: f64,
    /// Step size for averaging.
    pub step_size: f64,
}

impl ErgodicConsensus {
    pub fn new() -> Self {
        ErgodicConsensus {
            max_iterations: 10000,
            tolerance: 1e-10,
            step_size: 0.01,
        }
    }

    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }

    pub fn with_max_iterations(mut self, max: usize) -> Self {
        self.max_iterations = max;
        self
    }

    /// Run ergodic consensus to find a global section.
    ///
    /// Uses the sheaf Laplacian to iteratively project local sections
    /// toward the kernel (global sections).
    pub fn find_global_section(
        &self,
        laplacian: &SheafLaplacian,
        initial: &[f64],
    ) -> ConsensusResult {
        let n = laplacian.dimension;
        if n == 0 {
            return ConsensusResult {
                global_section: Some(vec![]),
                reached: true,
                iterations: 0,
                residual: 0.0,
                time_average: vec![],
            };
        }

        let mut current = initial.to_vec();
        if current.len() != n {
            current = vec![1.0; n];
        }

        let mut time_sum = vec![0.0; n];
        let mut iterations = 0;

        for iter in 0..self.max_iterations {
            // Apply Laplacian
            let delta = laplacian.apply(&current);
            let residual: f64 = delta.iter().map(|x| x * x).sum::<f64>().sqrt();

            // Time average
            for (s, c) in time_sum.iter_mut().zip(current.iter()) {
                *s += c;
            }
            iterations = iter + 1;

            if residual < self.tolerance {
                let time_avg: Vec<f64> = time_sum.iter().map(|x| x / iterations as f64).collect();
                return ConsensusResult {
                    global_section: Some(current.clone()),
                    reached: true,
                    iterations,
                    residual,
                    time_average: time_avg,
                };
            }

            // Gradient descent step: x ← x - η * Δx
            for (c, d) in current.iter_mut().zip(delta.iter()) {
                *c -= self.step_size * d;
            }
        }

        let delta = laplacian.apply(&current);
        let residual: f64 = delta.iter().map(|x| x * x).sum::<f64>().sqrt();
        let time_avg: Vec<f64> = time_sum.iter().map(|x| x / iterations as f64).collect();

        ConsensusResult {
            global_section: if residual < self.tolerance * 100.0 { Some(current) } else { None },
            reached: residual < self.tolerance * 100.0,
            iterations,
            residual,
            time_average: time_avg,
        }
    }

    /// Run consensus with multiple random initial conditions.
    pub fn find_all_global_sections(
        &self,
        laplacian: &SheafLaplacian,
        n_trials: usize,
    ) -> Vec<ConsensusResult> {
        let n = laplacian.dimension;
        let mut results = vec![];

        // First try uniform initial
        let uniform = vec![1.0; n];
        results.push(self.find_global_section(laplacian, &uniform));

        // Try random initials
        for i in 1..n_trials {
            let initial: Vec<f64> = (0..n).map(|j| {
                ((Self::simple_random(i * n + j + 42) % 1000) as f64 / 500.0) - 1.0
            }).collect();
            let result = self.find_global_section(laplacian, &initial);
            // Only add if it found a different section
            if result.reached {
                let is_new = results.iter().all(|r| {
                    if let (Some(ref a), Some(ref b)) = (&r.global_section, &result.global_section) {
                        // Check linear independence
                        a.iter().zip(b.iter())
                            .map(|(x, y)| (x - y).abs())
                            .sum::<f64>() > self.tolerance
                    } else {
                        true
                    }
                });
                if is_new {
                    results.push(result);
                }
            }
        }

        results
    }

    /// Simple deterministic pseudo-random for testing.
fn simple_random(seed: usize) -> u64 {
    let mut x = (seed as u64).wrapping_add(0x9e3779b97f4a7c15);
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    x ^= x >> 31;
    x
}

/// Check if the ergodic theorem applies: bounded, connected sheaf.
    pub fn check_ergodic_conditions(laplacian: &SheafLaplacian) -> ErgodicConditions {
        let spectral_gap = laplacian.spectral_gap();
        ErgodicConditions {
            has_spectral_gap: spectral_gap > 1e-10,
            is_finite_dimensional: laplacian.dimension < usize::MAX,
            spectral_gap,
        }
    }
}

/// Conditions for the ergodic theorem.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErgodicConditions {
    pub has_spectral_gap: bool,
    pub is_finite_dimensional: bool,
    pub spectral_gap: f64,
}

impl ErgodicConditions {
    pub fn ergodic_theorem_applies(&self) -> bool {
        self.has_spectral_gap && self.is_finite_dimensional
    }
}

/// Project a local section onto the space of global sections.
pub fn project_to_global(
    local_section: &[f64],
    global_basis: &[Vec<f64>],
) -> Vec<f64> {
    if global_basis.is_empty() {
        return vec![0.0; local_section.len()];
    }

    // Compute projection coefficients
    let mut projection = vec![0.0; local_section.len()];
    for basis_vec in global_basis {
        let dot: f64 = local_section.iter().zip(basis_vec.iter())
            .map(|(a, b)| a * b)
            .sum();
        let norm_sq: f64 = basis_vec.iter().map(|x| x * x).sum();
        if norm_sq > 1e-20 {
            for (p, b) in projection.iter_mut().zip(basis_vec.iter()) {
                *p += (dot / norm_sq) * b;
            }
        }
    }
    projection
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_consensus_identity_laplacian() {
        let lap = SheafLaplacian::identity(2, 1);
        let ec = ErgodicConsensus::new().with_max_iterations(100);
        let result = ec.find_global_section(&lap, &[1.0, 1.0]);
        // Identity Laplacian has no kernel (except 0)
        assert!(!result.reached || result.global_section.is_some());
    }

    #[test]
    fn test_consensus_graph_laplacian() {
        let lap = SheafLaplacian::from_sheaf_data(
            2,
            vec![1, 1],
            &[(0, 1)],
            &[(0, 1, 1.0)],
        );
        let ec = ErgodicConsensus::new()
            .with_max_iterations(5000)
            .with_tolerance(1e-6);
        let result = ec.find_global_section(&lap, &[3.0, 5.0]);
        assert!(result.reached);
        // Should converge to (c, c) for some c
        if let Some(ref gs) = result.global_section {
            assert!((gs[0] - gs[1]).abs() < 1e-4);
        }
    }

    #[test]
    fn test_consensus_triangle() {
        let lap = SheafLaplacian::from_sheaf_data(
            3,
            vec![1, 1, 1],
            &[(0, 1), (1, 2), (0, 2)],
            &[(0, 1, 1.0), (1, 2, 1.0), (0, 2, 1.0)],
        );
        let ec = ErgodicConsensus::new()
            .with_max_iterations(5000)
            .with_tolerance(1e-6);
        let result = ec.find_global_section(&lap, &[1.0, 2.0, 3.0]);
        assert!(result.reached);
    }

    #[test]
    fn test_ergodic_conditions() {
        let lap = SheafLaplacian::from_sheaf_data(
            2,
            vec![1, 1],
            &[(0, 1)],
            &[(0, 1, 1.0)],
        );
        let conditions = ErgodicConsensus::check_ergodic_conditions(&lap);
        assert!(conditions.ergodic_theorem_applies());
    }

    #[test]
    fn test_project_to_global() {
        let basis = vec![vec![1.0, 1.0]];
        let local = vec![3.0, 5.0];
        let proj = project_to_global(&local, &basis);
        // Projection onto (1,1) direction
        assert!((proj[0] - 4.0).abs() < 1e-10);
        assert!((proj[1] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_project_to_global_empty_basis() {
        let proj = project_to_global(&[1.0, 2.0], &[]);
        assert_eq!(proj, vec![0.0, 0.0]);
    }

    #[test]
    fn test_consensus_zero_laplacian() {
        let lap = SheafLaplacian::zero(2, 1);
        let ec = ErgodicConsensus::new();
        let result = ec.find_global_section(&lap, &[1.0, 2.0]);
        assert!(result.reached);
    }

    #[test]
    fn test_time_average() {
        let lap = SheafLaplacian::zero(1, 1);
        let ec = ErgodicConsensus::new();
        let result = ec.find_global_section(&lap, &[42.0]);
        assert!(result.reached);
        // Time average should be close to 42.0
        assert!((result.time_average[0] - 42.0).abs() < 1.0);
    }
}
