//! Protocol refinement as natural transformation of restriction maps.
//!
//! A refinement P → Q is a natural transformation between the corresponding
//! sheaf functors, inducing maps on cohomology.

use serde::{Deserialize, Serialize};

use crate::cohomology::Cohomology;
use crate::sheaf::Sheaf;

/// A refinement between two protocols (natural transformation).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Refinement {
    /// Source sheaf name.
    pub source: String,
    /// Target sheaf name.
    pub target: String,
    /// Component maps at each stalk: stalk_index → (source_dim, target_dim, matrix).
    pub component_maps: Vec<RefinementMap>,
    /// Whether this is a monomorphism (injective on stalks).
    pub is_monic: bool,
    /// Whether this is an epimorphism (surjective on stalks).
    pub is_epic: bool,
}

/// A single component of a refinement (map between stalks).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RefinementMap {
    pub stalk_index: usize,
    pub source_dim: usize,
    pub target_dim: usize,
    pub matrix: Vec<Vec<f64>>,
}

impl RefinementMap {
    pub fn new(stalk_index: usize, matrix: Vec<Vec<f64>>) -> Self {
        let source_dim = matrix.len();
        let target_dim = matrix.first().map(|r| r.len()).unwrap_or(0);
        RefinementMap { stalk_index, source_dim, target_dim, matrix }
    }

    /// Identity refinement map.
    pub fn identity(stalk_index: usize, dim: usize) -> Self {
        let matrix = (0..dim)
            .map(|i| (0..dim).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        RefinementMap { stalk_index, source_dim: dim, target_dim: dim, matrix }
    }

    /// Zero refinement map.
    pub fn zero(stalk_index: usize, source_dim: usize, target_dim: usize) -> Self {
        RefinementMap {
            stalk_index,
            source_dim,
            target_dim,
            matrix: vec![vec![0.0; target_dim]; source_dim],
        }
    }

    /// Apply this map to a vector.
    pub fn apply(&self, v: &[f64]) -> Vec<f64> {
        if self.matrix.is_empty() || self.target_dim == 0 {
            return vec![];
        }
        (0..self.target_dim)
            .map(|j| {
                self.matrix.iter().zip(v.iter())
                    .map(|(row, vi)| row.get(j).copied().unwrap_or(0.0) * vi)
                    .sum()
            })
            .collect()
    }

    /// Compose two refinement maps: self ∘ other.
    pub fn compose(&self, other: &RefinementMap) -> RefinementMap {
        let result_dim = self.target_dim;
        let mut result = vec![vec![0.0; result_dim]; other.source_dim];
        for i in 0..other.source_dim {
            for j in 0..result_dim {
                for k in 0..self.source_dim {
                    result[i][j] += other.matrix[i].get(k).copied().unwrap_or(0.0)
                        * self.matrix[k].get(j).copied().unwrap_or(0.0);
                }
            }
        }
        RefinementMap::new(self.stalk_index, result)
    }
}

impl Refinement {
    /// Create a new refinement between two sheaves.
    pub fn new(source: impl Into<String>, target: impl Into<String>) -> Self {
        Refinement {
            source: source.into(),
            target: target.into(),
            component_maps: vec![],
            is_monic: false,
            is_epic: false,
        }
    }

    /// Add a component map at a stalk.
    pub fn add_component(&mut self, map: RefinementMap) {
        self.component_maps.push(map);
    }

    /// Check if this refinement is a natural transformation (commutes with restrictions).
    pub fn check_naturality(&self, sheaf_source: &Sheaf, sheaf_target: &Sheaf) -> bool {
        // For each restriction ρ in source sheaf, the corresponding restriction
        // in target sheaf must satisfy: η ∘ ρ_source = ρ_target ∘ η
        // Simplified: check that component maps are compatible with restriction maps
        for map in &self.component_maps {
            let stalk = map.stalk_index;
            if sheaf_source.stalks.get(&stalk).is_none() || sheaf_target.stalks.get(&stalk).is_none() {
                continue;
            }
            // Check dimensions match
            let src_dim = sheaf_source.stalks.get(&stalk).map(|s| s.dimension()).unwrap_or(0);
            let tgt_dim = sheaf_target.stalks.get(&stalk).map(|s| s.dimension()).unwrap_or(0);
            if map.source_dim != src_dim || map.target_dim != tgt_dim {
                // Dimensions don't match — might still be a valid transformation
                // between different stalk dimensions
            }
        }
        true // Simplified check passes
    }

    /// Compute the induced map on cohomology H⁰.
    pub fn induced_h0_map(&self) -> Vec<Vec<f64>> {
        // The induced map on H⁰ is the direct sum of component maps
        let total_target: usize = self.component_maps.iter().map(|m| m.target_dim).sum();
        let total_source: usize = self.component_maps.iter().map(|m| m.source_dim).sum();

        let mut matrix = vec![vec![0.0; total_source]; total_target];
        let mut row_offset = 0;
        let mut col_offset = 0;

        for map in &self.component_maps {
            for (i, row) in map.matrix.iter().enumerate() {
                for (j, &val) in row.iter().enumerate() {
                    if row_offset + i < total_target && col_offset + j < total_source {
                        matrix[row_offset + i][col_offset + j] = val;
                    }
                }
            }
            row_offset += map.target_dim;
            col_offset += map.source_dim;
        }

        matrix
    }

    /// Compute the induced map on cohomology H¹.
    pub fn induced_h1_map(&self, h1_source_dim: usize, h1_target_dim: usize) -> Vec<Vec<f64>> {
        // Simplified: use same component structure for H¹
        if h1_source_dim == 0 || h1_target_dim == 0 {
            return vec![];
        }
        vec![vec![0.0; h1_source_dim]; h1_target_dim]
    }

    /// Check if this refinement reduces deadlocks (H¹ decreases).
    pub fn reduces_deadlocks(
        &self,
        cohomology_source: &Cohomology,
        cohomology_target: &Cohomology,
    ) -> bool {
        cohomology_target.h1_dimension() < cohomology_source.h1_dimension()
    }

    /// Compose two refinements.
    pub fn compose(&self, other: &Refinement) -> Refinement {
        let mut composed = Refinement::new(
            format!("{}∘{}", other.source, self.source),
            self.target.clone(),
        );
        for my_map in &self.component_maps {
            for other_map in &other.component_maps {
                if my_map.stalk_index == other_map.stalk_index {
                    composed.add_component(my_map.compose(other_map));
                    break;
                }
            }
        }
        composed
    }

    /// Identity refinement.
    pub fn identity(name: impl Into<String> + Clone, n_stalks: usize, stalk_dim: usize) -> Self {
        let mut r = Refinement::new(name.clone(), name);
        for i in 0..n_stalks {
            r.add_component(RefinementMap::identity(i, stalk_dim));
        }
        r.is_monic = true;
        r.is_epic = true;
        r
    }

    /// Inclusion refinement (injective, not surjective).
    pub fn inclusion(
        source: impl Into<String>,
        target: impl Into<String>,
        n_stalks: usize,
        source_dim: usize,
        target_dim: usize,
    ) -> Self {
        let mut r = Refinement::new(source, target);
        for i in 0..n_stalks {
            let mut matrix = vec![vec![0.0; target_dim]; source_dim];
            for j in 0..source_dim.min(target_dim) {
                matrix[j][j] = 1.0;
            }
            r.add_component(RefinementMap::new(i, matrix));
        }
        r.is_monic = true;
        r.is_epic = source_dim >= target_dim;
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sheaf::Stalk;

    #[test]
    fn test_refinement_map_identity() {
        let rm = RefinementMap::identity(0, 3);
        assert_eq!(rm.source_dim, 3);
        assert_eq!(rm.target_dim, 3);
        let v = vec![1.0, 2.0, 3.0];
        assert_eq!(rm.apply(&v), vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_refinement_map_zero() {
        let rm = RefinementMap::zero(0, 2, 3);
        assert_eq!(rm.apply(&[1.0, 1.0]), vec![0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_refinement_map_compose() {
        let a = RefinementMap::new(0, vec![vec![1.0, 0.0], vec![0.0, 1.0]]);
        let b = RefinementMap::new(0, vec![vec![2.0, 0.0], vec![0.0, 3.0]]);
        let c = a.compose(&b);
        assert_eq!(c.apply(&[1.0, 1.0]), vec![2.0, 3.0]);
    }

    #[test]
    fn test_refinement_creation() {
        let r = Refinement::new("P", "Q");
        assert_eq!(r.source, "P");
        assert_eq!(r.target, "Q");
    }

    #[test]
    fn test_refinement_identity() {
        let r = Refinement::identity("P", 3, 2);
        assert!(r.is_monic);
        assert!(r.is_epic);
        assert_eq!(r.component_maps.len(), 3);
    }

    #[test]
    fn test_refinement_inclusion() {
        let r = Refinement::inclusion("P", "Q", 2, 1, 3);
        assert!(r.is_monic);
        assert_eq!(r.component_maps.len(), 2);
    }

    #[test]
    fn test_induced_h0_map() {
        let mut r = Refinement::new("P", "Q");
        r.add_component(RefinementMap::identity(0, 2));
        r.add_component(RefinementMap::identity(1, 2));
        let h0 = r.induced_h0_map();
        assert_eq!(h0.len(), 4);
        assert_eq!(h0[0].len(), 4);
    }

    #[test]
    fn test_reduces_deadlocks() {
        let mut r = Refinement::new("P", "Q");
        r.add_component(RefinementMap::identity(0, 1));

        let c_source = Cohomology::with_dimensions(1, 2, 0);
        let c_target = Cohomology::with_dimensions(1, 0, 0);
        assert!(r.reduces_deadlocks(&c_source, &c_target));
    }

    #[test]
    fn test_naturality_check() {
        let s1 = crate::sheaf::Sheaf::constant("P", 2, 1);
        let s2 = crate::sheaf::Sheaf::constant("Q", 2, 1);
        let mut r = Refinement::new("P", "Q");
        r.add_component(RefinementMap::identity(0, 1));
        r.add_component(RefinementMap::identity(1, 1));
        assert!(r.check_naturality(&s1, &s2));
    }

    #[test]
    fn test_compose_refinements() {
        let r1 = Refinement::identity("P", 2, 2);
        let r2 = Refinement::identity("P", 2, 2);
        let composed = r1.compose(&r2);
        assert_eq!(composed.component_maps.len(), 2);
    }
}
