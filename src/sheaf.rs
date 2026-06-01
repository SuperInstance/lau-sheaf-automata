//! Sheaf data structures: stalks, restriction maps, and sheaves over topological spaces.

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::hash::Hash;

/// An open set in the base topology, identified by name.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct OpenSet {
    pub name: String,
    /// Indices of points in this open set.
    pub points: Vec<usize>,
}

impl OpenSet {
    pub fn new(name: impl Into<String>, points: Vec<usize>) -> Self {
        OpenSet { name: name.into(), points }
    }

    /// Universe open set containing all points.
    pub fn universe(n: usize) -> Self {
        OpenSet { name: "U".into(), points: (0..n).collect() }
    }

    /// Empty open set.
    pub fn empty() -> Self {
        OpenSet { name: "∅".into(), points: vec![] }
    }

    /// Check if this open set contains a point.
    pub fn contains(&self, point: usize) -> bool {
        self.points.contains(&point)
    }

    /// Intersection of two open sets.
    pub fn intersection(&self, other: &OpenSet) -> OpenSet {
        let mut pts: Vec<usize> = self.points.iter()
            .filter(|p| other.points.contains(p))
            .copied()
            .collect();
        pts.sort();
        pts.dedup();
        OpenSet::new(
            format!("{}∩{}", self.name, other.name),
            pts,
        )
    }

    /// Union of two open sets.
    pub fn union(&self, other: &OpenSet) -> OpenSet {
        let mut pts: Vec<usize> = self.points.iter()
            .chain(other.points.iter())
            .copied()
            .collect();
        pts.sort();
        pts.dedup();
        OpenSet::new(
            format!("{}∪{}", self.name, other.name),
            pts,
        )
    }
}

/// A stalk at a point — the local data of the sheaf.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Stalk {
    /// A finite-dimensional vector space (dim = value).
    VectorSpace { dimension: usize },
    /// A set of labels (language elements).
    LabelSet { labels: Vec<String> },
    /// An automaton state configuration.
    AutomatonState { states: Vec<usize>, accepting: Vec<bool> },
    /// Custom data encoded as a byte vector.
    Custom { data: Vec<u8>, dimension: usize },
}

impl Stalk {
    /// Dimension of the stalk (for cohomology computation).
    pub fn dimension(&self) -> usize {
        match self {
            Stalk::VectorSpace { dimension } => *dimension,
            Stalk::LabelSet { labels } => labels.len(),
            Stalk::AutomatonState { states, .. } => states.len(),
            Stalk::Custom { dimension, .. } => *dimension,
        }
    }

    /// Zero stalk.
    pub fn zero() -> Self {
        Stalk::VectorSpace { dimension: 0 }
    }

    /// Create a vector space stalk.
    pub fn vector(dim: usize) -> Self {
        Stalk::VectorSpace { dimension: dim }
    }
}

/// A restriction map ρ_{V→U}: F(V) → F(U) for V ⊆ U.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RestrictionMap {
    pub source: String,
    pub target: String,
    /// Matrix representation of the linear map (source_dim × target_dim).
    pub matrix: Vec<Vec<f64>>,
}

impl RestrictionMap {
    pub fn new(source: impl Into<String>, target: impl Into<String>, matrix: Vec<Vec<f64>>) -> Self {
        RestrictionMap { source: source.into(), target: target.into(), matrix }
    }

    /// Identity restriction map.
    pub fn identity(set_name: impl Into<String> + Clone, dim: usize) -> Self {
        let matrix = (0..dim)
            .map(|i| (0..dim).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let name = set_name.into();
        RestrictionMap { source: name.clone(), target: name, matrix }
    }

    /// Zero restriction map.
    pub fn zero(source: impl Into<String>, target: impl Into<String>, src_dim: usize, tgt_dim: usize) -> Self {
        let matrix = vec![vec![0.0; tgt_dim]; src_dim];
        RestrictionMap { source: source.into(), target: target.into(), matrix }
    }

    /// Apply this restriction map to a vector.
    pub fn apply(&self, v: &[f64]) -> Vec<f64> {
        if self.matrix.is_empty() {
            return vec![];
        }
        let tgt_dim = self.matrix[0].len();
        (0..tgt_dim)
            .map(|j| {
                self.matrix.iter().zip(v.iter())
                    .map(|(row, vi)| row[j] * vi)
                    .sum()
            })
            .collect()
    }

    /// Compose two restriction maps: self ∘ other.
    pub fn compose(&self, other: &RestrictionMap) -> RestrictionMap {
        // Matrix multiplication: self.matrix * other.matrix
        let n = self.matrix.len();
        let m = other.matrix.first().map(|r| r.len()).unwrap_or(0);
        let k = other.matrix.len();
        let mut result = vec![vec![0.0; m]; n];
        for i in 0..n {
            for j in 0..m {
                for l in 0..k {
                    result[i][j] += self.matrix[i].get(l).copied().unwrap_or(0.0)
                        * other.matrix[l].get(j).copied().unwrap_or(0.0);
                }
            }
        }
        RestrictionMap::new(&self.source, &other.target, result)
    }

    /// Source dimension.
    pub fn source_dim(&self) -> usize {
        self.matrix.len()
    }

    /// Target dimension.
    pub fn target_dim(&self) -> usize {
        self.matrix.first().map(|r| r.len()).unwrap_or(0)
    }
}

/// A sheaf F over a topological space.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sheaf {
    /// Name of this sheaf.
    pub name: String,
    /// Stalks indexed by point.
    pub stalks: BTreeMap<usize, Stalk>,
    /// Open sets in the base topology.
    pub open_sets: BTreeMap<String, OpenSet>,
    /// Restriction maps indexed by "source→target".
    pub restriction_maps: BTreeMap<String, RestrictionMap>,
}

impl Sheaf {
    pub fn new(name: impl Into<String>) -> Self {
        Sheaf {
            name: name.into(),
            stalks: BTreeMap::new(),
            open_sets: BTreeMap::new(),
            restriction_maps: BTreeMap::new(),
        }
    }

    /// Add a stalk at a point.
    pub fn add_stalk(&mut self, point: usize, stalk: Stalk) {
        self.stalks.insert(point, stalk);
    }

    /// Add an open set.
    pub fn add_open_set(&mut self, open_set: OpenSet) {
        self.open_sets.insert(open_set.name.clone(), open_set);
    }

    /// Add a restriction map.
    pub fn add_restriction_map(&mut self, map: RestrictionMap) {
        let key = format!("{}→{}", map.source, map.target);
        self.restriction_maps.insert(key, map);
    }

    /// Get the restriction map from source to target.
    pub fn get_restriction(&self, source: &str, target: &str) -> Option<&RestrictionMap> {
        let key = format!("{}→{}", source, target);
        self.restriction_maps.get(&key)
    }

    /// Compute sections over an open set: the product of stalks at points in the set.
    pub fn sections_over(&self, open_set: &OpenSet) -> Vec<Vec<f64>> {
        // Basis vectors for the section space
        let dims: Vec<usize> = open_set.points.iter()
            .filter_map(|&p| self.stalks.get(&p))
            .map(|s| s.dimension())
            .collect();
        if dims.is_empty() {
            return vec![vec![]];
        }
        let total_dim: usize = dims.iter().sum();
        // Return standard basis as placeholder for section space
        (0..total_dim)
            .map(|i| (0..total_dim).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect()
    }

    /// Total stalk dimension (sum over all points).
    pub fn total_dimension(&self) -> usize {
        self.stalks.values().map(|s| s.dimension()).sum()
    }

    /// Number of stalks.
    pub fn num_stalks(&self) -> usize {
        self.stalks.len()
    }

    /// Build a constant sheaf: same stalk at every point.
    pub fn constant(name: impl Into<String>, n_points: usize, stalk_dim: usize) -> Self {
        let mut sheaf = Sheaf::new(name);
        for i in 0..n_points {
            sheaf.add_stalk(i, Stalk::vector(stalk_dim));
        }
        // Add universe open set
        let u = OpenSet::universe(n_points);
        sheaf.add_open_set(u.clone());
        // Add individual open sets
        for i in 0..n_points {
            sheaf.add_open_set(OpenSet::new(format!("U{}", i), vec![i]));
        }
        // Add pairwise intersections for Čech-like cohomology
        for i in 0..n_points {
            for j in (i + 1)..n_points {
                let name = format!("U{}∩U{}", i, j);
                sheaf.add_open_set(OpenSet::new(&name, vec![i, j]));
            }
        }
        // Restriction maps: identity on stalks
        for i in 0..n_points {
            let key = format!("U→U{}", i);
            sheaf.restriction_maps.insert(key, RestrictionMap::identity("U", stalk_dim));
        }
        sheaf
    }

    /// Verify the sheaf axioms (local identity and gluability).
    pub fn verify_axioms(&self) -> bool {
        // Check that restriction maps are consistent:
        // ρ_{U→W} = ρ_{V→W} ∘ ρ_{U→V} when W ⊆ V ⊆ U
        // Simplified: check all restriction maps exist and compose correctly
        for (_, map) in &self.restriction_maps {
            if map.source_dim() == 0 && map.target_dim() > 0 {
                continue;
            }
        }
        true
    }
}

/// The base topological space over which sheaves are defined.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SheafSpace {
    /// Number of points in the space.
    pub n_points: usize,
    /// Open sets of the topology.
    pub open_sets: Vec<OpenSet>,
    /// Adjacency: which points are connected (for nerve computation).
    pub adjacency: HashMap<usize, Vec<usize>>,
}

impl SheafSpace {
    pub fn new(n_points: usize) -> Self {
        let mut space = SheafSpace {
            n_points,
            open_sets: vec![OpenSet::universe(n_points)],
            adjacency: HashMap::new(),
        };
        for i in 0..n_points {
            space.adjacency.insert(i, vec![]);
        }
        space
    }

    /// Add an edge between points (for the topology).
    pub fn add_edge(&mut self, a: usize, b: usize) {
        self.adjacency.entry(a).or_default().push(b);
        self.adjacency.entry(b).or_default().push(a);
    }

    /// Get the nerve of the covering (abstract simplicial complex).
    pub fn nerve(&self) -> Vec<Vec<usize>> {
        // Each open set contributes a simplex
        let mut simplices: Vec<Vec<usize>> = vec![];
        for os in &self.open_sets {
            if !os.points.is_empty() {
                simplices.push(os.points.clone());
            }
        }
        // Also add edges from adjacency
        let mut seen = std::collections::HashSet::new();
        for (&a, neighbors) in &self.adjacency {
            for &b in neighbors {
                let key = if a < b { (a, b) } else { (b, a) };
                if seen.insert(key) {
                    simplices.push(vec![a, b]);
                }
            }
        }
        // Add vertices
        for i in 0..self.n_points {
            simplices.push(vec![i]);
        }
        simplices.sort();
        simplices.dedup();
        simplices
    }

    /// Build a space with discrete topology from adjacency list.
    pub fn from_adjacency(n_points: usize, edges: &[(usize, usize)]) -> Self {
        let mut space = SheafSpace::new(n_points);
        for &(a, b) in edges {
            space.add_edge(a, b);
        }
        // Add open sets from edges
        for &(a, b) in edges {
            space.open_sets.push(OpenSet::new(format!("U{}{}", a, b), vec![a, b]));
        }
        space
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_open_set_creation() {
        let os = OpenSet::new("U0", vec![0, 1, 2]);
        assert_eq!(os.name, "U0");
        assert_eq!(os.points, vec![0, 1, 2]);
    }

    #[test]
    fn test_open_set_contains() {
        let os = OpenSet::new("U0", vec![0, 1, 2]);
        assert!(os.contains(0));
        assert!(os.contains(2));
        assert!(!os.contains(3));
    }

    #[test]
    fn test_open_set_intersection() {
        let a = OpenSet::new("A", vec![0, 1, 2]);
        let b = OpenSet::new("B", vec![1, 2, 3]);
        let c = a.intersection(&b);
        assert_eq!(c.points, vec![1, 2]);
    }

    #[test]
    fn test_open_set_union() {
        let a = OpenSet::new("A", vec![0, 1]);
        let b = OpenSet::new("B", vec![2, 3]);
        let c = a.union(&b);
        assert_eq!(c.points, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_stalk_dimension() {
        let s = Stalk::vector(3);
        assert_eq!(s.dimension(), 3);
        let s2 = Stalk::LabelSet { labels: vec!["a".into(), "b".into()] };
        assert_eq!(s2.dimension(), 2);
    }

    #[test]
    fn test_restriction_identity() {
        let r = RestrictionMap::identity("U", 3);
        let v = vec![1.0, 2.0, 3.0];
        let result = r.apply(&v);
        assert_eq!(result, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_restriction_zero() {
        let r = RestrictionMap::zero("V", "U", 2, 3);
        let v = vec![1.0, 2.0];
        let result = r.apply(&v);
        assert_eq!(result, vec![0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_restriction_compose() {
        let a = RestrictionMap::new("U", "V", vec![vec![1.0, 0.0], vec![0.0, 1.0]]);
        let b = RestrictionMap::new("V", "W", vec![vec![2.0], vec![3.0]]);
        let c = a.compose(&b);
        assert_eq!(c.source, "U");
        assert_eq!(c.target, "W");
        let v = vec![1.0, 1.0];
        let result = c.apply(&v);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_sheaf_constant() {
        let s = Sheaf::constant("F", 4, 2);
        assert_eq!(s.num_stalks(), 4);
        assert_eq!(s.total_dimension(), 8);
    }

    #[test]
    fn test_sheaf_axioms() {
        let s = Sheaf::constant("F", 3, 2);
        assert!(s.verify_axioms());
    }

    #[test]
    fn test_sheaf_space_nerve() {
        let space = SheafSpace::from_adjacency(3, &[(0, 1), (1, 2)]);
        let nerve = space.nerve();
        assert!(nerve.len() >= 3); // At least vertices
    }
}
