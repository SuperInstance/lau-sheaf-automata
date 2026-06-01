//! # lau-sheaf-automata
//!
//! Sheaf-theoretic protocol verification implementing Kimi's Theorem 2:
//! A protocol P is deadlock-free if and only if H¹(Sh(States); L(P)) = 0.
//!
//! Protocol composition is modeled as the cup product in sheaf cohomology,
//! and deadlock detection corresponds to finding non-zero obstruction classes
//! in H¹.

pub mod automata;
pub mod cohomology;
pub mod consensus;
pub mod cup_product;
pub mod laplacian;
pub mod protocol;
pub mod refinement;
pub mod sheaf;

pub use automata::{Dfa, Nfa};
pub use cohomology::{CochainComplex, Cohomology};
pub use consensus::{ConsensusResult, ErgodicConsensus};
pub use cup_product::CupProduct;
pub use laplacian::SheafLaplacian;
pub use protocol::{AgentConfig, Protocol, ProtocolSheaf};
pub use refinement::Refinement;
pub use sheaf::{OpenSet, RestrictionMap, Sheaf, SheafSpace, Stalk};

/// Result type alias for this crate.
pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

/// Deadlock detection result.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DeadlockResult {
    /// Whether a deadlock was detected (H¹ ≠ 0).
    pub has_deadlock: bool,
    /// Dimension of H¹ (number of independent obstruction classes).
    pub obstruction_dimension: usize,
    /// Description of obstruction classes.
    pub obstructions: Vec<String>,
}

impl DeadlockResult {
    /// Check deadlock freedom via Kimi's Theorem 2.
    pub fn from_cohomology(cohomology: &Cohomology) -> Self {
        let h1_dim = cohomology.h1_dimension();
        let obstructions: Vec<String> = (0..h1_dim)
            .map(|i| format!("obstruction_class_{}", i))
            .collect();
        DeadlockResult {
            has_deadlock: h1_dim > 0,
            obstruction_dimension: h1_dim,
            obstructions,
        }
    }
}
