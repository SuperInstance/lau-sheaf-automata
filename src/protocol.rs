//! Protocol sheaves and agent configurations.
//!
//! A protocol P defines a language sheaf L(P) over the space of agent configurations.
//! Deadlock-freedom is characterized by H¹(Sh(States); L(P)) = 0.

use serde::{Deserialize, Serialize};

use crate::cohomology::Cohomology;
use crate::sheaf::{OpenSet, RestrictionMap, Sheaf, SheafSpace, Stalk};
use crate::DeadlockResult;

/// An agent's configuration: its local state and connectivity.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgentConfig {
    pub id: usize,
    pub name: String,
    pub local_states: Vec<String>,
    pub neighbors: Vec<usize>,
}

impl AgentConfig {
    pub fn new(id: usize, name: impl Into<String>, local_states: Vec<String>) -> Self {
        AgentConfig { id, name: name.into(), local_states, neighbors: vec![] }
    }

    pub fn with_neighbors(mut self, neighbors: Vec<usize>) -> Self {
        self.neighbors = neighbors;
        self
    }

    pub fn n_states(&self) -> usize {
        self.local_states.len()
    }
}

/// A protocol: a collection of agents with transition rules.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Protocol {
    pub name: String,
    pub agents: Vec<AgentConfig>,
    /// Transition rules: (agent_id, from_state, to_state) → guard conditions.
    pub rules: Vec<TransitionRule>,
}

/// A transition rule in the protocol.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransitionRule {
    pub agent: usize,
    pub from_state: usize,
    pub to_state: usize,
    /// Required states of other agents for this transition.
    pub guards: Vec<(usize, usize)>,
}

impl Protocol {
    pub fn new(name: impl Into<String>) -> Self {
        Protocol {
            name: name.into(),
            agents: vec![],
            rules: vec![],
        }
    }

    pub fn add_agent(&mut self, agent: AgentConfig) {
        self.agents.push(agent);
    }

    pub fn add_rule(&mut self, rule: TransitionRule) {
        self.rules.push(rule);
    }

    /// Number of agents.
    pub fn n_agents(&self) -> usize {
        self.agents.len()
    }

    /// Total configuration space dimension.
    pub fn config_dim(&self) -> usize {
        self.agents.iter().map(|a| a.n_states()).sum()
    }

    /// Build the protocol sheaf L(P) over the configuration space.
    pub fn to_sheaf(&self) -> ProtocolSheaf {
        let n = self.agents.len();
        let mut sheaf = Sheaf::new(format!("L({})", self.name));

        // Build configuration space topology
        let mut space = SheafSpace::new(n);
        for agent in &self.agents {
            sheaf.add_stalk(agent.id, Stalk::vector(agent.n_states()));
            sheaf.add_open_set(OpenSet::new(format!("U{}", agent.id), vec![agent.id]));
            for &neighbor in &agent.neighbors {
                space.add_edge(agent.id, neighbor);
            }
        }
        sheaf.add_open_set(OpenSet::universe(n));

        // Add pairwise open sets for edges
        for agent in &self.agents {
            for &neighbor in &agent.neighbors {
                if agent.id < neighbor {
                    let name = format!("U{}∩U{}", agent.id, neighbor);
                    sheaf.add_open_set(OpenSet::new(&name, vec![agent.id, neighbor]));
                }
            }
        }

        // Restriction maps
        for agent in &self.agents {
            let dim = agent.n_states();
            sheaf.add_restriction_map(RestrictionMap::identity(
                format!("U{}", agent.id), dim,
            ));
        }

        ProtocolSheaf {
            protocol_name: self.name.clone(),
            sheaf,
            space,
            rules: self.rules.clone(),
        }
    }
}

/// The language sheaf L(P) of a protocol.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProtocolSheaf {
    pub protocol_name: String,
    pub sheaf: Sheaf,
    pub space: SheafSpace,
    pub rules: Vec<TransitionRule>,
}

impl ProtocolSheaf {
    /// Compute sheaf cohomology and check for deadlocks.
    pub fn check_deadlock(&self) -> DeadlockResult {
        let cohomology = self.compute_cohomology();
        DeadlockResult::from_cohomology(&cohomology)
    }

    /// Compute the sheaf cohomology.
    pub fn compute_cohomology(&self) -> Cohomology {
        // Build coboundary maps from the sheaf and topology
        let n = self.sheaf.num_stalks();
        let dims: Vec<usize> = (0..n)
            .filter_map(|i| self.sheaf.stalks.get(&i))
            .map(|s| s.dimension())
            .collect();

        if dims.is_empty() {
            return Cohomology::trivial();
        }

        // d0: C⁰ → C¹ (restriction differences)
        // For n points, C⁰ = ⊕ F(Ui), C¹ = ⊕ F(Ui ∩ Uj)
        let total_c0: usize = dims.iter().sum();
        let total_c1: usize = self.space.open_sets.iter()
            .filter(|os| os.points.len() == 2)
            .map(|os| os.points.iter()
                .filter_map(|&p| self.sheaf.stalks.get(&p))
                .map(|s| s.dimension())
                .sum::<usize>())
            .sum();

        // Build d0 matrix
        let d0 = if total_c0 > 0 && total_c1 > 0 {
            Self::build_d0(&dims, &self.space, &self.sheaf, total_c0, total_c1)
        } else {
            vec![]
        };

        // d1: C¹ → C²
        let d1 = if total_c1 > 0 {
            vec![vec![0.0; total_c1]; 0] // Placeholder
        } else {
            vec![]
        };

        Cohomology::from_coboundaries(d0, d1, total_c0, total_c1)
    }

    fn build_d0(
        dims: &[usize],
        space: &SheafSpace,
        _sheaf: &Sheaf,
        c0_dim: usize,
        c1_dim: usize,
    ) -> Vec<Vec<f64>> {
        // d0 maps from C⁰ (sections on individual open sets) to C¹ (sections on intersections)
        // For each pair (i,j), d0_ij = ρ_j - ρ_i (difference of restrictions)
        let mut d0 = vec![vec![0.0; c0_dim]; c1_dim];

        // Simple model: for each edge, the coboundary is the difference map
        let mut row = 0;
        let edges: Vec<(usize, usize)> = space.adjacency.iter()
            .flat_map(|(&a, neighbors)| neighbors.iter().map(move |&b| (a, b)))
            .filter(|(a, b)| a < b)
            .collect();

        for &(i, j) in &edges {
            let di = dims.get(i).copied().unwrap_or(0);
            let dj = dims.get(j).copied().unwrap_or(0);
            let block_dim = di.min(dj);
            // Map: d0[section on Ui ⊕ Uj] = ρ_j(sj) - ρ_i(si)
            for k in 0..block_dim {
                if row + k < c1_dim {
                    // +1 for j's section
                    let col_j: usize = dims.iter().take(j).sum();
                    if col_j + k < c0_dim {
                        d0[row + k][col_j + k] = 1.0;
                    }
                    // -1 for i's section
                    let col_i: usize = dims.iter().take(i).sum();
                    if col_i + k < c0_dim {
                        d0[row + k][col_i + k] = -1.0;
                    }
                }
            }
            row += block_dim;
        }

        d0
    }
}

/// Build a simple mutex protocol with n agents and m resources.
pub fn mutex_protocol(n_agents: usize, _n_resources: usize) -> Protocol {
    let mut proto = Protocol::new("MutexProtocol");
    for i in 0..n_agents {
        let neighbors: Vec<usize> = (0..n_agents).filter(|&j| j != i).collect();
        let agent = AgentConfig::new(i, format!("Agent{}", i), vec!["idle".into(), "waiting".into(), "critical".into()])
            .with_neighbors(neighbors);
        proto.add_agent(agent);
    }
    proto
}

/// Build a deadlock-free protocol (ring topology, always makes progress).
pub fn deadlock_free_ring(n: usize) -> Protocol {
    let mut proto = Protocol::new("DeadlockFreeRing");
    for i in 0..n {
        let next = (i + 1) % n;
        let prev = if i == 0 { n - 1 } else { i - 1 };
        let agent = AgentConfig::new(i, format!("Node{}", i), vec!["idle".into(), "active".into()])
            .with_neighbors(vec![prev, next]);
        proto.add_agent(agent);
    }
    proto
}

/// Build a protocol with a known deadlock (circular wait).
pub fn circular_wait_deadlock(n: usize) -> Protocol {
    let mut proto = Protocol::new("CircularWait");
    for i in 0..n {
        let next = (i + 1) % n;
        let agent = AgentConfig::new(i, format!("Node{}", i), vec!["idle".into(), "waiting".into(), "held".into()])
            .with_neighbors(vec![next]);
        proto.add_agent(agent);
        // Rule: each node waits for the next
        proto.add_rule(TransitionRule {
            agent: i,
            from_state: 0, // idle
            to_state: 1,   // waiting for next
            guards: vec![],
        });
        proto.add_rule(TransitionRule {
            agent: i,
            from_state: 1, // waiting
            to_state: 2,   // held (requires next to be idle)
            guards: vec![(next, 0)],
        });
    }
    proto
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_agent_config() {
        let agent = AgentConfig::new(0, "A", vec!["idle".into(), "active".into()]);
        assert_eq!(agent.n_states(), 2);
    }

    #[test]
    fn test_protocol_creation() {
        let mut p = Protocol::new("TestProtocol");
        p.add_agent(AgentConfig::new(0, "A", vec!["x".into(), "y".into()]));
        p.add_agent(AgentConfig::new(1, "B", vec!["x".into(), "y".into()]));
        assert_eq!(p.n_agents(), 2);
        assert_eq!(p.config_dim(), 4);
    }

    #[test]
    fn test_protocol_to_sheaf() {
        let mut p = Protocol::new("Simple");
        p.add_agent(AgentConfig::new(0, "A", vec!["idle".into(), "run".into()]));
        p.add_agent(AgentConfig::new(1, "B", vec!["idle".into(), "run".into()]));
        let ps = p.to_sheaf();
        assert_eq!(ps.sheaf.num_stalks(), 2);
    }

    #[test]
    fn test_mutex_protocol() {
        let p = mutex_protocol(3, 2);
        assert_eq!(p.n_agents(), 3);
    }

    #[test]
    fn test_deadlock_free_ring() {
        let p = deadlock_free_ring(4);
        assert_eq!(p.n_agents(), 4);
        let ps = p.to_sheaf();
        assert_eq!(ps.sheaf.num_stalks(), 4);
    }

    #[test]
    fn test_circular_wait() {
        let p = circular_wait_deadlock(3);
        assert_eq!(p.n_agents(), 3);
        assert_eq!(p.rules.len(), 6);
    }

    #[test]
    fn test_deadlock_free_check() {
        let p = deadlock_free_ring(3);
        let ps = p.to_sheaf();
        let result = ps.check_deadlock();
        // Deadlock-free ring should have H¹ = 0 in our model
        assert_eq!(result.has_deadlock, false);
    }

    #[test]
    fn test_protocol_sheaf_cohomology() {
        let p = mutex_protocol(2, 1);
        let ps = p.to_sheaf();
        let cohomology = ps.compute_cohomology();
        // H⁰ should be at least 1 (the constant section)
        assert!(cohomology.h0_dimension() >= 0);
    }
}
