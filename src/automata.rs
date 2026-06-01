//! Automata (DFA/NFA) → sheaf stalks conversion.
//!
//! Language sheaf L(P) maps agent configurations (automaton states) to stalks.

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

use crate::sheaf::{OpenSet, Sheaf, Stalk};

/// A deterministic finite automaton.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Dfa {
    pub states: Vec<String>,
    pub alphabet: Vec<String>,
    /// transitions[state_index][symbol_index] → target state index.
    pub transitions: BTreeMap<(usize, usize), usize>,
    pub initial: usize,
    pub accepting: Vec<bool>,
}

impl Dfa {
    pub fn new(
        states: Vec<String>,
        alphabet: Vec<String>,
        initial: usize,
        accepting: Vec<bool>,
    ) -> Self {
        Dfa { states, alphabet, transitions: BTreeMap::new(), initial, accepting }
    }

    pub fn add_transition(&mut self, from: usize, symbol: usize, to: usize) {
        self.transitions.insert((from, symbol), to);
    }

    /// Run the DFA on a word, returning final state index (or None if stuck).
    pub fn run(&self, word: &[usize]) -> Option<usize> {
        let mut current = self.initial;
        for &sym in word {
            current = self.transitions.get(&(current, sym)).copied()?;
        }
        Some(current)
    }

    /// Check if a word is accepted.
    pub fn accepts(&self, word: &[usize]) -> bool {
        self.run(word)
            .map(|s| self.accepting.get(s).copied().unwrap_or(false))
            .unwrap_or(false)
    }

    /// Number of states.
    pub fn n_states(&self) -> usize {
        self.states.len()
    }

    /// Convert this DFA into a sheaf stalk at a point.
    pub fn to_stalk(&self) -> Stalk {
        Stalk::AutomatonState {
            states: (0..self.states.len()).collect(),
            accepting: self.accepting.clone(),
        }
    }
}

/// A non-deterministic finite automaton.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Nfa {
    pub states: Vec<String>,
    pub alphabet: Vec<String>,
    /// transitions[(state, symbol)] → set of target states.
    pub transitions: BTreeMap<(usize, usize), BTreeSet<usize>>,
    pub initial: usize,
    pub accepting: Vec<bool>,
}

impl Nfa {
    pub fn new(
        states: Vec<String>,
        alphabet: Vec<String>,
        initial: usize,
        accepting: Vec<bool>,
    ) -> Self {
        Nfa { states, alphabet, transitions: BTreeMap::new(), initial, accepting }
    }

    pub fn add_transition(&mut self, from: usize, symbol: usize, to: usize) {
        self.transitions.entry((from, symbol)).or_default().insert(to);
    }

    /// Run the NFA, returning all reachable states.
    pub fn run(&self, word: &[usize]) -> BTreeSet<usize> {
        let mut current: BTreeSet<usize> = vec![self.initial].into_iter().collect();
        for &sym in word {
            let mut next = BTreeSet::new();
            for &s in &current {
                if let Some(targets) = self.transitions.get(&(s, sym)) {
                    next.extend(targets);
                }
            }
            current = next;
        }
        current
    }

    /// Check if a word is accepted.
    pub fn accepts(&self, word: &[usize]) -> bool {
        self.run(word).iter().any(|&s| self.accepting.get(s).copied().unwrap_or(false))
    }

    /// Subset construction: convert to DFA.
    pub fn to_dfa(&self) -> Dfa {
        let mut dfa_states: Vec<BTreeSet<usize>> = vec![];
        let mut dfa_names: Vec<String> = vec![];
        let mut dfa_accepting: Vec<bool> = vec![];
        let mut dfa_transitions: BTreeMap<(usize, usize), usize> = BTreeMap::new();
        let mut visited: BTreeMap<BTreeSet<usize>, usize> = BTreeMap::new();

        let initial_set: BTreeSet<usize> = vec![self.initial].into_iter().collect();
        dfa_states.push(initial_set.clone());
        dfa_names.push(format!("{:?}", initial_set));
        dfa_accepting.push(initial_set.iter().any(|&s| self.accepting.get(s).copied().unwrap_or(false)));
        visited.insert(initial_set, 0);

        let mut queue = vec![0usize];
        while let Some(dfa_state_idx) = queue.pop() {
            let nfa_state_set: BTreeSet<usize> = dfa_states[dfa_state_idx].clone();
            for (sym_idx, _sym) in self.alphabet.iter().enumerate() {
                let mut next: BTreeSet<usize> = BTreeSet::new();
                for &s in &nfa_state_set {
                    if let Some(targets) = self.transitions.get(&(s, sym_idx)) {
                        next.extend(targets);
                    }
                }
                if next.is_empty() {
                    continue;
                }
                if let Some(&existing) = visited.get(&next) {
                    dfa_transitions.insert((dfa_state_idx, sym_idx), existing);
                } else {
                    let new_idx = dfa_states.len();
                    dfa_states.push(next.clone());
                    dfa_names.push(format!("{:?}", next));
                    dfa_accepting.push(next.iter().any(|&s| self.accepting.get(s).copied().unwrap_or(false)));
                    visited.insert(next, new_idx);
                    dfa_transitions.insert((dfa_state_idx, sym_idx), new_idx);
                    queue.push(new_idx);
                }
            }
        }

        let initial = 0;
        Dfa {
            states: dfa_names,
            alphabet: self.alphabet.clone(),
            transitions: dfa_transitions,
            initial,
            accepting: dfa_accepting,
        }
    }

    /// Convert to sheaf stalk.
    pub fn to_stalk(&self) -> Stalk {
        Stalk::AutomatonState {
            states: (0..self.states.len()).collect(),
            accepting: self.accepting.clone(),
        }
    }
}

/// Build a language sheaf from a collection of automata indexed by agent.
pub fn build_language_sheaf(
    name: &str,
    automata: &[(usize, &Dfa)],
) -> Sheaf {
    let mut sheaf = Sheaf::new(name);
    for &(agent_idx, dfa) in automata {
        sheaf.add_stalk(agent_idx, dfa.to_stalk());
        sheaf.add_open_set(OpenSet::new(format!("U{}", agent_idx), vec![agent_idx]));
    }
    // Universe
    let n = automata.len();
    sheaf.add_open_set(OpenSet::universe(n));
    // Add restriction maps (identity for constant sheaf behavior)
    for &(agent_idx, dfa) in automata {
        let dim = dfa.n_states();
        sheaf.add_restriction_map(crate::sheaf::RestrictionMap::identity(
            format!("U{}", agent_idx), dim,
        ));
    }
    sheaf
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_dfa() -> Dfa {
        let mut dfa = Dfa::new(
            vec!["q0".into(), "q1".into()],
            vec!["a".into(), "b".into()],
            0,
            vec![false, true],
        );
        dfa.add_transition(0, 0, 1); // q0 --a--> q1
        dfa.add_transition(0, 1, 0); // q0 --b--> q0
        dfa.add_transition(1, 0, 0); // q1 --a--> q0
        dfa.add_transition(1, 1, 1); // q1 --b--> q1
        dfa
    }

    #[test]
    fn test_dfa_run() {
        let dfa = simple_dfa();
        assert_eq!(dfa.run(&[0]), Some(1)); // "a" → q1
        assert_eq!(dfa.run(&[1]), Some(0)); // "b" → q0
    }

    #[test]
    fn test_dfa_accepts() {
        let dfa = simple_dfa();
        assert!(dfa.accepts(&[0]));      // "a"
        assert!(!dfa.accepts(&[]));       // ε
        assert!(dfa.accepts(&[0, 1]));    // "ab"
    }

    #[test]
    fn test_dfa_to_stalk() {
        let dfa = simple_dfa();
        let stalk = dfa.to_stalk();
        assert_eq!(stalk.dimension(), 2);
    }

    #[test]
    fn test_nfa_run() {
        let mut nfa = Nfa::new(
            vec!["q0".into(), "q1".into()],
            vec!["a".into()],
            0,
            vec![false, true],
        );
        nfa.add_transition(0, 0, 0);
        nfa.add_transition(0, 0, 1);
        nfa.add_transition(1, 0, 1);

        let result = nfa.run(&[0]);
        assert!(result.contains(&0));
        assert!(result.contains(&1));
    }

    #[test]
    fn test_nfa_accepts() {
        let mut nfa = Nfa::new(
            vec!["q0".into(), "q1".into()],
            vec!["a".into()],
            0,
            vec![false, true],
        );
        nfa.add_transition(0, 0, 1);
        assert!(nfa.accepts(&[0]));
    }

    #[test]
    fn test_nfa_to_dfa() {
        let mut nfa = Nfa::new(
            vec!["q0".into(), "q1".into()],
            vec!["a".into()],
            0,
            vec![false, true],
        );
        nfa.add_transition(0, 0, 0);
        nfa.add_transition(0, 0, 1);
        nfa.add_transition(1, 0, 1);

        let dfa = nfa.to_dfa();
        assert!(dfa.accepts(&[0])); // "a" should be accepted
    }

    #[test]
    fn test_build_language_sheaf() {
        let dfa = simple_dfa();
        let sheaf = build_language_sheaf("L(P)", &[(0, &dfa), (1, &dfa)]);
        assert_eq!(sheaf.num_stalks(), 2);
        assert_eq!(sheaf.total_dimension(), 4);
    }
}
