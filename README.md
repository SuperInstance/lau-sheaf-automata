# lau-sheaf-automata

**Sheaf-theoretic protocol verification — deadlock-free iff H¹=0, composition as cup product.**

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)

---

## What This Does

This crate implements **Kimi's Theorem 2**: a protocol `P` is deadlock-free if and only if the first sheaf cohomology group vanishes:

```
P is deadlock-free  ⟺  H¹(Sh(States); L(P)) = 0
```

It provides:

- **Finite automata (DFA/NFA)** converted to sheaf stalks — language sheaves over agent configurations
- **Sheaf data structures** — stalks, restriction maps, open sets, sections
- **Sheaf cohomology computation** — H⁰ (global sections), H¹ (obstructions/deadlocks), H²
- **Protocol sheaves** — multi-agent protocols as sheaves over configuration spaces
- **Cup products** — protocol composition as `H¹ × H¹ → H²` in sheaf cohomology
- **Sheaf Laplacian** — efficient global section detection via spectral methods
- **Ergodic consensus** — distributed consensus as global section finding
- **Protocol refinement** — natural transformations between sheaf functors

The core workflow: define agents, build the protocol sheaf, compute cohomology, check if H¹ = 0.

---

## Key Idea

In a multi-agent protocol, each agent has a local state machine. The collection of all agent states forms a topological space (the **nerve** of the communication graph). A **sheaf** assigns data to each open set (agent neighborhood) with restriction maps for how data glues across overlaps.

| Sheaf Theory | Protocol Verification |
|---|---|
| Stalk at vertex `v` | Agent `v`'s local protocol state |
| Restriction map `ρ_{UV}` | Communication from agent `U` to neighbor `V` |
| Global section `s ∈ H⁰` | Consistent global protocol state |
| H¹ obstruction class | A deadlock — local states can't be patched globally |
| Cup product `∪: H¹ × H¹ → H²` | Composition of two protocols |
| Sheaf Laplacian `Δ_L` | Spectral tool for finding global sections |
| Refinement `P → Q` | Protocol specialization / abstraction |

If `H¹ ≠ 0`, there exist obstruction classes — local states that are individually valid but cannot be composed into a global consistent state. This is precisely a deadlock.

---

## Install

```toml
[dependencies]
lau-sheaf-automata = "0.1.0"
```

Requires Rust 2021 edition. Dependencies: `nalgebra`, `serde`.

---

## Quick Start

```rust
use lau_sheaf_automata::*;
use lau_sheaf_automata::sheaf::{OpenSet, Sheaf, Stalk, SheafSpace};

// --- Define agents ---
let agent1 = AgentConfig::new(0, "alice", vec!["idle".into(), "busy".into()])
    .with_neighbors(vec![1]);
let agent2 = AgentConfig::new(1, "bob", vec!["idle".into(), "busy".into()])
    .with_neighbors(vec![0]);

// --- Build protocol ---
let protocol = Protocol::new("handshake", vec![agent1, agent2]);

// --- Build protocol sheaf ---
let sheaf = protocol.build_sheaf();

// --- Compute cohomology ---
let cohomology = sheaf.compute_cohomology();
println!("H⁰ dimension (global sections): {}", cohomology.h0_dimension());
println!("H¹ dimension (deadlocks): {}", cohomology.h1_dimension());

// --- Check deadlock freedom ---
let deadlock_result = DeadlockResult::from_cohomology(&cohomology);
if deadlock_result.has_deadlock {
    println!("DEADLOCK DETECTED! {} obstruction classes",
        deadlock_result.obstruction_dimension);
} else {
    println!("Protocol is deadlock-free.");
}
```

### DFA → Sheaf Conversion

```rust
let mut dfa = Dfa::new(
    vec!["idle".into(), "request".into(), "done".into()],
    vec!["req".into(), "ack".into()],
    0,
    vec![false, false, true],
);
dfa.add_transition(0, 0, 1);  // idle --req--> request
dfa.add_transition(1, 1, 2);  // request --ack--> done

// Convert to sheaf stalk
let stalk = dfa.to_stalk();
```

### Protocol Composition via Cup Product

```rust
// Compose two protocols — their obstruction interaction lives in H²
let cup = CupProduct::new(h1_dim_a, h1_dim_b, h2_dim);
let combined_obstruction = cup.compute(&obstruction_a, &obstruction_b);
```

### Consensus via Sheaf Laplacian

```rust
let consensus = ErgodicConsensus::new()
    .with_tolerance(1e-8);
let result = consensus.compute(&laplacian, &initial_sections);
if result.reached {
    println!("Consensus reached in {} iterations", result.iterations);
}
```

---

## API Reference

### Core Types

| Type | Module | Description |
|---|---|---|
| `Dfa` | `automata` | Deterministic finite automaton |
| `Nfa` | `automata` | Nondeterministic finite automaton |
| `Stalk` | `sheaf` | Vector space attached to an open set |
| `RestrictionMap` | `sheaf` | Linear map between stalks on nested open sets |
| `OpenSet` | `sheaf` | Named subset of the base space |
| `Sheaf` | `sheaf` | Sheaf: stalks + restriction maps + topology |
| `SheafSpace` | `sheaf` | Topological base space with open cover |
| `CochainComplex` | `cohomology` | C⁰ → C¹ → C² complex |
| `Cohomology` | `cohomology` | Computed H⁰, H¹, H² groups |
| `AgentConfig` | `protocol` | Agent with local states and neighbors |
| `Protocol` | `protocol` | Multi-agent protocol with transition rules |
| `ProtocolSheaf` | `protocol` | Sheaf derived from a protocol |
| `CupProduct` | `cup_product` | Bilinear map H¹ × H¹ → H² |
| `SheafLaplacian` | `laplacian` | Spectral operator for global sections |
| `ErgodicConsensus` | `consensus` | Iterative consensus via Laplacian |
| `Refinement` | `refinement` | Natural transformation between sheaves |
| `DeadlockResult` | `lib` | Deadlock detection from cohomology |

### Key Methods

**Dfa / Nfa**
- `Dfa::new(states, alphabet, initial, accepting)` — create automaton
- `dfa.add_transition(from, symbol, to)` — add transition
- `dfa.run(&word)` — execute on input, get final state
- `dfa.accepts(&word)` — membership test
- `dfa.to_stalk()` — convert to sheaf stalk
- `Nfa::epsilon_closure(&state)` — ε-closure computation

**Sheaf**
- `Sheaf::new(space)` — create empty sheaf over a space
- `sheaf.add_stalk(open_set, stalk)` — attach data
- `sheaf.add_restriction(source, target, map)` — add restriction map
- `sheaf.section_at(&open_set)` — retrieve local section
- `sheaf.compute_cohomology()` → `Cohomology` — full cohomology computation

**Cohomology**
- `cohomology.h0_dimension()` — dim(H⁰) = number of global sections
- `cohomology.h1_dimension()` — dim(H¹) = number of obstruction classes
- `cohomology.h2_dimension()` — dim(H²) = cup product target
- `cohomology.h0_representatives()` — basis for global sections
- `cohomology.h1_representatives()` — basis for obstruction classes

**Protocol**
- `Protocol::new(name, agents)` — create protocol
- `protocol.build_sheaf()` → `ProtocolSheaf` — derive sheaf from agents
- `protocol.check_deadlock()` → `DeadlockResult` — via cohomology

**CupProduct**
- `CupProduct::new(h1_dim_a, h1_dim_b, h2_dim)` — create bilinear map
- `cup.compute(&class_a, &class_b)` → `Vec<f64>` — cup product of obstruction classes
- `cup.is_graded_commutative()` — verify `a ∪ b = (-1)^{pq} b ∪ a`

**SheafLaplacian**
- `SheafLaplacian::from_sheaf_data(n, stalk_dims, edges, restriction_norms)` — build matrix
- `laplacian.global_sections()` — compute ker(Δ_L) = H⁰
- `laplacian.spectrum()` — eigenvalues for spectral analysis

**ErgodicConsensus**
- `ErgodicConsensus::new()` — default consensus solver
- `consensus.compute(&laplacian, &initial)` → `ConsensusResult` — iterate to global section

**Refinement**
- `Refinement::new(source, target, component_maps)` — natural transformation
- `refinement.induced_map_on_cohomology(&h_source, &h_target)` — cohomology map
- `refinement.is_monic` / `is_epic` — injective/surjective check

---

## How It Works

### 1. Automata as Sheaf Stalks

Each agent's local automaton (DFA/NFA) defines a **stalk** — a vector space whose basis elements are the automaton states. Multiple agents give a product of stalks over the communication graph.

### 2. Communication as Restriction Maps

When agents communicate, their states must be compatible. Restriction maps encode this compatibility: the map `ρ_{UV}` specifies how agent `U`'s local state constrains neighbor `V`'s state. These maps form a **presheaf** functor from the nerve of the communication graph to vector spaces.

### 3. Cohomology = Deadlock Detection

The Čech cochain complex `C⁰ → C¹ → C²` is built from the sheaf:
- `C⁰` = assignments of states to individual agents
- `C¹` = assignments to pairwise overlaps (communication channels)
- `C²` = assignments to triple overlaps
- `d₀` = coboundary checking pairwise consistency
- `d₁` = coboundary checking triple consistency

Then:
- **H⁰ = ker(d₀)**: global sections — globally consistent state assignments
- **H¹ = ker(d₁)/im(d₀)**: obstruction classes — deadlocks

### 4. Cup Product = Protocol Composition

When two protocols `P` and `Q` are composed, their obstruction classes interact via the cup product `∪: H¹(P) × H¹(Q) → H²(P ⊗ Q)`. A non-zero result means the composed protocol has higher-order obstructions.

### 5. Sheaf Laplacian

The sheaf Laplacian `Δ_L = D - A_sheaf` is a discrete operator whose kernel equals H⁰ (global sections). It enables efficient spectral methods for finding consistent states without enumerating the full cochain complex.

### 6. Consensus via Laplacian Flow

Distributed consensus is equivalent to finding a global section. The Laplacian flow `dx/dt = -Δ_L x` converges to the nearest global section, providing an iterative distributed algorithm.

---

## The Math

### Sheaf Cohomology

Given a sheaf `F` on a topological space `X` with open cover `{Uᵢ}`:

```
Čech cochain complex:  C⁰ →[d⁰] C¹ →[d¹] C² → ...

C⁰ = Πᵢ F(Uᵢ)         (local sections)
C¹ = Πᵢ<ⱼ F(Uᵢ ∩ Uⱼ)   (pairwise compatibility data)
C² = Πᵢ<j<k F(Uᵢ ∩ Uⱼ ∩ Uₖ)  (triple compatibility)

H⁰(X; F) = ker(d⁰)     (global sections)
H¹(X; F) = ker(d¹)/im(d⁰)  (obstruction classes)
```

### Deadlock Theorem

```
P is deadlock-free ⟺ H¹(Sh(States); L(P)) = 0
```

**Proof sketch**: H¹ measures the failure of local data to patch into global data. If H¹ ≠ 0, there exist local states that are pairwise compatible but not globally consistent — a deadlock. If H¹ = 0, all pairwise compatibilities extend to global consistency.

### Cup Product

The cup product is the bilinear map:
```
∪: H^p(F) × H^q(G) → H^{p+q}(F ⊗ G)
```

Defined on cochains by:
```
(α ∪ β)(U₀...U_{p+q}) = α(U₀...U_p) · ρ(β(U_p...U_{p+q}))
```

It is graded-commutative: `α ∪ β = (-1)^{pq} β ∪ α`.

### Sheaf Laplacian

For a cellular sheaf on a graph with vertices `V` and edges `E`:
```
Δ_L = L_↓ + L_↑

where:
  (L_↓)ᵥ = Σ_{e:v∼w} ρ_{v,e}ᵀ ρ_{v,e}    (restriction norm²)
  (L_↑)ᵥ = Σ_{e:v∼w} ρ_{w,e}ᵀ ρ_{w,e}    (co-restriction norm²)
```

The kernel of `Δ_L` equals the global sections: `ker(Δ_L) = H⁰(X; F)`.

---

## License

MIT
