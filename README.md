# lau-sheaf-automata

**Sheaf-theoretic protocol verification — deadlock-free iff H¹ = 0, composition as cup product.**

A Rust library implementing Kimi's Theorem 2: a protocol P is deadlock-free if and only if H¹(Sh(States); L(P)) = 0. Protocol composition is modeled as the **cup product** in sheaf cohomology, and deadlock detection corresponds to finding non-zero **obstruction classes** in H¹.

[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)

---

## What This Does

In concurrent and distributed systems, **deadlock** — where agents cyclically wait on each other forever — is notoriously hard to detect. Traditional approaches (model checking, type systems) don't scale well to large systems with many agents.

This crate takes a **topological** approach:

1. Model a protocol's configuration space as a **sheaf** — a mathematical object that assigns local data (stalks) to regions of a space and connects them with restriction maps
2. Compute the **sheaf cohomology** H⁰ and H¹ using coboundary maps and linear algebra
3. **Deadlock ≡ H¹ ≠ 0**: a non-trivial cohomology class is an "obstruction" to deadlock freedom

Additionally:

- **Protocol composition** (running two protocols together) corresponds to the **cup product** ⌣: H¹ × H¹ → H²
- **Refinement** (making a protocol more specific) is a **natural transformation** between sheaf functors
- **Consensus** is detected via **ergodic theory** on the sheaf's state space
- The **sheaf Laplacian** provides spectral verification of global consistency

---

## Key Idea

A **cellular sheaf** F on a topological space X assigns:
- A **stalk** F(U) (a vector space) to each open set U
- A **restriction map** ρ_{V→U}: F(V) → F(U) for V ⊆ U

For a protocol P with agents {A₁, …, Aₙ}, we build the **language sheaf** L(P):
- Stalk at agent Aᵢ = its local state space
- Stalk at intersection Uᵢ ∩ Uⱼ = joint states visible to both agents
- Restriction maps = projections from joint to local states

The **Čech cohomology** of this sheaf is computed via coboundary maps:

> d⁰: C⁰ → C¹ (local → pairwise)  
> d¹: C¹ → C² (pairwise → triple intersections)

Then:
- H⁰ = ker(d⁰) = **global sections** (consistent configurations)
- H¹ = ker(d¹) / im(d⁰) = **obstructions** (deadlock classes)

**Kimi's Theorem 2**: P is deadlock-free ⟺ H¹ = 0.

---

## Install

```toml
[dependencies]
lau-sheaf-automata = "0.1"
```

---

## Quick Start

```rust
use lau_sheaf_automata::{
    Protocol, AgentConfig, Sheaf, Stalk, OpenSet, RestrictionMap,
    ConsensusResult, CupProduct, SheafLaplacian, DeadlockResult,
};

// --- Define a protocol with 3 agents ---
let mut protocol = Protocol::new("MutexProtocol");
for i in 0..3 {
    let neighbors: Vec<usize> = (0..3).filter(|&j| j != i).collect();
    let agent = AgentConfig::new(i, format!("Agent{}", i),
        vec!["idle".into(), "waiting".into(), "critical".into()])
        .with_neighbors(neighbors);
    protocol.add_agent(agent);
}

// --- Build the language sheaf L(P) ---
let protocol_sheaf = protocol.to_sheaf();

// --- Check for deadlocks ---
let result: DeadlockResult = protocol_sheaf.check_deadlock();
if result.has_deadlock {
    println!("DEADLOCK DETECTED! {} obstruction classes", result.obstruction_dimension);
    for obs in &result.obstructions {
        println!("  - {}", obs);
    }
} else {
    println!("Protocol is deadlock-free (H¹ = 0)");
}

// --- Protocol composition via cup product ---
// Composing two protocols: H¹(P₁) ⌣ H¹(P₂) → H²(P₁ ∘ P₂)

// --- Sheaf Laplacian for spectral verification ---
let lap = SheafLaplacian::from_sheaf_data(
    3,                    // 3 nodes
    vec![2, 2, 2],        // stalk dimensions
    &[(0, 1), (1, 2)],    // edges
    &[(0, 1, 1.0), (1, 2, 1.0)], // restriction weights
);
let kernel = lap.kernel();          // Global sections (H⁰)
let gap = lap.spectral_gap();       // Measures consensus speed
```

---

## API Reference

### `Sheaf` — Cellular Sheaf Data Structure

```rust
pub struct Sheaf {
    pub name: String,
    pub stalks: BTreeMap<usize, Stalk>,
    pub open_sets: BTreeMap<String, OpenSet>,
    pub restriction_maps: BTreeMap<String, RestrictionMap>,
}
```

| Method | Description |
|--------|-------------|
| `new(name)` | Create an empty sheaf |
| `add_stalk(point, stalk)` | Assign a stalk at a point |
| `add_open_set(open_set)` | Add an open set to the topology |
| `add_restriction_map(map)` | Add ρ_{V→U}: F(V) → F(U) |
| `get_restriction(source, target)` | Retrieve a restriction map |
| `sections_over(open_set)` | Compute section space basis |
| `total_dimension() → usize` | Sum of all stalk dimensions |
| `constant(name, n, dim)` | Build a constant sheaf (same stalk everywhere) |
| `verify_axioms() → bool` | Check consistency of restriction maps |

### `Stalk` — Local Data at a Point

```rust
pub enum Stalk {
    VectorSpace { dimension: usize },
    LabelSet { labels: Vec<String> },
    AutomatonState { states: Vec<usize>, accepting: Vec<bool> },
    Custom { data: Vec<u8>, dimension: usize },
}
```

### `OpenSet` — Region of the Base Space

| Method | Description |
|--------|-------------|
| `new(name, points)` | Create an open set |
| `universe(n)` | The full space |
| `intersection(other)` | Set intersection |
| `union(other)` | Set union |
| `contains(point)` | Membership test |

### `RestrictionMap` — Linear Map Between Stalks

| Method | Description |
|--------|-------------|
| `new(source, target, matrix)` | Create from a matrix |
| `identity(name, dim)` | Identity map |
| `zero(source, target, dims)` | Zero map |
| `apply(v)` | Apply to a vector |
| `compose(other)` | Compose two maps (matrix multiply) |

### `Protocol` / `ProtocolSheaf` — Protocol as Sheaf

| Method | Description |
|--------|-------------|
| `Protocol::new(name)` | Create empty protocol |
| `add_agent(agent)` | Add an agent |
| `add_rule(rule)` | Add a transition rule |
| `to_sheaf() → ProtocolSheaf` | Build the language sheaf L(P) |
| `ProtocolSheaf::check_deadlock() → DeadlockResult` | H¹ deadlock detection |
| `ProtocolSheaf::compute_cohomology() → Cohomology` | Full cohomology computation |

### `DeadlockResult` — Deadlock Detection Output

```rust
pub struct DeadlockResult {
    pub has_deadlock: bool,
    pub obstruction_dimension: usize,
    pub obstructions: Vec<String>,
}
```

### `Cohomology` — Sheaf Cohomology Groups

| Method | Description |
|--------|-------------|
| `h0_dimension() → usize` | dim H⁰ (global sections) |
| `h1_dimension() → usize` | dim H¹ (obstructions) |
| `is_deadlock_free() → bool` | H¹ = 0 |
| `euler_characteristic() → i64` | χ = dim H⁰ − dim H¹ |
| `betti_numbers() → (usize, usize, usize)` | (h⁰, h¹, h²) |
| `obstruction_classes() → Vec<Vec<f64>>` | Basis for H¹ |
| `trivial()` | Trivial cohomology (point space) |

### `CupProduct` — Protocol Composition

| Method | Description |
|--------|-------------|
| `compute(h1_a, h1_b) → Vec<Vec<f64>>` | Cup product H¹ × H¹ → H² |
| `is_associative() → bool` | Verify (α ⌣ β) ⌣ γ = α ⌣ (β ⌣ γ) |
| `is graded_commutative() → bool` | Verify α ⌣ β = (−1)^(pq) β ⌣ α |

### `Refinement` — Protocol Refinement (Natural Transformation)

| Method | Description |
|--------|-------------|
| `new(source, target, maps)` | Create refinement |
| `verify_naturality() → bool` | Check naturality squares commute |
| `induced_map_h0() → Vec<f64>` | Map on H⁰ (global sections) |
| `induced_map_h1() → Vec<f64>` | Map on H¹ (obstructions) |
| `is_safe_refinement() → bool` | Refinement preserves deadlock-freedom |

### `SheafLaplacian` — Spectral Protocol Verification

| Method | Description |
|--------|-------------|
| `from_sheaf_data(n, dims, edges, weights)` | Build from sheaf |
| `kernel() → Vec<Vec<f64>>` | Global sections (ker LΣ) |
| `kernel_dimension() → usize` | dim ker LΣ |
| `apply(v) → Vec<f64>` | LΣ · v |
| `eigenvalues() → Vec<f64>` | Eigenvalue spectrum |
| `spectral_gap() → f64` | Smallest non-zero eigenvalue |
| `is_global_section(v) → bool` | Is v in ker LΣ? |

### `ErgodicConsensus` — Consensus via Ergodic Theory

| Method | Description |
|--------|-------------|
| `new(n_states, tolerance)` | Create consensus checker |
| `check_consensus(trajectory) → ConsensusResult` | Check if trajectory converges |
| `mixing_time(trajectory) → usize` | Steps to reach consensus |

### Prebuilt Protocol Generators

| Function | Description |
|----------|-------------|
| `mutex_protocol(n_agents, n_resources)` | N-agent mutual exclusion |
| `deadlock_free_ring(n)` | Ring topology (always safe) |
| `circular_wait_deadlock(n)` | N-node circular wait (always deadlocked) |

---

## How It Works

### Architecture

```
Protocol (agents + rules)
  └─ to_sheaf() → ProtocolSheaf (language sheaf L(P))
       ├─ compute_cohomology() → Cohomology (H⁰, H¹)
       │    └─ check_deadlock() → DeadlockResult
       ├─ CupProduct (protocol composition ⌣: H¹ × H¹ → H²)
       ├─ Refinement (natural transformation between sheaf functors)
       ├─ SheafLaplacian (spectral verification)
       └─ ErgodicConsensus (trajectory-based consensus)

Foundation:
  Sheaf (stalks + open sets + restriction maps)
    └─ CochainComplex (d⁰, d¹ coboundary maps)
         └─ Cohomology (ker/im computation)
```

### Module Map

| Module | Contents |
|--------|----------|
| `sheaf` | `Sheaf`, `Stalk`, `OpenSet`, `RestrictionMap`, `SheafSpace` |
| `automata` | `Dfa`, `Nfa` — automaton structures |
| `cohomology` | `Cohomology`, `CochainComplex` — H⁰, H¹ computation |
| `protocol` | `Protocol`, `AgentConfig`, `ProtocolSheaf`, prebuilt protocols |
| `cup_product` | `CupProduct` — protocol composition in cohomology |
| `refinement` | `Refinement` — natural transformations between protocol sheaves |
| `consensus` | `ErgodicConsensus`, `ConsensusResult` — convergence detection |
| `laplacian` | `SheafLaplacian` — spectral analysis of protocol sheaves |

---

## The Math

### Sheaves

A **sheaf** F on a topological space X assigns algebraic data to regions of X in a way that is locally compatible:

1. **Stalks**: For each point x, a stalk F_x (e.g., a vector space)
2. **Restriction maps**: For V ⊆ U, a map ρ_{V→U}: F(U) → F(V)
3. **Gluing axiom**: Compatible local sections can be glued into a global section
4. **Identity axiom**: A section determined locally is determined uniquely

### Čech Cohomology

Given an open cover {Uᵢ} of X and a sheaf F, the **Čech cochain complex** is:

> C⁰ = ⊕ F(Uᵢ) — local sections  
> C¹ = ⊕ F(Uᵢ ∩ Uⱼ) — pairwise sections  
> C² = ⊕ F(Uᵢ ∩ Uⱼ ∩ Uₖ) — triple sections

With coboundary maps:
- d⁰: C⁰ → C¹ sends (sᵢ) to (ρⱼ(sᵢ) − ρᵢ(sⱼ)) on each intersection
- d¹: C¹ → C² defined analogously

The cohomology groups are:
- **H⁰(X, F) = ker(d⁰)** — global sections (consistent configurations)
- **H¹(X, F) = ker(d¹)/im(d⁰)** — obstruction classes (deadlock witnesses)

### Kimi's Theorem 2

> A protocol P is deadlock-free if and only if H¹(Sh(States); L(P)) = 0.

Intuition: H¹ measures the failure of local compatibility to produce global compatibility. A non-zero class in H¹ means there are local configurations that are pairwise compatible but have no global extension — exactly the situation in deadlock.

### Cup Product

The **cup product** ⌣: H^p × H^q → H^(p+q) composes cohomology classes. For protocols, this corresponds to **running two protocols simultaneously**:

- If P₁ has obstruction class α ∈ H¹ and P₂ has β ∈ H¹, then P₁ ∘ P₂ has class α ⌣ β ∈ H²
- The composed protocol can have deadlocks even when both sub-protocols are individually safe

### The Sheaf Laplacian

The **sheaf Laplacian** LΣ = (δ⁰)† δ⁰ is a positive semi-definite operator whose kernel is exactly H⁰ (global sections). Its spectral gap measures how quickly the protocol converges to a consistent state.

---

## Testing

```bash
cargo test
```

**75 tests** covering:

- Sheaf construction (stalks, open sets, restriction maps)
- Stalk operations and dimension calculations
- Open set intersection and union
- Restriction map composition and application
- Cochain complex construction and verification (d¹ ∘ d⁰ = 0)
- Cohomology computation (H⁰, H¹ dimensions, Euler characteristic)
- Kernel and image computation for coboundary maps
- Protocol construction and sheaf generation
- Deadlock detection (free ring vs. circular wait)
- Cup product computation
- Refinement naturality verification
- Sheaf Laplacian (kernel, spectral gap, global sections)
- Consensus detection via ergodic theory

---

## License

MIT
