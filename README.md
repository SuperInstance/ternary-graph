# ternary-graph

Graph algorithms operating on ternary-weighted edges (`-1`, `0`, `+1`).

---

## Why This Exists

Standard graph libraries treat edges as binary (present/absent) or weighted with positive reals. That works if your graph represents a road network or a social network where "friend" is the only relationship that matters.

But the world has more texture than that.

Not all relationships are positive. Social networks have allies **and** adversaries. Gene interactions can be activating or inhibiting. Financial correlations can be positive, negative, or absent. A standard graph library has no way to express that an edge *repels* rather than attracts.

Ternary-weighted graphs solve this: (+1) edges connect, (-1) edges separate, and (0) edges mean no relationship. Shortest paths handle negative-weight edges via Bellman-Ford. Community detection uses signed modularity that rewards positive within-group edges and penalizes negative ones. The signed Laplacian captures both attraction and repulsion in a single matrix.

This crate exists because binary graphs are a lie. Relationships have valence. This library models that.

---

> ⛏️ **DEEP CUT: Why Not Just Use Signed Weights?**  
> 
> You might ask: why restrict to {-1, 0, +1} when you could use the full real number line? The answer is: because continuous weights carry more information than is meaningful.  
> 
> A weight of -0.37 vs -0.42 doesn't capture anything useful about an adversarial relationship that -1 doesn't cover. But it *does* add floating-point drift, numerical instability, and algorithmic complexity to every graph operation.  
> 
> Ternary weights force you to ask: is this edge positive, negative, or absent? There's no room for spurious precision. The three values correspond to three distinct physical states — cooperation, neutral, competition — and those are the only states an agent needs to distinguish. Continuous weights would be a liability, not a feature.  
> 
> The real insight: restricting the value domain *increases* algorithmic robustness. Floyd-Warshall on ternary matrices doesn't have convergence issues. Signed modularity doesn't need threshold tuning. The ternary constraint is the feature, not the limitation.

---

## Core Concepts

- **`Ternary`** — Edge weight: `Neg` (−1, adversarial/inhibiting), `Zero` (0, no edge), `Pos` (+1, friendly/activating).
- **`TernaryGraph`** — Adjacency matrix representation with ternary weights. Supports both directed and undirected graphs.
- **Signed Laplacian** — `L = D − A`, where `D` uses absolute degrees. Captures the structure of positive and negative edges simultaneously.
- **Signed modularity** — Extends Newman-Girvan modularity to signed graphs: positive within-community edges increase modularity, negative ones decrease it.

## Quick Start

```toml
# Cargo.toml
[dependencies]
ternary-graph = "0.1"
```

```rust
use ternary_graph::*;

fn main() {
    // Build a signed social network
    let mut g = TernaryGraph::new(5, false);
    g.add_edge(0, 1, Ternary::Pos); // allies
    g.add_edge(1, 2, Ternary::Pos);
    g.add_edge(2, 3, Ternary::Neg); // adversaries
    g.add_edge(3, 4, Ternary::Pos);
    g.add_edge(0, 4, Ternary::Pos);

    // Shortest paths (handles negative weights via Bellman-Ford)
    let dist = shortest_paths(&g, 0);
    println!("Distances from node 0: {:?}", dist);

    // Community detection
    let communities = label_propagation(&g, 100);
    println!("Communities: {:?}", communities);

    // Spectral clustering into 2 groups
    let clusters = spectral_clustering(&g, 2);
    println!("Spectral clusters: {:?}", clusters);

    // Signed modularity of a partition
    let q = modularity(&g, &clusters);
    println!("Modularity: {:.4}", q);
}
```

## API Overview

### Graph Construction
- `TernaryGraph::new(n, directed)` — Create an `n`-vertex graph
- `g.add_edge(u, v, weight)` — Add a ternary-weighted edge
- `g.neighbors(v)` — Get neighbors with their edge weights
- `g.edge_count()`, `g.degree(v)` — Basic graph statistics

### Graph Matrices
- `g.laplacian()` — Standard Laplacian `L = D − A`
- `g.normalized_laplacian()` — Normalized Laplacian `D^{−1/2} L D^{−1/2}`
- `g.adjacency_f64()`, `g.degree_matrix()` — Raw matrix access

### Shortest Paths
- `shortest_paths(graph, source)` — Single-source via Bellman-Ford. Detects negative cycles and marks affected vertices as unreachable.
- `all_pairs_shortest_paths(graph)` — All-pairs via Floyd-Warshall.

### Community Detection
- `label_propagation(graph, max_iters)` — Weighted label propagation: positive edges attract, negative edges repel.
- `spectral_clustering(graph, k)` — Spectral partitioning using the signed Laplacian's Fiedler vector.
- `modularity(graph, communities)` — Signed modularity score for a given community assignment.
- `connected_components(graph)` — Components of the positive-weight subgraph.

## How It Works

**Bellman-Ford** relaxes all edges `n−1` times, then runs one additional pass to detect negative-weight cycles reachable from the source. Reachable vertices are propagated and marked as having undefined distance.

**Label propagation** initializes each vertex with a unique label, then iteratively updates each vertex's label to the most common label among its neighbors — weighted by edge sign (positive edges attract, negative edges repel). Convergence is typically fast (10-50 iterations for well-structured graphs).

**Spectral clustering** computes the signed Laplacian, sorts eigenvectors by eigenvalue, and uses the k smallest-eigenvalue eigenvectors to embed vertices. The smallest eigenvalue gives the Fiedler vector directly.

**Signed modularity** computes the difference between actual edge density within communities and expected edge density under a null model, where positive edges within communities increase the score and negative edges within communities penalize it.

License: MIT
