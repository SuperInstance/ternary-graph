# ternary-graph

Graph algorithms operating on ternary-weighted edges (`-1`, `0`, `+1`).

## Why This Exists

Not all relationships are positive. Social networks have allies and adversaries. Gene interactions can be activating or inhibiting. Financial correlations can be positive, negative, or absent. Standard graph libraries treat edges as binary (present/absent) or weighted with positive reals — they can't express that an edge *repels* rather than attracts. This crate provides graph algorithms that natively understand ternary edge weights: positive edges connect, negative edges separate, and zero edges mean no relationship. Shortest paths handle negative-weight edges via Bellman-Ford. Community detection uses signed modularity. Spectral clustering uses the signed Laplacian.

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

**Label propagation** iteratively assigns each vertex the label with the highest weighted vote from its neighbors (positive edges vote for the same label, negative edges vote against). Convergence occurs when no labels change.

**Spectral clustering** builds the signed Laplacian (`D` uses absolute-value degrees), shifts it to make all eigenvalues positive, then uses power iteration to find the Fiedler vector (second eigenvector). For `k=2`, vertices are split at the median; for larger `k`, a hash-based assignment across eigenvectors is used.

## Use Cases

1. **Social network analysis** — Detect communities in networks with friendships (positive) and rivalries (negative).
2. **Gene regulatory networks** — Identify functional modules where genes have activating (+1) or inhibiting (−1) interactions.
3. **Financial correlation graphs** — Cluster assets based on positive/negative/uncorrelated return relationships.
4. **Recommendation systems** — Model like/dislike/neutral relationships between users and items as a signed bipartite graph.

## Ecosystem

- [`ternary-clustering`](https://github.com/user/ternary-clustering) — Clustering algorithms for ternary data
- [`ternary-automata`](https://github.com/user/ternary-automata) — Cellular automata with ternary states
- [`ternary-projection`](https://github.com/user/ternary-projection) — Dimensionality reduction for ternary data

## License

MIT

## See Also
- **ternary-network** — related
- **ternary-topology** — related
- **ternary-mesh** — related
- **ternary-clustering** — related
- **ternary-path** — related
- **ternary-petri** — related

