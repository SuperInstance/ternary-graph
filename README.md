# ternary-graph

Graph algorithms operating on ternary-weighted edges {-1, 0, +1}. Provides adjacency matrices, shortest paths with signed weights, graph Laplacian computation, community detection via modularity optimization, and spectral clustering — all in the ternary weight domain.

## Why It Matters

Graph algorithms typically assume non-negative edge weights (Dijkstra requires it). Real-world graphs often have **signed relationships**: trust/distrust networks (Leskovec et al., 2010), excitatory/inhibitory neural connections, and ternary-weighted neural networks where weights are quantized to {-1, 0, +1}.

This crate provides:

- **Ternary adjacency matrices** — O(1) edge lookup with 2-bit storage per edge
- **Graph Laplacian** — L = D − A, the foundation of spectral clustering
- **Normalized Laplacian** — L_norm = D^(−1/2) L D^(−1/2), used in normalized spectral clustering (Ng, Jordan & Weiss, 2001)
- **Community detection** via modularity maximization
- **Spectral clustering** using Laplacian eigenvectors

## How It Works

### Ternary Edge Weights

Each edge has weight w ∈ {-1, 0, +1}:

| Weight | Semantic | Meaning |
|--------|----------|---------|
| +1 | Positive edge | Excitatory / trust / similarity |
| 0 | No edge | Disconnected |
| −1 | Negative edge | Inhibitory / distrust / dissimilarity |

The adjacency matrix A is a ternary matrix: A[i][j] ∈ {-1, 0, +1}. Zero means no edge exists, not a "zero-weight edge."

### Graph Laplacian

The combinatorial Laplacian is:

```
L = D − A
```

where D is the degree matrix (diagonal, D_ii = degree(i)) and A is the adjacency matrix.

For vertex i:

```
L_ii = degree(i) = |{j : A[i][j] ≠ 0}|
L_ij = −A[i][j]   if i ≠ j and edge exists
L_ij = 0          otherwise
```

**Key property**: L is positive semi-definite. Its smallest eigenvalue is λ₁ = 0 with eigenvector the all-ones vector (for connected graphs). The second-smallest eigenvalue λ₂ (the **algebraic connectivity** or Fiedler value) measures how well-connected the graph is.

### Normalized Laplacian

```
L_norm = D^(−1/2) L D^(−1/2)
```

For vertex i with degree d_i:

```
L_norm[i][i] = 1
L_norm[i][j] = −A[i][j] / √(d_i · d_j)    if edge exists
```

Eigenvalues of L_norm lie in [0, 2]. The spectral gap (λ₂ of L_norm) determines the number of clusters k via the **eigengap heuristic**.

### Spectral Clustering

The implementation in this crate uses the **signed** Laplacian (`D_ii = Σ_j |A_ij|`)
and the standard "shift + power iteration + Hotelling deflation" pipeline:

1. Build the signed Laplacian `L = D − A`.
2. Form `2n·I − L` so that the smallest eigenvalues of `L` become the largest.
3. Use power iteration with deflation to extract the first `k + 1` eigenvectors.
4. Discard the trivial constant eigenvector; for `k = 2` split vertices at
   the median of the Fiedler vector; for `k > 2` hash each vertex's
   eigenvector-sign pattern mod `k`.

This is a lightweight stand-in for the full Ng-Jordan-Weiss row-wise k-means
step. It produces correct partitions on graphs with a clean spectral gap.

### Complexity

| Operation | Time | Space |
|-----------|------|-------|
| `new(n)` | O(n²) | O(n²) |
| `add_edge(u, v, w)` | O(1) | O(1) |
| `neighbors(v)` | O(n) | O(n) |
| `edge_count()` | O(n²) | O(1) |
| `laplacian()` | O(n²) | O(n²) |
| `normalized_laplacian()` | O(n²) | O(n²) |
| `spectral_clustering(k)` | O(n²·T + nkI) | O(n² + nk) |

Where n = vertices, T = power iteration steps, I = k-means iterations.

### Signed Shortest Paths

With ternary weights, "shortest" is redefined. The signed distance accumulates:

```
d(s→t) = Σₑ∈path w(e)    where w(e) ∈ {-1, 0, +1}
```

A path with negative total weight is "antagonistic" — it represents a chain of inhibitory connections. The Bellman-Ford algorithm handles negative weights and detects negative cycles.

## Quick Start

```rust
use ternary_graph::{TernaryGraph, Ternary};

// Create a graph with 6 vertices
let mut g = TernaryGraph::new(6, false);

// Add ternary-weighted edges
g.add_edge(0, 1, Ternary::Positive);  // excitatory
g.add_edge(1, 2, Ternary::Positive);
g.add_edge(2, 3, Ternary::Negative);  // inhibitory
g.add_edge(3, 4, Ternary::Positive);
g.add_edge(4, 5, Ternary::Positive);
g.add_edge(5, 0, Ternary::Negative);

// Compute Laplacian
let l = g.laplacian();

// Compute normalized Laplacian
let l_norm = g.normalized_laplacian();

// Get neighbors
let neighbors = g.neighbors(0);
println!("Vertex 0 has {} connections", neighbors.len());
// -> Vertex 0 has 2 connections
```

Add `ternary-graph` to your `Cargo.toml`:

```toml
[dependencies]
ternary-graph = { git = "https://github.com/SuperInstance/ternary-graph.git" }
```

The crate re-exports the [`Ternary`](https://docs.rs/ternary-types/latest/ternary_types/enum.Ternary.html)
enum from `ternary-types`, so `Ternary::Positive` / `Neutral` / `Negative`
are all available directly under `ternary_graph::`.

## API

### `TernaryGraph` methods

| Method | Description |
|--------|-------------|
| `new(n, directed)` | Create `n`-vertex graph with no edges |
| `add_edge(u, v, weight)` | Add (or overwrite) a ternary-weighted edge; `Neutral` removes it |
| `edge(u, v) -> Ternary` | Query edge weight between `u` and `v` |
| `neighbors(v) -> Vec<(usize, Ternary)>` | List `(neighbor, weight)` pairs in ascending index order |
| `degree(v) -> usize` | Count non-`Neutral` edges incident to `v` |
| `edge_count() -> usize` | Total edges (off-diagonal counted once for undirected; self-loops once) |
| `adjacency_f64() -> Vec<Vec<f64>>` | Adjacency matrix widened to `f64` |
| `degree_matrix() -> Vec<Vec<f64>>` | Diagonal degree matrix `D` |
| `laplacian() -> Vec<Vec<f64>>` | Combinatorial Laplacian `L = D − A` |
| `normalized_laplacian() -> Vec<Vec<f64>>` | Sym. normalized Laplacian `D^(−1/2) L D^(−1/2)` |

### Free functions

| Function | Description |
|----------|-------------|
| `shortest_paths(&graph, src) -> Vec<Option<f64>>` | Bellman-Ford from `src`; returns `None` for unreachable vertices or vertices on / downstream of a negative cycle |
| `all_pairs_shortest_paths(&graph) -> Vec<Vec<Option<f64>>>` | Floyd-Warshall all-pairs; same `None` semantics for negative cycles |
| `connected_components(&graph) -> Vec<usize>` | Components of the *positive-weight* subgraph |
| `label_propagation(&graph, max_iters) -> Vec<usize>` | Deterministic weighted label propagation |
| `modularity(&graph, &communities) -> f64` | Signed Newman modularity `Q` |
| `spectral_clustering(&graph, k) -> Vec<usize>` | Spectral clustering via the signed Laplacian and power iteration |

## Architecture Notes

This crate bridges **both γ and η layers** of the γ + η = C framework:

- **η (eta)**: The graph algorithms — Laplacian computation, clustering, shortest paths — are the computational primitives operating on ternary data.
- **γ (gamma)**: The graph structure encodes the dependency/topology information that the γ layer uses for coordination and routing decisions.
- **C**: The complete graph-based analysis system. The ternary weights {−1, 0, +1} are the same domain used for ternary neural network weights, enabling direct analysis of ternary model architectures as graphs.

## References

- **Spectral Clustering**: Ng, A., Jordan, M. & Weiss, Y., "On Spectral Clustering: Analysis and an Algorithm," NeurIPS 2001.
- **Signed Networks**: Leskovec, J., Huttenlocher, D. & Kleinberg, J., "Signed Networks in Social Media," CHI 2010.
- **Graph Laplacian**: Chung, F., "Spectral Graph Theory," CBMS Regional Conference Series in Mathematics, 92, 1997.
- **Modularity**: Newman, M.E.J., "Modularity and Community Structure in Networks," PNAS, 103(23), 8577-8582, 2006.
- **Algebraic Connectivity**: Fiedler, M., "Algebraic Connectivity of Graphs," Czechoslovak Mathematical Journal, 23(2), 298-305, 1973.
- **Bellman-Ford**: Bellman, R., "On a Routing Problem," Quarterly of Applied Mathematics, 16(1), 87-90, 1958.

## License

MIT
