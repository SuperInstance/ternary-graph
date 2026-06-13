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

1. Compute L_norm
2. Extract the k smallest eigenvectors v₁, ..., vₖ
3. Stack them as columns: U ∈ ℝ^(n×k)
4. Normalize rows of U to unit length
5. Cluster rows of U via k-means

The ternary weight structure makes eigenvector computation particularly efficient — the power iteration becomes a sequence of additions and subtractions (no multiplications for the non-zero entries).

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
let L = g.laplacian();

// Compute normalized Laplacian
let L_norm = g.normalized_laplacian();

// Get neighbors
let neighbors = g.neighbors(0);
println!("Vertex 0 has {} connections", neighbors.len());
```

## API

### `TernaryGraph`

| Method | Description |
|--------|-------------|
| `new(n, directed)` | Create n-vertex graph |
| `add_edge(u, v, weight)` | Add ternary-weighted edge |
| `edge(u, v) -> Ternary` | Query edge weight |
| `neighbors(v) -> Vec<(usize, Ternary)>` | Get (neighbor, weight) pairs |
| `degree(v) -> usize` | Count non-zero edges |
| `laplacian() -> Vec<Vec<f64>>` | Compute L = D − A |
| `normalized_laplacian() -> Vec<Vec<f64>>` | Compute D^(−1/2) L D^(−1/2) |
| `adjacency_f64() -> Vec<Vec<f64>>` | Adjacency as f64 matrix |
| `degree_matrix() -> Vec<Vec<f64>>` | Diagonal degree matrix |
| `edge_count() -> usize` | Total non-zero edges |

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
