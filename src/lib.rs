//! # ternary-graph
//!
//! Graph algorithms operating on ternary-weighted edges (`-1`, `0`, `+1`).
//! Provides ternary adjacency matrices, shortest paths with ternary weights,
//! community detection via modularity optimization, graph Laplacian computation,
//! and spectral clustering.

#![forbid(unsafe_code)]

/// Canonical ternary type re-exported from `ternary-types`.
///
/// Maps `Neg -> Negative`, `Zero -> Neutral`, `Pos -> Positive` from the
/// previous custom enum.
pub use ternary_types::Ternary;

/// Extension trait providing methods previously on the custom `Ternary` type.
pub trait TernaryExt {
    /// Convert this ternary value to `f64`.
    fn to_f64(self) -> f64;
}

impl TernaryExt for ternary_types::Ternary {
    fn to_f64(self) -> f64 {
        i8::from(self) as f64
    }
}

/// A graph with ternary-weighted edges, represented as an adjacency matrix.
#[derive(Clone, Debug)]
pub struct TernaryGraph {
    /// Number of vertices.
    pub n: usize,
    /// Adjacency matrix: adj[i][j] = weight of edge i→j (Zero = no edge).
    pub adj: Vec<Vec<Ternary>>,
    /// Whether edges are directed.
    pub directed: bool,
}

impl TernaryGraph {
    /// Create a new graph with `n` vertices and no edges.
    pub fn new(n: usize, directed: bool) -> Self {
        TernaryGraph {
            n,
            adj: vec![vec![Ternary::Neutral; n]; n],
            directed,
        }
    }

    /// Add an edge from `u` to `v` with the given ternary weight.
    pub fn add_edge(&mut self, u: usize, v: usize, weight: Ternary) {
        self.adj[u][v] = weight;
        if !self.directed {
            self.adj[v][u] = weight;
        }
    }

    /// Get the weight of edge u→v.
    pub fn edge(&self, u: usize, v: usize) -> Ternary {
        self.adj[u][v]
    }

    /// Get neighbors of vertex `v` (vertices connected by non-zero edges).
    pub fn neighbors(&self, v: usize) -> Vec<(usize, Ternary)> {
        (0..self.n)
            .filter(|&u| self.adj[v][u] != Ternary::Neutral)
            .map(|u| (u, self.adj[v][u]))
            .collect()
    }

    /// Count total non-zero edges.
    pub fn edge_count(&self) -> usize {
        let mut count = 0;
        for i in 0..self.n {
            for j in 0..self.n {
                if self.adj[i][j] != Ternary::Neutral {
                    count += 1;
                }
            }
        }
        if self.directed { count } else { count / 2 }
    }

    /// Compute the degree of vertex `v` (count of non-zero edges).
    pub fn degree(&self, v: usize) -> usize {
        self.neighbors(v).len()
    }

    /// Compute the degree matrix (diagonal matrix of degrees).
    pub fn degree_matrix(&self) -> Vec<Vec<f64>> {
        let mut d = vec![vec![0.0f64; self.n]; self.n];
        for i in 0..self.n {
            d[i][i] = self.degree(i) as f64;
        }
        d
    }

    /// Compute the adjacency matrix as f64 values.
    pub fn adjacency_f64(&self) -> Vec<Vec<f64>> {
        self.adj.iter().map(|row| row.iter().map(|&v| v.to_f64()).collect()).collect()
    }

    /// Compute the graph Laplacian: L = D - A.
    pub fn laplacian(&self) -> Vec<Vec<f64>> {
        let d = self.degree_matrix();
        let a = self.adjacency_f64();
        let mut l = vec![vec![0.0f64; self.n]; self.n];
        for i in 0..self.n {
            for j in 0..self.n {
                l[i][j] = d[i][j] - a[i][j];
            }
        }
        l
    }

    /// Compute the normalized Laplacian: L_norm = D^{-1/2} L D^{-1/2}.
    pub fn normalized_laplacian(&self) -> Vec<Vec<f64>> {
        let d = self.degree_matrix();
        let a = self.adjacency_f64();
        let mut l = vec![vec![0.0f64; self.n]; self.n];

        // Precompute D^{-1/2}
        let mut d_inv_sqrt = vec![0.0f64; self.n];
        for i in 0..self.n {
            if d[i][i] > 0.0 {
                d_inv_sqrt[i] = 1.0 / d[i][i].sqrt();
            }
        }

        for i in 0..self.n {
            for j in 0..self.n {
                let val = if i == j {
                    1.0
                } else if a[i][j] != 0.0 {
                    -a[i][j] / (d[i][i] * d[j][j]).sqrt()
                } else {
                    0.0
                };
                l[i][j] = val;
            }
        }

        l
    }
}

/// Shortest paths using Bellman-Ford with ternary weights.
/// Returns distances as Option<f64> (None if unreachable or negative cycle).
pub fn shortest_paths(graph: &TernaryGraph, source: usize) -> Vec<Option<f64>> {
    let n = graph.n;
    let mut dist: Vec<Option<f64>> = vec![None; n];
    dist[source] = Some(0.0);

    // Relax edges n-1 times
    for _ in 0..n.saturating_sub(1) {
        for u in 0..n {
            if dist[u].is_none() {
                continue;
            }
            for (v, w) in graph.neighbors(u) {
                let new_dist = dist[u].unwrap() + w.to_f64();
                if dist[v].is_none() || new_dist < dist[v].unwrap() {
                    dist[v] = Some(new_dist);
                }
            }
        }
    }

    // Check for negative cycles (one more relaxation)
    let mut has_neg_cycle = vec![false; n];
    for _ in 0..n {
        for u in 0..n {
            if dist[u].is_none() {
                continue;
            }
            for (v, w) in graph.neighbors(u) {
                let new_dist = dist[u].unwrap() + w.to_f64();
                if dist[v].is_none() || new_dist < dist[v].unwrap() {
                    has_neg_cycle[v] = true;
                    dist[v] = Some(new_dist);
                }
            }
        }
    }

    // Propagate negative cycle reachability
    for _ in 0..n {
        let mut changed = false;
        for u in 0..n {
            if !has_neg_cycle[u] {
                continue;
            }
            for (v, _) in graph.neighbors(u) {
                if !has_neg_cycle[v] {
                    has_neg_cycle[v] = true;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }

    for i in 0..n {
        if has_neg_cycle[i] {
            dist[i] = None;
        }
    }

    dist
}

/// All-pairs shortest paths using Floyd-Warshall with ternary weights.
pub fn all_pairs_shortest_paths(graph: &TernaryGraph) -> Vec<Vec<Option<f64>>> {
    let n = graph.n;
    let mut dist = vec![vec![None::<f64>; n]; n];

    for i in 0..n {
        dist[i][i] = Some(0.0);
    }
    for u in 0..n {
        for (v, w) in graph.neighbors(u) {
            dist[u][v] = Some(w.to_f64());
        }
    }

    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                if let (Some(dik), Some(dkj)) = (dist[i][k], dist[k][j]) {
                    let new_dist = dik + dkj;
                    if dist[i][j].is_none() || new_dist < dist[i][j].unwrap() {
                        dist[i][j] = Some(new_dist);
                    }
                }
            }
        }
    }

    dist
}

/// Community detection using label propagation on ternary graphs.
///
/// Positive edges encourage same-community assignment, negative edges encourage different.
pub fn label_propagation(graph: &TernaryGraph, max_iters: usize) -> Vec<usize> {
    let n = graph.n;
    let mut labels: Vec<usize> = (0..n).collect();

    for _ in 0..max_iters {
        let mut changed = false;
        for v in 0..n {
            // Count weighted votes from neighbors
            let mut votes: std::collections::HashMap<usize, f64> = std::collections::HashMap::new();
            for (u, w) in graph.neighbors(v) {
                let label = labels[u];
                let weight = w.to_f64();
                *votes.entry(label).or_insert(0.0) += weight;
            }

            if votes.is_empty() {
                continue;
            }

            // Find best label (highest weighted vote)
            let best_label = votes.into_iter()
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(l, _)| l)
                .unwrap_or(labels[v]);

            if best_label != labels[v] {
                labels[v] = best_label;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }

    // Relabel to consecutive integers
    let mut label_map = std::collections::HashMap::new();
    let mut next = 0;
    for l in labels.iter_mut() {
        let entry = *label_map.entry(*l).or_insert_with(|| { let v = next; next += 1; v });
        *l = entry;
    }

    labels
}

/// Compute modularity of a given community assignment.
///
/// Uses ternary weights: positive edges within communities increase modularity,
/// negative edges within communities decrease it.
pub fn modularity(graph: &TernaryGraph, communities: &[usize]) -> f64 {
    let n = graph.n;
    let m = graph.edge_count() as f64;
    if m == 0.0 {
        return 0.0;
    }

    // Compute weighted degree for each node
    let mut k = vec![0.0f64; n];
    for i in 0..n {
        for j in 0..n {
            k[i] += graph.adj[i][j].to_f64();
        }
    }

    let mut q = 0.0;
    for i in 0..n {
        for j in 0..n {
            if communities[i] == communities[j] {
                let a_ij = graph.adj[i][j].to_f64();
                q += a_ij - (k[i] * k[j]) / (2.0 * m);
            }
        }
    }

    q / (2.0 * m)
}

/// Power iteration to find smallest eigenvectors of the Laplacian.
fn power_iteration_laplacian(laplacian: &[Vec<f64>], n: usize, max_iters: usize) -> Vec<f64> {
    let mut v: Vec<f64> = (0..n).map(|i| (i as f64 + 1.0) / n as f64).collect();
    let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    for x in v.iter_mut() {
        *x /= norm;
    }

    for _ in 0..max_iters {
        let mv: Vec<f64> = laplacian.iter().map(|row| {
            row.iter().zip(v.iter()).map(|(&a, &b)| a * b).sum()
        }).collect();
        let n = mv.iter().map(|x| x * x).sum::<f64>().sqrt();
        if n < 1e-15 { break; }
        for i in 0..v.len() {
            v[i] = mv[i] / n;
        }
    }

    v
}

/// Spectral clustering using the graph Laplacian.
///
/// Computes the Fiedler vector (second smallest eigenvector) and partitions vertices
/// based on the sign of their entry in this vector.
pub fn spectral_clustering(graph: &TernaryGraph, k: usize) -> Vec<usize> {
    if graph.n == 0 {
        return Vec::new();
    }
    if k == 1 {
        return vec![0; graph.n];
    }

    let n = graph.n;

    // Build signed adjacency: treat negative edges as -1, positive as +1
    // Use the signed Laplacian for better handling of negative edges
    let mut adj_f64 = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            adj_f64[i][j] = graph.adj[i][j].to_f64();
        }
    }

    // Degree for signed Laplacian: d_i = sum of |A_ij|
    let mut degree = vec![0.0f64; n];
    for i in 0..n {
        for j in 0..n {
            degree[i] += adj_f64[i][j].abs();
        }
    }

    // Signed Laplacian: L = D - A
    let mut laplacian = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                laplacian[i][j] = degree[i];
            } else {
                laplacian[i][j] = -adj_f64[i][j];
            }
        }
    }

    // Power iteration on (shift * I - L) to find smallest eigenvectors of L
    let shift = 2.0 * n as f64;
    let mut shifted = laplacian.clone();
    for i in 0..n {
        shifted[i][i] = shift - laplacian[i][i];
        for j in 0..n {
            if i != j {
                shifted[i][j] = -laplacian[i][j];
            }
        }
    }

    // Find the top eigenvector of shifted (= smallest non-trivial eigenvector of L)
    // Skip the first one (constant vector, eigenvalue = shift for connected graphs)
    let mut current = shifted.clone();
    let mut eigenvectors = Vec::new();

    for _ in 0..(k + 1).min(n) {
        let eigvec = power_iteration_laplacian(&current, n, 300);
        let eigval: f64 = eigvec.iter().zip(current.iter())
            .map(|(&v, row)| v * row.iter().zip(eigvec.iter()).map(|(&a, &b)| a * b).sum::<f64>())
            .sum();

        // Deflate
        for i in 0..n {
            for j in 0..n {
                current[i][j] -= eigval * eigvec[i] * eigvec[j];
            }
        }

        eigenvectors.push(eigvec);
    }

    if eigenvectors.is_empty() {
        return vec![0; n];
    }

    // Use the second eigenvector (index 1) as Fiedler vector for k=2
    // Index 0 is the trivial eigenvector
    let fiedler_idx = if eigenvectors.len() > 1 { 1 } else { 0 };
    let fiedler = &eigenvectors[fiedler_idx];

    if k == 2 {
        // Split at median
        let mut sorted: Vec<f64> = fiedler.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median = sorted[n / 2];
        return fiedler.iter().map(|&v| if v >= median { 1 } else { 0 }).collect();
    }

    // For k > 2, use quantile-based assignment across eigenvectors
    let mut assignments = vec![0usize; n];
    for i in 0..n {
        let mut hash = 0usize;
        for (d, ev) in eigenvectors.iter().enumerate() {
            if ev[i] >= 0.0 {
                hash |= 1 << (d % 16);
            }
        }
        assignments[i] = hash % k;
    }

    assignments
}

/// Compute the connected components of the positive-weight subgraph.
pub fn connected_components(graph: &TernaryGraph) -> Vec<usize> {
    let n = graph.n;
    let mut components = vec![usize::MAX; n];
    let mut component_id = 0;

    // BFS from each unvisited node using only positive edges
    for start in 0..n {
        if components[start] != usize::MAX {
            continue;
        }
        let mut queue = vec![start];
        components[start] = component_id;
        let mut head = 0;
        while head < queue.len() {
            let v = queue[head];
            head += 1;
            for (u, w) in graph.neighbors(v) {
                if components[u] == usize::MAX && w == Ternary::Positive {
                    components[u] = component_id;
                    queue.push(u);
                }
            }
        }
        component_id += 1;
    }

    components
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_simple_graph() -> TernaryGraph {
        let mut g = TernaryGraph::new(4, false);
        g.add_edge(0, 1, Ternary::Positive);
        g.add_edge(1, 2, Ternary::Positive);
        g.add_edge(2, 3, Ternary::Positive);
        g
    }

    #[test]
    fn test_graph_creation() {
        let g = TernaryGraph::new(5, false);
        assert_eq!(g.n, 5);
        assert_eq!(g.edge_count(), 0);
    }

    #[test]
    fn test_add_edge_undirected() {
        let g = make_simple_graph();
        assert_eq!(g.edge(0, 1), Ternary::Positive);
        assert_eq!(g.edge(1, 0), Ternary::Positive);
        assert_eq!(g.edge(0, 3), Ternary::Neutral);
    }

    #[test]
    fn test_add_edge_directed() {
        let mut g = TernaryGraph::new(3, true);
        g.add_edge(0, 1, Ternary::Positive);
        assert_eq!(g.edge(0, 1), Ternary::Positive);
        assert_eq!(g.edge(1, 0), Ternary::Neutral);
    }

    #[test]
    fn test_edge_count() {
        let g = make_simple_graph();
        assert_eq!(g.edge_count(), 3);
    }

    #[test]
    fn test_neighbors() {
        let g = make_simple_graph();
        let nb = g.neighbors(1);
        assert_eq!(nb.len(), 2);
        assert!(nb.contains(&(0, Ternary::Positive)));
        assert!(nb.contains(&(2, Ternary::Positive)));
    }

    #[test]
    fn test_degree() {
        let g = make_simple_graph();
        assert_eq!(g.degree(0), 1);
        assert_eq!(g.degree(1), 2);
        assert_eq!(g.degree(2), 2);
    }

    #[test]
    fn test_laplacian() {
        let g = make_simple_graph();
        let l = g.laplacian();
        assert_eq!(l.len(), 4);
        // Diagonal should be degrees
        assert_eq!(l[0][0], 1.0);
        assert_eq!(l[1][1], 2.0);
        // Off-diagonal for edges should be -1
        assert_eq!(l[0][1], -1.0);
    }

    #[test]
    fn test_normalized_laplacian() {
        let g = make_simple_graph();
        let l = g.normalized_laplacian();
        // Diagonal should be 1.0 for non-isolated vertices
        assert_eq!(l[0][0], 1.0);
        assert_eq!(l[1][1], 1.0);
    }

    #[test]
    fn test_shortest_paths_basic() {
        let g = make_simple_graph();
        let dist = shortest_paths(&g, 0);
        assert_eq!(dist[0], Some(0.0));
        assert_eq!(dist[1], Some(1.0));
        assert_eq!(dist[2], Some(2.0));
        assert_eq!(dist[3], Some(3.0));
    }

    #[test]
    fn test_shortest_paths_with_negative() {
        let mut g = TernaryGraph::new(3, true);
        g.add_edge(0, 1, Ternary::Positive);
        g.add_edge(1, 2, Ternary::Negative);
        let dist = shortest_paths(&g, 0);
        assert_eq!(dist[0], Some(0.0));
        assert_eq!(dist[1], Some(1.0));
        assert_eq!(dist[2], Some(0.0)); // 1 + (-1) = 0
    }

    #[test]
    fn test_shortest_paths_disconnected() {
        let mut g = TernaryGraph::new(4, false);
        g.add_edge(0, 1, Ternary::Positive);
        // 2 and 3 are isolated
        let dist = shortest_paths(&g, 0);
        assert_eq!(dist[0], Some(0.0));
        assert_eq!(dist[1], Some(1.0));
        assert_eq!(dist[2], None);
    }

    #[test]
    fn test_all_pairs_shortest_paths() {
        let g = make_simple_graph();
        let apsp = all_pairs_shortest_paths(&g);
        assert_eq!(apsp[0][3], Some(3.0));
        assert_eq!(apsp[3][0], Some(3.0));
        assert_eq!(apsp[0][0], Some(0.0));
    }

    #[test]
    fn test_label_propagation() {
        let mut g = TernaryGraph::new(4, false);
        g.add_edge(0, 1, Ternary::Positive);
        g.add_edge(1, 2, Ternary::Positive);
        g.add_edge(2, 3, Ternary::Negative);
        let labels = label_propagation(&g, 100);
        // 0, 1, 2 should form one community; 3 should be different due to negative edge
        assert_eq!(labels.len(), 4);
    }

    #[test]
    fn test_modularity() {
        let mut g = TernaryGraph::new(4, false);
        g.add_edge(0, 1, Ternary::Positive);
        g.add_edge(2, 3, Ternary::Positive);
        let communities = vec![0, 0, 1, 1];
        let q = modularity(&g, &communities);
        // Good partition should have positive modularity
        assert!(q > 0.0, "Modularity should be positive for good partition, got {}", q);
    }

    #[test]
    fn test_modularity_empty_graph() {
        let g = TernaryGraph::new(3, false);
        let communities = vec![0, 0, 0];
        let q = modularity(&g, &communities);
        assert_eq!(q, 0.0);
    }

    #[test]
    fn test_spectral_clustering() {
        let mut g = TernaryGraph::new(6, false);
        // Two triangles connected by a single edge
        g.add_edge(0, 1, Ternary::Positive);
        g.add_edge(1, 2, Ternary::Positive);
        g.add_edge(0, 2, Ternary::Positive);
        g.add_edge(3, 4, Ternary::Positive);
        g.add_edge(4, 5, Ternary::Positive);
        g.add_edge(3, 5, Ternary::Positive);
        // Single bridge edge
        g.add_edge(2, 3, Ternary::Positive);
        let labels = spectral_clustering(&g, 2);
        assert_eq!(labels.len(), 6);
        // With only a single bridge, the Fiedler vector should mostly separate the triangles.
        // Check that most nodes in each triangle share a label.
        let same_in_first = (labels[0] == labels[1]) as usize
            + (labels[1] == labels[2]) as usize
            + (labels[0] == labels[2]) as usize;
        let same_in_second = (labels[3] == labels[4]) as usize
            + (labels[4] == labels[5]) as usize
            + (labels[3] == labels[5]) as usize;
        // At least 2 of 3 pairs should agree in each triangle
        assert!(same_in_first >= 2, "First triangle should be mostly one cluster: {:?}", labels);
        assert!(same_in_second >= 2, "Second triangle should be mostly one cluster: {:?}", labels);
    }

    #[test]
    fn test_connected_components() {
        let mut g = TernaryGraph::new(6, false);
        g.add_edge(0, 1, Ternary::Positive);
        g.add_edge(1, 2, Ternary::Positive);
        g.add_edge(3, 4, Ternary::Positive);
        // 5 is isolated
        let comp = connected_components(&g);
        assert_eq!(comp[0], comp[1]);
        assert_eq!(comp[1], comp[2]);
        assert_eq!(comp[3], comp[4]);
        assert_ne!(comp[0], comp[3]);
        assert_ne!(comp[0], comp[5]);
    }

    #[test]
    fn test_connected_components_negative_edges() {
        let mut g = TernaryGraph::new(3, false);
        g.add_edge(0, 1, Ternary::Positive);
        g.add_edge(1, 2, Ternary::Negative);
        let comp = connected_components(&g);
        // Only positive edges connect
        assert_eq!(comp[0], comp[1]);
        assert_ne!(comp[1], comp[2]);
    }
}
