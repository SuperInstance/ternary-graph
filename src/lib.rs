//! # ternary-graph
//!
//! Graph algorithms over ternary-weighted edges (`-1`, `0`, `+1`): adjacency
//! matrices, shortest paths with signed weights (Bellman-Ford / Floyd-Warshall),
//! the combinatorial and normalized graph Laplacians, connected components,
//! label-propagation community detection, modularity scoring, and spectral
//! clustering.
//!
//! ## When to use this crate
//!
//! Use `ternary-graph` when your edges are signed and discrete — i.e. each
//! edge is one of `{+1, 0, -1}` rather than an arbitrary `f64`. This shows up
//! naturally in:
//!
//! - **Signed social networks** (trust `+1` / distrust `-1`, with `0` meaning
//!   "no relationship").
//! - **Ternary neural networks** whose weights are quantized to `{-1, 0, +1}`,
//!   analysed as a computational graph.
//! - **Excitatory / inhibitory neuronal or chemical-reaction graphs**.
//! - **Dependency / conflict graphs** where `+1` is "supports", `-1` is
//!   "blocks", and `0` is "no constraint".
//!
//! Because the edge alphabet is exactly `{+1, 0, -1}`, every algorithm in this
//! crate is built on integer additions and subtractions — no float weights are
//! introduced until they are required (Laplacian eigenvalues, modularity
//! normalization, etc.).
//!
//! ## Example
//!
//! ```
//! use ternary_graph::{TernaryGraph, Ternary};
//!
//! // 4-vertex undirected path  0 -- 1 -- 2 -- 3  (all excitatory).
//! let mut g = TernaryGraph::new(4, false);
//! g.add_edge(0, 1, Ternary::Positive);
//! g.add_edge(1, 2, Ternary::Positive);
//! g.add_edge(2, 3, Ternary::Positive);
//!
//! // Bellman-Ford handles the (potential) negative edges correctly.
//! let dist = ternary_graph::shortest_paths(&g, 0);
//! assert_eq!(dist[3], Some(3.0));
//! ```
//!
//! ## Edge semantics
//!
//! | Weight           | [`Ternary`] variant | Meaning                                |
//! |------------------|---------------------|----------------------------------------|
//! | `+1`             | [`Positive`]        | Excitatory / trust / attracts          |
//! | `0`              | [`Neutral`]         | No edge                                |
//! | `-1`             | [`Negative`]        | Inhibitory / distrust / repels         |
//!
//! [`Ternary`]: ternary_types::Ternary
//! [`Positive`]: ternary_types::Ternary::Positive
//! [`Negative`]: ternary_types::Ternary::Negative
//! [`Neutral`]: ternary_types::Ternary::Neutral

#![forbid(unsafe_code)]

/// Canonical ternary type re-exported from [`ternary_types`].
///
/// This is the same enum used across the SuperInstance ternary fleet:
/// `Negative = -1`, `Neutral = 0`, `Positive = +1`. Re-exporting it here means
/// callers only need a single `use` to bring both the graph type and its edge
/// weight type into scope.
pub use ternary_types::Ternary;

/// A graph whose edge weights are ternary (`-1`, `0`, `+1`), stored as a dense
/// adjacency matrix.
///
/// Storage is `O(n²)` regardless of how many edges are present, which makes this
/// structure ideal for small-to-medium dense signed graphs (e.g. trust matrices,
/// ternary-weighted network layers, dependency/conflict graphs). For very large
/// sparse graphs, prefer a CSR-style representation.
///
/// The adjacency value `Ternary::Neutral` (`0`) means "no edge"; the variants
/// `Ternary::Positive` and `Ternary::Negative` carry the `+1` / `-1` weight.
#[derive(Clone, Debug)]
pub struct TernaryGraph {
    /// Number of vertices. Vertices are indexed `0..n`.
    pub n: usize,
    /// Dense adjacency matrix: `adj[i][j]` is the weight of the edge `i -> j`
    /// (or `Ternary::Neutral` when no such edge exists).
    pub adj: Vec<Vec<Ternary>>,
    /// `true` for a directed graph; `false` for an undirected graph (where the
    /// adjacency matrix is kept symmetric by [`add_edge`](Self::add_edge)).
    pub directed: bool,
}

impl TernaryGraph {
    /// Create a new graph with `n` vertices and no edges.
    ///
    /// Every entry of the adjacency matrix starts as `Ternary::Neutral` (the
    /// "no edge" value). Runtime is `O(n²)` and the resulting matrix uses
    /// `O(n²)` memory.
    pub fn new(n: usize, directed: bool) -> Self {
        TernaryGraph {
            n,
            adj: vec![vec![Ternary::Neutral; n]; n],
            directed,
        }
    }

    /// Add (or overwrite) the edge `u -> v` with the given ternary `weight`.
    ///
    /// Passing `Ternary::Neutral` removes any existing edge, since `Neutral`
    /// represents "no edge" everywhere else in the API.
    ///
    /// For undirected graphs the mirror entry `v -> u` is set to the same
    /// weight, keeping the matrix symmetric. Self-loops (`u == v`) are
    /// permitted: they set the single diagonal entry and count as exactly one
    /// edge in [`edge_count`](Self::edge_count).
    ///
    /// # Panics
    ///
    /// Panics if `u >= n` or `v >= n`.
    pub fn add_edge(&mut self, u: usize, v: usize, weight: Ternary) {
        self.adj[u][v] = weight;
        if !self.directed && u != v {
            self.adj[v][u] = weight;
        }
    }

    /// Query the weight of the edge `u -> v`.
    ///
    /// Returns `Ternary::Neutral` when no edge is present (including the
    /// diagonal `u == v` if no self-loop was added).
    ///
    /// # Panics
    ///
    /// Panics if `u >= n` or `v >= n`.
    pub fn edge(&self, u: usize, v: usize) -> Ternary {
        self.adj[u][v]
    }

    /// List the `(neighbor, weight)` pairs of `v`, in ascending index order.
    ///
    /// Only non-`Neutral` entries are returned. A self-loop on `v` (if present)
    /// appears as the pair `(v, weight)`.
    ///
    /// # Panics
    ///
    /// Panics if `v >= n`.
    pub fn neighbors(&self, v: usize) -> Vec<(usize, Ternary)> {
        (0..self.n)
            .filter(|&u| self.adj[v][u] != Ternary::Neutral)
            .map(|u| (u, self.adj[v][u]))
            .collect()
    }

    /// Total number of edges in the graph.
    ///
    /// For a **directed** graph this is the number of non-`Neutral` cells in
    /// the adjacency matrix. For an **undirected** graph each off-diagonal
    /// edge is stored twice (as `(u,v)` and `(v,u)`) and is therefore counted
    /// once; self-loops are counted once.
    pub fn edge_count(&self) -> usize {
        let mut off_diagonal = 0usize;
        let mut self_loops = 0usize;
        for i in 0..self.n {
            for j in 0..self.n {
                if self.adj[i][j] != Ternary::Neutral {
                    if i == j {
                        self_loops += 1;
                    } else {
                        off_diagonal += 1;
                    }
                }
            }
        }
        if self.directed {
            off_diagonal + self_loops
        } else {
            // Each undirected edge is stored symmetrically as (i,j) and (j,i).
            self_loops + off_diagonal / 2
        }
    }

    /// Degree of `v` — the number of non-`Neutral` entries in row `v`,
    /// including a self-loop if one is present.
    ///
    /// # Panics
    ///
    /// Panics if `v >= n`.
    pub fn degree(&self, v: usize) -> usize {
        self.neighbors(v).len()
    }

    /// Return the diagonal degree matrix `D` as a dense `n × n` `f64` matrix.
    ///
    /// `D[i][i]` is the degree of vertex `i`; off-diagonal entries are zero.
    pub fn degree_matrix(&self) -> Vec<Vec<f64>> {
        let mut d = vec![vec![0.0f64; self.n]; self.n];
        for (i, row) in d.iter_mut().enumerate().take(self.n) {
            row[i] = self.degree(i) as f64;
        }
        d
    }

    /// Return a copy of the adjacency matrix with each ternary value widened
    /// to `f64` (`-1.0`, `0.0`, `+1.0`).
    pub fn adjacency_f64(&self) -> Vec<Vec<f64>> {
        self.adj
            .iter()
            .map(|row| row.iter().map(|&v| i8::from(v) as f64).collect())
            .collect()
    }

    /// Compute the combinatorial graph Laplacian `L = D − A`.
    ///
    /// `L[i][i]` is the degree of `i`, and `L[i][j] = -A[i][j]` for `i != j`.
    /// `L` is symmetric positive semi-definite for undirected graphs.
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

    /// Compute the symmetric normalized Laplacian
    /// `L_norm = D^(−1/2) · L · D^(−1/2)`.
    ///
    /// For vertices with non-zero degree `d_i`:
    ///
    /// ```text
    /// L_norm[i][i] = 1
    /// L_norm[i][j] = -A[i][j] / sqrt(d_i · d_j)   (i != j, edge present)
    /// ```
    ///
    /// For an **isolated** vertex (`d_i = 0`) the diagonal entry is left as
    /// `1.0`, matching the convention used by NetworkX's
    /// `normalized_laplacian_matrix`. Eigenvalues of `L_norm` lie in `[0, 2]`.
    pub fn normalized_laplacian(&self) -> Vec<Vec<f64>> {
        let d = self.degree_matrix();
        let a = self.adjacency_f64();
        let mut l = vec![vec![0.0f64; self.n]; self.n];

        // Precompute D^{-1/2}, using 0 for isolated vertices (their rows/columns
        // will produce only zero off-diagonals and the special diagonal below).
        let mut d_inv_sqrt = vec![0.0f64; self.n];
        for i in 0..self.n {
            if d[i][i] > 0.0 {
                d_inv_sqrt[i] = 1.0 / d[i][i].sqrt();
            }
        }

        for i in 0..self.n {
            for j in 0..self.n {
                l[i][j] = if i == j {
                    1.0
                } else if a[i][j] != 0.0 {
                    // Both endpoints necessarily have non-zero degree here,
                    // so d_inv_sqrt[i] and d_inv_sqrt[j] are both finite.
                    -a[i][j] * d_inv_sqrt[i] * d_inv_sqrt[j]
                } else {
                    0.0
                };
            }
        }

        l
    }
}

/// Single-source shortest paths via Bellman-Ford over ternary edge weights.
///
/// Distances are returned as `Option<f64>`:
/// - `Some(d)` — the shortest signed distance `source -> i` is `d`.
/// - `None`    — vertex `i` is unreachable from `source`, **or** `i` is
///   reachable from a negative cycle (in which case no finite shortest path
///   exists).
///
/// Because ternary edges may be `±1`, the algorithm handles negative weights
/// directly (unlike Dijkstra). Self-loops of weight `-1` are detected as
/// negative cycles.
///
/// # Panics
///
/// Panics if `source >= graph.n`.
pub fn shortest_paths(graph: &TernaryGraph, source: usize) -> Vec<Option<f64>> {
    let n = graph.n;
    let mut dist: Vec<Option<f64>> = vec![None; n];
    if n == 0 {
        return dist;
    }
    dist[source] = Some(0.0);

    // Relax edges n-1 times.
    for _ in 0..n.saturating_sub(1) {
        let mut changed = false;
        for u in 0..n {
            if dist[u].is_none() {
                continue;
            }
            let du = dist[u].unwrap();
            for (v, w) in graph.neighbors(u) {
                let new_dist = du + i8::from(w) as f64;
                if dist[v].is_none() || new_dist < dist[v].unwrap() {
                    dist[v] = Some(new_dist);
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }

    // Detect vertices whose distance can still decrease — these are on, or
    // reachable from, a negative cycle.
    let mut has_neg_cycle = vec![false; n];
    for _ in 0..n {
        let mut changed = false;
        for u in 0..n {
            if dist[u].is_none() {
                continue;
            }
            let du = dist[u].unwrap();
            for (v, w) in graph.neighbors(u) {
                let new_dist = du + i8::from(w) as f64;
                if dist[v].is_none() || new_dist < dist[v].unwrap() {
                    has_neg_cycle[v] = true;
                    dist[v] = Some(new_dist);
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }

    // Propagate negative-cycle reachability forward along every edge type.
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

/// All-pairs shortest paths via Floyd-Warshall over ternary edge weights.
///
/// Returns an `n × n` matrix where `dist[i][j] = Some(d)` is the shortest
/// signed distance `i -> j`, or `None` when `j` is unreachable from `i`.
/// The diagonal is `Some(0.0)` unless `i` lies on a negative cycle, in which
/// case it is `None` (consistent with [`shortest_paths`]).
pub fn all_pairs_shortest_paths(graph: &TernaryGraph) -> Vec<Vec<Option<f64>>> {
    let n = graph.n;
    let mut dist = vec![vec![None::<f64>; n]; n];

    for (i, row) in dist.iter_mut().enumerate().take(n) {
        row[i] = Some(0.0);
    }
    for (u, row) in dist.iter_mut().enumerate().take(n) {
        for (v, w) in graph.neighbors(u) {
            row[v] = Some(i8::from(w) as f64);
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

    // Floyd-Warshall detects a negative cycle when any dist[i][i] < 0 after
    // the main loop; vertices whose distance can be made arbitrarily negative
    // are exactly those reachable from such a cycle, which we expose as None.
    let mut on_neg_cycle = vec![false; n];
    for i in 0..n {
        if matches!(dist[i][i], Some(d) if d < 0.0) {
            on_neg_cycle[i] = true;
        }
    }
    // Propagate reachability forward along all edges.
    for _ in 0..n {
        let mut changed = false;
        for u in 0..n {
            if !on_neg_cycle[u] {
                continue;
            }
            for (v, _) in graph.neighbors(u) {
                if !on_neg_cycle[v] {
                    on_neg_cycle[v] = true;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    // Null out distances to/from every vertex on or downstream of a negative
    // cycle. We collect the affected indices first so we never hold two
    // mutable borrows of `dist` simultaneously.
    let affected: Vec<usize> = (0..n).filter(|&i| on_neg_cycle[i]).collect();
    for &i in &affected {
        // dist[i][j] = None  AND  dist[j][i] = None  for every j.
        for row in dist.iter_mut().take(n) {
            row[i] = None;
        }
        for cell in dist[i].iter_mut().take(n) {
            *cell = None;
        }
    }

    dist
}

/// Community detection via weighted label propagation.
///
/// Each vertex starts in its own community. On every iteration we re-label
/// vertex `v` with the label that has the largest summed edge weight from
/// `v`'s neighbors (`+1` for [`Ternary::Positive`] edges, `-1` for
/// [`Ternary::Negative`]). Positive edges therefore pull neighbours into the
/// same community while negative edges push them apart. Iteration stops after
/// `max_iters` rounds or once a full pass produces no relabeling.
///
/// Ties are broken by **smallest label id**, and vote tallying uses a
/// [`BTreeMap`](std::collections::BTreeMap) so the result is fully
/// **deterministic** for any given input graph.
///
/// The returned labels are relabeled to a contiguous `0..k` range.
pub fn label_propagation(graph: &TernaryGraph, max_iters: usize) -> Vec<usize> {
    let n = graph.n;
    let mut labels: Vec<usize> = (0..n).collect();

    for _ in 0..max_iters {
        let mut changed = false;
        for v in 0..n {
            // Tally weighted votes per neighbor-label. Using a BTreeMap makes
            // the iteration order deterministic, so ties resolve consistently.
            let mut votes: std::collections::BTreeMap<usize, f64> =
                std::collections::BTreeMap::new();
            for (u, w) in graph.neighbors(v) {
                let label = labels[u];
                let weight = i8::from(w) as f64;
                *votes.entry(label).or_insert(0.0) += weight;
            }

            if votes.is_empty() {
                continue;
            }

            // Pick the label with the highest weighted vote; ties broken by
            // smallest label id (BTreeMap iterates ascending, and we only
            // replace on a strictly-greater vote).
            let mut best_label = labels[v];
            let mut best_score = f64::NEG_INFINITY;
            for (&label, &score) in votes.iter() {
                if score > best_score {
                    best_score = score;
                    best_label = label;
                }
            }

            if best_label != labels[v] {
                labels[v] = best_label;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }

    // Relabel to consecutive integers starting at 0.
    let mut label_map = std::collections::HashMap::new();
    let mut next = 0;
    for l in labels.iter_mut() {
        let entry = *label_map.entry(*l).or_insert_with(|| {
            let v = next;
            next += 1;
            v
        });
        *l = entry;
    }

    labels
}

/// Compute the (signed) modularity `Q` of a community assignment.
///
/// Uses the standard Newman definition
///
/// ```text
/// Q = (1 / 2m) · Σ_{i,j : c_i = c_j} ( A[i][j] − k_i·k_j / (2m) )
/// ```
///
/// where `m` is [`edge_count`](TernaryGraph::edge_count), `A` is the
/// (signed) adjacency matrix, and `k_i = Σ_j A[i][j]` is the signed weighted
/// degree. Positive edges within a community raise `Q`; negative edges within
/// a community lower it. Returns `0.0` for an edgeless graph.
///
/// **Note on signed graphs.** Strictly modularity is defined for non-negative
/// weights; the formula above is the natural signed extension (Gómez,
/// Jensen & Arenas, *Phys. Rev. L* 2009) and is useful for ranking partitions
/// even when negative edges are present.
pub fn modularity(graph: &TernaryGraph, communities: &[usize]) -> f64 {
    let n = graph.n;
    let m = graph.edge_count() as f64;
    if m == 0.0 {
        return 0.0;
    }

    // Signed weighted degree k_i = sum_j A[i][j].
    let mut k = vec![0.0f64; n];
    for (i, ki) in k.iter_mut().enumerate().take(n) {
        for j in 0..n {
            *ki += i8::from(graph.adj[i][j]) as f64;
        }
    }

    let mut q = 0.0;
    for i in 0..n {
        for j in 0..n {
            if communities[i] == communities[j] {
                let a_ij = i8::from(graph.adj[i][j]) as f64;
                q += a_ij - (k[i] * k[j]) / (2.0 * m);
            }
        }
    }

    q / (2.0 * m)
}

/// Power iteration on `laplacian`, returning its top eigenvector (unit-norm).
///
/// Used by [`spectral_clustering`] after shifting the matrix so that the
/// smallest eigenvalue of the original Laplacian becomes the largest of the
/// shifted matrix. Convergence is faster when the top two eigenvalues are
/// well-separated. Returns the initial vector unchanged when the matrix is
/// all-zero (to avoid dividing by zero).
fn power_iteration_laplacian(laplacian: &[Vec<f64>], n: usize, max_iters: usize) -> Vec<f64> {
    let mut v: Vec<f64> = (0..n).map(|i| (i as f64 + 1.0) / n as f64).collect();
    let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm > 0.0 {
        for x in v.iter_mut() {
            *x /= norm;
        }
    }

    for _ in 0..max_iters {
        let mv: Vec<f64> = laplacian
            .iter()
            .map(|row| row.iter().zip(v.iter()).map(|(&a, &b)| a * b).sum())
            .collect();
        let new_norm = mv.iter().map(|x| x * x).sum::<f64>().sqrt();
        if new_norm < 1e-15 {
            break;
        }
        for (vi, mvi) in v.iter_mut().zip(mv.iter()) {
            *vi = *mvi / new_norm;
        }
    }

    v
}

/// Spectral clustering via the signed Laplacian.
///
/// Builds the **signed** Laplacian (`D − A` where `D_ii = Σ_j |A_ij|`), shifts
/// it by `2n·I` so that the smallest eigenvalues of `L` become the largest of
/// `2n·I − L`, and uses power iteration with Hotelling deflation to extract
/// the first `k + 1` eigenvectors. The first (trivial, constant) eigenvector
/// is discarded; the next ones form the embedding used for clustering.
///
/// - For `k = 2`, vertices are split at the median of the Fiedler vector
///   (the second-smallest eigenvector).
/// - For `k > 2`, each vertex is assigned to `sign(eigenvector)·hash mod k`
///   over the first `k + 1` non-trivial eigenvectors. This is a lightweight
///   stand-in for row-wise k-means; it gives reasonable partitions when the
///   spectral gap is clean but is **not** a full Ng-Jordan-Weiss pipeline.
/// - `k ≤ 1` returns every vertex in community `0`.
///
/// The result is deterministic for a given input graph.
///
/// # Complexity
///
/// `O(k · T · n²)` where `T` is the number of power-iteration steps (300).
pub fn spectral_clustering(graph: &TernaryGraph, k: usize) -> Vec<usize> {
    if graph.n == 0 {
        return Vec::new();
    }
    if k <= 1 {
        return vec![0; graph.n];
    }

    let n = graph.n;

    // Signed adjacency as f64 (-1.0 / 0.0 / +1.0).
    let mut adj_f64 = vec![vec![0.0f64; n]; n];
    for (i, row) in adj_f64.iter_mut().enumerate().take(n) {
        for (j, cell) in row.iter_mut().enumerate().take(n) {
            *cell = i8::from(graph.adj[i][j]) as f64;
        }
    }

    // Signed degree: d_i = sum_j |A_ij|.
    let mut degree = vec![0.0f64; n];
    for (i, di) in degree.iter_mut().enumerate().take(n) {
        for &v in adj_f64[i].iter().take(n) {
            *di += v.abs();
        }
    }

    // Signed Laplacian L = D - A.
    let mut laplacian = vec![vec![0.0f64; n]; n];
    for (i, row) in laplacian.iter_mut().enumerate().take(n) {
        for (j, cell) in row.iter_mut().enumerate().take(n) {
            *cell = if i == j { degree[i] } else { -adj_f64[i][j] };
        }
    }

    // Shift by 2n*I so power iteration on (2n*I - L) finds the *smallest*
    // eigenvectors of L.
    let shift = 2.0 * n as f64;
    let mut shifted = laplacian.clone();
    for (i, row) in shifted.iter_mut().enumerate().take(n) {
        for (j, cell) in row.iter_mut().enumerate().take(n) {
            *cell = if i == j {
                shift - laplacian[i][j]
            } else {
                -laplacian[i][j]
            };
        }
    }

    // Deflation-based extraction of the top (k+1) eigenvectors of `shifted`,
    // which are the smallest (k+1) eigenvectors of L.
    let mut current = shifted.clone();
    let mut eigenvectors = Vec::new();

    for _ in 0..(k + 1).min(n) {
        let eigvec = power_iteration_laplacian(&current, n, 300);
        let eigval: f64 = eigvec
            .iter()
            .zip(current.iter())
            .map(|(&v, row)| {
                v * row
                    .iter()
                    .zip(eigvec.iter())
                    .map(|(&a, &b)| a * b)
                    .sum::<f64>()
            })
            .sum();

        // Hotelling deflation: subtract the rank-1 component.
        for (i, row) in current.iter_mut().enumerate().take(n) {
            for (j, cell) in row.iter_mut().enumerate().take(n) {
                *cell -= eigval * eigvec[i] * eigvec[j];
            }
        }

        eigenvectors.push(eigvec);
    }

    if eigenvectors.is_empty() {
        return vec![0; n];
    }

    // The very first eigenvector of L for a connected graph is the constant
    // vector (eigenvalue 0); the Fiedler vector is the second-smallest.
    let fiedler_idx = if eigenvectors.len() > 1 { 1 } else { 0 };
    let fiedler = &eigenvectors[fiedler_idx];

    if k == 2 {
        // Split at the median of the Fiedler vector. Use total_cmp so the sort
        // is total even in the (theoretically impossible here) presence of NaN.
        let mut sorted: Vec<f64> = fiedler.to_vec();
        sorted.sort_by(f64::total_cmp);
        let median = sorted[n / 2];
        return fiedler
            .iter()
            .map(|&v| if v >= median { 1 } else { 0 })
            .collect();
    }

    // For k > 2: hash the sign pattern of the first k+1 eigenvectors and
    // reduce mod k. (Lightweight proxy for row-wise k-means.)
    let mut assignments = vec![0usize; n];
    for (i, slot) in assignments.iter_mut().enumerate().take(n) {
        let mut hash = 0usize;
        for (d, ev) in eigenvectors.iter().enumerate() {
            if ev[i] >= 0.0 {
                hash |= 1 << (d % 16);
            }
        }
        *slot = hash % k;
    }

    assignments
}

/// Connected components of the **positive-weight** subgraph.
///
/// Performs a BFS over edges whose weight is [`Ternary::Positive`] only —
/// neutral (no edge) and negative (inhibitory) edges do **not** connect their
/// endpoints in this view. Returns a vector of component ids in `0..k`,
/// assigned contiguously in BFS-start order (vertex 0 is always in
/// component 0). Isolated vertices each receive their own component id.
pub fn connected_components(graph: &TernaryGraph) -> Vec<usize> {
    let n = graph.n;
    let mut components = vec![usize::MAX; n];
    let mut component_id = 0;

    // BFS from each unvisited node using only positive edges.
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
        assert!(
            q > 0.0,
            "Modularity should be positive for good partition, got {}",
            q
        );
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
        assert!(
            same_in_first >= 2,
            "First triangle should be mostly one cluster: {:?}",
            labels
        );
        assert!(
            same_in_second >= 2,
            "Second triangle should be mostly one cluster: {:?}",
            labels
        );
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
