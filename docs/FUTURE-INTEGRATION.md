# Future Integration: ternary-graph

## Current State
Implements graph algorithms on ternary-weighted edges {-1, 0, +1}: `TernaryGraph` with adjacency matrix, shortest paths (handling negative edges via Bellman-Ford), community detection via modularity optimization, graph Laplacian computation, spectral clustering, and connected components.

## Integration Opportunities

### With ternary-cell / room-as-codespace
Rooms are vertices; passages are edges; edge weights encode passage type: +1 (preferred route), 0 (no connection), -1 (restricted/blocked). `shortest_path()` routes agents between rooms. `community_detection()` identifies room neighborhoods — groups of rooms that form natural zones. `spectral_clustering()` on the graph Laplacian produces higher-quality partitions than simple connectivity.

### With ternary-planning
`PlanGraph` from `ternary-planning` becomes a `TernaryGraph`. Task dependencies are +1 edges (strong), 0 edges (no dependency), -1 edges (conflict). `shortest_path()` finds the critical path through the plan. `connected_components()` identifies independent subplans that can execute in parallel.

### With ternary-markov
The transition matrix of a `TernaryMarkov` chain IS a weighted graph. `community_detection()` on this graph identifies metastable state clusters — sets of states that the chain visits frequently before transitioning to another cluster. The graph Laplacian's eigenvalues bound the mixing time from below.

## Potential in Mature Systems
In PLATO, the entire fleet is a `TernaryGraph`. Constructs are vertices; communication channels are edges weighted by reliability (+1 = reliable, 0 = offline, -1 = adversarial). `shortest_path()` routes messages through the most reliable path. When edges flip from +1 to -1 (construct compromise), `connected_components()` identifies the affected subgraph for isolation. The graph Laplacian's spectrum monitors fleet connectivity — a vanishing second eigenvalue indicates the fleet is fragmenting.

## Cross-Pollination Ideas
**Music × Graph:** Voice-leading graphs: vertices are chords, edges are ternary voice-leading distances (+1 smooth, 0 parallel, -1 leap). `shortest_path()` finds the smoothest chord progression between two harmonies. `community_detection()` identifies tonal regions (keys). The graph Laplacian's eigenvectors are the harmonic "modes." Connects to `ternary-music` and `flux-algebra-rs`.

**Social × Graph:** Agent interaction graphs with ternary trust: +1 (trusted), 0 (unknown), -1 (distrusted). `spectral_clustering()` identifies trust communities. `shortest_path()` with negative weights finds routes that avoid distrusted intermediaries — secure message routing.

## Dependencies for Next Steps
- Dynamic graph updates (edges flip as room conditions change)
- Sparse adjacency representation for large room fleets
- Integration with `ternary-consensus` for distributed graph algorithms
