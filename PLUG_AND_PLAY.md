# PLUG_AND_PLAY — ternary-graph

> *Integration guide.*

## Dependency

```toml
[dependencies]
ternary-graph = { git = "https://github.com/SuperInstance/ternary-graph.git" }
```

The crate is a library (no binary target) and has no Cargo features — every
algorithm is always available.

## Minimal example

```rust
use ternary_graph::{TernaryGraph, Ternary};

// Directed graph so the negative edge does not form a length-2 negative cycle.
let mut g = TernaryGraph::new(3, true);
g.add_edge(0, 1, Ternary::Positive);
g.add_edge(1, 2, Ternary::Negative);

let dist = ternary_graph::shortest_paths(&g, 0);
assert_eq!(dist[2], Some(0.0)); // 1 + (-1) = 0
```

> **Note**: in an *undirected* graph a negative edge forms a length-2 negative
> cycle (going back and forth costs `-2`), so `shortest_paths` will return
> `None` for any vertex on or reachable from that edge. Use a directed graph
> when you want negative edge weights without triggering negative-cycle
> detection.

## Compatibility

- **Rust edition**: 2021 (MSRV 1.75)
- **Targets**: every tier-1 Rust target (no platform-specific code; `#![forbid(unsafe_code)]`)
