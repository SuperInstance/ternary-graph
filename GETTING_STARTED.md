# Getting Started — ternary-graph

> *Estimated time: 5 minutes*

## Prerequisites

- **Rust 1.75+** (MSRV)
- Cargo (included with Rust)

## Installation

```toml
[dependencies]
ternary-graph = { git = "https://github.com/SuperInstance/ternary-graph.git" }
```

Or from source:

```bash
git clone https://github.com/SuperInstance/ternary-graph.git
cd ternary-graph
cargo build --release
cargo test
```

## Core Concept

This crate implements graph algorithms over edges weighted by the balanced
ternary digit set `{-1, 0, +1}`. `+1` is an excitatory / trust / attractive
edge; `-1` is inhibitory / distrust / repulsive; `0` (the `Ternary::Neutral`
variant) means "no edge".

## Quick Example

```rust
use ternary_graph::{TernaryGraph, Ternary};

let mut g = TernaryGraph::new(4, false);
g.add_edge(0, 1, Ternary::Positive);
g.add_edge(1, 2, Ternary::Positive);
g.add_edge(2, 3, Ternary::Positive);

let dist = ternary_graph::shortest_paths(&g, 0);
assert_eq!(dist[3], Some(3.0));
```

## Running Tests

```bash
cargo test
```

## Next Steps

- [README.md](./README.md) — Full API reference and mathematical background
- [ARCHITECTURE.md](./ARCHITECTURE.md) — Internal design
- [PLUG_AND_PLAY.md](./PLUG_AND_PLAY.md) — Integration
- [CONTRIBUTING.md](./CONTRIBUTING.md) — Contributing
