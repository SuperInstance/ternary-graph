# PLUG_AND_PLAY — Graph

> Graph algorithms on ternary-weighted edges

## 🚀 Quick Start

Add to your `Cargo.toml`:

```toml
[dependencies]
ternary-graph = { git = "https://github.com/SuperInstance/ternary-graph" }
```

Use in your code:

```rust
use ternary_graph::TernaryGraph;

let mut g = TernaryGraph::new(4);
g.add_edge(0, 1, TernaryGraph::Pos);
let path = g.shortest_path(0, 3);
```

## 📚 Available Documentation

| Document | Description |
|----------|-------------|
| `docs/FROM_BINARY.md` | Understanding ternary concepts as a binary programmer |
| `docs/MIGRATION.md` | Version migration guide |
| `docs/FUTURE-INTEGRATION.md` | Planned features and roadmap |

## 🔗 Integration

This crate is part of the [SuperInstance ternary fleet](https://github.com/SuperInstance). It uses the canonical `Ternary` type from `ternary-types` for cross-crate compatibility.

## 📄 License

MIT
