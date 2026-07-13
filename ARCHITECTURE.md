# Architecture — ternary-graph

> *Internal design and data flow.*

## Overview

`ternary-graph` is a pure-Rust library of graph algorithms whose edge weights
are restricted to the balanced ternary digit set `{-1, 0, +1}` (re-exported
from [`ternary-types`](https://github.com/SuperInstance/ternary-types) as
`ternary_graph::Ternary`).

The implementation lives in a single file, `src/lib.rs`, and is `#![forbid(unsafe_code)]`.

## Core Type

- **`TernaryGraph`** — Dense adjacency-matrix graph with `n` vertices. The
  matrix is `Vec<Vec<Ternary>>`; `Ternary::Neutral` represents "no edge".

## Public API Surface

### Methods on `TernaryGraph`

`new`, `add_edge`, `edge`, `neighbors`, `degree`, `edge_count`,
`adjacency_f64`, `degree_matrix`, `laplacian`, `normalized_laplacian`.

### Free functions

`shortest_paths` (Bellman-Ford, single source),
`all_pairs_shortest_paths` (Floyd-Warshall),
`connected_components` (BFS over positive-weight edges),
`label_propagation` (deterministic weighted label propagation),
`modularity` (signed Newman Q),
`spectral_clustering` (signed Laplacian + power iteration + Hotelling deflation).

## Ternary Edge Mapping

| Value | `Ternary` variant | Meaning |
|-------|-------------------|---------|
| `+1`  | `Positive`        | Attractive / excitatory / trust edge |
| `0`   | `Neutral`         | No edge |
| `-1`  | `Negative`        | Repulsive / inhibitory / distrust edge |

## Source Structure

- `src/lib.rs` — All algorithms and unit tests.
- `tests/` — None (all tests are in `src/lib.rs`).
- `docs/FUTURE-INTEGRATION.md` — Cross-crate integration ideas.

## Cross-Repo References

- [ternary-types](https://github.com/SuperInstance/ternary-types) — shared
  `Ternary` enum (`Negative`, `Neutral`, `Positive`) and `TernaryError`.
- [Full SuperInstance fleet](https://github.com/orgs/SuperInstance/repositories?q=ternary)
