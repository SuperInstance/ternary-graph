# Architecture — ternary-graph

> *Internal design and data flow.*

## Overview

This crate implements ternary {-1, 0, +1} semantics for the `graph` domain.
It is one of ~280 ternary crates in the SuperInstance fleet, all sharing Z₃ arithmetic
from [ternary-core](https://github.com/SuperInstance/ternary-core).

## Core Types

- **`TernaryGraph`**

## Key Functions

- `to_i8()`
- `to_f64()`
- `new()`
- `add_edge()`
- `edge()`
- `neighbors()`
- `edge_count()`
- `degree()`

## Ternary Mapping

| Value | Meaning |
|-------|---------|
| +1 | Attractive edge |
| 0  | No edge |
| -1 | Repulsive edge |

## Source Structure

1 Rust source file(s) in `src/`.
Language: Rust

## Cross-Repo References

- [ternary-core](https://github.com/SuperInstance/ternary-core) — shared Z₃ traits
- [ternary-types](https://github.com/SuperInstance/ternary-types) — type-level encodings
- [Full SuperInstance fleet](https://github.com/orgs/SuperInstance/repositories?q=ternary)
