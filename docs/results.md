# Results and evidence

| Source | Output under `results/theory/` | What is checked |
|---|---|---|
| `src/checks.py` | `exact_checks.*` | Seven finite invariant problems: 1,248 contexts, 46 observations, and 21 rational diagonal witnesses |
| `src/affine_profiles.py` | `affine_profiles.*` | Five finite point/direction searches and a rational attaining example; search maxima are lower bounds |
| `src/partition_reduction.py` | `partition_reduction.*` | 84 partition extraction cases and 364 normalized containment combinations |
| `src/affine_upper.py` | `affine_upper.*` | Six algebraic identities and finite adversarial checks of the three-dimensional affine bound |
| `src/affine_dimensions.py` | `affine_dimensions.*` | Dimension-dependent affine supports and 20 rational certificates |
| `src/affine_cells.py` | `affine_cells.*` | 30 interior-cell certificates and finite searches |
| `src/affine_two_level.py` | `affine_two_level.*` | Two-level domination checks and 203 witnesses, including ties and zero coefficients |
| `src/affine_certified.py` | `affine_certified.*` | Complete rational scalar covers and partition bounds: 12 cells, 12 budgets, 21,052 intervals |
| `src/affine_boundary.py` | `affine_boundary.*` | 25 boundary controls, eight partition controls, 23 budget roots, six transitions and exact algebraic identities |

The separate verifiers in `tests/` check the corresponding TOML inputs with Julia's `Rational{BigInt}` arithmetic. Files ending in `_julia.toml`, plus `julia_verification.toml`, record these checks. Failure-control records cover malformed input, invalid certificates and errors after otherwise successful mathematical checks.

## ReLU examples

`examples/relu.py` writes `results/relu.json`. It partitions the entire one-dimensional input interval at every ReLU crossing and uses exact rational arithmetic.

| Example | True minimum | Coarse merge | Refined merge |
|---|---:|---:|---:|
| Nonlinear suffix | `1/8` | Exact minimum `-1/8` | Exact minimum `1/8` |
| Affine readout | `9/4` | Witness value `-1/12` | Guaranteed lower bound `1/8` |

The nonlinear result checks the geometric cover and an attaining reachable point. The affine refined bound uses the written affine budget theorem; it is not an independently optimized abstract minimum. The generator budgets three and four count the diagonal summary, not every generator in the merged state. Refinement replaces an anchor. Adding an anchor while keeping the coarse one cannot remove its spurious state.

The examples establish a loss under a specified merge-before-query policy. Exact queries on separate regions can avoid this loss. They are constructed rational networks, not trained-model benchmarks.
