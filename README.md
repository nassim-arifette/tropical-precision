# Tropical Precision

Exact-arithmetic code for contextual precision and generator budgets in tropical abstract domains, with two constructed ReLU certification examples.

The accompanying manuscript by Nassim Arifette is available as a [PDF](paper/main.pdf), with [LaTeX sources, the figure and build instructions](paper/). This is the version intended for arXiv.

## Quick start

Python 3.11 or newer is sufficient for these read-only checks:

```sh
python -B examples/relu.py
python -B tests/run_tests.py
```

The ReLU example reports exact true/coarse/refined margins of `1/8`, `-1/8`, and `1/8` for the nonlinear classifier. For the affine classifier, it reports a true minimum of `9/4`, a spurious-point value of `-1/12`, and a theorem-derived refined lower bound of `1/8`.

The tests verify result provenance, complete JSON/TOML agreement, and rejection of altered certificates. They do not prove the universal mathematical theorems.

## Reproduce the results

Install Julia in addition to Python. The tested versions are Python 3.11.0 and Julia 1.12.7 on Windows; both use standard libraries only. No Python package installation, datasets, trained weights, GPU or external solver is needed.

```sh
python -B scripts/reproduce.py --certificates-only
python -B scripts/reproduce.py
```

The first command rechecks saved certificates in Julia. The second regenerates all nine mathematical studies, checks their certificates, and reruns the ReLU examples. Both refresh verification records in `results/`. Pass `--julia PATH` if Julia is not on PATH. Do not use Python's `-O` option or `PYTHONOPTIMIZE`.

See [reproduction instructions](docs/reproducibility.md) for expected outputs and [result descriptions](docs/results.md) for the scope of each computation.

## Repository structure

| Directory | Contents |
|---|---|
| `src/` | Exact Python computations and certificate generation |
| `examples/` | Runnable rational ReLU classifiers |
| `tests/` | Julia certificate verifiers and Python integrity tests |
| `scripts/` | Reproduction entry point |
| `results/` | Reference results and verification records |
| `docs/` | Reproduction instructions and interpretation of results |
| `paper/` | Manuscript PDF, LaTeX sources, figure and build instructions |

## Scope

The computations support explicit constructions, finite tests and certified scalar bounds. Universal results depend on the written mathematical arguments. The neural examples isolate a specified hidden-state merge; they do not measure performance on trained models or implement a general neural verifier.

## Citation and license

Citation metadata for the software is in [CITATION.cff](CITATION.cff). When reporting results, include the commit identifier used for reproduction.

The code and accompanying repository documentation are available under the [MIT license](LICENSE).
The manuscript, its LaTeX sources and its figure are covered by the separate [paper copyright notice](paper/LICENSE).
