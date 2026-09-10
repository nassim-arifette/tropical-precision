# Reproducing the results

## Requirements

- Python 3.11+; tested with 3.11.0.
- Julia; tested with 1.12.7. Required only to execute the separate certificate verifiers.
- Standard libraries only. Internet access is unnecessary after installing the runtimes.
- CPU execution; no GPU or proprietary solver. The validated host is Windows. Cross-platform reproduction has not yet been independently tested.

Run commands from the repository root. The entry points also work when called by absolute path from another working directory.

## Quick checks

```sh
python -B examples/relu.py
python -B tests/run_tests.py
```

The first command recomputes the ReLU examples using exact fractions and checks the retained `results/relu.json`. The second verifies all source hashes, mathematical JSON/TOML payloads, certificate counts and failure-control records, then tests deliberate corruptions on a temporary copy. Expect a final `All tests passed.` message. These checks take seconds on the development host.

## Separate certificate verification

```sh
python -B scripts/reproduce.py --certificates-only
```

This executes the Julia verifiers on the saved rational certificates, including malformed inputs and failed-command controls. It refreshes verifier records but skips the finite Python discovery grids. Output includes exact certificate counts and ends with `Reproduction completed successfully.`

## Full reproduction

```sh
python -B scripts/reproduce.py
```

This executes the nine Python studies in dependency order, runs the Julia verifiers, checks provenance, and regenerates the neural results. Allow several minutes on a desktop CPU; runtime varies with hardware. The boundary-budget study depends on the finite-budget results generated earlier in this command.

Result paths are relative to this repository, independent of the invoking shell's current directory. The command overwrites generated result files. Scientific payloads are deterministic; timestamps and runtime metadata may change. Text is UTF-8 with LF endings to preserve source hashes across Git checkouts.

For a custom Julia location:

```sh
python -B scripts/reproduce.py --julia /path/to/julia
```

Do not use `-O` or set `PYTHONOPTIMIZE`: the mathematical studies deliberately reject disabled assertions. The neural example uses explicit exceptions and remains checked under optimization.

## Interpreting evidence

Each study writes JSON results and TOML certificate inputs. Julia records successful verification only after its positive and negative checks finish. Missing, malformed or incomplete records are not successful verification. The Python provenance check compares complete mathematical payloads and the current executable sources. Written proofs belong to the accompanying paper and are not included in executable-source hashes.

Finite grids do not establish continuous optima. Complete scalar interval covers give bounds conditional on the paper's universal reductions. A successful run does not establish novelty, a universal proof, or trained-network performance. See [the result map](results.md).
