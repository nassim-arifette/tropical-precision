"""Exact, bounded checks of contextual formulas and dimension-sensitive loss.

Only the Python standard library is used. Finite-K cases enumerate EVERY
context subset and compute joined membership directly by normalized residuation.
They are compared with a separate coordinate-sector calculation. The continuous
diagonal controls carry explicit rational inclusion/affine-margin certificates;
they are not numerical proofs that a representative has the advertised slice.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Exact checks require assertions: remove -O and PYTHONOPTIMIZE")

import argparse
from datetime import datetime, timezone
from fractions import Fraction as Q
from itertools import product
import hashlib
import json
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]


def point(xs):
    if any(isinstance(x, float) for x in xs):
        raise TypeError("Use exact integer or rational inputs")
    return tuple(Q(x) for x in xs)


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), Q(0))


def linf(a, b):
    return max((abs(x - y) for x, y in zip(a, b)), default=Q(0))


def image(A, x):
    return tuple(dot(row, x) for row in A)


def signed(x):
    return tuple(x) + tuple(-t for t in x)


def signed_matrix(n):
    eye = [tuple(Q(int(i == j)) for j in range(n)) for i in range(n)]
    return tuple(eye + [tuple(-t for t in row) for row in eye])


def membership(G, z):
    """Direct normalized max-plus membership, no sector API used."""
    residuals = [min((Q(0),) + tuple(t - g for t, g in zip(z, v))) for v in G]
    return max(residuals) == 0 and all(
        max(v[j] + c for v, c in zip(G, residuals)) == z[j]
        for j in range(len(z))
    )


def missing_and_witnesses(A, G, K, x):
    """Independent form: sector inequalities in concrete source coordinates."""
    rows = ((Q(0),) * len(x),) + A
    augmented_G = [(Q(0),) + g for g in G]
    missing = []
    witnesses = {}
    for i, ai in enumerate(rows):
        normals = [tuple(aj[k] - ai[k] for k in range(len(x))) for aj in rows]
        covered = any(
            all(dot(normal, x) >= g[j] - g[i] for j, normal in enumerate(normals))
            for g in augmented_G
        )
        if not covered:
            missing.append(i)
            witnesses[i] = tuple(
                j for j, t in enumerate(K)
                if all(dot(normal, tuple(t[k] - x[k] for k in range(len(x)))) <= 0
                       for normal in normals)
            )
            assert witnesses[i], "Every sector fiber contains x itself"
    return tuple(missing), witnesses


def finite_case(name, A, G, K, observations):
    A, G, K = tuple(map(point, A)), tuple(map(point, G)), tuple(map(point, K))
    assert len(set(K)) == len(K)
    assert all(len(g) == len(A) for g in G)
    assert all(len(row) == len(K[0]) for row in A)
    images = tuple(image(A, x) for x in K)
    S = tuple(i for i, z in enumerate(images) if membership(G, z))
    assert S
    distances = tuple(tuple(linf(x, t) for t in K) for x in K)
    observed_values = [tuple(dot(point(w), x) for x in K) for w in observations]
    # Any assignment on a finite K is a continuous observation. These two
    # nonlinear assignments test (3) beyond distance or affine special cases.
    observed_values += [
        tuple(Q((i * i + 3 * i) % 7 - 3, 5) for i in range(len(K))),
        tuple(max(abs(x[j] - Q(j + 1, len(x) + 2)) for j in range(len(x))) for x in K),
    ]
    best_lip = Q(0)
    best_observation = [Q(0)] * len(observed_values)
    best_lip_context = ()
    contexts = 1 << len(K)
    for mask in range(contexts):
        T = tuple(i for i in range(len(K)) if mask & (1 << i))
        concrete = S + T
        generators = G + tuple(images[i] for i in T)
        joined = tuple(i for i, z in enumerate(images) if membership(generators, z))
        assert set(concrete) <= set(joined)
        loss = max(min(distances[i][j] for j in concrete) for i in joined)
        if loss > best_lip:
            best_lip, best_lip_context = loss, T
        for k, values in enumerate(observed_values):
            gap = min(values[j] for j in concrete) - min(values[i] for i in joined)
            best_observation[k] = max(best_observation[k], gap)

    formula_lip = Q(0)
    formula_observation = [Q(0)] * len(observed_values)
    for i, x in enumerate(K):
        missing, witnesses = missing_and_witnesses(A, G, K, x)
        assert bool(missing) == (i not in S)
        radius = min(
            [min(distances[i][j] for j in S)]
            + [max(distances[i][j] for j in witnesses[c]) for c in missing]
        )
        formula_lip = max(formula_lip, radius)
        for k, values in enumerate(observed_values):
            bound = min(
                [min(values[j] for j in S)]
                + [max(values[j] for j in witnesses[c]) for c in missing]
            ) - values[i]
            formula_observation[k] = max(formula_observation[k], bound)
    assert best_lip == formula_lip, (name, best_lip, formula_lip)
    assert best_observation == formula_observation, (name, best_observation, formula_observation)
    return dict(
        name=name, A=A, generators=G, K=K, concrete_indices=S,
        contexts_enumerated=contexts, observations_checked=len(observed_values),
        observation_values=observed_values,
        direct_lipschitz_loss=best_lip, sector_lipschitz_loss=formula_lip,
        direct_observation_losses=best_observation,
        sector_observation_losses=formula_observation,
        worst_lipschitz_context_indices=best_lip_context,
    )


def diagonal_generators(n, partition):
    return (
        signed((Q(0),) * n),
        signed((Q(1),) * n),
    ) + tuple((a,) * n + (-b,) * n for a, b in zip(partition, partition[1:]))


def reconstruction(G, T, z):
    family = tuple(G) + tuple(signed(t) for t in T)
    target = signed(z)
    coefficients = tuple(
        min((Q(0),) + tuple(a - b for a, b in zip(target, g))) for g in family
    )
    rebuilt = tuple(max(g[j] + c for g, c in zip(family, coefficients))
                    for j in range(len(target)))
    assert max(coefficients) == 0
    assert rebuilt == target
    return dict(coefficients=coefficients, reconstruction=rebuilt)


def interval_certificate(n, partition):
    G = diagonal_generators(n, partition)
    # Select the largest part of an anchor interval in the middle half.
    candidates = [(max(a, Q(1, 4)), min(b, Q(3, 4)), i + 2)
                  for i, (a, b) in enumerate(zip(partition, partition[1:]))]
    a, b, index = max(candidates, key=lambda x: x[1] - x[0])
    h, r = b - a, Q(1, 4)
    q = len(G) - 2
    assert Q(1, 2 * q) <= h <= Q(1, 2)
    z = (a, a, b) + (a,) * (n - 3)
    directions = (
        (1,) * n,
        (-1, 1, -1) + (-1,) * (n - 3),
        (1, -1, -1) + (-1,) * (n - 3),
    )
    T = tuple(tuple(z[j] + r * direction[j] for j in range(n)) for direction in directions)
    w = tuple(c / (4 - h) for c in (Q(1), Q(1), -(2 - h)) + (Q(0),) * (n - 3))
    assert sum(abs(t) for t in w) == 1
    assert all(0 <= c <= 1 for t in T for c in t)
    assert all(g <= t for g, t in zip(G[index], signed(z)))
    # Stronger than membership: explicitly verify the three fixed shifts and
    # a coefficient-zero anchor used by the universal lower-bound proof.
    rebuilt = tuple(max((G[index][j],) + tuple(signed(t)[j] - r for t in T))
                    for j in range(2 * n))
    assert rebuilt == signed(z)
    min_diagonal = min(Q(0), sum(w))
    concrete = min((min_diagonal,) + tuple(dot(w, t) for t in T))
    value_z = dot(w, z)
    gap = concrete - value_z
    assert gap == h / (4 * (4 - h))
    assert gap >= Q(1, 32 * q)
    # Several exact slice inclusions including every partition endpoint.
    sampled_diagonal = sorted(set(partition + tuple((a + b) / 2 for a, b in zip(partition, partition[1:]))))
    assert all(membership(G, signed((t,) * n)) for t in sampled_diagonal)
    return dict(
        dimension=n, generator_count=len(G), partition=partition, generators=G,
        anchor_index=index, middle_interval=(a, b), target=z, context=T,
        affine_coefficients=w, minimum_diagonal=min_diagonal,
        minimum_concrete_union=concrete, target_value=value_z, certified_gap=gap,
        proved_universal_lower_bound=Q(1, 32 * q),
        membership_certificate=reconstruction(G, T, z),
        slice_scope="Exact diagonal slice follows from equal positive/negative coordinate groups and the written reconstruction proof, not sampled inclusions.",
    )


def run():
    half = Q(1, 2)
    K1 = tuple((Q(i, 2),) for i in range(-2, 3))
    K01 = tuple((Q(i, 2),) for i in range(5))
    K2 = tuple(product((Q(0), half, Q(1)), repeat=2))
    affine2 = [(Q(a, 3), Q(b, 3)) for a, b in product((-1, 0, 1), repeat=2) if (a, b) != (0, 0)]
    cases = [
        finite_case("unsigned_line_no_interior_row", [(1,)], [(0,)], K1, [(1,), (-1,)]),
        finite_case("signed_line_off_lift_generator", [(1,), (-1,)], [(0, 0), (-1, -1)], K1, [(1,), (-1,)]),
        finite_case("ordered_exact_lift", [(1,), (2,)], [(0, 0), (1, 1)], K01, [(1,), (-1,)]),
        finite_case("rank_deficient_lift", [(1, 0)], [(0,)], tuple(product((0, half, 1), (0, 1))), affine2),
        finite_case("unsigned_plane", [(1, 0), (0, 1)], [(0, 0)], K2, affine2 + [(Q(-1, 3), Q(2, 3))]),
        finite_case("signed_diagonal_dimension_two", signed_matrix(2),
                    diagonal_generators(2, (Q(0), Q(1))), K2, affine2),
    ]
    z = point((Q(1, 4), Q(1, 4), Q(3, 4)))
    T = (point((half, half, 1)), point((0, half, half)), point((half, 0, half)))
    G = diagonal_generators(3, (Q(0), Q(1)))
    w = point((Q(1, 3), Q(1, 3), Q(-1, 3)))
    assert all(dot(w, t) == 0 for t in T)
    assert dot(w, z) == Q(-1, 12)
    reconstruction(G, T, z)
    cases.append(finite_case("signed_diagonal_dimension_three_affine_failure",
                            signed_matrix(3), G, ((0, 0, 0), (1, 1, 1), z) + T, [w]))
    assert cases[-1]["direct_observation_losses"][0] >= Q(1, 12)
    assert all(v == 0 for v in cases[-2]["direct_observation_losses"][:len(affine2)])

    certificates = []
    for n in (3, 4, 8):
        for q in (1, 2, 3, 7):
            uniform = tuple(Q(i, q) for i in range(q + 1))
            graded = tuple(Q(i * i, q * q) for i in range(q + 1))
            for partition in dict.fromkeys((uniform, graded)):
                certificates.append(interval_certificate(n, partition))
    return dict(
        schema="bounded-context-exact-controls-v1",
        complete=True, arithmetic="fractions.Fraction",
        scope="Exhaustive finite-context identities and explicit continuous-diagonal witnesses; no benchmark or universal numerical proof.",
        finite_cases=cases, continuous_diagonal_certificates=certificates,
        summary=dict(
            finite_cases=len(cases),
            contexts_enumerated=sum(c["contexts_enumerated"] for c in cases),
            fixed_observations_checked=sum(c["observations_checked"] for c in cases),
            diagonal_certificates=len(certificates), discrepancies=0,
        ),
    )


def jsonable(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, (tuple, list)):
        return [jsonable(y) for y in x]
    if isinstance(x, dict):
        return {k: jsonable(v) for k, v in x.items()}
    return x


def export_witnesses(report, destination):
    """Simple TOML interchange for a separately implemented Julia verifier."""
    lines = ['schema = "bounded-context-witnesses-v1"', "complete = true"]
    for case in report["continuous_diagonal_certificates"]:
        lines += ["", "[[witnesses]]"]
        for name in ("dimension", "generator_count"):
            lines.append(f'{name} = {case[name]}')
        for name in ("generators", "target", "context", "affine_coefficients"):
            lines.append(f'{name} = {json.dumps(jsonable(case[name]))}')
        for name in ("minimum_concrete_union", "target_value", "certified_gap"):
            lines.append(f'{name} = "{case[name]}"')
        lines.append('coefficients = ' + json.dumps(jsonable(case["membership_certificate"]["coefficients"])))
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/theory/exact_checks.json")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text('{"complete": false}\n', encoding="utf-8", newline="\n")
    args.output.with_suffix(".toml").write_text("complete = false\n", encoding="utf-8", newline="\n")
    sources = [Path(__file__)]
    metadata = dict(
        time_utc=datetime.now(timezone.utc).isoformat(),
        python=sys.version, platform=platform.platform(),
        source_sha256={str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        command="python -B src/checks.py",
    )
    try:
        report = {**run(), **metadata}
    except Exception as exc:
        report = dict(complete=False, error=repr(exc), **metadata)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
        args.output.with_suffix(".toml").write_text(
            'schema = "bounded-context-witnesses-v1"\ncomplete = false\n', encoding="utf-8", newline="\n"
        )
        raise
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(jsonable(report), indent=2) + "\n", encoding="utf-8", newline="\n")
    export_witnesses(report, args.output.with_suffix(".toml"))
    print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
