"""Exact adversarial instances for the diagonal partition-reduction theorem.

Certificates prove containment over the original generator list and an exact
continuous output slice. They do not prove the universal reduction theorem.
"""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

from checks import ROOT, jsonable, membership, signed


def anchor(g, n, clipped=True):
    a, b = max(g[:n]), min(-v for v in g[n:])
    return (max(Q(0), a), min(Q(1), b)) if clipped else (a, b)


def combination(G, target):
    weights = tuple(min((Q(0),) + tuple(x-y for x, y in zip(target, g))) for g in G)
    assert max(weights) == 0
    assert tuple(max(g[j]+w for g, w in zip(G, weights))
                 for j in range(len(target))) == target
    return weights


def input_exact_slice(G, n):
    """Sufficient whole-slice criterion, never a sampled exactness assertion."""
    assert (all(len(set(g[:n])) == 1 for g in G)
            or all(len(set(g[n:])) == 1 for g in G))
    lo = -min(max(g[n+j] for g in G) for j in range(n))
    hi = min(max(g[j] for g in G) for j in range(n))
    for g in G:
        a, b = anchor(g, n, clipped=False)
        a, b = max(a, lo), min(b, hi)
        assert a > b or (0 <= a <= b <= 1)


def reduce_case(name, G, n, inexact_point=None):
    G = tuple(tuple(Q(v) for v in g) for g in G)
    intervals = tuple(anchor(g, n) for g in G)
    assert sum(a < b for a, b in intervals) <= len(G)-2
    partition, sources = [Q(0)], []
    while partition[-1] < 1:
        reached = partition[-1]
        candidates = [(b, -k) for k, (a, b) in enumerate(intervals)
                      if a <= reached < b]
        assert candidates, "Uncovered interval"
        b, minus_k = max(candidates)
        sources.append(-minus_k)
        partition.append(b)
    output = (signed((Q(0),)*n), signed((Q(1),)*n)) + tuple(
        (a,)*n + (-b,)*n for a, b in zip(partition, partition[1:]))
    assert len(output) <= len(G)
    coefficients = tuple(combination(G, g) for g in output)
    if inexact_point is None:
        input_exact_slice(G, n)
    else:
        assert not (len(set(inexact_point)) == 1 and 0 <= inexact_point[0] <= 1)
        assert membership(G, signed(inexact_point))
    return dict(
        name=name, dimension=n, generators=G, partition=partition,
        source_indices=sources, output_generators=output, coefficients=coefficients,
        input_scope="exact_one_uniform_block" if inexact_point is None else "inexact_enclosure",
        inexact_point=() if inexact_point is None else inexact_point,
        inexact_coefficients=() if inexact_point is None else combination(G, signed(inexact_point)),
    )


def cases():
    for n in (1, 2, 3, 5):
        for q in (1, 2, 4):
            endpoints = (signed((Q(0),)*n), signed((Q(1),)*n))
            cells = tuple((Q(i, q), Q(i+1, q)) for i in range(q))
            overlaps = tuple((max(Q(0), a-Q(1, 4*q)),
                              min(Q(1), b+Q(1, 4*q))) for a, b in cells)
            canonical = tuple((a,)*n+(-b,)*n for a, b in cells)
            overlap = tuple((a,)*n+(-b,)*n for a, b in overlaps)
            variants = {
                "canonical": endpoints+canonical,
                "overlap": endpoints+overlap,
                "positive_uniform": endpoints+tuple(
                    (a,)*n+tuple(-b-Q(j, n+1) for j in range(n)) for a, b in overlaps),
                "negative_uniform": endpoints+tuple(
                    tuple(a-Q(j, n+1) for j in range(n))+(-b,)*n for a, b in overlaps),
                "remote_endpoints": ((Q(-1),)*n+(Q(0),)*n,
                                     (Q(2),)*n+(Q(-1),)*n)+canonical,
                "clutter": endpoints+overlap+(
                    signed((Q(1, 2),)*n), (Q(1),)*n+(Q(0),)*n,
                    overlap[0], (Q(1, 3),)*n+(Q(-2, 3),)*n),
            }
            for variant, G in variants.items():
                yield reduce_case(f"n{n}_q{q}_{variant}", G, n)
            point = (Q(-1),) if n == 1 else (Q(1, 4), Q(3, 4))+(Q(1, 4),)*(n-2)
            yield reduce_case(f"n{n}_q{q}_inexact", endpoints+canonical+(signed(point),), n, point)


def main():
    destination = ROOT/"results/theory/partition_reduction.json"
    interchange = destination.with_suffix(".toml")
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Invalidate both products before doing work, including interrupted runs.
    destination.write_text('{"complete": false}\n', encoding="utf-8", newline="\n")
    interchange.write_text("complete = false\n", encoding="utf-8", newline="\n")
    records = list(cases())
    sources = (Path(__file__), Path(__file__).with_name("checks.py"))
    report = dict(
        complete=True, cases=records,
        source_sha256={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sources},
        summary=dict(cases=len(records), exact_inputs=sum(c["input_scope"].startswith("exact") for c in records),
                     inexact_inputs=sum(c["input_scope"] == "inexact_enclosure" for c in records),
                     output_generator_certificates=sum(len(c["coefficients"]) for c in records)),
        scope="Finite exact inclusion and whole-slice certificates; not universal theorem verification.",
    )
    lines = ['schema = "diagonal-partition-certificates-v1"', "complete = true"]
    for case in records:
        lines.extend(["", "[[cases]]"])
        lines.extend(f"{key} = {json.dumps(jsonable(value))}" for key, value in case.items())
    interchange.write_text("\n".join(lines)+"\n", encoding="utf-8", newline="\n")
    destination.write_text(json.dumps(jsonable(report), indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps(report["summary"]))


if __name__ == "__main__":
    main()
