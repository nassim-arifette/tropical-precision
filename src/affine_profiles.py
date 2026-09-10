"""Targeted exact affine-loss discovery for diagonal anchor cubes in dimension 3.

Context optimization is continuous and exact, reduced to scalar breakpoints.
Only x and affine coefficient directions are sampled on explicitly finite grids.
Reported values are certified lower bounds, never claims of global optimality.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
import hashlib
import json
from pathlib import Path

from checks import ROOT, diagonal_generators, dot, export_witnesses, jsonable, reconstruction


def weight_shell(denominator):
    for a in range(-denominator, denominator + 1):
        for b in range(-denominator + abs(a), denominator - abs(a) + 1):
            c = denominator - abs(a) - abs(b)
            for signed_c in ((c, -c) if c else (0,)):
                yield a, b, signed_c


def sector_support(x, w, k, direction, D):
    """Integer numerator of max w·h in a signed sector of the continuous cube.

The denominators of x and w are D and W. For fixed s=|h_k| the other
coordinates maximize independently. The resulting concave piecewise-affine
function has breakpoints at distance-to-boundary values; all have denominator D.
"""
    radius = D - x[k] if direction > 0 else x[k]
    bounds = [D - x[j] if w[j] >= 0 else x[j] for j in range(3)]
    breaks = {0, radius}
    breaks.update(min(radius, bounds[j]) for j in range(3) if j != k)
    best, best_h = 0, (0, 0, 0)
    for s in breaks:
        h = tuple(direction * s if j == k
                  else (min(s, bounds[j]) if w[j] > 0
                        else -min(s, bounds[j]) if w[j] < 0 else 0)
                  for j in range(3))
        value = sum(w[j] * h[j] for j in range(3))
        if value > best:
            best, best_h = value, h
    return best, best_h


def profile(a, b, subdivisions, W):
    values = tuple(a + (b - a) * Q(i, subdivisions) for i in range(subdivisions + 1))
    from math import lcm
    D = lcm(*(t.denominator for t in values))
    grid = tuple(int(t * D) for t in values)
    weights = tuple(weight_shell(W))
    best, best_data, count = 0, None, 0
    for x in combinations_with_replacement(grid, 3):
        if x[0] == x[2]:
            continue
        missing = tuple((k, 1) for k in range(3) if x[k] > x[0]) + tuple(
            (k, -1) for k in range(3) if x[k] < x[2])
        for w in weights:
            count += 1
            deficit = min(0, sum(w)) * D - sum(w[j] * x[j] for j in range(3))
            if deficit <= best:
                continue
            supports = [sector_support(x, w, k, direction, D) for k, direction in missing]
            loss = min([deficit] + [v for v, _ in supports])
            if loss > best:
                best, best_data = loss, (x, w, missing, supports)
    assert best_data is not None
    x, w, missing, supports = best_data
    z = tuple(Q(v, D) for v in x)
    coeffs = tuple(Q(v, W) for v in w)
    T = tuple(tuple(Q(x[j] + h[j], D) for j in range(3)) for _, h in supports)
    partition = tuple(sorted({Q(0), a, b, Q(1)}))
    G = diagonal_generators(3, partition)
    concrete = min([min(Q(0), sum(coeffs))] + [dot(coeffs, t) for t in T])
    assert concrete - dot(coeffs, z) == Q(best, D * W)
    cert = reconstruction(G, T, z)
    return dict(
        interval=(a, b), coordinate_subdivisions=subdivisions,
        coordinate_denominator=D, weight_denominator=W, weight_directions=len(weights),
        point_direction_evaluations=count, certified_lower_bound=Q(best, D * W),
        target=z, affine_coefficients=coeffs, context=T, generators=G,
        concrete_minimum=concrete, target_value=dot(coeffs, z),
        reconstruction_certificate=cert,
        scope="Exact continuous-context optimum for each tested point/direction; finite-grid lower bound over points/directions."
    )


def family_witness(h):
    u, v = h * (1 - h), h * (2 - h)
    z = (u, u, v)
    T = ((1 - h, 1 - h, Q(1)), (Q(0), 2 * u, h), (2 * u, Q(0), h))
    w = (1 / (2 * (2 - h)), 1 / (2 * (2 - h)), -(1 - h) / (2 - h))
    G = diagonal_generators(3, (Q(0), Q(1)))
    assert all(0 <= t <= 1 for p in T for t in p)
    assert sum(abs(t) for t in w) == 1
    assert all(dot(w, t) == 0 for t in T)
    gap = h * (1 - h) ** 2 / (2 - h)
    assert dot(w, z) == -gap and sum(w) > 0
    return dict(
        parameter=h, target=z, context=T, affine_coefficients=w, generators=G,
        certified_gap=gap, reconstruction_certificate=reconstruction(G, T, z),
        scope="Analytic family valid for 0<h<1; this rational certificate does not prove a global affine optimum."
    )


def main():
    specs = [
        (Q(0), Q(1), 8, 24),
        (Q(0), Q(1, 4), 4, 24),
        (Q(1, 4), Q(1, 2), 4, 24),
        (Q(3, 8), Q(5, 8), 4, 24),
        (Q(7, 16), Q(9, 16), 4, 48),
    ]
    results = []
    destination = ROOT / "results/theory/affine_profiles.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text('{"complete": false}\n', encoding="utf-8", newline="\n")
    destination.with_suffix(".toml").write_text("complete = false\n", encoding="utf-8", newline="\n")
    try:
        for spec in specs:
            result = profile(*spec)
            results.append(result)
            print(json.dumps(jsonable({k: result[k] for k in
                                       ("interval", "certified_lower_bound", "target", "affine_coefficients", "point_direction_evaluations")})), flush=True)
        sources = [Path(__file__), Path(__file__).with_name("checks.py")]
        report = dict(
            complete=True, arithmetic="exact integers and fractions.Fraction",
            scientific_question="How do affine witness losses vary with anchor location and width?",
            results=results,
            source_sha256={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
            total_point_direction_evaluations=sum(r["point_direction_evaluations"] for r in results),
            analytic_family=family_witness(Q(3, 8)),
            scope="Targeted discovery, not global minimax optimization or a performance benchmark."
        )
        assert report["analytic_family"]["certified_gap"] == Q(75, 832) > Q(1, 12)
    except Exception as exc:
        destination.write_text(json.dumps(dict(complete=False, error=repr(exc))) + "\n", encoding="utf-8", newline="\n")
        destination.with_suffix(".toml").write_text(
            'schema = "bounded-context-witnesses-v1"\ncomplete = false\n', encoding="utf-8", newline="\n")
        raise
    destination.write_text(json.dumps(jsonable(report), indent=2) + "\n", encoding="utf-8", newline="\n")
    witnesses = []
    for item in results + [report["analytic_family"]]:
        is_family = "parameter" in item
        gap = item["certified_gap"] if is_family else item["certified_lower_bound"]
        witnesses.append(dict(
            dimension=3, generator_count=len(item["generators"]),
            generators=item["generators"], target=item["target"], context=item["context"],
            affine_coefficients=item["affine_coefficients"],
            minimum_concrete_union=Q(0) if is_family else item["concrete_minimum"],
            target_value=-gap if is_family else item["target_value"], certified_gap=gap,
            membership_certificate=item["reconstruction_certificate"],
        ))
    export_witnesses(dict(continuous_diagonal_certificates=witnesses), destination.with_suffix(".toml"))


if __name__ == "__main__":
    main()
