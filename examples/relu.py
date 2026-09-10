"""Exact, bounded regression of the two neural-interface examples.

Uses Fraction and full one-dimensional ReLU region splitting, not sampling.
The universal generator and budget claims still depend on the written proofs.
Default mode verifies the retained report; --write replaces it after a new run.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def require(value, message):
    if not value:
        raise ValueError(message)


def dot(a, b):
    require(len(a) == len(b), "Vector dimension mismatch")
    return sum((x * y for x, y in zip(a, b)), Q(0))


def signed(x):
    return tuple(x) + tuple(-v for v in x)


def combination(generators, coefficients):
    require(len(generators) == len(coefficients) and max(coefficients) == 0,
            "Combination is not normalized")
    return tuple(max(g[j] + c for g, c in zip(generators, coefficients))
                 for j in range(len(generators[0])))


def member(x, generators):
    coefficients = [min(Q(0), *(a - b for a, b in zip(x, g))) for g in generators]
    return max(coefficients) == 0 and combination(generators, coefficients) == tuple(x)


@dataclass(frozen=True)
class Piece:
    lo: Q
    hi: Q
    slope: tuple
    intercept: tuple

    def at(self, t):
        return tuple(a * t + b for a, b in zip(self.slope, self.intercept))


def layer(pieces, matrix, bias, relu=True):
    require(len(matrix) == len(bias), "Bias dimension mismatch")
    output = []
    for p in pieces:
        slope = tuple(dot(row, p.slope) for row in matrix)
        intercept = tuple(dot(row, p.intercept) + b for row, b in zip(matrix, bias))
        knots = {p.lo, p.hi}
        if relu:
            for a, b in zip(slope, intercept):
                if a and p.lo < -b / a < p.hi:
                    knots.add(-b / a)
        ordered = sorted(knots)
        for lo, hi in zip(ordered, ordered[1:]):
            active = [not relu or a * (lo + hi) / 2 + b >= 0
                      for a, b in zip(slope, intercept)]
            output.append(Piece(lo, hi,
                                tuple(a if keep else Q(0) for a, keep in zip(slope, active)),
                                tuple(b if keep else Q(0) for b, keep in zip(intercept, active))))
    for left, right in zip(output, output[1:]):
        require(left.hi == right.lo and left.at(left.hi) == right.at(right.lo),
                "Gap or discontinuity in exact region propagation")
    return output


def minimum(pieces, observation):
    return min(observation(p.at(t)) for p in pieces for t in (p.lo, p.hi))


def anchor_cover(generators, n):
    boxes = []
    for g in generators:
        lower = tuple(max(Q(0), g[j]) for j in range(n))
        upper = tuple(min(Q(1), -g[n+j]) for j in range(n))
        if all(l <= u for l, u in zip(lower, upper)):
            boxes.append((lower, upper))
    return boxes


def distance_to_box(x, box):
    lower, upper = box
    return max(max(l - v, v - u, Q(0)) for v, l, u in zip(x, lower, upper))


def nonlinear():
    zero, one = (Q(0),) * 2, (Q(1),) * 2
    p3 = [signed(zero), signed(one), (Q(0), Q(0), Q(-1), Q(-1))]
    p4 = p3[:2] + [(Q(0), Q(0), Q(-1, 2), Q(-1, 2)),
                    (Q(1, 2), Q(1, 2), Q(-1), Q(-1))]
    uv = [signed(zero), signed((Q(0), Q(3, 4))), (Q(0), Q(0), Q(0), Q(-3, 4))]
    uh = [signed(one), signed((Q(1, 4), Q(1))), (Q(1, 4), Q(1), Q(-1), Q(-1))]
    coarse = list(dict.fromkeys(p3 + uv + uh))
    fine = list(dict.fromkeys(p4 + uv + uh))
    require((len(coarse), len(fine)) == (7, 8), "Wrong full representation counts")
    target = (Q(1, 4), Q(3, 4))
    witness_g = p3 + [signed((Q(1, 4), Q(1))), signed((Q(0), Q(3, 4)))]
    coeff = [Q(-3, 4), Q(-3, 4), Q(0), Q(-1, 4), Q(-1, 4)]
    require(combination(witness_g, coeff) == signed(target), "Nonlinear join witness failed")
    require(member(signed(target), coarse) and not member(signed(target), fine),
            "Spurious state not removed")
    distance_bound = min(distance_to_box(target, box) for box in anchor_cover(fine, 2))
    require(distance_bound == Q(1, 4), "Refined anchor-cover distance is not 1/4")
    require(member(signed((Q(1, 2), Q(1, 2))), fine), "Missing attaining diagonal point")
    prefix = [Piece(Q(0), Q(3), (Q(1),), (Q(0),))]
    prefix = layer(prefix, [(Q(1),)] * 3, [Q(0), Q(-1), Q(-2)])
    prefix = layer(prefix, [(Q(0), Q(1), Q(-7, 4)), (Q(-3, 4), Q(7, 4), Q(-1))],
                   [Q(0), Q(3, 4)])
    expected = [(Q(0), Q(3, 4)), (Q(0), Q(0)), (Q(1), Q(1)), (Q(1, 4), Q(1))]
    require(len(prefix) == 3, "Prefix no longer has the three specified regions")
    for i, p in enumerate(prefix):
        require(p.lo == i and p.hi == i+1 and p.at(p.lo) == expected[i]
                and p.at(p.hi) == expected[i+1], "Wrong nonlinear-prefix segment")
    suffix = layer(prefix, [(Q(1), Q(0)), (Q(-1), Q(0)),
                            (Q(0), Q(1)), (Q(0), Q(-1))],
                   [Q(-1, 4), Q(1, 4), Q(-3, 4), Q(3, 4)])
    suffix = layer(suffix, [(Q(1), Q(1), Q(-1), Q(-1)),
                            (Q(0), Q(0), Q(1), Q(1))], [Q(0), Q(0)])
    suffix = layer(suffix, [(Q(1), Q(1)), (Q(0), Q(0))], [Q(0), Q(1, 8)], relu=False)
    true_min = minimum(suffix, lambda z: z[0]-z[1])
    require(true_min == Q(1, 8), "Wrong exact nonlinear classifier margin")
    return dict(architecture=[1, 3, 2, 4, 2, 2], input_interval=[Q(0), Q(3)],
                prefix_pieces=len(prefix), full_network_pieces=len(suffix),
                true_minimum=true_min, coarse_minimum=Q(-1, 8),
                refined_minimum=distance_bound-Q(1, 8), displayed_merge_generators=[7, 8],
                exact_suffix_evaluation=True, suffix_solver_benchmark=False)


def affine(matrix=None):
    rows = [[3, -6, 3, 2, -2, 0, 0, 0, 0], [3, -4, 1, -2, 2, 2, -2, 0, 0],
            [3, -5, 2, 0, 0, 2, -2, -2, 2], [3, -5, 2, 0, 0, 0, 0, 2, -2]]
    matrix = matrix or [[Q(v, 3) for v in row] for row in rows]
    contexts = [(Q(0), Q(2, 3), Q(1, 3), Q(1, 3)),
                (Q(2, 3), Q(0), Q(1, 3), Q(1, 3)),
                (Q(2, 3), Q(2, 3), Q(1), Q(1, 3)),
                (Q(2, 3), Q(2, 3), Q(1, 3), Q(1))]
    vertices = [(Q(0),) * 4, (Q(1),) * 4] + [t for t in contexts for _ in range(2)]
    pieces = layer([Piece(Q(0), Q(1), (Q(1),), (Q(0),))], [(Q(9),)] * 9,
                   [Q(-i) for i in range(9)])
    pieces = layer(pieces, matrix, [Q(0)] * 4)
    require(len(pieces) == 9, "Affine prefix does not have nine pieces")
    w = tuple(map(Q, (5, 4, -4, -4)))
    for i, p in enumerate(pieces):
        require(p.lo == Q(i, 9) and p.hi == Q(i+1, 9), "Wrong knot interval")
        require(p.at(p.lo) == vertices[i] and p.at(p.hi) == vertices[i+1], "Wrong affine-prefix vertex")
        for t in (p.lo, p.hi):
            z = p.at(t)
            require(all(0 <= v <= 1 for v in z), "Hidden state leaves cube")
            require(z[0]+z[1] == z[2]+z[3] and dot(w, z) == z[0], "Lost linear invariant")
    true_min = minimum(pieces, lambda z: Q(9, 4) + dot(w, z))
    true_max = -minimum(pieces, lambda z: -Q(9, 4) - dot(w, z))
    target = (Q(1, 3), Q(1, 3), Q(2, 3), Q(2, 3))
    anchor = (Q(0),)*4 + (Q(-1),)*4
    require(combination([anchor] + [signed(t) for t in contexts], [Q(0)]+[Q(-1, 3)]*4)
            == signed(target), "Affine join witness failed")
    witness_margin = Q(9, 4) + dot(w, target)
    require((true_min, true_max, witness_margin) == (Q(9, 4), Q(13, 4), Q(-1, 12)),
            "Wrong affine margin bounds")
    norm = sum(map(abs, w))
    lower = true_min - norm * Q(1, 8)
    require(lower == Q(1, 8), "Wrong theorem-derived affine lower bound")
    return dict(architecture=[1, 9, 4, 2], input_interval=[Q(0), Q(1)],
                affine_pieces=9, true_range=[true_min, true_max], witness_upper_bound=witness_margin,
                theorem_derived_refined_lower_bound=lower, exact_refined_minimum_computed=False,
                dependency="Written affine budget theorem at n=4, M=4; normalized loss 1/8",
                readout_l1=norm, selected_regions=5, separately_checked_transitions=4)


def serialize(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {k: serialize(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(v) for v in value]
    return value


def run():
    negative = []
    altered = [[Q(v, 3) for v in row] for row in
               [[3, -5, 3, 2, -2, 0, 0, 0, 0], [3, -4, 1, -2, 2, 2, -2, 0, 0],
                [3, -5, 2, 0, 0, 2, -2, -2, 2], [3, -5, 2, 0, 0, 0, 0, 2, -2]]]
    for name, check in (("altered_affine_weight", lambda: affine(altered)),
                        ("nonnormalized_witness", lambda: combination([(Q(0),)], [Q(-1)]))):
        try:
            check()
        except ValueError:
            negative.append(name)
        else:
            raise ValueError(f"Negative control unexpectedly passed: {name}")
    files = [Path(__file__).resolve()]
    return serialize(dict(schema="neural-interface-regression-v1", complete=True,
                          nonlinear=nonlinear(), affine=affine(), negative_controls=negative,
                          files_sha256={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                        for p in files},
                          scope="Exact example arithmetic and complete 1D region propagation; not a mechanized proof, general hull optimizer, or trained-network benchmark."))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    report = run()
    output = ROOT / "results/relu.json"
    if args.write:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    else:
        require(json.loads(output.read_text(encoding="utf-8")) == report,
                "Retained neural-interface report differs from current exact run")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
