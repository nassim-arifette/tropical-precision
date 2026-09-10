"""Exact dimension controls for the three-generator affine minimax.

Continuous sector supports are exact; finite point/direction grids are only
falsification checks. Rational attaining witnesses are separately verified in
Julia. The universal coefficient/Jensen bound is a written mathematical proof.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
import hashlib
import json

from checks import (ROOT, diagonal_generators, dot, export_witnesses, finite_case,
                    jsonable, reconstruction, signed_matrix)
from affine_profiles import sector_support as support_three


def weight_shell(n, total):
    """All integer vectors with l1 norm total, exactly once."""
    if n == 1:
        yield (total,)
        if total:
            yield (-total,)
        return
    for first in range(-total, total + 1):
        for rest in weight_shell(n - 1, total - abs(first)):
            yield (first,) + rest


def sector_support(x, w, k, sign, D):
    """Integer numerator of the exact continuous cube support."""
    n = len(x)
    assert len(w) == n and sign in (-1, 1)
    radius = D - x[k] if sign == 1 else x[k]
    bounds = tuple(D-x[j] if w[j] >= 0 else x[j] for j in range(n))
    breaks = {0, radius}
    breaks.update(min(radius, bounds[j]) for j in range(n) if j != k)
    best = 0
    for s in breaks:
        value = sign*w[k]*s + sum(abs(w[j])*min(s, bounds[j])
                                  for j in range(n) if j != k)
        best = max(best, value)
    return best


def sign_count_bound(p, k):
    assert p >= 2 and k >= 2
    return 1/(Q(4)+Q(1, p-1)+Q(1, k-1))


def dimension_bound(n):
    assert n != 3  # Its exact value is irrational, handled algebraically below.
    return Q(0) if n <= 2 else sign_count_bound(n//2, n-n//2)


def finite_checks(n, D, W):
    weights = tuple(weight_shell(n, W))
    assert len(weights) == len(set(weights))
    assert all(sum(abs(v) for v in w) == W for w in weights)
    pairs, positive, covered_extreme_checks = 0, 0, 0
    three_support_comparisons = 0
    best, maximizer = Q(0), None
    sign_counts = set()
    for x in combinations_with_replacement(range(D+1), n):
        missing = [(j, 1) for j in range(n) if x[j] > x[0]]
        missing += [(j, -1) for j in range(n) if x[j] < x[-1]]
        for w in weights:
            pairs += 1
            supports = [sector_support(x, w, j, sign, D) for j, sign in missing]
            if n == 3:
                previous = [support_three(x, w, j, sign, D)[0] for j, sign in missing]
                assert supports == previous
                three_support_comparisons += len(supports)
            base = min(0, sum(w))*D - sum(v*y for v, y in zip(w, x))
            numerator = min([base]+supports)
            delta = Q(numerator, D*W)
            if delta > best:
                best, maximizer = delta, dict(point=x, direction=w)
            if delta <= 0:
                continue
            positive += 1
            assert n >= 3
            assert all(2*abs(v) < W for v in w if v)
            for j, v in enumerate(w):
                if v > 0:
                    assert numerator <= (W-2*v)*x[j]
                    covered_extreme_checks += int(x[j] == x[-1])
                elif v < 0:
                    assert numerator <= (W+2*v)*(D-x[j])
                    covered_extreme_checks += int(x[j] == x[0])
            p, k = sum(v > 0 for v in w), sum(v < 0 for v in w)
            sign_counts.add((p, k))
            A, C = Q(sum(v for v in w if v > 0), W), Q(-sum(v for v in w if v < 0), W)
            aggregate = min(A,C)/(1+sum(Q(abs(v), W-2*abs(v)) for v in w if v))
            equalized = min(A,C)/(1+A/(1-2*A/p)+C/(1-2*C/k))
            assert delta <= aggregate <= equalized
            if min(p,k) == 1:
                c = A if p == 1 else C
                singleton = c*(1-2*c)/(2*(1-c)**2)
                assert delta <= singleton <= Q(1,8)
                if n == 3:
                    assert delta <= c*c*(1-2*c)/(1-c)**2
            else:
                assert equalized <= sign_count_bound(p,k)
            if n != 3:
                assert delta <= dimension_bound(n)
    return dict(dimension=n, coordinate_denominator=D, direction_denominator=W,
                weight_directions=len(weights), point_direction_pairs=pairs,
                positive_deficits=positive, covered_extreme_checks=covered_extreme_checks,
                sign_counts=sorted(sign_counts), grid_maximum=best, maximizer=maximizer,
                n3_support_comparisons=three_support_comparisons,
                scope="Finite falsification of continuous-support bounds; not a universal proof.")


def attaining_witness(p, k):
    n = p+k
    alpha, beta = Q(p-1,p), Q(k-1,k)
    delta = sign_count_bound(p,k)
    u, v = delta/alpha, 1-delta/beta
    assert 0 < u < Q(1,2) < v < 1 and v-u == 2*delta
    z = (u,)*p+(v,)*k
    w = (Q(1,2*p),)*p+(-Q(1,2*k),)*k
    T = tuple(tuple(Q(0) if j == i else 2*u if j < p else v-u
                    for j in range(n)) for i in range(p))
    T += tuple(tuple(u+1-v if j < p else Q(1) if j == i else 2*v-1
                     for j in range(n)) for i in range(p,n))
    assert all(0 <= t <= 1 for point in T for t in point)
    assert sum(abs(t) for t in w) == 1 and sum(w) == 0
    assert dot(w,z) == -delta and all(dot(w,t) == 0 for t in T)
    G = diagonal_generators(n,(Q(0),Q(1)))
    cert = reconstruction(G,T,z)
    expected_coefficients = (-v, -(1-u), Q(0))+(-u,)*p+(-(1-v),)*k
    assert cert["coefficients"] == expected_coefficients
    assert delta <= dimension_bound(n)
    if abs(p-k) <= 1:
        assert delta == dimension_bound(n)
    return dict(dimension=n, generator_count=3, positive_count=p, negative_count=k,
                balanced=abs(p-k)<=1, generators=G, target=z, context=T,
                affine_coefficients=w, minimum_concrete_union=Q(0),
                target_value=-delta, certified_gap=delta, dimension_upper_bound=dimension_bound(n),
                membership_certificate=cert)


def run():
    grids = [finite_checks(*spec) for spec in
             ((1,4,8),(2,4,8),(3,4,8),(4,3,8),(5,3,6),(6,2,4))]
    groups = {(p,k) for p in range(2,6) for k in range(2,6)}
    groups.update((n//2,n-n//2) for n in (*range(4,13),16,32))
    witnesses = [attaining_witness(p,k) for p,k in sorted(groups)]
    four = next(w for w in witnesses if w["dimension"] == 4)
    K = ((Q(0),)*4,(Q(1),)*4,four["target"])+four["context"]
    finite = finite_case("n4_attaining_affine_loss", signed_matrix(4), four["generators"],
                         K,[four["affine_coefficients"]])
    assert finite["contexts_enumerated"] == 128
    assert finite["direct_observation_losses"][0] == Q(1,6)
    assert finite["sector_observation_losses"][0] == Q(1,6)
    assert sum(g["covered_extreme_checks"] for g in grids) > 0
    return dict(complete=True, schema="affine-dimension-controls-v1",
                arithmetic="exact integers and fractions.Fraction", grids=grids,
                continuous_diagonal_certificates=witnesses, finite_context_check=finite,
                summary=dict(point_direction_pairs=sum(g["point_direction_pairs"] for g in grids),
                             positive_deficits=sum(g["positive_deficits"] for g in grids),
                             covered_extreme_checks=sum(g["covered_extreme_checks"] for g in grids),
                             rational_witnesses=len(witnesses), maximum_witness_dimension=32,
                             finite_contexts=128, discrepancies=0),
                scope="Bounded falsification checks and rational attainment; global upper bound is written mathematics.")


def main():
    output = ROOT/"results/theory/affine_dimensions.json"
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text('{"complete": false}\n',encoding="utf-8",newline="\n")
    output.with_suffix(".toml").write_text("complete = false\n",encoding="utf-8",newline="\n")
    report = run()
    sources = [Path(__file__), Path(__file__).with_name("checks.py"), Path(__file__).with_name("affine_profiles.py")]
    report["source_sha256"] = {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in sources}
    export_witnesses(report,output.with_suffix(".toml"))
    output.write_text(json.dumps(jsonable(report),indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps(report["summary"],indent=2))


if __name__ == "__main__":
    main()
