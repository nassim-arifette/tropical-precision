"""Constructive exact falsification of the all-cell two-level reduction.

Every positive tested configuration is mapped to an explicit two-level one
with at least its deficit. This does not globally optimize the scalar branches.
Saved witnesses are checked separately in Julia, including original deficits.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import lcm
from pathlib import Path
import hashlib
import json

from checks import ROOT, diagonal_generators, dot, jsonable, reconstruction
from affine_dimensions import sector_support, weight_shell
from affine_cells import point_deficit, interior, interior_cost

HALF = Q(1, 2)


def balanced_left(n, a, b):
    assert n >= 4 and 0 <= a <= b <= HALF
    p, k = n-2, 2
    alpha = Q(p-1, p)
    delta = min((b-a)/2, alpha*b/(1+2*alpha))
    u, v = max(a, b/(1+2*alpha)), b
    return delta, (u,)*p+(v,)*k, (Q(1, 2*p),)*p+(-Q(1, 4),)*k


def positive_singleton_left(n, a, b, A):
    assert n >= 3 and 0 <= a <= b <= HALF and 0 < A < HALF
    C, alpha = 1-A, 1-2*A
    B = A-C*(1-b)
    delta = max(Q(0), min(B-A*a, alpha*B/(A+alpha)))
    u, v = max(a, B/C), b
    return delta, (u,)+(v,)*(n-1), (A,)+(-C/(n-1),)*(n-1)


def negative_singleton_left(n, a, b, c):
    assert n >= 3 and 0 <= a <= b <= HALF and 0 < c < HALF
    p, A, gamma = n-1, 1-c, 1-2*c
    alpha = 1-2*A/p
    delta = max(Q(0), min(c*b-A*a, alpha*c*b/(alpha+A),
                          gamma*(c/A-a), c/(1+A/alpha+c/gamma)))
    u, v = max(a, delta/alpha), min(b, 1-delta/gamma)
    return delta, (u,)*p+(v,), (A/p,)*p+(-c,)


def central(n, a, b, A, k):
    p, C = n-k, 1-A
    assert p >= 1 and k >= 1 and a <= HALF <= b
    alpha, beta, m = 1-2*A/p, 1-2*C/k, min(A, C)
    assert 0 < A < 1 and alpha > 0 and beta > 0
    delta = max(Q(0), min(alpha/2, beta/2, m-A*a-C*(1-b),
                          (m-C*(1-b))/(1+A/alpha),
                          (m-A*a)/(1+C/beta),
                          m/(1+A/alpha+C/beta)))
    u, v = max(a, delta/alpha), 1-max(1-b, delta/beta)
    return delta, (u,)*p+(v,)*k, (A/p,)*p+(-C/k,)*k


def reduce_positive(x, w, a, b):
    n = len(x)
    assert n == len(w) and n >= 3 and all(a <= v <= b for v in x)
    norm = sum(abs(v) for v in w)
    assert 0 < norm <= 1
    normalized = tuple(v/norm for v in w)
    A = sum(v for v in normalized if v > 0)
    C = -sum(v for v in normalized if v < 0)
    p, k = sum(v > 0 for v in w), sum(v < 0 for v in w)
    assert A > 0 and C > 0 and all(abs(v) < HALF for v in normalized)
    u = sum(v*t for v, t in zip(normalized, x) if v > 0)/A
    v = -sum(v*t for v, t in zip(normalized, x) if v < 0)/C
    assert a <= u < v <= b
    reflected = u >= HALF and v > HALF
    if v <= HALF or reflected:
        aa, bb = (1-b, min(1-a, HALF)) if reflected else (a, min(b, HALF))
        pp, kk, AA, CC = (k, p, C, A) if reflected else (p, k, A, C)
        side = 'right_' if reflected else 'left_'
        if pp == 1:
            branch = side+'positive_singleton'
            delta, y, z = positive_singleton_left(n, aa, bb, AA)
        elif kk == 1:
            branch = side+'negative_singleton'
            delta, y, z = negative_singleton_left(n, aa, bb, CC)
        else:
            branch = side+'both'
            delta, y, z = balanced_left(n, aa, bb)
        if reflected:
            y, z = tuple(1-t for t in y), tuple(-t for t in z)
    else:
        branch = 'central'
        delta, y, z = central(n, a, b, A, k)
    assert delta > 0 and all(a <= t <= b for t in y)
    assert sum(abs(t) for t in z) == 1
    assert len(set(y)) == 2
    assert all(len({z[j] for j in range(n) if sign*z[j] > 0}) == 1 for sign in (-1, 1))
    assert all(y[j] == (min(y) if z[j] > 0 else max(y)) for j in range(n))
    return branch, delta, y, z


def sector_point(x, w, k, sign):
    """Explicit continuous-support maximizer, not a grid of context points."""
    radius = 1-x[k] if sign == 1 else x[k]
    bounds = tuple(1-t if v >= 0 else t for t, v in zip(x, w))
    nodes = {Q(0), radius}
    nodes.update(min(radius, bounds[j]) for j in range(len(x)) if j != k)
    candidates = []
    for s in sorted(nodes):
        h = tuple(sign*s if j == k else min(s, bounds[j]) if w[j] > 0
                  else -min(s, bounds[j]) if w[j] < 0 else Q(0) for j in range(len(x)))
        candidates.append((dot(w, h), tuple(t+v for t, v in zip(x, h))))
    return max(candidates, key=lambda pair: pair[0])[1]


def certificate(x, w, a, b):
    original = point_deficit(x, w)
    assert original > 0
    branch, delta, y, z = reduce_positive(x, w, a, b)
    assert point_deficit(y, z) == delta >= original
    missing = [(j, 1) for j in range(len(y)) if y[j] > min(y)]
    missing += [(j, -1) for j in range(len(y)) if y[j] < max(y)]
    T = tuple(dict.fromkeys(sector_point(y, z, j, sign) for j, sign in missing))
    G = diagonal_generators(len(x), tuple(sorted({Q(0), a, b, Q(1)})))
    minimum = min([min(Q(0), sum(z))]+[dot(z, t) for t in T])
    assert minimum-dot(z, y) == delta
    cert = reconstruction(G, T, y)
    return dict(dimension=len(x), interval=(a, b), branch=branch,
                original_target=x, original_weights=w, original_deficit=original,
                target=y, affine_coefficients=z, certified_gap=delta,
                generators=G, generator_count=len(G), context=T,
                minimum_concrete_union=minimum, target_value=dot(z, y),
                coefficients=cert['coefficients'])


def grid(n, a, b, W):
    nodes = tuple(a+(b-a)*Q(i, 3) for i in range(4))
    D = lcm(*(v.denominator for v in nodes))
    coordinates = tuple(int(v*D) for v in nodes)
    weights = tuple(weight_shell(n, W))
    counts = Counter()
    representatives = {}
    for xx in combinations_with_replacement(coordinates, n):
        missing = [(j, 1) for j in range(n) if xx[j] > xx[0]]
        missing += [(j, -1) for j in range(n) if xx[j] < xx[-1]]
        for ww in weights:
            counts['point_direction_pairs'] += 1
            numerator = min(0, sum(ww))*D-sum(t*v for t, v in zip(xx, ww))
            if numerator <= 0:
                continue
            numerator = min([numerator]+[sector_support(xx, ww, j, sign, D)
                                        for j, sign in missing])
            if numerator <= 0:
                continue
            x, w = tuple(Q(t, D) for t in xx), tuple(Q(v, W) for v in ww)
            original = Q(numerator, D*W)
            branch, bound, y, z = reduce_positive(x, w, a, b)
            actual = point_deficit(y, z)
            assert original <= actual == bound, (n, a, b, x, w, original, branch, bound, actual)
            if interior(a, b):
                assert bound <= interior_cost(n, a, b)
            counts['positive_deficits'] += 1
            counts['strict_improvements'] += int(bound > original)
            levels = len(set(x))
            unequal = any(len({abs(v) for v in w if sign*v > 0}) > 1 for sign in (-1, 1))
            zero = any(v == 0 for v in w)
            covered = any((v > 0 and t == max(x)) or (v < 0 and t == min(x)) for t, v in zip(x, w))
            counts['three_or_more_levels'] += int(levels >= 3)
            counts['unequal_sign_weights'] += int(unequal)
            counts['zero_coefficients'] += int(zero)
            counts['covered_extrema'] += int(covered)
            counts[branch] += 1
            key = (branch, levels >= 3, unequal, zero, covered)
            representatives.setdefault(key, (x, w, a, b))
    return (dict(dimension=n, interval=(a, b), coordinate_denominator=D,
                 direction_denominator=W, counts=dict(sorted(counts.items()))),
            list(representatives.values()))


def direct_branch_checks():
    """Test scalar feasibility/attainment at rational parameters, not optimality."""
    counts = Counter()
    specimens = []
    for n in (3, 4, 5, 8):
        for a, b in ((Q(0), Q(1, 10)), (Q(0), HALF), (Q(1, 4), HALF),
                     (Q(1, 10), Q(1, 5)), (Q(49, 100), HALF)):
            cases = [('positive', *positive_singleton_left(n, a, b, Q(j, 40))) for j in range(1, 20)]
            cases += [('negative', *negative_singleton_left(n, a, b, Q(j, 40))) for j in range(1, 20)]
            if n >= 4:
                cases += [('both', *balanced_left(n, a, b))]
            for branch, delta, x, w in cases:
                if delta > 0:
                    assert all(a <= v <= b for v in x)
                    assert point_deficit(x, w) == delta, (n, a, b, branch, delta, x, w)
                    counts['left_'+branch] += 1
        for a, b in ((Q(0), Q(1)), (Q(0), Q(3, 4)), (Q(1, 4), Q(3, 4)),
                     (Q(2, 5), Q(9, 10)), (HALF, Q(1)), (Q(0), HALF)):
            for k in range(1, n):
                for j in range(1, 20):
                    A = Q(j, 20)
                    if 1-2*A/(n-k) <= 0 or 1-2*(1-A)/k <= 0:
                        continue
                    delta, x, w = central(n, a, b, A, k)
                    if delta > 0:
                        assert point_deficit(x, w) == delta
                        counts['central'] += 1
    delta, x, w = negative_singleton_left(4, Q(0), Q(1, 10), Q(47, 100))
    assert delta == Q(4559, 176500) > Q(1, 40)
    assert x == (Q(141, 3530),)*3+(Q(1, 10),)
    specimens.append((x, w, Q(0), Q(1, 10)))
    specimens.append((tuple(1-v for v in x), tuple(-v for v in w), Q(9, 10), Q(1)))
    return dict(counts), specimens


def run():
    specifications = [(3, a, b, 14) for a, b in
                      ((Q(0), Q(1)), (Q(0), Q(1, 4)), (Q(0), HALF),
                       (Q(1, 4), HALF), (Q(1, 4), Q(3, 4)),
                       (HALF, Q(1)), (Q(3, 4), Q(1)))]
    specifications += [(4, a, b, 10) for a, b in
                       ((Q(0), Q(1)), (Q(0), HALF), (Q(1, 4), Q(3, 4)), (HALF, Q(1)))]
    specifications += [(n, Q(0), Q(1), 8) for n in (5, 6)]
    grids, records = [], []
    for spec in specifications:
        result, representatives = grid(*spec)
        grids.append(result)
        records.extend(certificate(*v) for v in representatives)
        print(json.dumps(jsonable(dict(dimension=spec[0], interval=spec[1:3],
                                       counts=result['counts']))), flush=True)
    branches, explicit = direct_branch_checks()
    records.extend(certificate(*v) for v in explicit)
    # Homogeneity controls also cover directions in the interior of the l1 ball.
    records.extend(certificate(tuple(Q(v) for v in c['original_target']),
                               tuple(Q(v)/3 for v in c['original_weights']),
                               *c['interval']) for c in list(records[:7]))
    totals = Counter()
    for g in grids:
        totals.update(g['counts'])
    for field in ('three_or_more_levels', 'unequal_sign_weights', 'zero_coefficients', 'covered_extrema',
                  'left_both', 'right_both', 'left_positive_singleton', 'right_positive_singleton',
                  'left_negative_singleton', 'right_negative_singleton', 'central'):
        assert totals[field] > 0, field
    return dict(complete=True, schema='affine-two-level-controls-v1',
                arithmetic='exact integers and fractions.Fraction', grids=grids,
                scalar_branch_attainment_checks=branches, witnesses=records,
                summary={**dict(totals), 'rational_certificates': len(records), 'discrepancies': 0},
                scope='Constructive domination for each positive tested configuration and rational branch attainment; not a global scalar optimization or a universal proof.')


def main():
    output = ROOT/'results/theory/affine_two_level.json'
    interchange = output.with_suffix('.toml')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text('{"complete": false}\n', encoding='utf-8', newline='\n')
    interchange.write_text('complete = false\n', encoding='utf-8', newline='\n')
    report = run()
    sources = [Path(__file__), Path(__file__).with_name('affine_cells.py'), Path(__file__).with_name('affine_dimensions.py'), Path(__file__).with_name('affine_profiles.py'), Path(__file__).with_name('checks.py')]
    report['source_sha256'] = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in sources}
    lines = ['schema = "affine-two-level-controls-v1"', 'complete = true']
    for case in report['witnesses']:
        lines += ['', '[[witnesses]]']
        lines += [f'{key} = {json.dumps(jsonable(value))}' for key, value in case.items()]
    interchange.write_text('\n'.join(lines)+'\n', encoding='utf-8', newline='\n')
    output.write_text(json.dumps(jsonable(report), indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(report['summary'], indent=2))


if __name__ == '__main__':
    main()
