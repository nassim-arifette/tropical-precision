"""Exact controls for interior affine cell costs, not a numerical limit proof.

Finite point/direction grids try to falsify the universal upper bounds. Rational
certificates attain the proved cell costs. Other cells in the complete diagonal
representative can have larger loss: these are pointwise cell certificates.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import lcm
from pathlib import Path
import hashlib
import json

from checks import ROOT, diagonal_generators, dot, jsonable, reconstruction
from affine_dimensions import sector_support, weight_shell


def interior(a, b):
    return b-a <= min(a, 1-b)


def interior_cost(n, a, b):
    assert 0 <= a <= b <= 1 and interior(a, b)
    if n <= 2:
        return Q(0)
    h = b-a
    return h*max(a, 1-b)/(2-h) if n == 3 else h/2


def global_upper_three(a, b):
    h = b-a
    return h*max(1-a, b)/(2+h)


def point_deficit(x, w):
    D = lcm(*(v.denominator for v in x))
    W = lcm(*(v.denominator for v in w))
    xx, ww = tuple(int(v*D) for v in x), tuple(int(v*W) for v in w)
    values = [min(0, sum(ww))*D-sum(u*v for u, v in zip(xx, ww))]
    values += [sector_support(xx, ww, j, 1, D) for j in range(len(x)) if xx[j] > min(xx)]
    values += [sector_support(xx, ww, j, -1, D) for j in range(len(x)) if xx[j] < max(xx)]
    return Q(min(values), D*W)


def grid(n, a, b, W):
    points = tuple(a+(b-a)*Q(i, 3) for i in range(4))
    D = lcm(*(v.denominator for v in points))
    coordinates = tuple(int(v*D) for v in points)
    weights = tuple(weight_shell(n, W))
    bound = interior_cost(n, a, b) if interior(a, b) else global_upper_three(a, b)
    pairs = positive = three_levels = unequal_weights = 0
    best = Q(0)
    for x in combinations_with_replacement(coordinates, n):
        missing = [(j, 1) for j in range(n) if x[j] > x[0]]
        missing += [(j, -1) for j in range(n) if x[j] < x[-1]]
        for w in weights:
            pairs += 1
            values = [min(0, sum(w))*D-sum(u*v for u, v in zip(x, w))]
            values += [sector_support(x, w, j, sign, D) for j, sign in missing]
            numerator = min(values)
            if numerator <= 0:
                continue
            delta = Q(numerator, D*W)
            positive += 1
            three_levels += int(len(set(x)) >= 3)
            unequal_weights += int(any(len({abs(v) for v in w if sign*v > 0}) > 1
                                       for sign in (-1, 1)))
            assert delta <= bound, (n, a, b, x, w, delta, bound)
            if n == 3:
                assert delta <= global_upper_three(a, b)
            best = max(best, delta)
    return dict(dimension=n, interval=(a,b), coordinate_denominator=D,
                direction_denominator=W, point_direction_pairs=pairs,
                positive_deficits=positive, positive_three_level_points=three_levels,
                positive_unequal_sign_weights=unequal_weights,
                interior_formula_applicable=interior(a,b), checked_upper_bound=bound,
                grid_maximum=best)


def witness(n, a, b):
    assert n >= 3 and 0 < b-a <= min(a,1-b)
    reflected = n == 3 and a > 1-b
    aa, bb = (1-b, 1-a) if reflected else (a,b)
    h = bb-aa
    delta = interior_cost(n, a, b)
    if n == 3:
        A, c = 1/(2-h), (1-h)/(2-h)
        s, r = delta/c, 1-bb
        x = (aa,aa,bb)
        w = (A/2,A/2,-c)
        T = ((aa+r,aa+r,bb+r),(aa-s,aa+s,bb-s),(aa+s,aa-s,bb-s))
    else:
        p, k = n-2, 2
        alpha, beta = Q(p-1,p), Q(k-1,k)
        s, r = delta/alpha, delta/beta
        x = (aa,)*p+(bb,)*k
        w = (Q(1,2*p),)*p+(-Q(1,2*k),)*k
        T = tuple(tuple(aa-s if j == i else aa+s if j < p else bb-s
                        for j in range(n)) for i in range(p))
        T += tuple(tuple(aa+r if j < p else bb+r if j == i else bb-r
                         for j in range(n)) for i in range(p,n))
    if reflected:
        x, T, w = tuple(1-v for v in x), tuple(tuple(1-v for v in t) for t in T), tuple(-v for v in w)
    assert all(0 <= v <= 1 for t in T for v in t)
    assert all(a <= v <= b for v in x)
    assert sum(abs(v) for v in w) == 1
    baseline = min(Q(0),sum(w))
    assert all(dot(w,t) == baseline for t in T)
    assert baseline-dot(w,x) == delta == point_deficit(x,w)
    G = diagonal_generators(n,tuple(sorted({Q(0),a,b,Q(1)})))
    cert = reconstruction(G,T,x)
    return dict(dimension=n, interval=(a,b), reflected=reflected,
                generator_count=len(G), generators=G, target=x, context=T,
                affine_coefficients=w, minimum_concrete_union=baseline,
                target_value=dot(w,x), certified_gap=delta,
                coefficients=cert['coefficients'])


def run():
    intervals = ((Q(0),Q(1)), (Q(0),Q(1,4)), (Q(0),Q(1,2)),
                 (Q(1,4),Q(1,2)), (Q(3,8),Q(5,8)), (Q(1,2),Q(3,4)),
                 (Q(1,10),Q(1,5)), (Q(49,100),Q(51,100)), (Q(3,4),Q(1)))
    grids = [grid(3,a,b,14) for a,b in intervals]
    grids += [grid(4,a,b,8) for a,b in (intervals[3],intervals[4],intervals[7])]
    cells = tuple(pair for pair in intervals if interior(*pair)) + ((Q(4,5),Q(9,10)),)
    witnesses = [witness(n,a,b) for n in (3,4,5,8,16) for a,b in cells]
    assert interior_cost(3,Q(1,4),Q(1,2)) == Q(1,14) > Q(1,16)
    assert interior_cost(3,Q(3,8),Q(5,8)) == Q(3,56) > Q(1,24)
    assert sum(g['positive_three_level_points'] for g in grids) > 0
    assert sum(g['positive_unequal_sign_weights'] for g in grids) > 0
    return dict(complete=True, schema='affine-cell-controls-v1',
                arithmetic='exact integers and fractions.Fraction', grids=grids,
                witnesses=witnesses,
                summary=dict(grids=len(grids), point_direction_pairs=sum(g['point_direction_pairs'] for g in grids),
                             positive_deficits=sum(g['positive_deficits'] for g in grids),
                             positive_three_level_points=sum(g['positive_three_level_points'] for g in grids),
                             positive_unequal_sign_weights=sum(g['positive_unequal_sign_weights'] for g in grids),
                             rational_witnesses=len(witnesses), discrepancies=0),
                scope='Finite upper-bound falsification and interior-cell attainers; the asymptotic constants have a written proof, not a numerical fit.')


def main():
    output = ROOT/'results/theory/affine_cells.json'
    interchange = output.with_suffix('.toml')
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text('{"complete": false}\n',encoding='utf-8',newline='\n')
    interchange.write_text('complete = false\n',encoding='utf-8',newline='\n')
    report = run()
    sources = [Path(__file__), Path(__file__).with_name('affine_dimensions.py'), Path(__file__).with_name('affine_profiles.py'), Path(__file__).with_name('checks.py')]
    report['source_sha256'] = {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    lines = ['schema = "affine-cell-controls-v1"','complete = true']
    for record in report['witnesses']:
        lines += ['', '[[witnesses]]']
        lines += [f'{k} = {json.dumps(jsonable(v))}' for k,v in record.items()]
    interchange.write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
    output.write_text(json.dumps(jsonable(report),indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(report['summary'],indent=2))


if __name__ == '__main__':
    main()
