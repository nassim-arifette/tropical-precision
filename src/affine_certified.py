"""Certified global affine cell and finite-budget bounds, conditional on D8/D9.

Floating arithmetic proposes knots only. Exact branch covers and rational
context witnesses establish the reported bounds; no optimizer success flag is
scientific evidence. Only standard libraries are needed.
"""
from fractions import Fraction as Q
from functools import lru_cache
from heapq import heappush, heappop
from math import sqrt
from itertools import combinations_with_replacement
from random import Random
from pathlib import Path
import hashlib
import json

from checks import ROOT, diagonal_generators, dot, jsonable, reconstruction
from affine_cells import interior, interior_cost, point_deficit, witness as interior_witness
from affine_dimensions import sector_support
import affine_two_level as two


def branches(n, a, b):
    result, constants = [], []
    for reflected, aa, bb in ((False, a, min(b, Q(1, 2))),
                              (True, 1-b, min(1-a, Q(1, 2)))):
        if aa >= bb:
            continue
        for kind in ('positive', 'negative'):
            result.append(dict(kind=kind, reflected=reflected, a=aa, b=bb,
                               lo=Q(0), hi=Q(1, 2)))
        if n >= 4:
            delta, x, w = two.balanced_left(n, aa, bb)
            if reflected:
                x, w = tuple(1-v for v in x), tuple(-v for v in w)
            constants.append((delta, x, w))
    if a < Q(1, 2) < b:
        for k in range(1, n):
            p = n-k
            result.append(dict(kind='central', reflected=False, a=a, b=b, k=k,
                               lo=max(Q(0), 1-Q(k, 2)), hi=min(Q(1), Q(p, 2))))
    return result, constants


def branch_at(n, spec, t):
    if t == spec['lo'] or t == spec['hi']:
        return Q(0), None, None
    if spec['kind'] == 'positive':
        result = two.positive_singleton_left(n, spec['a'], spec['b'], t)
    elif spec['kind'] == 'negative':
        result = two.negative_singleton_left(n, spec['a'], spec['b'], t)
    else:
        result = two.central(n, spec['a'], spec['b'], t, spec['k'])
    delta, x, w = result
    if spec['reflected']:
        x, w = tuple(1-v for v in x), tuple(-v for v in w)
    return delta, x, w


def context_certificate(n, a, b, x, w, delta):
    assert delta > 0 and point_deficit(x, w) == delta
    assert len(x) == len(w) == n and all(a <= v <= b for v in x)
    missing = [(j, 1) for j in range(n) if x[j] > min(x)]
    missing += [(j, -1) for j in range(n) if x[j] < max(x)]
    T = tuple(dict.fromkeys(two.sector_point(x, w, j, sign) for j, sign in missing))
    G = diagonal_generators(n, tuple(sorted({Q(0), a, b, Q(1)})))
    concrete = min([min(Q(0), sum(w))]+[dot(w, t) for t in T])
    assert concrete-dot(w, x) == delta
    assert all(0 <= v <= 1 for t in T for v in t)
    cert = reconstruction(G, T, x)
    return dict(dimension=n, interval=(a, b), generators=G, generator_count=len(G),
                target=x, affine_coefficients=w, context=T, coefficients=cert['coefficients'],
                minimum_concrete_union=concrete, target_value=dot(w, x), certified_gap=delta)


def certify_cell(n, a, b, tolerance=Q(1, 100000), *, shortcuts=True):
    assert isinstance(a, Q) and isinstance(b, Q) and isinstance(tolerance, Q)
    assert n >= 1 and 0 <= a <= b <= 1 and tolerance > 0
    if n <= 2 or a == b:
        return dict(dimension=n, interval=(a, b), method='zero', lower=Q(0), upper=Q(0))
    if shortcuts and interior(a, b):
        cost = interior_cost(n, a, b)
        return dict(dimension=n, interval=(a, b), method='interior', lower=cost, upper=cost,
                    witness=interior_witness(n, a, b))
    specs, constants = branches(n, a, b)
    best, best_x, best_w = max(constants, key=lambda row: row[0], default=(Q(0), None, None))
    heap, serial = [], 0
    evaluations = 0

    def insert(index, lo, hi):
        nonlocal best, best_x, best_w, serial, evaluations
        mid = (lo+hi)/2
        value, x, w = branch_at(n, specs[index], mid)
        evaluations += 1
        if value > best:
            best, best_x, best_w = value, x, w
        # D9's 2-Lipschitz bound: radius (hi-lo)/2 times constant 2.
        heappush(heap, (-(value+hi-lo), serial, index, lo, hi, value))
        serial += 1

    for index, spec in enumerate(specs):
        insert(index, spec['lo'], spec['hi'])
    while -heap[0][0] > best+tolerance:
        _, _, index, lo, hi, _ = heappop(heap)
        mid = (lo+hi)/2
        insert(index, lo, mid)
        insert(index, mid, hi)
    upper = max(best, -heap[0][0])
    leaves = [dict(branch=i, lo=lo, hi=hi, midpoint_value=value)
              for _, _, i, lo, hi, value in sorted(heap, key=lambda row: (row[2], row[3]))]
    result = dict(dimension=n, interval=(a, b), method='scalar', lower=best, upper=upper,
                  tolerance=tolerance, evaluations=evaluations, leaves=leaves)
    if best > 0:
        result['witness'] = context_certificate(n, a, b, best_x, best_w, best)
    assert 0 <= upper-best <= tolerance
    return result


def approximate_maximum(f, lo, hi):
    """Heuristic proposal only: a coarse scan followed by local golden search."""
    step = (hi-lo)/32
    nodes = [lo+i*step for i in range(33)]
    values = [f(x) for x in nodes]
    index = max(range(33), key=values.__getitem__)
    left, right = nodes[max(0, index-1)], nodes[min(32, index+1)]
    ratio = (sqrt(5)-1)/2
    c, d = right-ratio*(right-left), left+ratio*(right-left)
    fc, fd = f(c), f(d)
    for _ in range(32):
        if fc > fd:
            right, d, fd = d, c, fc
            c = right-ratio*(right-left); fc = f(c)
        else:
            left, c, fc = c, d, fd
            d = left+ratio*(right-left); fd = f(d)
    return max(max(values), fc, fd)


def approximate_left(n, a, b):
    h = b-a
    s = min(h, 2-sqrt(4-2*b))
    positive = s*(b-s)/(2-s)

    def negative(c):
        if c <= 0 or c >= .5:
            return 0.
        A, gamma = 1-c, 1-2*c
        alpha = 1-2*A/(n-1)
        return max(0., min(c*b-A*a, alpha*c*b/(alpha+A),
                           gamma*(c/A-a), c/(1+A/alpha+c/gamma)))

    cost = max(positive, approximate_maximum(negative, 0., .5))
    if n >= 4:
        cost = max(cost, min(h/2, (n-3)*b/(3*n-8)))
    return cost


@lru_cache(maxsize=100000)
def approximate_cell(n, a, b):
    if n <= 2 or a == b:
        return 0.
    h = b-a
    if h <= min(a, 1-b):
        return h*max(a, 1-b)/(2-h) if n == 3 else h/2
    costs = []
    if a <= .5:
        costs.append(approximate_left(n, a, min(b, .5)))
    if b >= .5:
        costs.append(approximate_left(n, 1-b, min(1-a, .5)))
    if a <= .5 <= b:
        for k in range(1, n):
            p = n-k
            lo, hi = max(0., 1-k/2), min(1., p/2)

            def central(A):
                if A <= lo or A >= hi:
                    return 0.
                C = 1-A
                alpha, beta, m = 1-2*A/p, 1-2*C/k, min(A, C)
                return max(0., min(alpha/2, beta/2, m-A*a-C*(1-b),
                                   (m-C*(1-b))/(1+A/alpha),
                                   (m-A*a)/(1+C/beta), m/(1+A/alpha+C/beta)))

            costs.append(approximate_maximum(central, lo, hi))
    return max(costs)


def propose_partition(n, q, denominator=10**9):
    if q == 1:
        return (Q(0), Q(1))
    if q == 2:
        return (Q(0), Q(1, 2), Q(1))  # D9's rigorous symmetry reduction.

    def greedy(delta):
        knots = [0.]
        for _ in range(q):
            a = knots[-1]
            if approximate_cell(n, a, 1.) <= delta:
                knots.append(1.)
                continue
            lo, hi = a, 1.
            for _ in range(27):
                mid = (lo+hi)/2
                if approximate_cell(n, a, mid) <= delta:
                    lo = mid
                else:
                    hi = mid
            knots.append(lo)
        return knots

    lo, hi = 0., .25
    for _ in range(27):
        mid = (lo+hi)/2
        if greedy(mid)[-1] == 1.:
            hi = mid
        else:
            lo = mid
    knots = greedy(hi)
    rational = tuple(Q(round(v*denominator), denominator) for v in knots)
    assert rational[0] == 0 and rational[-1] == 1
    assert all(a < b for a, b in zip(rational, rational[1:])), 'Proposal has empty cells'
    return rational


def balanced_partition(n,q):
    assert n>=4 and q>=2
    delta=1/(2*(Q(q+1)+Q(1,n-3)))
    boundary=delta*Q(3*n-8,n-3)
    return (Q(0),)+tuple(boundary+2*i*delta for i in range(q-1))+(Q(1),)


def certify_partition(n, partition, tolerance=Q(1, 100000)):
    assert partition[0] == 0 and partition[-1] == 1
    assert all(a < b for a, b in zip(partition, partition[1:]))
    cells = [certify_cell(n, a, b, tolerance) for a, b in zip(partition, partition[1:])]
    return dict(dimension=n, generator_budget=len(partition)+1, partition=partition,
                lower=min(c['lower'] for c in cells), upper=max(c['upper'] for c in cells),
                cells=cells, certificate_basis='D9 min-cell lower / max-cell upper bounds')


def m4_algebra_checks():
    # Rational polynomial identities used by the written irrational proof.
    for c in (Q(i, 11) for i in range(-11, 12)):
        assert (1-2*c)-(1-c)**2/2 == -(c*c+2*c-1)/2
        assert (2*c-6*c*c)*(1-c)+2*(c*c-2*c**3) == 2*c*(1-3*c+c*c)
        assert (1-3*c+c*c)-(2-5*c) == c*c+2*c-1
    # Isolate the positive root c*=sqrt(2)-1 and hence its squared loss.
    lo, hi = Q(4142135623730950, 10**16), Q(4142135623730951, 10**16)
    assert lo*lo+2*lo-1 < 0 < hi*hi+2*hi-1
    assert Q(2, 5) < lo < hi < Q(1, 2)
    assert Q(1, 24) < lo*lo/2
    return dict(identities=3, evaluations=69, root_interval=(lo, hi),
                loss_interval=(lo*lo/2, hi*hi/2))


def structural_checks():
    rng=Random(9026)
    supports=parameters=containments=0
    for n in (1,2,3,4,8):
        for _ in range(16):
            x=tuple(rng.randrange(13) for _ in range(n))
            y=tuple(rng.randrange(13) for _ in range(n))
            w=tuple(rng.randrange(-4,5) for _ in range(n))
            norm=sum(abs(v) for v in w)
            eta=max(abs(a-b) for a,b in zip(x,y))
            for j in range(n):
                for sign in (-1,1):
                    assert abs(sector_support(x,w,j,sign,12)-sector_support(y,w,j,sign,12))<=eta*norm
                    supports+=1
    for n in (3,4,8):
        for a,b in ((Q(0),Q(1,10)),(Q(1,10),Q(3,10)),(Q(1,4),Q(3,4)),(Q(0),Q(1))):
            specs,_=branches(n,a,b)
            for spec in specs:
                nodes=[spec['lo']+(spec['hi']-spec['lo'])*Q(i,16) for i in range(17)]
                for left,right in zip(nodes,nodes[1:]):
                    assert abs(branch_at(n,spec,left)[0]-branch_at(n,spec,right)[0])<=2*(right-left)
                    parameters+=1
    for q in range(1,5):
        partitions=[(0,)+knots+(4,) for knots in combinations_with_replacement(range(5),q-1)]
        for s in partitions:
            for t in partitions:
                assert any(t[i]<=s[i]<=s[i+1]<=t[i+1] for i in range(q))
                containments+=1
    return dict(support_lipschitz_pairs=supports,scalar_lipschitz_pairs=parameters,
                weak_partition_containments=containments,
                scope='Finite implementation controls; universal continuity and duality are written proofs.')


def run():
    algebra = m4_algebra_checks()
    structural = structural_checks()
    controls = []
    specifications = [(3,Q(0),Q(1)), (3,Q(0),Q(1,2)), (4,Q(0),Q(1)),
                      (4,Q(0),Q(1,2)), (8,Q(0),Q(1,2)),
                      (3,Q(1,4),Q(1,2)), (4,Q(1,4),Q(1,2)),
                      (4,Q(0),Q(1,10)), (4,Q(9,10),Q(1)),
                      (3,Q(0),Q(1,100)), (3,Q(1,3),Q(1,3)), (2,Q(0),Q(1))]
    for n,a,b in specifications:
        cell=certify_cell(n,a,b,Q(1,10000000),shortcuts=False)
        if interior(a,b) or n <= 2:
            expected=interior_cost(n,a,b) if n >= 3 else Q(0)
            assert cell['lower'] <= expected <= cell['upper']
        if n >= 4 and (a,b) == (Q(0),Q(1,2)):
            expected=Q(n-3,2*(3*n-8))
            assert cell['lower'] <= expected <= cell['upper']
        if n == 4 and (a,b) == (Q(0),Q(1)):
            assert cell['lower'] <= Q(1,6) <= cell['upper']
        if n == 3 and (a,b) == (Q(0),Q(1,2)):
            left,right=algebra['loss_interval']
            assert cell['lower'] <= right and left <= cell['upper']
        controls.append(cell)
        print(json.dumps(jsonable(dict(control=(n,a,b), lower=cell['lower'], upper=cell['upper'],
                                       evaluations=cell.get('evaluations',0)))),flush=True)
    records=[]
    for n,q in ((3,2),(3,3),(3,4),(3,6),(4,2),(4,3),(4,4),(4,6),(4,7),(4,8),(4,12),(8,3)):
        partition=(balanced_partition(n,q) if (n==4 and q<=6) or (n==8 and q==3)
                   else propose_partition(n,q))
        record=certify_partition(n,partition,Q(1,10000000))
        assert record['upper']-record['lower'] <= Q(3,10000000)
        if n==4 and q<=6:
            assert record['lower']==record['upper']==Q(1,2*(q+2))
        if n==4 and q>=7:
            assert record['lower']>Q(1,2*(q+2)), 'Must distinguish the failed balanced-law extrapolation'
        if n==8:
            assert record['lower']==record['upper']==Q(5,42)
        records.append(record)
        print(json.dumps(jsonable({k:record[k] for k in ('dimension','generator_budget','partition','lower','upper')})),flush=True)
    return dict(schema='affine-certified-minimax-v1', complete=True, algebra=algebra, structural_checks=structural,
                controls=controls, partitions=records,
                scope='Global cell and finite-budget intervals conditional on D8/D9; floating search only proposes rational knots.',
                summary=dict(cell_controls=len(controls), partitions=len(records),
                             maximum_global_interval_width=max(r['upper']-r['lower'] for r in records),
                             scalar_leaves=sum(len(c.get('leaves',[])) for c in controls)
                             +sum(len(c.get('leaves',[])) for r in records for c in r['cells'])))


def main():
    output=ROOT/'results/theory/affine_certified.json'
    interchange=output.with_suffix('.toml')
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text('{"complete": false}\n',encoding='utf-8',newline='\n')
    interchange.write_text('complete = false\n',encoding='utf-8',newline='\n')
    report=run()
    sources=[Path(__file__)]+[Path(__file__).with_name(x) for x in
            ('checks.py','affine_profiles.py','affine_dimensions.py','affine_cells.py','affine_two_level.py')]

    report['source_sha256']={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    # Julia reads explicit tables with rational strings using its standard
    # TOML parser; no JSON package or floating round-trip is required.
    lines=['schema = "affine-certified-minimax-v1"','complete = true']
    for cell in report['controls']:
        write_cell(lines,cell,'controls')
    for record in report['partitions']:
        lines += ['', '[[partitions]]']
        for key in ('dimension','generator_budget','partition','lower','upper'):
            lines.append(f'{key} = {json.dumps(jsonable(record[key]))}')
        for cell in record['cells']:
            write_cell(lines,cell,'partitions.cells')
    lines += ['', '[algebra]']
    lines += [f'{key} = {json.dumps(jsonable(value))}' for key,value in report['algebra'].items()]
    interchange.write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
    output.write_text(json.dumps(jsonable(report),indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(jsonable(report['summary']),indent=2))


def write_cell(lines,cell,table):
    lines += ['', '[['+table+']]']
    for key,value in cell.items():
        if key not in ('leaves','witness'):
            lines.append(f'{key} = {json.dumps(jsonable(value))}')
    if 'witness' in cell:
        lines += ['['+table+'.witness]']
        lines += [f'{key} = {json.dumps(jsonable(value))}' for key,value in cell['witness'].items()]
    for leaf in cell.get('leaves',[]):
        lines += ['[['+table+'.leaves]]']
        lines += [f'{key} = {json.dumps(jsonable(value))}' for key,value in leaf.items()]


if __name__ == '__main__':
    main()
