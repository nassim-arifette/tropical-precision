"""Exact checks of D10 boundary and all-budget laws; D9 certificates are controls.

Root isolation and algebra use rationals. Selected global partition controls
are certified independently of D10 by the older D8/D9 scalar-cover machinery.
Large-budget root checks are conditional on D10, not large geometric searches.
"""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib
import json

from checks import ROOT, jsonable
from affine_certified import certify_cell, certify_partition, write_cell


def B(d, c):
    return (1-2*c)*(2*(d+1)-d*c)/(2*(d+1)*(1-c)**2)


def F(d, c):
    return c*(d+2*c)*(1-2*c)/(2*(d+1)*(1-c)**2)


def P(d, q, c):
    return d+1+((q-5)*d-2)*c+(2*q-5)*(1-d)*c*c-4*(q-2)*c**3


def isolate(residual, lo, hi, *, increasing=False, steps=80):
    sign = 1 if increasing else -1
    assert sign*residual(lo) <= 0 <= sign*residual(hi)
    for _ in range(steps):
        mid = (lo+hi)/2
        value = sign*residual(mid)
        if value == 0:
            return mid, mid
        if value < 0:
            lo = mid
        else:
            hi = mid
    return lo, hi


def boundary_root(n, b):
    assert n >= 3 and 0 <= b <= Q(1, 2)
    d = n-3
    if not b:
        lo = hi = Q(1, 2)
    else:
        lo, hi = isolate(lambda c: B(d,c)-b, Q(2,5) if d == 0 else Q(3,8), Q(1,2))
    balanced = Q(d,3*d+1)*b
    return dict(dimension=n, endpoint=b, root_interval=(lo,hi),
                loss_interval=(max(balanced,F(d,hi)),max(balanced,F(d,lo))))


def budget_root(n, M):
    assert n >= 4 and M >= 4
    d,q = n-3,M-2
    lo,hi = isolate(lambda c: P(d,q,c),Q(3,8),Q(1,2))
    balanced = Q(1,2)/(q+1+Q(1,d))
    lower,upper = max(balanced,F(d,hi)),max(balanced,F(d,lo))
    return dict(dimension=n,generator_budget=M,root_interval=(lo,hi),
                balanced_value=balanced,loss_interval=(lower,upper),
                boundary_interval=(Q(1,2)-(q-2)*upper,Q(1,2)-(q-2)*lower),
                second_order_interval=(M*M*(lower-Q(1,2*M)),M*M*(upper-Q(1,2*M))))


def transition_root(n):
    d=n-3
    lo,hi=isolate(lambda c: 2*(3*d+1)*c*c+d*(4*d+1)*c-2*d*(d+1),
                  Q(0),Q(1,2),increasing=True)
    blo,bhi=B(d,hi),B(d,lo)
    kap=Q(d,3*d+1)
    mlo,mhi=4+(Q(1,2)/bhi-1)/kap,4+(Q(1,2)/blo-1)/kap
    assert mlo.numerator//mlo.denominator == mhi.numerator//mhi.denominator
    return dict(dimension=n,root_interval=(lo,hi),boundary_interval=(blo,bhi),
                budget_threshold_interval=(mlo,mhi),last_balanced_budget=mlo.numerator//mlo.denominator)


def algebra_checks():
    identities = [
        ('k_derivative',(2,2),lambda d,c:
         (d+4*c)*(2*(d+1)-d*c)+d*c*(d+2*c)
         -(2*d*(d+1)+8*(d+1)*c-2*d*c*c)),
        ('B_derivative',(1,2),lambda d,c:
         (-2*(2*(d+1)-d*c)-d*(1-2*c))*(1-c)
         +2*(1-2*c)*(2*(d+1)-d*c)+d+(d+4)*c),
        ('F_derivative',(1,3),lambda d,c:
         ((d+4*c)*(1-2*c)-2*c*(d+2*c))*(1-c)
         +2*c*(d+2*c)*(1-2*c)-d*(1-3*c)-4*c*(1-3*c+c*c)),
        ('cubic_expansion',(1,1,3),lambda d,q,c:
         (1-2*c)*(2*(d+1)+(q-3)*d*c+2*(q-2)*c*c)
         -(d+1)*(1-c)**2-P(d,q,c)),
        ('boundary_quadratic',(1,1,2),lambda d,b,c:
         2*(d+1)*b*(1-c)**2-(1-2*c)*(2*(d+1)-d*c)
         -(2*(b*(d+1)-d)*c*c+(5*d+4-4*b*(d+1))*c+2*(d+1)*(b-1))),
        ('switch_polynomial',(2,2),lambda d,c:
         (3*d+1)*c*(d+2*c)-d*(2*(d+1)-d*c)
         -(2*(3*d+1)*c*c+d*(4*d+1)*c-2*d*(d+1))),
        ('endpoint_gap',(2,),lambda d:(d+1)*(3*d+1)-d*(3*d+4)-1),
        ('dimension_four_threshold',(2,),lambda c:
         (Q(13,2)+4*c)*(20-41*c)-(48-84*c)+Q(41,2)*(8*c*c+5*c-4)),
        ('dimension_four_M9_root',(3,),lambda c:P(1,7,c)-(2-20*c**3)),
        ('dimension_four_M9_value',(3,),lambda c:
         c*(1+2*c)*(1-2*c)-(c-Q(2,5))+4*(c**3-Q(1,10))),
    ]
    records=[]
    for name,degrees,residual in identities:
        count=0
        for node in product(*(tuple(Q(i) for i in range(degree+1)) for degree in degrees)):
            assert residual(*node)==0,(name,node)
            count+=1
        records.append(dict(name=name,degree_bounds=degrees,evaluations=count))
    return records


def run():
    algebra=algebra_checks()
    boundaries=[]
    specs=[(n,b) for n in (3,4,5,8,16) for b in (Q(1,2),Q(1,4),Q(1,10),Q(1,100))]
    specs += [(4,Q(6,25)),(4,Q(239,1000)),(5,Q(99,1000)),(3,Q(0)),(4,Q(0))]
    for n,b in specs:
        record=boundary_root(n,b)
        cert=certify_cell(n,Q(0),b,Q(1,1000000),shortcuts=False)
        lo,hi=record['loss_interval']
        assert cert['lower']<=hi and lo<=cert['upper'],(record,cert['lower'],cert['upper'])
        record['certificates']=[cert]
        boundaries.append(record)
    transitions=[transition_root(n) for n in (4,5,6,8,16,100)]
    assert transitions[0]['last_balanced_budget']==8
    targets=[(4,7),(4,9),(5,18),(5,19),(6,31),(6,32),(8,70),(8,71)]
    targets += [(n,M) for n in (4,5,8,16,100) for M in (4,1000,1000000)]
    budgets=[]
    for n,M in targets:
        record=budget_root(n,M)
        record['certificates']=[]
        if (n,M) in targets[:8]:
            q=M-2
            blo,bhi=record['boundary_interval']
            b=Q(round((blo+bhi)/2*10**9),10**9)
            h=(1-2*b)/(q-2)
            assert 0<h<=b
            knots=(Q(0),)+tuple(b+i*h for i in range(q-1))+(Q(1),)
            cert=certify_partition(n,knots,Q(1,1000000))
            lo,hi=record['loss_interval']
            assert cert['lower']<=hi and lo<=cert['upper']
            assert cert['upper']-cert['lower']<Q(2,1000000)
            record['certificates']=[cert]
        budgets.append(record)
        print(json.dumps(jsonable({k:record[k] for k in ('dimension','generator_budget','loss_interval')})),flush=True)
    existing=json.loads((ROOT/'results/theory/affine_certified.json').read_text())
    comparisons=[]
    for cert in existing['partitions']:
        n,M=cert['dimension'],cert['generator_budget']
        if n<4:
            continue
        root=budget_root(n,M)
        lo,hi=root['loss_interval']
        assert Q(cert['lower'])<=hi and lo<=Q(cert['upper'])
        comparisons.append((n,M))
    return dict(schema='affine-boundary-budget-v1',complete=True,algebra=algebra,
                boundaries=boundaries,budgets=budgets,transitions=transitions,
                existing_D9_comparisons=comparisons,
                scope='D10 root and algebra checks with selected independent D8/D9 global certificate controls; not a machine proof of the universal theorem.')


def main():
    output=ROOT/'results/theory/affine_boundary.json'
    interchange=output.with_suffix('.toml')
    output.write_text('{"complete": false}\n',encoding='utf-8',newline='\n')
    interchange.write_text('complete = false\n',encoding='utf-8',newline='\n')
    report=run()
    sources=[Path(__file__)]+[Path(__file__).with_name(name) for name in
        ('checks.py','affine_profiles.py','affine_dimensions.py','affine_cells.py','affine_two_level.py','affine_certified.py')]

    sources += [ROOT/'results/theory/affine_certified.json']
    report['source_sha256']={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    lines=['schema = "affine-boundary-budget-v1"','complete = true']
    lines.append('existing_D9_comparisons = '+json.dumps(report['existing_D9_comparisons']))
    for family in ('algebra','boundaries','budgets','transitions'):
        for record in report[family]:
            lines += ['', '[['+family+']]']
            for key,value in record.items():
                if key!='certificates' or not value:
                    lines.append(f'{key} = {json.dumps(jsonable(value))}')
            for cert in record.get('certificates',[]):
                if family=='boundaries':
                    write_cell(lines,cert,family+'.certificates')
                else:
                    table=family+'.certificates'
                    lines += ['', '[['+table+']]']
                    for key,value in cert.items():
                        if key!='cells':lines.append(f'{key} = {json.dumps(jsonable(value))}')
                    for cell in cert['cells']:write_cell(lines,cell,table+'.cells')
    interchange.write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
    output.write_text(json.dumps(jsonable(report),indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(boundaries=len(report['boundaries']),budgets=len(report['budgets']),
                          global_partition_controls=sum(bool(r['certificates']) for r in report['budgets']),
                          algebra_evaluations=sum(r['evaluations'] for r in report['algebra']))))


if __name__=='__main__':
    main()
