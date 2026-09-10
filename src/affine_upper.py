"""Bounded exact checks of the sharp M=3,n=3 affine upper-bound proof.

Context supports are continuous and exact; point/direction checks are finite.
Polynomial identities are checked on unisolvent tensor grids with explicit
degree bounds, which establishes those identities, not the entire theorem.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement, product
from pathlib import Path
import hashlib
import json

from checks import ROOT
from affine_profiles import weight_shell, sector_support


def algebra_checks():
    identities = [
        ("weighted_denominators", (3, 3, 3),
         lambda a, b, c: ((a+b)*c+(a-b)**2)*c
         -(a+b)*(c*c-(a-b)**2)-(a+b+c)*(a-b)**2),
        ("mixed_cross_multiplication", (2, 3),
         lambda a, k: (1-k)**2*a*(1-2*a+k)
         -(1+k)**2*(a-k)*(1-2*a)
         -k*(8*a*a-(1+k)*(5+k)*a+(1+k)**2)),
        ("mixed_square_decomposition", (2, 4),
         lambda a, k: 8*a*a-(1+k)*(5+k)*a+(1+k)**2
         -8*(a-(1+k)*(5+k)/16)**2
         -(1+k)**2*(1-(5+k)**2/32)),
        ("positive_remainder", (2,),
         lambda k: 1-(5+k)**2/32-Q(7, 128)
         -(Q(1, 2)-k)*(Q(21, 2)+k)/32),
        ("derivative_numerator", (3,),
         lambda r: (2*r-3*r*r)*(1+r)-(r*r-r**3)
         -2*r*(1-r-r*r)),
        ("golden_value", (3,),
         lambda r: r*r*(1-r)-(5*r-3)*(1+r)
         +(r+3)*(r*r+r-1)),
    ]
    records = []
    for name, degrees, residual in identities:
        nodes = product(*(tuple(Q(i) for i in range(d+1)) for d in degrees))
        count = 0
        for node in nodes:
            assert residual(*node) == 0, (name, node)
            count += 1
        records.append(dict(name=name, degree_bounds=degrees, evaluations=count))
    return records


def finite_checks(D=8, W=12):
    counts = dict(all_same_sign=0, negative_min=0, negative_strict_max=0,
                  mixed=0, positive_defect=0)
    pairs, best, maximizer = 0, Q(0), None
    for x in combinations_with_replacement(range(D+1), 3):
        if x[0] == x[-1]:
            continue  # Diagonal points have no positive deficit, analytically.
        missing = [(k, 1) for k in range(3) if x[k] > min(x)]
        missing += [(k, -1) for k in range(3) if x[k] < max(x)]
        for w in weight_shell(W):
            pairs += 1
            values = [min(0, sum(w))*D-sum(a*b for a, b in zip(w, x))]
            values += [sector_support(x, w, k, sign, D)[0] for k, sign in missing]
            delta = Q(min(values), D*W)
            if delta > best:
                best, maximizer = delta, dict(point=x, direction=w)
            if all(v >= 0 for v in w) or all(v <= 0 for v in w):
                counts["all_same_sign"] += 1
                assert delta <= 0
                continue
            xx, ww = list(x), list(w)
            if sum(v < 0 for v in ww) == 2:
                xx, ww = [D-v for v in xx], [-v for v in ww]
            negative = next(j for j in range(3) if ww[j] < 0)
            c = Q(-ww[negative], W)
            A = 1-c
            upper = max(Q(0), c*c*(A-c)/(A*A))
            assert delta <= upper, (x, w, delta, upper)
            positive = [j for j in range(3) if j != negative]
            case = ("negative_min" if xx[negative] <= min(xx[j] for j in positive)
                    else "negative_strict_max" if xx[negative] > max(xx[j] for j in positive)
                    else "mixed")
            counts[case] += 1
            counts["positive_defect"] += int(delta > 0)
            if delta > 0:
                assert c < A
                if case == "mixed":
                    low = min(positive, key=lambda j: xx[j])
                    a, k = Q(ww[low], W), 1-2*c
                    assert 0 < k < a < Q(1, 2)
                    alpha, d = 1-2*a, a-k
                    bound = d*alpha*k/(a*(alpha+k))
                    assert delta <= bound <= upper
    assert sum(counts[k] for k in counts if k != "positive_defect") == pairs
    return dict(coordinate_denominator=D, direction_denominator=W,
                point_direction_pairs=pairs, case_counts=counts,
                grid_maximum=str(best), grid_maximizer=maximizer,
                scope="Finite exact checks of the global bound; this grid is not an optimum proof.")


def main():
    output = ROOT/"results/theory/affine_upper.json"
    interchange = output.with_suffix(".toml")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text('{"complete": false}\n', encoding="utf-8", newline="\n")
    interchange.write_text("complete = false\n", encoding="utf-8", newline="\n")
    algebra, finite = algebra_checks(), finite_checks()
    sources = [Path(__file__), Path(__file__).with_name("affine_profiles.py"), Path(__file__).with_name("checks.py")]
    report = dict(
        complete=True, arithmetic="exact integers and fractions.Fraction",
        algebra=algebra, finite_checks=finite,
        source_sha256={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sources},
        scope="Algebra identities plus finite sector-support checks; the global case proof is written mathematics.",
    )
    fields = dict(
        schema="affine-upper-algebra-v1", complete=True,
        identity_names=[c["name"] for c in algebra],
        degree_bounds=[list(c["degree_bounds"]) for c in algebra],
        evaluations=[c["evaluations"] for c in algebra],
        mixed_remainder_lower_bound="7/128", optimum_coefficients=["-3", "5"],
    )
    lines = [f"{key} = {json.dumps(value)}" for key, value in fields.items()]
    interchange.write_text("\n".join(lines)+"\n", encoding="utf-8", newline="\n")
    output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps(dict(algebra_identities=len(algebra), **finite)))


if __name__ == "__main__":
    main()
