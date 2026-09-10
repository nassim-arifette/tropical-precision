"""Read-only provenance and interchange checks for all nine exact studies.

This verifies coverage, current source/input bytes, mathematical export payloads,
and scoped verifier counts. It does not rerun the mathematical computations.
"""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import tomllib

if not __debug__:
    raise RuntimeError("Provenance checks require assertions: remove -O and PYTHONOPTIMIZE")

ROOT = Path(__file__).resolve().parents[1]
STUDY = "src/"
VERIFIERS = "tests/"
RESULTS = "results/theory/"


def spec(scripts, verifiers, schema, negative, receipt=None, extras=()):
    return dict(sources={STUDY+s+".py" for s in scripts} | set(extras),
                verifiers={VERIFIERS+s+".jl" for s in verifiers},
                schema=schema, negative=negative, receipt=receipt)


CORE = ["checks", "affine_profiles", "affine_dimensions", "affine_cells"]
WITNESS = ["verify_witnesses"]
CELLS = WITNESS+["verify_affine_cells"]
CERTIFIED = CELLS+["verify_affine_certified"]
SPECS = {
    "exact_checks": spec(["checks"], WITNESS,
                         "bounded-context-witnesses-v1", 7, "julia_verification"),
    "affine_profiles": spec(["checks","affine_profiles"], WITNESS, "bounded-context-witnesses-v1", 7),
    "partition_reduction": spec(["checks","partition_reduction"], WITNESS, "diagonal-partition-certificates-v1", 7),
    "affine_upper": spec(["checks","affine_profiles","affine_upper"], ["verify_affine_upper"], "affine-upper-algebra-v1", 5),
    "affine_dimensions": spec(["checks","affine_profiles","affine_dimensions"], WITNESS, "bounded-context-witnesses-v1", 7),
    "affine_cells": spec(CORE, CELLS, "affine-cell-controls-v1", 8),
    "affine_two_level": spec(CORE+["affine_two_level"], CELLS+["verify_affine_two_level"], "affine-two-level-controls-v1", 12),
    "affine_certified": spec(CORE+["affine_two_level","affine_certified"],
                             CERTIFIED, "affine-certified-minimax-v1", 19),
    "affine_boundary": spec(CORE+["affine_two_level","affine_certified","affine_boundary"],
                            CERTIFIED+["verify_affine_boundary"], "affine-boundary-budget-v1", 20,
                            extras=[RESULTS+"affine_certified.json"]),
}


def digest(path):
    data = path.read_bytes()
    if path.suffix in (".py", ".jl", ".md", ".json", ".toml") and b"\r\n" in data:
        raise ValueError(f"Noncanonical newlines: {path.relative_to(ROOT)}")
    return hashlib.sha256(data).hexdigest()


def hashes_current(recorded, expected):
    assert set(recorded) == set(expected), "Incomplete source-hash coverage"
    for file, value in recorded.items():
        assert digest(ROOT/file) == value, f"Stale receipt source: {file}"


def legacy_witness(case):
    names = ("dimension", "generator_count", "generators", "target", "context",
             "affine_coefficients", "minimum_concrete_union", "target_value", "certified_gap")
    result = {name:case[name] for name in names}
    result["coefficients"] = case["membership_certificate"]["coefficients"]
    return result


def mathematical_export(stem, data):
    result = dict(schema=SPECS[stem]["schema"], complete=data["complete"])
    if stem in ("exact_checks","affine_dimensions"):
        result["witnesses"] = [legacy_witness(c) for c in data["continuous_diagonal_certificates"]]
    elif stem == "affine_profiles":
        records = []
        for item in data["results"]+[data["analytic_family"]]:
            family = "parameter" in item
            gap = item["certified_gap"] if family else item["certified_lower_bound"]
            record = {k:item[k] for k in ("generators","target","context","affine_coefficients")}
            record.update(dimension=3, generator_count=len(item["generators"]),
                          minimum_concrete_union="0" if family else item["concrete_minimum"],
                          target_value=str(-Fraction(gap)) if family else item["target_value"],
                          certified_gap=gap, coefficients=item["reconstruction_certificate"]["coefficients"])
            records.append(record)
        result["witnesses"] = records
    elif stem == "partition_reduction":
        result["cases"] = data["cases"]
    elif stem == "affine_upper":
        result.update(identity_names=[r["name"] for r in data["algebra"]],
                      degree_bounds=[r["degree_bounds"] for r in data["algebra"]],
                      evaluations=[r["evaluations"] for r in data["algebra"]],
                      mixed_remainder_lower_bound="7/128", optimum_coefficients=["-3","5"])
    elif stem in ("affine_cells","affine_two_level"):
        result["witnesses"] = data["witnesses"]
    elif stem == "affine_certified":
        result.update(controls=data["controls"], algebra=data["algebra"],
                      partitions=[{k:v for k,v in r.items() if k!="certificate_basis"}
                                  for r in data["partitions"]])
    elif stem == "affine_boundary":
        result.update({k:data[k] for k in ("algebra","boundaries","budgets","transitions","existing_D9_comparisons")})
    else:
        raise ValueError(stem)
    return result


def scalar_leaves(cell):
    return len(cell.get("leaves", []))


def verify(stem, verifier_receipt=None):
    configuration = SPECS[stem]
    directory = ROOT/RESULTS
    data = json.loads((directory/(stem+".json")).read_text(encoding="utf-8"))
    export_path = directory/(stem+".toml")
    export = tomllib.loads(export_path.read_text(encoding="utf-8"))
    receipt_name = configuration["receipt"] or stem+"_julia"
    if verifier_receipt is not None:
        assert verifier_receipt == receipt_name, "Unexpected verifier receipt"
    receipt = tomllib.loads((directory/(receipt_name+".toml")).read_text(encoding="utf-8"))
    assert data["complete"] is True and export["complete"] is True and receipt["complete"] is True, "Incomplete receipt"
    hashes_current(data["source_sha256"], configuration["sources"])
    assert export == mathematical_export(stem,data), f"Python/Julia mathematical payload differs: {stem}"
    assert digest(export_path) == receipt["input_sha256"], "Stale Julia input"
    if "source_sha256" in receipt:
        hashes_current(receipt["source_sha256"], configuration["verifiers"])
    else:
        assert len(configuration["verifiers"]) == 1, "Missing Julia dependency hashes"
        assert digest(ROOT/next(iter(configuration["verifiers"]))) == receipt["verifier_sha256"], "Stale Julia verifier"
    negative_key = "malformed_witnesses_rejected" if configuration["verifiers"] == set(VERIFIERS+s+".jl" for s in WITNESS) else "malformed_inputs_rejected"
    assert receipt[negative_key] == configuration["negative"], "Negative controls missing"
    if "witnesses" in export:
        assert len(export["witnesses"]) == receipt["witnesses_verified"]
        assert [Fraction(r["certified_gap"]) for r in export["witnesses"]] == [
            Fraction(s.replace("//","/")) for s in receipt["gaps"]], "Verified gap inventory differs"
    elif stem == "partition_reduction":
        assert len(export["cases"]) == receipt["partitions_verified"]
        assert sum(len(c["coefficients"]) for c in export["cases"]) == receipt["generator_certificates_verified"]
    elif stem == "affine_upper":
        assert len(export["identity_names"]) == receipt["identities_verified"] == 6
        assert sum(export["evaluations"]) == receipt["polynomial_evaluations"] == 102
    elif stem == "affine_certified":
        assert len(export["controls"]) == receipt["cell_controls_verified"]
        assert len(export["partitions"]) == receipt["partitions_verified"]
        leaves = sum(map(scalar_leaves, export["controls"]))
        leaves += sum(scalar_leaves(c) for r in export["partitions"] for c in r["cells"])
        assert leaves == receipt["scalar_intervals_verified"]
        assert export["algebra"]["evaluations"] == receipt["algebra_evaluations"] == 69
    elif stem == "affine_boundary":
        for family,key in (("boundaries","boundary_cells_verified"),("budgets","budget_roots_verified"),
                           ("transitions","transition_roots_verified")):
            assert len(export[family]) == receipt[key]
        assert sum(len(r["certificates"]) for r in export["budgets"]) == receipt["global_partition_controls_verified"]
        leaves = sum(scalar_leaves(c) for r in export["boundaries"] for c in r["certificates"])
        leaves += sum(scalar_leaves(c) for r in export["budgets"] for p in r["certificates"] for c in p["cells"])
        assert leaves == receipt["scalar_intervals_verified"]
        assert sum(r["evaluations"] for r in export["algebra"]) == receipt["algebra_evaluations"] == 74
    return receipt


def verify_failures():
    directory = ROOT/RESULTS
    legacy = tomllib.loads((directory/"verifier_failures.toml").read_text(encoding="utf-8"))
    assert legacy["complete"] is True
    assert legacy["verifier_sha256"] == digest(ROOT/VERIFIERS/"verify_witnesses.jl")
    assert legacy["runner_sha256"] == digest(ROOT/VERIFIERS/"check_verifier_failures.jl")
    assert legacy["input_sha256"] == digest(directory/"exact_checks.toml")
    assert set(legacy["checks_passed"]) == {"missing_input","malformed_toml","invalid_certificate","failed_negative_control"}
    certified = tomllib.loads((directory/"affine_certified_failures.toml").read_text(encoding="utf-8"))
    assert certified["complete"] is True
    hashes_current(certified["source_sha256"], {VERIFIERS+s+".jl" for s in CERTIFIED+["check_affine_certified_failures"]})
    assert certified["input_sha256"] == digest(directory/"affine_certified.toml")
    assert set(certified["checks_passed"]) == {"missing_input","malformed_toml","invalid_certificate","reporting_after_all_checks"}
    return 8


def verify_all():
    receipts = {stem:verify(stem) for stem in SPECS}
    failure_count = verify_failures()
    print(f"Current provenance verified: {len(receipts)} studies, full mathematical interchange payloads, "
          f"source/input coverage, and {failure_count} command-failure controls.")
    return receipts


if __name__ == "__main__":
    verify_all()
