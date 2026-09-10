"""Reject changed scientific payloads even when completion/source flags survive.

Only a checked temporary copy is mutated. These are provenance controls, not
new mathematical experiments or tests of the universal theorems.
"""
import json
from pathlib import Path
import shutil
import tempfile

import verify_provenance as provenance


def main():
    source = provenance.ROOT
    work = (source/".build").resolve()
    work.mkdir(exist_ok=True)
    passed = []
    with tempfile.TemporaryDirectory(prefix="provenance-controls-", dir=work) as directory:
        copied = Path(directory).resolve()
        assert copied.is_relative_to(work)
        for folder in ("src", "tests", "results/theory"):
            shutil.copytree(source/folder, copied/folder, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        provenance.ROOT = copied
        try:
            provenance.verify_all()
            mutations = {
                "exact_checks": lambda d: d["continuous_diagonal_certificates"][0]["target"].__setitem__(0,"999"),
                "affine_profiles": lambda d: d["results"][0]["target"].__setitem__(0,"999"),
                "partition_reduction": lambda d: d["cases"][0]["partition"].__setitem__(0,"999"),
                "affine_upper": lambda d: d["algebra"][0].__setitem__("evaluations",999),
                "affine_dimensions": lambda d: d["continuous_diagonal_certificates"][0]["target"].__setitem__(0,"999"),
                "affine_cells": lambda d: d["witnesses"][0]["target"].__setitem__(0,"999"),
                "affine_two_level": lambda d: d["witnesses"][0]["target"].__setitem__(0,"999"),
                "affine_certified": lambda d: d["controls"][0].__setitem__("upper","999"),
                "affine_boundary": lambda d: d["budgets"][0]["loss_interval"].__setitem__(0,"999"),
            }
            for stem,change in list(mutations.items())+[
                ("affine_boundary",lambda d:d.__setitem__("source_sha256",{})),
                ("affine_certified",lambda d:d.__setitem__("complete",False)),
            ]:
                path = copied/provenance.RESULTS/(stem+".json")
                original = path.read_bytes()
                data = json.loads(original)
                change(data)
                path.write_text(json.dumps(data),encoding="utf-8",newline="\n")
                rejected = False
                try:
                    provenance.verify(stem)
                except AssertionError:
                    rejected = True
                finally:
                    path.write_bytes(original)
                assert rejected, f"Accepted changed evidence: {stem}"
                passed.append(stem)
            provenance.verify_all()
        finally:
            provenance.ROOT = source
    print(f"Rejected {len(passed)} changed payload/completion/coverage records; original evidence untouched.")


if __name__ == "__main__":
    main()
