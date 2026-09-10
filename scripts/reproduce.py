"""Regenerate exact results and verify their certificates using Python and Julia."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--julia", default=shutil.which("julia"))
    parser.add_argument("--certificates-only", action="store_true")
    args = parser.parse_args()
    if sys.version_info < (3, 11) or not __debug__:
        parser.error("Use Python 3.11+ without -O or PYTHONOPTIMIZE")
    if not args.julia:
        parser.error("Julia was not found; provide --julia PATH")
    root = Path(__file__).resolve().parents[1]
    (root / "results/theory").mkdir(parents=True, exist_ok=True)
    if not args.certificates_only:
        for name in ("checks", "affine_profiles", "partition_reduction", "affine_upper",
                     "affine_dimensions", "affine_cells", "affine_two_level",
                     "affine_certified", "affine_boundary"):
            print(f"Running {name}", flush=True)
            subprocess.run([sys.executable, "-B", f"src/{name}.py"], cwd=root, check=True)
    subprocess.run([args.julia, "--startup-file=no", "--threads=1", "tests/run_all.jl"], cwd=root, check=True)
    for script in ("tests/verify_provenance.py", "tests/test_provenance.py"):
        subprocess.run([sys.executable, "-B", script], cwd=root, check=True)
    neural = [sys.executable, "-B", "examples/relu.py"]
    if not args.certificates_only:
        neural.append("--write")
    subprocess.run(neural, cwd=root, check=True)
    print("Reproduction completed successfully.")


if __name__ == "__main__":
    main()
