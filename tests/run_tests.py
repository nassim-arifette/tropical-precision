"""Check retained results and failure controls without rerunning discovery grids."""
import ast
from pathlib import Path
import subprocess
import sys


def main():
    if sys.version_info < (3, 11) or not __debug__:
        raise RuntimeError("Use Python 3.11+ without -O or PYTHONOPTIMIZE")
    root = Path(__file__).resolve().parents[1]
    for directory in ("src", "tests", "scripts", "examples"):
        for path in (root / directory).glob("*.py"):
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for script in ("tests/verify_provenance.py", "tests/test_provenance.py", "examples/relu.py"):
        subprocess.run([sys.executable, "-B", str(root / script)], cwd=root, check=True)
    result = subprocess.run([sys.executable, "-O", "-B", str(root / "tests/verify_provenance.py")],
                            cwd=root, capture_output=True, text=True)
    if result.returncode == 0 or "require assertions" not in result.stderr:
        raise ValueError("Disabled-assertion mode was not rejected correctly")
    print("All tests passed.")


if __name__ == "__main__":
    main()
