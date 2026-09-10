"""Compile the canonical manuscript and optionally package the same paper for arXiv."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arxiv", action="store_true",
                        help="also write the ignored arxiv/arxiv-source.zip and output/pdf/arxiv.pdf")
    args = parser.parse_args()
    engine = shutil.which("pdflatex")
    if engine is None:
        raise RuntimeError("Install a TeX distribution with PDFLaTeX and add pdflatex to PATH.")

    sources = sorted(PAPER.rglob("*.tex"))
    if not sources or not (PAPER / "main.tex").is_file():
        raise RuntimeError("The manuscript source paper/main.tex is missing.")
    for source in sources:
        for name in re.findall(r"\\input\{([^}]+)\}", source.read_text(encoding="utf-8")):
            included = (PAPER / name).with_suffix(".tex").resolve()
            if not included.is_relative_to(PAPER.resolve()) or not included.is_file():
                raise RuntimeError(f"Missing or external manuscript input: {name}")

    work = ROOT / ".build/paper"
    work.mkdir(parents=True, exist_ok=True)
    for run in range(1, 4):
        transcript = work / f"pdflatex-{run}.txt"
        with transcript.open("w", encoding="utf-8") as stream:
            result = subprocess.run(
                [engine, "-interaction=nonstopmode", "-halt-on-error", "-file-line-error",
                 "-no-shell-escape", "-recorder", f"-output-directory={work}", "main.tex"],
                cwd=PAPER, stdout=stream, stderr=subprocess.STDOUT,
            )
        if result.returncode:
            raise RuntimeError(f"PDFLaTeX failed; see {transcript}")
    log = (work / "main.log").read_text(encoding="utf-8", errors="replace")
    if re.search(r"undefined|multiply defined|destination with the same identifier|Overfull \\[hv]box", log):
        raise RuntimeError(f"Resolve reference or layout warnings in {work / 'main.log'}")
    match = re.search(r"Output written on .*?\((\d+) pages?", log, flags=re.S)
    if match is None or not (work / "main.pdf").is_file():
        raise RuntimeError("PDFLaTeX did not report a completed PDF.")

    destination = PAPER / "main.pdf"
    shutil.copyfile(work / "main.pdf", destination)
    report = {"pdf": str(destination), "pages": int(match.group(1)), "source_files": len(sources)}
    if args.arxiv:
        upload = ROOT / "arxiv/arxiv-source.zip"
        upload.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(upload, "w", zipfile.ZIP_DEFLATED) as bundle:
            for source in sources:
                bundle.write(source, source.relative_to(PAPER).as_posix())
        preview = ROOT / "output/pdf/arxiv.pdf"
        preview.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(destination, preview)
        report.update(arxiv_sources=str(upload), arxiv_pdf=str(preview))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError) as error:
        print(f"Paper build failed: {error}", file=sys.stderr)
        sys.exit(1)
