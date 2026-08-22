#!/usr/bin/env python3
"""Check internal consistency of the public verification release.

This is a packaging/audit utility, not part of the mathematical proof. It
checks that both C checkers contain exactly the tensors in final_tensors.json,
that the committed chart polynomials are regenerated exactly, and that the
leaf CSV is reproduced byte-for-byte from the readable certificate.
"""
from __future__ import annotations

import base64
import importlib.util
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def load_expected() -> tuple[dict[tuple[int, ...], int], dict[tuple[int, ...], int]]:
    data = json.loads((ROOT / "data/final_tensors.json").read_text(encoding="utf-8"))
    cubic = {tuple(i): int(v) for i, v in data["cubic"]["entries"]}
    quartic = {tuple(i): int(v) for i, v in data["quartic"]["entries"]}
    return cubic, quartic


def function_body(source: str, function_name: str) -> str:
    marker = f"static int64_t {function_name}("
    start = source.find(marker)
    if start < 0:
        fail(f"could not find {function_name} in C source")
    brace = source.find("{", start)
    depth = 0
    for position in range(brace, len(source)):
        if source[position] == "{":
            depth += 1
        elif source[position] == "}":
            depth -= 1
            if depth == 0:
                return source[brace + 1:position]
    fail(f"unterminated {function_name} in C source")
    raise AssertionError


def parse_c_tensor(path: Path, function_name: str, degree: int) -> dict[tuple[int, ...], int]:
    body = function_body(path.read_text(encoding="utf-8"), function_name)
    pattern = re.compile(r"case\s+\d+\s*:\s*return\s+(-?\d+)\s*;\s*/\*\s*(\d{%d})\s*\*/" % degree)
    entries: dict[tuple[int, ...], int] = {}
    for value, indices in pattern.findall(body):
        key = tuple(int(c) for c in indices)
        entries[key] = int(value)
    if not entries:
        fail(f"found no entries for {function_name} in {path}")
    return entries


def load_certificate_module():
    path = ROOT / "tools/certificate_text_tools.py"
    spec = importlib.util.spec_from_file_location("certificate_text_tools", path)
    if spec is None or spec.loader is None:
        fail("could not load certificate_text_tools.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    required = [
        ROOT / "verify_grunbaum_tensors.c",
        ROOT / "verification_manifest.md",
        ROOT / "grunbaum_leaf_boxes.csv",
    ]
    for path in required:
        if not path.is_file():
            fail(f"missing manuscript-listed file: {path.name}")

    cubic, quartic = load_expected()
    for checker in [ROOT / "verify_grunbaum_tensors.c",
                    ROOT / "checkers/check_grunbaum_text.c"]:
        if parse_c_tensor(checker, "A5", 3) != cubic:
            fail(f"cubic tensor mismatch in {checker}")
        if parse_c_tensor(checker, "B4", 4) != quartic:
            fail(f"quartic tensor mismatch in {checker}")
        print(f"OK tensors: {checker.relative_to(ROOT)}")

    subprocess.run(
        [sys.executable, str(ROOT / "python/generate_final_chart_polynomials.py"), "--check"],
        cwd=ROOT,
        check=True,
    )

    cert_tools = load_certificate_module()
    certificate = cert_tools.read_text(ROOT / "certificate/grunbaum_certificate.txt")
    summary = cert_tools.validate_tree(certificate)
    if summary["nodes"] != 77713 or summary["leaves"] != 38857 or summary["max_depth"] != 60:
        fail(f"unexpected certificate statistics: {summary}")
    print("OK readable certificate structure")

    with tempfile.TemporaryDirectory() as directory:
        regenerated = Path(directory) / "grunbaum_leaf_boxes.csv"
        binary = Path(directory) / "grunbaum_certificate.bin"
        cert_tools.write_leaf_csv(certificate, regenerated)
        cert_tools.write_binary(certificate, binary)
        if regenerated.read_bytes() != (ROOT / "grunbaum_leaf_boxes.csv").read_bytes():
            fail("grunbaum_leaf_boxes.csv is not reproduced by the readable certificate")

        c_source = (ROOT / "verify_grunbaum_tensors.c").read_text(encoding="utf-8")
        start = c_source.find("const char CERTIFICATE_BASE64[] =")
        if start < 0:
            fail("embedded certificate declaration not found")
        encoded_region = c_source[start:]
        chunks = re.findall(r'"([A-Za-z0-9+/=]+)"', encoded_region)
        if not chunks:
            fail("embedded base64 certificate not found")
        embedded = base64.b64decode("".join(chunks), validate=True)
        if embedded != binary.read_bytes():
            fail("embedded C certificate differs from the readable certificate")
    print("OK leaf CSV byte-for-byte regeneration")
    print("OK embedded certificate equals readable certificate")

    prohibited = ["B5_SORTED", "SIMPLE_B_SORTED", "defwall_exact_polynomials"]
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.name == "check_release_consistency.py":
            continue
        if path.suffix.lower() not in {".py", ".c", ".h", ".json", ".md", ".txt"}:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        for term in prohibited:
            if term in text:
                fail(f"obsolete identifier {term!r} remains in {path.relative_to(ROOT)}")
    print("OK no obsolete quartic identifiers")
    print("RELEASE CONSISTENCY CHECK PASSED")


if __name__ == "__main__":
    main()
