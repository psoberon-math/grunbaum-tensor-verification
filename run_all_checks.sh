#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
ROOT="$(cd "$(dirname "$0")" && pwd)"
TMPDIR_CHECK="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_CHECK"' EXIT
CC_BIN="${CC:-cc}"

cd "$ROOT"
python3 python/generate_final_chart_polynomials.py --check
python3 tools/check_release_consistency.py

"$CC_BIN" -O3 -std=gnu11 verify_grunbaum_tensors.c \
  -o "$TMPDIR_CHECK/verify_grunbaum_tensors"
"$TMPDIR_CHECK/verify_grunbaum_tensors"

"$CC_BIN" -O3 -std=gnu11 checkers/check_grunbaum_text.c \
  -o "$TMPDIR_CHECK/check_grunbaum_text"
"$TMPDIR_CHECK/check_grunbaum_text" certificate/grunbaum_certificate.txt
