#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
TMPDIR_CHECK="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_CHECK"' EXIT
CC_BIN="${CC:-cc}"
"$CC_BIN" -O3 -std=gnu11 "$ROOT/verify_grunbaum_tensors.c" \
  -o "$TMPDIR_CHECK/verify_grunbaum_tensors"
"$TMPDIR_CHECK/verify_grunbaum_tensors"
