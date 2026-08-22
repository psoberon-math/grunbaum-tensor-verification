# Verification manifest

**Release purpose.** This repository supplies the exact finite certificate for
the tensor theorem used in the four-hyperplane counterexample. The certified
finite statement is that the five explicit tensor contractions have no common
zero on the normalized quaternion chart `[-1,1]^6`.

The analytic perturbation argument, the signed-permutation reduction to one
chart, and the passage from `SO(4)` to `O(4)` are proved in the manuscript and
are outside the scope of these programs.

## 1. Single tensor data set

The sole machine-readable tensor definition in this release is

```text
data/final_tensors.json
```

It contains five times the manuscript cubic `A` and the manuscript quartic
`B`. Multiplication of `A` by five clears its only denominator and does not
change any zero set.

The quartic has exactly the following nine nonzero symmetric entries:

```text
B0013=-2  B0033=-4  B0133=3
B0223=-1  B1123=-1  B1133=-1
B1223=-1  B1233=1   B2333=-6
```

Omitted symmetric entries are zero. No exploratory, rationalized, or
alternative quartic is included in this public release.

Both C checkers contain these same entries. The script
`tools/check_release_consistency.py` compares their hard-coded tensors with
`data/final_tensors.json`, regenerates the chart polynomials, and regenerates
the CSV from the text certificate.

## 2. Files named in the manuscript

### `verify_grunbaum_tensors.c`

This is the authoritative one-file checker. It contains the complete
subdivision certificate in an embedded encoding. Starting from the displayed
`5A` and `B`, it:

1. constructs the quaternion-chart frame entries;
2. derives the four cubic and one quartic six-variable polynomials;
3. converts them exactly to a common degree-four tensor-product Bernstein
   representation;
4. creates outward fixed-point interval enclosures of all root control data;
5. replays every midpoint de Casteljau subdivision; and
6. verifies every leaf witness using signed `__int128` integer arithmetic.

It uses no floating-point arithmetic and no third-party library.

Compile and run:

```bash
cc -O3 -std=gnu11 verify_grunbaum_tensors.c -o verify_grunbaum_tensors
./verify_grunbaum_tensors
```

A successful run exits with status `0` and ends with:

```text
CERTIFIED
  nodes:      77713
  leaves:     38857
  max depth:  60
  minimum exact integer leaf margin: 121251517
The five tensor contractions have no common zero on [-1,1]^6.
```

The checker requires GCC or Clang support for signed `__int128`.

### `grunbaum_leaf_boxes.csv`

This is an explicit list of all 38,857 terminal boxes and their integer
separating covectors. For each coordinate, `(k,d)` represents the interval

```text
[k/2^d, (k+1)/2^d]
```

in the Bernstein cube `[0,1]^6`. The exact format is specified in
`certificate/LEAF_BOXES_FORMAT.md`. The boxes have disjoint interiors and
exact total dyadic volume one.

The CSV is redundant with `certificate/grunbaum_certificate.txt`; it is
included so that independent implementations need not parse the preorder
subdivision tree. It does not depend logically on trusting the supplied C or
Python parsers: the rows themselves are the proposed boxes and witnesses.
An independent checker must still reconstruct the exact Bernstein control
vectors and verify the listed inequalities.

### `verification_manifest.md`

This file records the trust boundary, commands, file roles, and release
consistency. Cryptographic hashes for every retained file are in the root
`MANIFEST.sha256`.

## 3. Readable certificate route

The proof data are also separated from the checking program:

```text
certificate/grunbaum_certificate.txt
checkers/check_grunbaum_text.c
```

Compile and run:

```bash
cc -O3 -std=gnu11 checkers/check_grunbaum_text.c -o check_grunbaum_text
./check_grunbaum_text certificate/grunbaum_certificate.txt
```

This checker uses the same exact arithmetic as the self-contained checker and
returns the same node count, leaf count, maximum depth, and minimum exact
integer margin. The certificate format and root-normalization convention are
specified in `certificate/CERTIFICATE_FORMAT.md`.

## 4. Exact polynomial data

`python/generate_final_chart_polynomials.py` derives the final chart
polynomials from `data/final_tensors.json` using only the Python standard
library. The committed output is

```text
data/final_chart_polynomials.json
```

Check exact reproduction with:

```bash
python3 python/generate_final_chart_polynomials.py --check
```

The five term counts are

```text
316, 316, 316, 316, 832.
```

## 5. Optional Python files

The Python programs are supplementary. They are not needed for the exact C
certificate.

| File | Role |
|---|---|
| `python/generate_final_chart_polynomials.py` | Final-only exact tensor-to-chart generator; standard library only. |
| `python/bernstein_tools.py` | Floating Bernstein conversion, de Casteljau subdivision, and separator-search utilities. |
| `python/generate_certificate.py` | Untrusted producer of candidate subdivision trees and quantized witnesses. |
| `python/verify_python_replay.py` | Independent fixed-point interval Bernstein replay. |
| `python/python_replay_output.txt` | Recorded completed output of the cleaned replay, if present. |
| `python/certificate_generation_output.txt` | Recorded producer output, if present. |

Install the recorded NumPy dependency with:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

Run the independent replay:

```bash
PYTHONUNBUFFERED=1 python3 python/verify_python_replay.py \
  --maxnodes 1000000 --report 20000 --maxdepth 120
```

The producer may use floating-point heuristics. Its output supports no theorem
until an exact checker accepts it.

## 6. Trust boundary

The trusted finite computation consists of:

- the C source of an exact checker;
- the exact tensor constants `5A` and `B`;
- integer/rational power-to-Bernstein conversion;
- directed midpoint de Casteljau subdivision;
- the complete tree and integer leaf witnesses; and
- signed `__int128` verification of every strict leaf inequality.

The following are not trusted:

- floating-point separator search;
- split-selection heuristics;
- Python/NumPy as used by the certificate producer; and
- saved terminal transcripts.

They may be wrong without invalidating a certificate accepted by either exact
C checker.

## 7. Recorded checks for this cleaned release

The repository was assembled and checked on Linux x86_64 with GCC 14.2.0.
The following checks were completed:

```text
self-contained C checker: CERTIFIED
readable-certificate C checker: CERTIFIED
text certificate: 77713 nodes, 38857 leaves, maximum depth 60
CSV regenerated from text certificate: byte-for-byte identical
chart JSON regenerated from final_tensors.json: byte-for-byte identical
both C tensor tables compared with final_tensors.json: identical
cleaned exact Python replay: CERTIFIED (141699 nodes, 70850 leaves)
cleaned certificate producer reproduced the readable certificate byte-for-byte
```

The recorded C output is `expected_output.txt`.

## 8. Hash policy

Adding a file does not change the SHA-256 hash of an unchanged existing file.
It does, however, make an exhaustive manifest incomplete. Therefore
`MANIFEST.sha256` was generated only after the final repository contents were
assembled. Any later addition, deletion, rename, or edit requires regenerating
the root manifest and creating a new tagged release.

Verify all retained files with:

Linux:

```bash
sha256sum -c MANIFEST.sha256
```

macOS:

```bash
shasum -a 256 -c MANIFEST.sha256
```

## 9. Key SHA-256 values

```text
b144a05028d09522606fcfc499488b16e24616fba23de08f878fa2367e16fe1d  verify_grunbaum_tensors.c
eec1e986e5cbebedc6cb58d39455f93da2fc1deb1f06d921c1584d7ed10e9274  grunbaum_leaf_boxes.csv
66d91ccd5ded4c9042a37b5ec5abf9d03fe1cf8f813b7aee8580bf488930ac6d  certificate/grunbaum_certificate.txt
aea0b93fc39a89cd2feae60f73d733aa0372f6b401064593600011ce402709fb  data/final_tensors.json
184f0ac9684418345e4795388e9e901ff0d7283d594a53bce380aac982e177a4  data/final_chart_polynomials.json
```

The compact binary reconstructed from the readable certificate, and the binary
regenerated independently by `python/generate_certificate.py`, both have
SHA-256

```text
66381f1454d7611abc10660e62cf33e617a2ff2dddeb032346688c75752f2ceb.
```

The exhaustive root `MANIFEST.sha256` is the authoritative inventory for the
release as a whole.
