# Verification files for the four-hyperplane tensor obstruction

This repository contains the exact finite verification used in the paper
*Four hyperplanes do not always equipartition a mass in* `R^4`.

The finite statement checked here is:

> The five explicit contractions determined by the cubic tensor `A` and the
> quartic tensor `B` have no common zero on the normalized quaternion chart
> `[-1,1]^6`.

The analytic perturbation argument and the symmetry reduction from `O(4)` to
this chart are mathematical arguments in the manuscript and are not checked by
these programs.

## Fastest verification: one self-contained C file

The three files named in the manuscript are at the repository root:

- `verify_grunbaum_tensors.c` — self-contained exact checker;
- `verification_manifest.md` — detailed trust boundary, commands, and file inventory;
- `grunbaum_leaf_boxes.csv` — explicit terminal boxes and separating witnesses.

Compile and run the authoritative checker with GCC or Clang:

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

The checker requires compiler support for signed `__int128`. It uses no
floating-point arithmetic and no third-party library.

## Readable proof data

The certificate is also supplied separately from the checker:

- `certificate/grunbaum_certificate.txt` gives the full subdivision tree and
  one exact integer witness at each leaf;
- `grunbaum_leaf_boxes.csv` lists all terminal boxes explicitly in dyadic
  coordinates, together with their witnesses;
- `checkers/check_grunbaum_text.c` verifies the readable certificate directly;
- `certificate/CERTIFICATE_FORMAT.md` and
  `certificate/LEAF_BOXES_FORMAT.md` specify the formats and exact arithmetic.

To run the checker for the readable certificate:

```bash
cc -O3 -std=gnu11 checkers/check_grunbaum_text.c -o check_grunbaum_text
./check_grunbaum_text certificate/grunbaum_certificate.txt
```

## Exact tensor and polynomial data

`data/final_tensors.json` is the sole machine-readable tensor definition in
this release. It contains `5A` and the sparse quartic `B` printed in the
manuscript. There is no alternative or exploratory quartic in this repository.

The exact chart polynomials are in `data/final_chart_polynomials.json`. They
can be regenerated and compared byte for byte with:

```bash
python3 python/generate_final_chart_polynomials.py --check
```

## Optional Python reproduction

The Python programs document how the subdivision was found and provide an
independent rigorous replay. They are not needed for the C certificate.
Install the recorded NumPy version with:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

See `python/README.md` and `verification_manifest.md` for commands and the
precise trust boundary.

## Hashes

After the repository contents are finalized, every retained file except the
manifest itself is listed in `MANIFEST.sha256`.

Linux:

```bash
sha256sum -c MANIFEST.sha256
```

macOS:

```bash
shasum -a 256 -c MANIFEST.sha256
```
