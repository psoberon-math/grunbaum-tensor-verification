# Optional Python reproduction files

None of these files is required for the authoritative C verification.

## Files

- `generate_final_chart_polynomials.py` derives the exact one-chart polynomial
  data from `data/final_tensors.json`. It contains no alternative quartic.
- `bernstein_tools.py` contains floating Bernstein conversion, subdivision,
  and separator-search utilities.
- `generate_certificate.py` is the untrusted floating-point producer for a
  subdivision tree and integer leaf witnesses.
- `verify_python_replay.py` is an independent interval Bernstein replay.

## Regenerate and check the exact chart data

```bash
python3 python/generate_final_chart_polynomials.py
python3 python/generate_final_chart_polynomials.py --check
```

## Run the independent exact replay

```bash
PYTHONUNBUFFERED=1 python3 python/verify_python_replay.py \
  --maxnodes 1000000 --report 20000 --maxdepth 120
```

The expected final conclusion is

```text
CERTIFIED: exact interval Bernstein separation excludes a common zero on [-1,1]^6
```

## Regenerate a candidate serialized certificate

```bash
PYTHONUNBUFFERED=1 python3 python/generate_certificate.py \
  --output certificate/grunbaum_certificate.generated.bin \
  --maxnodes 1000000 --maxdepth 120 --report 10000
```

This program is not trusted: it uses floating-point heuristics to choose
splits and separators. Its output has mathematical force only after an exact
checker accepts it.
