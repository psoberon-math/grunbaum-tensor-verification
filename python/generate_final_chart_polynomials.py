#!/usr/bin/env python3
"""Generate the five exact polynomial contractions on the final p0=r0=1 chart.

The only tensor input is ``data/final_tensors.json``.  It contains five times
A (to clear the denominator in the manuscript cubic) and the sparse integer
quartic B printed in the manuscript.  No exploratory or alternative quartic
is present in this public generator.

The output is ``data/final_chart_polynomials.json``.  Its five components are
four contractions of 5A and one contraction of B, in variables
(x1,x2,x3,y1,y2,y3) corresponding to p=(1,x1,x2,x3), r=(1,y1,y2,y3).
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

NVAR_MASTER = 8
Exponent = Tuple[int, ...]
Polynomial = Dict[Exponent, int]


def clean(poly: Mapping[Exponent, int]) -> Polynomial:
    return {m: int(c) for m, c in poly.items() if c}


def add(*polys: Mapping[Exponent, int]) -> Polynomial:
    out: defaultdict[Exponent, int] = defaultdict(int)
    for poly in polys:
        for monomial, coefficient in poly.items():
            out[monomial] += int(coefficient)
    return clean(out)


def scale(poly: Mapping[Exponent, int], scalar: int) -> Polynomial:
    return clean({m: scalar * c for m, c in poly.items()})


def multiply(left: Mapping[Exponent, int], right: Mapping[Exponent, int]) -> Polynomial:
    out: defaultdict[Exponent, int] = defaultdict(int)
    for a, ca in left.items():
        for b, cb in right.items():
            out[tuple(x + y for x, y in zip(a, b))] += ca * cb
    return clean(out)


def variable(index: int) -> Polynomial:
    monomial = [0] * NVAR_MASTER
    monomial[index] = 1
    return {tuple(monomial): 1}


P = [variable(i) for i in range(4)]
R = [variable(4 + i) for i in range(4)]


def quaternion_multiply(a: Sequence[Polynomial], b: Sequence[Polynomial]) -> List[Polynomial]:
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return [
        add(multiply(aw, bw), scale(multiply(ax, bx), -1),
            scale(multiply(ay, by), -1), scale(multiply(az, bz), -1)),
        add(multiply(aw, bx), multiply(ax, bw), multiply(ay, bz),
            scale(multiply(az, by), -1)),
        add(multiply(aw, by), scale(multiply(ax, bz), -1), multiply(ay, bw),
            multiply(az, bx)),
        add(multiply(aw, bz), multiply(ax, by), scale(multiply(ay, bx), -1),
            multiply(az, bw)),
    ]


def spin_matrix() -> List[List[Polynomial]]:
    """Return M[row][column] for columns p*e_j*conj(r)."""
    conjugate_r = [R[0], scale(R[1], -1), scale(R[2], -1), scale(R[3], -1)]
    one = {(0,) * NVAR_MASTER: 1}
    zero: Polynomial = {}
    columns: List[List[Polynomial]] = []
    for j in range(4):
        basis = [zero, zero, zero, zero]
        basis[j] = one
        columns.append(quaternion_multiply(quaternion_multiply(P, basis), conjugate_r))
    return [[columns[column][row] for column in range(4)] for row in range(4)]


def read_tensor_entries(path: Path) -> tuple[dict[tuple[int, ...], int], dict[tuple[int, ...], int]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("format") != "grunbaum-final-tensors-v1":
        raise ValueError(f"unexpected tensor format in {path}")
    cubic = {tuple(indices): int(value) for indices, value in data["cubic"]["entries"]}
    quartic = {tuple(indices): int(value) for indices, value in data["quartic"]["entries"]}
    if any(tuple(sorted(k)) != k or len(k) != 3 for k in cubic):
        raise ValueError("malformed cubic tensor entry")
    if any(tuple(sorted(k)) != k or len(k) != 4 for k in quartic):
        raise ValueError("malformed quartic tensor entry")
    return cubic, quartic


def tensor_evaluate(entries: Mapping[tuple[int, ...], int], vectors: Sequence[Sequence[Polynomial]]) -> Polynomial:
    degree = len(vectors)
    out: Polynomial = {}
    for indices in itertools.product(range(4), repeat=degree):
        coefficient = entries.get(tuple(sorted(indices)), 0)
        if coefficient == 0:
            continue
        term: Polynomial = {(0,) * NVAR_MASTER: coefficient}
        for vector, index in zip(vectors, indices):
            term = multiply(term, vector[index])
        out = add(out, term)
    return out


def chart_p0_r0(poly: Mapping[Exponent, int]) -> Polynomial:
    """Substitute p0=r0=1, retaining p1,p2,p3,r1,r2,r3 in this order."""
    keep = (1, 2, 3, 5, 6, 7)
    out: defaultdict[Exponent, int] = defaultdict(int)
    for monomial, coefficient in poly.items():
        out[tuple(monomial[i] for i in keep)] += coefficient
    return clean(out)


def encode(poly: Mapping[Exponent, int]) -> list[list[object]]:
    return [[list(monomial), coefficient] for monomial, coefficient in sorted(poly.items())]


def generate(tensor_path: Path) -> dict[str, object]:
    cubic, quartic = read_tensor_entries(tensor_path)
    matrix = spin_matrix()
    columns = [[matrix[row][column] for row in range(4)] for column in range(4)]

    components: list[Polynomial] = []
    for omit in range(4):
        selected = [columns[j] for j in range(4) if j != omit]
        components.append(tensor_evaluate(cubic, selected))
    components.append(tensor_evaluate(quartic, columns))
    charts = [chart_p0_r0(poly) for poly in components]

    counts = [len(poly) for poly in charts]
    expected = [316, 316, 316, 316, 832]
    if counts != expected:
        raise AssertionError(f"unexpected chart term counts {counts}; expected {expected}")

    return {
        "format": "grunbaum-final-chart-polynomials-v1",
        "source_tensor_file": "final_tensors.json",
        "chart": "p=(1,x1,x2,x3), r=(1,y1,y2,y3)",
        "variables": ["x1", "x2", "x3", "y1", "y2", "y3"],
        "components": [
            "5A(q1,q2,q3)",
            "5A(q0,q2,q3)",
            "5A(q0,q1,q3)",
            "5A(q0,q1,q2)",
            "B(q0,q1,q2,q3)",
        ],
        "term_counts": counts,
        "polynomials": [encode(poly) for poly in charts],
    }


def canonical_json(data: object) -> str:
    return json.dumps(data, ensure_ascii=False, indent=2) + "\n"


def main() -> None:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tensors", type=Path, default=here.parent / "data" / "final_tensors.json")
    parser.add_argument("--output", type=Path, default=here.parent / "data" / "final_chart_polynomials.json")
    parser.add_argument("--check", action="store_true", help="fail if the committed output differs")
    args = parser.parse_args()

    rendered = canonical_json(generate(args.tensors))
    if args.check:
        current = args.output.read_text(encoding="utf-8")
        if current != rendered:
            raise SystemExit(f"ERROR: {args.output} is not reproduced by {args.tensors}")
        print(f"VERIFIED: {args.output} is reproduced exactly")
        return

    args.output.write_text(rendered, encoding="utf-8", newline="\n")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
