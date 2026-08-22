#!/usr/bin/env python3
"""Generate a serialized Bernstein subdivision certificate.

This is the untrusted producer. It may use floating point to choose split
axes and candidate separating directions. Every emitted leaf direction is
quantized to a common signed-integer scale. The separate C checker derives
the five chart polynomials from the displayed tensors and verifies every
subdivision and every leaf inequality using integer arithmetic only.
"""
from __future__ import annotations

import argparse
import json
import math
import struct
import time
from pathlib import Path

import numpy as np

from bernstein_tools import T, to_bernstein, split_half, midpoint_value, split_score

S = 1 << 52
EPS = np.finfo(float).eps
WEIGHT_SCALE = 1 << 52
LEAF = 6
MAGIC = b"G4BCERT1"


def exact_root(raw):
    t6 = np.rint(6 * T).astype(np.int64)
    assert np.max(abs(6 * T - t6)) == 0
    roots = []
    for poly in raw:
        a = np.zeros((5,) * 6, dtype=np.int64)
        for m, c in poly:
            a[tuple(m)] += int(c)
        for axis in range(6):
            a = np.moveaxis(np.tensordot(t6, a, axes=(1, axis)), 0, axis)
        roots.append(a)
    mids, los, his = [], [], []
    for a, poly in zip(roots, raw):
        den = int(np.max(np.abs(a)))
        assert den > 0
        flat = a.ravel()
        lo = np.fromiter(((int(v) * S) // den for v in flat),
                         dtype=np.int64, count=flat.size)
        hi = np.fromiter((-((-int(v) * S) // den) for v in flat),
                         dtype=np.int64, count=flat.size)
        los.append(lo.reshape(a.shape))
        his.append(hi.reshape(a.shape))
        mids.append(to_bernstein(poly))
    return mids, los, his


L = np.zeros((5, 5), dtype=np.int64)
R = np.zeros((5, 5), dtype=np.int64)
for i in range(5):
    for k in range(i + 1):
        L[i, k] = (1 << (4 - i)) * math.comb(i, k)
for i in range(5):
    for k in range(i, 5):
        R[i, k] = (1 << i) * math.comb(4 - i, k - i)
assert np.all(L.sum(1) == 16) and np.all(R.sum(1) == 16)


def divlo(x):
    return x // 16


def divhi(x):
    return -((-x) // 16)


def split_enclosure(arr, axis, lower):
    x = np.moveaxis(arr, axis, -1)
    assert int(np.max(np.abs(x))) < (1 << 58)
    left = np.tensordot(x, L.T, axes=(-1, 0))
    right = np.tensordot(x, R.T, axes=(-1, 0))
    op = divlo if lower else divhi
    return np.moveaxis(op(left), -1, axis), np.moveaxis(op(right), -1, axis)


def candidate_direction(polys, iterations=12):
    for i, p in enumerate(polys):
        if float(np.min(p)) > 0:
            u = np.zeros(5)
            u[i] = 1
            return float(np.min(p)), u
        if float(np.max(p)) < 0:
            u = np.zeros(5)
            u[i] = -1
            return float(-np.max(p)), u
    y = np.asarray([midpoint_value(p) for p in polys])
    best = (-np.inf, None)
    if np.linalg.norm(y) == 0:
        return best
    for _ in range(iterations + 1):
        combo = sum(z * p for z, p in zip(y, polys))
        idx = int(np.argmin(combo))
        mdot = float(combo.flat[idx])
        margin = mdot / np.linalg.norm(y)
        if margin > best[0]:
            best = (margin, y.copy())
        ix = np.unravel_index(idx, combo.shape)
        v = np.asarray([p[ix] for p in polys])
        d = v - y
        dd = float(d @ d)
        if dd == 0:
            break
        gamma = float(np.clip(-(y @ d) / dd, 0, 1))
        if gamma == 0:
            break
        y = y + gamma * d
        if np.linalg.norm(y) == 0:
            break
    return best


def rigorous_float_lower(u, los, his):
    terms = []
    for z, lo, hi in zip(u, los, his):
        a = lo if z >= 0 else hi
        terms.append(float(z) * (a.astype(np.float64) / S))
    val = sum(terms)
    ab = sum(abs(x) for x in terms)
    lower = val - 64 * EPS * ab - np.nextafter(0.0, 1.0)
    return float(np.min(lower))


def quantize_direction(u):
    m = float(np.max(np.abs(u)))
    if not (m > 0):
        raise ValueError("zero separating direction")
    w = np.rint((u / m) * WEIGHT_SCALE).astype(np.int64)
    if not np.any(w):
        raise ValueError("direction vanished during quantization")
    nz = [abs(int(x)) for x in w if x]
    shift = min((x & -x).bit_length() - 1 for x in nz)
    if shift:
        w = w >> shift
    return tuple(int(x) for x in w)


def longdouble_lower(weights, los, his):
    total = None
    for w, lo, hi in zip(weights, los, his):
        a = lo if w >= 0 else hi
        term = np.longdouble(w) * a.astype(np.longdouble)
        total = term if total is None else total + term
    return float(np.min(total))


def write_certificate(path: Path, tokens: bytearray, witnesses):
    with path.open("wb") as f:
        f.write(MAGIC)
        f.write(struct.pack("<QQQ", len(tokens), len(witnesses), WEIGHT_SCALE))
        f.write(tokens)
        for w in witnesses:
            f.write(struct.pack("<5q", *w))


def run(source: Path, output: Path, maxnodes: int, maxdepth: int, report: int):
    raw = json.loads(source.read_text(encoding="utf-8"))["polynomials"]
    ps, lo, hi = exact_root(raw)
    stack = [(np.zeros(6, dtype=np.int16), ps, lo, hi)]
    tokens = bytearray()
    witnesses = []
    nodes = leaves = deep = failed = scalar = 0
    worst_float = float("inf")
    worst_integer_ld = float("inf")
    started = time.time()

    while stack and nodes < maxnodes:
        ds, p, l, h = stack.pop()
        nodes += 1
        margin, u = candidate_direction(p, iterations=12)
        certified = False
        weights = None
        em = None

        if int(ds.sum()) >= 60:
            scalar_margins = []
            for i in range(5):
                scalar_margins.append(max(int(np.min(l[i])), -int(np.max(h[i]))))
            i = int(np.argmax(scalar_margins))
            if scalar_margins[i] > 0:
                sign = 1 if int(np.min(l[i])) > 0 else -1
                weights = [0] * 5
                weights[i] = sign
                weights = tuple(weights)
                em = scalar_margins[i] / S
                certified = True
                scalar += 1

        if margin > 0 and not certified:
            em = rigorous_float_lower(u, l, h)
            if em > 0:
                weights = quantize_direction(u)
                qlower = longdouble_lower(weights, l, h)
                if qlower > 0:
                    certified = True
                    worst_integer_ld = min(worst_integer_ld, qlower)

        if not certified and margin > 0:
            margin, u = candidate_direction(p, iterations=48)
            if margin > 0:
                em = rigorous_float_lower(u, l, h)
                if em > 0:
                    weights = quantize_direction(u)
                    qlower = longdouble_lower(weights, l, h)
                    if qlower > 0:
                        certified = True
                        worst_integer_ld = min(worst_integer_ld, qlower)

        if certified:
            tokens.append(LEAF)
            witnesses.append(weights)
            leaves += 1
            worst_float = min(worst_float, float(em))
            continue

        if margin > 0:
            failed += 1
        if int(ds.sum()) >= maxdepth:
            deep += 1
            raise RuntimeError(f"unresolved box at depth vector {ds.tolist()}")

        axis = int(np.argmax([split_score(p, i) for i in range(6)]))
        tokens.append(axis)
        fp = [split_half(z, axis) for z in p]
        il = [split_enclosure(z, axis, True) for z in l]
        ih = [split_enclosure(z, axis, False) for z in h]
        nd = ds.copy()
        nd[axis] += 1
        stack.append((nd.copy(), [z[1] for z in fp], [z[1] for z in il], [z[1] for z in ih]))
        stack.append((nd, [z[0] for z in fp], [z[0] for z in il], [z[0] for z in ih]))

        if report and nodes % report == 0:
            print("PROGRESS", nodes, "stack", len(stack), "leaves", leaves,
                  "failed", failed, "scalar", scalar, "depth", int(ds.sum()),
                  "worst_float", worst_float,
                  "sec", round(time.time() - started, 1), flush=True)

    if stack or deep:
        raise RuntimeError("certificate generation incomplete")
    if nodes != 2 * leaves - 1:
        raise RuntimeError("tree is not full binary")
    write_certificate(output, tokens, witnesses)
    print("FINAL nodes", nodes, "leaves", leaves, "scalar", scalar,
          "failed", failed, "worst_float", worst_float,
          "worst_integer_longdouble", worst_integer_ld,
          "bytes", output.stat().st_size,
          "seconds", time.time() - started, flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    here = Path(__file__).resolve().parent
    ap.add_argument("--source", type=Path, default=here.parent / "data" / "final_chart_polynomials.json")
    ap.add_argument("--output", type=Path, default=here.parent / "certificate" / "grunbaum_certificate.generated.bin")
    ap.add_argument("--maxnodes", type=int, default=1_000_000)
    ap.add_argument("--maxdepth", type=int, default=120)
    ap.add_argument("--report", type=int, default=10_000)
    args = ap.parse_args()
    run(args.source, args.output, args.maxnodes, args.maxdepth, args.report)
