#!/usr/bin/env python3
"""Floating Bernstein subdivision utilities for the final chart polynomials.

Bernstein convex-hull bounds preserve polynomial cancellation far better than
naive interval arithmetic.  This script is still *not* a proof because it
uses binary floating point without recorded outward rounding.  If it closes
all boxes with a robust coefficient margin, the identical dyadic algorithm
can be rerun over integers/rationals to produce a checkable certificate.
"""

import json
import math
import argparse
from pathlib import Path

import numpy as np


DEG = 4


def power_to_bernstein_matrix(n=DEG):
    # p(x)=sum_e a_e x^e, x=2t-1; return b_i in degree-n Bernstein basis.
    t = np.zeros((n+1, n+1))
    for i in range(n+1):
        for e in range(n+1):
            z = 0.
            for k in range(min(i, e)+1):
                z += (math.comb(e, k)*2**k*(-1)**(e-k)
                      *math.comb(i, k)/math.comb(n, k))
            t[i, e] = z
    return t


T = power_to_bernstein_matrix()


def decode_power(raw):
    a = np.zeros((DEG+1,)*6)
    for m, c in raw:
        a[tuple(m)] += c
    return a


def transform_axis(a, axis):
    z = np.tensordot(T, a, axes=(1, axis))
    return np.moveaxis(z, 0, axis)


def to_bernstein(raw):
    a = decode_power(raw)
    for axis in range(6):
        a = transform_axis(a, axis)
    s = np.max(abs(a))
    return a/s


def split_half(a, axis):
    x = np.moveaxis(a, axis, -1)
    work = x.copy()
    left = np.empty_like(x)
    right = np.empty_like(x)
    left[..., 0] = work[..., 0]
    right[..., DEG] = work[..., DEG]
    for k in range(1, DEG+1):
        work[..., :DEG+1-k] = (
            work[..., :DEG+1-k]+work[..., 1:DEG+2-k])/2
        left[..., k] = work[..., 0]
        right[..., DEG-k] = work[..., DEG-k]
    return np.moveaxis(left, -1, axis), np.moveaxis(right, -1, axis)


def excluded(polys, tol=0.):
    for p in polys:
        if np.min(p) > tol or np.max(p) < -tol:
            return True
    return False


MIDWEIGHT = np.asarray([math.comb(DEG, i)/2**DEG for i in range(DEG+1)])


def midpoint_value(p):
    z = p
    # Contract one axis at a time; all axes use the same midpoint weights.
    for _ in range(6):
        z = np.tensordot(MIDWEIGHT, z, axes=(0, 0))
    return float(z)


def separating_margin(polys, iterations=4):
    """A valid vector-valued Bernstein separation lower bound.

    With all components represented in the same degree-4 tensor Bernstein
    basis, the image lies in the convex hull of its vector control points.
    The direction of the box-midpoint value is an inexpensive candidate
    separating functional; positive margin certifies that zero is absent
    even when every scalar component interval contains zero.
    """
    y = np.asarray([midpoint_value(p) for p in polys])
    n = np.linalg.norm(y)
    if n == 0:
        return -np.inf
    best = -np.inf
    for _ in range(iterations+1):
        combo = sum(z*p for z, p in zip(y, polys))
        flat_index = int(np.argmin(combo))
        mdot = float(combo.flat[flat_index])
        best = max(best, mdot/np.linalg.norm(y))
        # Gilbert/Frank-Wolfe step toward the most adverse control point.
        idx = np.unravel_index(flat_index, combo.shape)
        v = np.asarray([p[idx] for p in polys])
        d = v-y
        dd = float(np.dot(d, d))
        if dd == 0:
            break
        gamma = float(np.clip(-np.dot(y, d)/dd, 0., 1.))
        if gamma == 0:
            break
        y = y+gamma*d
        if np.linalg.norm(y) == 0:
            break
    return best


def split_score(polys, axis):
    # Bernstein adjacent differences are derivative control coefficients up
    # to the common degree; their current-box magnitude is a direct measure
    # of how much subdivision along this coordinate can shrink the hull.
    return sum(float(np.max(abs(np.diff(p, axis=axis)))) for p in polys)


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Run a non-rigorous floating Bernstein probe.")
    parser.add_argument("--source", type=Path,
                        default=here.parent / "data" / "final_chart_polynomials.json")
    parser.add_argument("--max-nodes", type=int, default=100000)
    parser.add_argument("--max-depth", type=int, default=60)
    parser.add_argument("--tolerance", type=float, default=1.0e-12)
    args = parser.parse_args()

    with args.source.open(encoding="utf-8") as f:
        polys = [to_bernstein(p) for p in json.load(f)["polynomials"]]
    stack = [(np.zeros(6, dtype=np.int16), polys)]
    nodes = removed = deep = 0
    worst_margin = np.inf
    while stack and nodes < args.max_nodes:
        depths, ps = stack.pop()
        nodes += 1
        margins = [max(float(np.min(p)), float(-np.max(p))) for p in ps]
        margin = max(max(margins), separating_margin(ps))
        if margin > args.tolerance:
            removed += 1
            worst_margin = min(worst_margin, margin)
            continue
        if int(np.sum(depths)) >= args.max_depth:
            deep += 1
            continue
        axis = int(np.argmax([split_score(ps, i) for i in range(6)]))
        lr = [split_half(p, axis) for p in ps]
        child_depth = depths.copy()
        child_depth[axis] += 1
        stack.append((child_depth.copy(), [z[1] for z in lr]))
        stack.append((child_depth, [z[0] for z in lr]))
        if nodes % 1000 == 0:
            print(nodes, "stack", len(stack), "removed", removed,
                  "deep", deep, "depth", int(np.sum(depths)), flush=True)
    print("FINAL nodes", nodes, "stack", len(stack), "removed", removed,
          "deep", deep, "worst_positive_margin", worst_margin,
          "tolerance", args.tolerance)


if __name__ == "__main__":
    main()
