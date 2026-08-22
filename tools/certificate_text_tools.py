#!/usr/bin/env python3
"""Convert and inspect the Grünbaum Bernstein subdivision certificate.

Text format (G4BCERT-TEXT-1):

    G4BCERT-TEXT-1 nodes=<N> leaves=<L> witness_scale=<S>
    split x1
    split y2
    leaf w0 w1 w2 w3 w4
    ...

Node records are in preorder, with the left child before the right child.
A split is always at the midpoint of the indicated coordinate.  Coordinate
order is x1,x2,x3,y1,y2,y3.  The root box is [0,1]^6 in Bernstein parameter
coordinates, corresponding to [-1,1]^6 under z -> 2z-1.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterable, Iterator, TextIO

BINARY_MAGIC = b"G4BCERT1"
TEXT_MAGIC = "G4BCERT-TEXT-1"
AXES = ("x1", "x2", "x3", "y1", "y2", "y3")
AXIS_TO_INDEX = {name: i for i, name in enumerate(AXES)}
LEAF_TOKEN = 6
EXPECTED_SCALE = 1 << 52


@dataclass(frozen=True)
class Certificate:
    tokens: bytes
    witnesses: tuple[tuple[int, int, int, int, int], ...]
    witness_scale: int

    @property
    def nodes(self) -> int:
        return len(self.tokens)

    @property
    def leaves(self) -> int:
        return len(self.witnesses)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_binary(path: Path) -> Certificate:
    data = path.read_bytes()
    if len(data) < 32 or data[:8] != BINARY_MAGIC:
        raise ValueError("bad binary certificate magic")
    nodes, leaves, scale = struct.unpack_from("<QQQ", data, 8)
    if nodes != 2 * leaves - 1:
        raise ValueError("binary header does not describe a full binary tree")
    token_start = 32
    token_end = token_start + nodes
    witness_end = token_end + leaves * 5 * 8
    if witness_end != len(data):
        raise ValueError("binary certificate length does not match its header")
    tokens = data[token_start:token_end]
    if any(token > LEAF_TOKEN for token in tokens):
        raise ValueError("binary certificate contains an invalid node token")
    witnesses = tuple(
        struct.unpack_from("<5q", data, token_end + 40 * i)
        for i in range(leaves)
    )
    cert = Certificate(tokens=tokens, witnesses=witnesses, witness_scale=scale)
    validate_tree(cert)
    return cert


def write_binary(cert: Certificate, path: Path) -> None:
    validate_tree(cert)
    with path.open("wb") as f:
        f.write(BINARY_MAGIC)
        f.write(struct.pack("<QQQ", cert.nodes, cert.leaves, cert.witness_scale))
        f.write(cert.tokens)
        for witness in cert.witnesses:
            f.write(struct.pack("<5q", *witness))


def _parse_header(line: str) -> tuple[int, int, int]:
    fields = line.strip().split()
    if len(fields) != 4 or fields[0] != TEXT_MAGIC:
        raise ValueError("bad text certificate header")
    values: dict[str, int] = {}
    for field in fields[1:]:
        if "=" not in field:
            raise ValueError(f"malformed header field: {field}")
        key, raw = field.split("=", 1)
        values[key] = int(raw)
    if set(values) != {"nodes", "leaves", "witness_scale"}:
        raise ValueError("text header must specify nodes, leaves, and witness_scale")
    return values["nodes"], values["leaves"], values["witness_scale"]


def read_text(path: Path) -> Certificate:
    with path.open("r", encoding="utf-8") as f:
        header = f.readline()
        if not header:
            raise ValueError("empty text certificate")
        nodes, leaves, scale = _parse_header(header)
        tokens = bytearray()
        witnesses: list[tuple[int, int, int, int, int]] = []
        for line_number, line in enumerate(f, start=2):
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            if fields[0] == "split":
                if len(fields) != 2 or fields[1] not in AXIS_TO_INDEX:
                    raise ValueError(f"line {line_number}: malformed split record")
                tokens.append(AXIS_TO_INDEX[fields[1]])
            elif fields[0] == "leaf":
                if len(fields) != 6:
                    raise ValueError(f"line {line_number}: leaf needs five weights")
                witness = tuple(int(x) for x in fields[1:])
                if not any(witness):
                    raise ValueError(f"line {line_number}: zero leaf witness")
                tokens.append(LEAF_TOKEN)
                witnesses.append(witness)  # type: ignore[arg-type]
            else:
                raise ValueError(f"line {line_number}: unknown record {fields[0]!r}")
    cert = Certificate(bytes(tokens), tuple(witnesses), scale)
    if cert.nodes != nodes or cert.leaves != leaves:
        raise ValueError(
            f"text header counts ({nodes}, {leaves}) do not match data "
            f"({cert.nodes}, {cert.leaves})"
        )
    validate_tree(cert)
    return cert


def write_text(cert: Certificate, path: Path) -> None:
    validate_tree(cert)
    leaf_index = 0
    with path.open("w", encoding="utf-8", newline="\n") as f:
        f.write(
            f"{TEXT_MAGIC} nodes={cert.nodes} leaves={cert.leaves} "
            f"witness_scale={cert.witness_scale}\n"
        )
        for token in cert.tokens:
            if token == LEAF_TOKEN:
                witness = cert.witnesses[leaf_index]
                leaf_index += 1
                f.write("leaf " + " ".join(str(x) for x in witness) + "\n")
            else:
                f.write(f"split {AXES[token]}\n")
    if leaf_index != cert.leaves:
        raise AssertionError("internal leaf count mismatch")


def validate_tree(cert: Certificate) -> dict[str, object]:
    if cert.nodes != 2 * cert.leaves - 1:
        raise ValueError("certificate is not a full binary tree by node count")
    if cert.witness_scale != EXPECTED_SCALE:
        raise ValueError(
            f"unexpected witness scale {cert.witness_scale}; expected {EXPECTED_SCALE}"
        )

    token_pos = 0
    leaf_pos = 0
    max_depth = 0
    axis_splits = [0] * 6
    depth_hist: dict[int, int] = {}

    # Iterative preorder parser. Each stack entry is a node depth still to read.
    stack = [0]
    while stack:
        depth = stack.pop()
        if token_pos >= cert.nodes:
            raise ValueError("tree ends before all expected children are present")
        token = cert.tokens[token_pos]
        token_pos += 1
        max_depth = max(max_depth, depth)
        if token == LEAF_TOKEN:
            if leaf_pos >= cert.leaves:
                raise ValueError("tree contains more leaves than witnesses")
            if not any(cert.witnesses[leaf_pos]):
                raise ValueError(f"leaf {leaf_pos} has the zero witness")
            leaf_pos += 1
            depth_hist[depth] = depth_hist.get(depth, 0) + 1
        elif 0 <= token < 6:
            axis_splits[token] += 1
            # Push right then left, so the left child is consumed next.
            stack.append(depth + 1)
            stack.append(depth + 1)
        else:
            raise ValueError(f"invalid token {token} at preorder node {token_pos - 1}")

    if token_pos != cert.nodes:
        raise ValueError("tree has trailing node records after the root is complete")
    if leaf_pos != cert.leaves:
        raise ValueError("number of leaf records does not match witnesses")

    return {
        "nodes": cert.nodes,
        "leaves": cert.leaves,
        "internal_nodes": cert.nodes - cert.leaves,
        "max_depth": max_depth,
        "axis_splits": dict(zip(AXES, axis_splits)),
        "leaf_depth_histogram": dict(sorted(depth_hist.items())),
    }


def print_summary(cert: Certificate, path: Path | None = None) -> None:
    summary = validate_tree(cert)
    if path is not None:
        print(f"file: {path}")
        print(f"sha256: {sha256(path)}")
    print(f"format: {TEXT_MAGIC if path and path.suffix != '.bin' else 'G4BCERT1'}")
    print(f"nodes: {summary['nodes']}")
    print(f"leaves: {summary['leaves']}")
    print(f"internal nodes: {summary['internal_nodes']}")
    print(f"max depth: {summary['max_depth']}")
    print(f"witness scale: {cert.witness_scale}")
    print("splits by coordinate:")
    for axis, count in summary["axis_splits"].items():  # type: ignore[union-attr]
        print(f"  {axis}: {count}")
    print("leaf depths:")
    histogram = summary["leaf_depth_histogram"]  # type: ignore[assignment]
    for depth, count in histogram.items():  # type: ignore[union-attr]
        print(f"  {depth}: {count}")



def write_leaf_csv(cert: Certificate, path: Path) -> None:
    """Write explicit dyadic leaf boxes and witnesses as CSV.

    For each coordinate, (k,d) denotes the interval
    [k/2^d,(k+1)/2^d] in the Bernstein parameter cube [0,1]^6.
    """
    validate_tree(cert)
    token_pos = 0
    leaf_pos = 0
    node_id = 0
    volume_numerators: dict[int, int] = {}

    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        header = ["leaf_id", "node_id", "depth"]
        for axis in AXES:
            header.extend([f"{axis}_k", f"{axis}_d"])
        header.extend(["w0", "w1", "w2", "w3", "w4"])
        writer.writerow(header)

        # Stack entries: (depth, cells), cells are six (k,d) pairs.
        stack: list[tuple[int, tuple[tuple[int, int], ...]]] = [
            (0, tuple((0, 0) for _ in AXES))
        ]
        while stack:
            depth, cells = stack.pop()
            if token_pos >= cert.nodes:
                raise ValueError("tree ended while producing leaf CSV")
            current_node_id = node_id
            node_id += 1
            token = cert.tokens[token_pos]
            token_pos += 1
            if token == LEAF_TOKEN:
                witness = cert.witnesses[leaf_pos]
                row: list[int] = [leaf_pos, current_node_id, depth]
                for k, d in cells:
                    row.extend([k, d])
                row.extend(witness)
                writer.writerow(row)
                total_depth = sum(d for _, d in cells)
                volume_numerators[total_depth] = volume_numerators.get(total_depth, 0) + 1
                leaf_pos += 1
            else:
                axis = token
                left = list(cells)
                right = list(cells)
                k, d = cells[axis]
                left[axis] = (2 * k, d + 1)
                right[axis] = (2 * k + 1, d + 1)
                # Preorder: process left first.
                stack.append((depth + 1, tuple(right)))
                stack.append((depth + 1, tuple(left)))

    if token_pos != cert.nodes or leaf_pos != cert.leaves:
        raise ValueError("tree was not consumed exactly while producing leaf CSV")

    # Exact cover-volume check using a common dyadic denominator.
    max_exp = max(volume_numerators, default=0)
    scaled_sum = sum(count << (max_exp - exp) for exp, count in volume_numerators.items())
    if scaled_sum != 1 << max_exp:
        raise ValueError("leaf boxes do not have total volume one")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    export = sub.add_parser("export", help="convert binary certificate to text")
    export.add_argument("binary", type=Path)
    export.add_argument("text", type=Path)

    imp = sub.add_parser("import", help="convert text certificate to binary")
    imp.add_argument("text", type=Path)
    imp.add_argument("binary", type=Path)

    inspect = sub.add_parser("inspect", help="validate and summarize a certificate")
    inspect.add_argument("certificate", type=Path)

    roundtrip = sub.add_parser("roundtrip", help="prove text conversion is lossless")
    roundtrip.add_argument("binary", type=Path)
    roundtrip.add_argument("text", type=Path)
    roundtrip.add_argument("rebuilt_binary", type=Path)

    leaves = sub.add_parser("leaves", help="write explicit dyadic leaf boxes as CSV")
    leaves.add_argument("certificate", type=Path)
    leaves.add_argument("csv", type=Path)

    args = parser.parse_args()

    if args.command == "export":
        cert = read_binary(args.binary)
        write_text(cert, args.text)
        print_summary(cert, args.text)
    elif args.command == "import":
        cert = read_text(args.text)
        write_binary(cert, args.binary)
        print_summary(cert, args.binary)
    elif args.command == "inspect":
        path: Path = args.certificate
        cert = read_binary(path) if path.suffix == ".bin" else read_text(path)
        print_summary(cert, path)
    elif args.command == "leaves":
        path: Path = args.certificate
        cert = read_binary(path) if path.suffix == ".bin" else read_text(path)
        write_leaf_csv(cert, args.csv)
        print(f"wrote {cert.leaves} leaf boxes to {args.csv}")
        print(f"sha256: {sha256(args.csv)}")
    elif args.command == "roundtrip":
        original = read_binary(args.binary)
        write_text(original, args.text)
        parsed = read_text(args.text)
        if parsed != original:
            raise SystemExit("ERROR: text representation changed the certificate")
        write_binary(parsed, args.rebuilt_binary)
        original_hash = sha256(args.binary)
        rebuilt_hash = sha256(args.rebuilt_binary)
        print(f"original sha256: {original_hash}")
        print(f"rebuilt  sha256: {rebuilt_hash}")
        if original_hash != rebuilt_hash:
            raise SystemExit("ERROR: rebuilt binary differs from original")
        print("ROUNDTRIP VERIFIED: text certificate is lossless")


if __name__ == "__main__":
    main()
