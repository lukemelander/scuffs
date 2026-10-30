#!/usr/bin/env python3
"""Compare an unstructured SCUFFS solution (data.vtk) with a structured v1.0
solution (data.dat, Tecplot block format) on the same grid.

The structured-to-SU2 converter numbers cells in the same order v1 writes
them (j outer, i inner), so the two files line up cell by cell.

    tools/compare_v1.py v1/data.dat data.vtk
"""
import re
import sys


def read_v1(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    for k, line in enumerate(lines):
        if "ZONE" in line.upper():
            break
    else:
        sys.exit(f"{path}: no ZONE line")
    nx = int(re.search(r"\bI\s*=\s*(\d+)", lines[k], re.I).group(1))
    ny = int(re.search(r"\bJ\s*=\s*(\d+)", lines[k], re.I).group(1))
    vals = [float(v) for v in " ".join(lines[k + 1:]).split()]
    nn, nc = nx * ny, (nx - 1) * (ny - 1)
    out, pos = {}, 2 * nn
    for name in ("rho", "u", "v", "p", "mach"):
        out[name] = vals[pos:pos + nc]
        pos += nc
    return out


def read_vtk(path):
    with open(path) as fh:
        tokens = fh.read().split()
    out, ncell = {}, None
    k = 0
    while k < len(tokens):
        if tokens[k] == "CELL_DATA":
            ncell = int(tokens[k + 1])
        if tokens[k] == "SCALARS":
            name = tokens[k + 1]
            # SCALARS name type ncomp LOOKUP_TABLE default <values>
            start = k + 6
            out[name] = [float(v) for v in tokens[start:start + ncell]]
            k = start + ncell
            continue
        k += 1
    return out


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    a, b = read_v1(sys.argv[1]), read_vtk(sys.argv[2])
    if len(a["rho"]) != len(b["rho"]):
        sys.exit(f"different cell counts: {len(a['rho'])} vs {len(b['rho'])}")
    print(f"{'var':6s} {'max |diff|':>12s} {'max rel diff':>14s}  (relative to max |value|)")
    for name in ("rho", "u", "v", "p", "mach"):
        x, y = a[name], b[name]
        scale = max(abs(v) for v in x) or 1.0
        d = max(abs(p - q) for p, q in zip(x, y))
        print(f"{name:6s} {d:12.3e} {d / scale:14.3e}")


if __name__ == "__main__":
    main()
