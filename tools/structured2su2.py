#!/usr/bin/env python3
"""Convert a structured SCUFFS grid (Tecplot ASCII, nodes only) to SU2.

The structured grid is nx by ny nodes, written in block format with i
varying fastest:

    <title line>
    ZONE I= <nx> J= <ny> ...
    <nx*ny x-coordinates>
    <nx*ny y-coordinates>

Each structured cell becomes a quad. The four sides of the grid become
boundary markers, named with --imin/--imax/--jmin/--jmax. Sides given the
same name are merged into one marker. A side can be split into two markers
at a node index with --split, e.g. for the flat plate:

    tools/structured2su2.py grids/plate.dat grids/plate.su2 \\
        --imin inlet --imax outlet --jmax top \\
        --split jmin:10:symmetry:plate

Node and edge numbering in the .su2 file is 0-based, as SU2 expects.
"""
import argparse
import re
import sys

SIDES = ("imin", "imax", "jmin", "jmax")


def read_structured(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    zone_line = None
    for k, line in enumerate(lines):
        if "ZONE" in line.upper():
            zone_line = k
            break
    if zone_line is None:
        sys.exit(f"{path}: no ZONE line found")
    zone = lines[zone_line]
    mi = re.search(r"\bI\s*=\s*(\d+)", zone, re.IGNORECASE)
    mj = re.search(r"\bJ\s*=\s*(\d+)", zone, re.IGNORECASE)
    if not (mi and mj):
        sys.exit(f"{path}: could not read I= and J= from the ZONE line")
    nx, ny = int(mi.group(1)), int(mj.group(1))
    values = " ".join(lines[zone_line + 1:]).split()
    nn = nx * ny
    if len(values) < 2 * nn:
        sys.exit(f"{path}: expected {2 * nn} coordinates, found {len(values)}")
    x = [float(v) for v in values[:nn]]
    y = [float(v) for v in values[nn:2 * nn]]
    return nx, ny, x, y


def side_edges(side, nx, ny):
    """Boundary edges of one side as (position along the side, node a, node b),
    1-based structured indices; position k means the edge from node k to k+1."""
    nid = lambda i, j: (j - 1) * nx + (i - 1)          # 0-based node id
    if side == "imin":
        return [(j, nid(1, j), nid(1, j + 1)) for j in range(1, ny)]
    if side == "imax":
        return [(j, nid(nx, j), nid(nx, j + 1)) for j in range(1, ny)]
    if side == "jmin":
        return [(i, nid(i, 1), nid(i + 1, 1)) for i in range(1, nx)]
    if side == "jmax":
        return [(i, nid(i, ny), nid(i + 1, ny)) for i in range(1, nx)]
    raise ValueError(side)


def parse_split(text):
    try:
        side, node, before, after = text.split(":")
        node = int(node)
    except ValueError:
        sys.exit(f"--split {text}: expected SIDE:NODE:NAME_BEFORE:NAME_AFTER")
    if side not in SIDES:
        sys.exit(f"--split {text}: side must be one of {', '.join(SIDES)}")
    return side, node, before, after


def main():
    ap = argparse.ArgumentParser(
        description="Convert a structured SCUFFS grid (.dat) to an SU2 mesh (.su2).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example (flat plate):\n"
               "  tools/structured2su2.py grids/plate.dat grids/plate.su2 \\\n"
               "      --imin inlet --imax outlet --jmax top --split jmin:10:symmetry:plate")
    ap.add_argument("input", help="structured grid (.dat)")
    ap.add_argument("output", help="SU2 mesh to write (.su2)")
    for side in SIDES:
        ap.add_argument(f"--{side}", default=side, metavar="NAME",
                        help=f"marker name for the {side} side (default: {side})")
    ap.add_argument("--triangulate", action="store_true",
                    help="split every quad into two triangles, alternating the diagonal in a "
                         "checkerboard pattern (for testing on triangle meshes)")
    ap.add_argument("--split", action="append", default=[], metavar="SIDE:NODE:BEFORE:AFTER",
                    help="split a side into two markers: edges starting at node index < NODE "
                         "(1-based, along the side) get BEFORE, the rest get AFTER. Can be repeated.")
    args = ap.parse_args()

    nx, ny, x, y = read_structured(args.input)

    splits = {}
    for text in args.split:
        side, node, before, after = parse_split(text)
        length = (ny if side in ("imin", "imax") else nx)
        if not 2 <= node <= length - 1:
            sys.exit(f"--split {text}: NODE must be between 2 and {length - 1} for {side}")
        splits[side] = (node, before, after)

    # quads, counter-clockwise in (i, j) index space
    nid = lambda i, j: (j - 1) * nx + (i - 1)
    cells = []
    for j in range(1, ny):
        for i in range(1, nx):
            a, b, c, d = nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)
            if not args.triangulate:
                cells.append((a, b, c, d))
            elif (i + j) % 2 == 0:          # diagonal a-c
                cells += [(a, b, c), (a, c, d)]
            else:                           # diagonal b-d
                cells += [(a, b, d), (b, c, d)]

    # boundary edges grouped by marker name (in first-seen order)
    markers = {}
    for side in SIDES:
        for k, a, b in side_edges(side, nx, ny):
            if side in splits:
                node, before, after = splits[side]
                name = before if k < node else after
            else:
                name = getattr(args, side)
            markers.setdefault(name, []).append((a, b))

    with open(args.output, "w") as out:
        out.write(f"% Converted from {args.input} ({nx} x {ny} nodes) by structured2su2.py\n")
        out.write("NDIME= 2\n")
        out.write(f"NELEM= {len(cells)}\n")
        for e, cell in enumerate(cells):
            etype = 9 if len(cell) == 4 else 5
            out.write(f"{etype} {' '.join(map(str, cell))} {e}\n")
        out.write(f"NPOIN= {nx * ny}\n")
        for k in range(nx * ny):
            out.write(f"{x[k]:.15e} {y[k]:.15e} {k}\n")
        out.write(f"NMARK= {len(markers)}\n")
        for name, edges in markers.items():
            out.write(f"MARKER_TAG= {name}\n")
            out.write(f"MARKER_ELEMS= {len(edges)}\n")
            for a, b in edges:
                out.write(f"3 {a} {b}\n")

    summary = ", ".join(f"{name} ({len(edges)})" for name, edges in markers.items())
    kind = "triangles" if args.triangulate else "quads"
    print(f"{args.input} -> {args.output}: {nx} x {ny} nodes, {len(cells)} {kind}, markers: {summary}")


if __name__ == "__main__":
    main()
