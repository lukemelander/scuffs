#!/usr/bin/env python3
"""Compare flat-plate boundary-layer profiles from a SCUFFS solution with the
Blasius solution.

    tools/blasius.py examples/plate_viscous_implicit/data.vtk \\
        --input examples/plate_viscous_implicit/input.inp --x 0.3 0.6 0.9

For each station x it takes the column of cells closest to x and prints u/U
against two similarity variables, x' = x - x_le being the distance from the
leading edge:

  eta   = y sqrt(U / (nu_inf x'))                   (plain Blasius scaling)
  eta_H = sqrt(U / (nu_inf x')) * int_0^y rho/rho_inf dy   (Howarth transform)

The Howarth transform accounts for the lower density near a hot (adiabatic)
wall in compressible flow; with constant viscosity it is not exact, so
expect the CFD to sit between the two comparisons at Mach 1.5.

Freestream values come from the input file (p_inf, rho_inf, M_inf); the
viscosity is Sutherland's law at the freestream temperature, as in SCUFFS.
"""
import argparse
import math
import sys

GAM, RAIR = 1.4, 287.05


def blasius(eta_max=10.0, n=4000):
    """f'(eta) of the Blasius solution (f''' + f f''/2 = 0) on a uniform grid."""
    h = eta_max / n
    y = [0.0, 0.0, 0.332057336215]  # f, f', f''
    etas, fp = [0.0], [0.0]

    def rhs(s):
        return [s[1], s[2], -0.5 * s[0] * s[2]]

    for i in range(n):
        k1 = rhs(y)
        k2 = rhs([y[j] + 0.5 * h * k1[j] for j in range(3)])
        k3 = rhs([y[j] + 0.5 * h * k2[j] for j in range(3)])
        k4 = rhs([y[j] + h * k3[j] for j in range(3)])
        y = [y[j] + h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) for j in range(3)]
        etas.append((i + 1) * h)
        fp.append(y[1])
    return etas, fp


def interp(x, xs, ys):
    if x >= xs[-1]:
        return ys[-1]
    i = max(0, min(int(x / (xs[1] - xs[0])), len(xs) - 2))
    t = (x - xs[i]) / (xs[i + 1] - xs[i])
    return ys[i] + t * (ys[i + 1] - ys[i])


def read_vtk(path):
    with open(path) as fh:
        tok = fh.read().split()
    k = tok.index("POINTS")
    npts = int(tok[k + 1])
    pts = [(float(tok[k + 3 + 3 * i]), float(tok[k + 4 + 3 * i])) for i in range(npts)]
    k = tok.index("CELLS")
    ncell = int(tok[k + 1])
    pos = k + 3
    centres = []
    for _ in range(ncell):
        nv = int(tok[pos])
        ids = [int(t) for t in tok[pos + 1:pos + 1 + nv]]
        centres.append((sum(pts[i][0] for i in ids) / nv, sum(pts[i][1] for i in ids) / nv))
        pos += nv + 1
    fields = {}
    k = 0
    while True:
        try:
            k = tok.index("SCALARS", k)
        except ValueError:
            break
        name = tok[k + 1]
        start = k + 6
        fields[name] = [float(t) for t in tok[start:start + ncell]]
        k = start
    return centres, fields


def read_freestream(path):
    with open(path) as fh:
        lines = [l for l in fh if l.strip()]
    p, rho, mach = (float(v.replace("d", "e").replace("D", "e")) for v in lines[1].split()[:3])
    T = p / (rho * RAIR)
    mu = 1.458e-6 * T ** 1.5 / (T + 110.4)
    U = mach * math.sqrt(GAM * p / rho)
    return rho, U, mu


def main():
    ap = argparse.ArgumentParser(description="Compare plate boundary-layer profiles with Blasius.")
    ap.add_argument("vtk", help="solution file (data.vtk)")
    ap.add_argument("--input", required=True, help="input file of the run (for the freestream)")
    ap.add_argument("--x", type=float, nargs="+", default=[0.3, 0.6, 0.9], help="stations")
    ap.add_argument("--xle", type=float, default=0.09, help="leading-edge x (default 0.09)")
    ap.add_argument("--csv", help="also write the profiles to this CSV file")
    args = ap.parse_args()

    rho_inf, U, mu = read_freestream(args.input)
    nu = mu / rho_inf
    centres, f = read_vtk(args.vtk)
    for name in ("u", "rho"):
        if name not in f:
            sys.exit(f"{args.vtk}: no '{name}' field")
    etas, fp = blasius()
    rows = []

    print(f"Freestream: U = {U:.2f} m/s, nu = {nu:.4e} m^2/s, leading edge x = {args.xle}")
    for xs in args.x:
        xp = xs - args.xle
        if xp <= 0:
            print(f"x = {xs}: upstream of the leading edge, skipped")
            continue
        # column of cells closest to the station
        xc = min((c[0] for c in centres), key=lambda v: abs(v - xs))
        col = sorted((c[1], f["u"][i], f["rho"][i]) for i, c in enumerate(centres)
                     if abs(c[0] - xc) < 1e-9 * max(1.0, abs(xc)))
        scale = math.sqrt(U / (nu * xp))
        print(f"\nx = {xc:.4f} (x' = {xp:.4f}, Re_x = {U * xp / nu:.3e}), {len(col)} cells")
        print(f"{'eta':>8} {'u/U':>8} {'Blasius':>8} | {'eta_H':>8} {'Blasius':>8}")
        int_rho, yprev, rprev = 0.0, 0.0, col[0][2] if col else 1.0
        err_plain, err_howarth = 0.0, 0.0
        for y, u, r in col:
            # Howarth: integrate rho/rho_inf from the wall to the cell centre
            int_rho += 0.5 * (r + rprev) / rho_inf * (y - yprev) if yprev > 0 else r / rho_inf * y
            yprev, rprev = y, r
            eta, eta_h = y * scale, int_rho * scale
            if eta > 9:
                break
            ub, ubh = interp(eta, etas, fp), interp(eta_h, etas, fp)
            err_plain = max(err_plain, abs(u / U - ub))
            err_howarth = max(err_howarth, abs(u / U - ubh))
            print(f"{eta:8.3f} {u / U:8.4f} {ub:8.4f} | {eta_h:8.3f} {ubh:8.4f}")
            rows.append((xs, y, eta, eta_h, u / U, ub, ubh))
        print(f"max |u/U - Blasius|: plain scaling {err_plain:.3f}, Howarth {err_howarth:.3f}")

    if args.csv:
        with open(args.csv, "w") as out:
            out.write("x,y,eta,eta_howarth,u_over_U,blasius_eta,blasius_eta_howarth\n")
            for r in rows:
                out.write(",".join(f"{v:.6g}" for v in r) + "\n")
        print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()
