# fv2d: 2D finite-volume solver for compressible flow

A first-order, cell-centred finite-volume solver for the 2D Euler and
thin-layer Navier–Stokes equations on structured (i, j) grids, written in
Fortran 90. It was written for a graduate course on finite-volume methods for
hypersonic CFD and later cleaned up and debugged.

- Inviscid fluxes: modified Steger–Warming (MSW) flux-vector splitting with
  pressure-gradient blending
- Viscous fluxes: thin-layer Navier–Stokes in the j direction (Sutherland's
  law for viscosity)
- Time integration: explicit (forward Euler), or implicit (backward Euler with
  point-Jacobi relaxation) with an optional CFL ramp for steady-state runs
- Boundary conditions: supersonic inflow, characteristic outflow, slip wall,
  mirror, adiabatic no-slip wall
- Output: Tecplot ASCII files that also open in ParaView

## Repository layout

```
.
├── Makefile        build rules (gfortran by default)
├── run.sh          runs an example case
├── src/            Fortran sources (main.f90, solver.f90, util.f90)
├── grids/          structured grids (Tecplot ASCII, nodes only)
└── examples/       one directory per case: input.inp + case.conf
```

## Building

Requires a Fortran compiler (gfortran by default) and `make`.

```
make            # optimised build -> bin/fv2d
make debug      # -O0 with bounds checking and backtraces
make clean
make FC=ifx FFLAGS=-O3   # different compiler or flags
```

## Running

Run an example from the top level:

```
./run.sh examples/cylinder_m6_implicit
```

`run.sh` builds the solver if needed, runs it in the case directory, and
writes `data.dat`, `residual.dat` and `run.log` there. Extra arguments are
passed on to the solver.

The solver can also be run directly:

```
bin/fv2d --input path/to/input.inp --grid path/to/grid.dat
```

| Option | Meaning |
|---|---|
| `-i`, `--input FILE` | Input file (default `input.inp`) |
| `-g`, `--grid FILE` | Grid file (default: picked from the case type on line 3 of the input file, e.g. 4 -> `plate.dat`, looked up in the current directory) |
| `-h`, `--help` | Show usage |

`--input=FILE` and `--grid=FILE` also work.

## Examples

| Case | Grid | Description |
|---|---|---|
| `ramp_explicit` | `ramp.dat` | Mach 2.5 inviscid flow in a channel with a ramp, explicit |
| `ramp_implicit` | `ramp.dat` | Same, implicit with CFL ramp |
| `ramp_aligned_implicit` | `rampb.dat` | Same, on a grid aligned with the ramp shock |
| `cylinder_m6_implicit` | `cyl.dat` | Mach 6 inviscid flow over a cylinder (bow shock), implicit |
| `plate_viscous_implicit` | `plate.dat` | Mach 1.5 laminar flat-plate boundary layer, implicit |
| `plate_viscous_explicit` | `plate.dat` | Same, explicit (slow: the time step is set by the near-wall viscous limit) |
| `shock_tube` | `block.dat` | Shock tube in a closed box (gas at rest, lower half at 2x pressure and density) |

Each case directory has `input.inp` (solver settings) and `case.conf`, which
names the grid file in `grids/` that `run.sh` passes with `--grid`.

## Input file

The input file has six lines, read with list-directed reads, so anything after
the values on a line is treated as a comment.

```
1 1 1              !ivisc, impl, cfl_ramp
1000.0 0.0071 1.5  !p_inf rho_inf M_inf
4 5000 500         !case type, iterations, save interval
0.4d0              !CFL
2,1,4,5            !BCs at i=nx, i=1, j=ny, j=1
1                  !initial condition
```

| Line | Values | Meaning |
|---|---|---|
| 1 | `ivisc` | 0 = inviscid (Euler), 1 = viscous (thin-layer NS, j direction only) |
|   | `impl` | 0 = explicit, 1 = implicit (point Jacobi, `kmax`+1 = 10 sweeps, set in `src/util.f90`) |
|   | `cfl_ramp` | 0 = constant CFL, 1 = built-in CFL ramp (see below) |
| 2 | `p_inf rho_inf M_inf` | Freestream pressure [Pa], density [kg/m^3] and Mach number. The flow is in +x. Viscosity uses Sutherland's law at the freestream temperature and is held constant |
| 3 | `stype` | Case type (also picks the default grid): 0 = block, 1 = ramp, 2 = shock-aligned ramp, 3 = cylinder, 4 = flat plate |
|   | `iter` | Number of iterations |
|   | `isave` | Write `data.dat` every `isave` iterations (it is also written at the end of the run) |
| 4 | `CFL` | CFL number (starting value if the ramp is on). Explicit viscous runs also include the viscous stability limit |
| 5 | `bc(1) bc(2) bc(3) bc(4)` | Boundary type at i=nx (right), i=1 (left), j=ny (top), j=1 (bottom) |
| 6 | `switch` | Initial condition: 1 = uniform freestream, 2 = shock tube (lower half of j at 2x p and rho, at rest), 3 = freestream with a 1% bump at i=(nx-1)/2 |

Boundary types:

| Code | Type |
|---|---|
| 1 | Supersonic inflow (ghost cell held at freestream) |
| 2 | Outflow (characteristic; only outgoing waves are taken from the interior) |
| 3 | Inviscid (slip) wall / symmetry plane |
| 4 | Mirror (ghost cell = interior cell) |
| 5 | Adiabatic no-slip wall (intended for j=1; its inviscid flux is a slip wall) |

The case type is still read when `--grid` is given, because it also selects
case-specific setup. For the flat plate (`stype = 4`), bottom faces with
i < `iplate` (= 11, in `src/util.f90`) are slip walls, so the plate starts at
x = 0.1 on `plate.dat`.

With `cfl_ramp = 1` the CFL is set at these iterations (meant for implicit
runs; edit `cfl_ramp` in `src/util.f90` to change it):

| Iteration | 1 | 20 | 30 | 40 | 500 | 700 | 1000 | 1250 | 1500 |
|---|---|---|---|---|---|---|---|---|---|
| CFL | 0.4 | 0.8 | 1.6 | 3.2 | 6.4 | 12.8 | 25.6 | 50 | 100 |

## Grid files

Tecplot ASCII, node coordinates only, in block format:

```
<title line>
ZONE I= <nx> J= <ny> ...
<nx*ny x-coordinates, i fastest>
<nx*ny y-coordinates, i fastest>
```

## Output

| File | Contents |
|---|---|
| `data.dat` | Tecplot ASCII: X, Y at nodes; rho, u, v, p and Mach cell-centred |
| `residual.dat` | Iteration and RMS change of the four conserved variables, each normalised by its freestream value |
| `checkdf.dat` | Grid-metric debug output |

Each iteration also prints `iteration, physical time, CFL, dt, residual` to
the screen.

To view a solution in ParaView: open `data.dat` with the Tecplot reader, click
Apply, switch the view to 2D and colour by `m`, `p` or `rho`.

## Notes

- The cylinder case reproduces the bow-shock stand-off from Billig's
  correlation (about 0.44 R at Mach 6).
- The implicit Jacobians use an unblended Steger–Warming-type split. The
  blended MSW Jacobians are inconsistent with the flux at shocks and made the
  implicit scheme unstable above CFL ~20.
- An MPI-parallel version is in progress.
