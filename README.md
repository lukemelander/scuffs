# SCUFFS

**S**hock-capturing **C**ompressible **U**nstructured **F**inite-volume **F**low **S**olver

SCUFFS is a first-order, cell-centred finite-volume solver for the 2D Euler
equations on unstructured meshes of triangles and quads, written in Fortran
90. It was written for a graduate course on finite-volume methods for
hypersonic CFD on structured grids, and has since been cleaned up, debugged
and moved to unstructured meshes.

- Inviscid fluxes: modified Steger–Warming (MSW) flux-vector splitting with
  pressure-gradient blending
- Time integration: explicit (forward Euler), or implicit (backward Euler with
  point-Jacobi relaxation) with an optional CFL ramp for steady-state runs
- Boundary conditions by named mesh marker: supersonic inflow,
  characteristic outflow, slip wall, mirror, adiabatic wall
- Meshes: SU2 format; a converter turns structured grids into SU2 meshes
- Output: legacy VTK files that open in ParaView

Viscous terms (full Navier–Stokes) and second-order fluxes are in progress;
the structured v1.0 release (git tag `v1.0`) has thin-layer viscous terms.

## Repository layout

```
.
├── Makefile        build rules (gfortran by default)
├── run.sh          runs an example case
├── src/            Fortran sources (mesh, util, solver, main)
├── grids/          SU2 meshes (*.su2) and the original structured grids (*.dat)
├── tools/          grid converter and comparison scripts (Python 3)
└── examples/       one directory per case: input.inp + case.conf
```

## Building

Requires a Fortran compiler (gfortran by default) and `make`.

```
make            # optimised build -> bin/scuffs
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
writes `data.vtk`, `residual.dat` and `run.log` there. Extra arguments are
passed on to the solver.

The solver can also be run directly:

```
bin/scuffs --grid path/to/mesh.su2 --input path/to/input.inp
bin/scuffs --grid path/to/mesh.su2 --mesh-check
```

| Option | Meaning |
|---|---|
| `-g`, `--grid FILE` | SU2 mesh (required) |
| `-i`, `--input FILE` | Input file (default `input.inp`) |
| `--mesh-check` | Read and check the mesh, write `mesh.vtk` and `mesh_boundary.vtk`, then stop |
| `-h`, `--help` | Show usage |

`--input=FILE` and `--grid=FILE` also work.

## Examples

| Case | Mesh | Description |
|---|---|---|
| `ramp_explicit` | `ramp.su2` | Mach 2.5 inviscid flow in a channel with a ramp, explicit |
| `ramp_implicit` | `ramp.su2` | Same, implicit with CFL ramp |
| `ramp_aligned_implicit` | `rampb.su2` | Same, on a grid aligned with the ramp shock |
| `cylinder_m6_implicit` | `cyl.su2` | Mach 6 inviscid flow over a cylinder (bow shock), implicit |
| `shock_tube` | `block.su2` | Shock tube in a closed box (gas at rest, lower half at 2x pressure and density) |
| `plate_viscous_implicit` | `plate.su2` | Mach 1.5 laminar flat plate, implicit (needs the viscous terms; not runnable yet) |
| `plate_viscous_explicit` | `plate.su2` | Same, explicit (not runnable yet) |

Each case directory has `input.inp` (solver settings) and `case.conf`, which
names the mesh in `grids/` that `run.sh` passes with `--grid`.

## Input file

The input file is read with list-directed reads, so anything after the values
on a line is treated as a comment.

```
0 1 1              !ivisc, impl, cfl_ramp
1000.0 0.0071 6.0  !p_inf rho_inf M_inf
2001 500           !iterations, save interval
0.4d0              !CFL
3                  !number of boundary markers
farfield 1         !marker name, BC type
outlet   2
wall     3
1                  !initial condition
```

| Line | Values | Meaning |
|---|---|---|
| 1 | `ivisc` | 0 = inviscid (Euler). 1 = viscous (not available yet on unstructured meshes) |
|   | `impl` | 0 = explicit, 1 = implicit (point Jacobi, `kmax`+1 = 10 sweeps, set in `src/util.f90`) |
|   | `cfl_ramp` | 0 = constant CFL, 1 = built-in CFL ramp (see below) |
| 2 | `p_inf rho_inf M_inf` | Freestream pressure [Pa], density [kg/m^3] and Mach number. The flow is in +x |
| 3 | `iter` | Number of iterations |
|   | `isave` | Write `data.vtk` every `isave` iterations (it is also written at the end of the run) |
| 4 | `CFL` | CFL number (starting value if the ramp is on) |
| 5 | `nbc` | Number of boundary-condition lines that follow |
| 6.. | `name type` | One line per mesh marker: the marker name from the `.su2` file and its BC type. Every marker in the mesh needs a line |
| last | `switch` | Initial condition: 1 = uniform freestream; 2 = shock tube (gas at rest, cells below the middle of the mesh's y range at 2x p and rho); 3 = freestream with a 1% bump in the cells at the middle of the x range |

Boundary types:

| Code | Type |
|---|---|
| 1 | Supersonic inflow (ghost cell held at freestream) |
| 2 | Outflow (characteristic; only outgoing waves are taken from the interior) |
| 3 | Inviscid (slip) wall / symmetry plane |
| 4 | Mirror (ghost cell = interior cell) |
| 5 | Adiabatic no-slip wall (a slip wall for the inviscid flux) |

With `cfl_ramp = 1` the CFL is set at these iterations (meant for implicit
runs; edit `cfl_ramp` in `src/util.f90` to change it):

| Iteration | 1 | 20 | 30 | 40 | 500 | 700 | 1000 | 1250 | 1500 |
|---|---|---|---|---|---|---|---|---|---|
| CFL | 0.4 | 0.8 | 1.6 | 3.2 | 6.4 | 12.8 | 25.6 | 50 | 100 |

## Meshes

Meshes use the SU2 native format (2D): triangles (type 5) and quads (type 9),
with boundary edges grouped into named markers. Gmsh can export SU2 directly.

The structured grids from v1.0 are converted with:

```
tools/convert_grids.sh                          # all example grids
tools/structured2su2.py grids/plate.dat grids/plate.su2 \
    --imin inlet --imax outlet --jmax top --split jmin:10:symmetry:plate
```

`--imin/--imax/--jmin/--jmax` name the four sides (sides with the same name
become one marker); `--split SIDE:NODE:BEFORE:AFTER` splits a side into two
markers at a node index.

## Output

| File | Contents |
|---|---|
| `data.vtk` | Legacy VTK: cell-centred rho, u, v, p, Mach and a velocity vector |
| `residual.dat` | Iteration and RMS change of the four conserved variables, each normalised by its freestream value |

Each iteration also prints `iteration, physical time, CFL, dt, residual` to
the screen.

To view a solution in ParaView: open `data.vtk`, click Apply, switch the view
to 2D and colour by `mach`, `p` or `rho`.

`tools/compare_v1.py v1/data.dat data.vtk` compares a solution with the
structured v1.0 code on the same grid (cells line up one to one).

## Notes

- The cylinder case reproduces the bow-shock stand-off from Billig's
  correlation (about 0.44 R at Mach 6).
- The implicit Jacobians use an unblended Steger–Warming-type split. The
  blended MSW Jacobians are inconsistent with the flux at shocks and made the
  implicit scheme unstable above CFL ~20.
- An MPI-parallel version (pSCUFFS) is planned.
