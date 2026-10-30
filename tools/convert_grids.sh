#!/usr/bin/env bash
# Convert the structured example grids in grids/ to SU2 meshes, with marker
# names that match each case's boundaries.
#   tools/convert_grids.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONV="python3 $ROOT/tools/structured2su2.py"
G="$ROOT/grids"

# Channel with a ramp (and the shock-aligned version): inflow left, outflow
# right, walls top and bottom
$CONV "$G/ramp.dat"  "$G/ramp.su2"  --imin inlet --imax outlet --jmin wall --jmax wall
$CONV "$G/rampb.dat" "$G/rampb.su2" --imin inlet --imax outlet --jmin wall --jmax wall

# Cylinder (front half): i runs around the body, j outward. Both i ends are
# outflow, the outer boundary is inflow, j=1 is the cylinder wall
$CONV "$G/cyl.dat" "$G/cyl.su2" --imin outlet --imax outlet --jmin wall --jmax farfield
# Same mesh split into triangles (for testing on triangle meshes)
$CONV "$G/cyl.dat" "$G/cyl_tri.su2" --imin outlet --imax outlet --jmin wall --jmax farfield --triangulate

# Flat plate: the bottom is a symmetry plane up to node 10 (x = 0.09), then
# the plate; the top is a mirror boundary
$CONV "$G/plate.dat" "$G/plate.su2" --imin inlet --imax outlet --jmax top \
      --split jmin:10:symmetry:plate

# Closed box for the shock tube: mirror ends, walls top and bottom
$CONV "$G/block.dat" "$G/block.su2" --imin ends --imax ends --jmin wall --jmax wall
