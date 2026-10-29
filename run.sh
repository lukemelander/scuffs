#!/usr/bin/env bash
# Run an example case:
#   ./run.sh examples/<case> [extra solver options]
# Builds the solver if needed, then runs it inside the case directory so the
# outputs (data.dat, residual.dat, run.log) are written next to its input.inp.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
BIN="$ROOT/bin/scuffs"

if [ $# -lt 1 ]; then
  echo "usage: $0 examples/<case> [solver options]"
  echo "cases:"; ls -1 "$ROOT/examples"
  exit 1
fi

CASE="$(cd "$1" && pwd)"; shift
[ -f "$CASE/input.inp" ] || { echo "no input.inp in $CASE"; exit 1; }
[ -f "$CASE/case.conf" ] || { echo "no case.conf in $CASE"; exit 1; }

# case.conf sets GRID (a file name in grids/)
GRID=""
source "$CASE/case.conf"
[ -n "$GRID" ] || { echo "GRID not set in $CASE/case.conf"; exit 1; }

[ -x "$BIN" ] || make -C "$ROOT"

cd "$CASE"
"$BIN" --input input.inp --grid "$ROOT/grids/$GRID" "$@" | tee run.log
