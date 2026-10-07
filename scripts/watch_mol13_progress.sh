#!/usr/bin/env bash
# Live monitor for mol13_pose1_monomer equil/production (run in a separate terminal).
# Usage: watch -n 30 bash /path/to/urease_md/scripts/watch_mol13_progress.sh
set -euo pipefail
SYS=mol13_pose1_monomer
BASE=/path/to/urease_md
CURRENT_OUT=""
CURRENT_SEG="(none)"

# Last matching file wins (later segments overwrite).
for N in 4 5; do
  OUT=${BASE}/simulations/${SYS}/equil_seg${N}/equil_seg${N}.out
  if [[ -f "$OUT" ]]; then
    CURRENT_OUT=$OUT
    CURRENT_SEG=seg${N}
  fi
done
for N in 1 2 3; do
  OUT=${BASE}/simulations/${SYS}/production_rep${N}/production_rep${N}.out
  if [[ -f "$OUT" ]]; then
    CURRENT_OUT=$OUT
    CURRENT_SEG=rep${N}
  fi
done

echo "=== mol13 === $(date) === Active: ${CURRENT_SEG} ==="
if [[ -z "${CURRENT_OUT}" ]]; then
  echo "No equil_seg4/5 or production_rep .out found yet."
  exit 0
fi

NSTEP=$(grep -E "^ NSTEP" "${CURRENT_OUT}" 2>/dev/null | tail -1 | awk '{print $3}')
NS=$(awk -v n="${NSTEP:-0}" 'BEGIN {printf "%.2f", n * 0.004 / 1000}')
echo "Progress: NSTEP=${NSTEP}  ${NS}/10 ns (equil) or scale for prod"
tail -15 "${CURRENT_OUT}" 2>/dev/null | grep -E "NSTEP|TEMP|Density|NMR" || true
NAN=$(grep -c "NaN" "${CURRENT_OUT}" 2>/dev/null || true)
if [[ "${NAN}" == "0" ]]; then
  echo "NaN: CLEAN"
else
  echo "NaN: ${NAN} FOUND"
fi
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi --query-gpu=index,utilization.gpu,temperature.gpu,power.draw \
    --format=csv,noheader 2>/dev/null | grep -v "^0," || nvidia-smi --query-gpu=index,utilization.gpu,temperature.gpu,power.draw --format=csv,noheader
fi
