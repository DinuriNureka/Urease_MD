#!/bin/bash
# scripts/run_mmpbsa.sh
# -------------------------------------------------------
# MM-GBSA binding free energy for one system.
#
# Usage:
#     bash scripts/run_mmpbsa.sh mol01_pose01
#
# Inputs:
#     systems/<name>/solvated_hmr.parm7
#     systems/<name>/prod_rep01/prod.nc (or simulations/<name>/prod/prod_*.nc legacy)
#
# Outputs:
#     analysis/mmpbsa/results/<name>/mmpbsa.dat
# -------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config.env"

NAME=${1:?"Usage: bash run_mmpbsa.sh <name>"}
SYS="${PROJECT_ROOT}/systems/${NAME}"
if [[ -f "${SYS}/solvated_hmr.prmtop" ]]; then
    PARM="${SYS}/solvated_hmr.prmtop"
elif [[ -f "${SYS}/solvated_hmr.parm7" ]]; then
    PARM="${SYS}/solvated_hmr.parm7"
else
    PARM="${SYS}/solvated.parm7"
fi
TRAJ="${PROJECT_ROOT}/simulations/${NAME}/prod"
OUT="${PROJECT_ROOT}/analysis/mmpbsa/results/${NAME}"

log() { echo "[$(date '+%H:%M:%S')] [$NAME] $1"; }

mkdir -p "${OUT}"

REIMAGED="${TRAJ}/${PROD_TRAJ_REIMAGED:-production_reimaged.nc}"
PRE="${PROJECT_ROOT}/scripts/preprocess_trajectory.sh"
if [[ ! -f "${REIMAGED}" ]]; then
  log "Pre-processing trajectory (autoimage) for MM-GBSA..."
  shopt -s nullglob
  if [[ -f "${TRAJ}/prod.nc" ]]; then
    _PROD_NC=("${TRAJ}/prod.nc")
  else
    _PROD_NC=("${TRAJ}"/prod_*.nc)
  fi
  shopt -u nullglob
  [[ ${#_PROD_NC[@]} -gt 0 ]] || {
    log "[FAIL] No trajectory under ${TRAJ} — run production or preprocess manually."
    exit 1
  }
  bash "${PRE}" "${PARM}" "${REIMAGED}" "${_PROD_NC[@]}"
else
  log "Using existing re-imaged trajectory: ${REIMAGED}"
fi

log "Generating stripped topologies..."
ante-MMPBSA.py -p "${PARM}" -c "${OUT}/complex.parm7" -s ':WAT:Na+:Cl-'
ante-MMPBSA.py -p "${PARM}" -c "${OUT}/receptor.parm7" -s ':WAT:Na+:Cl-:LIG'
ante-MMPBSA.py -p "${PARM}" -c "${OUT}/ligand.parm7" -s ':WAT:Na+:Cl-' -l ':LIG'

log "Running MMPBSA.py..."
MMPBSA.py -O \
    -i "${PROJECT_ROOT}/analysis/mmpbsa/mmpbsa.in" \
    -o "${OUT}/mmpbsa.dat" \
    -sp "${PARM}" \
    -cp "${OUT}/complex.parm7" \
    -rp "${OUT}/receptor.parm7" \
    -lp "${OUT}/ligand.parm7" \
    -y "${REIMAGED}"

log "[DONE] Results: ${OUT}/mmpbsa.dat"
grep "DELTA TOTAL" "${OUT}/mmpbsa.dat" || true
