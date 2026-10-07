#!/usr/bin/env bash
# Re-image periodic trajectories with cpptraj autoimage (mandatory before analysis).
# Usage:
#   bash scripts/preprocess_trajectory.sh <prmtop> <out_traj.nc> <raw1.nc> [raw2.nc ...]
set -euo pipefail
_SCRIPT_ROOT="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=/dev/null
source "${_SCRIPT_ROOT}/../config.env" 2>/dev/null || true

PRMTOP="${1:?prmtop}"
OUT_TRAJ="${2:?output trajectory}"
shift 2
RAW_TRAJS=("$@")
[[ ${#RAW_TRAJS[@]} -ge 1 ]] || {
  echo "Usage: bash scripts/preprocess_trajectory.sh <prmtop> <out_traj.nc> <raw.nc> [more.nc ...]" >&2
  exit 2
}

CPPTRAJ_CMD="${CPPTRAJ:-}"
if [[ -z "${CPPTRAJ_CMD}" ]]; then
  CPPTRAJ_CMD="$(command -v cpptraj 2>/dev/null || true)"
fi
if [[ -z "${CPPTRAJ_CMD}" && -n "${AMBER_HOME:-}" ]]; then
  CPPTRAJ_CMD="${AMBER_HOME}/bin/cpptraj"
fi
[[ -x "${CPPTRAJ_CMD}" ]] || {
  echo "[ERROR] preprocess_trajectory.sh: cpptraj not found (set CPPTRAJ or AMBER_HOME)." >&2
  exit 1
}

INP="$(mktemp)"
trap 'rm -f "${INP}"' EXIT

{
  echo "parm ${PRMTOP}"
  for t in "${RAW_TRAJS[@]}"; do
    echo "trajin ${t}"
  done
  echo "autoimage"
  echo "trajout ${OUT_TRAJ} netcdf"
  echo "run"
} > "${INP}"

"${CPPTRAJ_CMD}" -i "${INP}"

if [[ ! -f "${OUT_TRAJ}" ]]; then
  echo "[ERROR] preprocess_trajectory.sh: output not created: ${OUT_TRAJ}" >&2
  exit 1
fi
echo "[OK] Re-imaged trajectory written to ${OUT_TRAJ}"
