#!/usr/bin/env bash
# Usage: bash scripts/fix_trajectory.sh <system_name>
# e.g.:  bash scripts/fix_trajectory.sh mol1_pose2_monomer
# Fixes periodic boundary artifacts and aligns trajectory to reference structure
# (autoimage + backbone RMS fit) for all production segment trajectories per replicate.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${SCRIPT_DIR}/../config.env"

SYS="${1:?Usage: bash scripts/fix_trajectory.sh <system_name>}"

PRMTOP="${PROJECT_ROOT}/systems/${SYS}/solvated_hmr.prmtop"
SIM="${PROJECT_ROOT}/simulations/${SYS}"

[[ -f "${PRMTOP}" ]] || {
  echo "[FAIL] Missing topology: ${PRMTOP}" >&2
  exit 1
}
[[ -d "${SIM}" ]] || {
  echo "[FAIL] Missing simulations dir: ${SIM}" >&2
  exit 1
}

CPPTRAJ_CMD="${CPPTRAJ:-}"
if [[ -z "${CPPTRAJ_CMD}" ]]; then
  CPPTRAJ_CMD="$(command -v cpptraj 2>/dev/null || true)"
fi
if [[ -z "${CPPTRAJ_CMD}" && -n "${AMBER_HOME:-}" ]]; then
  CPPTRAJ_CMD="${AMBER_HOME}/bin/cpptraj"
fi
[[ -x "${CPPTRAJ_CMD}" ]] || {
  echo "[FAIL] cpptraj not found (set CPPTRAJ or AMBER_HOME)." >&2
  exit 1
}

# Last protein/solute residue index before solvent (for backbone RMS mask :1-NRES&!@H=)
NRES="$(python3 "${SCRIPT_DIR}/_prmtop_last_solute_residue.py" "${PRMTOP}")" || {
  echo "[FAIL] Could not determine NRES from ${PRMTOP}" >&2
  exit 1
}

echo "[INFO] System=${SYS}  prmtop=${PRMTOP}  NRES=${NRES} (mask :1-${NRES}&!@H=)"
echo "[INFO] cpptraj=${CPPTRAJ_CMD}"

shopt -s nullglob

# Discover replicate indices from prod_rep##_seg01 directories
declare -A SEEN_REP=()
for segdir in "${SIM}"/prod_rep*_seg01/; do
  [[ -d "${segdir}" ]] || continue
  base="$(basename "${segdir}")" # prod_rep01_seg01
  if [[ "${base}" =~ ^prod_rep([0-9]+)_seg01$ ]]; then
    rep="${BASH_REMATCH[1]}"
    SEEN_REP["${rep}"]=1
  fi
done

if [[ ${#SEEN_REP[@]} -eq 0 ]]; then
  echo "[WARN] No prod_rep##_seg01 directories under ${SIM} — nothing to process."
  exit 0
fi

mapfile -t REP_LIST < <(printf '%s\n' "${!SEEN_REP[@]}" | sort -n)

for rep in "${REP_LIST[@]}"; do
  # Ordered segment NC files: prod_rep01_seg01.nc, prod_rep01_seg02.nc, ...
  mapfile -t NC_FILES < <(find "${SIM}" -maxdepth 2 -type f \
    -regex ".*/prod_rep${rep}_seg[0-9][0-9]*/prod_rep${rep}_seg[0-9][0-9]*\\.nc\$" | sort -V)
  if [[ ${#NC_FILES[@]} -eq 0 ]]; then
    echo "[WARN] Replicate ${rep}: no .nc files found — skip"
    continue
  fi

  out_dir="${SIM}/prod_rep${rep}"
  mkdir -p "${out_dir}"

  imaged="${out_dir}/traj_rep${rep}_imaged.nc"
  protein="${out_dir}/traj_rep${rep}_protein.nc"
  rmsd="${out_dir}/rmsd_rep${rep}.dat"
  clog="${out_dir}/fix_traj_rep${rep}.cpptraj.log"

  inp="$(mktemp)"
  {
    echo "parm ${PRMTOP}"
    for nc in "${NC_FILES[@]}"; do
      echo "trajin ${nc}"
    done
    echo "autoimage"
    echo "rms first :1-${NRES}"'&!@H= out '"${rmsd}"
    echo "trajout ${imaged} netcdf"
    echo "strip :WAT,Na+,Cl-"
    echo "trajout ${protein} netcdf"
    echo "run"
  } > "${inp}"

  echo "=== Replicate ${rep}: ${#NC_FILES[@]} segment file(s) → ${imaged} / ${protein} ==="
  "${CPPTRAJ_CMD}" -i "${inp}" 2>&1 | tee "${clog}"
  rm -f "${inp}"

  # Frame counts from log (cpptraj version-dependent strings)
  in_frames="$(grep -E 'contains [0-9]+ frames|Read [0-9]+ frames' "${clog}" | tail -1 | grep -oE '[0-9]+' | tail -1 || true)"
  out_frames="$(grep -E 'Written.*frames|Wrote [0-9]+ frames' "${clog}" | tail -1 | grep -oE '[0-9]+' | tail -1 || true)"
  if [[ -z "${in_frames}" ]]; then
    in_frames="(see ${clog})"
  fi
  if [[ -z "${out_frames}" ]]; then
    out_frames="${in_frames}"
  fi
  echo "[SUMMARY] rep=${rep}  segments=${#NC_FILES[@]}  input_frames≈${in_frames}  output_imaged+protein_frames≈${out_frames}"
done

echo "[OK] fix_trajectory.sh finished for ${SYS}"
