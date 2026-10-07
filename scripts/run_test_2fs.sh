#!/usr/bin/env bash
# Usage: bash scripts/run_test_2fs.sh <system_name>
# Runs a single 10 ns production test at 2 fs without HMR (no segments, no replicates).
# Uses the STANDARD (non-HMR) topology: systems/<name>/solvated.prmtop
# Output: simulations/<name>/test_2fs/test_2fs.{out,rst7,nc,mdinfo}
#
# Prerequisites: solvated.prmtop and a completed equil_seg*.rst7 under simulations/<name>/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${SCRIPT_DIR}/../config.env"

SYS="${1:?Usage: bash scripts/run_test_2fs.sh <system_name>}"

PRMTOP="${PROJECT_ROOT}/systems/${SYS}/solvated.prmtop"
MDIN="${PROJECT_ROOT}/simulations/input_files/production_2fs_nohmr.in"
OUTDIR="${PROJECT_ROOT}/simulations/${SYS}/test_2fs"
DISANG_TARGET="${PROJECT_ROOT}/ni_ni.RST"

export CUDA_DEVICE_ORDER="${CUDA_DEVICE_ORDER:-PCI_BUS_ID}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-1}"

PMEMD="${PMEMD_CUDA_BIN:-/opt/software/AMBER/pmemd24_cuda128/bin/pmemd.cuda}"
[[ -x "${PMEMD}" ]] || {
  echo "[FAIL] Not executable: ${PMEMD}" >&2
  exit 1
}

[[ -f "${PRMTOP}" ]] || {
  echo "[FAIL] Missing non-HMR topology: ${PRMTOP}" >&2
  echo "       Generate with tLeap (solvated.prmtop) before running this test." >&2
  exit 1
}
[[ -f "${MDIN}" ]] || {
  echo "[FAIL] Missing ${MDIN}" >&2
  exit 1
}
[[ -f "${DISANG_TARGET}" ]] || {
  echo "[FAIL] Missing Ni–Ni restraint file: ${DISANG_TARGET}" >&2
  exit 1
}

mkdir -p "${OUTDIR}"
# Short symlink for documentation / tools that expect restraints next to run (optional)
ln -sfn "${DISANG_TARGET}" "${OUTDIR}/ni_ni.RST"

# Pick last equil_segNN.rst7 by largest numeric NN (equil_seg01 … equil_seg10; skips equil_seg3c etc.)
LAST_RST=""
MAX_N=-1
shopt -s nullglob
for f in "${PROJECT_ROOT}/simulations/${SYS}"/equil_seg*/equil_seg*.rst7; do
  [[ -f "${f}" ]] || continue
  dir="$(basename "$(dirname "${f}")")"
  num="${dir#equil_seg}"
  if [[ "${num}" =~ ^[0-9]+$ ]]; then
    if (( 10#num > MAX_N )); then
      MAX_N=$((10#num))
      LAST_RST="${f}"
    fi
  fi
done
shopt -u nullglob

[[ -n "${LAST_RST}" ]] || {
  echo "[FAIL] No equil_seg*/equil_seg*.rst7 found under simulations/${SYS}/" >&2
  exit 1
}

echo "=== 2 fs NO-HMR production test ==="
echo "system:       ${SYS}"
echo "topology:     ${PRMTOP}"
echo "mdin:         ${MDIN}"
echo "start rst7:   ${LAST_RST} (equil_seg index ${MAX_N})"
echo "output dir:   ${OUTDIR}"
echo "PMEMD:        ${PMEMD}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"
echo

OUT_BASE="${OUTDIR}/test_2fs"
"${PMEMD}" -O \
  -i "${MDIN}" \
  -o "${OUT_BASE}.out" \
  -p "${PRMTOP}" \
  -c "${LAST_RST}" \
  -r "${OUT_BASE}.rst7" \
  -x "${OUT_BASE}.nc" \
  -ref "${LAST_RST}" \
  -inf "${OUT_BASE}.mdinfo"

echo
echo "--- Performance (ns/day) ---"
if [[ -f "${OUT_BASE}.mdinfo" ]]; then
  grep 'ns/day' "${OUT_BASE}.mdinfo" | tail -5 || true
else
  grep 'ns/day' "${OUT_BASE}.out" | tail -5 || true
fi

echo "[OK] Outputs: ${OUT_BASE}.{out,rst7,nc,mdinfo}"
