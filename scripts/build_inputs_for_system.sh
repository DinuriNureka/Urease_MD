#!/usr/bin/env bash
# scripts/build_inputs_for_system.sh <SYS>
# Build per-system restraints/ni_ni.RST and input_files/{equil_seg,production}.in
set -euo pipefail
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/post_equil_lib.sh
source "${_SCRIPT_DIR}/post_equil_lib.sh"
source "${_SCRIPT_DIR}/../config.env"

SYS=${1:?"Usage: bash scripts/build_inputs_for_system.sh <system_name>"}

PRMTOP="${PROJECT_ROOT}/systems/${SYS}/solvated_hmr.prmtop"
SIM="${PROJECT_ROOT}/simulations/${SYS}"
EQUIL3_RST="${SIM}/equil3/equil3.rst7"
RST_OUT="${SIM}/restraints/ni_ni.RST"
MDIN_EQ="${SIM}/input_files/equil_seg.in"
MDIN_PR="${SIM}/input_files/production.in"
MASK_TMP="${SIM}/restraints/.ni_mask_cpptraj.txt"

[[ -f "${PRMTOP}" ]] || {
  echo "[FAIL] ${SYS}: missing ${PRMTOP}"
  exit 1
}
[[ -f "${EQUIL3_RST}" ]] || {
  echo "[FAIL] ${SYS}: equil3 not complete — run run_equilibration.sh first (${EQUIL3_RST} missing)"
  exit 1
}

if ! post_equil_extract_ni_atoms "${PRMTOP}" "${EQUIL3_RST}" "${MASK_TMP}"; then
  echo "[FAIL] ${SYS}: could not extract two NI atom indices"
  exit 1
fi

post_equil_write_ni_ni_rst "${RST_OUT}" "${NI1_ATOM}" "${NI2_ATOM}"
DISANG_ABS="$(post_equil_disang_short_abs "${SYS}" "${RST_OUT}")"
post_equil_write_mdin_pair "${DISANG_ABS}" "${MDIN_EQ}" "${MDIN_PR}"
post_equil_ensure_short_links "${SYS}" || exit 1

# Verification
if ! grep -q 'mcbarint=100' "${MDIN_EQ}"; then
  echo "[FAIL] ${SYS}: mcbarint=100 missing in equil_seg.in"
  exit 1
fi
if ! grep -q 'ig=-1' "${MDIN_EQ}" || ! grep -q 'ig=-1' "${MDIN_PR}"; then
  echo "[FAIL] ${SYS}: ig=-1 missing in one or both mdin files"
  exit 1
fi
DISANG_EXPECT="${PROJECT_ROOT}/nmr_disang/${SYS}.RST"
if ! grep -qF "DISANG=${DISANG_EXPECT}" "${MDIN_EQ}"; then
  echo "[FAIL] ${SYS}: equil_seg.in DISANG= must use short path ${DISANG_EXPECT}"
  exit 1
fi
[[ -L "${DISANG_EXPECT}" || -f "${DISANG_EXPECT}" ]] || {
  echo "[FAIL] ${SYS}: DISANG symlink missing: ${DISANG_EXPECT}"
  exit 1
}
if ! grep -qE 'iat=[0-9]+,[0-9]+,' "${RST_OUT}"; then
  echo "[FAIL] ${SYS}: iat= line malformed in ni_ni.RST"
  exit 1
fi

echo "[PASS] ${SYS} inputs built (NI iat=${NI1_ATOM},${NI2_ATOM} DISANG=${DISANG_EXPECT})"
