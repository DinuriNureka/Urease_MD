#!/usr/bin/env bash
# scripts/run_system_full.sh <SYS>
# Full post-equil3 pipeline: equil_seg1–5 (20 ns each) + production_rep1–3 (20 ns each)
set -euo pipefail
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/post_equil_lib.sh
source "${_SCRIPT_DIR}/post_equil_lib.sh"
source "${_SCRIPT_DIR}/../config.env"

export PATH="/usr/bin:${AMBER_HOME}/bin:/usr/local/bin:${PATH}"
export CUDA_DEVICE_ORDER="${CUDA_DEVICE_ORDER:-PCI_BUS_ID}"
export CUDA_VISIBLE_DEVICES=1

PMEMD="${PMEMD_CUDA_BIN:-/opt/software/AMBER/pmemd24_cuda128/bin/pmemd.cuda}"
[[ -x "${PMEMD}" ]] || {
  echo "[FAIL] Not executable: ${PMEMD}"
  exit 1
}

SYS=${1:?"Usage: bash scripts/run_system_full.sh <system_name>"}

LOG="${PROJECT_ROOT}/logs/run_${SYS}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

echo "=== run_system_full.sh ${SYS} ==="
echo "LOG: ${LOG}"
echo "PMEMD: ${PMEMD}  CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"

bash "${_SCRIPT_DIR}/build_inputs_for_system.sh" "${SYS}"

SIM="${PROJECT_ROOT}/simulations/${SYS}"
# pmemd -i/-p/-c use short paths via post_equil_lib.sh (Fortran 80-char limit)

declare -A SEG_STAT SEG_NAN
declare -A REP_STAT REP_NAN

run_equil_seg() {
  local n=$1 crd_in=$2
  local dir="${SIM}/equil_seg${n}"
  local pre="equil_seg${n}"
  if [[ -s "${dir}/${pre}.rst7" ]]; then
    echo "seg${n}: already done, skipping"
    SEG_STAT[${n}]="SKIP"
    SEG_NAN[${n}]="—"
    return 0
  fi
  echo "--- seg${n}: starting from ${crd_in} ---"
  set +e
  post_equil_run_pmemd_segment "${dir}" "${MDIN_EQ}" "${PRMTOP}" "${crd_in}" "${pre}" "${PMEMD}"
  local rc=$?
  set -e
  if [[ "${rc}" -ne 0 ]]; then
    echo "[FAIL] seg${n}: pmemd exit code ${rc}"
    exit 1
  fi
  if ! post_equil_verify_mdout "${dir}/${pre}.out" 2500000; then
    echo "[FAIL] seg${n}: mdout verification failed"
    exit 1
  fi
  SEG_STAT[${n}]="DONE"
  SEG_NAN[${n}]="0"
}

run_prod_rep() {
  local n=$1 crd_in=$2
  local dir="${SIM}/production_rep${n}"
  local pre="production_rep${n}"
  if [[ -s "${dir}/${pre}.rst7" ]]; then
    echo "rep${n}: already done, skipping"
    REP_STAT[${n}]="SKIP"
    REP_NAN[${n}]="—"
    return 0
  fi
  echo "--- rep${n}: starting from ${crd_in} ---"
  set +e
  post_equil_run_pmemd_segment "${dir}" "${SYS}" prod "${crd_in}" "${pre}" "${PMEMD}"
  local rc=$?
  set -e
  if [[ "${rc}" -ne 0 ]]; then
    echo "[FAIL] rep${n}: pmemd exit code ${rc}"
    exit 1
  fi
  if ! post_equil_verify_mdout "${dir}/${pre}.out" 5000000; then
    echo "[FAIL] rep${n}: mdout verification failed"
    exit 1
  fi
  REP_STAT[${n}]="DONE"
  REP_NAN[${n}]="0"
}

# seg1 ← equil3
run_equil_seg 1 "${SIM}/equil3/equil3.rst7"
# seg2–5 chain
run_equil_seg 2 "${SIM}/equil_seg1/equil_seg1.rst7"
run_equil_seg 3 "${SIM}/equil_seg2/equil_seg2.rst7"
run_equil_seg 4 "${SIM}/equil_seg3/equil_seg3.rst7"
run_equil_seg 5 "${SIM}/equil_seg4/equil_seg4.rst7"

SEG5_RST="${SIM}/equil_seg5/equil_seg5.rst7"
[[ -s "${SEG5_RST}" ]] || {
  echo "[FAIL] Missing ${SEG5_RST} — cannot start production"
  exit 1
}

run_prod_rep 1 "${SEG5_RST}"
run_prod_rep 2 "${SEG5_RST}"
run_prod_rep 3 "${SEG5_RST}"

fmt_seg() {
  local n=$1
  local st="${SEG_STAT[${n}]:-?}"
  if [[ "${st}" == "DONE" ]]; then
    echo "seg${n}:  DONE  (NaN=${SEG_NAN[${n}]}, NSTEP=2500000)"
  elif [[ "${st}" == "SKIP" ]]; then
    echo "seg${n}:  SKIP  (existing rst7)"
  else
    echo "seg${n}:  ${st}"
  fi
}

fmt_rep() {
  local n=$1
  local st="${REP_STAT[${n}]:-?}"
  if [[ "${st}" == "DONE" ]]; then
    echo "rep${n}:  DONE  (NaN=${REP_NAN[${n}]}, NSTEP=5000000)"
  elif [[ "${st}" == "SKIP" ]]; then
    echo "rep${n}:  SKIP  (existing rst7)"
  else
    echo "rep${n}:  ${st}"
  fi
}

echo ""
echo "=== ${SYS} COMPLETE ==="
fmt_seg 1
fmt_seg 2
fmt_seg 3
fmt_seg 4
fmt_seg 5
fmt_rep 1
fmt_rep 2
fmt_rep 3
echo "Total production trajectory: 60 ns (3 × 20 ns)"
