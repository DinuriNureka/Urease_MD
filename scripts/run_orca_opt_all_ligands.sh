#!/usr/bin/env bash
# Run ORCA geometry optimizations for 36 ligands with bounded parallelism.
#
# Uses pre-generated setups in:
#   params/inhibitors/orca_opt/molNN/
# where each directory contains:
#   - molNN.xyz
#   - molNN.inp
#
# Parallel policy:
#   - launch up to 8 concurrent ORCA jobs
#   - each job uses nprocs 6 and maxcore 3000 (from input files)
#
# Outputs per job:
#   - params/inhibitors/orca_opt/molNN/molNN.log
#   - on success, copy gbw to resp_results/molNN/resp.gbw
#
# Usage (from project root):
#   bash scripts/run_orca_opt_all_ligands.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/config.env"

ORCA_BIN="${ORCA_PATH:-$(command -v orca || true)}"
MAP_CSV="${ROOT}/params/inhibitors/orca_opt/ligand_name_mapping.csv"
OPT_ROOT="${ROOT}/params/inhibitors/orca_opt"
RESULTS_ROOT="${ROOT}/resp_results"
MAX_JOBS="${ORCA_OPT_MAX_JOBS:-8}"

log_prefix="[orca_opt_all]"
log_msg() { echo "[$(date '+%F %T')] ${log_prefix} $*"; }

[[ -n "${ORCA_BIN}" && -x "${ORCA_BIN}" ]] || {
  echo "${log_prefix} ERROR: ORCA executable not found. ORCA_PATH='${ORCA_PATH:-}'" >&2
  exit 2
}
[[ -f "${MAP_CSV}" ]] || {
  echo "${log_prefix} ERROR: mapping CSV missing: ${MAP_CSV}" >&2
  exit 2
}

mkdir -p "${RESULTS_ROOT}"

FAILFILE="${OPT_ROOT}/._orca_fail_ids.$$"
OKFILE="${OPT_ROOT}/._orca_ok_ids.$$"
rm -f "${FAILFILE}" "${OKFILE}"

run_one() {
  local mol_id="$1"
  local workdir="${OPT_ROOT}/${mol_id}"
  local inp="${workdir}/${mol_id}.inp"
  local log="${workdir}/${mol_id}.log"
  local gbw="${workdir}/${mol_id}.gbw"

  if [[ ! -f "${inp}" ]]; then
    log_msg "${mol_id}: missing input ${inp}"
    echo "${mol_id}" >> "${FAILFILE}"
    return
  fi

  log_msg "${mol_id}: started"
  (
    cd "${workdir}" || exit 1
    "${ORCA_BIN}" "${mol_id}.inp" > "${mol_id}.log" 2>&1
  )
  local rc=$?

  if [[ ${rc} -ne 0 ]]; then
    log_msg "${mol_id}: FAIL exit=${rc}"
    echo "${mol_id}" >> "${FAILFILE}"
    return
  fi

  if ! grep -q 'TERMINATED NORMALLY' "${log}"; then
    log_msg "${mol_id}: FAIL no normal termination banner"
    echo "${mol_id}" >> "${FAILFILE}"
    return
  fi

  if [[ ! -f "${gbw}" ]]; then
    log_msg "${mol_id}: FAIL missing ${gbw}"
    echo "${mol_id}" >> "${FAILFILE}"
    return
  fi

  mkdir -p "${RESULTS_ROOT}/${mol_id}"
  cp -f "${gbw}" "${RESULTS_ROOT}/${mol_id}/resp.gbw"
  log_msg "${mol_id}: OK"
  echo "${mol_id}" >> "${OKFILE}"
}

# Read mol_id list from mapping CSV (skip header)
mapfile -t MOL_IDS < <(tail -n +2 "${MAP_CSV}" | awk -F',' '{print $1}')

log_msg "ORCA binary: ${ORCA_BIN}"
log_msg "Launching ${#MOL_IDS[@]} jobs with max parallel=${MAX_JOBS}"

for mol_id in "${MOL_IDS[@]}"; do
  while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
    wait -n || true
  done
  run_one "${mol_id}" &
done

wait || true

ok=0
fail=0
[[ -f "${OKFILE}" ]] && ok=$(sort -u "${OKFILE}" | wc -l | tr -d ' ')
[[ -f "${FAILFILE}" ]] && fail=$(sort -u "${FAILFILE}" | wc -l | tr -d ' ')

rm -f "${OKFILE}" "${FAILFILE}"
log_msg "finished: succeeded=${ok} failed=${fail} total=${#MOL_IDS[@]}"

[[ "${fail}" -eq 0 ]]
