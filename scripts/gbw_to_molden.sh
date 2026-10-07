#!/usr/bin/env bash
# scripts/gbw_to_molden.sh
#
# Convert ORCA resp.gbw → resp.molden for Multiwfn (orca_2mkl).
# ORCA 6 writes the molden file as resp.molden.input; we rename to resp.molden.
#
# Usage:
#   bash scripts/gbw_to_molden.sh [resp_results_dir]
# Env:
#   ORCA_2MKL  default: /opt/software/ORCA/build/orca_2mkl

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RESULTS_DIR="${1:-${ROOT}/resp_results}"
ORCA_2MKL="${ORCA_2MKL:-/opt/software/ORCA/build/orca_2mkl}"

log_prefix="[gbw_to_molden]"
log_msg() { echo "[$(date '+%F %T')] ${log_prefix} $*"; }

[[ -x "${ORCA_2MKL}" ]] || { echo "${log_prefix} ERROR: orca_2mkl not executable: ${ORCA_2MKL}" >&2; exit 2; }
[[ -d "${RESULTS_DIR}" ]] || { echo "${log_prefix} ERROR: not a directory: ${RESULTS_DIR}" >&2; exit 2; }

ok=0
skip=0
fail=0

for job_dir in "${RESULTS_DIR}"/mol*/; do
  [[ -d "${job_dir}" ]] || continue
  mol_id="$(basename "${job_dir}")"
  if [[ ! -f "${job_dir}/resp.gbw" ]]; then
    log_msg "${mol_id}: skip (no resp.gbw)"
    skip=$((skip + 1))
    continue
  fi
  (
    cd "${job_dir}" || exit 1
    "${ORCA_2MKL}" resp -molden >/dev/null 2>&1
    if [[ -f resp.molden.input ]]; then
      mv -f resp.molden.input resp.molden
    elif [[ -f resp.molden ]]; then
      :
    else
      echo "${mol_id}: no molden output" >&2
      exit 1
    fi
  ) && {
    log_msg "${mol_id}: wrote resp.molden"
    ok=$((ok + 1))
  } || {
    log_msg "${mol_id}: FAILED"
    fail=$((fail + 1))
  }
done

log_msg "done: ok=${ok} skip=${skip} fail=${fail}"
[[ "${fail}" -eq 0 ]] || exit 1
