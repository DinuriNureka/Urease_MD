#!/usr/bin/env bash
# Apply Multiwfn charges to gaff2_params/molNN/molNN.mol2 using
# resp_results/molNN/resp_multiwfn.chg
#
# Usage:
#   bash scripts/apply_multiwfn_charges_batch.sh
#
# Env:
#   GAFF2_DIR  default: gaff2_params

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
GAFF2_DIR="${GAFF2_DIR:-${ROOT}/gaff2_params}"
RESULTS_DIR="${ROOT}/resp_results"
PY="${ROOT}/scripts/apply_multiwfn_charges.py"

ok=0
fail=0

for job_dir in "${RESULTS_DIR}"/mol*/; do
  [[ -d "${job_dir}" ]] || continue
  mol_id="$(basename "${job_dir}")"
  chg="${job_dir}/resp_multiwfn.chg"
  mol_in="${GAFF2_DIR}/${mol_id}/${mol_id}.mol2"
  mol_out="${GAFF2_DIR}/${mol_id}/${mol_id}.mol2"
  if [[ ! -f "${chg}" ]]; then
    echo "SKIP ${mol_id}: no ${chg}" >&2
    fail=$((fail + 1))
    continue
  fi
  if [[ ! -f "${mol_in}" ]]; then
    echo "SKIP ${mol_id}: no ${mol_in}" >&2
    fail=$((fail + 1))
    continue
  fi
  tmp="${mol_out}.tmp"
  if python3 "${PY}" --mol2 "${mol_in}" --chg "${chg}" --out "${tmp}"; then
    mv -f "${tmp}" "${mol_out}"
    ok=$((ok + 1))
  else
    rm -f "${tmp}"
    fail=$((fail + 1))
  fi
done

echo "apply_multiwfn_charges_batch: ok=${ok} fail=${fail}"
[[ "${fail}" -eq 0 ]] || exit 1
