#!/usr/bin/env bash
set -uo pipefail

: <<'DOC'
run_antechamber_resp.sh

Purpose:
  Build GAFF2 mol2/frcmod for ligands using ORCA 6 output. AmberTools accepts
  ORCA output as -fi orcout (not orcaout). Antechamber -c resp requires Gaussian
  gout/gesp; this script therefore:
    1) antechamber -fi orcout -c mul  (GAFF2 types + coordinates from ORCA)
    2) apply_pc_chelpg_to_mol2.py     (replace charges with ORCA CHELPG from resp.pc_chelpg)
    3) parmchk2

Inputs:
  - resp_results/rmsd_check.csv
  - resp_inputs/job_manifest.csv
  - resp_results/molNN/resp.out and resp.pc_chelpg

Outputs:
  - resp_results/molNN/molNN.mol2
  - resp_results/molNN/molNN.frcmod
  - gaff2_params/molNN/molNN.mol2
  - gaff2_params/molNN/molNN.frcmod
  - gaff2_params/antechamber.log

Example:
  bash run_antechamber_resp.sh
DOC

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR=resp_results
INPUTS_DIR=resp_inputs
PARAMS_DIR=gaff2_params
RMSD_CSV=resp_results/rmsd_check.csv
MANIFEST=resp_inputs/job_manifest.csv
LOG=${PARAMS_DIR}/antechamber.log

log_msg() {
  local log_prefix="[run_antechamber_resp]"
  local msg="$1"
  echo "[$(date '+%F %T')] ${msg}" | tee -a "${LOG}"
}

main() {
  local log_prefix="[main]"
  mkdir -p "${PARAMS_DIR}"
  : > "${LOG}"

  [[ -f "${MANIFEST}" ]] || { echo "${log_prefix} ERROR: missing ${MANIFEST}" >&2; exit 2; }
  [[ -f "${RMSD_CSV}" ]] || { echo "${log_prefix} ERROR: missing ${RMSD_CSV}" >&2; exit 2; }

  declare -A CHARGES
  while IFS=, read -r mol_id original_filename xyz_file net_charge multiplicity n_atoms validation_status; do
    [[ "${mol_id}" == "mol_id" ]] && continue
    CHARGES["${mol_id}"]="${net_charge}"
  done < "${MANIFEST}"

  local summary_file
  summary_file="$(mktemp)"
  local ok_count=0
  local total_count=0

  while IFS=, read -r mol_id original_filename n_atoms n_heavy_atoms rmsd_angstrom geometry_status scf_energy_Ha homo_eV lumo_eV gap_eV; do
    [[ "${mol_id}" == "mol_id" ]] && continue
    [[ "${geometry_status}" != "PASS" && "${geometry_status}" != "WARN" ]] && continue
    total_count=$((total_count + 1))

    local nc
    nc="${CHARGES[${mol_id}]:-0}"
    local mol2_ok="NO" frcmod_ok="NO" miss="NA" status="FAILED"
    local job_dir="${RESULTS_DIR}/${mol_id}"
    local dst_dir="${PARAMS_DIR}/${mol_id}"
    mkdir -p "${dst_dir}"

    if [[ ! -d "${job_dir}" ]]; then
      log_msg "${mol_id}: missing job directory ${job_dir}"
      echo "${mol_id},${nc},${mol2_ok},${frcmod_ok},${miss},${status}" >> "${summary_file}"
      continue
    fi

    (
      cd "${job_dir}" || exit 1
      antechamber \
        -i resp.out \
        -fi orcout \
        -o "${mol_id}_mul.mol2" \
        -fo mol2 \
        -c mul \
        -s 2 \
        -nc "${nc}" \
        -m 1 \
        -at gaff2 \
        -rn MOL \
        -pf y \
        >> "../../${LOG}" 2>&1
    ) || true

    if [[ -s "${job_dir}/${mol_id}_mul.mol2" && -f "${job_dir}/resp.pc_chelpg" ]]; then
      python3 "${SCRIPT_ROOT}/apply_pc_chelpg_to_mol2.py" \
        --mol2 "${job_dir}/${mol_id}_mul.mol2" \
        --pc "${job_dir}/resp.pc_chelpg" \
        --out "${job_dir}/${mol_id}.mol2" \
        >> "${LOG}" 2>&1 || true
      rm -f "${job_dir}/${mol_id}_mul.mol2"
    fi

    if [[ -s "${job_dir}/${mol_id}.mol2" ]]; then
      mol2_ok="YES"
    fi

    if [[ "${mol2_ok}" == "YES" ]]; then
      (
        cd "${job_dir}" || exit 1
        parmchk2 \
          -i "${mol_id}.mol2" \
          -f mol2 \
          -o "${mol_id}.frcmod" \
          -s gaff2 \
          >> "../../${LOG}" 2>&1
      ) || true
      if [[ -s "${job_dir}/${mol_id}.frcmod" ]]; then
        frcmod_ok="YES"
        miss="$(grep -c "ATTN, need" "${job_dir}/${mol_id}.frcmod" || true)"
      fi
    fi

    if [[ "${mol2_ok}" == "YES" && "${frcmod_ok}" == "YES" ]]; then
      cp -f "${job_dir}/${mol_id}.mol2" "${dst_dir}/"
      cp -f "${job_dir}/${mol_id}.frcmod" "${dst_dir}/"
      status="OK"
      ok_count=$((ok_count + 1))
      if [[ "${miss}" != "NA" && "${miss}" -gt 0 ]]; then
        log_msg "WARNING: ${mol_id}.frcmod has ${miss} missing parameters — inspect before tLeap"
      fi
    else
      status="FAILED"
    fi

    log_msg "${mol_id} | charge=${nc} | mol2=${mol2_ok} | frcmod=${frcmod_ok} | missing_params=${miss} | status=${status}"
    echo "${mol_id},${nc},${mol2_ok},${frcmod_ok},${miss},${status}" >> "${summary_file}"
  done < "${RMSD_CSV}"

  echo "mol_id | net_charge | mol2 | frcmod | missing_params | status"
  echo "-------------------------------------------------------------"
  while IFS=, read -r mol_id nc mol2_ok frcmod_ok miss status; do
    echo "${mol_id} | ${nc} | ${mol2_ok} | ${frcmod_ok} | ${miss} | ${status}"
  done < "${summary_file}"
  rm -f "${summary_file}"
  echo
  echo "${ok_count}/36 molecules successfully parameterized"
}

main "$@"

# Note: True antechamber RESP (-c resp) needs Gaussian gout/gesp; ORCA users get
# CHELPG ESP-fitted charges via resp.pc_chelpg (see apply_pc_chelpg_to_mol2.py).

