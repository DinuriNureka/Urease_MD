#!/usr/bin/env bash
set -uo pipefail

: <<'DOC'
run_resp_batch.sh

Purpose:
  Run ORCA 6.1.1 single-point calculations with CHELPG (for antechamber RESP)
  and electric-property output (ESP grid → resp.vpot when supported).

Inputs:
  - resp_inputs/job_manifest.csv
  - Optional flags: --manifest, --orca, --retry-failed, --dry-run

Outputs:
  - resp_results/molNN/{resp.inp,resp.out,resp.vpot?,DONE|FAILED}
  - resp_results/resp_batch.log

Example:
  bash run_resp_batch.sh
  bash run_resp_batch.sh --dry-run
  bash run_resp_batch.sh --retry-failed
DOC

ORCA=/opt/software/ORCA/build/orca
MAX_PARALLEL=6
NPROCS=8
MAXCORE=3500
MANIFEST=resp_inputs/job_manifest.csv
RESULTS_DIR=resp_results
LOG=${RESULTS_DIR}/resp_batch.log
RETRY_FAILED=0
DRY_RUN=0

log_msg() {
  local log_prefix="[run_resp_batch]"
  local msg="$1"
  echo "[$(date '+%F %T')] ${msg}" | tee -a "${LOG}"
}

parse_args() {
  local log_prefix="[parse_args]"
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --manifest)
        MANIFEST="${2:-}"
        shift 2
        ;;
      --orca)
        ORCA="${2:-}"
        shift 2
        ;;
      --retry-failed)
        RETRY_FAILED=1
        shift
        ;;
      --dry-run)
        DRY_RUN=1
        shift
        ;;
      *)
        echo "${log_prefix} ERROR: Unknown argument: $1" >&2
        exit 2
        ;;
    esac
  done
}

# Write ORCA input for a single molecule (same template as production jobs).
# Use out_path=- to print to stdout (dry-run).
write_resp_inp() {
  local out_path="$1"
  local mol_id="$2"
  local charge="$3"
  local mult="$4"
  # Note (ORCA 6.1.1): "GTTEN" and "Print[P_ESP]" are rejected by the parser.
  # %elprop + PrintLevel 3 still produces resp.vpot (grid) alongside resp.out.
  _resp_inp_heredoc() {
    cat <<EOF
! wB97X-D4 def2-SVP def2/J RIJCOSX TightSCF CHELPG
%maxcore ${MAXCORE}
%pal nprocs ${NPROCS} end
%elprop
  Dipole true
  Quadrupole true
end
%output
  PrintLevel 3
  Print[P_Hirshfeld] 1
  Print[P_Mulliken] 1
end
* xyzfile ${charge} ${mult} ${mol_id}.xyz
EOF
  }
  if [[ "${out_path}" == "-" ]]; then
    _resp_inp_heredoc
  else
    _resp_inp_heredoc > "${out_path}"
  fi
}

# Print the first resp.inp that would be generated (first manifest row).
dry_run_print_first_inp() {
  local log_prefix="[dry_run]"
  [[ -f "${MANIFEST}" ]] || {
    echo "${log_prefix} ERROR: manifest not found: ${MANIFEST}" >&2
    exit 2
  }
  local mol_id="" charge="" mult=""
  while IFS=, read -r mol_id _ _ charge mult _ val_status; do
    [[ "${mol_id}" == "mol_id" ]] && continue
    [[ "${val_status}" == "FAIL_FORMAT" ]] && continue
    break
  done < "${MANIFEST}"
  if [[ -z "${mol_id}" ]]; then
    echo "${log_prefix} ERROR: no valid rows in ${MANIFEST}" >&2
    exit 2
  fi
  echo "=== DRY RUN: showing first resp.inp that would be generated (${mol_id}) ==="
  write_resp_inp "-" "${mol_id}" "${charge}" "${mult}"
}

run_job() {
  local log_prefix="[run_job]"
  local mol_id="$1"
  local charge="$2"
  local mult="$3"
  local xyz_src="resp_inputs/${mol_id}.xyz"
  local job_dir="${RESULTS_DIR}/${mol_id}"
  local out_file="${job_dir}/resp.out"
  local inp_file="${job_dir}/resp.inp"

  mkdir -p "${job_dir}"
  cp -f "${xyz_src}" "${job_dir}/${mol_id}.xyz"

  if [[ -f "${job_dir}/DONE" ]]; then
    log_msg "${log_prefix} ${mol_id}: already done, skipping"
    return 0
  fi
  if [[ -f "${job_dir}/FAILED" && "${RETRY_FAILED}" -ne 1 ]]; then
    log_msg "${log_prefix} ${mol_id}: previously failed, skipping (use --retry-failed to retry)"
    return 0
  fi
  rm -f "${job_dir}/FAILED"

  write_resp_inp "${inp_file}" "${mol_id}" "${charge}" "${mult}"

  (
    cd "${job_dir}" || exit 1
    "${ORCA}" resp.inp > resp.out 2>&1
  ) || true

  if grep -q "FINAL SINGLE POINT ENERGY" "${out_file}"; then
    local e
    e="$(grep "FINAL SINGLE POINT ENERGY" "${out_file}" | tail -1 | awk '{print $5}')"
    touch "${job_dir}/DONE"
    rm -f "${job_dir}/FAILED"
    log_msg "${log_prefix} ${mol_id}: SCF CONVERGED, E = ${e} Ha"
    return 0
  fi

  touch "${job_dir}/FAILED"
  rm -f "${job_dir}/DONE"
  log_msg "${log_prefix} ${mol_id}: SCF FAILED — inspect ${job_dir}/resp.out"
  return 1
}

main() {
  local log_prefix="[main]"
  parse_args "$@"

  if [[ "${DRY_RUN}" -eq 1 ]]; then
    dry_run_print_first_inp
    exit 0
  fi

  mkdir -p "${RESULTS_DIR}"
  : > "${LOG}"

  [[ -f "${MANIFEST}" ]] || { echo "${log_prefix} ERROR: manifest not found: ${MANIFEST}" >&2; exit 2; }
  [[ -x "${ORCA}" ]] || { echo "${log_prefix} ERROR: ORCA not executable: ${ORCA}" >&2; exit 2; }

  local pids=()
  local total=0 done_n=0 failed_n=0 skipped_n=0

  while IFS=, read -r mol_id orig_file xyz_file charge mult n_atoms val_status; do
    [[ "${mol_id}" == "mol_id" ]] && continue
    total=$((total + 1))
    if [[ "${val_status}" == "FAIL_FORMAT" ]]; then
      skipped_n=$((skipped_n + 1))
      continue
    fi

    (
      run_job "${mol_id}" "${charge}" "${mult}"
    ) < /dev/null &
    pids+=($!)

    if [[ ${#pids[@]} -ge ${MAX_PARALLEL} ]]; then
      wait "${pids[0]}" || true
      pids=("${pids[@]:1}")
    fi
  done < "${MANIFEST}"
  wait || true

  done_n="$(ls -1 "${RESULTS_DIR}"/mol*/DONE 2>/dev/null | wc -l | tr -d ' ')"
  failed_n="$(ls -1 "${RESULTS_DIR}"/mol*/FAILED 2>/dev/null | wc -l | tr -d ' ')"

  {
    echo "============================================================"
    echo "RESP BATCH COMPLETE"
    echo "Total:   ${total}"
    echo "Done:    ${done_n}"
    echo "Failed:  ${failed_n}"
    echo "Skipped: ${skipped_n}"
    echo "============================================================"
  } | tee -a "${LOG}"
}

main "$@"
