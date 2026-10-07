#!/usr/bin/env bash
# Ni-center ORCA chain (sequential, unattended):
#   small_opt_wB97XD4 -> small_freq_wB97XD4 -> large_mk_wB97XD4
#
# Default: wait for an already-running small_opt (manual nohup) to finish, verify
# success, then run freq and large_mk. Use after: nohup orca small_opt_wB97XD4.inp ...
#
# --from-scratch: run small_opt first (blocking), then freq, then large_mk.
#   Do not use while another small_opt is already running.
#
# Environment:
#   NI_CENTER_WAIT_PID  — optional explicit PID to wait on (main orca driver)
#
# Logs: ${PROJECT_ROOT}/logs/ni_center_chain.log (tee)
#
# Usage:
#   bash scripts/run_ni_center_chain.sh
#   bash scripts/run_ni_center_chain.sh --from-scratch
#   NI_CENTER_WAIT_PID=620247 bash scripts/run_ni_center_chain.sh
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "${ROOT}/config.env"

LOG="${ROOT}/logs/ni_center_chain.log"
mkdir -p "${ROOT}/logs"
# When launched from run_all_qm.sh, RUN_ALL_QM_CHILD=1 — skip internal tee; parent tees to NI + this log.
if [[ "${RUN_ALL_QM_CHILD:-0}" != "1" ]]; then
  exec > >(tee -a "${LOG}") 2>&1
fi

log() { echo "[$(date '+%F %T')] $*"; }

status_line() {
  # Optional: written when Ni-center is driven from run_all_qm.sh (NI_STATUS_FILE set).
  [[ -n "${NI_STATUS_FILE:-}" ]] || return 0
  echo "$1" >> "${NI_STATUS_FILE}"
}

ORCA_BIN="${ORCA_PATH:-}"
[[ -n "${ORCA_BIN}" && -x "${ORCA_BIN}" ]] || { log "ERROR: ORCA not executable (ORCA_PATH='${ORCA_PATH:-}')"; exit 2; }

ORCA_DIR="$(dirname "${ORCA_BIN}")"
export PATH="${ORCA_DIR}:${PATH}"
export LD_LIBRARY_PATH="${ORCA_DIR}:${LD_LIBRARY_PATH:-}"

NI_CALC_DIR="${ROOT}/params/ni_center/orca_calcs"
SMALL_OPT_DIR="${NI_CALC_DIR}/small_opt"
SMALL_FREQ_DIR="${NI_CALC_DIR}/small_freq"
LARGE_MK_DIR="${NI_CALC_DIR}/large_mk"

SMALL_OPT_INP="small_opt_wB97XD4.inp"
SMALL_OPT_OUT="small_opt_wB97XD4.out"
SMALL_FREQ_INP="small_freq_wB97XD4.inp"
SMALL_FREQ_OUT="small_freq_wB97XD4.out"
LARGE_MK_INP="large_mk_wB97XD4.inp"
LARGE_MK_OUT="large_mk_wB97XD4.out"

FROM_SCRATCH=0
for arg in "$@"; do
  case "$arg" in
    --from-scratch) FROM_SCRATCH=1 ;;
    -h|--help)
      echo "Usage: $0 [--from-scratch]"
      exit 0
      ;;
  esac
done

small_opt_failed_in_log() {
  local out="$1"
  [[ -f "$out" ]] || return 1
  grep -q "ORCA finished by error termination" "$out" 2>/dev/null \
    || grep -q "error termination in" "$out" 2>/dev/null
}

small_opt_succeeded_in_log() {
  local out="$1"
  [[ -f "$out" ]] || return 1
  grep -q "HURRAY" "$out" 2>/dev/null \
    || grep -q "THE OPTIMIZATION HAS CONVERGED" "$out" 2>/dev/null \
    || grep -q "OPTIMIZATION RUN DONE" "$out" 2>/dev/null
}

detect_small_opt_orca_pid() {
  # Main ORCA driver (not orca_leanscf_mpi workers)
  pgrep -f 'orca small_opt_wB97XD4.inp' 2>/dev/null | head -1 || true
}

wait_for_small_opt_completion() {
  local out="${SMALL_OPT_DIR}/${SMALL_OPT_OUT}"
  local pid="${NI_CENTER_WAIT_PID:-}"
  if [[ -z "${pid}" ]]; then
    pid="$(detect_small_opt_orca_pid)"
  fi

  log "Waiting for small_opt to finish..."
  if [[ -n "${pid}" ]]; then
    log "Monitoring ORCA driver PID: ${pid} (file: ${out})"
  else
    log "No orca small_opt_wB97XD4.inp PID found; polling ${out} only"
  fi

  local spin=0
  while true; do
    if small_opt_failed_in_log "${out}"; then
      log "ERROR: small_opt failed (see ${out})"
      return 1
    fi
    if small_opt_succeeded_in_log "${out}"; then
      log "small_opt success markers present in ${out}"
      # Ensure driver has exited so GBW is fully written
      if [[ -n "${pid}" ]] && kill -0 "${pid}" 2>/dev/null; then
        log "Waiting for PID ${pid} to exit..."
        while kill -0 "${pid}" 2>/dev/null; do sleep 15; done
      fi
      sleep 3
      return 0
    fi

    if [[ -n "${pid}" ]]; then
      if ! kill -0 "${pid}" 2>/dev/null; then
        log "PID ${pid} no longer running; verifying output..."
        sleep 5
        if small_opt_failed_in_log "${out}"; then
          log "ERROR: small_opt failed after process exit"
          return 1
        fi
        if small_opt_succeeded_in_log "${out}"; then
          log "small_opt completed successfully"
          return 0
        fi
        log "ERROR: process exited without success markers in ${out}"
        return 1
      fi
    else
      # No PID: poll file until we see completion or stall forever — user should start from --from-scratch
      if [[ -f "${out}" ]] && [[ $(stat -c%s "${out}" 2>/dev/null || echo 0) -gt 0 ]]; then
        :
      fi
    fi

    spin=$((spin + 1))
    if (( spin % 10 == 0 )); then
      log "... still waiting (small_opt) ..."
    fi
    sleep 30
  done
}

run_small_opt_blocking() {
  log "[NI] step 1: running ${SMALL_OPT_INP} (blocking)"
  cd "${SMALL_OPT_DIR}"
  if pgrep -f 'orca small_opt_wB97XD4.inp' >/dev/null 2>&1; then
    log "ERROR: another small_opt_wB97XD4.inp is already running. Do not use --from-scratch."
    log "Hint: run this script without --from-scratch to wait for the current job."
    exit 1
  fi
  "${ORCA_BIN}" "${SMALL_OPT_INP}" > "${SMALL_OPT_OUT}" 2>&1
  if small_opt_failed_in_log "${SMALL_OPT_DIR}/${SMALL_OPT_OUT}"; then
    log "ERROR: small_opt failed"
    exit 1
  fi
  if ! small_opt_succeeded_in_log "${SMALL_OPT_DIR}/${SMALL_OPT_OUT}"; then
    log "ERROR: small_opt finished without convergence markers"
    exit 1
  fi
  if [[ ! -f "small_opt.gbw" && ! -f "small_opt_wB97XD4.gbw" ]]; then
    log "ERROR: small_opt missing GBW output"
    status_line "small_opt:FAIL"
    exit 1
  fi
  log "[NI] small_opt OK"
  status_line "small_opt:OK"
}

copy_opt_to_small_freq() {
  log "[NI] copying optimized geometry and GBW to small_freq/"
  cd "${SMALL_FREQ_DIR}"
  if [[ -f "${SMALL_OPT_DIR}/small_opt_wB97XD4.xyz" ]]; then
    cp -f "${SMALL_OPT_DIR}/small_opt_wB97XD4.xyz" ./small_freq.xyz
    log "Using small_opt_wB97XD4.xyz -> small_freq/small_freq.xyz"
  else
    cp -f "${SMALL_OPT_DIR}/small_opt.xyz" ./small_freq.xyz
    log "Using small_opt.xyz -> small_freq/small_freq.xyz"
  fi
  if [[ -f "${SMALL_OPT_DIR}/small_opt.gbw" ]]; then
    cp -f "${SMALL_OPT_DIR}/small_opt.gbw" ./small_opt.gbw
  elif [[ -f "${SMALL_OPT_DIR}/small_opt_wB97XD4.gbw" ]]; then
    cp -f "${SMALL_OPT_DIR}/small_opt_wB97XD4.gbw" ./small_opt.gbw
    log "Copied small_opt_wB97XD4.gbw -> small_freq/small_opt.gbw"
  else
    log "ERROR: no GBW in small_opt to copy"
    exit 1
  fi
}

freq_succeeded() {
  local dir="$1"
  [[ -f "${dir}/small_freq.hess" || -f "${dir}/small_freq_wB97XD4.hess" ]] || return 1
  local out="${dir}/${SMALL_FREQ_OUT}"
  [[ -f "${out}" ]] || return 1
  grep -q "ORCA TERMINATED NORMALLY" "${out}" 2>/dev/null \
    || grep -q "VIBRATIONAL FREQUENCIES" "${out}" 2>/dev/null
}

run_small_freq() {
  log "[NI] step 2: ${SMALL_FREQ_INP}"
  cd "${SMALL_FREQ_DIR}"
  "${ORCA_BIN}" "${SMALL_FREQ_INP}" > "${SMALL_FREQ_OUT}" 2>&1
  if ! freq_succeeded "${SMALL_FREQ_DIR}"; then
    log "ERROR: small_freq failed (check ${SMALL_FREQ_DIR}/${SMALL_FREQ_OUT} and .hess)"
    exit 1
  fi
  log "[NI] small_freq OK"
}

large_mk_succeeded() {
  local dir="$1"
  local out="${dir}/${LARGE_MK_OUT}"
  [[ -f "${out}" ]] || return 1
  grep -q "SCF CONVERGED" "${out}" && grep -q "ORCA TERMINATED NORMALLY" "${out}"
}

run_large_mk() {
  log "[NI] step 3: ${LARGE_MK_INP}"
  cd "${LARGE_MK_DIR}"
  "${ORCA_BIN}" "${LARGE_MK_INP}" > "${LARGE_MK_OUT}" 2>&1
  if ! large_mk_succeeded "${LARGE_MK_DIR}"; then
    log "ERROR: large_mk failed (check ${LARGE_MK_DIR}/${LARGE_MK_OUT})"
    status_line "large_mk:FAIL"
    exit 1
  fi
  log "[NI] large_mk OK (SCF converged, normal termination)"
  status_line "large_mk:OK"
}

# --- main ---
log "======== Ni-center chain start (FROM_SCRATCH=${FROM_SCRATCH}) ========"

if [[ "${FROM_SCRATCH}" -eq 1 ]]; then
  run_small_opt_blocking
else
  if [[ -f "${SMALL_OPT_DIR}/${SMALL_OPT_OUT}" ]] && small_opt_succeeded_in_log "${SMALL_OPT_DIR}/${SMALL_OPT_OUT}"; then
    log "small_opt already completed successfully; skipping wait."
  else
    wait_for_small_opt_completion || { status_line "small_opt:FAIL"; exit 1; }
  fi
  if ! small_opt_succeeded_in_log "${SMALL_OPT_DIR}/${SMALL_OPT_OUT}"; then
    log "ERROR: small_opt output does not show success"
    status_line "small_opt:FAIL"
    exit 1
  fi
  if [[ ! -f "${SMALL_OPT_DIR}/small_opt.gbw" && ! -f "${SMALL_OPT_DIR}/small_opt_wB97XD4.gbw" ]]; then
    log "ERROR: missing GBW after small_opt"
    status_line "small_opt:FAIL"
    exit 1
  fi
  status_line "small_opt:OK"
fi

copy_opt_to_small_freq
run_small_freq
run_large_mk

log "======== Ni-center chain finished successfully ========"
exit 0
