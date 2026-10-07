#!/usr/bin/env bash
# Master QM orchestration:
# - Ni-center ORCA workflow (sequential): small_opt -> small_freq -> large_mk
# - Ligand ORCA workflow (parallel batch): 6 jobs at a time (6 cores/job)
# Both tracks run concurrently and share total machine cores.
#
# Core budgeting (requested):
# - Ni-center track: 12 cores (maxcore 4000 MB each job)
# - Ligand track:    6 jobs * 6 cores = 36 cores
# Total: 48 cores, leaving headroom on 56 cores.
#
# Logging:
# - Master log
# - Ni-center track log
# - Ligand track log
# - Per-job ORCA logs remain in their own directories
#
# Usage:
#   bash scripts/run_all_qm.sh
#   SKIP_LIGAND=1 bash scripts/run_all_qm.sh   # Ni-center ORCA only (no ligand batch)
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "${ROOT}/config.env"

STAMP="$(date '+%Y-%m-%d_%H%M%S')"
RUN_DIR="${ROOT}/logs/qm_run_all_${STAMP}"
mkdir -p "${RUN_DIR}"

MASTER_LOG="${RUN_DIR}/master.log"
NI_LOG="${RUN_DIR}/ni_center_track.log"
LIG_LOG="${RUN_DIR}/ligands_orca_track.log"

exec > >(tee -a "${MASTER_LOG}") 2>&1

log() { echo "[$(date '+%F %T')] $*"; }

MAP_CSV="${ROOT}/params/inhibitors/orca_opt/ligand_name_mapping.csv"
ORCA_BIN="${ORCA_PATH:-$(command -v orca || true)}"

[[ -f "${MAP_CSV}" ]] || { log "ERROR: missing mapping CSV: ${MAP_CSV}"; exit 2; }
[[ -n "${ORCA_BIN}" && -x "${ORCA_BIN}" ]] || { log "ERROR: ORCA executable not found (ORCA_PATH='${ORCA_PATH:-}')."; exit 2; }

log "ORCA binary: ${ORCA_BIN}"
log "Mapping CSV: ${MAP_CSV}"
log "Run dir: ${RUN_DIR}"

ORCA_DIR="$(dirname "${ORCA_BIN}")"
export PATH="${ORCA_DIR}:${PATH}"
export LD_LIBRARY_PATH="${ORCA_DIR}:${LD_LIBRARY_PATH:-}"

MAX_LIG_JOBS="${ORCA_LIG_OPT_MAX_JOBS:-6}"   # 6 at a time
NI_CORES="${ORCA_NI_NPROCS:-12}"            # 12 cores for Ni-center track
NI_MAXCORE_MB="${ORCA_NI_MAXCORE_MB:-4000}" # 4000 MB per ORCA job

NI_CALC_DIR="${ROOT}/params/ni_center/orca_calcs"

ni_track() {
  (
    set -euo pipefail
    # One copy to ni_center_track.log; master log still receives this subshell's stdout.
    exec > >(tee -a "${NI_LOG}") 2>&1

    log "[NI] starting Ni-center track (wB97X-D4) in background"
    TMP_TAG="run_all_qm_${STAMP}"

    # Ensure clean status artifacts
    NI_STATUS_FILE="${RUN_DIR}/ni_status_${TMP_TAG}.txt"
    : > "${NI_STATUS_FILE}"

    # Delegates to scripts/run_ni_center_chain.sh --from-scratch (small_opt -> small_freq -> large_mk).
    # Also mirrors ORCA progress into logs/ni_center_chain.log.
    export NI_STATUS_FILE
    export RUN_ALL_QM_CHILD=1
    bash "${ROOT}/scripts/run_ni_center_chain.sh" --from-scratch 2>&1 | tee -a "${ROOT}/logs/ni_center_chain.log"

    log "[NI] Ni-center track finished successfully"
  )
}

ligand_track() {
  (
    set -euo pipefail
    log "[LIG] starting ligand Opt+Freq batch track in background"

    # Run exactly the existing per-ligand batch runner, but cap concurrency to 6.
    # Each ORCA job uses %pal nprocs 6 inside the per-ligand inputs.
    export ORCA_OPT_MAX_JOBS="${MAX_LIG_JOBS}"
    bash "${ROOT}/scripts/run_orca_opt_all_ligands.sh"

    log "[LIG] ligand track finished (script exit code=$?)"
  ) 2>&1 | tee -a "${LIG_LOG}"
}

if [[ "${SKIP_LIGAND:-0}" == "1" ]]; then
  log "SKIP_LIGAND=1 — Ni-center track only (no ligand batch)"
  ni_track
  NI_EXIT=$?
  log "Ni-center track exit code: ${NI_EXIT}"
  LIG_PID=""
else
  log "Launching both tracks concurrently..."

  ni_track &
  NI_PID=$!

  ligand_track &
  LIG_PID=$!

  wait "${NI_PID}" || true
  wait "${LIG_PID}" || true
fi

log "Run finished. Collecting success/failure..."

REPORT_FILE="${RUN_DIR}/success_report.txt"
: > "${REPORT_FILE}"

echo "Ni-center track (small_opt -> small_freq -> large_mk):" >> "${REPORT_FILE}"
for f in "${RUN_DIR}"/ni_status_run_all_qm_${STAMP}.txt; do
  [[ -f "$f" ]] || continue
  cat "$f" >> "${REPORT_FILE}"
done

echo "" >> "${REPORT_FILE}"
echo "Ligand Opt+Freq track (per molNN):" >> "${REPORT_FILE}"

ok=()
fail=()

while IFS=, read -r mol_id _rest1 _rest2 _rest3 _rest4 _rest5 _rest6; do
  [[ "${mol_id}" == "mol_id" ]] && continue
  lig_log="${ROOT}/params/inhibitors/orca_opt/${mol_id}/${mol_id}.log"
  lig_gbw="${ROOT}/resp_results/${mol_id}/resp.gbw"

  if [[ -f "${lig_gbw}" ]] && [[ -f "${lig_log}" ]] && grep -q "TERMINATED NORMALLY" "${lig_log}"; then
    ok+=("${mol_id}")
  else
    fail+=("${mol_id}")
  fi
done < <(tail -n +2 "${MAP_CSV}")

echo "succeeded: ${#ok} / ${#ok[@]} (expected 36)" >> "${REPORT_FILE}"
echo "failed: ${#fail} (failed mol_ids follow)" >> "${REPORT_FILE}"
if [[ ${#fail[@]} -gt 0 ]]; then
  printf 'fail_ids: %s\n' "${fail[*]}" >> "${REPORT_FILE}"
fi

log "Success report written to: ${REPORT_FILE}"
cat "${REPORT_FILE}"

exit 0

