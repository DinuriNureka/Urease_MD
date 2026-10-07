#!/usr/bin/env bash
# scripts/run_multiwfn_esp_map.sh
#
# Generate colored ESP map on vdW surface (Multiwfn) → esp_map.png per molecule.
# May require X11; use xvfb-run if headless:  xvfb-run -a bash scripts/run_multiwfn_esp_map.sh
#
# Usage:
#   export MULTIWFN=/path/to/Multiwfn
#   bash scripts/run_multiwfn_esp_map.sh [resp_results_dir]
#
# Env:
#   MULTIWFN              Multiwfn executable
#   MULTIWFN_SKIP_ESP_MAP  set to 1 to skip (optional)

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RESULTS_DIR="${1:-${ROOT}/resp_results}"
MULTIWFN_BIN="${MULTIWFN:-Multiwfn}"

log_prefix="[run_multiwfn_esp_map]"
log_msg() { echo "[$(date '+%F %T')] ${log_prefix} $*"; }

if [[ "${MULTIWFN_SKIP_ESP_MAP:-0}" == "1" ]]; then
  log_msg "MULTIWFN_SKIP_ESP_MAP=1 — skipping ESP map generation."
  exit 0
fi

if ! command -v "${MULTIWFN_BIN}" >/dev/null 2>&1; then
  echo "${log_prefix} ERROR: Multiwfn not found. Set MULTIWFN." >&2
  exit 2
fi

RUNNER=( )
if command -v xvfb-run >/dev/null 2>&1; then
  RUNNER=(xvfb-run -a)
  log_msg "using xvfb-run for headless rendering"
fi

ok=0
skip=0
fail=0

for job_dir in "${RESULTS_DIR}"/mol*/; do
  [[ -d "${job_dir}" ]] || continue
  mol_id="$(basename "${job_dir}")"
  if [[ ! -f "${job_dir}/resp.molden" ]]; then
    skip=$((skip + 1))
    continue
  fi
  (
    cd "${job_dir}" || exit 1
    # Space function (5) -> Colored ESP on surface (12) -> vdW (1) -> map ESP (2) -> save (-1) -> exit
    "${RUNNER[@]}" "${MULTIWFN_BIN}" resp.molden <<'EOF' >esp_map_run.log 2>&1
5
12
1
2
-1
0
EOF
    # Multiwfn may write a default PNG name; collect the newest image
    shopt -s nullglob
    imgs=( *.png *.bmp )
    if [[ ${#imgs[@]} -gt 0 ]]; then
      latest="$(ls -t "${imgs[@]}" 2>/dev/null | head -1)"
      [[ -n "${latest}" ]] && cp -f "${latest}" esp_map.png
    fi
    if [[ -f esp_map.png ]]; then
      exit 0
    fi
    exit 1
  ) && {
    log_msg "${mol_id}: esp_map.png"
    ok=$((ok + 1))
  } || {
    log_msg "${mol_id}: ESP map FAILED (see ${job_dir}/esp_map_run.log)"
    fail=$((fail + 1))
  }
done

log_msg "done: ok=${ok} skip=${skip} fail=${fail}"
# ESP maps are best-effort (GUI/DISLIN); do not fail the script.
exit 0
