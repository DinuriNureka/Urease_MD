#!/usr/bin/env bash
# scripts/run_multiwfn_resp.sh
#
# Run Multiwfn RESP fitting non-interactively on resp.molden in each molNN dir.
# Writes resp_multiwfn_resp.log, resp_multiwfn.chg, resp_multiwfn_metrics.txt
#
# Usage:
#   export MULTIWFN=/path/to/Multiwfn   # required
#   bash scripts/run_multiwfn_resp.sh [resp_results_dir]
#
# Env:
#   MULTIWFN  Multiwfn executable (default: Multiwfn in PATH)

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RESULTS_DIR="${1:-${ROOT}/resp_results}"
MULTIWFN_BIN="${MULTIWFN:-Multiwfn}"

log_prefix="[run_multiwfn_resp]"
log_msg() { echo "[$(date '+%F %T')] ${log_prefix} $*"; }

if ! command -v "${MULTIWFN_BIN}" >/dev/null 2>&1; then
  echo "${log_prefix} ERROR: Multiwfn not found. Set MULTIWFN to the binary path." >&2
  exit 2
fi

ok=0
skip=0
fail=0

for job_dir in "${RESULTS_DIR}"/mol*/; do
  [[ -d "${job_dir}" ]] || continue
  mol_id="$(basename "${job_dir}")"
  if [[ ! -f "${job_dir}/resp.molden" ]]; then
    log_msg "${mol_id}: skip (no resp.molden — run scripts/gbw_to_molden.sh first)"
    skip=$((skip + 1))
    continue
  fi

  (
    cd "${job_dir}" || exit 1
    # Multiwfn menu (user-specified): Population (7) -> RESP (18) -> ...
    "${MULTIWFN_BIN}" resp.molden <<'EOF' >resp_multiwfn_resp.log 2>&1
7
18
2
1
y
0
0
EOF
    python3 "${ROOT}/scripts/parse_multiwfn_resp.py" resp_multiwfn_resp.log \
      --out-chg resp_multiwfn.chg \
      --out-metrics resp_multiwfn_metrics.txt
  ) && {
    log_msg "${mol_id}: OK"
    ok=$((ok + 1))
  } || {
    log_msg "${mol_id}: FAILED"
    fail=$((fail + 1))
  }
done

log_msg "done: ok=${ok} skip=${skip} fail=${fail}"
exit 0
