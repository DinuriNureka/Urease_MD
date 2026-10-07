#!/usr/bin/env bash
# Extract RESP fitting RMSE from Multiwfn logs for all compounds
# Usage: bash scripts/check_resp_rmse.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

echo "RESP Fitting RMSE Summary"
echo "========================="
echo ""

RESP_DIR="${ROOT}/resp_results"
THRESH=0.01

if [[ ! -d "${RESP_DIR}" ]]; then
    echo "ERROR: ${RESP_DIR} not found"
    exit 1
fi

mapfile -t LOGS < <(find "${RESP_DIR}" -maxdepth 2 -type f -name 'resp_multiwfn_resp.log' | sort -V)

if [[ ${#LOGS[@]} -eq 0 ]]; then
    echo "No resp_multiwfn_resp.log files found under ${RESP_DIR}"
    exit 0
fi

printf "%-8s %-12s %-12s %s\n" "mol_id" "RMSE" "RRMSE" "flag"
echo "-------- ------------ ------------ ----"

FLAGGED=()
for f in "${LOGS[@]}"; do
    mol_id="$(basename "$(dirname "$f")")"
    line="$(grep -E '^\s*RMSE:' "$f" | tail -1 || true)"
    if [[ -z "${line}" ]]; then
        printf "%-8s %-12s %-12s %s\n" "${mol_id}" "N/A" "N/A" "no RMSE line"
        continue
    fi
    rmse="$(echo "${line}" | awk '{print $2}')"
    rrmse="$(echo "${line}" | awk '{print $4}')"
    flag="OK"
    if awk -v r="${rmse}" -v t="${THRESH}" 'BEGIN { exit (r > t) ? 0 : 1 }'; then
        flag="REVIEW"
        FLAGGED+=("${mol_id} (${rmse})")
    fi
    printf "%-8s %-12s %-12s %s\n" "${mol_id}" "${rmse}" "${rrmse}" "${flag}"
done

echo ""
echo "Project mean RMSE: ~0.004 e (previously approved by Sir)"
echo "Flag any compound with RMSE > ${THRESH} e for review"
if [[ ${#FLAGGED[@]} -gt 0 ]]; then
    echo ""
    echo "REVIEW list: ${FLAGGED[*]}"
fi
