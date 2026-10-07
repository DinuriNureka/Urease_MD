#!/bin/bash
# scripts/verify_phase1.sh
# -------------------------------------------------------
# Phase 1 completion checklist.
# Verifies all outputs from structure preparation.
#
# Usage:
#     bash scripts/verify_phase1.sh
#
# Exit codes:
#   0 — all checks passed
#   1 — one or more failed
# -------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config.env"

PASS=0
FAIL=0

check() {
    local desc=$1
    local result=$2
    if [ "$result" = "true" ]; then
        printf "  [PASS] %s\n" "$desc"
        PASS=$((PASS + 1))
    else
        printf "  [FAIL] %s\n" "$desc"
        FAIL=$((FAIL + 1))
    fi
}

echo "=== Phase 1 Completion Checklist ==="

PREP="${PROJECT_ROOT}/structures/prepared/${PDB_ID}_prepared.pdb"

check "Prepared PDB exists" \
    "$([ -s "$PREP" ] && echo true || echo false)"

CHAINS=$(grep "^ATOM" "$PREP" 2>/dev/null | awk '{print $5}' | sort -u | wc -l || echo 0)
check "24 chains present (found: ${CHAINS})" \
    "$([ "$CHAINS" -ge 20 ] && echo true || echo false)"

NI_COUNT=$(grep -c "NI" "$PREP" 2>/dev/null || echo 0)
check "24 Ni atoms present (found: ${NI_COUNT})" \
    "$([ "$NI_COUNT" -ge 24 ] && echo true || echo false)"

HIS_COUNT=$(grep "^ATOM" "$PREP" 2>/dev/null | awk '{print $4}' | grep -c "^HIS$" || true)
HIS_COUNT=${HIS_COUNT:-0}
HIS_COUNT=$(echo "$HIS_COUNT" | tr -d '[:space:]')
check "No unassigned HIS residues (found: ${HIS_COUNT})" \
    "$([ "$HIS_COUNT" -eq 0 ] && echo true || echo false)"

check "Renumbering map exists" \
    "$([ -f "${PROJECT_ROOT}/structures/cleaned/${PDB_ID}_clean_renum.txt" ] && echo true || echo false)"

echo ""
echo "Results: ${PASS} passed, ${FAIL} failed"
if [ "$FAIL" -eq 0 ]; then
    echo "[OK] Phase 1 complete. Proceed to Phase 2."
    exit 0
else
    echo "[BLOCKED] Fix failures before proceeding."
    exit 1
fi
