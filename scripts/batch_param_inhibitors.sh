#!/bin/bash
# scripts/batch_param_inhibitors.sh
# -------------------------------------------------------
# GAFF2 parameterization for all inhibitors.
#
# Usage:
#     bash scripts/batch_param_inhibitors.sh
#
# Reads charges from params/inhibitors/charges.txt (optional):
#     mol01  0
#     aha   -1
#
# Outputs:
#     params/inhibitors/mol2/*_gaff2.mol2
#     params/inhibitors/frcmod/*.frcmod
# -------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config.env"

STRUCT_DIR="${PROJECT_ROOT}/params/inhibitors/structures"
MOL2_DIR="${PROJECT_ROOT}/params/inhibitors/mol2"
FRCMOD_DIR="${PROJECT_ROOT}/params/inhibitors/frcmod"
LOG_DIR_INH="${PROJECT_ROOT}/params/inhibitors/logs"
CHARGE_FILE="${PROJECT_ROOT}/params/inhibitors/charges.txt"

mkdir -p "${MOL2_DIR}" "${FRCMOD_DIR}" "${LOG_DIR_INH}"

log() { echo "[$(date '+%H:%M:%S')] $1"; }

get_charge() {
    local name=$1
    if [ -f "${CHARGE_FILE}" ]; then
        local charge
        charge=$(grep "^${name}[[:space:]]" "${CHARGE_FILE}" 2>/dev/null | awk '{print $2}')
        [ -n "${charge}" ] && echo "${charge}" && return
    fi
    echo "0"
}

SUCCESS=0; FAILED=0; WARN=0

log "=== Batch GAFF2 Parameterization ==="

for mol2 in "${STRUCT_DIR}"/*.mol2; do
    [ -f "$mol2" ] || continue
    NAME=$(basename "$mol2" .mol2)
    NC=$(get_charge "$NAME")

    log "Processing ${NAME} (charge=${NC})..."

    antechamber \
        -i "${mol2}" -fi mol2 \
        -o "${MOL2_DIR}/${NAME}_gaff2.mol2" -fo mol2 \
        -c bcc -at gaff2 -nc "${NC}" -rn LIG -s 0 \
        > "${LOG_DIR_INH}/${NAME}_ac.log" 2>&1

    if [ $? -ne 0 ]; then
        log "  [FAIL] antechamber failed for ${NAME}"
        FAILED=$((FAILED + 1))
        continue
    fi

    parmchk2 \
        -i "${MOL2_DIR}/${NAME}_gaff2.mol2" -f mol2 \
        -o "${FRCMOD_DIR}/${NAME}.frcmod" -s gaff2 \
        > "${LOG_DIR_INH}/${NAME}_pc.log" 2>&1

    if [ $? -ne 0 ]; then
        log "  [FAIL] parmchk2 failed for ${NAME}"
        FAILED=$((FAILED + 1))
        continue
    fi

    if grep -q "ATTN" "${FRCMOD_DIR}/${NAME}.frcmod" 2>/dev/null; then
        log "  [WARN] ${NAME}: has ATTN warnings in frcmod"
        WARN=$((WARN + 1))
    fi

    SUCCESS=$((SUCCESS + 1))

    rm -f ANTECHAMBER_*.AC* ATOMTYPE.INF sqm.* PREP.INF NEWPDB.PDB
done

log ""
log "=== Summary ==="
log "  Success: ${SUCCESS}"
log "  Failed:  ${FAILED}"
log "  ATTN warnings: ${WARN}"

if [ ${WARN} -gt 0 ]; then
    log ""
    log "Molecules needing manual review:"
    for frc in "${FRCMOD_DIR}"/*.frcmod; do
        if grep -q "ATTN" "$frc" 2>/dev/null; then
            name=$(basename "$frc" .frcmod)
            log "  ${name}: $(grep -c 'ATTN' "$frc") parameters"
        fi
    done
fi

[ ${FAILED} -eq 0 ] && exit 0 || exit 1
