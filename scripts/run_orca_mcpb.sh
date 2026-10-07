#!/bin/bash
# scripts/run_orca_mcpb.sh
# -------------------------------------------------------
# Run ORCA calculations for MCPB.py parameterization:
#   1. FM geometry optimization (mult=5, high-spin quintet)
#   2. BS-DFT frequency / Hessian (mult=1, BrokenSym 2,2)
#   3. BS-DFT ESP for RESP charges (mult=1, BrokenSym 2,2)
#
# The broken-symmetry approach captures the antiferromagnetic
# coupling between the two Ni2+ (d8) centers in urease.
#
# Usage:
#     bash scripts/run_orca_mcpb.sh
#
# Expected runtime: 6-24 hours total (dominated by Freq)
# -------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config.env"

ORCA="${ORCA_PATH}"
ORCA_DIR="$(dirname "${ORCA_PATH}")"
export PATH="${ORCA_DIR}:${PATH}"
export LD_LIBRARY_PATH="${ORCA_DIR}:${LD_LIBRARY_PATH:-}"

CALC_DIR="${PROJECT_ROOT}/params/ni_center/orca_calcs"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_DIR}/orca_mcpb.log"; }

# === Step 1: FM Geometry Optimization ===
log "=== STEP 1: FM geometry optimization (mult=5) ==="
cd "${CALC_DIR}/small_opt"

${ORCA} small_opt.inp > small_opt.out 2>&1

if grep -q "OPTIMIZATION HAS CONVERGED" small_opt.out; then
    log "[OK] FM optimization converged"
    ENERGY=$(grep "FINAL SINGLE POINT ENERGY" small_opt.out | tail -1 | awk '{print $5}')
    log "  Final energy: ${ENERGY} Eh"
else
    log "[FAIL] FM optimization did not converge. Check small_opt.out"
    exit 1
fi

# === Step 2: BS-DFT Frequency calculation (Hessian) ===
log "=== STEP 2: BS-DFT frequency calculation ==="
cd "${CALC_DIR}/small_freq"

# Copy optimized geometry and FM orbitals (BS starting point)
cp "${CALC_DIR}/small_opt/small_opt.xyz" small_freq.xyz
cp "${CALC_DIR}/small_opt/small_opt.gbw" small_opt.gbw

${ORCA} small_freq.inp > small_freq.out 2>&1

if [ -f "small_freq.hess" ]; then
    log "[OK] Hessian file generated: small_freq.hess"

    # Verify BS solution via <S^2>
    S2=$(grep "Expectation value of <S\*\*2>" small_freq.out | tail -1 | awk '{print $NF}')
    log "  <S^2> = ${S2:-N/A} (expect ~2.0 for BS-DFT AFM state)"

    N_IMAG=$(grep -c "***imaginary mode***" small_freq.out || true)
    if [ "$N_IMAG" -eq 0 ]; then
        log "  No imaginary frequencies (clean minimum)"
    elif [ "$N_IMAG" -le 3 ]; then
        log "  WARNING: ${N_IMAG} small imaginary modes (usually OK for parameterization)"
    else
        log "  WARNING: ${N_IMAG} imaginary modes — check geometry"
    fi
else
    log "[FAIL] No .hess file produced. Check small_freq.out"
    exit 1
fi

# === Step 3: BS-DFT Large model ESP ===
log "=== STEP 3: BS-DFT ESP for RESP charges ==="
cd "${CALC_DIR}/large_mk"

${ORCA} large_mk.inp > large_mk.out 2>&1

if grep -q "SCF CONVERGED" large_mk.out; then
    log "[OK] BS-DFT SCF converged"

    S2=$(grep "Expectation value of <S\*\*2>" large_mk.out | tail -1 | awk '{print $NF}')
    log "  <S^2> = ${S2:-N/A} (expect ~2.0 for BS-DFT AFM state)"

    ${ORCA_2MKL} large_mk -molden 2>/dev/null || true
    if [ -f "large_mk.molden.input" ]; then
        log "[OK] Molden file generated for Multiwfn RESP fitting"
    else
        log "[WARNING] Molden file not generated — run orca_2mkl manually"
    fi
else
    log "[FAIL] BS-DFT SCF did not converge. Check large_mk.out"
    exit 1
fi

log ""
log "=== ALL ORCA CALCULATIONS COMPLETE ==="
log "BS-DFT workflow preserves antiferromagnetic coupling."
log "Next: Convert Hessian to fchk format (Step 2.6)"
log "      Run RESP fitting with Multiwfn (Step 2.7)"
