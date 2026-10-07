#!/bin/bash
# scripts/check_dependencies.sh
# -------------------------------------------------------
# Verify all required software is installed and accessible.
# Run BEFORE starting any phase.
#
# Exit codes:
#   0 — all dependencies found
#   1 — one or more missing
# -------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config.env"

MISSING=0

check_cmd() {
    local name=$1
    local cmd=$2
    if command -v "$cmd" &>/dev/null; then
        local ver
        ver=$(timeout 5 "$cmd" --version 2>&1 | head -1 || echo "version unknown")
        printf "  %-20s %-10s %s\n" "$name" "[OK]" "$ver"
    else
        printf "  %-20s %-10s\n" "$name" "[MISSING]"
        MISSING=$((MISSING + 1))
    fi
}

check_file() {
    local name=$1
    local path=$2
    if [ -x "$path" ]; then
        printf "  %-20s %-10s %s\n" "$name" "[OK]" "$path"
    else
        printf "  %-20s %-10s expected at %s\n" "$name" "[MISSING]" "$path"
        MISSING=$((MISSING + 1))
    fi
}

echo "=== Software Dependencies ==="
echo ""
echo "--- AmberTools ---"
check_cmd "pdb4amber"     pdb4amber
check_cmd "tleap"         tleap
check_cmd "cpptraj"       cpptraj
check_cmd "antechamber"   antechamber
check_cmd "parmchk2"      parmchk2
check_cmd "MCPB.py"       MCPB.py
check_cmd "MMPBSA.py"     MMPBSA.py
check_cmd "parmed"        parmed
check_cmd "reduce"        reduce
check_cmd "ambmask"       ambmask
check_cmd "ante-MMPBSA.py" ante-MMPBSA.py

echo ""
echo "--- AMBER (licensed) ---"
check_cmd "pmemd.cuda"    pmemd.cuda
check_file "pmemd.cuda (project)" "${PMEMD_CUDA_BIN:-}"

echo ""
echo "--- QM Software ---"
check_file "ORCA 6.x"     "${ORCA_PATH}"
check_file "Multiwfn"     "${MULTIWFN_PATH}"

echo ""
echo "--- Utilities ---"
check_cmd "OpenBabel"     obabel
check_cmd "Python 3"      python3

echo ""
echo "--- GPU ---"
if nvidia-smi &>/dev/null; then
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
    GPU_MEM=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader | head -1)
    printf "  %-20s %-10s %s (%s)\n" "NVIDIA GPU" "[OK]" "$GPU_NAME" "$GPU_MEM"
else
    printf "  %-20s %-10s\n" "NVIDIA GPU" "[MISSING]"
    MISSING=$((MISSING + 1))
fi

echo ""
if [ "$MISSING" -eq 0 ]; then
    echo "[OK] All dependencies found."
    exit 0
else
    echo "[FAIL] ${MISSING} dependencies missing. Install them before proceeding."
    exit 1
fi
