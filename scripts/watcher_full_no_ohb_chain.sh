#!/usr/bin/env bash
# Wait for tighter small_opt HURRAY, then freq → CHELPG → hess2fchk → MCPB 2–4.
# ORCA inputs use nprocs 32 / %maxcore 1800 (see params/ni_center_no_ohb/orca_calcs/*/*.inp).
set -euo pipefail

ROOT="/path/to/urease_md"
ORCA="${ORCA:-/opt/software/ORCA/build/orca}"
# Open MPI default slots can be 28 (cores) while PAL nprocs=32 — env below avoids startup abort on NumFreq/CHELPG.
OPT_OUT="${ROOT}/params/ni_center_no_ohb/orca_calcs/small_opt/small_opt_wb97xd4.out"
OPT_XYZ="${ROOT}/params/ni_center_no_ohb/orca_calcs/small_opt/small_opt_wb97xd4.xyz"
FREQ_DIR="${ROOT}/params/ni_center_no_ohb/orca_calcs/small_freq"
LARGE_DIR="${ROOT}/params/ni_center_no_ohb/orca_calcs/large_mk"
MCPB_DIR="${ROOT}/params/ni_center_no_ohb/mcpb_input"

export AMBERHOME="${AMBERHOME:-/usr/local}"
export AMBER_HOME="${AMBER_HOME:-/usr/local}"

echo "Watcher started at $(date)"

if grep -q "HURRAY" "${OPT_OUT}" 2>/dev/null; then
    echo "=== HURRAY already present in ${OPT_OUT} — skipping optimization wait ==="
else
    while ! grep -q "HURRAY" "${OPT_OUT}" 2>/dev/null; do
        CYCLES=$(grep -c "GEOMETRY OPTIMIZATION CYCLE" "${OPT_OUT}" 2>/dev/null || echo 0)
        echo "$(date +%H:%M) — cycle ${CYCLES}"
        sleep 60
    done
    echo "=== HURRAY at $(date) ==="
fi

NI_COUNT=$(grep -cE '^[[:space:]]*Ni[[:space:]]' "${OPT_XYZ}" 2>/dev/null || echo 0)
echo "Ni atoms in optimized xyz: ${NI_COUNT}"
if [[ "${NI_COUNT}" -lt 2 ]]; then
    echo "[FATAL] Fewer than 2 Ni in ${OPT_XYZ}"
    exit 1
fi

if [[ ! -f "${FREQ_DIR}/small_freq.xyz" ]]; then
    cp -a "${OPT_XYZ}" "${FREQ_DIR}/small_freq.xyz"
    echo "Geometry copied to freq (small_freq.xyz was missing)"
else
    echo "small_freq.xyz present — not overwriting (NumFreq may be running or geometry already staged)"
fi

echo "=== Waiting for NumFreq: ${FREQ_OUT} (ORCA TERMINATED NORMALLY) ==="
while true; do
    if grep -q "ORCA TERMINATED NORMALLY" "${FREQ_OUT}" 2>/dev/null; then
        echo "=== NumFreq completed $(date) ==="
        break
    fi
    if grep -qE "ORCA finished by error termination|error termination in SCF gradient" "${FREQ_OUT}" 2>/dev/null; then
        echo "=== NumFreq FAILED (see ${FREQ_OUT}) ==="
        tail -30 "${FREQ_OUT}" || true
        exit 1
    fi
    DISP=$(grep -ciE "DISPLACED GEOMETRY|displaced geometry" "${FREQ_OUT}" 2>/dev/null || echo 0)
    echo "$(date +%H:%M) — NumFreq pending or running (${DISP} displacement line(s) so far)"
    sleep 120
done

NIMAG=$(grep -c "\*\*\*imaginary mode\*\*\*" "${FREQ_OUT}" 2>/dev/null || true)
NIMAG=${NIMAG:-0}
echo "Imaginary frequencies: ${NIMAG}"
if [[ "${NIMAG}" -gt 0 ]]; then
    grep "\*\*\*imaginary mode\*\*\*" "${FREQ_OUT}" || true
    echo '[WARN] Imaginary modes present — review coordinates'
fi

echo "=== CHELPG starting $(date) ==="
cd "${LARGE_DIR}"
env OMPI_MCA_rmaps_base_oversubscribe=1 OMPI_MCA_hwloc_base_binding_policy=none \
    "${ORCA}" large_mk.inp > large_mk.out 2>&1
if grep -q "ORCA TERMINATED NORMALLY" large_mk.out; then
    echo "=== CHELPG done $(date) ==="
else
    echo "=== CHELPG FAILED ===" && tail -20 large_mk.out && exit 1
fi

echo "=== Converting hess to fchk ==="
cd "${ROOT}"
python3 scripts/orca_hess_to_fchk.py \
    "${FREQ_DIR}/small_freq.hess" \
    "${MCPB_DIR}/urease_no_ohb_small_opt.fchk"

python3 - << PY
from pathlib import Path
import re
p = Path("${MCPB_DIR}/urease_no_ohb_small_opt.fchk")
text = p.read_text()
m = re.search(r"Atomic numbers\s+I\s+N=\s*(\d+)\s*\n((?:\s*\d+)+)", text)
if not m:
    print("Could not parse atomic numbers from fchk")
else:
    nums = [int(x) for x in m.group(2).split()]
    z28 = nums.count(28)
    print(f"Z=28 (Ni) count in fchk: {z28} (expect 2)")
PY

echo "=== MCPB.py step 2 ==="
cd "${MCPB_DIR}"
python3 "${ROOT}/scripts/run_mcpb_patched.py" -i urease_no_ohb.in -s 2 2>&1 | tee mcpb_step2.log

echo "=== MCPB.py step 3 ==="
python3 "${ROOT}/scripts/run_mcpb_patched.py" -i urease_no_ohb.in -s 3 2>&1 | tee mcpb_step3.log

echo "=== MCPB.py step 4 ==="
python3 "${ROOT}/scripts/run_mcpb_patched.py" -i urease_no_ohb.in -s 4 2>&1 | tee mcpb_step4.log

echo "=== Checking frcmod ==="
if grep "Y4" "${MCPB_DIR}/urease_no_ohb_mcpbpy.frcmod" 2>/dev/null; then
    echo "WARNING: Y4 still present"
else
    echo "OK: No Y4"
fi

echo "=== Bond / Seminario excerpt from step 2 ==="
grep -E "bond|Bond|angle|M1|M2|Ni|length|Angstrom|Y[0-9]-M" "${MCPB_DIR}/mcpb_step2.log" | head -50 || true

echo "=== ALL DONE at $(date) ==="
