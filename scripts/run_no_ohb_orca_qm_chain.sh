#!/usr/bin/env bash
# After small_opt completes with HURRAY: copy xyz → NumFreq → hess→fchk → MCPB 2–4.
# Edit ORCA path and nprocs if needed. Do not run until opt has zero issues in .out.

set -euo pipefail
ROOT="/path/to/urease_md"
ORCA="${ORCA:-/opt/software/ORCA/build/orca}"
OPT_DIR="${ROOT}/params/ni_center_no_ohb/orca_calcs/small_opt"
FREQ_DIR="${ROOT}/params/ni_center_no_ohb/orca_calcs/small_freq"
LARGE_DIR="${ROOT}/params/ni_center_no_ohb/orca_calcs/large_mk"
MCPB_DIR="${ROOT}/params/ni_center_no_ohb/mcpb_input"
OPT_OUT="${OPT_DIR}/small_opt_wb97xd4.out"
OPT_XYZ="${OPT_DIR}/small_opt_wb97xd4.xyz"

export AMBERHOME="${AMBERHOME:-/usr/local}"
export AMBER_HOME="${AMBER_HOME:-/usr/local}"

grep -q "HURRAY" "${OPT_OUT}" || { echo "Optimization not converged: ${OPT_OUT}"; exit 1; }
NI_IN_XYZ=$(grep -cE '^[[:space:]]*Ni[[:space:]]' "${OPT_XYZ}" || true)
if [[ "${NI_IN_XYZ}" -lt 2 ]]; then
  echo "[FATAL] Fewer than 2 Ni in ${OPT_XYZ}"
  exit 1
fi

cp -a "${OPT_XYZ}" "${FREQ_DIR}/small_freq.xyz"
echo "Copied optimized geometry to ${FREQ_DIR}/small_freq.xyz"

cd "${FREQ_DIR}"
"${ORCA}" small_freq.inp > small_freq.out 2>&1
grep -q "ORCA TERMINATED NORMALLY" small_freq.out || { tail -40 small_freq.out; exit 1; }

python3 "${ROOT}/scripts/orca_hess_to_fchk.py" \
  "${FREQ_DIR}/small_freq.hess" \
  "${MCPB_DIR}/urease_no_ohb_small_opt.fchk"

cd "${MCPB_DIR}"
python3 "${ROOT}/scripts/run_mcpb_patched.py" -i urease_no_ohb.in -s 2 2>&1 | tee mcpb_step2.log
python3 "${ROOT}/scripts/run_mcpb_patched.py" -i urease_no_ohb.in -s 3 2>&1 | tee mcpb_step3.log
python3 "${ROOT}/scripts/run_mcpb_patched.py" -i urease_no_ohb.in -s 4 2>&1 | tee mcpb_step4.log

echo "Done. Validate frcmod bond lengths and imaginary modes before tleap."
