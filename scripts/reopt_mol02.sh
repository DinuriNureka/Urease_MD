#!/usr/bin/env bash
# Re-run ORCA Opt+Freq for mol02 with tighter thresholds (see mol02_reopt.inp).
# Uses the same MPI layout as params/inhibitors/orca_opt/molNN/ jobs (nprocs 6, maxcore 3000).
#
# Usage (from project root):
#   bash scripts/reopt_mol02.sh
#
# What it does:
#   1. Copies the current optimized geometry from params/.../mol02.xyz into resp_results/mol02/mol02.xyz
#   2. Runs ORCA on resp_results/mol02/mol02_reopt.inp → mol02_reopt.log, mol02_reopt.xyz, etc.
#   3. Parses frequencies; if there are no imaginary modes, replaces mol02.xyz (params + resp_results) with the new geometry
#
# Environment: ORCA_PATH from config.env (default project config)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/config.env"

ORCA_BIN="${ORCA_PATH:-$(command -v orca || true)}"
WORKDIR="${ROOT}/resp_results/mol02"
SRC_XYZ="${ROOT}/params/inhibitors/orca_opt/mol02/mol02.xyz"
INP="${WORKDIR}/mol02_reopt.inp"
LOG="${WORKDIR}/mol02_reopt.log"
OUT_XYZ="${WORKDIR}/mol02_reopt.xyz"
CANON_XYZ="${SRC_XYZ}"
RESP_XYZ="${WORKDIR}/mol02.xyz"

log() { echo "[reopt_mol02] $(date '+%F %T') $*"; }

[[ -x "${ORCA_BIN}" ]] || { echo "ERROR: ORCA not executable: ${ORCA_BIN}" >&2; exit 2; }
[[ -f "${INP}" ]] || { echo "ERROR: missing ${INP}" >&2; exit 2; }
[[ -f "${SRC_XYZ}" ]] || { echo "ERROR: missing ${SRC_XYZ}" >&2; exit 2; }

mkdir -p "${WORKDIR}"
cp -f "${SRC_XYZ}" "${RESP_XYZ}"
log "Starting geometry copied to ${RESP_XYZ} (from params)"

log "Running: ${ORCA_BIN} mol02_reopt.inp (cwd=${WORKDIR})"
(
  cd "${WORKDIR}" || exit 1
  "${ORCA_BIN}" mol02_reopt.inp > "${LOG}" 2>&1
)
rc=$?
if [[ ${rc} -ne 0 ]]; then
  log "ERROR: ORCA exited with code ${rc}"
  exit "${rc}"
fi

if ! grep -q 'TERMINATED NORMALLY' "${LOG}"; then
  log "ERROR: ${LOG} does not contain TERMINATED NORMALLY — not updating mol02.xyz"
  exit 3
fi

export REOPT_LOG="${LOG}"
export REOPT_OUT_XYZ="${OUT_XYZ}"
export REOPT_CANON_XYZ="${CANON_XYZ}"
export REOPT_RESP_XYZ="${RESP_XYZ}"

python3 - <<'PY'
import os
import re
import sys
from pathlib import Path

log_path = Path(os.environ["REOPT_LOG"])
out_xyz = Path(os.environ["REOPT_OUT_XYZ"])
canon_xyz = Path(os.environ["REOPT_CANON_XYZ"])
resp_xyz = Path(os.environ["REOPT_RESP_XYZ"])

text = log_path.read_text(errors="replace")
if "VIBRATIONAL FREQUENCIES" in text:
    text = text.split("VIBRATIONAL FREQUENCIES", 1)[1]

# Frequency lines: "     6:    -18.48 cm**-1  ***imaginary mode***"
pat = re.compile(
    r"^\s*(\d+)\s*:\s*([+-]?\d+\.?\d*)\s*cm\*\*-1\s*(.*)$",
    re.MULTILINE,
)
entries = []
for m in pat.finditer(text):
    rest = m.group(3).lower()
    val = float(m.group(2))
    imag = "imaginary" in rest or val < 0.0
    entries.append((m.group(1), val, imag))

has_imag = any(im for _, _, im in entries) or "***imaginary mode***" in text.lower()

pos_vib = [v for _, v, im in entries if not im and v > 1.0]
lowest = min(pos_vib) if pos_vib else float("nan")

print()
print("--- Frequency summary (mol02_reopt) ---")
n_imag = sum(1 for _, _, im in entries if im)
print(f"Lines parsed in frequency table: {len(entries)}")
print(f"Count with imaginary flag / negative: {n_imag}")
if entries:
    print("Sample (index : cm^-1, imaginary):", entries[:8], "...")
print(f"Lowest real vibrational frequency (>1 cm^-1): {lowest:.4f} cm^-1")
print("----------------------------------------")
print()

if has_imag:
    print("DECISION: Imaginary frequency still present — NOT overwriting mol02.xyz")
    sys.exit(4)

if not out_xyz.is_file():
    print(f"ERROR: expected optimized geometry at {out_xyz}", file=sys.stderr)
    sys.exit(5)

body = out_xyz.read_text()
canon_xyz.write_text(body)
resp_xyz.write_text(body)
print(f"Updated: {canon_xyz}")
print(f"Updated: {resp_xyz}")
print("DECISION: No imaginary modes — mol02.xyz replaced with reoptimized geometry.")
sys.exit(0)
PY
status=$?
exit "${status}"
