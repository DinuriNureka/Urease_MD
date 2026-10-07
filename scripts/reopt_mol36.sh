#!/usr/bin/env bash
# Re-run ORCA Opt+Freq for mol36 (thiourea) from planar NH2 guess (mol36_planar.xyz).
# Same MPI as params/inhibitors/orca_opt/molNN/ (nprocs 6, maxcore 3000).
#
# Usage (from project root):
#   bash scripts/reopt_mol36.sh
#
# After a clean frequency run, replaces mol36.xyz in params and resp_results with
# mol36_reopt.xyz and writes resp_results/mol36/REDO_MULTIWFN_RESP.txt as a reminder
# to rerun Multiwfn RESP for this ligand only.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/config.env"

ORCA_BIN="${ORCA_PATH:-$(command -v orca || true)}"
WORKDIR="${ROOT}/resp_results/mol36"
PLANAR="${WORKDIR}/mol36_planar.xyz"
INP="${WORKDIR}/mol36_reopt.inp"
LOG="${WORKDIR}/mol36_reopt.log"
OUT_XYZ="${WORKDIR}/mol36_reopt.xyz"
CANON_XYZ="${ROOT}/params/inhibitors/orca_opt/mol36/mol36.xyz"
RESP_XYZ="${WORKDIR}/mol36.xyz"
FLAG="${WORKDIR}/REDO_MULTIWFN_RESP.txt"

log() { echo "[reopt_mol36] $(date '+%F %T') $*"; }

[[ -x "${ORCA_BIN}" ]] || { echo "ERROR: ORCA not executable: ${ORCA_BIN}" >&2; exit 2; }
[[ -f "${INP}" ]] || { echo "ERROR: missing ${INP}" >&2; exit 2; }
[[ -f "${PLANAR}" ]] || { echo "ERROR: missing ${PLANAR}" >&2; exit 2; }

mkdir -p "${WORKDIR}"

log "Running: ${ORCA_BIN} mol36_reopt.inp (cwd=${WORKDIR})"
(
  cd "${WORKDIR}" || exit 1
  "${ORCA_BIN}" mol36_reopt.inp > "${LOG}" 2>&1
)
rc=$?
if [[ ${rc} -ne 0 ]]; then
  log "ERROR: ORCA exited with code ${rc}"
  exit "${rc}"
fi

if ! grep -q 'TERMINATED NORMALLY' "${LOG}"; then
  log "ERROR: ${LOG} does not contain TERMINATED NORMALLY — not updating mol36.xyz"
  exit 3
fi

export REOPT_LOG="${LOG}"
export REOPT_OUT_XYZ="${OUT_XYZ}"
export REOPT_CANON_XYZ="${CANON_XYZ}"
export REOPT_RESP_XYZ="${RESP_XYZ}"
export REOPT_FLAG="${FLAG}"

python3 - <<'PY'
import os
import re
import sys
from pathlib import Path
from datetime import datetime, timezone

log_path = Path(os.environ["REOPT_LOG"])
out_xyz = Path(os.environ["REOPT_OUT_XYZ"])
canon_xyz = Path(os.environ["REOPT_CANON_XYZ"])
resp_xyz = Path(os.environ["REOPT_RESP_XYZ"])
flag_path = Path(os.environ["REOPT_FLAG"])

text = log_path.read_text(errors="replace")
if "VIBRATIONAL FREQUENCIES" in text:
    text = text.split("VIBRATIONAL FREQUENCIES", 1)[1]

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
print("--- Frequency summary (mol36_reopt) ---")
n_imag = sum(1 for _, _, im in entries if im)
print(f"Lines parsed in frequency table: {len(entries)}")
print(f"Count with imaginary flag / negative: {n_imag}")
if entries:
    print("Sample (index : cm^-1, imaginary):", entries[:8], "...")
print(f"Lowest real vibrational frequency (>1 cm^-1): {lowest:.4f} cm^-1")
print("----------------------------------------")
print()

if has_imag:
    print("DECISION: Imaginary frequency still present — NOT overwriting mol36.xyz")
    sys.exit(4)

if not out_xyz.is_file():
    print(f"ERROR: expected optimized geometry at {out_xyz}", file=sys.stderr)
    sys.exit(5)

body = out_xyz.read_text()
canon_xyz.write_text(body)
resp_xyz.write_text(body)
print(f"Updated: {canon_xyz}")
print(f"Updated: {resp_xyz}")

ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
flag_path.write_text(
    f"Geometry updated for mol36 (thiourea) from mol36_reopt at {ts}.\n"
    "Rerun Multiwfn RESP for this ligand only, then update charges in mol36.mol2:\n"
    "  - scripts/gbw_to_molden.sh (or only mol36 dir)\n"
    "  - scripts/run_multiwfn_resp.sh (or only mol36)\n"
    "  - scripts/apply_multiwfn_charges_batch.sh (or apply_multiwfn_charges.py for mol36)\n"
    "  - scripts/collect_multiwfn_rmse.py\n"
)
print(f"Wrote: {flag_path}")
print("DECISION: No imaginary modes — mol36.xyz replaced; Multiwfn RESP must be redone for mol36.")
sys.exit(0)
PY
status=$?
exit "${status}"
