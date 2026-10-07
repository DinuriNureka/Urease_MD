#!/usr/bin/env bash
set -euo pipefail

# Build and compare 9 Å buffer: truncated oct (solvateOct) vs rectangular (solvateBox).
# Writes only under params/box_comparison/ (does not modify simulations/ or systems/).

PROJECT_ROOT="/path/to/urease_md"
cd "${PROJECT_ROOT}"

BOX_DIR="${PROJECT_ROOT}/params/box_comparison"
REPORT="${BOX_DIR}/box_comparison_9A_report.txt"

TLEAP_OCT="${BOX_DIR}/tleap_oct_9A.in"
TLEAP_RECT="${BOX_DIR}/tleap_rect_9A.in"

OCT_TOP="${BOX_DIR}/oct_9A.prmtop"
OCT_CRD="${BOX_DIR}/oct_9A.inpcrd"
RECT_TOP="${BOX_DIR}/rect_9A.prmtop"
RECT_CRD="${BOX_DIR}/rect_9A.inpcrd"

OCT_LOG="${BOX_DIR}/tleap_oct_9A.log"
RECT_LOG="${BOX_DIR}/tleap_rect_9A.log"

require_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing command: $1" >&2; exit 2; }; }
require_cmd tleap

AMBPDB="$(command -v ambpdb || command -v bashambpdb || true)"
[[ -n "${AMBPDB}" ]] || { echo "Missing ambpdb/bashambpdb on PATH." >&2; exit 2; }

pdb_from_topcrd() {
  local ptop="$1" crd="$2" outpdb="$3"
  "${AMBPDB}" -p "${ptop}" -c "${crd}" > "${outpdb}"
}

run_tleap() {
  local in="$1" outlog="$2"
  echo "Running: tleap -f ${in}"
  tleap -f "${in}" > "${outlog}" 2>&1
  if grep -qiE "FATAL|fatal" "${outlog}"; then
    echo "tleap reported FATAL in ${outlog}" >&2
    exit 1
  fi
}

count_wat_from_tleap_log() {
  local lg="$1"
  awk 'toupper($1)=="WAT" && $2 ~ /^[0-9]+$/ {n=$2} END{print n+0}' "$lg"
}

natom_from_prmtop() {
  python3 - "$1" <<'PY'
import sys
from pathlib import Path
ptop = Path(sys.argv[1])
lines = ptop.read_text(errors="ignore").splitlines()
for i, ln in enumerate(lines):
    if ln.strip().startswith("%FLAG POINTERS"):
        j = i + 2
        ints = []
        while j < len(lines) and not lines[j].startswith("%FLAG"):
            ints.extend(lines[j].split())
            if len(ints) >= 1:
                break
            j += 1
        if ints:
            print(ints[0])
            sys.exit(0)
print("NA")
PY
}

# Box: a b c alpha beta gamma from prmtop+inpcrd (Å, degrees) via parmed
box_summary() {
  python3 - "$1" "$2" <<'PY'
import sys
try:
    import parmed as pmd
except ImportError:
    print("NA (parmed not available)")
    sys.exit(0)
import numpy as np
p = pmd.load_file(sys.argv[1], sys.argv[2])
if getattr(p, "box", None) is None:
    print("NA")
else:
    b = p.box
    if hasattr(b, "a"):
        print(
            f"a={b.a:.4f} Å  b={b.b:.4f} Å  c={b.c:.4f} Å  "
            f"alpha={b.alpha:.4f}°  beta={b.beta:.4f}°  gamma={b.gamma:.4f}°"
        )
    elif isinstance(b, np.ndarray) and b.size >= 6:
        a, b_, c, al, be, ga = b.flat[:6]
        print(
            f"a={float(a):.4f} Å  b={float(b_):.4f} Å  c={float(c):.4f} Å  "
            f"alpha={float(al):.4f}°  beta={float(be):.4f}°  gamma={float(ga):.4f}°"
        )
    else:
        print(str(b))
PY
}

echo "=== tleap: octahedral 9 Å ==="
run_tleap "${TLEAP_OCT}" "${OCT_LOG}"
[[ -f "${OCT_TOP}" && -f "${OCT_CRD}" ]] || { echo "Missing ${OCT_TOP} or ${OCT_CRD}" >&2; exit 1; }

echo "=== tleap: rectangular 9 Å ==="
run_tleap "${TLEAP_RECT}" "${RECT_LOG}"
[[ -f "${RECT_TOP}" && -f "${RECT_CRD}" ]] || { echo "Missing ${RECT_TOP} or ${RECT_CRD}" >&2; exit 1; }

echo "=== ambpdb: full-system PDBs ==="
pdb_from_topcrd "${OCT_TOP}" "${OCT_CRD}" "${BOX_DIR}/oct_9A_full.pdb"
pdb_from_topcrd "${RECT_TOP}" "${RECT_CRD}" "${BOX_DIR}/rect_9A_full.pdb"
echo "Wrote ${BOX_DIR}/oct_9A_full.pdb"
echo "Wrote ${BOX_DIR}/rect_9A_full.pdb"

oct_nat=$(natom_from_prmtop "${OCT_TOP}")
rect_nat=$(natom_from_prmtop "${RECT_TOP}")

oct_wat=$(count_wat_from_tleap_log "${OCT_LOG}")
rect_wat=$(count_wat_from_tleap_log "${RECT_LOG}")

oct_vol=$(grep -E "Volume:.*A\^3" "${OCT_LOG}" | tail -1 | sed 's/^[[:space:]]*//' || true)
rect_vol=$(grep -E "Volume:.*A\^3" "${RECT_LOG}" | tail -1 | sed 's/^[[:space:]]*//' || true)

oct_box=$(box_summary "${OCT_TOP}" "${OCT_CRD}")
rect_box=$(box_summary "${RECT_TOP}" "${RECT_CRD}")

diff=$((rect_nat - oct_nat))
pct="NA"
if [[ "${oct_nat}" =~ ^[0-9]+$ && "${rect_nat}" =~ ^[0-9]+$ && "${oct_nat}" -gt 0 ]]; then
  pct=$(python3 - <<PY
rect=${rect_nat}
oct=${oct_nat}
print(f"{(rect-oct)/oct*100.0:.4f}")
PY
)
fi

{
  echo "9 Å buffer box comparison (mol6 urease monomer + LIG, ff19SB + GAFF2, TIP3P)"
  echo "Date: $(date -Is)"
  echo
  echo "Inputs:"
  echo "  Octahedral:  ${TLEAP_OCT}"
  echo "  Rectangular: ${TLEAP_RECT}"
  echo
  echo "=== Truncated octahedron (oct_9A) ==="
  echo "  Files: ${OCT_TOP}"
  echo "         ${OCT_CRD}"
  echo "  PDB:   ${BOX_DIR}/oct_9A_full.pdb"
  echo "  Log:   ${OCT_LOG}"
  echo "  NATOM (prmtop POINTERS): ${oct_nat}"
  echo "  WAT residues (tleap log): ${oct_wat}"
  echo "  tleap volume line: ${oct_vol}"
  echo "  Box (parmed): ${oct_box}"
  echo
  echo "=== Rectangular (rect_9A) ==="
  echo "  Files: ${RECT_TOP}"
  echo "         ${RECT_CRD}"
  echo "  PDB:   ${BOX_DIR}/rect_9A_full.pdb"
  echo "  Log:   ${RECT_LOG}"
  echo "  NATOM (prmtop POINTERS): ${rect_nat}"
  echo "  WAT residues (tleap log): ${rect_wat}"
  echo "  tleap volume line: ${rect_vol}"
  echo "  Box (parmed): ${rect_box}"
  echo
  echo "=== Difference (rect - oct) ==="
  echo "  NATOM difference: ${diff}"
  echo "  Water difference: $((rect_wat - oct_wat))"
  echo "  Percent atom count difference (rect vs oct): (rect-oct)/oct*100 = ${pct}%"
  echo "    (positive => rectangular has more atoms; negative => oct has more)"
} > "${REPORT}"

echo "Wrote ${REPORT}"
sed -n '1,80p' "${REPORT}"
