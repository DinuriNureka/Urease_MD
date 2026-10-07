#!/usr/bin/env bash
set -euo pipefail

# Compare truncated octahedral (solvateOct) vs rectangular (solvateBox) water boxes
# for the mol6_pose2_monomer build recipe, without touching any simulation outputs.
#
# Outputs:
#   params/box_comparison/oct/solvated.{prmtop,inpcrd,full.pdb}
#   params/box_comparison/rect/solvated.{prmtop,inpcrd,full.pdb}
#   params/box_comparison/box_comparison_report.txt

PROJECT_ROOT="/path/to/urease_md"
cd "${PROJECT_ROOT}"

BOX_DIR="${PROJECT_ROOT}/params/box_comparison"
OCT_DIR="${BOX_DIR}/oct"
RECT_DIR="${BOX_DIR}/rect"
REPORT="${BOX_DIR}/box_comparison_report.txt"

TLEAP_OCT="${BOX_DIR}/tleap_octahedral.in"
TLEAP_RECT="${BOX_DIR}/tleap_rectangular.in"

mkdir -p "${OCT_DIR}" "${RECT_DIR}"

require_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing command: $1" >&2; exit 2; }; }
require_cmd tleap

# Prefer ambpdb but allow bashambpdb for compatibility with user instructions.
AMBPDB="$(command -v ambpdb || command -v bashambpdb || true)"
[[ -n "${AMBPDB}" ]] || { echo "Missing ambpdb/bashambpdb on PATH." >&2; exit 2; }

run_tleap() {
  local in="$1" outlog="$2"
  echo "Running: tleap -f ${in}"
  tleap -f "${in}" > "${outlog}" 2>&1
  if grep -qiE "FATAL|fatal|Error:" "${outlog}"; then
    echo "tleap reported fatal errors in ${outlog}" >&2
    exit 1
  fi
}

pdb_from_topcrd() {
  local ptop="$1" crd="$2" outpdb="$3"
  "${AMBPDB}" -p "${ptop}" -c "${crd}" > "${outpdb}"
}

count_atoms_in_pdb() {
  # Count ATOM+HETATM records. (Used as a quick cross-check only.)
  grep -Ec '^(ATOM|HETATM)' "$1"
}

count_wat_from_tleap_log() {
  # Prefer tleap's own water residue count (PDB residue IDs can wrap at 9999).
  local lg="$1"
  # Example line: "        WAT     36924"
  awk 'toupper($1)=="WAT" && $2 ~ /^[0-9]+$/ {n=$2} END{print n+0}' "$lg"
}

count_total_residues_added_from_tleap_log() {
  # Example: "  Added 36936 residues."
  local lg="$1"
  awk 'tolower($1)=="added" && $2 ~ /^[0-9]+$/ && tolower($3) ~ /^resid/ {n=$2} END{print n+0}' "$lg"
}

volume_line_from_tleap_log() {
  local lg="$1"
  grep -E "Volume:.*A\^3" "$lg" | tail -1 | sed 's/^[[:space:]]*//'
}

box_from_prmtop() {
  # Quick summary: NATOM + box type (if present) and last BOX_DIMENSIONS line if any.
  local ptop="$1"
  python3 - "$ptop" <<'PY'
import sys
from pathlib import Path
ptop = Path(sys.argv[1])
txt = ptop.read_text(errors="ignore").splitlines()
nat = None
for ln in txt:
    if ln.strip().startswith("%FLAG POINTERS"):
        break
for i, ln in enumerate(txt):
    if ln.strip().startswith("%FLAG POINTERS"):
        # Next %FORMAT line then data lines with integers; NATOM is first integer.
        j = i + 2
        ints = []
        while j < len(txt) and not txt[j].startswith("%FLAG"):
            ints.extend(txt[j].split())
            if len(ints) >= 1:
                break
            j += 1
        if ints:
            nat = ints[0]
        break
print(f"NATOM={nat or 'NA'}")
PY
}

echo "=== Building octahedral box ==="
run_tleap "${TLEAP_OCT}" "${OCT_DIR}/tleap.log"
[[ -f "${OCT_DIR}/solvated.prmtop" && -f "${OCT_DIR}/solvated.inpcrd" ]] || { echo "Oct build missing outputs." >&2; exit 1; }
pdb_from_topcrd "${OCT_DIR}/solvated.prmtop" "${OCT_DIR}/solvated.inpcrd" "${OCT_DIR}/solvated_full.pdb"

echo "=== Building rectangular box ==="
run_tleap "${TLEAP_RECT}" "${RECT_DIR}/tleap.log"
[[ -f "${RECT_DIR}/solvated.prmtop" && -f "${RECT_DIR}/solvated.inpcrd" ]] || { echo "Rect build missing outputs." >&2; exit 1; }
pdb_from_topcrd "${RECT_DIR}/solvated.prmtop" "${RECT_DIR}/solvated.inpcrd" "${RECT_DIR}/solvated_full.pdb"

oct_atoms=$(box_from_prmtop "${OCT_DIR}/solvated.prmtop" | awk -F= '{print $2}')
rect_atoms=$(box_from_prmtop "${RECT_DIR}/solvated.prmtop" | awk -F= '{print $2}')

oct_atoms_pdb=$(count_atoms_in_pdb "${OCT_DIR}/solvated_full.pdb")
rect_atoms_pdb=$(count_atoms_in_pdb "${RECT_DIR}/solvated_full.pdb")

oct_wat=$(count_wat_from_tleap_log "${OCT_DIR}/tleap.log")
rect_wat=$(count_wat_from_tleap_log "${RECT_DIR}/tleap.log")

oct_added=$(count_total_residues_added_from_tleap_log "${OCT_DIR}/tleap.log")
rect_added=$(count_total_residues_added_from_tleap_log "${RECT_DIR}/tleap.log")

diff_atoms=$((rect_atoms - oct_atoms))
diff_wat=$((rect_wat - oct_wat))

oct_atom_only=$(grep -c '^ATOM' "${OCT_DIR}/solvated_full.pdb" || true)
rect_atom_only=$(grep -c '^ATOM' "${RECT_DIR}/solvated_full.pdb" || true)
oct_hetatm=$(grep -c '^HETATM' "${OCT_DIR}/solvated_full.pdb" || true)
rect_hetatm=$(grep -c '^HETATM' "${RECT_DIR}/solvated_full.pdb" || true)

pct_atoms_rect_vs_oct="NA"
pct_wat_rect_vs_oct="NA"
if [[ "${oct_atoms}" -gt 0 ]]; then
  pct_atoms_rect_vs_oct=$(python3 - <<PY
rect=${rect_atoms}
oct=${oct_atoms}
print(f"{(rect-oct)/oct*100.0:.2f}")
PY
)
fi
if [[ "${oct_wat}" -gt 0 ]]; then
  pct_wat_rect_vs_oct=$(python3 - <<PY
rect=${rect_wat}
oct=${oct_wat}
print(f"{(rect-oct)/oct*100.0:.2f}")
PY
)
fi

{
  echo "Box comparison report (mol6 monomer build recipe)"
  echo "Date: $(date -Is)"
  echo
  echo "Inputs:"
  echo "  Octahedral:   ${TLEAP_OCT}"
  echo "  Rectangular:  ${TLEAP_RECT}"
  echo
  echo "Outputs:"
  echo "  Octahedral dir:  ${OCT_DIR}"
  echo "  Rectangular dir: ${RECT_DIR}"
  echo
  echo "Counts (from PDB; includes grep -c '^ATOM' as requested):"
  echo "  Octahedral:  grep -c '^ATOM' = ${oct_atom_only} ; grep -Ec '^(ATOM|HETATM)' = ${oct_atoms_pdb} (HETATM lines: ${oct_hetatm})"
  echo "  Rectangular: grep -c '^ATOM' = ${rect_atom_only} ; grep -Ec '^(ATOM|HETATM)' = ${rect_atoms_pdb} (HETATM lines: ${rect_hetatm})"
  echo
  echo "Counts (from prmtop NATOM):"
  echo "  Octahedral atoms:  ${oct_atoms}"
  echo "  Rectangular atoms: ${rect_atoms}"
  echo "  Difference (rect - oct): ${diff_atoms}"
  echo
  echo "Water molecules (from tleap log WAT residue count):"
  echo "  Octahedral waters:  ${oct_wat}"
  echo "  Rectangular waters: ${rect_wat}"
  echo "  Difference (rect - oct): ${diff_wat}"
  echo
  echo "Residues added by solvation (from tleap log):"
  echo "  Octahedral Added:  ${oct_added}"
  echo "  Rectangular Added: ${rect_added}"
  echo
  echo "tleap-reported solvation volumes:"
  echo "  Octahedral:  $(volume_line_from_tleap_log "${OCT_DIR}/tleap.log")"
  echo "  Rectangular: $(volume_line_from_tleap_log "${RECT_DIR}/tleap.log")"
  echo
  echo "Relative size (rect vs oct):"
  echo "  Atoms:  (rect-oct)/oct*100 = ${pct_atoms_rect_vs_oct}%  (negative means oct is larger)"
  echo "  Waters: (rect-oct)/oct*100 = ${pct_wat_rect_vs_oct}%  (negative means oct has more waters)"
  echo
  echo "Interpretation:"
  echo "  Positive (rect-oct) means rectangular is larger; negative means truncated oct is larger."
  echo "  For this mol6 monomer recipe at 12.0 Å buffer, compare the tleap volumes and NATOM above."
  echo
  echo "Sanity:"
  echo "  Octahedral:  NATOM=${oct_atoms}"
  echo "  Rectangular: NATOM=${rect_atoms}"
  echo
  echo "tleap logs:"
  echo "  Octahedral:  ${OCT_DIR}/tleap.log"
  echo "  Rectangular: ${RECT_DIR}/tleap.log"
} > "${REPORT}"

echo "Wrote report: ${REPORT}"
sed -n '1,40p' "${REPORT}"

