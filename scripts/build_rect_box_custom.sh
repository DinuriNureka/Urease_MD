#!/usr/bin/env bash
set -euo pipefail
# Build mol6 monomer + LIG in a fixed orthorhombic box (tleap setBox centers + solvatebox 0).
# 1) Export dry PDB for extent. 2) Compute clearances. 3) Full solvated prmtop/inpcrd + PDB.
#
# Usage:
#   bash scripts/build_rect_box_custom.sh [A] [B] [C]
# Default A B C = 100 110 108 (see docs: 100 Å on Z is too small for this solute)
#
# Outputs under params/box_comparison/:
#   rect_<A>x<B>x<C>.{prmtop,inpcrd,_full.pdb}, tleap log, dry PDB, report snippet

PROJECT_ROOT="/path/to/urease_md"
cd "${PROJECT_ROOT}"
BOX_DIR="${PROJECT_ROOT}/params/box_comparison"
SCRIPT_ROOT="${PROJECT_ROOT}/scripts"

TA="${1:-100}"
TB="${2:-110}"
TC="${3:-108}"
TAG="${TA}x${TB}x${TC}"

DRY_PDB="${BOX_DIR}/mol6_dry_for_box_${TAG}.pdb"
DRY_TLEAP="${BOX_DIR}/tleap_dry_export_${TAG}.in"
FINAL_TLEAP="${BOX_DIR}/tleap_rect_${TAG}.in"
FINAL_LOG="${BOX_DIR}/tleap_rect_${TAG}.log"
OUT_TOP="${BOX_DIR}/rect_${TAG}.prmtop"
OUT_CRD="${BOX_DIR}/rect_${TAG}.inpcrd"
OUT_PDB="${BOX_DIR}/rect_${TAG}_full.pdb"
REPORT="${BOX_DIR}/box_rect_${TAG}_report.txt"

require_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing: $1" >&2; exit 2; }; }
require_cmd tleap
AMBPDB="$(command -v ambpdb || command -v bashambpdb || true)"
[[ -n "${AMBPDB}" ]] || { echo "Missing ambpdb" >&2; exit 2; }

sed "s/mol6_dry_for_box_100x110x100\\.pdb/mol6_dry_for_box_${TAG}.pdb/" \
  "${BOX_DIR}/tleap_rect_100x110x100_dry_export.in" > "${DRY_TLEAP}"

echo "=== Dry export (for clearance) ==="
tleap -f "${DRY_TLEAP}" > "${BOX_DIR}/tleap_dry_export_${TAG}.log" 2>&1
grep -qi "Fatal Error" "${BOX_DIR}/tleap_dry_export_${TAG}.log" && {
  echo "tleap dry export failed" >&2
  tail -20 "${BOX_DIR}/tleap_dry_export_${TAG}.log"
  exit 1
}

CLEAR="$("${SCRIPT_ROOT}/compute_clearance_for_box.py" "${DRY_PDB}" "${TA}" "${TB}" "${TC}" 2>/dev/null)" || {
  echo "[FAIL] Solute does not fit in ${TA} x ${TB} x ${TC} Å (axis-aligned bounding box)." >&2
  "${SCRIPT_ROOT}/compute_clearance_for_box.py" "${DRY_PDB}" "${TA}" "${TB}" "${TC}" || true
  exit 1
}

echo "setBox centers: ${CLEAR}"

cat > "${FINAL_TLEAP}" <<EOF
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.tip3p

addAtomTypes {
        { "M1"  "Ni" "sp3" }
        { "M2"  "Ni" "sp3" }
        { "Y1"  "O" "sp3" }
        { "Y2"  "N" "sp3" }
        { "Y3"  "N" "sp3" }
        { "Y4"  "O" "sp3" }
        { "Y5"  "N" "sp3" }
        { "Y6"  "N" "sp3" }
        { "Y7"  "O" "sp3" }
        { "Y8"  "O" "sp3" }
}

KCX = loadmol2 ${PROJECT_ROOT}/params/ni_center/mcpb_input/KCX_full_stopgap.mol2
OHB = loadmol2 ${PROJECT_ROOT}/params/ni_center/mcpb_input/OHB.mol2
loadamberparams ${PROJECT_ROOT}/params/ni_center/mcpb_input/KCX.frcmod
loadamberparams ${PROJECT_ROOT}/params/ni_center/mcpb_input/OHB.frcmod
loadamberparams ${PROJECT_ROOT}/params/ni_center/mcpb_input/urease_mcpbpy.frcmod
loadamberparams ${PROJECT_ROOT}/systems/mol6_pose2_full/provisional_autopatch_optB.frcmod

LIG = loadmol2 ${PROJECT_ROOT}/systems/mol6_pose2_full/mol6_pose2_gaff2.mol2
loadamberparams ${PROJECT_ROOT}/systems/mol6_pose2_full/mol6_pose2.frcmod

protein = loadpdb ${PROJECT_ROOT}/structures/monomer/monomer_BQ_apo.pdb
mol = combine { protein LIG }

check mol

source ${PROJECT_ROOT}/params/ni_center/monomer/bonds_monomer.txt
source ${PROJECT_ROOT}/params/ni_center/monomer/set_patch_monomer.txt

check mol

center mol
setBox mol centers ${CLEAR}

solvatebox mol TIP3PBOX 0.0

addionsrand mol Na+ 0
addionsrand mol Cl- 0

saveamberparm mol ${OUT_TOP} ${OUT_CRD}
quit
EOF

echo "=== Full solvated build ${TAG} Å ==="
tleap -f "${FINAL_TLEAP}" > "${FINAL_LOG}" 2>&1
grep -qi "Fatal Error" "${FINAL_LOG}" && {
  echo "tleap final build failed" >&2
  tail -30 "${FINAL_LOG}"
  exit 1
}

"${AMBPDB}" -p "${OUT_TOP}" -c "${OUT_CRD}" > "${OUT_PDB}"

BOX_LINE="$(tail -1 "${OUT_CRD}")"
{
  echo "Rectangular box build (mol6 monomer + LIG)"
  echo "Target edges (Å): ${TA} x ${TB} x ${TC}"
  echo "setBox centers: ${CLEAR}"
  echo "Last line of inpcrd (a b c alpha beta gamma):"
  echo "  ${BOX_LINE}"
  echo "Files:"
  echo "  ${OUT_TOP}"
  echo "  ${OUT_CRD}"
  echo "  ${OUT_PDB}"
  echo "Dry reference PDB: ${DRY_PDB}"
} > "${REPORT}"

echo "Wrote ${REPORT}"
cat "${REPORT}"
