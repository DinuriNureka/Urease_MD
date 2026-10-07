#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# build_all_no_ohb_systems.sh
#
# Build all 36 monomer (chains B+Q) solvated AMBER systems using the *no-OHB*
# Ni-center parameter set.
#
# Input list:
#   params/inhibitors/orca_opt/ligand_name_mapping.csv
#
# Output per system:
#   systems/<compound>_pose<pose>_monomer/solvated_hmr.prmtop
#   systems/<compound>_pose<pose>_monomer/solvated_hmr.rst7
#
# Logs:
#   logs/build_no_ohb_<sys>.log
#   logs/build_no_ohb_summary.txt
#
# IMPORTANT:
# - This script does NOT rename any existing files.
# - It skips systems that already have solvated_hmr.prmtop > 10 MB.
# - It will NOT build systems whose selected ORCA pose has imaginary frequencies.
# -----------------------------------------------------------------------------
set -euo pipefail

PROJECT_ROOT="/path/to/urease_md"
cd "${PROJECT_ROOT}"

# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config.env"

require_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "[FAIL] command not found: $1" >&2; exit 2; }; }
require_file() { [[ -f "$1" ]] || { echo "[FAIL] missing file: $1" >&2; exit 2; }; }

log_ts() { date '+%F %T'; }

CPPTRAJ_BIN="${CPPTRAJ:-/usr/local/bin/cpptraj}"
if [[ ! -x "${CPPTRAJ_BIN}" ]] && command -v cpptraj >/dev/null 2>&1; then
  CPPTRAJ_BIN="$(command -v cpptraj)"
fi

TLEAP_BIN="${TLEAP_BIN:-$(command -v tleap)}"

ORCA_MAP_CSV="${PROJECT_ROOT}/params/inhibitors/orca_opt/ligand_name_mapping.csv"
MONOMER_APO_PDB="${MONOMER_APO_PDB:-${PROJECT_ROOT}/structures/monomer/monomer_BQ_apo_no_ohb.pdb}"

# no-OHB Ni-center assets (CRITICAL)
KCX_MOL2="${PROJECT_ROOT}/params/ni_center_no_ohb/mcpb_input/KCX_full_stopgap.mol2"
KCX_FRCMOD="${PROJECT_ROOT}/params/ni_center_no_ohb/mcpb_input/KCX.frcmod"
UREASE_FRCMOD="${PROJECT_ROOT}/params/ni_center_no_ohb/mcpb_input/urease_no_ohb_mcpbpy.frcmod"
MON_BONDS="${PROJECT_ROOT}/params/ni_center_no_ohb/monomer/bonds_monomer_no_ohb.txt"
MON_PATCH="${PROJECT_ROOT}/params/ni_center_no_ohb/monomer/set_patch_monomer_no_ohb.txt"

# Same naming as existing working systems: octahedron with 12 Å buffer
WATER_BOX="TIP3PBOX"
BOX_BUFFER="9.0"

LOG_DIR="${PROJECT_ROOT}/logs"
mkdir -p "${LOG_DIR}"
SUMMARY="${LOG_DIR}/build_no_ohb_summary.txt"

require_cmd python3
require_cmd "${TLEAP_BIN}"
require_cmd "${CPPTRAJ_BIN}"
require_file "${ORCA_MAP_CSV}"
require_file "${MONOMER_APO_PDB}"
require_file "${KCX_MOL2}"
require_file "${KCX_FRCMOD}"
require_file "${UREASE_FRCMOD}"
require_file "${MON_BONDS}"
require_file "${MON_PATCH}"

PROVISIONAL_AUTOPATCH_FRCMOD="${PROVISIONAL_AUTOPATCH_FRCMOD:-${PROJECT_ROOT}/systems/mol6_pose2_full/provisional_autopatch_optB.frcmod}"
require_file "${PROVISIONAL_AUTOPATCH_FRCMOD}"

stat_size_bytes() { stat -c%s "$1" 2>/dev/null || echo 0; }

natoms_from_prmtop() {
  python3 - "$1" <<'PY'
import sys
try:
    import parmed as pmd
except Exception:
    print("0")
    raise SystemExit(0)
p = pmd.load_file(sys.argv[1])
print(len(p.atoms))
PY
}

orca_pose_is_clean_minimum() {
  # Return success (0) only if:
  # - ORCA log for mol_id shows normal termination
  # - ORCA log contains no explicit "***imaginary mode***" marker
  # - selected pose xyz exists either in workdir OR in params/inhibitors/structures/
  #
  # Note: workdir is the authoritative mapping for mol_id → ORCA directory (e.g. aha uses mol01/).
  local mol_id="$1" compound="$2" pose="$3" workdir="$4" selected_xyz="$5"
  local wd="${PROJECT_ROOT}/${workdir}"
  local log=""
  local xyz_wd="${wd}/${compound}_pose${pose}_lig.xyz"
  local xyz_struct="${PROJECT_ROOT}/params/inhibitors/structures/${selected_xyz}"
  local xyz_opt1="${wd}/${mol_id}_pose${pose}.xyz"
  local xyz_opt2="${wd}/${mol_id}_pose${pose}_28core.xyz"

  # ORCA log: try canonical, then pose-specific fallbacks (used by thiourea pose5, etc.)
  if [[ -f "${wd}/${mol_id}.log" ]]; then
    log="${wd}/${mol_id}.log"
  elif [[ -f "${wd}/${mol_id}_pose${pose}.log" ]]; then
    log="${wd}/${mol_id}_pose${pose}.log"
  elif [[ -f "${wd}/${mol_id}_pose${pose}_28core.log" ]]; then
    log="${wd}/${mol_id}_pose${pose}_28core.log"
  elif [[ -f "${wd}/${compound}_pose${pose}.log" ]]; then
    log="${wd}/${compound}_pose${pose}.log"
  fi

  [[ -n "${log}" && -f "${log}" ]] || return 1
  grep -q "ORCA TERMINATED NORMALLY" "${log}" || return 1
  if grep -qi "\\*\\*\\*imaginary mode\\*\\*\\*" "${log}"; then
    return 1
  fi
  [[ -f "${xyz_wd}" || -f "${xyz_struct}" || -f "${xyz_opt1}" || -f "${xyz_opt2}" ]] || return 1
  return 0
}

nini_distance_A() {
  # Print Ni–Ni distance (Å) or empty string.
  local prmtop="$1" rst7="$2"
  local tmp
  tmp="$(mktemp)"
  # Preferred per prompt (may fail if mask selects multiple atoms in some topologies)
  "${CPPTRAJ_BIN}" -p "${prmtop}" >/dev/null 2>&1 <<EOF || true
trajin ${rst7}
distance ni_ni :NI@NI :NI@NI out ${tmp}
run
EOF
  local d
  d="$(awk 'NF>=2 && $1 ~ /^[0-9]+$/ {print $2; exit}' "${tmp}" 2>/dev/null || true)"
  if [[ -z "${d}" ]]; then
    "${CPPTRAJ_BIN}" -p "${prmtop}" >/dev/null 2>&1 <<EOF || true
trajin ${rst7}
distance ni_ni :570@NI :571@NI out ${tmp}
run
EOF
    d="$(awk 'NF>=2 && $1 ~ /^[0-9]+$/ {print $2; exit}' "${tmp}" 2>/dev/null || true)"
  fi
  rm -f "${tmp}"
  printf '%s' "${d}"
}

has_no_ohb_residue() {
  local prmtop="$1"
  local out
  out="$("${CPPTRAJ_BIN}" 2>/dev/null <<EOF || true
parm ${prmtop}
resinfo :OHB
quit
EOF
)"
  # Pass if cpptraj indicates 0 residues selected.
  if echo "${out}" | grep -qiE "0 residues selected|No residues selected"; then
    return 0
  fi
  # Conservative fallback: fail if OHB appears as a selected residue.
  if echo "${out}" | grep -q "OHB"; then
    return 1
  fi
  return 0
}

resolve_inhibitor_params() {
  # Echo "mol2_path frcmod_path" or return non-zero.
  local compound="$1"
  local mol2="${PROJECT_ROOT}/params/inhibitors/mol2/${compound}_gaff2.mol2"
  local frc1="${PROJECT_ROOT}/params/inhibitors/frcmod/${compound}.frcmod"
  [[ -f "${mol2}" ]] || return 1
  [[ -f "${frc1}" ]] || return 1
  printf '%s %s\n' "${mol2}" "${frc1}"
}

build_one() {
  local mol_id="$1" compound="$2" pose="$3" selected_xyz="$4" net_charge="$5" mult="$6" workdir="$7"
  local sys="${compound}_pose${pose}_monomer"
  local sdir="${PROJECT_ROOT}/systems/${sys}"
  local logf="${LOG_DIR}/build_no_ohb_${sys}.log"

  mkdir -p "${sdir}"

  # Preflight checks outside per-system log redirection so failures are unambiguous
  # (and so callers piping stdout/stderr can see progress).
  local log="${PROJECT_ROOT}/${workdir}/${mol_id}.log"
  local xyz_wd="${PROJECT_ROOT}/${workdir}/${compound}_pose${pose}_lig.xyz"
  local xyz_struct="${PROJECT_ROOT}/params/inhibitors/structures/${selected_xyz}"

  if ! orca_pose_is_clean_minimum "${mol_id}" "${compound}" "${pose}" "${workdir}" "${selected_xyz}"; then
    {
      echo "[$(log_ts)] build_no_ohb: ${sys}"
      echo "  mol_id=${mol_id} compound=${compound} pose=${pose} net_charge=${net_charge} mult=${mult}"
      echo "  workdir=${workdir}"
      echo "  solvation=solvatebox ${WATER_BOX} ${BOX_BUFFER}"
      echo ""
      echo "[FAIL] ORCA pose not accepted for ${compound} (mol_id=${mol_id}, pose=${pose}, workdir=${workdir}, xyz=${selected_xyz})"
      echo "       log_exists=$([[ -f \"${log}\" ]] && echo yes || echo no) log=${log}"
      if [[ -f "${log}" ]]; then
        echo "       terminated_normally=$(grep -q \"ORCA TERMINATED NORMALLY\" \"${log}\" && echo yes || echo no)"
        echo "       imaginary_mode_lines=$(grep -ci \"\\*\\*\\*imaginary mode\\*\\*\\*\" \"${log}\" 2>/dev/null || echo 0)"
      fi
      echo "       xyz_in_workdir=$([[ -f \"${xyz_wd}\" ]] && echo yes || echo no) path=${xyz_wd}"
      echo "       xyz_in_structures=$([[ -f \"${xyz_struct}\" ]] && echo yes || echo no) path=${xyz_struct}"
    } > "${logf}" 2>&1
    echo -e "${sys}\tFAIL\t-\torca_pose_not_clean" >> "${SUMMARY}"
    return 1
  fi

  local imol ifrc
  if ! read -r imol ifrc < <(resolve_inhibitor_params "${compound}"); then
    {
      echo "[$(log_ts)] build_no_ohb: ${sys}"
      echo "  mol_id=${mol_id} compound=${compound} pose=${pose} net_charge=${net_charge} mult=${mult}"
      echo "  workdir=${workdir}"
      echo "  solvation=solvatebox ${WATER_BOX} ${BOX_BUFFER}"
      echo ""
      echo "[FAIL] Missing inhibitor params for ${compound}:"
      echo "       ${PROJECT_ROOT}/params/inhibitors/mol2/${compound}_gaff2.mol2"
      echo "       ${PROJECT_ROOT}/params/inhibitors/frcmod/${compound}.frcmod"
    } > "${logf}" 2>&1
    echo -e "${sys}\tFAIL\t-\tmissing_ligand_params" >> "${SUMMARY}"
    return 1
  fi

  echo "[$(log_ts)] (${sys}) Running tleap" >&2

  {
    echo "[$(log_ts)] build_no_ohb: ${sys}"
    echo "  mol_id=${mol_id} compound=${compound} pose=${pose} net_charge=${net_charge} mult=${mult}"
    echo "  workdir=${workdir}"
    echo "  solvation=solvatebox ${WATER_BOX} ${BOX_BUFFER}"
    echo ""

    local tleap_in="${sdir}/tleap_no_ohb_${sys}.in"
    local tleap_log="${sdir}/tleap_no_ohb_${sys}.log"

    cat > "${tleap_in}" <<EOF
# Generated by build_all_no_ohb_systems.sh — ${sys}
source ${PROTEIN_FF}
source ${LIGAND_FF}
source ${WATER_FF}

addAtomTypes {
        { "M1"  "Ni" "sp3" }
        { "M2"  "Ni" "sp3" }
        { "Y1"  "O"  "sp3" }
        { "Y2"  "N"  "sp3" }
        { "Y3"  "N"  "sp3" }
        { "Y4"  "O"  "sp3" }
        { "Y5"  "N"  "sp3" }
        { "Y6"  "N"  "sp3" }
        { "Y7"  "O"  "sp3" }
        { "Y8"  "O"  "sp3" }
}

# no-OHB Ni center parameters (NO OHB LOADS)
KCX = loadmol2 ${KCX_MOL2}
loadamberparams ${KCX_FRCMOD}
loadamberparams ${UREASE_FRCMOD}

# Provisional protein-only parameters (validated reference build)
loadamberparams ${PROVISIONAL_AUTOPATCH_FRCMOD}

# Inhibitor parameters (docked coordinates in mol2; RESP charges in mol2)
loadamberparams ${ifrc}
LIG = loadmol2 ${imol}

protein = loadpdb ${MONOMER_APO_PDB}
mol = combine { protein LIG }

# Apply no-OHB monomer metal bonds + patch charges/types
source ${MON_BONDS}
source ${MON_PATCH}

check mol
solvatebox mol ${WATER_BOX} ${BOX_BUFFER}

# Neutralize only (match existing systems)
addionsrand mol Na+ 0
addionsrand mol Cl- 0

saveamberparm mol ${sdir}/monomer_${compound}.prmtop ${sdir}/monomer_${compound}.inpcrd
quit
EOF

    echo "[$(log_ts)] Running tleap"
    if ! "${TLEAP_BIN}" -f "${tleap_in}" > "${tleap_log}" 2>&1; then
      echo "[FAIL] tleap non-zero exit (see ${tleap_log})"
      echo -e "${sys}\tFAIL\t-\ttleap_failed" >> "${SUMMARY}"
      return 1
    fi
    if grep -qiE "FATAL|Could not|could not" "${tleap_log}"; then
      echo "[FAIL] tleap fatal content (see ${tleap_log})"
      echo -e "${sys}\tFAIL\t-\ttleap_fatal" >> "${SUMMARY}"
      return 1
    fi

    local ptop="${sdir}/monomer_${compound}.prmtop"
    local pcrd="${sdir}/monomer_${compound}.inpcrd"
    [[ -f "${ptop}" && -f "${pcrd}" ]] || {
      echo "[FAIL] Missing tleap outputs: ${ptop} / ${pcrd}"
      echo -e "${sys}\tFAIL\t-\tmissing_tleap_outputs" >> "${SUMMARY}"
      return 1
    }

    echo "[$(log_ts)] Running cpptraj fixatomorder"
    CPPTRAJ="${CPPTRAJ_BIN}" bash "${PROJECT_ROOT}/scripts/fix_atomorder_amber.sh" "${ptop}" "${pcrd}" || {
      echo "[FAIL] fixatomorder failed"
      echo -e "${sys}\tFAIL\t-\tfixatomorder_failed" >> "${SUMMARY}"
      return 1
    }

    echo "[$(log_ts)] Applying HMR (parmed)"
    python3 - <<PY
import parmed as pmd
p = pmd.load_file("${ptop}", "${pcrd}")
pmd.tools.HMassRepartition(p).execute()
p.save("${sdir}/solvated_hmr.prmtop", overwrite=True)
p.save("${sdir}/solvated_hmr.rst7", overwrite=True)
PY

    local natoms
    natoms="$(natoms_from_prmtop "${sdir}/solvated_hmr.prmtop")"
    echo "[OK] Built ${sys} natoms=${natoms}"
    echo -e "${sys}\tPASS\t${natoms}\t-" >> "${SUMMARY}"
    return 0
  } > "${logf}" 2>&1
}

verify_one() {
  local sys="$1"
  local sdir="${PROJECT_ROOT}/systems/${sys}"
  local top="${sdir}/solvated_hmr.prmtop"
  local rst="${sdir}/solvated_hmr.rst7"

  local ok=1
  [[ -f "${top}" && $(stat_size_bytes "${top}") -gt $((10 * 1024 * 1024)) ]] || ok=0
  [[ -f "${rst}" && $(stat_size_bytes "${rst}") -gt $((1 * 1024 * 1024)) ]] || ok=0

  local d
  d="$(nini_distance_A "${top}" "${rst}")"
  if [[ -z "${d}" ]]; then ok=0; fi
  if [[ -n "${d}" ]]; then
    python3 - <<PY || ok=0
d = float("${d}")
raise SystemExit(0 if (3.0 <= d <= 4.5) else 1)
PY
  fi

  has_no_ohb_residue "${top}" || ok=0

  if (( ok )); then
    echo "PASS ${sys} Ni-Ni=${d} Å"
    return 0
  fi
  echo "FAIL ${sys} Ni-Ni=${d:-NA} Å"
  return 1
}

echo -e "system\tstatus\tnatoms\tnote" > "${SUMMARY}"

PASS=0
FAIL=0
SKIP=0

while IFS=, read -r mol_id _orig compound sel_xyz pose net_charge multiplicity workdir; do
  [[ "${mol_id}" == "mol_id" ]] && continue

  mol_id="$(echo "${mol_id}" | tr -d '\r')"
  compound="$(echo "${compound}" | tr -d '\r')"
  sel_xyz="$(echo "${sel_xyz}" | tr -d '\r')"
  pose="$(echo "${pose}" | tr -d '\r')"
  net_charge="$(echo "${net_charge}" | tr -d '\r')"
  multiplicity="$(echo "${multiplicity}" | tr -d '\r')"
  workdir="$(echo "${workdir}" | tr -d '\r')"

  sys="${compound}_pose${pose}_monomer"
  if build_one "${mol_id}" "${compound}" "${pose}" "${sel_xyz}" "${net_charge}" "${multiplicity}" "${workdir}"; then
    PASS=$((PASS + 1))
  else
    FAIL=$((FAIL + 1))
  fi
done < "${ORCA_MAP_CSV}"

echo "" >> "${SUMMARY}"
echo -e "TOTAL\tPASS=${PASS}\tFAIL=${FAIL}\tSKIP=${SKIP}" >> "${SUMMARY}"

echo "[INFO] Build finished. Summary: ${SUMMARY}"
echo "[INFO] Starting verification loop (36 systems)..."

VPASS=0
VFAIL=0
while IFS=, read -r _mol_id _orig compound _sel_xyz pose _net _mult _workdir; do
  [[ "${_mol_id}" == "mol_id" ]] && continue
  compound="$(echo "${compound}" | tr -d '\r')"
  pose="$(echo "${pose}" | tr -d '\r')"
  sys="${compound}_pose${pose}_monomer"
  if verify_one "${sys}"; then
    VPASS=$((VPASS + 1))
  else
    VFAIL=$((VFAIL + 1))
  fi
done < "${ORCA_MAP_CSV}"

echo "[INFO] Verification summary: PASS=${VPASS} FAIL=${VFAIL}"
exit $(( VFAIL == 0 ? 0 : 1 ))

