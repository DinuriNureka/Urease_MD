#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# verify_system_ready.sh
#
# Standalone safety verification for a built system directory.
#
# Usage:
#   bash scripts/verify_system_ready.sh <system_dir>
#
# Exit codes:
#   0 = all PASS
#   1 = any FAIL
# -----------------------------------------------------------------------------
set -uo pipefail

SDIR="${1:?Usage: bash scripts/verify_system_ready.sh <system_dir>}"
PROJECT_ROOT="/path/to/urease_md"
MAP_CSV="${PROJECT_ROOT}/params/inhibitors/orca_opt/ligand_name_mapping.csv"

# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config.env"
set +e  # Allow check commands to return non-zero without exiting

CPPTRAJ="${AMBER_HOME}/bin/cpptraj"
[[ -x "${CPPTRAJ}" ]] || CPPTRAJ="$(command -v cpptraj || true)"
CPPTRAJ_OK=1
if [[ ! -x "${CPPTRAJ}" ]]; then
  CPPTRAJ_OK=0
fi

PRMTOP="${SDIR}/solvated_hmr.prmtop"
RST7="${SDIR}/solvated_hmr.rst7"

FAILS=0
passfail() { local ok="$1"; shift; printf '%s\n' "$([[ "${ok}" -eq 1 ]] && echo "PASS" || echo "FAIL")" "$*"; }

echo "=== SYSTEM: $(basename "${SDIR}") ==="

# CHECK 1 — No OHB in topology
OHB_COUNT="NA"
if [[ "${CPPTRAJ_OK}" -eq 1 && -s "${PRMTOP}" ]]; then
  _ohb_out="$(printf 'parm %s\nresinfo :OHB\nquit\n' "${PRMTOP}" | "${CPPTRAJ}" 2>/dev/null)"
  if echo "${_ohb_out}" | grep -q " OHB "; then
    OHB_COUNT=1
  else
    OHB_COUNT=0
  fi
fi
if [[ "${OHB_COUNT}" == "0" ]]; then
  echo "CHECK1 [PASS] OHB residues in topology: ${OHB_COUNT}"
else
  echo "CHECK1 [FAIL] OHB residues in topology: ${OHB_COUNT}"
  FAILS=$((FAILS + 1))
fi

# CHECK 2 — Ni-Ni distance 3.0–4.5 Å
NI_NI="NA"
NI1=""
NI2=""
if [[ -s "${PRMTOP}" && -s "${RST7}" ]]; then
  NI_RES="$(python3 -c "
import parmed as pmd, sys
p = pmd.load_file('${SDIR}/solvated_hmr.prmtop')
ni = [str(i+1) for i,r in enumerate(p.residues) if r.name=='NI']
print(ni[0], ni[1]) if len(ni)>=2 else sys.exit(1)
" 2>/dev/null)"
  NI1="$(echo "${NI_RES}" | cut -d' ' -f1)"
  NI2="$(echo "${NI_RES}" | cut -d' ' -f2)"

  if [[ -n "${NI1}" && -n "${NI2}" ]]; then
    TMPDAT="$(mktemp /tmp/nini_XXXXXX.dat)"
    printf 'parm %s\ntrajin %s\ndistance NiNi :%s@NI :%s@NI out %s\nrun\nquit\n' \
      "${SDIR}/solvated_hmr.prmtop" \
      "${SDIR}/solvated_hmr.rst7" \
      "${NI1}" "${NI2}" "${TMPDAT}" | \
      "${CPPTRAJ}" >/dev/null 2>/dev/null
    NI_NI="$(awk 'NR==2{print $2}' "${TMPDAT}" 2>/dev/null || true)"
    rm -f "${TMPDAT}"
  fi
fi

OK="$(awk -v d="${NI_NI:-0}" 'BEGIN{print (d+0>=3.0 && d+0<=4.5)?"yes":"no"}')"
if [[ "${OK}" == "yes" ]]; then
  echo "CHECK2 [PASS] Ni-Ni distance: ${NI_NI} Å"
else
  echo "CHECK2 [FAIL] Ni-Ni distance: ${NI_NI} Å"
  FAILS=$((FAILS+1))
fi

# CHECK 3 — No-OHB frcmod used in tleap log
TLEAP_LOG="${SDIR}/tleap_no_ohb_$(basename ${SDIR}).log"
CHECK3_OK=0
if [[ -f "${TLEAP_LOG}" ]]; then
  grep -q "urease_no_ohb_mcpbpy.frcmod" "${TLEAP_LOG}"
  [[ $? -eq 0 ]] && CHECK3_OK=1
else
  # Also try any tleap_pipeline_*.log as fallback
  TLEAP_LOG="$(ls "${SDIR}"/tleap_pipeline_*.log 2>/dev/null | tail -1)"
  if [[ -n "${TLEAP_LOG}" ]]; then
    grep -q "urease_no_ohb_mcpbpy.frcmod" "${TLEAP_LOG}"
    [[ $? -eq 0 ]] && CHECK3_OK=1
  fi
fi
if [[ "${CHECK3_OK}" -eq 1 ]]; then
  echo "CHECK3 [PASS] no-OHB frcmod in tleap log"
else
  echo "CHECK3 [FAIL] no-OHB frcmod in tleap log"
  FAILS=$((FAILS + 1))
fi

# CHECK 4 — No dodecamer reference in system dir
HITS="$(
  grep -rl "dodecamer\\|12-mer\\|12mer\\|full_dodecamer" "${SDIR}" \
    --include="*.in" --include="*.log" --include="*.sh" 2>/dev/null || true
)"
if [[ -z "${HITS}" ]]; then
  echo "CHECK4 [PASS] No dodecamer references"
else
  echo "CHECK4 [FAIL] No dodecamer references"
  echo "${HITS}" | sed 's/^/  - /'
  FAILS=$((FAILS + 1))
fi

# CHECK 5 — Ligand ORCA: zero imaginary frequencies, normal termination
SYS_NAME="$(basename "${SDIR}")"
COMPOUND="$(echo "${SYS_NAME}" | sed -E 's/_pose[0-9]+_monomer$//')"
WORKDIR="$(awk -F',' -v c="${COMPOUND}" 'NR>1 && $3==c {print $8; exit}' "${MAP_CSV}" 2>/dev/null | tr -d '\r')"

ORCA_OK=0
if [[ -n "${WORKDIR}" && -d "${PROJECT_ROOT}/${WORKDIR}" ]]; then
  LOG_CAND=""
  # Prefer molXX.log based on mapping (mol_id column)
  MOL_ID="$(awk -F',' -v c="${COMPOUND}" 'NR>1 && $3==c {print $1; exit}' "${MAP_CSV}" 2>/dev/null | tr -d '\r')"
  if [[ -n "${MOL_ID}" && -f "${PROJECT_ROOT}/${WORKDIR}/${MOL_ID}.log" ]]; then
    LOG_CAND="${PROJECT_ROOT}/${WORKDIR}/${MOL_ID}.log"
  else
    # Fallback: first *.log in workdir
    LOG_CAND="$(ls "${PROJECT_ROOT}/${WORKDIR}"/*.log 2>/dev/null | head -1 || true)"
  fi

  if [[ -n "${LOG_CAND}" && -f "${LOG_CAND}" ]]; then
    grep -q "ORCA TERMINATED NORMALLY" "${LOG_CAND}"
    if [[ $? -eq 0 ]]; then
      IMAG_LINES="$(grep -ci "imaginary mode" "${LOG_CAND}" 2>/dev/null)"
      if [[ "${IMAG_LINES}" -eq 0 ]]; then
        ORCA_OK=1
      else
        grep -qE "Number of imaginary.*0|Total number of imaginary perturbations.*0" "${LOG_CAND}" 2>/dev/null
        [[ $? -eq 0 ]] && ORCA_OK=1
      fi
    fi
  fi
fi

if [[ "${ORCA_OK}" -eq 1 ]]; then
  echo "CHECK5 [PASS] ORCA terminated normally, 0 imaginary modes"
else
  echo "CHECK5 [FAIL] ORCA terminated normally, 0 imaginary modes"
  FAILS=$((FAILS + 1))
fi

if [[ "${FAILS}" -eq 0 ]]; then
  echo "OVERALL: READY"
  exit 0
fi
echo "OVERALL: NOT READY (${FAILS} checks failed)"
exit 1

