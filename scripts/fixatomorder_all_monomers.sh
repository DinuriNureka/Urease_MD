#!/usr/bin/env bash
# Batch fixatomorder for all systems/*_monomer/ (Ni/OHB contiguous molecule fix).
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1
# shellcheck source=/dev/null
source config.env
export PATH="/usr/bin:${AMBER_HOME}/bin:${PATH}"
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/openmpi/lib:${LD_LIBRARY_PATH}"

LOGFILE="logs/fixatomorder_all36_$(date +%Y%m%d_%H%M%S).log"
mkdir -p logs

CPPTRAJ="${CPPTRAJ:-/usr/local/bin/cpptraj}"
[[ -x "${CPPTRAJ}" ]] || CPPTRAJ="$(command -v cpptraj)"

echo "=== fixatomorder rebuild for all monomer systems ===" | tee "${LOGFILE}"
echo "Started: $(date)" | tee -a "${LOGFILE}"
echo "cpptraj: ${CPPTRAJ}" | tee -a "${LOGFILE}"

PASS=0
FAIL=0
SKIP=0

ni_dist_awk() {
  # Standard data file: header line then frame data
  awk 'NF>=2 && $1 ~ /^[0-9]+$/ {print $2; exit}' "$1" 2>/dev/null || true
}

for SYS_DIR in systems/*_monomer/; do
  NAME=$(basename "${SYS_DIR}")

  PRMTOP=$(find "${SYS_DIR}" -maxdepth 1 -name 'monomer_*.prmtop' ! -name '*BAD*' ! -name '*backup*' ! -name '*.pre_fixatomorder' 2>/dev/null | sort | head -1)
  INPCRD=$(find "${SYS_DIR}" -maxdepth 1 -name 'monomer_*.inpcrd' ! -name '*BAD*' ! -name '*backup*' ! -name '*.pre_fixatomorder' 2>/dev/null | sort | head -1)

  if [[ -z "${PRMTOP}" || -z "${INPCRD}" ]]; then
    echo "SKIP ${NAME}: no prmtop/inpcrd found" | tee -a "${LOGFILE}"
    SKIP=$((SKIP + 1))
    continue
  fi

  echo -n "${NAME}: " | tee -a "${LOGFILE}"

  # Before fixatomorder, Ni/OHB sit in a tail segment; tLeap residue IDs are often 808/809 (not 570/571).
  PRE_DAT=$(mktemp)
  "${CPPTRAJ}" -p "${PRMTOP}" >>"${LOGFILE}" 2>&1 <<EOF || true
trajin ${INPCRD}
distance ni_ni :808@NI :809@NI
run
writedata ${PRE_DAT} ni_ni
EOF
  PRE_DIST=$(ni_dist_awk "${PRE_DAT}")
  rm -f "${PRE_DAT}"

  PREFIX="/tmp/fixatom_${NAME}_$$"
  if ! "${CPPTRAJ}" >>"${LOGFILE}" 2>&1 <<EOF
parm ${PRMTOP}
trajin ${INPCRD}
fixatomorder outprefix ${PREFIX}
trajout ${PREFIX}.rst7 restart
run
EOF
  then
    echo "FAIL — cpptraj fixatomorder error" | tee -a "${LOGFILE}"
    FAIL=$((FAIL + 1))
    rm -f "${PREFIX}"* 2>/dev/null || true
    continue
  fi

  FIXED_PRMTOP_ACTUAL="${PREFIX}.$(basename "${PRMTOP}")"
  FIXED_RST7="${PREFIX}.rst7"

  if [[ ! -f "${FIXED_PRMTOP_ACTUAL}" || ! -f "${FIXED_RST7}" ]]; then
    echo "FAIL — missing output (${FIXED_PRMTOP_ACTUAL})" | tee -a "${LOGFILE}"
    FAIL=$((FAIL + 1))
    rm -f "${PREFIX}"* 2>/dev/null || true
    continue
  fi

  POST_DAT=$(mktemp)
  "${CPPTRAJ}" -p "${FIXED_PRMTOP_ACTUAL}" >>"${LOGFILE}" 2>&1 <<EOF || true
trajin ${FIXED_RST7}
distance ni_ni :570@NI :571@NI
run
writedata ${POST_DAT} ni_ni
EOF
  POST_DIST=$(ni_dist_awk "${POST_DAT}")
  rm -f "${POST_DAT}"

  OK=$(python3 -c "
d = float('${POST_DIST:-nan}') if '${POST_DIST:-}' else float('nan')
print('YES' if (d == d and 2.5 < d < 5.0) else 'NO')
" 2>/dev/null || echo "NO")

  if [[ "${OK}" == "YES" ]]; then
    cp -a "${PRMTOP}" "${PRMTOP}.pre_fixatomorder"
    cp -a "${INPCRD}" "${INPCRD}.pre_fixatomorder"
    mv -f "${FIXED_PRMTOP_ACTUAL}" "${PRMTOP}"
    mv -f "${FIXED_RST7}" "${INPCRD}"
    echo "OK  Ni-Ni: ${PRE_DIST:-?} -> ${POST_DIST} A" | tee -a "${LOGFILE}"
    PASS=$((PASS + 1))
  else
    echo "FAIL  Ni-Ni post-fix: ${POST_DIST:-empty} (expected 2.5-5.0 A)" | tee -a "${LOGFILE}"
    FAIL=$((FAIL + 1))
  fi

  rm -f "${PREFIX}"* 2>/dev/null || true
done

echo "" | tee -a "${LOGFILE}"
echo "==============================" | tee -a "${LOGFILE}"
echo "PASS: ${PASS}  FAIL: ${FAIL}  SKIP: ${SKIP}" | tee -a "${LOGFILE}"
echo "Finished: $(date)" | tee -a "${LOGFILE}"
