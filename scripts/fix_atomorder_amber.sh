#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# fix_atomorder_amber.sh — reorder atoms after tLeap saveamberparm (Ni/OHB fix)
#
# tLeap can leave one "molecule" with non-contiguous atom ranges (protein tail
# + KCX/OHB). cpptraj fixatomorder makes molecules contiguous; pmemd.MPI and
# consistent analysis require matching coordinates, so we trajout a restart.
#
# Usage:
#   bash scripts/fix_atomorder_amber.sh <prmtop> <inpcrd_or_rst7>
#
# Replaces both files in place. Requires cpptraj (override with CPPTRAJ=...).
# -----------------------------------------------------------------------------
set -euo pipefail

PRMTOP="${1:?Usage: $0 <prmtop> <coords>}"
CRD="${2:?}"

CPPTRAJ_EXE="${CPPTRAJ:-/usr/local/bin/cpptraj}"
if [[ ! -x "${CPPTRAJ_EXE}" ]] && command -v cpptraj >/dev/null 2>&1; then
  CPPTRAJ_EXE="$(command -v cpptraj)"
fi
[[ -x "${CPPTRAJ_EXE}" ]] || {
  echo "fix_atomorder_amber: cpptraj not found (set CPPTRAJ=...)" >&2
  exit 1
}

TMPD="$(mktemp -d)"
cleanup() { rm -rf "${TMPD}"; }
trap cleanup EXIT

BASE="$(basename "${PRMTOP}")"
PREFIX="${TMPD}/fo"

"${CPPTRAJ_EXE}" <<EOF
parm ${PRMTOP}
trajin ${CRD}
fixatomorder outprefix ${PREFIX}
trajout ${TMPD}/fixed.rst7 restart
run
EOF

OUTTOP="${PREFIX}.${BASE}"
[[ -f "${OUTTOP}" ]] || {
  echo "fix_atomorder_amber: expected topology ${OUTTOP}" >&2
  exit 1
}
[[ -f "${TMPD}/fixed.rst7" ]] || {
  echo "fix_atomorder_amber: missing reordered coordinates" >&2
  exit 1
}

mv -f "${OUTTOP}" "${PRMTOP}"
mv -f "${TMPD}/fixed.rst7" "${CRD}"

exit 0
