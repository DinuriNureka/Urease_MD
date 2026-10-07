#!/usr/bin/env bash
# Ni–Ni distance on final equil_free restart vs crystal-derived tolerance.
# Usage: validate_ni_geometry.sh <system_name> <prmtop> <equil_free_dir>
set -euo pipefail

SCRIPT_ROOT="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=/dev/null
source "${SCRIPT_ROOT}/../config.env"

NAME="${1:?system name}"
PRMTOP="${2:?prmtop}"
EQUIL_DIR="${3:?equil_free directory}"

RST="${EQUIL_DIR}/equil_free.rst7"
[[ -f "${RST}" ]] || {
  echo "[GATE FAIL] Missing equil restart: ${RST}" >&2
  exit 1
}
[[ -f "${PRMTOP}" ]] || {
  echo "[GATE FAIL] Missing prmtop: ${PRMTOP}" >&2
  exit 1
}

NI1="${EQUIL_NI1:-570}"
NI2="${EQUIL_NI2:-571}"

CPP="$(command -v cpptraj || true)"
if [[ -z "${CPP}" && -n "${AMBERHOME:-}" ]]; then
  CPP="${AMBERHOME}/bin/cpptraj"
fi
if [[ -z "${CPP}" && -n "${AMBER_HOME:-}" ]]; then
  CPP="${AMBER_HOME}/bin/cpptraj"
fi
[[ -x "${CPP}" ]] || {
  echo "[GATE FAIL] cpptraj not found (PATH or AMBERHOME/AMBER_HOME)." >&2
  exit 1
}

INP="$(mktemp)"
DIST="$(mktemp)"
trap 'rm -f "${INP}" "${DIST}"' EXIT

cat >"${INP}" <<EOF
parm ${PRMTOP}
trajin ${RST}
autoimage
distance dNi :${NI1}@NI :${NI2}@NI out ${DIST}
run
EOF

if ! "${CPP}" -i "${INP}" >/dev/null 2>&1; then
  echo "[GATE FAIL] cpptraj failed for Ni–Ni geometry check (${NAME})." >&2
  exit 1
fi

d=""
while read -r line; do
  [[ "${line}" =~ ^# ]] && continue
  [[ -z "${line// }" ]] && continue
  d=$(echo "${line}" | awk '{print $2}')
done < "${DIST}"

if [[ -z "${d}" ]]; then
  echo "[GATE FAIL] No Ni–Ni distance parsed for ${NAME}." >&2
  exit 1
fi

if ! python3 -c "import sys; d=float(sys.argv[1]); sys.exit(0 if 3.21<=d<=3.81 else 1)" "${d}"; then
  echo "[GATE FAIL] Ni–Ni distance ${d} Å out of crystal range [3.21, 3.81] Å (ref 3.51 Å) for ${NAME}." >&2
  exit 1
fi

echo "[GATE PASS] Ni–Ni geometry OK (${d} Å) for ${NAME}."
exit 0
