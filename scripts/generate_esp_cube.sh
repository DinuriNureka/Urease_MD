#!/usr/bin/env bash
# Generate ESP cube from ORCA .gbw for visualization (Multiwfn).
# Usage: bash scripts/generate_esp_cube.sh <path_to_orca_output_dir>
# Requires: Multiwfn, orca_2mkl (see config.env ORCA_2MKL / ORCA_PATH)
set -euo pipefail
SCRIPT_ROOT="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=/dev/null
source "${SCRIPT_ROOT}/../config.env"

COMPOUND_DIR="${1:?Usage: $0 <path_to_orca_output_dir>}"
cd "${COMPOUND_DIR}"

GBW=$(find . -maxdepth 6 -name "*.gbw" -print -quit)
if [[ -z "${GBW}" ]]; then
  echo "ERROR: No .gbw file found under ${COMPOUND_DIR}" >&2
  exit 1
fi

GBW_DIR=$(dirname "${GBW}")
BASE=$(basename "${GBW%.gbw}")
cd "${GBW_DIR}"

ORCA2MKL="${ORCA_2MKL:-}"
[[ -x "${ORCA2MKL}" ]] || ORCA2MKL="$(command -v orca_2mkl 2>/dev/null || true)"
if [[ -z "${ORCA2MKL}" || ! -x "${ORCA2MKL}" ]]; then
  echo "ERROR: orca_2mkl not found (set ORCA_2MKL or install)." >&2
  exit 1
fi

echo "Processing: ${GBW}"
echo "Converting .gbw → molden..."
"${ORCA2MKL}" "${BASE}" -molden

MOLDEN="${BASE}.molden.input"
if [[ ! -f "${MOLDEN}" ]]; then
  echo "ERROR: Expected molden file not found: $(pwd)/${MOLDEN}" >&2
  exit 1
fi

MULTIWFN_CMD="${MULTIWFN_PATH:-$(command -v Multiwfn 2>/dev/null || true)}"
if [[ -n "${MULTIWFN_CMD}" && -d "${MULTIWFN_CMD}" ]]; then
  MULTIWFN_CMD="${MULTIWFN_CMD%/}/Multiwfn"
fi
if [[ -z "${MULTIWFN_CMD}" || ! -x "${MULTIWFN_CMD}" ]]; then
  echo "ERROR: Multiwfn not found (set MULTIWFN_PATH to the binary)." >&2
  exit 1
fi

echo "Generating ESP grid + Gaussian cube via Multiwfn (2026.x: main 5 → ESP 12 → medium grid 2 → export cube 2)..."
cat <<'MWFN' | "${MULTIWFN_CMD}" "${MOLDEN}"
5
12
2
2
0
q
MWFN

echo ""
echo "Done. Default output name: totesp.cub (see $(pwd)/)"
echo "View with:"
echo "  vmd totesp.cub"
echo "  jmol totesp.cub  # console: isosurface vdw map totesp.cub"
