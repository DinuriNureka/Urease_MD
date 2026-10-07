#!/usr/bin/env bash
# Generate ESP surface cube files for all compounds
# For visual QC of RESP charge fitting
# Usage: bash scripts/batch_esp_check.sh [single_mol_id]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/config.env"

cd "${ROOT}"

ORCA_DIR="${ROOT}/params/inhibitors/orca_opt"
ORCA_2MKL_BIN="${ORCA_2MKL:-${ORCA_PATH%/*}/orca_2mkl}"
MULTIWFN_BIN="${MULTIWFN_PATH}"
SINGLE="${1:-}"

if [[ -n "${SINGLE}" ]]; then
    MOLS=("${SINGLE}")
else
    mapfile -t MOLS < <(find "${ORCA_DIR}" -maxdepth 1 -type d -name 'mol[0-9][0-9]' -printf '%f\n' | sort -V)
fi

echo "ESP Cube Generation for RESP Surface Visualization"
echo "==================================================="
echo "Compounds to process: ${#MOLS[@]}"
echo ""

FAILED=()
SKIPPED=()
SUCCESS=()

for MOL in "${MOLS[@]}"; do
    DIR="${ORCA_DIR}/${MOL}"
    GBW="${DIR}/${MOL}.gbw"

    echo "--- ${MOL} ---"

    if [[ ! -d "${DIR}" ]]; then
        echo "  SKIP: directory missing"
        SKIPPED+=("${MOL}")
        continue
    fi

    if [[ ! -f "${GBW}" ]]; then
        echo "  SKIP: No .gbw file"
        SKIPPED+=("${MOL}")
        continue
    fi

    LOG="${DIR}/${MOL}.log"
    if [[ -f "${LOG}" ]]; then
        if ! grep -q "ORCA TERMINATED NORMALLY" "${LOG}"; then
            echo "  SKIP: ORCA did not terminate normally"
            SKIPPED+=("${MOL}")
            continue
        fi

        if grep -q '\*\*\*imaginary mode\*\*\*' "${LOG}"; then
            echo "  WARNING: Has imaginary vibrational mode(s)!"
            grep '\*\*\*imaginary mode\*\*\*' "${LOG}" || true
        fi
    fi

    pushd "${DIR}" >/dev/null

    MOLDEN=""
    if [[ -f "${MOL}.molden.input" ]]; then
        MOLDEN="${MOL}.molden.input"
    elif [[ -f "${MOL}.molden" ]]; then
        MOLDEN="${MOL}.molden"
    fi

    if [[ -z "${MOLDEN}" ]]; then
        echo "  Converting .gbw → .molden..."
        if ! "${ORCA_2MKL_BIN}" "${MOL}" -molden 2>/dev/null; then
            echo "  FAIL: orca_2mkl failed"
            FAILED+=("${MOL}")
            popd >/dev/null || true
            continue
        fi
        if [[ -f "${MOL}.molden.input" ]]; then
            MOLDEN="${MOL}.molden.input"
        elif [[ -f "${MOL}.molden" ]]; then
            MOLDEN="${MOL}.molden"
        else
            echo "  FAIL: orca_2mkl produced no .molden file"
            FAILED+=("${MOL}")
            popd >/dev/null || true
            continue
        fi
    fi

    echo "  Generating ESP cube via Multiwfn (${MOLDEN})..."
    # Multiwfn 2026.x: grid (5) → ESP (12) → medium grid (2) → export cube (2); writes totesp.cub
    if ! "${MULTIWFN_BIN}" "${MOLDEN}" <<'MWFN' >/dev/null 2>&1
5
12
2
2
0
q
MWFN
    then
        echo "  FAIL: Multiwfn ESP export failed"
        FAILED+=("${MOL}")
        popd >/dev/null || true
        continue
    fi

    CUBES=$(find . -maxdepth 1 -name '*.cub' 2>/dev/null | wc -l)
    if [[ "${CUBES}" -gt 0 ]]; then
        echo "  OK: ${CUBES} cube file(s) generated"
        SUCCESS+=("${MOL}")
    else
        echo "  FAIL: No cube files generated"
        FAILED+=("${MOL}")
    fi

    popd >/dev/null
done

echo ""
echo "==================================================="
echo "SUMMARY"
echo "  Success: ${#SUCCESS[@]} — ${SUCCESS[*]:-none}"
echo "  Skipped: ${#SKIPPED[@]} — ${SKIPPED[*]:-none}"
echo "  Failed:  ${#FAILED[@]} — ${FAILED[*]:-none}"
echo ""
echo "To view ESP surfaces:"
echo "  VMD:  vmd <mol>.cub → Isosurface, color by Volume"
echo "  Jmol: jmol <mol>.cub → isosurface vdw map <mol>_esp.cub"
echo ""
echo "Look for:"
echo "  - Smooth charge distribution (no spikes)"
echo "  - Negative regions (red) on electronegative atoms (O, N)"
echo "  - Positive regions (blue) on H, electropositive groups"
echo "  - Any unusual charge hotspots near the pharmacophore"
