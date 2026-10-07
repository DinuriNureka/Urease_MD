#!/usr/bin/env bash
# One-line status per system (equil3, seg1–5, rep1–3)
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

printf "%-35s %s\n" "SYSTEM" "STATUS"
printf "%-35s %s\n" "------" "------"

while IFS= read -r SYS; do
  [[ -z "${SYS}" ]] && continue
  DIR="${ROOT}/simulations/${SYS}"

  if [[ -f "${DIR}/equil3/equil3.rst7" ]]; then E3="E3:✓"; else E3="E3:-"; fi

  if [[ -f "${DIR}/crashtest/patched_test/patched_test.rst7" ]] || [[ -f "${DIR}/equil_seg1/equil_seg1.rst7" ]]; then
    S1="S1:✓"
  else
    S1="S1:-"
  fi

  if [[ -f "${DIR}/equil_seg2/equil_seg2.rst7" ]]; then S2="S2:✓"; else S2="S2:-"; fi
  if [[ -f "${DIR}/equil_seg3/equil_seg3.rst7" ]]; then S3="S3:✓"; else S3="S3:-"; fi
  if [[ -f "${DIR}/equil_seg4/equil_seg4.rst7" ]]; then S4="S4:✓"; else S4="S4:-"; fi
  if [[ -f "${DIR}/equil_seg5/equil_seg5.rst7" ]]; then S5="S5:✓"; else S5="S5:-"; fi

  if [[ -f "${DIR}/production_rep1/production_rep1.rst7" ]]; then R1="R1:✓"; else R1="R1:-"; fi
  if [[ -f "${DIR}/production_rep2/production_rep2.rst7" ]]; then R2="R2:✓"; else R2="R2:-"; fi
  if [[ -f "${DIR}/production_rep3/production_rep3.rst7" ]]; then R3="R3:✓"; else R3="R3:-"; fi

  printf "%-35s %s %s %s %s %s %s %s %s %s\n" "${SYS}" "${E3}" "${S1}" "${S2}" "${S3}" "${S4}" "${S5}" "${R1}" "${R2}" "${R3}"
done < <(find "${ROOT}/systems" -mindepth 2 -maxdepth 2 -name 'solvated_hmr.prmtop' -printf '%h\n' | sed "s|^${ROOT}/systems/||" | sort -u)
