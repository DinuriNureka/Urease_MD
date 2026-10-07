#!/usr/bin/env bash
# Regenerate systems/*/ni_ni.RST for every directory that has solvated_hmr.prmtop.
# Usage (from project root):
#   source config.env
#   bash scripts/gen_ni_ni_rst_all.sh
set -euo pipefail
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../config.env
source "${_SCRIPT_DIR}/../config.env"
export PATH="/usr/bin:${AMBER_HOME}/bin:/usr/local/bin:${PATH}"

parse_ni_pair() {
  local prmtop="$1" out
  out=$(printf 'parm %s\natominfo :570@NI\natominfo :571@NI\natominfo @NI\nquit\n' "$prmtop" | cpptraj 2>/dev/null) || true
  echo "$out" | awk '
    function flush_allni(   i, min1, min2, v) {
      if (nni == 2) { print allni[1], allni[2]; return }
      if (nni < 2) return
      min1 = allni[1] + 0
      for (i = 2; i <= nni; i++)
        if (allni[i] + 0 < min1) min1 = allni[i] + 0
      min2 = 0
      for (i = 1; i <= nni; i++) {
        v = allni[i] + 0
        if (v == min1) continue
        if (min2 == 0 || v < min2) min2 = v
      }
      if (min1 && min2) print min1, min2
    }
    /atominfo :570@NI/ { sec = 570; next }
    /atominfo :571@NI/ { sec = 571; next }
    /atominfo @NI/ { sec = 999; next }
    sec == 570 && $1 ~ /^[0-9]+$/ && $2 == "NI" { a570 = $1; sec = 0; next }
    sec == 571 && $1 ~ /^[0-9]+$/ && $2 == "NI" { a571 = $1; sec = 0; next }
    sec == 999 && $1 ~ /^[0-9]+$/ && $2 == "NI" { allni[++nni] = $1; next }
    END {
      if (a570 && a571) print a570, a571
      else flush_allni()
    }
  '
}

for d in "${PROJECT_ROOT}/systems"/*/; do
  prm="${d}solvated_hmr.prmtop"
  [[ -f "${prm}" ]] || continue
  name=$(basename "${d}")
  pair=$(parse_ni_pair "${prm}")
  if [[ -z "${pair}" ]] || [[ "${pair}" != *" "* ]]; then
    echo "[FAIL] ${name}: could not determine two NI atom indices" >&2
    continue
  fi
  read -r i1 i2 <<<"${pair}"
  cat >"${d}ni_ni.RST" <<EOF
&rst
  iat=${i1},${i2},
  r1=2.5, r2=3.3, r3=3.7, r4=4.5,
  rk2=50.0, rk3=50.0,
/
EOF
  echo "Wrote ${d}ni_ni.RST (iat=${i1},${i2})"
done
