#!/usr/bin/env bash
# Queue all systems with solvated_hmr.prmtop serially on GPU 1 (mol13 skipped — use run_mol13_now.sh)
set -euo pipefail
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${_SCRIPT_DIR}/.."
cd "${ROOT}"

while IFS= read -r SYS; do
  [[ -z "${SYS}" ]] && continue
  if [[ "${SYS}" == "mol13_pose1_monomer" ]]; then
    echo "Skipping mol13_pose1_monomer (run scripts/run_mol13_now.sh separately)"
    continue
  fi

  echo "=============================="
  echo "Starting: ${SYS}  ($(date))"
  echo "=============================="

  set +e
  bash "${_SCRIPT_DIR}/run_system_full.sh" "${SYS}"
  STATUS=$?
  set -e

  if [[ "${STATUS}" -ne 0 ]]; then
    echo "[QUEUE ERROR] ${SYS} failed with exit ${STATUS} — continuing to next system"
  else
    echo "[QUEUE OK] ${SYS} complete"
  fi
done < <(find systems -mindepth 2 -maxdepth 2 -name 'solvated_hmr.prmtop' -printf '%h\n' | sed 's|^systems/||' | sort -u)

echo "=============================="
echo "Queue finished: $(date)"
echo "=============================="
