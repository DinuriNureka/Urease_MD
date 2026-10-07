#!/usr/bin/env bash
# Queue mol1, mol13, mol14, mol19, mol33: wait for mol6 equilibration to finish, then
# run equilibration + production sequentially on CUDA_VISIBLE_DEVICES=1.
set -euo pipefail
cd /path/to/urease_md
source config.env
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES=1
export PATH="/usr/bin:${AMBER_HOME}/bin:${PATH}"
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/openmpi/lib:${LD_LIBRARY_PATH:-}"

LOG="logs/queue_5_$(date +%Y%m%d_%H%M%S).log"

echo "=== Queue 5 priority compounds ===" | tee "${LOG}"
echo "Started: $(date)" | tee -a "${LOG}"
echo "Waiting for mol6_pose2_monomer equilibration to finish..." | tee -a "${LOG}"

# Wait until no process whose command line contains mol6_pose2_monomer (equilibration driver + engines)
while pgrep -f "mol6_pose2_monomer" >/dev/null 2>&1; do
  sleep 60
done
echo "mol6 done: $(date)" | tee -a "${LOG}"

for MOL in mol1 mol13 mol14 mol19 mol33; do
  SYS_DIR=$(ls -d systems/"${MOL}"_pose*_monomer 2>/dev/null | head -1)
  if [[ -z "${SYS_DIR}" ]]; then
    echo "SKIP ${MOL}: no system dir" | tee -a "${LOG}"
    continue
  fi
  NAME=$(basename "${SYS_DIR}")

  echo "" | tee -a "${LOG}"
  echo "============================================" | tee -a "${LOG}"
  echo "Starting ${MOL} (${NAME}): $(date)" | tee -a "${LOG}"
  echo "============================================" | tee -a "${LOG}"

  echo "Equilibration..." | tee -a "${LOG}"
  set +e
  bash scripts/run_equilibration.sh "${NAME}" >>"${LOG}" 2>&1
  EQUIL_EXIT=$?
  set -e
  echo "Equilibration exit: ${EQUIL_EXIT}" | tee -a "${LOG}"

  if [[ ${EQUIL_EXIT} -ne 0 ]]; then
    echo "FAILED ${MOL} equilibration — skipping production" | tee -a "${LOG}"
    continue
  fi

  echo "Production..." | tee -a "${LOG}"
  set +e
  bash scripts/run_production.sh "${NAME}" >>"${LOG}" 2>&1
  PROD_EXIT=$?
  set -e
  echo "Production exit: ${PROD_EXIT}" | tee -a "${LOG}"

  echo "Completed ${MOL}: $(date) equil=${EQUIL_EXIT} prod=${PROD_EXIT}" | tee -a "${LOG}"
done

echo "" | tee -a "${LOG}"
echo "=== Queue finished: $(date) ===" | tee -a "${LOG}"
