#!/usr/bin/env bash
# Run equilibration + production for five inhibitors given by CANONICAL chemistry names.
#
# Naming (see params/inhibitors/COMPOUND_NAMING.txt):
#   - "mol1" … "mol33" = numbered synthesis series (NOT sequential ORCA folders mol01–mol36).
#   - Standards at end of full screen: aha, dmj, thiourea (not in this batch).
# Full-project order when listing compounds: mol1 → mol2 → … → mol33, then aha, dmj, thiourea.
#
# This batch: canonical compounds in numeric order (1, 13, 14, 19, 33) — not file/ORCA order.
# ORCA batch folders for reference (ligand_name_mapping.csv mol_id → compound):
#   mol1  → params/.../orca_opt/mol13
#   mol13 → params/.../orca_opt/mol06
#   mol14 → params/.../orca_opt/mol07
#   mol19 → params/.../orca_opt/mol12
#   mol33 → params/.../orca_opt/mol28
#
set -euo pipefail
cd /path/to/urease_md
source config.env
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES=1
export PATH="/usr/bin:${AMBER_HOME}/bin:${PATH}"
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/openmpi/lib:${LD_LIBRARY_PATH:-}"

# Canonical chemistry IDs, sorted by compound number (not mol01–mol36 batch order).
CANONICAL_COMPOUNDS=(mol1 mol13 mol14 mol19 mol33)

LOG="batch_5_compounds.log"
{
  echo "Batch run started: $(date)"
  echo "Compounds (canonical chemistry order, mol1–mol33 series): ${CANONICAL_COMPOUNDS[*]}"
  echo "These are NOT AHA/DMJ/thiourea and NOT ordered as orca_opt/mol01..mol36."
} > "${LOG}"

for MOL in "${CANONICAL_COMPOUNDS[@]}"; do
  echo "" >> "${LOG}"
  echo "============================================" >> "${LOG}"
  echo "Starting canonical ${MOL}: $(date)" >> "${LOG}"
  echo "============================================" >> "${LOG}"

  # System directory: systems/<compound>_pose<k>_monomer (same tag as in ligand_name_mapping selected_pose)
  SYS_DIR=$(ls -d systems/${MOL}_pose*_monomer 2>/dev/null | sort -V | head -1)
  if [[ -z "${SYS_DIR}" ]]; then
    SYS_DIR=$(ls -d "systems/${MOL}_monomer" 2>/dev/null | head -1)
  fi

  if [[ -z "${SYS_DIR}" ]]; then
    echo "SKIP ${MOL}: no system directory found (expected systems/${MOL}_pose*_monomer)" >> "${LOG}"
    continue
  fi

  NAME=$(basename "${SYS_DIR}")
  echo "  → simulation/system basename: ${NAME}" >> "${LOG}"

  echo "Running equilibration for ${NAME} (canonical ${MOL})..." >> "${LOG}"
  bash scripts/run_equilibration.sh "${NAME}" >> "${LOG}" 2>&1
  EQUIL_EXIT=$?
  echo "Equilibration exit: ${EQUIL_EXIT}" >> "${LOG}"

  if [[ ${EQUIL_EXIT} -ne 0 ]]; then
    echo "FAILED equilibration for ${NAME}, skipping production" >> "${LOG}"
    continue
  fi

  echo "Running production for ${NAME}..." >> "${LOG}"
  bash scripts/run_production.sh "${NAME}" >> "${LOG}" 2>&1
  PROD_EXIT=$?
  echo "Production exit: ${PROD_EXIT}" >> "${LOG}"

  echo "Completed canonical ${MOL} (${NAME}): $(date), equil=${EQUIL_EXIT}, prod=${PROD_EXIT}" >> "${LOG}"
done

echo "" >> "${LOG}"
echo "Batch run finished: $(date)" >> "${LOG}"
