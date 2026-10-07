#!/bin/bash
# Sequential queue for 5 chosen systems
cd /path/to/urease_md
source config.env

SYSTEMS=(
  "mol13_pose1_monomer"
  "mol14_pose2_monomer"
  "mol19_pose2_monomer"
  "mol1_pose2_monomer"
  "mol33_pose3_monomer"
)

for SYS in "${SYSTEMS[@]}"; do
  echo "[$(date)] Starting ${SYS}"
  bash scripts/run_equilibration.sh ${SYS} \
    > logs/equil_${SYS}.log 2>&1
  echo "[$(date)] Finished ${SYS} exit: $?"
done
