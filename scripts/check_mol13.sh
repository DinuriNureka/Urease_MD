#!/bin/bash
cd /path/to/urease_md
SIM=simulations/mol13_pose1_monomer

echo "=== mol13 equil_free @ $(date +%H:%M:%S) ==="

if [ -f ${SIM}/equil_free/equil_free.mdinfo ]; then
  # mdinfo is updated live during the run
  echo "--- mdinfo (live) ---"
  cat ${SIM}/equil_free/equil_free.mdinfo
fi

echo ""
if [ -f ${SIM}/equil_free/equil_free.out ]; then
  PROGRESS=$(grep -E "NSTEP\s*=" ${SIM}/equil_free/equil_free.out | tail -1 | grep -oE "NSTEP =\s*[0-9]+" | grep -oE "[0-9]+")
  NAN_COUNT=$(grep -c NaN ${SIM}/equil_free/equil_free.out)
  LAST_TEMP=$(grep "TEMP(K)" ${SIM}/equil_free/equil_free.out | tail -1 | grep -oE "TEMP\(K\) =\s*[0-9.]+" | grep -oE "[0-9.]+$")
  LAST_DENSITY=$(grep "Density" ${SIM}/equil_free/equil_free.out | tail -1 | grep -oE "[0-9]+\.[0-9]+")
  PCT=$(awk "BEGIN {printf \"%.2f\", $PROGRESS / 25000000 * 100}")
  REMAINING=$((25000000 - PROGRESS))
  echo "--- mdout state ---"
  echo "NSTEP:    $PROGRESS / 25,000,000  (${PCT}%, ${REMAINING} to go)"
  echo "TEMP:     $LAST_TEMP K"
  echo "Density:  $LAST_DENSITY g/cm³"
  echo "NaN:      $NAN_COUNT $([ "$NAN_COUNT" -eq 0 ] && echo '✅' || echo '❌ STOP NOW')"
fi

echo ""
pgrep -af "pmemd.cuda.*mol13" > /dev/null && echo "Process: ALIVE ✅" || echo "Process: DEAD ❌"

echo ""
echo "--- GPU ---"
nvidia-smi --query-gpu=index,utilization.gpu,memory.used,temperature.gpu,power.draw,fan.speed --format=csv,noheader | grep -E "P100"
