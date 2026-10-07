#!/usr/bin/env bash
# One-off test: mol6_pose2_monomer from gpu_bridge through short equil_free GPU test.
# Not used by the main pipeline; run manually per project instructions.
set -euo pipefail
cd /path/to/urease_md
source config.env
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES=1
export PATH="/usr/bin:${AMBER_HOME}/bin:${PATH}"
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/openmpi/lib:${LD_LIBRARY_PATH:-}"

SIM=simulations/mol6_pose2_monomer
SYS=systems/mol6_pose2_monomer
PRMTOP=${SYS}/monomer_mol6.prmtop
INP=simulations/input_files

echo "=========================================="
echo "  Testing pipeline: mol6_pose2_monomer"
echo "=========================================="

mkdir -p "${SIM}/gpu_bridge"
"${PMEMD_CUDA_BIN:-/opt/software/AMBER/pmemd24_cuda128/bin/pmemd.cuda}" -O \
  -i "${INP}/gpu_bridge.in" \
  -o "${SIM}/gpu_bridge/gpu_bridge.out" \
  -p "${PRMTOP}" \
  -c "${SIM}/equil_cpu_npt/equil_cpu_npt.rst7" \
  -r "${SIM}/gpu_bridge/gpu_bridge.rst7" \
  -inf "${SIM}/gpu_bridge/gpu_bridge.mdinfo" 2>&1
echo "gpu_bridge exit=$?"

mkdir -p "${SIM}/equil1"
mpirun --oversubscribe -np 50 "${AMBER_HOME}/bin/pmemd.MPI" -O \
  -i "${INP}/equil1_cpu.in" \
  -o "${SIM}/equil1/equil1.out" \
  -p "${PRMTOP}" \
  -c "${SIM}/gpu_bridge/gpu_bridge.rst7" \
  -r "${SIM}/equil1/equil1.rst7" \
  -x "${SIM}/equil1/equil1.nc" \
  -ref "${SIM}/gpu_bridge/gpu_bridge.rst7" \
  -inf "${SIM}/equil1/equil1.mdinfo" 2>&1
echo "equil1 exit=$?"
grep "Density" "${SIM}/equil1/equil1.out" | tail -2

mkdir -p "${SIM}/equil2"
mpirun --oversubscribe -np 50 "${AMBER_HOME}/bin/pmemd.MPI" -O \
  -i "${INP}/equil2_cpu.in" \
  -o "${SIM}/equil2/equil2.out" \
  -p "${PRMTOP}" \
  -c "${SIM}/equil1/equil1.rst7" \
  -r "${SIM}/equil2/equil2.rst7" \
  -x "${SIM}/equil2/equil2.nc" \
  -ref "${SIM}/equil1/equil1.rst7" \
  -inf "${SIM}/equil2/equil2.mdinfo" 2>&1
echo "equil2 exit=$?"
grep "Density" "${SIM}/equil2/equil2.out" | tail -2

mkdir -p "${SIM}/equil3"
mpirun --oversubscribe -np 50 "${AMBER_HOME}/bin/pmemd.MPI" -O \
  -i "${INP}/equil3_cpu.in" \
  -o "${SIM}/equil3/equil3.out" \
  -p "${PRMTOP}" \
  -c "${SIM}/equil2/equil2.rst7" \
  -r "${SIM}/equil3/equil3.rst7" \
  -x "${SIM}/equil3/equil3.nc" \
  -ref "${SIM}/equil2/equil2.rst7" \
  -inf "${SIM}/equil3/equil3.mdinfo" 2>&1
echo "equil3 exit=$?"
grep "Density" "${SIM}/equil3/equil3.out" | tail -2

mkdir -p "${SIM}/gpu_bridge2"
"${PMEMD_CUDA_BIN:-/opt/software/AMBER/pmemd24_cuda128/bin/pmemd.cuda}" -O \
  -i "${INP}/gpu_bridge.in" \
  -o "${SIM}/gpu_bridge2/gpu_bridge2.out" \
  -p "${PRMTOP}" \
  -c "${SIM}/equil3/equil3.rst7" \
  -r "${SIM}/gpu_bridge2/gpu_bridge2.rst7" \
  -inf "${SIM}/gpu_bridge2/gpu_bridge2.mdinfo" 2>&1
echo "gpu_bridge2 exit=$?"

cat > /tmp/equil_free_test.in << 'MDIN'
equil_free GPU test - 2 ps only
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=1000, dt=0.002,
  ntt=11, temp0=300.0, tautp=2.0,
  ntb=2, ntp=1, barostat=2,
  pres0=1.0,
  cut=9.0,
  ntc=2, ntf=2,
  ntpr=500, ntwx=0, ntwr=1000,
  ntxo=2, ioutfm=1,
  ig=-1,
  ntr=1,
  restraintmask=':NI,KCX,OHB&!@H=',
  restraint_wt=1.0,
  watnam='WAT', owtnm='O',
 /
MDIN

"${PMEMD_CUDA_BIN:-/opt/software/AMBER/pmemd24_cuda128/bin/pmemd.cuda}" -O \
  -i /tmp/equil_free_test.in \
  -o /tmp/equil_free_test.out \
  -p "${PRMTOP}" \
  -c "${SIM}/gpu_bridge2/gpu_bridge2.rst7" \
  -r /tmp/equil_free_test.rst7 \
  -ref "${SIM}/gpu_bridge2/gpu_bridge2.rst7" \
  -inf /tmp/equil_free_test.mdinfo 2>&1
echo "equil_free_test exit=$?"
grep "Density" /tmp/equil_free_test.out | tail -2
grep "ns/day" /tmp/equil_free_test.out | tail -1

echo ""
echo "=========================================="
echo "  PIPELINE TEST SUMMARY"
echo "=========================================="
echo "gpu_bridge:      $([ -f ${SIM}/gpu_bridge/gpu_bridge.rst7 ] && echo 'OK' || echo 'FAILED')"
echo "equil1:          $([ -f ${SIM}/equil1/equil1.rst7 ] && echo 'OK' || echo 'FAILED')"
echo "equil2:          $([ -f ${SIM}/equil2/equil2.rst7 ] && echo 'OK' || echo 'FAILED')"
echo "equil3:          $([ -f ${SIM}/equil3/equil3.rst7 ] && echo 'OK' || echo 'FAILED')"
echo "gpu_bridge2:     $([ -f ${SIM}/gpu_bridge2/gpu_bridge2.rst7 ] && echo 'OK' || echo 'FAILED')"
echo "equil_free_test: $([ -f /tmp/equil_free_test.rst7 ] && echo 'OK' || echo 'FAILED')"
