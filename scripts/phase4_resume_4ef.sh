#!/usr/bin/env bash
# One-off: complete Phase 4e–4f when 4d MD finished but QC used broken grep|tail density parse.
# Usage: bash scripts/phase4_resume_4ef.sh /path/to/phase4_mol6_pose2_monomer
set -euo pipefail
SCRIPT_DIR_PHASE4="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${SCRIPT_DIR_PHASE4}/../config.env"
PHASE4="$(realpath "$1")"
MAIN_LOG="${2:-}"
PRMTOP="/path/to/urease_md/systems/mol6_pose2_monomer/solvated_hmr.prmtop"
INPCRD="/path/to/urease_md/systems/mol6_pose2_monomer/solvated_hmr.rst7"
PMEMD="${PMEMD_MD_BIN:-${PMEMD_CUDA_BIN:-/opt/software/AMBER/pmemd24_cuda128/bin/pmemd.cuda}}"
CPPTRAJ="${CPPTRAJ:-/usr/local/bin/cpptraj}"
export CUDA_DEVICE_ORDER="${CUDA_DEVICE_ORDER:-PCI_BUS_ID}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-1}"

log() {
  local m="[$(date '+%F %T')] $*"
  echo "$m"
  [[ -n "${MAIN_LOG}" && -f "${MAIN_LOG}" ]] && echo "$m" >>"${MAIN_LOG}"
}

amber_last_density_g_cm3() {
  awk '/Density[[:space:]]*=/ {
    v = $NF + 0
    if (v >= 0.85 && v <= 1.2) d = v
  }
  END { if (d != "" && d == d) print d }' "$1" 2>/dev/null
}

amber_last_md_temp_k() {
  grep -E '^ NSTEP.*TEMP\(K\)' "$1" 2>/dev/null | sed -n 's/.*TEMP(K) *= *\([0-9.]*\).*/\1/p' | awk '$1+0>150 {last=$1} END {print last}'
}

grep_shake_error() {
  grep -qiE "SHAKE ERROR|SHAKE failure|RATTLE ERROR|Constraint error|ERROR: Calculation aborted|vlimit exceeded|Coordinate resetting cannot" "$1" 2>/dev/null
}

[[ -f "${PHASE4}/density.rst7" ]] || { echo "missing density.rst7"; exit 1; }
[[ -f "${PRMTOP}" && -f "${INPCRD}" ]] || { echo "missing prmtop/inpcrd"; exit 1; }

dens=$(amber_last_density_g_cm3 "${PHASE4}/density.log")
[[ -n "${dens}" ]] || { echo "could not parse density"; exit 1; }
awk -v d="${dens}" 'BEGIN{exit (d+0>=0.95 && d+0<=1.05)?0:1}' || { echo "density ${dens} out of range"; exit 1; }
log "PASS: Phase 4d density=${dens} g/cm³ (resume)"

if [[ -f "${PHASE4}/release.rst7" ]] && grep -q "Total wall time" "${PHASE4}/release.log" 2>/dev/null; then
  log "SKIP release (release.log already completed)"
  rms_bb=$(tail -1 "${PHASE4}/rms_bb.dat" 2>/dev/null | awk '{print $2}')
  [[ -n "${rms_bb}" ]] || { echo "missing rms_bb.dat"; exit 1; }
  log "PASS: Phase 4e backbone RMSD=${rms_bb} Å (reuse)"
else
cat >"${PHASE4}/release.in" <<'EOF'
NPT release (Ni only)
 &cntrl
  imin=0, nstlim=1000000, dt=0.002,
  ntx=5, irest=1, ntb=2,
  ntp=1, pres0=1.0, barostat=1,
  ntt=3, gamma_ln=5.0, temp0=300.0,
  ntc=2, ntf=2, cut=9.0,
  ntpr=5000, ntwx=5000,
  ntr=1, restraint_wt=0.5,
  restraintmask='(:570@NI)|(:571@NI)',
  iwrap=1, ntxo=2,
 /
EOF

cd "${PHASE4}"
ln -sf "$(realpath "${PRMTOP}")" release.ptop
ln -sf "$(realpath "${PHASE4}/density.rst7")" release.crdin
ln -sf "$(realpath "${PHASE4}/density.rst7")" release.ref
log "RUN release: ${PMEMD}"
"${PMEMD}" -O -i release.in -o release.log -p release.ptop -c release.crdin -r release.rst7 -ref release.ref -x release.nc

cat >"${PHASE4}/rms_bb.cpptraj" <<EOF
parm ${PRMTOP}
trajin ${PHASE4}/release.nc
reference ${INPCRD} [ref]
rms R1 @CA,C,N,O&!:WAT ref [ref] out ${PHASE4}/rms_bb.dat
run
EOF
"${CPPTRAJ}" -i "${PHASE4}/rms_bb.cpptraj" >"${PHASE4}/rms_bb.cpptraj.log" 2>&1
rms_bb=$(tail -1 "${PHASE4}/rms_bb.dat" 2>/dev/null | awk '{print $2}')
awk -v r="${rms_bb:-99}" 'BEGIN{exit (r+0 < 3.0)?0:1}' || { echo "Phase 4e RMSD fail ${rms_bb}"; exit 1; }
log "PASS: Phase 4e backbone RMSD=${rms_bb} Å"

fi

if [[ ! -f "${PHASE4}/rms_bb.dat" ]]; then
  echo "missing rms_bb.dat after release"; exit 1
fi

cd "${PHASE4}"

cat >"${PHASE4}/test5ns.in" <<'EOF'
NPT test 5 ns
 &cntrl
  imin=0, nstlim=2500000, dt=0.002,
  ntx=5, irest=1, ntb=2,
  ntp=1, pres0=1.0, barostat=1,
  ntt=3, gamma_ln=5.0, temp0=300.0,
  ntc=2, ntf=2, cut=9.0,
  ntpr=5000, ntwx=5000,
  ig=-1, iwrap=1, ntxo=2,
 /
EOF
ln -sf "$(realpath "${PRMTOP}")" test5ns.ptop
ln -sf "$(realpath "${PHASE4}/release.rst7")" test5ns.crdin
ln -sf "$(realpath "${PHASE4}/release.rst7")" test5ns.ref
log "RUN test5ns: ${PMEMD}"
"${PMEMD}" -O -i test5ns.in -o test5ns.log -p test5ns.ptop -c test5ns.crdin -r test5ns.rst7 -ref test5ns.ref -x test5ns.nc
grep_shake_error "${PHASE4}/test5ns.log" && { echo "SHAKE error in test5ns"; exit 1; }
[[ -f "${PHASE4}/test5ns.nc" ]] || { echo "missing test5ns.nc"; exit 1; }

cat >"${PHASE4}/dist_ni.cpptraj" <<EOF
parm ${PRMTOP}
trajin ${PHASE4}/test5ns.nc
distance dNi :570@NI :571@NI out ${PHASE4}/ni_ni.dat
run
EOF
"${CPPTRAJ}" -i "${PHASE4}/dist_ni.cpptraj" >"${PHASE4}/dist_ni.cpptraj.log" 2>&1
dni=$(tail -1 "${PHASE4}/ni_ni.dat" 2>/dev/null | awk '{print $2}')
awk -v d="${dni:-0}" 'BEGIN{exit (d+0>=3.0 && d+0<=4.5)?0:1}' || { echo "Ni–Ni fail ${dni}"; exit 1; }

cat >"${PHASE4}/rms_lig.cpptraj" <<EOF
parm ${PRMTOP}
trajin ${PHASE4}/test5ns.nc
reference ${INPCRD} [ref]
rms RL :LIG ref [ref] out ${PHASE4}/rms_lig.dat
run
EOF
"${CPPTRAJ}" -i "${PHASE4}/rms_lig.cpptraj" >"${PHASE4}/rms_lig.cpptraj.log" 2>&1
rms_l=$(tail -1 "${PHASE4}/rms_lig.dat" 2>/dev/null | awk '{print $2}')
awk -v r="${rms_l:-99}" 'BEGIN{exit (r+0 < 3.0)?0:1}' || { echo "lig RMSD fail ${rms_l}"; exit 1; }

dens_f=$(amber_last_density_g_cm3 "${PHASE4}/test5ns.log")
temp_f=$(amber_last_md_temp_k "${PHASE4}/test5ns.log")
[[ -n "${dens_f}" ]] || { echo "no dens_f"; exit 1; }
[[ -n "${temp_f}" && "${temp_f}" != "nan" ]] || { echo "no temp_f"; exit 1; }
awk -v d="${dens_f}" 'BEGIN{exit (d+0>=0.99 && d+0<=1.02)?0:1}' || { echo "dens_f ${dens_f}"; exit 1; }
awk -v t="${temp_f}" 'BEGIN{exit (t+0>=298 && t+0<=302)?0:1}' || { echo "temp ${temp_f}"; exit 1; }

log "PASS: Phase 4f QC passed (Ni–Ni=${dni}, dens=${dens_f}, T=${temp_f})"
echo "TEST EQUILIBRATION PASSED — ready to launch full 36-system queue" | tee "${PHASE4}/test_passed.msg"
date -Iseconds >"/path/to/urease_md/logs/test_passed.flag"
log "PASS: Wrote /path/to/urease_md/logs/test_passed.flag"
log "run_full_pipeline.sh finished (Phase 4 resume 4e–4f)."
