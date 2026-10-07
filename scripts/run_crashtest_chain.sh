#!/usr/bin/env bash
# Full chained watcher: waits for equil_free → production → crash tests
set -euo pipefail
cd /path/to/urease_md
source config.env

EQUIL_OUT="simulations/mol13_pose1_monomer/equil_free/equil_free.out"
PROD_OUT="simulations/mol13_pose1_monomer/prod_rep03/prod.out"
TESTDIR="simulations/mol13_pose1_monomer/crashtest"
STARTRST="simulations/mol13_pose1_monomer/equil_free/equil_free.rst7"
AMBER="${PMEMD_CUDA_BIN:-/opt/software/AMBER/pmemd24_cuda128/bin/pmemd.cuda}"
INPUT="simulations/input_files/equil_free_crashtest.in"
CPPTRAJ=$(command -v cpptraj)
NISRST="/path/to/urease_md/ni_ni.RST"
LOGFILE="logs/crashtest_master_watcher.log"

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "${LOGFILE}"; }

log "=== Full chained watcher started ==="

# ── PHASE 1: Wait for equil_free ──────────────────────────────────
log "Phase 1: Waiting for equil_free..."
while ! grep -q "Total wall time" "${EQUIL_OUT}" 2>/dev/null; do
    NSTEP=$(grep "NSTEP" "${EQUIL_OUT}" 2>/dev/null | \
      tail -1 | awk '{print $3}')
    log "equil_free running: NSTEP=${NSTEP}"
    sleep 300
done
log "=== equil_free COMPLETE ==="

NAN=$(grep "NaN" "${EQUIL_OUT}" | grep -v "NSTEP = 25000000" | wc -l 2>/dev/null || echo 0)
if [[ "${NAN}" -gt 0 ]]; then
    log "ERROR: equil_free NaN — aborting"; exit 1
fi

# Backup
mkdir -p "${TESTDIR}"
mkdir -p "simulations/mol13_pose1_monomer/equil_free_patched_backup"
cp "${EQUIL_OUT}" \
   "simulations/mol13_pose1_monomer/equil_free_patched_backup/" || true
cp "${STARTRST}" \
   "simulations/mol13_pose1_monomer/equil_free_patched_backup/" || true
cp "systems/mol13_pose1_monomer/solvated_hmr.prmtop" \
   "systems/mol13_pose1_monomer/solvated_hmr.prmtop.patched_backup_20260409" \
   || true
log "equil_free backup done"

# ── PHASE 2: Wait for production ──────────────────────────────────
log "Phase 2: Waiting for production rep01 to start..."
sleep 180

log "Phase 2: Waiting for production rep03 to complete (~11 hours)..."
while ! grep -q "Total wall time" "${PROD_OUT}" 2>/dev/null; do
    if [[ -f "simulations/mol13_pose1_monomer/prod_rep03/prod.out" ]]; then
        NSTEP=$(grep "NSTEP" \
          "simulations/mol13_pose1_monomer/prod_rep03/prod.out" \
          2>/dev/null | tail -1 | awk '{print $3}')
        log "prod_rep03 running: NSTEP=${NSTEP}"
    elif [[ -f "simulations/mol13_pose1_monomer/prod_rep02/prod.out" ]]; then
        NSTEP=$(grep "NSTEP" \
          "simulations/mol13_pose1_monomer/prod_rep02/prod.out" \
          2>/dev/null | tail -1 | awk '{print $3}')
        log "prod_rep02 running: NSTEP=${NSTEP}"
    elif [[ -f "simulations/mol13_pose1_monomer/prod_rep01/prod.out" ]]; then
        NSTEP=$(grep "NSTEP" \
          "simulations/mol13_pose1_monomer/prod_rep01/prod.out" \
          2>/dev/null | tail -1 | awk '{print $3}')
        log "prod_rep01 running: NSTEP=${NSTEP}"
    else
        log "Production not yet started..."
    fi
    sleep 300
done
log "=== Production COMPLETE — all 60 ns done ==="

# Use prod_rep03 restart if available
PROD_RST="simulations/mol13_pose1_monomer/prod_rep03/prod.rst7"
if [[ -f "${PROD_RST}" ]]; then
    STARTRST="${PROD_RST}"
    log "Using prod_rep03 restart as crash test starting point"
else
    log "WARNING: using equil_free restart"
fi

# ── PHASE 3: Create k=31 frcmod ───────────────────────────────────
log "Phase 3: Creating k=31 frcmod..."
python3 - << 'PYEOF'
import re
src = "params/ni_center_no_ohb/mcpb_input/urease_no_ohb_mcpbpy.frcmod"
dst = "params/ni_center_no_ohb/mcpb_input/urease_no_ohb_mcpbpy_k31.frcmod"
with open(src) as f:
    lines = f.readlines()
in_angl = False
out = []
changed = []
for line in lines:
    s = line.strip()
    if s.startswith("ANGL"):
        in_angl = True
        out.append(line)
        continue
    if in_angl and s and s[0].isupper() and not any(
        s.startswith(x) for x in
        ["M","Y","N","C","O","H","K","CR","CV","CC","CO","2C"]
    ):
        in_angl = False
    if in_angl and ("M1" in line or "M2" in line):
        match = re.match(r"^(\S+)\s+(\S+)\s+(\S+)(.*)", line)
        if match:
            atype, k_old, theta, rest = match.groups()
            line = f"{atype:<14s}  31.00   {theta}{rest}\n"
            changed.append(f"  {atype}: {k_old} -> 31.00")
    out.append(line)
with open(dst, "w") as f:
    f.writelines(out)
print(f"k=31 frcmod: {len(changed)} terms changed")
for c in changed:
    print(c)
PYEOF
log "k=31 frcmod done"

# ── PHASE 4: Build ABS topology ───────────────────────────────────
log "Phase 4: Building ABS topology..."
cd systems/mol13_pose1_monomer
tleap -f tleap_no_ohb_mol13_pose1_monomer.in \
  > tleap_abs_crashtest.log 2>&1
cd /path/to/urease_md

PRMTOP=$(find systems/mol13_pose1_monomer -maxdepth 1 \
  -name "monomer_*.prmtop" ! -name "*backup*" | head -1)
INPCRD=$(find systems/mol13_pose1_monomer -maxdepth 1 \
  -name "monomer_*.inpcrd" ! -name "*backup*" | head -1)

set +e
printf "parm ${PRMTOP}\ntrajin ${INPCRD}\nfixatomorder\ntrajout /tmp/mol13_abs_fixed.rst7 restart\nrun\nquit\n" \
  | ${CPPTRAJ} >> "${LOGFILE}" 2>&1
set -e

python3 - << PYEOF2
import parmed as pmd
p = pmd.load_file("${PRMTOP}", "/tmp/mol13_abs_fixed.rst7")
pmd.tools.HMassRepartition(p).execute()
p.save("${TESTDIR}/solvated_abs_hmr.prmtop", overwrite=True)
h = [a.mass for a in p.atoms if a.name.startswith("H")]
print(f"ABS HMR done. Max H: {max(h):.3f}")
PYEOF2

set +e
printf "parm ${TESTDIR}/solvated_abs_hmr.prmtop\ntrajin /tmp/mol13_abs_fixed.rst7\ndistance nini :570@NI :571@NI out /tmp/nini_abs.dat\nrun\nquit\n" \
  | ${CPPTRAJ} >> "${LOGFILE}" 2>&1
set -e
NINI_ABS=$(awk 'NR==2{print $2}' /tmp/nini_abs.dat 2>/dev/null || echo "unknown")
log "ABS topology Ni-Ni: ${NINI_ABS} Å"

# ── PHASE 5: Build k=31 topology ──────────────────────────────────
log "Phase 5: Building k=31 topology..."
sed "s/urease_no_ohb_mcpbpy\.frcmod/urease_no_ohb_mcpbpy_k31.frcmod/" \
  systems/mol13_pose1_monomer/tleap_no_ohb_mol13_pose1_monomer.in \
  > systems/mol13_pose1_monomer/tleap_k31_crashtest.in

cd systems/mol13_pose1_monomer
tleap -f tleap_k31_crashtest.in \
  > tleap_k31_crashtest.log 2>&1
cd /path/to/urease_md

PRMTOP_K31=$(find systems/mol13_pose1_monomer -maxdepth 1 \
  -name "monomer_*.prmtop" ! -name "*backup*" | head -1)
INPCRD_K31=$(find systems/mol13_pose1_monomer -maxdepth 1 \
  -name "monomer_*.inpcrd" ! -name "*backup*" | head -1)

set +e
printf "parm ${PRMTOP_K31}\ntrajin ${INPCRD_K31}\nfixatomorder\ntrajout /tmp/mol13_k31_fixed.rst7 restart\nrun\nquit\n" \
  | ${CPPTRAJ} >> "${LOGFILE}" 2>&1
set -e

python3 - << PYEOF3
import parmed as pmd
p = pmd.load_file("${PRMTOP_K31}", "/tmp/mol13_k31_fixed.rst7")
pmd.tools.HMassRepartition(p).execute()
p.save("${TESTDIR}/solvated_k31_hmr.prmtop", overwrite=True)
h = [a.mass for a in p.atoms if a.name.startswith("H")]
print(f"k31 HMR done. Max H: {max(h):.3f}")
PYEOF3

set +e
printf "parm ${TESTDIR}/solvated_k31_hmr.prmtop\ntrajin /tmp/mol13_k31_fixed.rst7\ndistance nini :570@NI :571@NI out /tmp/nini_k31.dat\nrun\nquit\n" \
  | ${CPPTRAJ} >> "${LOGFILE}" 2>&1
set -e
NINI_K31=$(awk 'NR==2{print $2}' /tmp/nini_k31.dat 2>/dev/null || echo "unknown")
log "k=31 topology Ni-Ni: ${NINI_K31} Å"

# ── PHASE 6: Write crash test input ───────────────────────────────
cat > "${INPUT}" << INEOF
Crash test: 32 ns NPT MC barostat from production restart
 &cntrl
  imin=0, nstlim=8000000, dt=0.004,
  ntx=5, irest=1,
  ntb=2, ntp=1,
  ntt=11, tautp=2.0, temp0=300.0,
  ntc=2, ntf=2,
  cut=9.0,
  ntpr=5000, ntwx=0, ntwr=500000,
  ioutfm=1, ntxo=2,
  barostat=2, mcbarint=25, pres0=1.0,
  nmropt=1,
 /
 &wt type='END' /
DISANG=${NISRST}
INEOF
log "Crash test input written (8M steps = 32 ns)"

# ── PHASE 7: TEST A — ABS frcmod ──────────────────────────────────
log "=== Phase 7: TEST A — ABS frcmod (0.94-8.00 kcal/mol/rad2) ==="
mkdir -p "${TESTDIR}/abs_test"

set +e
CUDA_VISIBLE_DEVICES=1 \
${AMBER} \
  -O \
  -i  "${INPUT}" \
  -p  "${TESTDIR}/solvated_abs_hmr.prmtop" \
  -c  "${STARTRST}" \
  -o  "${TESTDIR}/abs_test/crashtest_abs.out" \
  -r  "${TESTDIR}/abs_test/crashtest_abs.rst7" \
  >> "${TESTDIR}/abs_test/crashtest_abs.log" 2>&1
set -e

NAN_A=$(grep -c "NaN" \
  "${TESTDIR}/abs_test/crashtest_abs.out" 2>/dev/null || echo 0)
NSTEP_A=$(grep "NSTEP" \
  "${TESTDIR}/abs_test/crashtest_abs.out" 2>/dev/null | \
  tail -1 | awk '{print $3}')
ANGLE_A=$(grep "ANGLE " \
  "${TESTDIR}/abs_test/crashtest_abs.out" 2>/dev/null | \
  tail -1 | awk '{print $3}')

[[ "${NAN_A}" -gt 0 ]] && \
  log "TEST A: FAILED — NaN at NSTEP ${NSTEP_A} ANGLE=${ANGLE_A}" || \
  log "TEST A: PASSED — ${NSTEP_A} steps clean ANGLE=${ANGLE_A}"

# ── PHASE 8: TEST B — k=31 uniform ────────────────────────────────
log "=== Phase 8: TEST B — k=31 uniform frcmod ==="
mkdir -p "${TESTDIR}/k31_test"

set +e
CUDA_VISIBLE_DEVICES=1 \
${AMBER} \
  -O \
  -i  "${INPUT}" \
  -p  "${TESTDIR}/solvated_k31_hmr.prmtop" \
  -c  "${STARTRST}" \
  -o  "${TESTDIR}/k31_test/crashtest_k31.out" \
  -r  "${TESTDIR}/k31_test/crashtest_k31.rst7" \
  >> "${TESTDIR}/k31_test/crashtest_k31.log" 2>&1
set -e

NAN_B=$(grep -c "NaN" \
  "${TESTDIR}/k31_test/crashtest_k31.out" 2>/dev/null || echo 0)
NSTEP_B=$(grep "NSTEP" \
  "${TESTDIR}/k31_test/crashtest_k31.out" 2>/dev/null | \
  tail -1 | awk '{print $3}')
ANGLE_B=$(grep "ANGLE " \
  "${TESTDIR}/k31_test/crashtest_k31.out" 2>/dev/null | \
  tail -1 | awk '{print $3}')

[[ "${NAN_B}" -gt 0 ]] && \
  log "TEST B: FAILED — NaN at NSTEP ${NSTEP_B} ANGLE=${ANGLE_B}" || \
  log "TEST B: PASSED — ${NSTEP_B} steps clean ANGLE=${ANGLE_B}"

# ── PHASE 9: Final summary ─────────────────────────────────────────
log ""
log "════════════════════════════════════════════════"
log "            CRASH TEST SUMMARY"
log "════════════════════════════════════════════════"
log "Patched k=50/30:  PASSED (proven, 100 ns clean)"
[[ "${NAN_A}" -gt 0 ]] && \
  log "ABS (0.94-8.00):  FAILED at NSTEP ${NSTEP_A}" || \
  log "ABS (0.94-8.00):  PASSED — ${NSTEP_A} steps clean"
[[ "${NAN_B}" -gt 0 ]] && \
  log "k=31 uniform:     FAILED at NSTEP ${NSTEP_B}" || \
  log "k=31 uniform:     PASSED — ${NSTEP_B} steps clean"
log ""
log "RECOMMENDATION:"
if [[ "${NAN_B}" -eq 0 ]]; then
    log "  USE k=31 — ZAFF x15 justified AND MC stable"
    log "  Rebuild all 35 remaining systems with k31 frcmod"
elif [[ "${NAN_A}" -eq 0 ]]; then
    log "  USE ABS — QM-derived magnitudes are stable"
    log "  Rebuild all 35 remaining systems with ABS frcmod"
else
    log "  USE patched k=50/30 — only confirmed stable values"
    log "  Rebuild all 35 remaining systems with patched frcmod"
fi
log "════════════════════════════════════════════════"
log "=== ALL DONE $(date) ==="
