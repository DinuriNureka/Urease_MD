#!/usr/bin/env bash
# Top-up mol13 from seg08 rst7 -> seg09 -> seg10 to reach 100 ns
source config.env
export CUDA_VISIBLE_DEVICES=1
export PATH="/usr/bin:${AMBER_HOME}/bin:/usr/local/bin:${PATH}"

BASEDIR="${PROJECT_ROOT}"
SYS="mol13_pose1_monomer"
PMEMD="${PMEMD_CUDA_BIN:-/opt/software/AMBER/pmemd24_cuda128/bin/pmemd.cuda}"
PRMTOP="${BASEDIR}/systems/${SYS}/solvated_hmr.prmtop"
EQUIL_NSTLIM=2500000

# Find NI restraint symlink
RST_FILE=""
for p in "${BASEDIR}/rst/${SYS}.RST" \
          "${BASEDIR}/simulations/${SYS}/restraints/ni_ni.RST"; do
  if [[ -f "$p" ]]; then RST_FILE="$p"; break; fi
done
[[ -z "$RST_FILE" ]] && { echo "[FAIL] No RST file for ${SYS}"; exit 1; }
echo "[INFO] RST file: ${RST_FILE} ($(wc -c < "$RST_FILE") bytes)"

# Function to write mdin with DISANG embedded (mirrors run_resilient_queue.sh)
write_mdin() {
  local out=$1 disang=$2
  mkdir -p "$(dirname "$out")"
  cat > "$out" << MDIN
Topup mol13 — 10 ns NPT (ntt=11, MC barostat, nmropt=1)
 &cntrl
  imin=0, nstlim=${EQUIL_NSTLIM}, dt=0.004
  ntx=5, irest=1,
  ntb=2, ntp=1,
  ntt=11, tautp=2.0, temp0=300.0,
  ntc=2, ntf=2,
  cut=9.0,
  fswitch=8.0,
  ntpr=5000, ntwx=5000, ntwr=5000,
  ioutfm=1, ntxo=2,
  barostat=2, mcbarint=100, pres0=1.0,
  nmropt=1,
  ntr=0,
  ig=-1,
  iwrap=1,
 /
 &wt type='END' /
DISANG=${disang}
MDIN
}

# Start from seg08 rst7
PREV_RST="${BASEDIR}/simulations/${SYS}/equil_seg07/equil_seg07.rst7"
[[ -f "$PREV_RST" ]] || { echo "[FAIL] seg08 rst7 missing"; exit 1; }

for SEG in equil_seg08 equil_seg09 equil_seg10; do
  SEGDIR="${BASEDIR}/simulations/${SYS}/${SEG}"
  mkdir -p "${SEGDIR}"

  # Skip if already complete
  if [[ -f "${SEGDIR}/${SEG}.rst7" ]]; then
    SIZE=$(stat -c%s "${SEGDIR}/${SEG}.rst7" 2>/dev/null || echo 0)
    if [[ "$SIZE" -gt 1000000 ]]; then
      echo "[SKIP] ${SEG} already has rst7 (${SIZE} bytes)"
      PREV_RST="${SEGDIR}/${SEG}.rst7"
      continue
    fi
  fi

  # Write mdin with DISANG embedded
  MDIN="/tmp/mol13_topup_${SEG}.in"
  write_mdin "${MDIN}" "${RST_FILE}"
  echo "[INFO] mdin written: ${MDIN}"
  echo "[INFO] DISANG line: $(grep DISANG ${MDIN})"

  echo "[RUN] ${SEG} from ${PREV_RST}"
  cd "${SEGDIR}"

  "${PMEMD}" -O \
    -i "${MDIN}" \
    -o "${SEG}.out" \
    -p "${PRMTOP}" \
    -c "${PREV_RST}" \
    -r "${SEG}.rst7" \
    -x "${SEG}.nc"

  RC=$?
  NAN_COUNT=$(grep -ac NaN "${SEG}.out" 2>/dev/null || echo 0)
  echo "[INFO] ${SEG} exit=${RC} nan=${NAN_COUNT}"

  cd "${BASEDIR}"

  if [[ ! -f "${SEGDIR}/${SEG}.rst7" ]]; then
    echo "[FAIL] ${SEG} rst7 not written — aborting"
    exit 1
  fi

  PREV_RST="${SEGDIR}/${SEG}.rst7"
  echo "[DONE] ${SEG}"
done

echo "[TOPUP COMPLETE] mol13 seg09+seg10 done. Check TIME in seg10 should be ~100803 ps"
# NOTE: patch applied below — use short mdin path
