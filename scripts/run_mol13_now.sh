#!/usr/bin/env bash
# scripts/run_mol13_now.sh
# mol13: seg1 = pre-existing crashtest/patched_test; run seg2–5 + production_rep1–3
set -euo pipefail
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/post_equil_lib.sh
source "${_SCRIPT_DIR}/post_equil_lib.sh"
source "${_SCRIPT_DIR}/../config.env"

export PATH="/usr/bin:${AMBER_HOME}/bin:/usr/local/bin:${PATH}"
export CUDA_DEVICE_ORDER="${CUDA_DEVICE_ORDER:-PCI_BUS_ID}"
export CUDA_VISIBLE_DEVICES=1

PMEMD="${PMEMD_CUDA_BIN:-/opt/software/AMBER/pmemd24_cuda128/bin/pmemd.cuda}"
[[ -x "${PMEMD}" ]] || {
  echo "[FAIL] Not executable: ${PMEMD}"
  exit 1
}

SYS="mol13_pose1_monomer"
SIM="${PROJECT_ROOT}/simulations/${SYS}"
PRMTOP="${PROJECT_ROOT}/systems/${SYS}/solvated_hmr.prmtop"
SEG1_RST="${SIM}/crashtest/patched_test/patched_test.rst7"
RST_OUT="${SIM}/restraints/ni_ni.RST"
MDIN_EQ="${SIM}/input_files/equil_seg.in"
MDIN_PR="${SIM}/input_files/production.in"
MASK_TMP="${SIM}/restraints/.ni_mask_cpptraj.txt"
# pmemd -i/-p/-c use short symlinks via post_equil_lib.sh (Fortran 80-char limit)

LOG="${PROJECT_ROOT}/logs/mol13_run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

echo "=== run_mol13_now.sh ==="
echo "LOG: ${LOG}"
echo "PMEMD: ${PMEMD}  CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"

[[ -f "${SEG1_RST}" ]] || {
  echo "[FAIL] seg1 restart missing: ${SEG1_RST}"
  exit 1
}
[[ -f "${PRMTOP}" ]] || {
  echo "[FAIL] Missing topology ${PRMTOP}"
  exit 1
}

# Coordinates for @NI mask (any frame matching this prmtop)
COORD_MASK=""
for cand in \
  "${SIM}/gpu_bridge2/gpu_bridge2.rst7" \
  "${SEG1_RST}" \
  "${SIM}/equil3/equil3.rst7"; do
  if [[ -f "${cand}" ]]; then
    COORD_MASK="${cand}"
    break
  fi
done
[[ -n "${COORD_MASK}" ]] || {
  echo "[FAIL] No rst7 found for cpptraj mask extraction"
  exit 1
}
echo "Using coordinates for NI mask: ${COORD_MASK}"

if ! post_equil_extract_ni_atoms "${PRMTOP}" "${COORD_MASK}" "${MASK_TMP}"; then
  exit 1
fi

post_equil_write_ni_ni_rst "${RST_OUT}" "${NI1_ATOM}" "${NI2_ATOM}"
DISANG_ABS="$(post_equil_disang_short_abs "${SYS}" "${RST_OUT}")"
post_equil_write_mdin_pair "${DISANG_ABS}" "${MDIN_EQ}" "${MDIN_PR}"
post_equil_ensure_short_links "${SYS}" || exit 1

echo "Built ${RST_OUT} (iat=${NI1_ATOM},${NI2_ATOM}) and mdin + short m/p/r symlink layout."

SEG3C_DIR="${SIM}/equil_seg3c"
SEG3C_PRE="equil_seg3c"
SEG2_RST="${SIM}/equil_seg2/equil_seg2.rst7"

declare -A SEG_STAT SEG_NAN
declare -A REP_STAT REP_NAN

run_equil_seg() {
  local n=$1 crd_in=$2
  local dir="${SIM}/equil_seg${n}"
  local pre="equil_seg${n}"
  if [[ -s "${dir}/${pre}.rst7" ]]; then
    echo "seg${n}: already done, skipping"
    SEG_STAT[${n}]="SKIP"
    SEG_NAN[${n}]="—"
    return 0
  fi
  echo "--- seg${n}: starting from ${crd_in} ---"
  set +e
  post_equil_run_pmemd_segment "${dir}" "${SYS}" equil "${crd_in}" "${pre}" "${PMEMD}"
  local rc=$?
  set -e
  if [[ "${rc}" -ne 0 ]]; then
    echo "[FAIL] seg${n}: pmemd exit code ${rc}"
    exit 1
  fi
  if ! post_equil_verify_mdout "${dir}/${pre}.out" 2500000; then
    echo "[FAIL] seg${n}: mdout verification failed"
    exit 1
  fi
  SEG_STAT[${n}]="DONE"
  SEG_NAN[${n}]="0"
}

# Re-run seg3 into equil_seg3c/ from clean seg2.rst7 (original equil_seg3/ is preserved).
run_seg3c_from_seg2() {
  local logf="${PROJECT_ROOT}/logs/mol13_seg3c_worker.log"
  local pidf="${SEG3C_DIR}/seg3c_pmemd.pid"

  [[ -s "${SEG2_RST}" ]] || {
    echo "[FAIL] seg3c needs clean restart: ${SEG2_RST}"
    exit 1
  }

  if [[ -f "${SEG3C_DIR}/${SEG3C_PRE}.out" ]] && [[ -s "${SEG3C_DIR}/${SEG3C_PRE}.rst7" ]]; then
    if post_equil_verify_mdout "${SEG3C_DIR}/${SEG3C_PRE}.out" 2500000; then
      echo "seg3c: already complete, skipping"
      return 0
    fi
    echo "[FAIL] ${SEG3C_DIR} has partial outputs; remove equil_seg3c/ and re-run or wait for a running job to finish."
    exit 1
  fi

  mkdir -p "${SEG3C_DIR}"
  post_equil_ensure_short_links "${SYS}" "${SEG2_RST}" || exit 1

  local need_launch=1
  if [[ -f "${pidf}" ]]; then
    local opid
    opid=$(tr -d ' \n' <"${pidf}" || true)
    if [[ -n "${opid}" ]] && kill -0 "${opid}" 2>/dev/null; then
      need_launch=0
      echo "seg3c: pmemd already running (pid ${opid}), waiting for completion..."
    else
      rm -f "${pidf}"
    fi
  fi

  if [[ "${need_launch}" -eq 1 ]]; then
    echo "--- seg3c: launching from ${SEG2_RST} (background, log ${logf}) ---"
    nohup bash -c "
set -euo pipefail
export CUDA_VISIBLE_DEVICES=1
cd \"${SEG3C_DIR}\"
\"${PMEMD}\" -O \\
  -i \"${PROJECT_ROOT}/m/${SYS}_eq.in\" \\
  -p \"${PROJECT_ROOT}/p/${SYS}.prmtop\" \\
  -c \"${PROJECT_ROOT}/r/${SYS}_in.rst7\" \\
  -o \"${SEG3C_PRE}.out\" \\
  -r \"${SEG3C_PRE}.rst7\" \\
  -x \"${SEG3C_PRE}.nc\" &
PMEMD_PID=\$!
echo \"\${PMEMD_PID}\" > seg3c_pmemd.pid
wait \"\${PMEMD_PID}\"
rc=\$?
rm -f seg3c_pmemd.pid
if [[ -f mdinfo ]]; then
  mv -f mdinfo \"${SEG3C_PRE}.mdinfo\"
fi
exit \"\$rc\"
" >> "${logf}" 2>&1 &
    echo "seg3c launcher PID: $!"
    sleep 5
  fi

  echo "Waiting for seg3c to complete..."
  while true; do
    if post_equil_verify_mdout "${SEG3C_DIR}/${SEG3C_PRE}.out" 2500000 2>/dev/null; then
      break
    fi
    local live=0
    if [[ -f "${pidf}" ]]; then
      local opid2
      opid2=$(tr -d ' \n' <"${pidf}" || true)
      if [[ -n "${opid2}" ]] && kill -0 "${opid2}" 2>/dev/null; then
        live=1
      fi
    fi
    if [[ "${live}" -eq 0 ]]; then
      sleep 15
      if post_equil_verify_mdout "${SEG3C_DIR}/${SEG3C_PRE}.out" 2500000 2>/dev/null; then
        break
      fi
      echo "[FAIL] seg3c process died before mdout showed NSTEP=2500000 (see ${SEG3C_DIR}/${SEG3C_PRE}.out and ${logf})"
      exit 1
    fi
    sleep 60
  done

  echo "seg3c complete, proceeding to seg4"
}

run_prod_rep() {
  local n=$1 crd_in=$2
  local dir="${SIM}/production_rep${n}"
  local pre="production_rep${n}"
  if [[ -s "${dir}/${pre}.rst7" ]]; then
    echo "rep${n}: already done, skipping"
    REP_STAT[${n}]="SKIP"
    REP_NAN[${n}]="—"
    return 0
  fi
  echo "--- rep${n}: starting from ${crd_in} ---"
  set +e
  post_equil_run_pmemd_segment "${dir}" "${SYS}" prod "${crd_in}" "${pre}" "${PMEMD}"
  local rc=$?
  set -e
  if [[ "${rc}" -ne 0 ]]; then
    echo "[FAIL] rep${n}: pmemd exit code ${rc}"
    exit 1
  fi
  if ! post_equil_verify_mdout "${dir}/${pre}.out" 5000000; then
    echo "[FAIL] rep${n}: mdout verification failed"
    exit 1
  fi
  REP_STAT[${n}]="DONE"
  REP_NAN[${n}]="0"
}

# seg2 from seg1 (patched_test.rst7); do not move patched_test.
# Do NOT re-run into equil_seg3/ — that tree is kept for records; continuation uses clean equil_seg3c/ from seg2.
run_equil_seg 2 "${SEG1_RST}"
run_seg3c_from_seg2
run_equil_seg 4 "${SIM}/equil_seg3c/equil_seg3c.rst7"
run_equil_seg 5 "${SIM}/equil_seg4/equil_seg4.rst7"

SEG5_RST="${SIM}/equil_seg5/equil_seg5.rst7"
[[ -s "${SEG5_RST}" ]] || {
  echo "[FAIL] Missing ${SEG5_RST} — cannot start production"
  exit 1
}

run_prod_rep 1 "${SEG5_RST}"
run_prod_rep 2 "${SEG5_RST}"
run_prod_rep 3 "${SEG5_RST}"

fmt_seg() {
  local n=$1
  local st="${SEG_STAT[${n}]:-?}"
  if [[ "${st}" == "DONE" ]]; then
    echo "seg${n}:  DONE  (NaN=${SEG_NAN[${n}]}, NSTEP=2500000)"
  elif [[ "${st}" == "SKIP" ]]; then
    echo "seg${n}:  SKIP  (existing rst7)"
  else
    echo "seg${n}:  ${st}"
  fi
}

fmt_rep() {
  local n=$1
  local st="${REP_STAT[${n}]:-?}"
  if [[ "${st}" == "DONE" ]]; then
    echo "rep${n}:  DONE  (NaN=${REP_NAN[${n}]}, NSTEP=5000000)"
  elif [[ "${st}" == "SKIP" ]]; then
    echo "rep${n}:  SKIP  (existing rst7)"
  else
    echo "rep${n}:  ${st}"
  fi
}

echo ""
echo "=== ${SYS} COMPLETE ==="
echo "seg1:  DONE (pre-existing patched_test)"
echo "seg3:  KEEP (original equil_seg3/ trajectory preserved — not used after seg3c)"
fmt_seg 2
echo "seg3c: DONE (clean 2.5M-step NPT from seg2.rst7 → equil_seg3c.rst7)"
fmt_seg 4
fmt_seg 5
fmt_rep 1
fmt_rep 2
fmt_rep 3
echo "Total production trajectory: 60 ns (3 × 20 ns)"
