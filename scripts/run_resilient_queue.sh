#!/usr/bin/env bash
# scripts/run_resilient_queue.sh
# Sequential resilient MD queue on GPU 1 (NaN-aware retries, resumable skips).
# Per system:
#   Equil: 10 × 10 ns = 100 ns  →  simulations/<sys>/equil_seg01 … equil_seg10
#   Prod:  3 replicates × 2 × 10 ns = 60 ns → prod_rep01_seg01, prod_rep01_seg02, … prod_rep03_seg02
# All replicates start from equil_seg10.rst7. Do not run concurrent pmemd on the same GPU.
#
# Directory layout (each segment: <name>.{out,rst7,nc}):
#   simulations/<sys>/equil_seg01/equil_seg01.{out,rst7,nc}
#   …
#   simulations/<sys>/equil_seg10/equil_seg10.{out,rst7,nc}
#   simulations/<sys>/prod_rep01_seg01/prod_rep01_seg01.{out,rst7,nc}
#   simulations/<sys>/prod_rep01_seg02/prod_rep01_seg02.{out,rst7,nc}
#   … (prod_rep02_*, prod_rep03_*)
#
# Launch (when ready):
#   nohup bash scripts/run_resilient_queue.sh > logs/resilient_queue.log 2>&1 &
set -euo pipefail
set -o pipefail

_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/post_equil_lib.sh
source "${_SCRIPT_DIR}/post_equil_lib.sh"
source "${_SCRIPT_DIR}/../config.env"

export PATH="/usr/bin:${AMBER_HOME}/bin:/usr/local/bin:${PATH}"
export CUDA_DEVICE_ORDER="${CUDA_DEVICE_ORDER:-PCI_BUS_ID}"
export CUDA_VISIBLE_DEVICES=1

# Segment protocol (must match simulations/input_files/equil_free.in and production.in nstlim)
EQUIL_SEGS=10
EQUIL_NSTLIM=2500000
PROD_REPS=3
PROD_SEGS_PER_REP=2
# ~98% of EQUIL_NSTLIM — near-complete trajectory recovery threshold
NEAR_COMPLETE_NSTEP=2450000

PMEMD="${PMEMD_CUDA_BIN:-/opt/software/AMBER/pmemd24_cuda128/bin/pmemd.cuda}"
[[ -x "${PMEMD}" ]] || {
  echo "[FAIL] Not executable: ${PMEMD}" >&2
  exit 1
}

QUEUE=(
  mol1_pose2_monomer
  mol13_pose1_monomer
  mol14_pose2_monomer
  mol19_pose1_monomer
  mol33_pose1_monomer
)

BASEDIR="${PROJECT_ROOT}"
MQUEUE="${BASEDIR}/m"
PQUEUE="${BASEDIR}/p"
INP="${BASEDIR}/simulations/input_files"
mkdir -p "${MQUEUE}" "${PQUEUE}" "${BASEDIR}/logs"

TS="$(date +%Y%m%d_%H%M%S)"
MAIN_LOG="${BASEDIR}/logs/resilient_queue_${TS}.log"
exec > >(tee -a "${MAIN_LOG}") 2>&1

echo "=== run_resilient_queue.sh === ${TS}"
echo "MAIN_LOG=${MAIN_LOG}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"
echo "PMEMD=${PMEMD}"
echo "EQUIL_SEGS=${EQUIL_SEGS} × 10 ns | PROD_REPS=${PROD_REPS} × ${PROD_SEGS_PER_REP} × 10 ns"

resolve_ni_rst() {
  local sys=$1 p
  for p in \
    "${BASEDIR}/rst/${sys}.RST" \
    "${BASEDIR}/simulations/${sys}/restraints/ni_ni.RST"; do
    if [[ -f "${p}" ]]; then
      # Do not use readlink -f: symlinks under rst/ would resolve to long
      # systems/.../ni_ni.RST paths (>80 chars) and pmemd truncates DISANG= at 80.
      echo "${p}"
      return 0
    fi
  done
  echo "[FAIL] No ni_ni.RST for ${sys}" >&2
  return 1
}

ensure_prmtop_link() {
  local sys=$1
  local prm="${BASEDIR}/systems/${sys}/solvated_hmr.prmtop"
  [[ -f "${prm}" ]] || {
    echo "[FAIL] Missing ${prm}" >&2
    return 1
  }
  ln -sf "$(readlink -f "${prm}" 2>/dev/null || echo "${prm}")" "${PQUEUE}/${sys}.prmtop"
  post_equil__assert_path_short "prmtop" "${PQUEUE}/${sys}.prmtop" || return 1
}

# Unique short mdin link per segment (avoids path-length / collision issues)
link_mdin_short() {
  local sys=$1 tag=$2 real_mdin=$3
  local s="${MQUEUE}/${sys}_${tag}.in"
  ln -sf "$(readlink -f "${real_mdin}" 2>/dev/null || echo "${real_mdin}")" "${s}"
  post_equil__assert_path_short "mdin -i" "${s}" || return 1
  printf '%s' "${s}"
}

link_crd_short() {
  local sys=$1 real_crd=$2
  post_equil_short_crd_path "${sys}" "${real_crd}"
}

# Mdin body matches equil_free / production protocol (DISANG set to absolute path).
write_mdin_segment() {
  local out=$1 disang_abs=$2
  mkdir -p "$(dirname "${out}")"
  local mdin
  mdin=$'Resilient queue — 10 ns NPT (ntt=11, MC barostat, nmropt=1)\n'
  mdin+=$' &cntrl\n'
  mdin+="  imin=0, nstlim=${EQUIL_NSTLIM}, dt=0.004"$'\n'
  mdin+=$'  ntx=5, irest=1,\n'
  mdin+=$'  ntb=2, ntp=1,\n'
  mdin+=$'  ntt=11, tautp=2.0, temp0=300.0,\n'
  mdin+=$'  ntc=2, ntf=2,\n'
  mdin+=$'  cut=9.0,\n'
  mdin+=$'  fswitch=8.0,\n'
  mdin+=$'  ntpr=5000, ntwx=5000, ntwr=5000,\n'
  mdin+=$'  ioutfm=1, ntxo=2,\n'
  mdin+=$'  barostat=2, mcbarint=100, pres0=1.0,\n'
  mdin+=$'  nmropt=1,\n'
  mdin+=$'  ntr=0,\n'
  mdin+=$'  ig=-1,\n'
  mdin+=$'  iwrap=1,\n'
  mdin+=$' /\n'
  mdin+=$' &wt type=\'END\' /\n'
  mdin+="DISANG=${disang_abs}"$'\n'
  printf '%s' "${mdin}" >"${out}"
}

cpptraj_nini_ok() {
  local prmtop=$1 rst7=$2
  local cpptraj tdat tlog
  cpptraj="$(post_equil__cpptraj_bin)"
  tdat="$(mktemp /tmp/rq_ni_XXXXXX.dat)"
  tlog="$(mktemp /tmp/rq_cpp_XXXXXX.log)"
  "${cpptraj}" -p "${prmtop}" <<EOF >"${tlog}" 2>&1
trajin ${rst7}
distance rq_ni :${NI1_RES:-570}@NI :${NI2_RES:-571}@NI
run
writedata ${tdat} rq_ni
quit
EOF
  local v
  v=$(awk 'NR==2{print $2}' "${tdat}" 2>/dev/null || echo nan)
  rm -f "${tdat}" "${tlog}"
  [[ "${v}" != "nan" ]] && [[ -n "${v}" ]] && [[ "${v}" != "NaN" ]]
}

mdout_nan_count() {
  local out=$1 n
  [[ -f "${out}" ]] || {
    echo 999
    return 0
  }
  n=$(grep -c 'NaN' "${out}" 2>/dev/null || true)
  echo "${n:-0}"
}

last_clean_nstep() {
  local out=$1 v
  [[ -f "${out}" ]] || {
    echo 0
    return 0
  }
  v=$(grep 'NSTEP' "${out}" | grep -v NaN | tail -1 | sed -n 's/.*NSTEP[[:space:]]*=[[:space:]]*\([0-9][0-9]*\).*/\1/p')
  [[ -n "${v}" ]] && echo "${v}" || echo 0
}

newest_rst7_in_dir() {
  local d=$1
  ls -t "${d}"/*.rst7 2>/dev/null | head -1 || true
}

run_pmemd_once() {
  local run_dir=$1 mdin_short=$2 prm_short=$3 crd_short=$4 prefix=$5
  mkdir -p "${run_dir}"
  (
    cd "${run_dir}" || exit 1
    "${PMEMD}" -O \
      -i "${mdin_short}" \
      -p "${prm_short}" \
      -c "${crd_short}" \
      -o "${prefix}.out" \
      -r "${prefix}.rst7" \
      -x "${prefix}.nc"
    rc=$?
    [[ -f mdinfo ]] && mv -f mdinfo "${prefix}.mdinfo" 2>/dev/null || true
    exit "${rc}"
  )
}

# Returns: 0 success, 1 hard fail, 2 caller should retry (clear segment outputs, adjust crd per caller)
run_segment() {
  local sys=$1 phase_label=$2 run_dir=$3 prefix=$4 crd_initial=$5 disang_abs=$6 prm_short=$7
  local mdin_real="${run_dir}/${prefix}.in"
  local mdout="${run_dir}/${prefix}.out"
  local prm_real="${BASEDIR}/systems/${sys}/solvated_hmr.prmtop"
  local mdin_short crd_short crd_current
  local early_retry=0 rerun_from_prev=0

  write_mdin_segment "${mdin_real}" "${disang_abs}"
  mdin_short="$(link_mdin_short "${sys}" "${phase_label//\//_}" "${mdin_real}")" || return 1

  crd_current="${crd_initial}"

  while true; do
    crd_short="$(link_crd_short "${sys}" "${crd_current}")" || return 1
    echo "[RUN] ${sys} ${phase_label} crd=${crd_current} early_retry=${early_retry} rerun_from_prev=${rerun_from_prev}"

    set +e
    run_pmemd_once "${run_dir}" "${mdin_short}" "${prm_short}" "${crd_short}" "${prefix}"
    local rc=$?
    set -e

    local nan_c
    nan_c=$(mdout_nan_count "${mdout}")

    if [[ "${nan_c}" -eq 0 ]] && [[ "${rc}" -eq 0 ]]; then
      echo "[OK] ${sys} ${phase_label} completed (${EQUIL_NSTLIM} steps × dt=0.004 ps)"
      return 0
    fi

    if [[ "${nan_c}" -eq 0 ]] && [[ "${rc}" -ne 0 ]]; then
      echo "[WARN] ${sys} ${phase_label} pmemd rc=${rc}, NaN count 0 — treating as failure"
    fi

    echo "[NAN] ${sys} ${phase_label} nan_lines=${nan_c} pmemd_rc=${rc}"

    local last_clean
    last_clean=$(last_clean_nstep "${mdout}")
    last_clean="${last_clean:-0}"
    echo "[INFO] last_clean_nstep=${last_clean}"

    if [[ "${last_clean}" -ge "${NEAR_COMPLETE_NSTEP}" ]]; then
      local cand
      cand=$(newest_rst7_in_dir "${run_dir}")
      if [[ -n "${cand}" ]] && cpptraj_nini_ok "${prm_real}" "${cand}"; then
        echo "[RECOVER] ${sys} ${phase_label} accept near-complete rst7 ${cand}"
        return 0
      fi
      echo "[RECOVER] ${sys} ${phase_label} near-complete rst7 missing or cpptraj NiNi invalid — rerun-from-prev path"
      rerun_from_prev=$((rerun_from_prev + 1))
      if [[ "${rerun_from_prev}" -le 2 ]]; then
        echo "[RERUN] ${sys} ${phase_label} signal caller: restart segment from upstream rst7 (${rerun_from_prev}/2)"
        return 2
      fi
      echo "[FAIL] ${sys} ${phase_label} after near-complete recovery attempts"
      return 1
    fi

    if [[ "${early_retry}" -lt 3 ]]; then
      early_retry=$((early_retry + 1))
      echo "[EARLY CRASH] ${sys} ${phase_label} at NSTEP=${last_clean} — retry ${early_retry}/3 from same crd"
      rm -f "${mdout}" "${run_dir}/${prefix}.rst7" "${run_dir}/${prefix}.nc" 2>/dev/null || true
      crd_current="${crd_initial}"
      continue
    fi

    rerun_from_prev=$((rerun_from_prev + 1))
    if [[ "${rerun_from_prev}" -le 2 ]]; then
      echo "[RERUN] ${sys} ${phase_label} signal caller: full segment rerun from upstream (${rerun_from_prev}/2)"
      return 2
    fi

    echo "[FAIL] ${sys} ${phase_label} NaN recovery exhausted"
    return 1
  done
}

rst7_skip_ok() {
  local prm=$1 rst=$2
  [[ -s "${rst}" ]] || return 1
  cpptraj_nini_ok "${prm}" "${rst}"
}

# First equil_free segment: start after standard equil pipeline (prefer equil3).
find_starting_equil_rst() {
  local sim=$1 p
  for p in \
    "${sim}/equil3/equil3.rst7" \
    "${sim}/gpu_bridge2/gpu_bridge2.rst7" \
    "${sim}/equil2/equil2.rst7" \
    "${sim}/equil1/equil1.rst7"; do
    if [[ -f "${p}" ]]; then
      printf '%s' "${p}"
      return 0
    fi
  done
  return 1
}

run_one_system() {
  local sys=$1
  local sim="${BASEDIR}/simulations/${sys}"
  local prm="${BASEDIR}/systems/${sys}/solvated_hmr.prmtop"
  local disang prm_short
  disang="$(resolve_ni_rst "${sys}")" || return 1

  echo "========== SYSTEM ${sys} =========="
  [[ -f "${prm}" ]] || {
    echo "[SKIP] missing prmtop ${prm}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${sys}" "0/${EQUIL_SEGS}" "0/${PROD_REPS}" "0/$((PROD_REPS * PROD_SEGS_PER_REP))" "SKIP_NO_PRM" "0" >>"${BASEDIR}/logs/resilient_queue_summary.tsv"
    return 1
  }
  ensure_prmtop_link "${sys}" || return 1
  prm_short="${PQUEUE}/${sys}.prmtop"

  local crashes_recovered=0 eq_done=0 pr_rep_done=0 pr_seg_done=0 status=OK
  local es seg_name segdir seg_rst crd_in r prev_pad
  local final_equil_rst="${sim}/equil_seg10/equil_seg10.rst7"

  prev_rst=""
  for es in $(seq -w 1 "${EQUIL_SEGS}"); do
    seg_name="equil_seg${es}"
    segdir="${sim}/${seg_name}"
    seg_rst="${segdir}/${seg_name}.rst7"
    if rst7_skip_ok "${prm}" "${seg_rst}"; then
      echo "[SKIP] ${sys} ${seg_name} (valid rst7)"
      prev_rst="${seg_rst}"
      eq_done=$((eq_done + 1))
      continue
    fi
    mkdir -p "${segdir}"
    if [[ $((10#${es})) -eq 1 ]]; then
      crd_in="$(find_starting_equil_rst "${sim}")" || true
      if [[ -z "${crd_in}" ]]; then
        echo "[FAIL] ${sys}: no equil3 / gpu_bridge2 / equil2 / equil1.rst7 for ${seg_name}"
        status=FAIL
        break
      fi
    else
      prev_pad=$(printf '%02d' $((10#${es} - 1)))
      crd_in="${sim}/equil_seg${prev_pad}/equil_seg${prev_pad}.rst7"
    fi

    while true; do
      set +e
      run_segment "${sys}" "${seg_name}" "${segdir}" "${seg_name}" "${crd_in}" "${disang}" "${prm_short}"
      r=$?
      set -e
      if [[ "${r}" -eq 0 ]]; then
        eq_done=$((eq_done + 1))
        prev_rst="${seg_rst}"
        break
      fi
      if [[ "${r}" -eq 2 ]]; then
        if [[ $((10#${es})) -lt 2 ]]; then
          echo "[FAIL] ${sys} ${seg_name} cannot rewind before first equil segment"
          status=FAIL
          break 2
        fi
        prev_pad=$(printf '%02d' $((10#${es} - 1)))
        crd_in="${sim}/equil_seg${prev_pad}/equil_seg${prev_pad}.rst7"
        crashes_recovered=$((crashes_recovered + 1))
        echo "[RECHAIN] ${sys} ${seg_name} crd_in=${crd_in}"
        rm -f "${segdir}/${seg_name}.out" "${seg_rst}" "${segdir}/${seg_name}.nc" 2>/dev/null || true
        continue
      fi
      status=FAIL
      break 2
    done
    [[ "${status}" == "OK" ]] || break
  done

  if [[ "${status}" != "OK" ]]; then
    echo "[ABORT] ${sys} equilibration incomplete"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${sys}" "${eq_done}/${EQUIL_SEGS}" "${pr_rep_done}/${PROD_REPS}" "${pr_seg_done}/$((PROD_REPS * PROD_SEGS_PER_REP))" "${status}" "${crashes_recovered}" >>"${BASEDIR}/logs/resilient_queue_summary.tsv"
    return 1
  fi

  if [[ ! -s "${final_equil_rst}" ]]; then
    echo "[FAIL] ${sys} missing ${final_equil_rst}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${sys}" "${eq_done}/${EQUIL_SEGS}" "${pr_rep_done}/${PROD_REPS}" "${pr_seg_done}/$((PROD_REPS * PROD_SEGS_PER_REP))" "FAIL" "${crashes_recovered}" >>"${BASEDIR}/logs/resilient_queue_summary.tsv"
    return 1
  fi

  local rep crd_prod s01 s02 d01 d02 r01 r02

  for rep in $(seq -w 1 "${PROD_REPS}"); do
    s01="prod_rep${rep}_seg01"
    s02="prod_rep${rep}_seg02"
    d01="${sim}/${s01}"
    d02="${sim}/${s02}"
    r01="${d01}/${s01}.rst7"
    r02="${d02}/${s02}.rst7"

    while true; do
      local rep_restart=0

      if rst7_skip_ok "${prm}" "${r01}"; then
        echo "[SKIP] ${sys} ${s01} (valid rst7)"
        pr_seg_done=$((pr_seg_done + 1))
      else
        mkdir -p "${d01}"
        crd_prod="${final_equil_rst}"
        while true; do
          set +e
          run_segment "${sys}" "${s01}" "${d01}" "${s01}" "${crd_prod}" "${disang}" "${prm_short}"
          r=$?
          set -e
          if [[ "${r}" -eq 0 ]]; then
            pr_seg_done=$((pr_seg_done + 1))
            break
          fi
          if [[ "${r}" -eq 2 ]]; then
            crashes_recovered=$((crashes_recovered + 1))
            echo "[RECHAIN] ${sys} ${s01} reset to ${final_equil_rst}"
            rm -f "${d01}/${s01}.out" "${r01}" "${d01}/${s01}.nc" 2>/dev/null || true
            crd_prod="${final_equil_rst}"
            continue
          fi
          status=FAIL
          break 3
        done
      fi
      [[ "${status}" == "OK" ]] || break

      if rst7_skip_ok "${prm}" "${r02}"; then
        echo "[SKIP] ${sys} ${s02} (valid rst7)"
        pr_seg_done=$((pr_seg_done + 1))
      else
        mkdir -p "${d02}"
        crd_prod="${r01}"
        while true; do
          set +e
          run_segment "${sys}" "${s02}" "${d02}" "${s02}" "${crd_prod}" "${disang}" "${prm_short}"
          r=$?
          set -e
          if [[ "${r}" -eq 0 ]]; then
            pr_seg_done=$((pr_seg_done + 1))
            break
          fi
          if [[ "${r}" -eq 2 ]]; then
            crashes_recovered=$((crashes_recovered + 1))
            echo "[RECHAIN] ${sys} ${s02} restart full replicate from ${final_equil_rst}"
            rm -f "${d01}/${s01}.out" "${r01}" "${d01}/${s01}.nc" "${d02}/${s02}.out" "${r02}" "${d02}/${s02}.nc" 2>/dev/null || true
            # Undo seg01 count (skip or completed) before re-running the pair
            pr_seg_done=$((pr_seg_done - 1))
            rep_restart=1
            break
          fi
          status=FAIL
          break 3
        done
      fi

      if [[ "${rep_restart}" -eq 1 ]]; then
        continue
      fi
      pr_rep_done=$((pr_rep_done + 1))
      break
    done
    [[ "${status}" == "OK" ]] || break
  done

  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${sys}" "${eq_done}/${EQUIL_SEGS}" "${pr_rep_done}/${PROD_REPS}" "${pr_seg_done}/$((PROD_REPS * PROD_SEGS_PER_REP))" "${status}" "${crashes_recovered}" >>"${BASEDIR}/logs/resilient_queue_summary.tsv"
  [[ "${status}" == "OK" ]]
}

: >"${BASEDIR}/logs/resilient_queue_summary.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "system" "equil_segs_done" "prod_reps_done" "prod_segs_done" "status" "crashes_recovered" | tee -a "${BASEDIR}/logs/resilient_queue_summary.tsv"

for sys in "${QUEUE[@]}"; do
  SYSLOG="${BASEDIR}/logs/${sys}_resilient_${TS}.log"
  set +e
  run_one_system "${sys}" 2>&1 | tee -a "${SYSLOG}"
  rsys=$?
  set -e
  if [[ "${rsys}" -ne 0 ]]; then
    echo "[QUEUE] non-zero exit for ${sys} (continuing queue)"
  fi
done

echo ""
echo "=== SUMMARY (logs/resilient_queue_summary.tsv) ==="
column -t -s $'\t' "${BASEDIR}/logs/resilient_queue_summary.tsv" 2>/dev/null || cat "${BASEDIR}/logs/resilient_queue_summary.tsv"
echo "=== DONE main log: ${MAIN_LOG} ==="
echo "READY TO LAUNCH — queue not started by this script."
echo "Launch: nohup bash scripts/run_resilient_queue.sh > logs/resilient_queue.log 2>&1 &"
