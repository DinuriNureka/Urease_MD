#!/bin/bash
# scripts/run_equilibration.sh
# -------------------------------------------------------
# Sequential equilibration for one system.
#
# Usage:
#     bash scripts/run_equilibration.sh mol01_pose01
#
# Runs: MIN1 -> MIN2 -> HEAT_NVT_CPU -> HEAT -> EQUIL_CPU_NPT (CPU pmemd) -> EQUIL1 -> EQUIL2 -> EQUIL3 -> EQUIL_FREE
# Minimization: sander. CPU NPT density step: pmemd (ntt=3). GPU: pmemd.cuda (ntt=11). CUDA_DEVICE_ORDER/CUDA_VISIBLE_DEVICES from config.env.
# -------------------------------------------------------
set -euo pipefail
_SCRIPT_ROOT="$(cd "$(dirname "$0")" && pwd)"
source "${_SCRIPT_ROOT}/../config.env"
export PATH="/usr/bin:${AMBER_HOME}/bin:/usr/local/bin:${PATH}"
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/openmpi/lib:${LD_LIBRARY_PATH:-}"
export CUDA_DEVICE_ORDER="${CUDA_DEVICE_ORDER:-PCI_BUS_ID}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-1}"
# shellcheck source=scripts/amber_protocol_gates.sh
source "${_SCRIPT_ROOT}/amber_protocol_gates.sh"

NAME=${1:?"Usage: bash run_equilibration.sh <system_name>"}
SDIR="${PROJECT_ROOT}/systems/${NAME}"
SIM="${PROJECT_ROOT}/simulations/${NAME}"
INP="${PROJECT_ROOT}/simulations/input_files"

INHIB="$(echo "${NAME}" | sed -E 's/_pose[0-9]*//; s/_monomer$//')"

# === PARAMETER SAFETY GATE ===
echo "[GATE] Verifying system parameters before MD launch..."
if ! bash "${_SCRIPT_ROOT}/verify_system_ready.sh" "${SDIR}"; then
  echo "[ABORT] System failed safety verification. Fix issues above."
  echo "[ABORT] Will not launch MD on an unverified system."
  exit 1
fi
echo "[GATE] Verification passed. Proceeding with MD."
# ==============================

PMEMD_GPU_BIN="${PMEMD_CUDA_BIN:-/opt/software/AMBER/pmemd24_cuda128/bin/pmemd.cuda}"
SANDER_BIN="${SANDER_BIN:-${AMBER_HOME}/bin/sander}"
[[ -x "${SANDER_BIN}" ]] || SANDER_BIN="/usr/local/bin/sander"
PMEMD_CPU_BIN="${AMBER_HOME}/bin/pmemd"
[[ -x "${PMEMD_CPU_BIN}" ]] || PMEMD_CPU_BIN="${SANDER_BIN}"

# Topology/coordinates selection:
# Require `solvated_hmr.prmtop` + `solvated_hmr.rst7` and use them for ALL stages.
# This is the single authoritative topology used across the full equilibration pipeline.
PARM="${SDIR}/solvated_hmr.prmtop"
RST="${SDIR}/solvated_hmr.rst7"
[[ -f "${PARM}" ]] || { echo "[FAIL] Missing required topology: ${PARM}" >&2; exit 1; }
[[ -f "${RST}"  ]] || { echo "[FAIL] Missing required coordinates: ${RST}" >&2; exit 1; }

# Protocol gates (exit 1 if topology or mdin do not match project rules)
gate_equilibration_md_inputs "${INP}"

log() { echo "[$(date '+%H:%M:%S')] [$NAME] $1"; }

# Last instantaneous density (g/cm³) from pmemd/sander mdout NSTEP blocks.
read_last_equil_cpu_density() {
  local mdout=$1
  grep "Density" "${mdout}" | grep -v "^|" | awk '{print $3}' | grep -E '^[0-9]' | tail -1
}

# After equil_cpu_npt: inform only; do not abort (GPU may still work if slightly low).
warn_equil_cpu_npt_density() {
  local mdout=$1
  [[ -f "${mdout}" ]] || {
    log "[WARN] Missing mdout for density check: ${mdout}"
    return 0
  }
  local D
  D=$(read_last_equil_cpu_density "${mdout}")
  [[ -n "${D}" ]] || {
    log "[WARN] Could not parse Density from ${mdout}"
    return 0
  }
  export GATE_D="${D}"
  python3 - << 'PYEOF'
import os
d = float(os.environ["GATE_D"])
if d < 0.98:
    print(f"[WARN] Density after equil_cpu_npt: {d:.4f} g/cm³ (< 0.98) — continuing to GPU; monitor equil1")
else:
    print(f"[OK] Density after equil_cpu_npt: {d:.4f} g/cm³ (>= 0.98)")
PYEOF
}

HEAT_RST="${SIM}/heat/heat.rst7"
EQUIL_NPT_RST="${SIM}/equil_cpu_npt/equil_cpu_npt.rst7"
EQUIL_NPT_MDOUT="${SIM}/equil_cpu_npt/equil_cpu_npt.out"
GPU_BRIDGE_RST="${SIM}/gpu_bridge/gpu_bridge.rst7"
GPU_BRIDGE2_RST="${SIM}/gpu_bridge2/gpu_bridge2.rst7"
EQUIL1_RST="${SIM}/equil1/equil1.rst7"
EQUIL2_RST="${SIM}/equil2/equil2.rst7"
EQUIL3_RST="${SIM}/equil3/equil3.rst7"
EQUIL_FREE_RST="${SIM}/equil_free/equil_free.rst7"

# Resume: pick earliest incomplete stage (requires equil_cpu_npt before GPU bridge).
RESUME_MODE="full"
if [[ "${EQUIL_FORCE_FULL:-0}" != "1" ]]; then
  if [[ -f "${HEAT_RST}" ]] && [[ ! -s "${EQUIL_NPT_RST}" ]]; then
    log "heat.rst7 present but no equil_cpu_npt.rst7 — running equil_cpu_npt only, then downstream stages"
    RESUME_MODE="from_equil_cpu_npt_only"
  elif [[ -s "${EQUIL_NPT_RST}" ]] && [[ ! -f "${GPU_BRIDGE_RST}" ]]; then
    log "equil_cpu_npt done but no gpu_bridge.rst7 — resuming from gpu_bridge"
    RESUME_MODE="from_gpu_bridge"
  elif [[ -f "${GPU_BRIDGE_RST}" ]] && [[ ! -f "${EQUIL1_RST}" ]]; then
    log "gpu_bridge done but no equil1.rst7 — resuming from equil1 (GPU)"
    RESUME_MODE="from_equil1"
  elif [[ -f "${EQUIL1_RST}" ]] && [[ ! -f "${EQUIL2_RST}" ]]; then
    RESUME_MODE="from_equil2"
  elif [[ -f "${EQUIL2_RST}" ]] && [[ ! -f "${EQUIL3_RST}" ]]; then
    RESUME_MODE="from_equil3"
  elif [[ -f "${EQUIL3_RST}" ]] && [[ ! -f "${GPU_BRIDGE2_RST}" ]]; then
    log "equil3 done but no gpu_bridge2.rst7 — resuming from gpu_bridge2"
    RESUME_MODE="from_gpu_bridge2"
  elif [[ -f "${GPU_BRIDGE2_RST}" ]] && [[ ! -f "${EQUIL_FREE_RST}" ]]; then
    log "gpu_bridge2 done but no equil_free.rst7 — resuming from equil_free"
    RESUME_MODE="from_equil_free_only"
  fi
fi

# kind: min = sander; cpu = pmemd; gpu = pmemd.cuda; cpu_mpi = pmemd.MPI + mpirun
run_stage() {
    local stage=$1
    local input=$2
    local prev_rst=$3
    local kind=$4

    local dir="${SIM}/${stage}"
    mkdir -p "${dir}"

    case "${kind}" in
      min)
        [[ -x "${SANDER_BIN}" ]] || { echo "[FAIL] Missing executable: ${SANDER_BIN}" >&2; exit 1; }
        log "--- ${stage} (min: $(basename "${SANDER_BIN}")) ---"
        "${SANDER_BIN}" -O \
            -i "${input}" \
            -o "${dir}/${stage}.out" \
            -p "${PARM}" \
            -c "${prev_rst}" \
            -r "${dir}/${stage}.rst7" \
            -ref "${prev_rst}" \
            -inf "${dir}/${stage}.mdinfo"
        ;;
      cpu)
        [[ -x "${PMEMD_CPU_BIN}" ]] || { echo "[FAIL] Missing executable: ${PMEMD_CPU_BIN}" >&2; exit 1; }
        log "--- ${stage} (cpu: $(basename "${PMEMD_CPU_BIN}")) ---"
        "${PMEMD_CPU_BIN}" -O \
            -i "${input}" \
            -o "${dir}/${stage}.out" \
            -p "${PARM}" \
            -c "${prev_rst}" \
            -r "${dir}/${stage}.rst7" \
            -x "${dir}/${stage}.nc" \
            -ref "${prev_rst}" \
            -inf "${dir}/${stage}.mdinfo"
        ;;
      gpu)
        export CUDA_DEVICE_ORDER="${CUDA_DEVICE_ORDER:-PCI_BUS_ID}"
        export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-1}"
        [[ -x "${PMEMD_GPU_BIN}" ]] || { echo "[FAIL] Missing executable: ${PMEMD_GPU_BIN}" >&2; exit 1; }
        log "--- ${stage} (gpu: $(basename "${PMEMD_GPU_BIN}")) ---"
        "${PMEMD_GPU_BIN}" -O \
            -i "${input}" \
            -o "${dir}/${stage}.out" \
            -p "${PARM}" \
            -c "${prev_rst}" \
            -r "${dir}/${stage}.rst7" \
            -x "${dir}/${stage}.nc" \
            -ref "${prev_rst}" \
            -inf "${dir}/${stage}.mdinfo"
        ;;
      cpu_mpi)
        local PMEMD_MPI="${AMBER_HOME}/bin/pmemd.MPI"
        [[ -x "${PMEMD_MPI}" ]] || { echo "[FAIL] Missing ${PMEMD_MPI}" >&2; exit 1; }
        command -v mpirun &>/dev/null || { echo "[FAIL] mpirun not in PATH" >&2; exit 1; }
        local ORCA_RUNNING
        ORCA_RUNNING="$(pgrep -c orca 2>/dev/null || true)"
        [[ -z "${ORCA_RUNNING}" ]] && ORCA_RUNNING=0
        local TC="${TOTAL_CORES:-${N_CPU_CORES:-56}}"
        local NPT_CORES
        if [[ "${ORCA_RUNNING}" -gt 0 ]]; then
            NPT_CORES="${MPI_CORES_NPT:-48}"
        else
            NPT_CORES=$((TC - 2))
        fi
        [[ "${NPT_CORES}" -lt 1 ]] && NPT_CORES=1
        log "--- ${stage} (cpu_mpi: ${NPT_CORES} ranks, $(basename "${PMEMD_MPI}"), pgrep orca=${ORCA_RUNNING}) ---"
        mpirun --oversubscribe -np "${NPT_CORES}" "${PMEMD_MPI}" -O \
            -i "${input}" \
            -o "${dir}/${stage}.out" \
            -p "${PARM}" \
            -c "${prev_rst}" \
            -r "${dir}/${stage}.rst7" \
            -x "${dir}/${stage}.nc" \
            -ref "${prev_rst}" \
            -inf "${dir}/${stage}.mdinfo"
        ;;
      *) echo "[FAIL] Unknown engine kind: ${kind}" >&2; exit 1 ;;
    esac

    log "[OK] ${stage}"
}

# CPU NPT density convergence — use MPI pmemd when available (see config.env ORCA_CORES / MPI_CORES_NPT).
run_equil_cpu_npt() {
    local stage=equil_cpu_npt
    local input="${INP}/equil_cpu_npt.in"
    local prev_rst="${HEAT_RST}"
    local dir="${SIM}/${stage}"
    mkdir -p "${dir}"

    local PMEMD_MPI="${AMBER_HOME}/bin/pmemd.MPI"
    local TC="${TOTAL_CORES:-${N_CPU_CORES:-56}}"
    local ORCA_RUNNING
    ORCA_RUNNING="$(pgrep -c orca 2>/dev/null || true)"
    [[ -z "${ORCA_RUNNING}" ]] && ORCA_RUNNING=0

    if [[ -x "${PMEMD_MPI}" ]] && command -v mpirun &>/dev/null; then
        local NPT_CORES
        if [[ "${ORCA_RUNNING}" -gt 0 ]]; then
            NPT_CORES="${MPI_CORES_NPT:-48}"
        else
            NPT_CORES=$((TC - 2))
        fi
        [[ "${NPT_CORES}" -lt 1 ]] && NPT_CORES=1
        log "--- ${stage} (cpu: mpi ${NPT_CORES} ranks, $(basename "${PMEMD_MPI}")) ---"
        log "ORCA processes (pgrep -c orca): ${ORCA_RUNNING} — NPT cores=${NPT_CORES}"
        mpirun --oversubscribe -np "${NPT_CORES}" "${PMEMD_MPI}" -O \
            -i "${input}" \
            -o "${dir}/${stage}.out" \
            -p "${PARM}" \
            -c "${prev_rst}" \
            -r "${dir}/${stage}.rst7" \
            -x "${dir}/${stage}.nc" \
            -ref "${prev_rst}" \
            -inf "${dir}/${stage}.mdinfo"
    else
        log "pmemd.MPI not found or mpirun missing — single-core CPU NPT ($(basename "${PMEMD_CPU_BIN}"))"
        [[ -x "${PMEMD_CPU_BIN}" ]] || { echo "[FAIL] Missing executable: ${PMEMD_CPU_BIN}" >&2; exit 1; }
        log "--- ${stage} (cpu: $(basename "${PMEMD_CPU_BIN}")) ---"
        "${PMEMD_CPU_BIN}" -O \
            -i "${input}" \
            -o "${dir}/${stage}.out" \
            -p "${PARM}" \
            -c "${prev_rst}" \
            -r "${dir}/${stage}.rst7" \
            -x "${dir}/${stage}.nc" \
            -ref "${prev_rst}" \
            -inf "${dir}/${stage}.mdinfo"
    fi

    log "[OK] ${stage}"
}

log "=== Starting equilibration (resume mode: ${RESUME_MODE}) ==="

# Post–CPU-NPT GPU path: bridge → equil1–3 (GPU) → bridge2 → equil_free (GPU).
run_post_npt_stages() {
    run_stage gpu_bridge  "${INP}/gpu_bridge.in"   "${EQUIL_NPT_RST}"              gpu
    run_stage equil1      "${INP}/equil1.in"        "${GPU_BRIDGE_RST}"             gpu
    run_stage equil2      "${INP}/equil2.in"        "${EQUIL1_RST}"                 gpu
    run_stage equil3      "${INP}/equil3.in"        "${EQUIL2_RST}"                 gpu
    run_stage gpu_bridge2 "${INP}/gpu_bridge.in"    "${EQUIL3_RST}"                 gpu
    run_stage equil_free  "${INP}/equil_free.in"    "${GPU_BRIDGE2_RST}"            gpu
}

case "${RESUME_MODE}" in
  full)
    # Default: always run min1 from ${RST}. Set EQUIL_SKIP_MIN1=1 to skip min1 when
    # ${SIM}/min1/min1.rst7 already exists (e.g. good Ni geometry; restart from min2 only).
    if [[ "${EQUIL_SKIP_MIN1:-0}" != "1" ]]; then
      run_stage min1           "${INP}/min1.in"           "${RST}"                    gpu
    else
      [[ -s "${SIM}/min1/min1.rst7" ]] || {
        echo "[FAIL] EQUIL_SKIP_MIN1=1 requires existing ${SIM}/min1/min1.rst7" >&2
        exit 1
      }
      log "--- min1 SKIPPED (EQUIL_SKIP_MIN1=1), using ${SIM}/min1/min1.rst7 ---"
    fi
    run_stage min2           "${INP}/min2.in"           "${SIM}/min1/min1.rst7"     gpu
    run_stage heat_nvt_cpu   "${INP}/heat_nvt_cpu.in"   "${SIM}/min2/min2.rst7"     cpu_mpi
    run_stage heat           "${INP}/heat.in"           "${SIM}/heat_nvt_cpu/heat_nvt_cpu.rst7" gpu
    run_equil_cpu_npt
    warn_equil_cpu_npt_density "${EQUIL_NPT_MDOUT}"
    run_post_npt_stages
    ;;
  from_equil_cpu_npt_only)
    run_equil_cpu_npt
    warn_equil_cpu_npt_density "${EQUIL_NPT_MDOUT}"
    run_post_npt_stages
    ;;
  from_gpu_bridge)
    [[ -f "${EQUIL_NPT_MDOUT}" ]] && warn_equil_cpu_npt_density "${EQUIL_NPT_MDOUT}"
    run_post_npt_stages
    ;;
  from_equil1)
    run_stage equil1      "${INP}/equil1.in"        "${GPU_BRIDGE_RST}"             gpu
    run_stage equil2      "${INP}/equil2.in"        "${EQUIL1_RST}"                 gpu
    run_stage equil3      "${INP}/equil3.in"        "${EQUIL2_RST}"                 gpu
    run_stage gpu_bridge2 "${INP}/gpu_bridge.in"    "${EQUIL3_RST}"                 gpu
    run_stage equil_free  "${INP}/equil_free.in"    "${GPU_BRIDGE2_RST}"            gpu
    ;;
  from_equil2)
    run_stage equil2      "${INP}/equil2.in"        "${EQUIL1_RST}"                 gpu
    run_stage equil3      "${INP}/equil3.in"        "${EQUIL2_RST}"                 gpu
    run_stage gpu_bridge2 "${INP}/gpu_bridge.in"    "${EQUIL3_RST}"                 gpu
    run_stage equil_free  "${INP}/equil_free.in"    "${GPU_BRIDGE2_RST}"            gpu
    ;;
  from_equil3)
    run_stage equil3      "${INP}/equil3.in"        "${EQUIL2_RST}"                 gpu
    run_stage gpu_bridge2 "${INP}/gpu_bridge.in"    "${EQUIL3_RST}"                 gpu
    run_stage equil_free  "${INP}/equil_free.in"    "${GPU_BRIDGE2_RST}"            gpu
    ;;
  from_gpu_bridge2)
    run_stage gpu_bridge2 "${INP}/gpu_bridge.in"    "${EQUIL3_RST}"                 gpu
    run_stage equil_free  "${INP}/equil_free.in"    "${GPU_BRIDGE2_RST}"            gpu
    ;;
  from_equil_free_only)
    run_stage equil_free  "${INP}/equil_free.in"    "${GPU_BRIDGE2_RST}"            gpu
    ;;
  *)
    echo "[FAIL] Unknown RESUME_MODE: ${RESUME_MODE}" >&2
    exit 1
    ;;
esac

log "=== Equilibration COMPLETE ==="
log "Final restart: ${SIM}/equil_free/equil_free.rst7"

# --- Post-equilibration gates (order: Ni geometry → convergence → production protocol → launch) ---
EQUIL_FREE_DIR="${SIM}/equil_free"
export EQUIL_SYSTEM_NAME="${NAME}"

"${_SCRIPT_ROOT}/validate_ni_geometry.sh" "${NAME}" "${PARM}" "${EQUIL_FREE_DIR}" || exit 1

gate_equil_convergence "${SDIR}" "${INHIB}" "${PARM}" "${EQUIL_FREE_DIR}" "${NAME}" || exit 1

PROD_PARM="${SDIR}/solvated_hmr.prmtop"
[[ -f "${PROD_PARM}" ]] || { echo "[FAIL] Missing required topology for production: ${PROD_PARM}" >&2; exit 1; }
gate_production_protocol "${NAME}" "${SDIR}" "${INP}/production.in" "${PROD_PARM}" || exit 1

if [[ "${LAUNCH_PRODUCTION_AFTER_EQUIL:-1}" == "1" ]]; then
  echo "All gates passed. Launching production with HMR topology."
  "${_SCRIPT_ROOT}/run_production.sh" "${NAME}"
else
  log "LAUNCH_PRODUCTION_AFTER_EQUIL=0 — skipping run_production.sh (gates passed)."
fi
