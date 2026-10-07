#!/bin/bash
# scripts/run_production.sh
# -------------------------------------------------------
# Production MD: N independent replicates from the same equilibrated restart.
#
# USAGE:
#   bash scripts/run_production.sh <sim_name> [n_replicates]
#
# Arguments:
#   sim_name      Basename only, e.g. mol6_pose2_monomer
#   n_replicates  Optional; default from config NREP (typically 3)
#
# Each replicate: 20 ns (single segment; production.in dt=0.002, nstlim=10000000), same
#   -p systems/<name>/solvated_hmr.prmtop
#   -c simulations/<name>/equil_free/equil_free.rst7
# Output per replicate under systems/<name>/prod_repNN/ (prod.rst7, prod.nc, prod.out).
# ig=-1 in production.in gives a different seed each run → independent trajectories.
# Replicates run sequentially on one GPU.
#
# Legacy simulations/<name>/prod/ with prod_01.nc is not used by this script.
# -------------------------------------------------------
set -euo pipefail
_SCRIPT_ROOT="$(cd "$(dirname "$0")" && pwd)"
source "${_SCRIPT_ROOT}/../config.env"
export PATH="/usr/bin:${AMBER_HOME}/bin:/usr/local/bin:${PATH}"
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/openmpi/lib:${LD_LIBRARY_PATH:-}"
export CUDA_DEVICE_ORDER="${CUDA_DEVICE_ORDER:-PCI_BUS_ID}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-1}"
PMEMD_CUDA_BIN="${PMEMD_CUDA_BIN:-/opt/software/AMBER/pmemd24_cuda128/bin/pmemd.cuda}"
[[ -x "${PMEMD_CUDA_BIN}" ]] || {
  echo "[FAIL] Not executable: ${PMEMD_CUDA_BIN}" >&2
  exit 1
}
# shellcheck source=scripts/amber_protocol_gates.sh
source "${_SCRIPT_ROOT}/amber_protocol_gates.sh"

NAME=${1:?"Usage: bash run_production.sh <name> [n_replicates]"}
NREP=${2:-${NREP:-3}}

SDIR="${PROJECT_ROOT}/systems/${NAME}"
INP="${PROJECT_ROOT}/simulations/input_files/production.in"
START_RST="${PROJECT_ROOT}/simulations/${NAME}/equil_free/equil_free.rst7"

# Prefer the HMR topology produced by the builder/equilibration step.
if [[ -f "${SDIR}/solvated_hmr.prmtop" ]]; then
    PARM="${SDIR}/solvated_hmr.prmtop"
elif [[ -f "${SDIR}/solvated_hmr.parm7" ]]; then
    PARM="${SDIR}/solvated_hmr.parm7"
else
    PARM="${SDIR}/solvated.parm7"
fi

gate_production_protocol "${NAME}" "${SDIR}" "${INP}" "${PARM}"

PRMTOP="${PARM}"
NS_PER_REP="${PROD_NS_PER_REP:-20}"

[[ -f "${START_RST}" ]] || {
  echo "[FAIL] Missing equilibrated restart: ${START_RST}" >&2
  exit 1
}

log() { echo "[$(date '+%H:%M:%S')] [$NAME] $1"; }

echo "Production: ${NREP} replicates × ${NS_PER_REP} ns = $((NREP * NS_PER_REP)) ns total"
echo "Topology:   ${PRMTOP}"
echo "Starting from: ${START_RST}"
log "=== Production: ${NREP} × ${NS_PER_REP} ns (sequential replicates) ==="

for REP in $(seq -w 1 "${NREP}"); do
    REP_DIR="${SDIR}/prod_rep${REP}"
    mkdir -p "${REP_DIR}"
    OUT="${REP_DIR}/prod"

    if [[ -f "${OUT}.rst7" ]] && [[ -f "${OUT}.nc" ]]; then
        log "Rep ${REP}: ${OUT}.rst7 + ${OUT}.nc exist — skipping"
        continue
    fi

    log "Rep ${REP}: starting (output ${REP_DIR}/)"

    if ! "${PMEMD_CUDA_BIN}" -O \
        -i "${INP}" \
        -o "${OUT}.out" \
        -p "${PARM}" \
        -c "${START_RST}" \
        -r "${OUT}.rst7" \
        -x "${OUT}.nc" \
        -e "${OUT}.mden" \
        -ref "${START_RST}" \
        -inf "${OUT}.mdinfo"
    then
        log "[FAIL] Replicate ${REP}"
        exit 1
    fi

    NS_DAY=$(grep "ns/day" "${OUT}.out" | tail -1 | awk '{print $4}')
    log "[OK] Rep ${REP}: ${NS_DAY} ns/day"

    PROD_REI="${REP_DIR}/${PROD_TRAJ_REIMAGED:-production_reimaged.nc}"
    log "Post-processing rep ${REP}: autoimage → ${PROD_REI}"
    bash "${_SCRIPT_ROOT}/preprocess_trajectory.sh" "${PARM}" "${PROD_REI}" "${OUT}.nc"
done

log "=== Production COMPLETE: $((NREP * NS_PER_REP)) ns over ${NREP} replicates ==="
