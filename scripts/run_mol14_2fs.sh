#!/usr/bin/env bash
# scripts/run_mol14_2fs.sh
# -------------------------------------------------------
# Run mol14_pose2_monomer at 2 fs WITHOUT HMR — supervisor validation test.
#
# Protocol:
#   Min1 (sander) → Min2 (sander) →
#   Heat (pmemd.MPI 54 ranks, heat_nvt_cpu.in, 1 ps NVT 0→50 K) →
#   Equil (pmemd.cuda, equil_free_2fs.in, 100 ns NPT) →
#   Production 3 × 20 ns (pmemd.cuda, production_2fs_single.in)
#
# Output: simulations/mol14_pose2_monomer_2fs/
#   min1/  min2/  heat/  equil/  prod_rep1/  prod_rep2/  prod_rep3/
#
# TOPOLOGY : simulations/mol14_pose2_monomer/solvated.prmtop  (NO HMR)
# COORDS   : simulations/mol14_pose2_monomer/solvated.inpcrd
#
# RESUMABLE: Each stage checks for its output rst7 before running.
#            If equil/equil.rst7 exists, all equilibration is skipped.
#
# Usage:
#   nohup bash scripts/run_mol14_2fs.sh > logs/mol14_2fs_run.log 2>&1 &
#   echo "PID: $!"
# -------------------------------------------------------
set -euo pipefail

_SCRIPT_ROOT="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${_SCRIPT_ROOT}/.." && pwd)"
source "${ROOT}/config.env"

# ── Configuration ────────────────────────────────────────────────────────────
export CUDA_VISIBLE_DEVICES=1
export CUDA_DEVICE_ORDER="${CUDA_DEVICE_ORDER:-PCI_BUS_ID}"
export PATH="${AMBER_HOME}/bin:/usr/local/bin:${PATH}"

SYS="mol14_pose2_monomer"
PRMTOP="${ROOT}/simulations/${SYS}/solvated.prmtop"
INPCRD="${ROOT}/simulations/${SYS}/solvated.inpcrd"
OUTDIR="${ROOT}/simulations/mol14_pose2_monomer_2fs"
INP="${ROOT}/simulations/input_files"
DISANG="${ROOT}/ni_ni.RST"

MPI_RANKS="${MPI_RANKS_HEAT:-54}"
PMEMD="${PMEMD_CUDA_BIN:-/opt/software/AMBER/pmemd24_cuda128/bin/pmemd.cuda}"
PMEMD_MPI="${AMBER_HOME}/bin/pmemd.MPI"
SANDER="${SANDER_BIN:-${AMBER_HOME}/bin/sander}"
[[ -x "${SANDER}" ]] || SANDER="/usr/local/bin/sander"

# ── Internal log setup ───────────────────────────────────────────────────────
mkdir -p "${ROOT}/logs"
LOGFILE="${ROOT}/logs/mol14_2fs_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOGFILE}") 2>&1

log() { echo "[$(date '+%H:%M:%S')] [mol14_2fs] $1"; }
nan_check() {
  local outfile="$1" stage="$2"
  local n
  n=$(grep -ac "NaN" "${outfile}" 2>/dev/null || echo 0)
  [[ "${n}" -gt 0 ]] && log "[WARN] NaN detected: ${n} occurrences in ${stage} output (${outfile})"
}
nsday() {
  local mdinfo="$1"
  grep -a "ns/day\|Elapsed\|Total wall time" "${mdinfo}" 2>/dev/null | tail -5 \
    || log "  (performance info not available in ${mdinfo})"
}

# ── Pre-flight checks ────────────────────────────────────────────────────────
log "=== mol14_pose2_monomer 2fs no-HMR test pipeline ==="
log "Topology : ${PRMTOP}"
log "Coords   : ${INPCRD}"
log "Output   : ${OUTDIR}"
log "GPU      : CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"
log "GPU bin  : ${PMEMD}"
log "MPI bin  : ${PMEMD_MPI} (${MPI_RANKS} ranks for heat)"
log "Log      : ${LOGFILE}"

[[ -f "${PRMTOP}" ]] || {
  log "[FAIL] Non-HMR topology not found: ${PRMTOP}"
  log "       Run first: bash scripts/build_nonhmr_prmtop.sh ${SYS}"
  exit 1
}
[[ -f "${INPCRD}" ]] || {
  log "[FAIL] Starting coordinates not found: ${INPCRD}"
  log "       Run first: bash scripts/build_nonhmr_prmtop.sh ${SYS}"
  exit 1
}
[[ -x "${PMEMD}" ]]     || { log "[FAIL] Not executable: ${PMEMD}"     >&2; exit 1; }
[[ -x "${SANDER}" ]]    || { log "[FAIL] sander not found"              >&2; exit 1; }
[[ -x "${PMEMD_MPI}" ]] || { log "[FAIL] pmemd.MPI not found: ${PMEMD_MPI}" >&2; exit 1; }
command -v mpirun &>/dev/null || { log "[FAIL] mpirun not in PATH" >&2; exit 1; }
[[ -L "${DISANG}" || -f "${DISANG}" ]] || { log "[FAIL] ni_ni.RST not found: ${DISANG}" >&2; exit 1; }
for f in min1.in min2.in heat_nvt_cpu.in equil_free_2fs.in production_2fs_single.in; do
  [[ -f "${INP}/${f}" ]] || { log "[FAIL] Missing input file: ${INP}/${f}" >&2; exit 1; }
done

# Create output subdirectories.
mkdir -p "${OUTDIR}/min1" "${OUTDIR}/min2" "${OUTDIR}/heat" "${OUTDIR}/equil" \
         "${OUTDIR}/prod_rep1" "${OUTDIR}/prod_rep2" "${OUTDIR}/prod_rep3"

# ── Stage helpers ─────────────────────────────────────────────────────────────
run_min() {
  # Usage: run_min <label> <mdin> <rst_in> <ref> <outdir>
  local label="$1" mdin="$2" rst_in="$3" ref="$4" odir="$5"
  local rst_out="${odir}/${label}.rst7"
  local out="${odir}/${label}.out"

  if [[ -f "${rst_out}" ]]; then
    log "[SKIP] ${label} — ${rst_out} exists"
    return 0
  fi
  log "--- ${label} (sander minimization) ---"
  "${SANDER}" -O \
    -i  "${mdin}" \
    -o  "${out}" \
    -p  "${PRMTOP}" \
    -c  "${rst_in}" \
    -r  "${rst_out}" \
    -ref "${ref}" \
    -x  "${odir}/${label}.nc" \
    -inf "${odir}/${label}.mdinfo"
  log "${label} done."
  nan_check "${out}" "${label}"
}

run_mpi() {
  # Usage: run_mpi <label> <mdin> <rst_in> <ref> <outdir> <np>
  local label="$1" mdin="$2" rst_in="$3" ref="$4" odir="$5" np="$6"
  local rst_out="${odir}/${label}.rst7"
  local out="${odir}/${label}.out"

  if [[ -f "${rst_out}" ]]; then
    log "[SKIP] ${label} — ${rst_out} exists"
    return 0
  fi
  log "--- ${label} (pmemd.MPI, ${np} ranks) ---"
  mpirun --oversubscribe -np "${np}" "${PMEMD_MPI}" -O \
    -i  "${mdin}" \
    -o  "${out}" \
    -p  "${PRMTOP}" \
    -c  "${rst_in}" \
    -r  "${rst_out}" \
    -ref "${ref}" \
    -x  "${odir}/${label}.nc" \
    -inf "${odir}/${label}.mdinfo"
  log "${label} done."
  nan_check "${out}" "${label}"
}

run_gpu() {
  # Usage: run_gpu <label> <mdin> <rst_in> <ref_or_empty> <outdir>
  local label="$1" mdin="$2" rst_in="$3" ref="${4:-}" odir="$5"
  local rst_out="${odir}/${label}.rst7"
  local out="${odir}/${label}.out"

  if [[ -f "${rst_out}" ]]; then
    log "[SKIP] ${label} — ${rst_out} exists"
    return 0
  fi
  log "--- ${label} (pmemd.cuda, CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}) ---"
  local args=(
    -O
    -i  "${mdin}"
    -o  "${out}"
    -p  "${PRMTOP}"
    -c  "${rst_in}"
    -r  "${rst_out}"
    -x  "${odir}/${label}.nc"
    -inf "${odir}/${label}.mdinfo"
  )
  [[ -n "${ref}" ]] && args+=(-ref "${ref}")
  "${PMEMD}" "${args[@]}"
  log "${label} done."
  nan_check "${out}" "${label}"
}

# ── Resumption check ─────────────────────────────────────────────────────────
EQUIL_RST="${OUTDIR}/equil/equil.rst7"
if [[ -f "${EQUIL_RST}" ]]; then
  log "[RESUME] equil/equil.rst7 found — skipping all equilibration."
else

  # ── Stage 1: Minimization 1 — solvent + ions relaxed, heavy restraints on solute ──
  run_min "min1" "${INP}/min1.in" "${INPCRD}" "${INPCRD}" "${OUTDIR}/min1"

  # ── Stage 2: Minimization 2 — active site restrained ─────────────────────────
  run_min "min2" "${INP}/min2.in" "${OUTDIR}/min1/min1.rst7" "${OUTDIR}/min1/min1.rst7" "${OUTDIR}/min2"

  # ── Stage 3: Short NVT CPU warmup 0→50 K, 1 ps (pmemd.MPI, 54 ranks) ─────────
  run_mpi "heat" "${INP}/heat_nvt_cpu.in" \
    "${OUTDIR}/min2/min2.rst7" \
    "${OUTDIR}/min2/min2.rst7" \
    "${OUTDIR}/heat" \
    "${MPI_RANKS}"

  # ── Stage 4: 100 ns free NPT equilibration (pmemd.cuda, single run) ──────────
  log "--- equil (pmemd.cuda, 100 ns, equil_free_2fs.in) ---"
  [[ -f "${EQUIL_RST}" ]] || {
    "${PMEMD}" -O \
      -i  "${INP}/equil_free_2fs.in" \
      -o  "${OUTDIR}/equil/equil.out" \
      -p  "${PRMTOP}" \
      -c  "${OUTDIR}/heat/heat.rst7" \
      -r  "${EQUIL_RST}" \
      -x  "${OUTDIR}/equil/equil.nc" \
      -inf "${OUTDIR}/equil/equil.mdinfo"
    log "equil done."
    nan_check "${OUTDIR}/equil/equil.out" "equil"
  }

  log ""
  log "=== Equilibration performance (ns/day) ==="
  nsday "${OUTDIR}/equil/equil.mdinfo"
  log ""

fi  # end equilibration block

# ── Stages 5–7: Production 3 × 20 ns ────────────────────────────────────────
log "=== Production: 3 replicates × 20 ns, dt=0.002, ig=-1 ==="
log "All replicates start from the same equil.rst7 with random seeds."

for REP in 1 2 3; do
  REPDIR="${OUTDIR}/prod_rep${REP}"
  RST_OUT="${REPDIR}/prod_rep${REP}.rst7"
  OUT="${REPDIR}/prod_rep${REP}.out"

  if [[ -f "${RST_OUT}" ]]; then
    log "[SKIP] prod_rep${REP} — ${RST_OUT} exists"
    continue
  fi

  log "--- prod_rep${REP} (pmemd.cuda, 20 ns) ---"
  "${PMEMD}" -O \
    -i  "${INP}/production_2fs_single.in" \
    -o  "${OUT}" \
    -p  "${PRMTOP}" \
    -c  "${EQUIL_RST}" \
    -r  "${RST_OUT}" \
    -x  "${REPDIR}/prod_rep${REP}.nc" \
    -inf "${REPDIR}/prod_rep${REP}.mdinfo"

  log "prod_rep${REP} done."
  nan_check "${OUT}" "prod_rep${REP}"

  log "  prod_rep${REP} performance (ns/day):"
  nsday "${REPDIR}/prod_rep${REP}.mdinfo"
  log ""
done

log ""
log "=== run_mol14_2fs.sh: ALL STAGES COMPLETE ==="
log "System    : ${SYS} (no HMR, dt=0.002)"
log "Topology  : ${PRMTOP}"
log "Equil RST : ${EQUIL_RST}"
log "Production: ${OUTDIR}/prod_rep{1,2,3}/"
log "Full log  : ${LOGFILE}"
log ""
