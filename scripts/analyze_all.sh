#!/bin/bash
# scripts/analyze_all.sh
# -------------------------------------------------------
# cpptraj structural analysis for all completed systems.
#
# Usage:
#     bash scripts/analyze_all.sh
#
# Outputs:
#     analysis/results/rmsd/<name>_rmsd.dat
#     analysis/results/rmsf/<name>_rmsf.dat
#     analysis/results/rmsd/<name>_lig.dat
#     analysis/results/ni_dist/<name>_nini.dat
#     analysis/results/hbonds/<name>_hb_avg.dat
# -------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config.env"

log() { echo "[$(date '+%H:%M:%S')] $1"; }

ANALYZED=0
SKIPPED=0

for sysdir in "${PROJECT_ROOT}"/systems/*/; do
    NAME=$(basename "${sysdir}")
    # New layout: systems/<name>/prod_rep01/prod.nc; legacy: simulations/<name>/prod/prod_*.nc
    if [[ -f "${sysdir}/prod_rep01/prod.nc" ]]; then
        dir="${sysdir}/prod_rep01"
    elif [[ -d "${PROJECT_ROOT}/simulations/${NAME}/prod" ]]; then
        dir="${PROJECT_ROOT}/simulations/${NAME}/prod"
    else
        continue
    fi
    SYS="${PROJECT_ROOT}/systems/${NAME}"
    if [[ -f "${SYS}/solvated_hmr.prmtop" ]]; then
        PARM="${SYS}/solvated_hmr.prmtop"
    elif [[ -f "${SYS}/solvated_hmr.parm7" ]]; then
        PARM="${SYS}/solvated_hmr.parm7"
    elif [[ -f "${SYS}/solvated.parm7" ]]; then
        PARM="${SYS}/solvated.parm7"
    else
        PARM="${SYS}/solvated_hmr.parm7"
    fi

    # Monomer vs dodecamer residue numbering:
    # - Monomer heterodimer (mol6_pose2_monomer): Ni at 570/571 (cpptraj resinfo :NI).
    # - Full dodecamer topology uses Ni residues at 601/602.
    if [[ "${SYSTEM_TYPE:-}" == "monomer_heterodimer" ]]; then
        NI1_RES=570
        NI2_RES=571
    else
        NI1_RES=601
        NI2_RES=602
    fi

    [ -f "${PARM}" ] || { SKIPPED=$((SKIPPED + 1)); continue; }
    ls "${dir}"/prod_*.nc &>/dev/null || { SKIPPED=$((SKIPPED + 1)); continue; }

    log "Analyzing: ${NAME}"

    # Prefer trajectory already re-imaged by preprocess_trajectory.sh; else raw prod_*.nc + autoimage.
    REI="${dir}/${PROD_TRAJ_REIMAGED:-production_reimaged.nc}"
    if [[ -f "${REI}" ]]; then
      TRAJ_BLOCK="trajin ${REI}"
    elif [[ -f "${dir}/prod.nc" ]]; then
      TRAJ_BLOCK="trajin ${dir}/prod.nc
autoimage"
    else
      TRAJ_BLOCK="trajin ${dir}/prod_*.nc
autoimage"
    fi

    cpptraj -p "${PARM}" << EOF
${TRAJ_BLOCK}
rms fit @CA,C,N first
rms RMSD @CA,C,N first out ${PROJECT_ROOT}/analysis/results/rmsd/${NAME}_rmsd.dat
atomicfluct RMSF @CA out ${PROJECT_ROOT}/analysis/results/rmsf/${NAME}_rmsf.dat byres
rms LIG :LIG&!@H= first out ${PROJECT_ROOT}/analysis/results/rmsd/${NAME}_lig.dat
distance NiNi :${NI1_RES}@NI :${NI2_RES}@NI out ${PROJECT_ROOT}/analysis/results/ni_dist/${NAME}_nini.dat
distance FLAP :322@CA :${NI1_RES}@NI out ${PROJECT_ROOT}/analysis/results/rmsd/${NAME}_flap.dat
hbond HB :LIG avgout ${PROJECT_ROOT}/analysis/results/hbonds/${NAME}_hb_avg.dat
run
quit
EOF

    ANALYZED=$((ANALYZED + 1))
done

log ""
log "=== Analysis Summary ==="
log "  Analyzed: ${ANALYZED}"
log "  Skipped:  ${SKIPPED}"
