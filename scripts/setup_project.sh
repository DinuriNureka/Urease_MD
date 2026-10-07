#!/bin/bash
# scripts/setup_project.sh
# -------------------------------------------------------
# Creates the complete project directory tree.
# Run ONCE at the start of the project.
#
# Usage:
#     bash scripts/setup_project.sh /path/to/urease_md
# -------------------------------------------------------
set -euo pipefail

PROJECT=${1:?"Usage: bash setup_project.sh /path/to/project"}

dirs=(
    # Documentation
    "docs"

    # Structures pipeline
    "structures/raw"
    "structures/cleaned"
    "structures/prepared"

    # Parameters — Ni center (MCPB.py)
    "params/ni_center/model_cluster/residue_mol2"
    "params/ni_center/mcpb_input"
    "params/ni_center/orca_calcs/small_opt"
    "params/ni_center/orca_calcs/small_freq"
    "params/ni_center/orca_calcs/large_mk"
    "params/ni_center/resp_charges"
    "params/ni_center/final"

    # Parameters — Inhibitors (GAFF2)
    "params/inhibitors/structures"
    "params/inhibitors/mol2"
    "params/inhibitors/frcmod"
    "params/inhibitors/logs"

    # Assembled systems (78 total — created per-system)
    "systems"

    # Simulation runs
    "simulations/input_files"
    "simulations/scripts"

    # Analysis
    "analysis/cpptraj_scripts"
    "analysis/mmpbsa/strip_topologies"
    "analysis/mmpbsa/results"
    "analysis/results/rmsd"
    "analysis/results/rmsf"
    "analysis/results/ni_dist"
    "analysis/results/hbonds"
    "analysis/results/binding_energies"
    "analysis/scripts"

    # Logs
    "logs"

    # Scripts
    "scripts"
)

for d in "${dirs[@]}"; do
    mkdir -p "${PROJECT}/${d}"
done

echo "[OK] Project tree created at: ${PROJECT}"
echo "[NEXT] Copy config.env into ${PROJECT}/ and edit paths."
