#!/bin/bash
# Wrapper for ORCA on this machine — sets MPI oversubscribe flags
# Usage: bash scripts/run_orca.sh path/to/input.inp path/to/output.out
ORCA=/opt/software/ORCA/build/orca
INP="${1:?Usage: run_orca.sh input.inp output.out}"
OUT="${2:?Usage: run_orca.sh input.inp output.out}"
env OMPI_MCA_rmaps_base_oversubscribe=1 \
    OMPI_MCA_hwloc_base_binding_policy=none \
    ${ORCA} "${INP}" > "${OUT}" 2>&1
