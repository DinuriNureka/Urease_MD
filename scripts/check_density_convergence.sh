#!/usr/bin/env bash
# Check density convergence from an equil_cpu_npt mdout (.out)
# Usage: bash scripts/check_density_convergence.sh <path_to_equil_cpu_npt.out>
set -euo pipefail

OUTFILE="${1:?Usage: $0 <equil_cpu_npt.out>}"

if [[ ! -f "${OUTFILE}" ]]; then
    echo "ERROR: file not found: ${OUTFILE}" >&2
    exit 1
fi

echo "Density convergence analysis"
echo "============================="
echo ""

awk '
function parse_density(line,   rest) {
    if (index(line, "Density") == 0) return ""
    if (index(line, "=") == 0) return ""
    sub(/^.*Density[[:space:]]*=[[:space:]]*/, "", line)
    sub(/[[:space:]].*$/, "", line)
    return line
}
/^ NSTEP =/ {
    step = ""
    tps = ""
    for (i = 1; i <= NF; i++) {
        if ($i == "NSTEP") step = $(i+2)
        if ($i == "TIME(PS)") tps = $(i+2)
    }
}
/Density/ && /=/ && $0 !~ /^\|/ {
    dens = parse_density($0)
    if (dens != "" && step != "") {
        printf "Step %7s  Time %6.1f ps  Density %8.4f", step, tps+0, dens+0
        if ((dens+0) >= 0.98) printf "  ** ABOVE 0.98 **"
        printf "\n"
    }
}
' "${OUTFILE}"

echo ""
echo "First step at or above 0.98 g/cm³:"
awk '
function parse_density(line) {
    if (index(line, "Density") == 0) return ""
    sub(/^.*Density[[:space:]]*=[[:space:]]*/, "", line)
    sub(/[[:space:]].*$/, "", line)
    return line
}
/^ NSTEP =/ {
    for (i = 1; i <= NF; i++) {
        if ($i == "NSTEP") step = $(i+2)
        if ($i == "TIME(PS)") tps = $(i+2)
    }
}
/Density/ && /=/ && $0 !~ /^\|/ {
    dens = parse_density($0)
    if (dens != "" && (dens+0) >= 0.98 && !found) {
        printf "  Step %s (%.1f ps) — density %.4f\n", step, tps+0, dens+0
        found = 1
    }
}
END {
    if (!found) print "  Not reached in this file (or threshold not met)"
}
' "${OUTFILE}"

echo ""
echo "Recommendation (after reviewing the table above):"
echo "  If density ≥ 0.98 g/cm³ by ~50 ps, consider EQUIL_CPU_NPT_NSTLIM=30000 (60 ps) in config.env + equil_cpu_npt.in"
echo "  If density ≥ 0.98 g/cm³ by ~30 ps, consider 20000 (40 ps)"
echo "  Keep 50000 (100 ps) if convergence is slow"
echo ""
echo "Note: update simulations/input_files/equil_cpu_npt.in nstlim together with EQUIL_CPU_NPT_NSTLIM in config.env."
