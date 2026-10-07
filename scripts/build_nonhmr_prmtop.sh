#!/usr/bin/env bash
# scripts/build_nonhmr_prmtop.sh
# -------------------------------------------------------
# Build a NON-HMR solvated topology for a given monomer system.
# Mirrors the original pipeline (tLeap + cpptraj fixatomorder) but
# skips the ParmEd HMassRepartition step.
#
# Usage:
#   bash scripts/build_nonhmr_prmtop.sh [system_name]
#   Default system_name: mol14_pose2_monomer
#
# Inputs:
#   systems/<name>/tleap_pipeline_*.in  (latest tleap input)
#   All source files referenced inside that tleap input
#
# Outputs (both locations):
#   systems/<name>/solvated.prmtop        ← canonical topology (no HMR)
#   systems/<name>/solvated.inpcrd        ← matching starting coordinates
#   simulations/<name>/solvated.prmtop    ← symlink to above
#   simulations/<name>/solvated.inpcrd    ← symlink to above
#
# SAFETY: Does NOT overwrite solvated_hmr.prmtop or solvated_hmr.rst7.
# -------------------------------------------------------
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "${ROOT}/config.env"

SYS="${1:-mol14_pose2_monomer}"
SYSDIR="${ROOT}/systems/${SYS}"
SIMDIR="${ROOT}/simulations/${SYS}"
CPPTRAJ_EXE="${CPPTRAJ:-/usr/local/bin/cpptraj}"

log() { echo "[$(date '+%H:%M:%S')] [build_nonhmr_prmtop] $1"; }
die() { log "[FAIL] $1" >&2; exit 1; }

# ── Sanity checks ──────────────────────────────────────────────────────────
[[ -d "${SYSDIR}" ]] || die "System directory not found: ${SYSDIR}"
[[ -d "${SIMDIR}" ]] || die "Simulation directory not found: ${SIMDIR}"

OUT_PRMTOP="${SYSDIR}/solvated.prmtop"
OUT_INPCRD="${SYSDIR}/solvated.inpcrd"

if [[ -f "${OUT_PRMTOP}" && -f "${OUT_INPCRD}" ]]; then
  log "[SKIP] ${OUT_PRMTOP} already exists — delete it to rebuild."
  exit 0
fi

# ── Find the canonical tLeap input for this system ─────────────────────────
# Use the most recently generated pipeline tleap input.
TLEAP_ORIG=$(ls -t "${SYSDIR}"/tleap_pipeline_*.in 2>/dev/null | head -1 || true)
[[ -n "${TLEAP_ORIG}" ]] || die "No tleap_pipeline_*.in found in ${SYSDIR}"
log "Using tLeap template: ${TLEAP_ORIG}"

# ── Determine inhibitor name from system name ──────────────────────────────
INHIB="$(echo "${SYS}" | sed -E 's/_pose[0-9]*//; s/_monomer$//')"
log "Inhibitor: ${INHIB}"

# ── Generate a new tLeap input saving to solvated.prmtop / solvated.inpcrd ─
TLEAP_NEW="${SYSDIR}/tleap_nonhmr_$(date +%Y%m%d_%H%M%S).in"
TLEAP_LOG="${SYSDIR}/tleap_nonhmr_$(date +%Y%m%d_%H%M%S).log"

# Copy the original tLeap input, replacing the final saveamberparm destination.
# Original saves to monomer_<inhib>.prmtop / monomer_<inhib>.inpcrd; we redirect
# to a temp name so the originals are untouched.
TMPBASE="${SYSDIR}/solvated_nonhmr_tmp"

sed "s|saveamberparm mol .*|saveamberparm mol ${TMPBASE}.prmtop ${TMPBASE}.inpcrd|g" \
  "${TLEAP_ORIG}" > "${TLEAP_NEW}"

# Also remove any savepdb line to skip writing a PDB (not needed here).
sed -i '/^savepdb /d' "${TLEAP_NEW}"

log "Running tLeap (this re-solvates the system with a fresh random ion placement)..."
tleap -f "${TLEAP_NEW}" > "${TLEAP_LOG}" 2>&1 || {
  log "tLeap log tail:"
  tail -20 "${TLEAP_LOG}"
  die "tLeap failed — check ${TLEAP_LOG}"
}
# Check for genuine errors: tleap reports "Errors = N" in exit line.
# "Warnings" are expected (646 is normal for MCPB/monomer systems).
if ! grep -q "Errors = 0" "${TLEAP_LOG}"; then
  grep -iE "^Error|could not|FATAL" "${TLEAP_LOG}" | head -10
  die "tLeap reported errors (Errors > 0) — check ${TLEAP_LOG}"
fi
log "tLeap OK. Atoms: $(head -2 "${TMPBASE}.inpcrd" | tail -1 | awk '{print $1}')"

[[ -f "${TMPBASE}.prmtop" ]] || die "tLeap did not produce ${TMPBASE}.prmtop"
[[ -f "${TMPBASE}.inpcrd" ]] || die "tLeap did not produce ${TMPBASE}.inpcrd"

# ── cpptraj fixatomorder ───────────────────────────────────────────────────
log "Running cpptraj fixatomorder (makes Ni/OHB molecule atoms contiguous)..."
[[ -x "${CPPTRAJ_EXE}" ]] || CPPTRAJ_EXE="$(command -v cpptraj 2>/dev/null)" \
  || die "cpptraj not found"

TMPDIR_FIX="$(mktemp -d)"
cleanup_tmp() { rm -rf "${TMPDIR_FIX}"; }
trap cleanup_tmp EXIT

OUTPREFIX="${TMPDIR_FIX}/fo"
TMPBASE_BASENAME="$(basename "${TMPBASE}.prmtop")"

"${CPPTRAJ_EXE}" <<EOF
parm ${TMPBASE}.prmtop
trajin ${TMPBASE}.inpcrd
fixatomorder outprefix ${OUTPREFIX}
trajout ${TMPDIR_FIX}/fixed.rst7 restart
run
EOF

OUTTOP="${OUTPREFIX}.${TMPBASE_BASENAME}"
[[ -f "${OUTTOP}" ]] || die "cpptraj fixatomorder did not produce ${OUTTOP}"
[[ -f "${TMPDIR_FIX}/fixed.rst7" ]] || die "cpptraj fixatomorder did not produce fixed.rst7"

# Move fixed topology and coordinates into place.
cp -a "${OUTTOP}" "${OUT_PRMTOP}"
cp -a "${TMPDIR_FIX}/fixed.rst7" "${OUT_INPCRD}"
rm -f "${TMPBASE}.prmtop" "${TMPBASE}.inpcrd"
log "fixatomorder done. Output: ${OUT_PRMTOP}"

# ── Verify Ni–Ni distance (~3.5 Å expected in active site) ─────────────────
log "Verifying Ni–Ni distance in active site..."
NI_DIST=$("${CPPTRAJ_EXE}" <<'CPPTRAJ_EOF' 2>/dev/null | grep "^Ni-Ni" | awk '{print $2}'
parm PRMTOP_PLACEHOLDER
trajin INPCRD_PLACEHOLDER
distance Ni-Ni :NI@NI1 :NI@NI2 out /dev/null
run
CPPTRAJ_EOF
) || true

# Run again with real paths via here-doc substitution.
NI_DIST=$("${CPPTRAJ_EXE}" <<EOF 2>&1 | grep "Ni-Ni" | tail -1 | awk '{print $NF}'
parm ${OUT_PRMTOP}
trajin ${OUT_INPCRD}
distance Ni-Ni :NI@NI1 :NI@NI2 out /dev/null
run
EOF
) || NI_DIST=""

if [[ -n "${NI_DIST}" ]]; then
  log "Ni–Ni distance: ${NI_DIST} Å (expected ~3.5 Å)"
  python3 - <<PYEOF
d = float("${NI_DIST}")
if not (2.5 <= d <= 5.0):
    print(f"[WARN] Ni-Ni distance {d:.2f} Å is outside expected 2.5–5.0 Å range — verify topology")
else:
    print(f"[OK]   Ni-Ni distance {d:.2f} Å is within expected range")
PYEOF
else
  log "[WARN] Could not parse Ni–Ni distance — verify manually with cpptraj"
fi

# ── Atom count consistency check ───────────────────────────────────────────
NATOM_NONHMR=$(python3 -c "
import re
with open('${OUT_PRMTOP}') as f:
    txt = f.read()
m = re.search(r'%FLAG POINTERS.*?%FORMAT\([^)]+\)\s+(.*?)%FLAG', txt, re.DOTALL)
if m:
    nums = m.group(1).split()
    print(nums[0])
" 2>/dev/null || echo "unknown")
log "solvated.prmtop atom count: ${NATOM_NONHMR}"

# ── Create symlinks in simulations/<name>/ ─────────────────────────────────
for F in solvated.prmtop solvated.inpcrd; do
  TARGET="${SYSDIR}/${F}"
  LINK="${SIMDIR}/${F}"
  if [[ -L "${LINK}" ]]; then
    rm -f "${LINK}"
  elif [[ -f "${LINK}" ]]; then
    log "[WARN] ${LINK} exists and is not a symlink — skipping symlink creation"
    continue
  fi
  ln -s "${TARGET}" "${LINK}"
  log "Symlink: ${LINK} → ${TARGET}"
done

log ""
log "=== build_nonhmr_prmtop: DONE ==="
log "Non-HMR topology : ${OUT_PRMTOP}"
log "Coordinates       : ${OUT_INPCRD}"
log "Symlinks          : ${SIMDIR}/solvated.{prmtop,inpcrd}"
log ""
log "Next step: bash scripts/run_mol14_2fs.sh"
