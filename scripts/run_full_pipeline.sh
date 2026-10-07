#!/usr/bin/env bash
# =============================================================================
# run_full_pipeline.sh — Ni-center ORCA → MCPB → monomer builds → test MD → queue
# =============================================================================
# PROJECT_ROOT must be: /path/to/urease_md
#
# Phases:
#   1 — Wait for Ni-center ORCA (small_opt / small_freq / large_mk) completion
#   2 — Hessian → fchk, MCPB.py -s 2 & -s 3, frcmod QC vs B3LYP baseline, patch charges
#   3 — Build 36 monomer systems from ligand_name_mapping.csv (tleap + HMR)
#   4 — GPU test equilibration on mol6_pose2 monomer (pmemd.cuda)
#   5 — If Phase 4 passes: sequential production for all 36 systems (optional long run)
#
# Usage (review before launch):
#   nohup bash scripts/run_full_pipeline.sh > logs/pipeline_stdout.log 2>&1 &
#
# Environment (optional overrides):
#   CUDA_VISIBLE_DEVICES   default 1 (Tesla P100-PCIE-16GB). Use 0 only for the display K2200 (4 GB).
#   PIPELINE_NI_R1         cpptraj mask for first Ni  (default :808@NI — monomer prmtop numbering)
#   PIPELINE_NI_R2         cpptraj mask for second Ni (default :809@NI)
#   MONOMER_APO_PDB        apo monomer (chains B+Q) for combine { protein LIG }
#   PIPELINE_LIGAND_REF_MOL2  reference docked mol2 for centroid QC (default mol6_gaff2.mol2)
#   PIPELINE_LIGAND_CENTROID_MAX_DIST  max |Δcentroid| vs ref (Å, default 15)
#   PIPELINE_TEST_SYSTEM   default mol6_pose2_monomer
#   RUN_PHASE5             default 1; set 0 to stop after Phase 4 flag (still runs 1–4)
#   PIPELINE_ONLY_PHASE4   default 0; set 1 to skip Phases 1–3 and run GPU test only (expects solvated_hmr.*)
#   MAPPING_CSV            default params/inhibitors/orca_opt/ligand_name_mapping.csv
#
# =============================================================================
set -euo pipefail

PROJECT_ROOT="/path/to/urease_md"
cd "${PROJECT_ROOT}"

# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config.env"

STAMP="$(date '+%Y-%m-%d_%H%M%S')"
LOG_DIR="${PROJECT_ROOT}/logs"
MAIN_LOG="${LOG_DIR}/full_pipeline_${STAMP}.log"
mkdir -p "${LOG_DIR}"
exec > >(tee -a "${MAIN_LOG}") 2>&1

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-1}"
export PATH="${AMBER_HOME:-/usr/local}/bin:${PATH}"

PMEMD_CUDA_BIN="${PMEMD_CUDA_BIN:-${AMBER_HOME:-/usr/local}/bin/${PMEMD_CUDA:-pmemd.cuda}}"
# Optional override for all pmemd dynamics (defaults to GPU). Export PMEMD_MD_BIN to force CPU pmemd for debugging.
PMEMD_MD_BIN="${PMEMD_MD_BIN:-${PMEMD_CUDA_BIN}}"
CPPTRAJ_BIN="${CPPTRAJ:-/usr/local/bin/cpptraj}"
MCPB_PY="$(command -v MCPB.py || true)"
MCPB_RUNNER="${PROJECT_ROOT}/scripts/run_mcpb_patched.py"

RUN_PHASE5="${RUN_PHASE5:-1}"
PIPELINE_TEST_SYSTEM="${PIPELINE_TEST_SYSTEM:-mol6_pose2_monomer}"
# cpptraj masks: residue numbers for the two Ni ions (dinuclear site). Defaults match monomer
# mol6_pose2_monomer / config.env EQUIL_NI1–2; override for dodecamer builds (e.g. :570@NI :571@NI).
PIPELINE_NI_R1="${PIPELINE_NI_R1:-:808@NI}"
PIPELINE_NI_R2="${PIPELINE_NI_R2:-:809@NI}"
MAPPING_CSV="${MAPPING_CSV:-${PROJECT_ROOT}/params/inhibitors/orca_opt/ligand_name_mapping.csv}"
MONOMER_APO_PDB="${MONOMER_APO_PDB:-${PROJECT_ROOT}/structures/monomer/monomer_BQ_apo.pdb}"
PIPELINE_LIGAND_REF_MOL2="${PIPELINE_LIGAND_REF_MOL2:-${PROJECT_ROOT}/params/inhibitors/mol2/mol6_gaff2.mol2}"
PIPELINE_LIGAND_CENTROID_MAX_DIST="${PIPELINE_LIGAND_CENTROID_MAX_DIST:-15}"

SMALL_OPT_DIR="${PROJECT_ROOT}/params/ni_center_no_ohb/orca_calcs/small_opt"
SMALL_FREQ_DIR="${PROJECT_ROOT}/params/ni_center_no_ohb/orca_calcs/small_freq"
LARGE_MK_DIR="${PROJECT_ROOT}/params/ni_center_no_ohb/orca_calcs/large_mk"
MCPB_DIR="${PROJECT_ROOT}/params/ni_center_no_ohb/mcpb_input"
MON_PATCH="${PROJECT_ROOT}/params/ni_center_no_ohb/monomer/set_patch_monomer_no_ohb.txt"
MON_BONDS="${PROJECT_ROOT}/params/ni_center_no_ohb/monomer/bonds_monomer_no_ohb.txt"
PROVISIONAL_AUTOPATCH="${PROVISIONAL_AUTOPATCH:-${PROJECT_ROOT}/systems/mol6_pose2_full/provisional_autopatch_optB.frcmod}"

SIM_ROOT="${PROJECT_ROOT}/simulations/full_pipeline_${STAMP}"
PHASE4_DIR="${SIM_ROOT}/phase4_${PIPELINE_TEST_SYSTEM}"
TEST_FLAG="${LOG_DIR}/test_passed.flag"

POLL_SEC=$((5 * 60))
PHASE1_TIMEOUT_SEC=$((48 * 3600))

log() { echo "[$(date '+%F %T')] $*"; }
die() { log "FAIL: $*"; exit 1; }
pass() { log "PASS: $*"; }
fail() { log "FAIL: $*"; exit 1; }

require_file() { [[ -f "$1" ]] || die "missing file: $1"; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "command not found: $1"; }

# Mol2 centroid (Å) from TRIPOS ATOM block — for docked-pose consistency checks
mol2_centroid_xyz() {
  python3 - "$1" <<'PY'
import sys
from pathlib import Path
p = Path(sys.argv[1])
lines = p.read_text().splitlines()
in_atom = False
sx = sy = sz = 0.0
n = 0
for line in lines:
    if "@<TRIPOS>ATOM" in line:
        in_atom = True
        continue
    if in_atom:
        if line.startswith("@") or not line.strip():
            break
        parts = line.split()
        if len(parts) >= 4:
            sx += float(parts[2])
            sy += float(parts[3])
            sz += float(parts[4])
            n += 1
if n == 0:
    print("nan nan nan")
    sys.exit(1)
print(f"{sx/n:.6f} {sy/n:.6f} {sz/n:.6f}")
PY
}

# Return true (0) if |centroid - ref_centroid| <= max_dist (Å)
ligand_centroid_ok_vs_ref() {
  local mol2="$1"
  local ref_xyz="$2"
  local max_d="$3"
  local xyz
  xyz=$(mol2_centroid_xyz "${mol2}")
  [[ "${xyz}" != *nan* ]] || return 1
  python3 <<PY
import math
xyz = list(map(float, "${xyz}".split()))
ref = list(map(float, "${ref_xyz}".split()))
mx = float("${max_d}")
d = math.sqrt(sum((a - b) ** 2 for a, b in zip(xyz, ref)))
raise SystemExit(0 if d <= mx else 1)
PY
}

# ----------------------------------------------------------------------------- phase 1
phase1_wait_orca() {
  log "======== PHASE 1: Ni-center ORCA completion (poll every ${POLL_SEC}s, max ${PHASE1_TIMEOUT_SEC}s) ========"
  local small_opt_out="${SMALL_OPT_DIR}/small_opt_wB97XD4.out"
  local large_mk_out="${LARGE_MK_DIR}/large_mk_wB97XD4.out"
  [[ -f "${large_mk_out}" ]] || large_mk_out="${LARGE_MK_DIR}/large_mk.out"

  local t0
  t0=$(date +%s)
  while true; do
    local now elapsed
    now=$(date +%s)
    elapsed=$((now - t0))
    if (( elapsed > PHASE1_TIMEOUT_SEC )); then
      die "Phase 1 timeout (48 h) waiting for Ni-center ORCA."
    fi

    local ok_opt ok_hess ok_chelpg
    ok_opt=0
    ok_hess=0
    ok_chelpg=0

    if [[ -f "${small_opt_out}" ]] && grep -q "OPTIMIZATION HAS CONVERGED\|HURRAY" "${small_opt_out}" 2>/dev/null; then
      ok_opt=1
    fi

    if [[ -f "${SMALL_FREQ_DIR}/small_freq_wB97XD4.hess" ]]; then
      ok_hess=1
    fi

    if [[ -f "${large_mk_out}" ]] && grep -qi "CHELPG" "${large_mk_out}" 2>/dev/null; then
      ok_chelpg=1
    fi

    local orca_busy=0
    if pgrep -f 'orca.*small_opt_wB97XD4\.inp' >/dev/null 2>&1; then orca_busy=1; fi
    if pgrep -f 'orca.*small_freq_wB97XD4\.inp' >/dev/null 2>&1; then orca_busy=1; fi
    if pgrep -f 'orca.*large_mk_wB97XD4\.inp' >/dev/null 2>&1; then orca_busy=1; fi

    if (( ok_opt && ok_hess && ok_chelpg && ! orca_busy )); then
      pass "Phase 1 — small_opt converged, small_freq .hess present, large_mk contains CHELPG, ORCA drivers idle."
      return 0
    fi

    log "Phase 1 pending — opt=${ok_opt} hess=${ok_hess} chelpg=${ok_chelpg} orca_busy=${orca_busy} (${elapsed}s elapsed)"
    sleep "${POLL_SEC}"
  done
}

# ----------------------------------------------------------------------------- phase 2
phase2_mcpb_rebuild() {
  log "======== PHASE 2: MCPB.py Seminario + charge fitting ========"
  require_cmd python3
  [[ -n "${MCPB_PY}" && -f "${MCPB_RUNNER}" ]] || die "MCPB.py not in PATH or missing ${MCPB_RUNNER}."
  export MCPB_PY

  local hess_file="${SMALL_FREQ_DIR}/small_freq_wB97XD4.hess"
  require_file "${hess_file}"

  local fchk_out="${MCPB_DIR}/urease_small_fc.fchk"
  log "Converting Hessian: ${hess_file} → ${fchk_out}"
  python3 "${PROJECT_ROOT}/scripts/orca_hess_to_fchk.py" "${hess_file}" "${fchk_out}"

  # Baseline = pre-existing parameters (B3LYP or prior MCPB) — snapshot before this run overwrites urease_no_ohb_mcpbpy.frcmod
  local baseline="${MCPB_DIR}/urease_mcpbpy_b3lyp_baseline.frcmod"
  if [[ ! -f "${baseline}" ]]; then
    if [[ -f "${MCPB_DIR}/urease_no_ohb_mcpbpy.frcmod" ]]; then
      log "Saving baseline frcmod (current urease_no_ohb_mcpbpy.frcmod) to ${baseline}"
      cp -a "${MCPB_DIR}/urease_no_ohb_mcpbpy.frcmod" "${baseline}"
    elif [[ -f "${MCPB_DIR}/urease_mcpbpy_pre.frcmod" ]]; then
      cp -a "${MCPB_DIR}/urease_mcpbpy_pre.frcmod" "${baseline}"
      log "Using urease_mcpbpy_pre.frcmod as baseline."
    else
      log "WARN: No baseline frcmod for bond QC — create ${baseline} manually or accept skip."
    fi
  fi

  # Gaussian-like ESP log for MCPB Step 3 from ORCA .vpot (ωB97X-D4 large_mk only)
  local vpot="${LARGE_MK_DIR}/large_mk_wB97XD4.vpot"
  require_file "${vpot}"
  log "Building Gaussian-compatible ESP log from ${vpot}"
  python3 "${MCPB_DIR}/orca_vpot_to_gaussian_log.py" "${vpot}" "${MCPB_DIR}/urease_large_mk.log"

  (
    cd "${MCPB_DIR}"
    log "Running MCPB.py -s 2 (Seminario)…"
    python3 "${MCPB_RUNNER}" -i urease.in -s 2 2>&1 | tee "${LOG_DIR}/mcpb_step2_${STAMP}.log"
  )
  (
    cd "${MCPB_DIR}"
    log "Running MCPB.py -s 3 (charges)…"
    python3 "${MCPB_RUNNER}" -i urease.in -s 3 2>&1 | tee "${LOG_DIR}/mcpb_step3_${STAMP}.log"
  )

  [[ -f "${MCPB_DIR}/urease_no_ohb_mcpbpy.frcmod" ]] || die "Phase 2: urease_no_ohb_mcpbpy.frcmod not produced."

  if [[ -f "${baseline}" ]]; then
    log "QC: comparing new frcmod BOND terms to baseline ${baseline}"
    python3 - "${baseline}" "${MCPB_DIR}/urease_no_ohb_mcpbpy.frcmod" <<'PY'
import sys, re
from pathlib import Path

def load_bonds(path):
    d = {}
    text = Path(path).read_text().splitlines()
    in_bond = False
    for line in text:
        if line.strip().startswith("BOND"):
            in_bond = True
            continue
        if line.strip().startswith("ANGL"):
            break
        if not in_bond or not line.strip() or line.strip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        # AMBER frcmod: "M1-Y4  k  r0" or "C -Y1  k  r0" (split breaks "C -Y1")
        key = parts[0]
        if "M1" not in key and "M2" not in key:
            continue
        k, r0 = float(parts[1]), float(parts[2])
        d[key] = (k, r0)
    return d

bpath, npath = sys.argv[1], sys.argv[2]
old = load_bonds(bpath)
new = load_bonds(npath)
if not old:
    print("WARN: baseline has no BOND entries — skip QC.")
    sys.exit(0)
flag = False
for key in sorted(set(old) & set(new)):
    ko, ro = old[key]
    kn, rn = new[key]
    dr = abs(rn - ro)
    dk = abs(kn - ko) / max(abs(ko), 1e-6)
    if dr > 0.1:
        print(f"FLAG bond {key}: r0 {ro:.4f} → {rn:.4f} (Δ={dr:.4f} Å)")
        flag = True
    if dk > 0.5:
        print(f"FLAG bond {key}: k {ko:.4f} → {kn:.4f} (rel Δ={dk:.2%})")
        flag = True
if flag:
    print("QC bond comparison: FLAG (see above).")
    sys.exit(2)
print("QC bond comparison: no bond length Δ>0.1 Å nor k rel Δ>50%.")
PY
    local qc=$?
    if [[ "${qc}" -eq 2 ]]; then
      die "Phase 2 QC: Ni–ligand bond/force constant drift exceeds thresholds (see log)."
    fi
  fi

  phase2_update_monomer_patch_charges
  pass "Phase 2 complete — MCPB frcmod and monomer patch updated."
}

phase2_update_monomer_patch_charges() {
  require_file "${MON_PATCH}"
  local ni1="${MCPB_DIR}/NI1.mol2"
  local ni2="${MCPB_DIR}/NI2.mol2"
  [[ -f "${ni1}" && -f "${ni2}" ]] || {
    log "WARN: NI1.mol2 / NI2.mol2 missing — skip automated patch charge update."
    return 0
  }
  log "Updating Ni charges in ${MON_PATCH} from MCPB mol2 files."
  python3 - "${ni1}" "${ni2}" "${MON_PATCH}" <<'PY'
import sys
from pathlib import Path

def ni_charge_from_mol2(path):
    lines = Path(path).read_text().splitlines()
    in_atom = False
    for line in lines:
        if "@<TRIPOS>ATOM" in line:
            in_atom = True
            continue
        if in_atom:
            if line.startswith("@") or not line.strip():
                break
            p = line.split()
            if len(p) >= 9 and (p[1] == "NI" or p[5].upper() == "NI"):
                return float(p[8])
    raise ValueError(f"No NI atom in {path}")

ni1, ni2, patch_path = sys.argv[1:4]
c1 = ni_charge_from_mol2(ni1)
c2 = ni_charge_from_mol2(ni2)
text = Path(patch_path).read_text().splitlines()
out = []
for line in text:
    if "mol.1046@NI" in line and "charge" in line:
        parts = line.split()
        if len(parts) >= 4:
            parts[-1] = f"{c1:.6f}"
            line = " ".join(parts)
    elif "mol.1047@NI" in line and "charge" in line:
        parts = line.split()
        if len(parts) >= 4:
            parts[-1] = f"{c2:.6f}"
            line = " ".join(parts)
    out.append(line)
Path(patch_path).write_text("\n".join(out) + "\n")
print(f"Updated mol.1046@NI charge → {c1:.6f}, mol.1047@NI → {c2:.6f}")
PY
}

# ----------------------------------------------------------------------------- phase 3
phase3_build_monomers() {
  log "======== PHASE 3: batch monomer builds (36 ligands) ========"
  [[ -f "${MAPPING_CSV}" ]] || die "Missing ${MAPPING_CSV}"
  require_cmd tleap
  require_cmd python3
  [[ -x "${CPPTRAJ_BIN}" ]] || die "cpptraj not executable: ${CPPTRAJ_BIN}"
  require_file "${MONOMER_APO_PDB}"

  local ref_centroid=""
  if [[ -f "${PIPELINE_LIGAND_REF_MOL2}" ]]; then
    ref_centroid=$(mol2_centroid_xyz "${PIPELINE_LIGAND_REF_MOL2}")
    log "Ligand centroid QC: reference ${PIPELINE_LIGAND_REF_MOL2} → (${ref_centroid}) Å; max Δ = ${PIPELINE_LIGAND_CENTROID_MAX_DIST} Å"
  else
    log "WARN: PIPELINE_LIGAND_REF_MOL2 missing — skipping ligand centroid vs reference check."
  fi

  local ok=0 fail=0
  local -a fail_reasons=()

  local ref_natoms=0
  local ref_pt="${PROJECT_ROOT}/systems/mol6_pose2_monomer/monomer_mol6.prmtop"
  if [[ -f "${ref_pt}" ]]; then
    ref_natoms=$(python3 - "${ref_pt}" <<'PY'
import sys
try:
    import parmed as pmd
    print(len(pmd.load_file(sys.argv[1]).atoms))
except Exception:
    print(0)
PY
)
  fi
  [[ "${ref_natoms}" =~ ^[0-9]+$ ]] || ref_natoms=0

  # CSV: mol_id, original_filename, compound, selected_structure_xyz, selected_pose, net_charge, multiplicity, workdir
  while IFS=, read -r mol_id _orig compound _xyz pose _netc _mult _workdir; do
    [[ "${mol_id}" == "mol_id" ]] && continue
    compound=$(echo "${compound}" | tr -d '\r')
    pose=$(echo "${pose}" | tr -d '\r')
    local sys_name="${compound}_pose${pose}_monomer"
    local sdir="${PROJECT_ROOT}/systems/${sys_name}"
    mkdir -p "${sdir}"

    local imol="${PROJECT_ROOT}/params/inhibitors/mol2/${compound}_gaff2.mol2"
    local pose_full="${PROJECT_ROOT}/systems/${compound}_pose${pose}_full"
    local ifrc_pose="${pose_full}/${compound}_pose${pose}.frcmod"
    local prov_pose="${pose_full}/provisional_autopatch_optB.frcmod"
    local ifrc="${PROJECT_ROOT}/params/inhibitors/frcmod/${compound}.frcmod"
    [[ -f "${ifrc_pose}" ]] && ifrc="${ifrc_pose}"
    local prov_use="${PROVISIONAL_AUTOPATCH}"
    [[ -f "${prov_pose}" ]] && prov_use="${prov_pose}"
    if [[ ! -f "${imol}" || ! -f "${ifrc}" ]]; then
      ((fail++)) || true
      fail_reasons+=("${sys_name}: missing mol2 or frcmod for ${compound}")
      continue
    fi

    if [[ -n "${ref_centroid}" ]]; then
      local cx
      cx=$(mol2_centroid_xyz "${imol}")
      log "  ${sys_name}: ligand mol2 centroid (${cx})"
      if ! ligand_centroid_ok_vs_ref "${imol}" "${ref_centroid}" "${PIPELINE_LIGAND_CENTROID_MAX_DIST}"; then
        ((fail++)) || true
        fail_reasons+=("${sys_name}: ligand centroid too far from reference (docked pose check — see ${imol})")
        continue
      fi
    fi

    local tleap_in="${sdir}/tleap_pipeline_${STAMP}.in"
    local tleap_log="${sdir}/tleap_pipeline_${STAMP}.log"

    cat > "${tleap_in}" <<EOF
# Generated by run_full_pipeline.sh — ${sys_name}
# Ligand coordinates come entirely from the docked GAFF2 mol2 (not gas-phase isolated coords).
# Assembly: apo monomer PDB + mol2 → combine { protein LIG }
source ${PROTEIN_FF}
source ${LIGAND_FF}
source ${WATER_FF}

addAtomTypes {
        { "M1"  "Ni" "sp3" }
        { "M2"  "Ni" "sp3" }
        { "Y1"  "O" "sp3" }
        { "Y2"  "N" "sp3" }
        { "Y3"  "N" "sp3" }
        { "Y4"  "O" "sp3" }
        { "Y5"  "N" "sp3" }
        { "Y6"  "N" "sp3" }
        { "Y7"  "O" "sp3" }
        { "Y8"  "O" "sp3" }
}

KCX = loadmol2 ${MCPB_DIR}/KCX_full_stopgap.mol2
loadamberparams ${MCPB_DIR}/KCX.frcmod
loadamberparams ${MCPB_DIR}/urease_no_ohb_mcpbpy.frcmod
loadamberparams ${ifrc}
loadamberparams ${prov_use}
LIG = loadmol2 ${imol}
protein = loadpdb ${MONOMER_APO_PDB}
mol = combine { protein LIG }
check mol

source ${MON_BONDS}
source ${MON_PATCH}

check mol
solvateoct mol ${WATER_BOX} ${BOX_BUFFER}
addionsrand mol Na+ 0
addionsrand mol Cl- 0

saveamberparm mol ${sdir}/monomer_${compound}.prmtop ${sdir}/monomer_${compound}.inpcrd
quit
EOF

    if ! tleap -f "${tleap_in}" > "${tleap_log}" 2>&1; then
      ((fail++)) || true
      fail_reasons+=("${sys_name}: tleap exit non-zero")
      continue
    fi
    if grep -qi "FATAL\|fatal\|Could not\|could not" "${tleap_log}"; then
      ((fail++)) || true
      fail_reasons+=("${sys_name}: tleap FATAL in log")
      continue
    fi

    local ptop="${sdir}/monomer_${compound}.prmtop"
    local pcrd="${sdir}/monomer_${compound}.inpcrd"
    local sz
    sz=$(stat -c%s "${ptop}" 2>/dev/null || echo 0)
    if [[ "${sz}" -le 1048576 ]]; then
      ((fail++)) || true
      fail_reasons+=("${sys_name}: prmtop too small (${sz} bytes)")
      continue
    fi
    sz=$(stat -c%s "${pcrd}" 2>/dev/null || echo 0)
    if [[ "${sz}" -le 1048576 ]]; then
      ((fail++)) || true
      fail_reasons+=("${sys_name}: inpcrd too small (${sz} bytes)")
      continue
    fi

    local qsum natoms
    if ! read -r qsum natoms < <(python3 - "${ptop}" <<'PY'
import sys
try:
    import parmed as pmd
except ImportError:
    print("nan nan")
    sys.exit(1)
p = pmd.load_file(sys.argv[1])
q = sum(p.parm_data['CHARGE'])
print(f"{q:.8f} {len(p.atoms)}")
PY
); then
      ((fail++)) || true
      fail_reasons+=("${sys_name}: parmed/charge parse failed")
      continue
    fi
    if ! python3 <<PY
q = float("${qsum}")
raise SystemExit(0 if abs(q - round(q)) < 1.0 else 1)
PY
    then
      ((fail++)) || true
      fail_reasons+=("${sys_name}: total charge not near integer (${qsum})")
      continue
    fi
    if [[ "${ref_natoms}" -gt 0 ]]; then
      if ! python3 <<PY
n, ref = int("${natoms}"), int("${ref_natoms}")
lo, hi = int(ref * 0.9), int(ref * 1.1) + 1
raise SystemExit(0 if lo <= n <= hi else 1)
PY
      then
        ((fail++)) || true
        fail_reasons+=("${sys_name}: atom count not within ±10% of mol6 reference (${natoms} vs ref ${ref_natoms})")
        continue
      fi
    fi

    python3 -c "
import parmed as pmd
# Must load coordinates with topology — prmtop alone yields wrong natoms for rst7 save
p = pmd.load_file('${ptop}', '${pcrd}')
pmd.tools.HMassRepartition(p).execute()
p.save('${sdir}/solvated_hmr.prmtop', overwrite=True)
p.save('${sdir}/solvated_hmr.rst7', overwrite=True)
" || {
      ((fail++)) || true
      fail_reasons+=("${sys_name}: HMR failed")
      continue
    }

    ((ok++)) || true
    pass "Built ${sys_name}"
  done < <(tail -n +2 "${MAPPING_CSV}")

  log "Phase 3 summary: ${ok}/36 succeeded, ${fail} failed."
  if [[ "${fail}" -gt 0 ]]; then
    for r in "${fail_reasons[@]}"; do log "  - ${r}"; done
  fi
  if [[ "${fail}" -gt 3 ]]; then
    die "Phase 3: more than 3 systems failed — stopping."
  fi
  if [[ "${ok}" -lt 1 ]]; then
    die "Phase 3: no systems built successfully."
  fi
}

# ----------------------------------------------------------------------------- MD helpers (phase 4/5)
# pmemd_run executable stage prmtop crd ref [traj.nc]
# Run from PHASE4_DIR with short symlink basenames — Fortran OPEN() can fail on long full paths (>~80 chars).
pmemd_run() {
  local exe="$1" stage="$2" prmtop="$3" crd="$4" ref="$5"
  local traj=""
  [[ $# -ge 6 ]] && traj="$6"
  local log_out="${PHASE4_DIR}/${stage}.log"
  local rst_out="${PHASE4_DIR}/${stage}.rst7"
  local mdin="${PHASE4_DIR}/${stage}.in"
  [[ -f "${mdin}" ]] || die "pmemd_run: missing ${mdin}"
  mkdir -p "${PHASE4_DIR}"
  local ps="${stage}.ptop"
  local cs="${stage}.crdin"
  local rs="${stage}.ref"
  ln -sf "$(realpath "${prmtop}")" "${PHASE4_DIR}/${ps}"
  ln -sf "$(realpath "${crd}")" "${PHASE4_DIR}/${cs}"
  ln -sf "$(realpath "${ref}")" "${PHASE4_DIR}/${rs}"
  # Redirect pmemd stdout to /dev/null so command substitution $(pmemd_run ...) only captures the restart path
  if [[ -n "${traj}" ]]; then
    echo "[$(date '+%F %T')] RUN ${stage}: ${exe} (trajectory $(basename "${traj}"))" >&2
    ( cd "${PHASE4_DIR}" && "${exe}" -O -i "${stage}.in" -o "${stage}.log" \
      -p "${ps}" -c "${cs}" -r "${stage}.rst7" -ref "${rs}" -x "$(basename "${traj}")" ) >/dev/null \
      || die "pmemd ${stage} exited with error (see ${log_out})"
  else
    echo "[$(date '+%F %T')] RUN ${stage}: ${exe}" >&2
    ( cd "${PHASE4_DIR}" && "${exe}" -O -i "${stage}.in" -o "${stage}.log" \
      -p "${ps}" -c "${cs}" -r "${stage}.rst7" -ref "${rs}" ) >/dev/null \
      || die "pmemd ${stage} exited with error (see ${log_out})"
  fi
  printf '%s\n' "${rst_out}"
}

energy_last() {
  local lg="$1"
  # Prefer Etot / EAMBER; else last NSTEP ENERGY line (pmemd 24 minimization)
  local e
  e=$(grep "Etot" "${lg}" 2>/dev/null | tail -1 | awk -F'Etot' '{print $2}' | tr '=,' ' ' | awk '{print $1}')
  if [[ -n "${e}" ]]; then
    printf '%s\n' "${e}"
    return
  fi
  e=$(grep "EAMBER" "${lg}" 2>/dev/null | tail -1 | awk '{print $3}')
  if [[ -n "${e}" ]]; then
    printf '%s\n' "${e}"
    return
  fi
  grep -E '^[[:space:]]+[0-9]+[[:space:]]+-[0-9]' "${lg}" 2>/dev/null | tail -1 | awk '{print $2}'
}

# Last *instantaneous* temperature from pmemd MDOUT. AMBER prints a final RMS fluctuation block that
# repeats NSTEP/TEMP(K) lines with TEMP≈0–1 K; tailing `grep TEMP` then `awk '{print $NF}'` wrongly returns PRESS (0.0).
amber_last_md_temp_k() {
  local lg="$1"
  grep -E '^ NSTEP.*TEMP\(K\)' "${lg}" 2>/dev/null | sed -n 's/.*TEMP(K) *= *\([0-9.]*\).*/\1/p' | awk '$1+0>150 {last=$1} END {print last}'
}

# Last plausible mass density (g/cm³) from mdout. pmemd can print a bogus ~0 value in a post-run block; skip those.
amber_last_density_g_cm3() {
  awk '/Density[[:space:]]*=/ {
    v = $NF + 0
    if (v >= 0.85 && v <= 1.2) d = v
  }
  END { if (d != "" && d == d) print d }' "$1" 2>/dev/null
}

grep_shake_error() {
  grep -qiE "SHAKE ERROR|SHAKE failure|RATTLE ERROR|Constraint error|ERROR: Calculation aborted|vlimit exceeded|Coordinate resetting cannot" "$1" 2>/dev/null
}

# Phase 5: QC last 10 ns of 100 ns equilibration (ntwx=2500 → 10000 frames; use frames 9001–10000)
phase5_qc_eq100ns_last10ns() {
  local wdir="$1"
  local prmtop="$2"
  local refcrd="$3"
  local sys_label="$4"
  local nc="${wdir}/eq100ns.nc"
  local cpp="${wdir}/qc_eq100ns_last10ns.cpptraj"
  [[ -f "${nc}" ]] || {
    log "FAIL ${sys_label}: missing ${nc}"
    return 1
  }

  cat > "${cpp}" <<EOF
parm ${prmtop}
trajin ${nc} 9001 10000
reference ${refcrd} [ref]
rms BB @CA,C,N,O&!:WAT ref [ref] out ${wdir}/qc_eq_bb.dat
distance dNi ${PIPELINE_NI_R1} ${PIPELINE_NI_R2} out ${wdir}/qc_eq_nini.dat
rms LIG :LIG ref [ref] out ${wdir}/qc_eq_lig.dat
run
EOF
  "${CPPTRAJ_BIN}" -i "${cpp}" > "${wdir}/qc_eq100ns.cpptraj.log" 2>&1 || {
    log "FAIL ${sys_label}: cpptraj QC crashed"
    return 1
  }

  local max_bb max_lig dmin dmax
  max_bb=$(awk 'NR>1 && $2+0==$2 {if($2>m||m=="") m=$2} END{print m+0}' "${wdir}/qc_eq_bb.dat")
  max_lig=$(awk 'NR>1 && $2+0==$2 {if($2>m||m=="") m=$2} END{print m+0}' "${wdir}/qc_eq_lig.dat")
  dmin=$(awk 'NR>1 && $2+0==$2 {v=$2; if(min==""||v<min) min=v} END{print min+0}' "${wdir}/qc_eq_nini.dat")
  dmax=$(awk 'NR>1 && $2+0==$2 {v=$2; if(max==""||v>max) max=v} END{print max+0}' "${wdir}/qc_eq_nini.dat")

  log "Phase5 QC ${sys_label} (last 10 ns): max backbone RMSD=${max_bb} Å (limit ${PIPELINE_EQ_BB_RMS_MAX}), Ni–Ni [${dmin}, ${dmax}] Å, max ligand RMSD=${max_lig} Å (limit ${PIPELINE_EQ_LIG_RMS_MAX})"

  awk -v m="${max_bb}" -v lim="${PIPELINE_EQ_BB_RMS_MAX}" 'BEGIN{exit !(m+0 <= lim+0)}' || {
    log "FAIL ${sys_label}: backbone RMSD max ${max_bb} > ${PIPELINE_EQ_BB_RMS_MAX}"
    return 1
  }
  awk -v lo="${dmin}" -v hi="${dmax}" 'BEGIN{exit !(lo+0 >= 3.0 && hi+0 <= 4.5)}' || {
    log "FAIL ${sys_label}: Ni–Ni distance out of 3.0–4.5 Å over last 10 ns (min=${dmin}, max=${dmax})"
    return 1
  }
  awk -v m="${max_lig}" -v lim="${PIPELINE_EQ_LIG_RMS_MAX}" 'BEGIN{exit !(m+0 <= lim+0)}' || {
    log "FAIL ${sys_label}: ligand RMSD max ${max_lig} > ${PIPELINE_EQ_LIG_RMS_MAX} (pocket)"
    return 1
  }
  pass "Phase5 equilibration QC passed for ${sys_label}"
  return 0
}

# ----------------------------------------------------------------------------- phase 4
phase4_test_gpu() {
  log "======== PHASE 4: test equilibration (${PIPELINE_TEST_SYSTEM}) ========"
  [[ -x "${PMEMD_CUDA_BIN}" ]] || die "pmemd.cuda not executable: ${PMEMD_CUDA_BIN}"
  [[ -x "${PMEMD_MD_BIN}" ]] || die "pmemd MD binary not executable: ${PMEMD_MD_BIN}"
  require_cmd "${CPPTRAJ_BIN}"

  local tdir="${PROJECT_ROOT}/systems/${PIPELINE_TEST_SYSTEM}"
  local compound="mol6"
  local prmtop="${tdir}/solvated_hmr.prmtop"
  local inpcrd="${tdir}/solvated_hmr.rst7"

  if [[ ! -f "${prmtop}" || ! -f "${inpcrd}" ]]; then
    log "Building ${PIPELINE_TEST_SYSTEM} (apo PDB + docked mol6 mol2 → combine)…"
    compound="mol6"
    local sdir="${tdir}"
    mkdir -p "${sdir}"
    require_file "${MONOMER_APO_PDB}"
    require_file "${PROJECT_ROOT}/params/inhibitors/mol2/mol6_gaff2.mol2"
    local tleap_in="${sdir}/tleap_phase4_${STAMP}.in"
    cat > "${tleap_in}" <<EOF
source ${PROTEIN_FF}
source ${LIGAND_FF}
source ${WATER_FF}
addAtomTypes {
        { "M1"  "Ni" "sp3" }
        { "M2"  "Ni" "sp3" }
        { "Y1"  "O" "sp3" }
        { "Y2"  "N" "sp3" }
        { "Y3"  "N" "sp3" }
        { "Y4"  "O" "sp3" }
        { "Y5"  "N" "sp3" }
        { "Y6"  "N" "sp3" }
        { "Y7"  "O" "sp3" }
        { "Y8"  "O" "sp3" }
}
KCX = loadmol2 ${MCPB_DIR}/KCX_full_stopgap.mol2
loadamberparams ${MCPB_DIR}/KCX.frcmod
loadamberparams ${MCPB_DIR}/urease_no_ohb_mcpbpy.frcmod
loadamberparams ${PROJECT_ROOT}/systems/mol6_pose2_full/mol6_pose2.frcmod
loadamberparams ${PROJECT_ROOT}/systems/mol6_pose2_full/provisional_autopatch_optB.frcmod
LIG = loadmol2 ${PROJECT_ROOT}/params/inhibitors/mol2/mol6_gaff2.mol2
protein = loadpdb ${MONOMER_APO_PDB}
mol = combine { protein LIG }
source ${MON_BONDS}
source ${MON_PATCH}
check mol
solvateoct mol ${WATER_BOX} ${BOX_BUFFER}
addionsrand mol Na+ 0
addionsrand mol Cl- 0
saveamberparm mol ${sdir}/monomer_mol6.prmtop ${sdir}/monomer_mol6.inpcrd
quit
EOF
    tleap -f "${tleap_in}" > "${sdir}/tleap_phase4_${STAMP}.log" 2>&1
    CPPTRAJ="${CPPTRAJ_BIN}" bash "${PROJECT_ROOT}/scripts/fix_atomorder_amber.sh" \
      "${sdir}/monomer_mol6.prmtop" "${sdir}/monomer_mol6.inpcrd" || die "Phase 4: fixatomorder failed"
    python3 -c "
import parmed as pmd
p = pmd.load_file('${sdir}/monomer_mol6.prmtop', '${sdir}/monomer_mol6.inpcrd')
pmd.tools.HMassRepartition(p).execute()
p.save('${prmtop}', overwrite=True)
p.save('${inpcrd}', overwrite=True)
"
  fi

  mkdir -p "${PHASE4_DIR}"

  # Backbone heavy atoms + both Ni + bridging hydroxide (OHB); OHB O–H is vulnerable during 0→300 K heat.
  local mask_heavy="(@CA,C,N,O&!:WAT,Na+,Cl-)|${PIPELINE_NI_R1}|${PIPELINE_NI_R2}|:OHB"
  local mask_ni_only="(${PIPELINE_NI_R1})|(${PIPELINE_NI_R2})"

  log "Phase 4: minimization ${PMEMD_CUDA_BIN}; MD ${PMEMD_MD_BIN} (CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES})."

  # 4a restrained min
  cat > "${PHASE4_DIR}/min1.in" <<EOF
Minimization (restrained)
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntx=1, ntb=1, cut=${PME_CUTOFF:-9.0},
  ntc=2, ntf=2,
  ntpr=500,
  ntr=1,
  restraint_wt=10.0,
  restraintmask='${mask_heavy}',
 /
EOF
  local rst1
  rst1=$(pmemd_run "${PMEMD_CUDA_BIN}" "min1" "${prmtop}" "${inpcrd}" "${inpcrd}")
  local e1
  e1=$(energy_last "${PHASE4_DIR}/min1.log")
  awk -v e="${e1}" 'BEGIN{if (e+0 < 0 && e==e) exit 0; exit 1}' || fail "Phase 4a: final energy invalid (${e1})"
  pass "Phase 4a min1 Etot=${e1}"

  # 4b unrestrained min
  cat > "${PHASE4_DIR}/min2.in" <<EOF
Minimization (unrestrained)
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntx=1, ntb=1, cut=${PME_CUTOFF:-9.0},
  ntc=2, ntf=2,
  ntpr=500,
 /
EOF
  local rst2
  rst2=$(pmemd_run "${PMEMD_CUDA_BIN}" "min2" "${prmtop}" "${rst1}" "${rst1}")
  local e2
  e2=$(energy_last "${PHASE4_DIR}/min2.log")
  # Unrestrained min need not be more negative than restrained min (restraint release can raise E).
  awk -v e="${e2}" 'BEGIN{if (e+0 < 0 && e==e) exit 0; exit 1}' || fail "Phase 4b: final energy invalid (${e2})"
  pass "Phase 4b min2 Etot=${e2} (min1 was ${e1})"

  # 4b2 O–H distance on bridging OHB (expect ~0.96 Å after min); catch bad MCPB/geometry before heating.
  local cpp_ohb="${PHASE4_DIR}/ohb_dist.cpptraj"
  cat > "${cpp_ohb}" <<EOF
parm ${prmtop}
trajin ${rst2}
distance ohb_bond :OHB@O :OHB@H out ${PHASE4_DIR}/ohb_dist.dat
run
EOF
  "${CPPTRAJ_BIN}" -i "${cpp_ohb}" > "${PHASE4_DIR}/ohb_dist.cpptraj.log" 2>&1 \
    || fail "Phase 4b2: cpptraj OHB distance failed (see ${PHASE4_DIR}/ohb_dist.cpptraj.log)"
  local ohb_d
  ohb_d=$(awk '!/^#/ && NF >= 2 && $1 ~ /^[0-9]+$/ {d = $2} END { print d }' "${PHASE4_DIR}/ohb_dist.dat")
  [[ -n "${ohb_d}" && "${ohb_d}" != "nan" ]] || fail "Phase 4b2: could not parse OHB O–H distance from ohb_dist.dat"
  awk -v d="${ohb_d}" 'BEGIN { if (d + 0 >= 0.8 && d + 0 <= 1.2) exit 0; exit 1 }' \
    || fail "Phase 4b2: OHB O–H distance ${ohb_d} Å not in 0.8–1.2 Å (inspect MCPB / minimization)"
  pass "Phase 4b2 OHB O–H bond=${ohb_d} Å (ok for heating)"

  # 4c heating — 2 fs timestep for 1 ns (cold start). 4d density uses 2 fs on GPU: pmemd.cuda can halt with
  # "Periodic box dimensions have changed too much" when jumping from NVT heat to NPT at 4 fs (HMR); 2 fs avoids it.
  # Release / longer NPT can remain at 4 fs once the box is equilibrated.
  cat > "${PHASE4_DIR}/heat.in" <<EOF
Heating NVT (2 fs, 1 ns; restraints: backbone + Ni + OHB)
 &cntrl
  imin=0, nstlim=500000, dt=0.002,
  ntx=1, irest=0, ntb=1,
  ntt=3, gamma_ln=5.0, tempi=0.0, temp0=300.0,
  ntc=2, ntf=2, cut=${PME_CUTOFF:-9.0},
  ntpr=500, ntwx=0,
  ntr=1, restraint_wt=10.0,
  restraintmask='${mask_heavy}',
 /
EOF
  local r_heat
  r_heat=$(pmemd_run "${PMEMD_MD_BIN}" "heat" "${prmtop}" "${rst2}" "${rst2}")
  local temp_heat
  temp_heat=$(amber_last_md_temp_k "${PHASE4_DIR}/heat.log")
  [[ -n "${temp_heat}" && "${temp_heat}" != "nan" ]] || fail "Phase 4c: could not read final T from heat.log"
  awk -v t="${temp_heat}" 'BEGIN{if (t+0>=295 && t+0<=305) exit 0; exit 1}' || fail "Phase 4c: final T ${temp_heat} not in 295–305 K"
  grep_shake_error "${PHASE4_DIR}/heat.log" && fail "Phase 4c: SHAKE/related error in heat.log"
  pass "Phase 4c heat final T=${temp_heat} K"

  # 4d NPT density (Berendsen) — 2 fs, 1e6 steps = 2 ns (same time as 500k×4 fs). GPU pmemd halts if dt=4 fs here
  # ("Periodic box dimensions have changed too much" after NVT→NPT); 2 fs avoids it.
  cat > "${PHASE4_DIR}/density.in" <<EOF
NPT equilibration
 &cntrl
  imin=0, nstlim=1000000, dt=0.002,
  ntx=5, irest=1, ntb=2,
  ntp=1, pres0=1.0, barostat=1,
  ntt=3, gamma_ln=5.0, temp0=300.0,
  ntc=2, ntf=2, cut=${PME_CUTOFF:-9.0},
  ntpr=5000,
  ntr=1, restraint_wt=2.0,
  restraintmask='${mask_heavy}',
 /
EOF
  local r_den
  r_den=$(pmemd_run "${PMEMD_MD_BIN}" "density" "${prmtop}" "${r_heat}" "${r_heat}")
  local dens
  dens=$(amber_last_density_g_cm3 "${PHASE4_DIR}/density.log")
  [[ -n "${dens}" ]] || fail "Phase 4d: could not read density from density.log"
  awk -v d="${dens}" 'BEGIN{if (d+0>=0.95 && d+0<=1.05) exit 0; exit 1}' || fail "Phase 4d: density ${dens} not in 0.95–1.05 g/cm³"
  pass "Phase 4d density=${dens} g/cm³"

  # 4e light restraints Ni only (NetCDF for backbone RMSD)
  local nc_rel="${PHASE4_DIR}/release.nc"
  cat > "${PHASE4_DIR}/release.in" <<EOF
NPT release (Ni only) — 2 fs (same GPU grid issue as density when stepping to 4 fs right after NPT)
 &cntrl
  imin=0, nstlim=1000000, dt=0.002,
  ntx=5, irest=1, ntb=2,
  ntp=1, pres0=1.0, barostat=1,
  ntt=3, gamma_ln=5.0, temp0=300.0,
  ntc=2, ntf=2, cut=${PME_CUTOFF:-9.0},
  ntpr=5000, ntwx=5000,
  ntr=1, restraint_wt=0.5,
  restraintmask='${mask_ni_only}',
  iwrap=1, ntxo=2,
 /
EOF
  local r_rel
  r_rel=$(pmemd_run "${PMEMD_MD_BIN}" "release" "${prmtop}" "${r_den}" "${r_den}" "${nc_rel}")

  # backbone RMSD vs start
  local cpp_rms="${PHASE4_DIR}/rms_bb.cpptraj"
  cat > "${cpp_rms}" <<EOF
parm ${prmtop}
trajin ${nc_rel}
reference ${inpcrd} [ref]
rms R1 @CA,C,N,O&!:WAT ref [ref] out ${PHASE4_DIR}/rms_bb.dat
run
EOF
  "${CPPTRAJ_BIN}" -i "${cpp_rms}" > "${PHASE4_DIR}/rms_bb.cpptraj.log" 2>&1
  local rms_bb
  rms_bb=$(tail -1 "${PHASE4_DIR}/rms_bb.dat" 2>/dev/null | awk '{print $2}')
  awk -v r="${rms_bb:-99}" 'BEGIN{if (r+0 < 3.0) exit 0; exit 1}' || fail "Phase 4e: backbone RMSD ${rms_bb} Å ≥ 3.0"
  pass "Phase 4e backbone RMSD=${rms_bb} Å"

  # 4f 5 ns unrestrained NPT test (NetCDF trajectory for cpptraj)
  local nc="${PHASE4_DIR}/test5ns.nc"
  cat > "${PHASE4_DIR}/test5ns.in" <<EOF
NPT test 5 ns — 2 fs, 2.5e6 steps = 5 ns (was 1.25e6×4 fs)
 &cntrl
  imin=0, nstlim=2500000, dt=0.002,
  ntx=5, irest=1, ntb=2,
  ntp=1, pres0=1.0, barostat=1,
  ntt=3, gamma_ln=5.0, temp0=300.0,
  ntc=2, ntf=2, cut=${PME_CUTOFF:-9.0},
  ntpr=5000, ntwx=5000,
  ig=-1, iwrap=1, ntxo=2,
 /
EOF
  local r5
  r5=$(pmemd_run "${PMEMD_MD_BIN}" "test5ns" "${prmtop}" "${r_rel}" "${r_rel}" "${nc}")
  grep_shake_error "${PHASE4_DIR}/test5ns.log" && fail "Phase 4f: SHAKE/RATTLE or fatal error in test5ns.log"
  [[ -f "${nc}" ]] || fail "Phase 4f: ${nc} not written — cannot run trajectory QC."

  # Ni–Ni distance + ligand RMSD (full test5ns trajectory)
  local cpp_d="${PHASE4_DIR}/dist_ni.cpptraj"
  cat > "${cpp_d}" <<EOF
parm ${prmtop}
trajin ${nc}
distance dNi ${PIPELINE_NI_R1} ${PIPELINE_NI_R2} out ${PHASE4_DIR}/ni_ni.dat
run
EOF
  "${CPPTRAJ_BIN}" -i "${cpp_d}" > "${PHASE4_DIR}/dist_ni.cpptraj.log" 2>&1
  local dni
  dni=$(tail -1 "${PHASE4_DIR}/ni_ni.dat" 2>/dev/null | awk '{print $2}')
  awk -v d="${dni:-0}" 'BEGIN{if (d+0>=3.0 && d+0<=4.5) exit 0; exit 1}' || fail "Phase 4f: Ni–Ni distance ${dni} Å not in 3.0–4.5"

  local cpp_l="${PHASE4_DIR}/rms_lig.cpptraj"
  cat > "${cpp_l}" <<EOF
parm ${prmtop}
trajin ${nc}
reference ${inpcrd} [ref]
rms RL :LIG ref [ref] out ${PHASE4_DIR}/rms_lig.dat
run
EOF
  "${CPPTRAJ_BIN}" -i "${cpp_l}" > "${PHASE4_DIR}/rms_lig.cpptraj.log" 2>&1
  local rms_l
  rms_l=$(tail -1 "${PHASE4_DIR}/rms_lig.dat" 2>/dev/null | awk '{print $2}')
  awk -v r="${rms_l:-99}" 'BEGIN{if (r+0 < 3.0) exit 0; exit 1}' || fail "Phase 4f: ligand RMSD ${rms_l} Å ≥ 3.0"

  local dens_f temp_f
  dens_f=$(amber_last_density_g_cm3 "${PHASE4_DIR}/test5ns.log")
  [[ -n "${dens_f}" ]] || fail "Phase 4f: could not read density from test5ns.log"
  temp_f=$(amber_last_md_temp_k "${PHASE4_DIR}/test5ns.log")
  [[ -n "${temp_f}" && "${temp_f}" != "nan" ]] || fail "Phase 4f: could not read final T from test5ns.log"
  awk -v d="${dens_f}" 'BEGIN{if (d+0>=0.99 && d+0<=1.02) exit 0; exit 1}' || fail "Phase 4f: density ${dens_f} not in 0.99–1.02"
  awk -v t="${temp_f}" 'BEGIN{if (t+0>=298 && t+0<=302) exit 0; exit 1}' || fail "Phase 4f: T ${temp_f} not in 298–302 K"

  pass "Phase 4f QC passed (Ni–Ni=${dni}, dens=${dens_f}, T=${temp_f})"

  echo "TEST EQUILIBRATION PASSED — ready to launch full 36-system queue" | tee "${LOG_DIR}/test_passed.msg"
  date -Iseconds > "${TEST_FLAG}"
  pass "Wrote ${TEST_FLAG}"
}

# ----------------------------------------------------------------------------- phase 5
phase5_production_queue() {
  log "======== PHASE 5: full production queue (36 systems, sequential) ========"
  [[ "${RUN_PHASE5}" == "1" ]] || { log "RUN_PHASE5!=1 — skip Phase 5."; return 0; }
  [[ -f "${TEST_FLAG}" ]] || die "Phase 5: ${TEST_FLAG} missing — Phase 4 did not pass."

  local prod_root="${PROJECT_ROOT}/simulations/production_queue_${STAMP}"
  mkdir -p "${prod_root}"

  while IFS=, read -r mol_id _orig compound _xyz pose _netc _mult _workdir; do
    [[ "${mol_id}" == "mol_id" ]] && continue
    compound=$(echo "${compound}" | tr -d '\r')
    pose=$(echo "${pose}" | tr -d '\r')
    local sys_name="${compound}_pose${pose}_monomer"
    local sdir="${PROJECT_ROOT}/systems/${sys_name}"
    local prmtop="${sdir}/solvated_hmr.prmtop"
    local inpcrd="${sdir}/solvated_hmr.rst7"
    [[ -f "${prmtop}" && -f "${inpcrd}" ]] || {
      log "SKIP ${sys_name}: missing solvated_hmr topology."
      continue
    }

    local wdir="${prod_root}/${sys_name}"
    mkdir -p "${wdir}"
    log "Production run for ${sys_name} — workdir ${wdir}"

    local mask_heavy="(@CA,C,N,O&!:WAT,Na+,Cl-)|${PIPELINE_NI_R1}|${PIPELINE_NI_R2}|:OHB"
    local mask_ni_only="(${PIPELINE_NI_R1})|(${PIPELINE_NI_R2})"

    cat > "${wdir}/min1.in" <<EOF
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500, ntx=1, ntb=1, cut=${PME_CUTOFF:-9.0},
  ntc=2, ntf=2, ntpr=500, ntr=1, restraint_wt=10.0,
  restraintmask='${mask_heavy}',
 /
EOF
    "${PMEMD_CUDA_BIN}" -O -i "${wdir}/min1.in" -o "${wdir}/min1.log" -p "${prmtop}" -c "${inpcrd}" -r "${wdir}/min1.rst7" -ref "${inpcrd}"

    cat > "${wdir}/min2.in" <<EOF
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500, ntx=1, ntb=1, cut=${PME_CUTOFF:-9.0}, ntc=2, ntf=2, ntpr=500,
 /
EOF
    "${PMEMD_CUDA_BIN}" -O -i "${wdir}/min2.in" -o "${wdir}/min2.log" -p "${prmtop}" -c "${wdir}/min1.rst7" -r "${wdir}/min2.rst7" -ref "${wdir}/min1.rst7"

    cat > "${wdir}/heat.in" <<EOF
 &cntrl
  imin=0, nstlim=500000, dt=0.002, ntx=1, irest=0, ntb=1, ntt=3, gamma_ln=5.0,
  tempi=0.0, temp0=300.0, ntc=2, ntf=2, cut=${PME_CUTOFF:-9.0}, ntpr=500,
  ntr=1, restraint_wt=10.0, restraintmask='${mask_heavy}',
 /
EOF
    "${PMEMD_MD_BIN}" -O -i "${wdir}/heat.in" -o "${wdir}/heat.log" -p "${prmtop}" -c "${wdir}/min2.rst7" -r "${wdir}/heat.rst7" -ref "${wdir}/min2.rst7"

    # NPT density: 2 fs / 1e6 steps = 2 ns (matches Phase 4d; GPU-safe after NVT heat)
    cat > "${wdir}/density.in" <<EOF
 &cntrl
  imin=0, nstlim=1000000, dt=0.002, ntx=5, irest=1, ntb=2,
  ntp=1, barostat=1, pres0=1.0, ntt=3, gamma_ln=5.0, temp0=300.0,
  ntc=2, ntf=2, cut=${PME_CUTOFF:-9.0}, ntpr=5000,
  ntr=1, restraint_wt=2.0, restraintmask='${mask_heavy}',
 /
EOF
    "${PMEMD_MD_BIN}" -O -i "${wdir}/density.in" -o "${wdir}/density.log" -p "${prmtop}" -c "${wdir}/heat.rst7" -r "${wdir}/density.rst7" -ref "${wdir}/heat.rst7"

    # Release NPT: 2 fs on GPU (same rationale as Phase 4e)
    cat > "${wdir}/release.in" <<EOF
 &cntrl
  imin=0, nstlim=1000000, dt=0.002, ntx=5, irest=1, ntb=2,
  ntp=1, barostat=1, pres0=1.0, ntt=3, gamma_ln=5.0, temp0=300.0,
  ntc=2, ntf=2, cut=${PME_CUTOFF:-9.0}, ntpr=5000, ntwx=5000,
  ntr=1, restraint_wt=0.5, restraintmask='${mask_ni_only}',
  iwrap=1, ntxo=2,
 /
EOF
    "${PMEMD_MD_BIN}" -O -i "${wdir}/release.in" -o "${wdir}/release.log" -p "${prmtop}" -c "${wdir}/density.rst7" -r "${wdir}/release.rst7" -ref "${wdir}/density.rst7"

    cat > "${wdir}/eq100ns.in" <<EOF
 &cntrl
  imin=0, nstlim=25000000, dt=0.004, ntx=5, irest=1, ntb=2,
  ntp=2, barostat=2, pres0=1.0, ntt=3, gamma_ln=5.0, temp0=300.0,
  ntc=2, ntf=2, cut=${PME_CUTOFF:-9.0}, ntpr=5000, ntwx=2500, ntwr=250000,
  ig=-1, iwrap=1, ntxo=2,
 /
EOF
    log "Starting 100 ns equilibration for ${sys_name}…"
    "${PMEMD_MD_BIN}" -O -i "${wdir}/eq100ns.in" -o "${wdir}/eq100ns.log" \
      -p "${prmtop}" -c "${wdir}/release.rst7" -r "${wdir}/eq100ns.rst7" \
      -x "${wdir}/eq100ns.nc" -ref "${wdir}/release.rst7"

    if ! phase5_qc_eq100ns_last10ns "${wdir}" "${prmtop}" "${inpcrd}" "${sys_name}"; then
      log "SKIP production replicates for ${sys_name} — equilibration QC failed."
      continue
    fi

    local rep
    for rep in 1 2 3; do
      cat > "${wdir}/prod20_rep${rep}.in" <<EOF
 &cntrl
  imin=0, nstlim=5000000, dt=0.004, ntx=5, irest=1, ntb=2,
  ntp=1, barostat=1, pres0=1.0, ntt=3, gamma_ln=5.0, temp0=300.0,
  ntc=2, ntf=2, cut=${PME_CUTOFF:-9.0}, ntpr=5000, ntwx=2500,
  ig=-1, iwrap=1, ntxo=2,
 /
EOF
      "${PMEMD_MD_BIN}" -O -i "${wdir}/prod20_rep${rep}.in" -o "${wdir}/prod20_rep${rep}.log" \
        -p "${prmtop}" -c "${wdir}/eq100ns.rst7" -r "${wdir}/prod20_rep${rep}.rst7" \
        -x "${wdir}/prod20_rep${rep}.nc" -ref "${wdir}/eq100ns.rst7"
    done

    log "Finished ${sys_name} — trajectories under ${wdir}"
  done < <(tail -n +2 "${MAPPING_CSV}")

  pass "Phase 5 production queue completed (see ${prod_root})."
}

# ----------------------------------------------------------------------------- main
main() {
  log "run_full_pipeline.sh start — STAMP=${STAMP}"
  log "MAIN_LOG=${MAIN_LOG}"
  log "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"

  if [[ "${PIPELINE_ONLY_PHASE4:-0}" == "1" ]]; then
    log "PIPELINE_ONLY_PHASE4=1 — skipping Phases 1–3, running GPU test (Phase 4) only."
    phase4_test_gpu
    log "run_full_pipeline.sh finished (Phase 4 only)."
    return 0
  fi

  phase1_wait_orca
  phase2_mcpb_rebuild
  phase3_build_monomers

  log "======== PHASE 3b: post-build QC (before GPU Phase 4) ========"
  export PROJECT_ROOT
  export CPPTRAJ="${CPPTRAJ:-${CPPTRAJ_BIN}}"
  export PIPELINE_NI_R1 PIPELINE_NI_R2
  if ! python3 "${PROJECT_ROOT}/scripts/phase3_post_build_qc.py"; then
    die "Phase 3 post-build QC failed — fix systems/logs before re-running."
  fi
  pass "Phase 3b post-build QC passed — starting Phase 4 GPU test."

  phase4_test_gpu
  phase5_production_queue

  log "run_full_pipeline.sh finished successfully."
}

main "$@"
