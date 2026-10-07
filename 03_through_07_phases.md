# Phases 3–7: Inhibitor Parameterization Through Analysis

---

## Phase 3: Inhibitor Parameterization (GAFF2 + AM1-BCC)

**Goal:** Parameterize all 36 inhibitor molecules.  
**Time estimate:** ~1 day (automated)

---

### Step 3.1 — Extract ligand structures from GFN2-xTB complexes

**DO:**

```python
#!/usr/bin/env python3
"""
scripts/extract_ligands.py
---------------------------------------------------------------------------
Extract ligand atoms from GFN2-xTB complex XYZ files.
Ligand = first N atoms. Residue block = last 257 atoms.

Usage:
    python3 scripts/extract_ligands.py /path/to/xyz_dir params/inhibitors/structures/

Inputs:
    Directory of *_pose*.xyz files from GFN2-xTB optimization.

Outputs:
    One .xyz per ligand in the output directory.
---------------------------------------------------------------------------
"""

import os
import sys
import glob
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

RESIDUE_ATOMS = 257  # Atoms in the DFT active site cluster


def extract_ligand(xyz_file: Path, output_dir: Path) -> int:
    """Extract ligand, return atom count. Returns 0 on failure."""
    lines = xyz_file.read_text().splitlines()
    total = int(lines[0].strip())
    n_lig = total - RESIDUE_ATOMS

    if n_lig <= 0:
        log.warning(f"  {xyz_file.name}: {total} atoms, expected > {RESIDUE_ATOMS}")
        return 0

    out = output_dir / f"{xyz_file.stem}_lig.xyz"
    with open(out, "w") as f:
        f.write(f"{n_lig}\n")
        f.write(f"Ligand from {xyz_file.name}\n")
        for i in range(2, 2 + n_lig):
            f.write(lines[i] + "\n")

    return n_lig


def main() -> int:
    if len(sys.argv) < 3:
        print(f"Usage: python3 {sys.argv[0]} <xyz_dir> <output_dir>")
        return 1

    xyz_dir = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    xyz_files = sorted(xyz_dir.glob("*.xyz"))
    log.info(f"Found {len(xyz_files)} XYZ files in {xyz_dir}")

    success = 0
    for xf in xyz_files:
        n = extract_ligand(xf, out_dir)
        if n > 0:
            log.info(f"  {xf.stem}: {n} atoms")
            success += 1

    log.info(f"Extracted {success}/{len(xyz_files)} ligands to {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**CHECK:**

```bash
ls params/inhibitors/structures/*_lig.xyz | wc -l
# Expected: 36 (one per unique inhibitor) or 78 (one per pose)
```

**WHY:** The GFN2-xTB complex files contain both the ligand and the 257-atom active site cluster. We need the isolated ligand structure for antechamber. The atom count arithmetic (`total − 257`) reliably separates them because the DFT input files were constructed with the ligand first.

---

### Step 3.2 — Convert to mol2 format

**DO:**

```bash
source config.env
cd "${PROJECT_ROOT}/params/inhibitors/structures"

for xyz in *_lig.xyz; do
    name="${xyz%_lig.xyz}"
    obabel "${xyz}" -O "${name}.mol2" --gen3d 2>/dev/null
    echo "  Converted: ${name}"
done
```

**CHECK:**

```bash
# Verify all mol2 files have atoms
for mol2 in *.mol2; do
    n=$(grep -c "^" "$mol2")
    echo "${mol2}: ${n} lines"
done
```

---

### Step 3.3 — Batch antechamber + parmchk2

**DO:**

```bash
#!/bin/bash
# scripts/batch_param_inhibitors.sh
# -------------------------------------------------------
# GAFF2 parameterization for all inhibitors.
#
# Usage:
#     bash scripts/batch_param_inhibitors.sh
#
# Reads charges from params/inhibitors/charges.txt (optional):
#     mol01  0
#     aha   -1
#
# Outputs:
#     params/inhibitors/mol2/*_gaff2.mol2
#     params/inhibitors/frcmod/*.frcmod
# -------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config.env"

STRUCT_DIR="${PROJECT_ROOT}/params/inhibitors/structures"
MOL2_DIR="${PROJECT_ROOT}/params/inhibitors/mol2"
FRCMOD_DIR="${PROJECT_ROOT}/params/inhibitors/frcmod"
LOG_DIR_INH="${PROJECT_ROOT}/params/inhibitors/logs"
CHARGE_FILE="${PROJECT_ROOT}/params/inhibitors/charges.txt"

mkdir -p "${MOL2_DIR}" "${FRCMOD_DIR}" "${LOG_DIR_INH}"

log() { echo "[$(date '+%H:%M:%S')] $1"; }

get_charge() {
    local name=$1
    if [ -f "${CHARGE_FILE}" ]; then
        local charge
        charge=$(grep "^${name}[[:space:]]" "${CHARGE_FILE}" 2>/dev/null | awk '{print $2}')
        [ -n "${charge}" ] && echo "${charge}" && return
    fi
    echo "0"
}

SUCCESS=0; FAILED=0; WARN=0

log "=== Batch GAFF2 Parameterization ==="

for mol2 in "${STRUCT_DIR}"/*.mol2; do
    [ -f "$mol2" ] || continue
    NAME=$(basename "$mol2" .mol2)
    NC=$(get_charge "$NAME")

    log "Processing ${NAME} (charge=${NC})..."

    # --- antechamber ---
    antechamber \
        -i "${mol2}" -fi mol2 \
        -o "${MOL2_DIR}/${NAME}_gaff2.mol2" -fo mol2 \
        -c bcc -at gaff2 -nc "${NC}" -rn LIG -s 0 \
        > "${LOG_DIR_INH}/${NAME}_ac.log" 2>&1

    if [ $? -ne 0 ]; then
        log "  [FAIL] antechamber failed for ${NAME}"
        FAILED=$((FAILED + 1))
        continue
    fi

    # --- parmchk2 ---
    parmchk2 \
        -i "${MOL2_DIR}/${NAME}_gaff2.mol2" -f mol2 \
        -o "${FRCMOD_DIR}/${NAME}.frcmod" -s gaff2 \
        > "${LOG_DIR_INH}/${NAME}_pc.log" 2>&1

    if [ $? -ne 0 ]; then
        log "  [FAIL] parmchk2 failed for ${NAME}"
        FAILED=$((FAILED + 1))
        continue
    fi

    # --- Check for ATTN warnings ---
    if grep -q "ATTN" "${FRCMOD_DIR}/${NAME}.frcmod" 2>/dev/null; then
        log "  [WARN] ${NAME}: has ATTN warnings in frcmod"
        WARN=$((WARN + 1))
    fi

    SUCCESS=$((SUCCESS + 1))

    # Clean antechamber temp files
    rm -f ANTECHAMBER_*.AC* ATOMTYPE.INF sqm.* PREP.INF NEWPDB.PDB
done

log ""
log "=== Summary ==="
log "  Success: ${SUCCESS}"
log "  Failed:  ${FAILED}"
log "  ATTN warnings: ${WARN}"

if [ ${WARN} -gt 0 ]; then
    log ""
    log "Molecules needing manual review:"
    for frc in "${FRCMOD_DIR}"/*.frcmod; do
        if grep -q "ATTN" "$frc" 2>/dev/null; then
            name=$(basename "$frc" .frcmod)
            log "  ${name}: $(grep -c 'ATTN' "$frc") parameters"
        fi
    done
fi

[ ${FAILED} -eq 0 ] && exit 0 || exit 1
```

**CHECK:**

```bash
# 1. Count successes
ls params/inhibitors/mol2/*_gaff2.mol2 | wc -l
# Expected: 36

# 2. Verify charge sums
for mol2 in params/inhibitors/mol2/*_gaff2.mol2; do
    name=$(basename "$mol2" _gaff2.mol2)
    charge=$(awk '/@<TRIPOS>ATOM/,/@<TRIPOS>/' "$mol2" | \
        grep -v "@" | awk '{sum+=$NF} END {printf "%.3f", sum}')
    echo "  ${name}: ${charge}"
done

# 3. Quick tleap test (catches missing parameters)
for mol2 in params/inhibitors/mol2/*_gaff2.mol2; do
    name=$(basename "$mol2" _gaff2.mol2)
    frc="params/inhibitors/frcmod/${name}.frcmod"
    result=$(tleap -f - <<< "
source leaprc.gaff2
LIG = loadmol2 ${mol2}
loadamberparams ${frc}
check LIG
quit" 2>&1 | grep -ci "error")
    [ "${result}" -eq 0 ] && echo "  ${name}: OK" || echo "  ${name}: ERRORS"
done
```

**WHY:** GAFF2 is the current standard for organic drug-like molecules in AMBER. AM1-BCC charges approximate RESP quality at ~1000× lower cost. The batch script processes all 36 inhibitors identically, with charge overrides from `charges.txt` for ionisable molecules. The ATTN-warning check catches parameters assigned by loose analogy that may need manual refinement.

---

## Phase 4: System Assembly

**Goal:** Build 78 solvated AMBER topologies (one per inhibitor–pose combination).  
**Time estimate:** 2–3 days (scriptable)

---

### Step 4.1 — Place each inhibitor into the dodecamer

Use the Kabsch alignment script from the original guide to superimpose each GFN2-xTB optimized pose onto one active site via the Ni atom positions.

**DO:**

```bash
# For each of the 78 pose combinations:
source config.env

python3 scripts/place_inhibitor.py \
    --protein "${PROJECT_ROOT}/structures/prepared/${PDB_ID}_prepared.pdb" \
    --complex "/path/to/gfn2xtb/mol01_pose01.xyz" \
    --output "${PROJECT_ROOT}/systems/mol01_pose01/complex.pdb" \
    --target_chain B --ni1 601 --ni2 602
```

**CHECK:**

```bash
# Verify inhibitor was placed (LIG residue present)
grep "LIG" systems/mol01_pose01/complex.pdb | wc -l
# Expected: >0

# Verify Ni atoms still present
grep "NI" systems/mol01_pose01/complex.pdb | wc -l
# Expected: 24
```

### Step 4.2 — Build topology with tleap

**DO:**

```bash
#!/bin/bash
# scripts/build_system.sh
# -------------------------------------------------------
# Build one complete solvated AMBER system.
#
# Usage:
#     bash scripts/build_system.sh mol01_pose01
#
# Inputs:
#     systems/<name>/complex.pdb
#     params/ni_center/final/*
#     params/inhibitors/mol2/<inhib>_gaff2.mol2
#     params/inhibitors/frcmod/<inhib>.frcmod
#
# Outputs:
#     systems/<name>/solvated.parm7
#     systems/<name>/solvated.rst7
#     systems/<name>/solvated_hmr.parm7  (HMR for 4 fs)
# -------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config.env"

NAME=${1:?"Usage: bash build_system.sh <system_name>"}
INHIB=$(echo "${NAME}" | sed 's/_pose[0-9]*//')

SDIR="${PROJECT_ROOT}/systems/${NAME}"
PDIR="${PROJECT_ROOT}/params/ni_center/final"
IMOL="${PROJECT_ROOT}/params/inhibitors/mol2/${INHIB}_gaff2.mol2"
IFRC="${PROJECT_ROOT}/params/inhibitors/frcmod/${INHIB}.frcmod"
BOND_CMDS="${PROJECT_ROOT}/params/ni_center/final/bond_commands.txt"

log() { echo "[$(date '+%H:%M:%S')] $1"; }

# Validate inputs
for f in "${SDIR}/complex.pdb" "${PDIR}/urease.frcmod" "${PDIR}/urease.lib" \
         "${IMOL}" "${IFRC}" "${BOND_CMDS}"; do
    [ -f "$f" ] || { log "[FAIL] Missing: $f"; exit 1; }
done

log "Building: ${NAME}"

# --- Generate tleap input ---
cat > "${SDIR}/tleap.in" << TLEAPEOF
# ============================================================
# tleap: ${NAME}
# Full dodecamer + one inhibitor in TIP3P
# Generated by build_system.sh — do not edit manually
# ============================================================

# Force fields
source ${PROTEIN_FF}
source ${LIGAND_FF}
source ${WATER_FF}

# Ni center parameters
loadamberparams ${PDIR}/urease.frcmod
loadoff ${PDIR}/urease.lib

# Inhibitor parameters
loadamberparams ${IFRC}
LIG = loadmol2 ${IMOL}

# Load complex
mol = loadpdb ${SDIR}/complex.pdb

# Metal-ligand bonds (all 12 active sites)
$(cat "${BOND_CMDS}")

# Check before solvation
check mol

# Solvate (truncated octahedron — saves ~29% water vs cubic)
solvateoct mol ${WATER_BOX} ${BOX_BUFFER}

# Neutralize
addions mol Na+ 0
addions mol Cl- 0

# Physiological salt (150 mM NaCl)
# N_ions ≈ 0.0187 × ${SALT_CONC} × N_water
# Calculated after solvation; using ~140 as estimate for ~50K waters
addionsrand mol Na+ 140 Cl- 140

# Save
saveamberparm mol ${SDIR}/solvated.parm7 ${SDIR}/solvated.rst7
savepdb mol ${SDIR}/solvated.pdb

quit
TLEAPEOF

# --- Run tleap ---
tleap -f "${SDIR}/tleap.in" > "${SDIR}/tleap.log" 2>&1

# --- Check for errors ---
if grep -qi "error\|fatal\|could not" "${SDIR}/tleap.log"; then
    log "[FAIL] tleap errors — check ${SDIR}/tleap.log"
    grep -i "error\|fatal" "${SDIR}/tleap.log"
    exit 1
fi

# --- Apply HMR ---
log "Applying HMR..."
python3 -c "
import parmed as pmd
p = pmd.load_file('${SDIR}/solvated.parm7')
pmd.tools.HMassRepartition(p).execute()
p.save('${SDIR}/solvated_hmr.parm7', overwrite=True)
"

# --- Verify ---
CHARGE=$(cpptraj -p "${SDIR}/solvated.parm7" -c "${SDIR}/solvated.rst7" <<< "charge" 2>&1 | \
    grep -i "charge" | awk '{print $NF}')
NATOMS=$(cpptraj -p "${SDIR}/solvated.parm7" <<< "parminfo" 2>&1 | \
    grep "Atoms" | awk '{print $2}')

log "[OK] ${NAME}: ${NATOMS} atoms, charge = ${CHARGE}"
log "  Topology: ${SDIR}/solvated.parm7"
log "  HMR:      ${SDIR}/solvated_hmr.parm7"
log "  Coords:   ${SDIR}/solvated.rst7"
```

**CHECK:**

```bash
# Total charge near zero
cpptraj -p systems/mol01_pose01/solvated.parm7 <<< "charge" 2>&1 | grep charge
# Expected: magnitude < 0.01

# System size reasonable
cpptraj -p systems/mol01_pose01/solvated.parm7 <<< "parminfo" 2>&1 | grep Atoms
# Expected: ~200,000 atoms

# Visual inspection
echo "pymol systems/mol01_pose01/solvated.pdb"
echo "Check: Ni coordination, ligand placement, box shape"
```

**WHY:** Each system gets its own topology because the inhibitor identity differs. The truncated octahedral box saves ~29% solvent compared to cubic, which matters enormously for a dodecamer (~50,000 fewer water molecules × 78 systems = significant storage and compute savings). HMR is applied immediately so both standard and HMR topologies are available. The bond commands file (`bond_commands.txt`) is generated once by `generate_bond_commands.py` and reused for all 78 systems.

---

## Phase 5: Equilibration

**Goal:** Multi-stage relaxation from minimization through NPT equilibration.  
**Time estimate:** ~3 hours per system on GPU

---

### Step 5.1 — Write all AMBER input files

Save these in `simulations/input_files/`. Every file has a title line, inline comments on non-obvious parameters, and matches config.env settings.

**MIN1 — Solvent-only minimization** (`min1.in`):

```
Minimization stage 1: solvent relaxation only
 500 kcal/mol/A2 on all solute atoms; water+ions free
 &cntrl
   imin=1, maxcyc=5000, ncyc=2500,   ! 2500 SD then 2500 CG
   ntb=1,                              ! PBC, constant volume
   ntr=1,                              ! Positional restraints ON
   restraintmask='!(:WAT,Na+,Cl-)',    ! Restrain everything except solvent
   restraint_wt=500.0,                 ! Very strong: hold solute rigid
   cut=9.0,                            ! PME cutoff
   ntpr=100,                           ! Energy output every 100 steps
 /
```

**MIN2 — Full system minimization** (`min2.in`):

```
Minimization stage 2: full system, no restraints
 &cntrl
   imin=1, maxcyc=10000, ncyc=5000,
   ntb=1,
   ntr=0,                              ! No restraints
   cut=9.0,
   ntpr=100,
 /
```

**HEAT — NVT 0→300 K** (`heat.in`):

```
NVT heating: 0 to 300 K over 100 ps
 Backbone + metal center restrained at 10 kcal/mol/A2
 &cntrl
   imin=0, irest=0, ntx=1,            ! New simulation, no prior velocities
   nstlim=50000, dt=0.002,            ! 100 ps at 2 fs/step
   ntc=2, ntf=2,                       ! SHAKE on H-bonds
   ntb=1, ntp=0,                       ! NVT (constant volume during heating)
   ntt=3, gamma_ln=1.0,               ! Langevin thermostat, gentle coupling
   tempi=0.0, temp0=300.0,            ! Start at 0 K, target 300 K
   ntr=1,                              ! Restraints ON
   restraintmask='(@CA,C,N,O&!:WAT,Na+,Cl-)|(:NI1,NI2,KCX,OHB&!@H=)',
   restraint_wt=10.0,
   nmropt=1,                           ! Enable temperature ramping
   ig=-1,                              ! Random seed (ESSENTIAL)
   cut=9.0,
   ntpr=500, ntwx=500, ntwr=5000,
   ioutfm=1,                           ! NetCDF trajectory
 /
 &wt type='TEMP0', istep1=0, istep2=50000,
     value1=0.0, value2=300.0, /       ! Linear ramp 0→300 K
 &wt type='END' /
```

**EQUIL1–3 — NPT with decreasing restraints** (`equil1.in`, `equil2.in`, `equil3.in`):

```
NPT equilibration stage N: restraint_wt = {10.0, 1.0, 0.1}
 Berendsen barostat for fast density convergence
 &cntrl
   imin=0, irest=1, ntx=5,            ! Continue from previous
   nstlim=100000, dt=0.002,           ! 200 ps
   ntc=2, ntf=2,
   ntb=2, ntp=1, pres0=1.0, taup=1.0, ! NPT, 1 atm
   barostat=1,                         ! Berendsen (equilibration only)
   ntt=3, gamma_ln=5.0, temp0=300.0,  ! Strong Langevin for stability
   ntr=1,
   restraintmask='(@CA,C,N,O&!:WAT,Na+,Cl-)|(:NI1,NI2,KCX,OHB&!@H=)',
   restraint_wt=REPLACE_WT,           ! 10.0 → 1.0 → 0.1 across stages
   ig=-1,
   cut=9.0,
   ntpr=500, ntwx=500, ntwr=10000,
   ioutfm=1,
 /
```

Generate all three:

```bash
for wt in 10.0 1.0 0.1; do
    stage=$( [ "$wt" = "10.0" ] && echo 1 || ([ "$wt" = "1.0" ] && echo 2 || echo 3) )
    sed "s/REPLACE_WT/${wt}/" equil_template.in > "equil${stage}.in"
done
```

**EQUIL_FREE — Near-production equilibration** (`equil_free.in`):

```
Near-production equilibration: metal center restraints only (2 ns)
 Monte Carlo barostat (correct ensemble); gentle Langevin
 &cntrl
   imin=0, irest=1, ntx=5,
   nstlim=1000000, dt=0.002,          ! 2 ns
   ntc=2, ntf=2,
   ntb=2, ntp=1, pres0=1.0,
   barostat=2,                         ! Monte Carlo (production-grade)
   ntt=3, gamma_ln=1.0, temp0=300.0,  ! Production-level coupling
   ntr=1,
   restraintmask=':NI1,NI2,KCX,OHB&!@H=',  ! Metal center only
   restraint_wt=1.0,
   ig=-1,
   cut=9.0,
   ntpr=5000, ntwx=5000, ntwr=50000,
   ioutfm=1,
 /
```

**CHECK:**

```bash
# Verify all input files parse correctly (no syntax errors)
for inp in simulations/input_files/*.in; do
    name=$(basename "$inp")
    # Quick parse test: AMBER will error on malformed &cntrl
    echo "  ${name}: $(grep -c '&cntrl' "$inp") namelist block(s)"
done

# Verify restraint masks select the right atoms
ambmask -p systems/mol01_pose01/solvated.parm7 \
    -c systems/mol01_pose01/solvated.rst7 \
    -find '(@CA,C,N,O&!:WAT,Na+,Cl-)|(:NI1,NI2,KCX,OHB&!@H=)' \
    -prnlev 1 2>&1 | tail -3
# Expected: ~24,000 atoms (backbone) + ~100 atoms (metal centers)
```

**WHY for each parameter:**

| Parameter | Choice | Scientific reason |
|-----------|--------|-------------------|
| `barostat=1` (equil) → `barostat=2` (prod) | Berendsen→MC | Berendsen converges density fast but doesn't produce correct NPT ensemble. MC barostat does, and has better GPU performance. |
| `gamma_ln=5.0` (equil) → `1.0` (prod) | Strong→gentle | 5.0 ps⁻¹ provides strong temperature control during density equilibration. 1.0 ps⁻¹ preserves natural dynamics in production. |
| `ig=-1` | Random seed | Generates different initial velocities each run. ESSENTIAL for independent replicates and proper ensemble sampling. Using a fixed seed produces identical trajectories. |
| `nmropt=1` + `&wt TEMP0` | Temperature ramp | Linear 0→300 K prevents thermal shock that could distort the protein. |
| Metal center restraints kept through EQUIL_FREE | Safety | The bonded model should maintain geometry, but 1.0 kcal/mol/Å² light restraints provide a safety net during initial equilibration. Remove for production if stable. |

---

### Step 5.2 — Run equilibration

**DO:**

```bash
#!/bin/bash
# scripts/run_equilibration.sh
# -------------------------------------------------------
# Sequential equilibration for one system.
#
# Usage:
#     bash scripts/run_equilibration.sh mol01_pose01
#
# Runs: MIN1 → MIN2 → HEAT → EQUIL1 → EQUIL2 → EQUIL3 → EQUIL_FREE
# -------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config.env"

NAME=${1:?"Usage: bash run_equilibration.sh <system_name>"}
SDIR="${PROJECT_ROOT}/systems/${NAME}"
SIM="${PROJECT_ROOT}/simulations/${NAME}"
INP="${PROJECT_ROOT}/simulations/input_files"
PARM="${SDIR}/solvated.parm7"       # 2 fs for equilibration
RST="${SDIR}/solvated.rst7"

log() { echo "[$(date '+%H:%M:%S')] [$NAME] $1"; }

run_stage() {
    local stage=$1 input=$2 prev_rst=$3
    local dir="${SIM}/${stage}"
    mkdir -p "${dir}"

    log "--- ${stage} ---"
    ${PMEMD_CUDA} -O \
        -i "${input}" \
        -o "${dir}/${stage}.out" \
        -p "${PARM}" \
        -c "${prev_rst}" \
        -r "${dir}/${stage}.rst7" \
        -x "${dir}/${stage}.nc" \
        -ref "${prev_rst}" \
        -inf "${dir}/${stage}.mdinfo"

    if [ $? -ne 0 ]; then
        log "[FAIL] ${stage}"
        exit 1
    fi
    log "[OK] ${stage}"
}

log "=== Starting equilibration ==="

run_stage min1       "${INP}/min1.in"       "${RST}"
run_stage min2       "${INP}/min2.in"       "${SIM}/min1/min1.rst7"
run_stage heat       "${INP}/heat.in"       "${SIM}/min2/min2.rst7"
run_stage equil1     "${INP}/equil1.in"     "${SIM}/heat/heat.rst7"
run_stage equil2     "${INP}/equil2.in"     "${SIM}/equil1/equil1.rst7"
run_stage equil3     "${INP}/equil3.in"     "${SIM}/equil2/equil2.rst7"
run_stage equil_free "${INP}/equil_free.in" "${SIM}/equil3/equil3.rst7"

log "=== Equilibration COMPLETE ==="
log "Final restart: ${SIM}/equil_free/equil_free.rst7"
```

**CHECK:**

```bash
# 1. All stages completed
for stage in min1 min2 heat equil1 equil2 equil3 equil_free; do
    [ -f "simulations/mol01_pose01/${stage}/${stage}.rst7" ] \
        && echo "  ${stage}: OK" || echo "  ${stage}: MISSING"
done

# 2. Temperature stable at 300 K (from last equil stage)
grep "TEMP(K)" simulations/mol01_pose01/equil_free/equil_free.out | tail -5
# Expected: 300 ± 5 K

# 3. Density converged
grep "Density" simulations/mol01_pose01/equil_free/equil_free.out | tail -5
# Expected: ~1.01-1.04 g/cm³

# 4. Ni-Ni distance stable
cpptraj -p systems/mol01_pose01/solvated.parm7 << EOF
trajin simulations/mol01_pose01/equil_free/equil_free.nc
distance NiNi :601@NI :602@NI out ni_check.dat
run
quit
EOF
tail -5 ni_check.dat
# Expected: ~3.5 ± 0.3 Å
```

---

## Phase 6: Production MD

**Goal:** Unrestrained NPT at 300 K, 1 atm using HMR topology (4 fs timestep).

---

### Step 6.1 — Production input file

**`simulations/input_files/production.in`:**

```
Production MD: NPT 300K, 1 atm, 4 fs timestep (HMR)
 Monte Carlo barostat, Langevin thermostat, no restraints
 &cntrl
   imin=0, irest=1, ntx=5,
   nstlim=2500000, dt=0.004,          ! 10 ns per segment
   ntc=2, ntf=2,                       ! SHAKE on H-bonds
   ntb=2, ntp=1, pres0=1.0,
   barostat=2,                         ! Monte Carlo (correct NPT)
   ntt=3, gamma_ln=1.0, temp0=300.0,
   ig=-1,                              ! Random seed
   ntr=0,                              ! NO restraints
   cut=9.0,
   ntpr=5000, ntwx=5000, ntwe=5000,   ! Output every 20 ps
   ntwr=-25000,                        ! Overwrite restart every 100 ps
   ioutfm=1, iwrap=1,                 ! NetCDF, wrap coords
 /
```

### Step 6.2 — Run production

**DO:**

```bash
#!/bin/bash
# scripts/run_production.sh
# -------------------------------------------------------
# Production MD for one system (10 ns segments).
#
# Usage:
#     bash scripts/run_production.sh mol01_pose01 [n_segments]
#     Default: 1 segment = 10 ns
# -------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config.env"

NAME=${1:?"Usage: bash run_production.sh <name> [n_segments]"}
NSEG=${2:-1}

SDIR="${PROJECT_ROOT}/systems/${NAME}"
PROD="${PROJECT_ROOT}/simulations/${NAME}/prod"
INP="${PROJECT_ROOT}/simulations/input_files/production.in"
PARM="${SDIR}/solvated_hmr.parm7"      # HMR topology for 4 fs

mkdir -p "${PROD}"

log() { echo "[$(date '+%H:%M:%S')] [$NAME] $1"; }

log "=== Production: ${NSEG} × 10 ns ==="

for SEG in $(seq -w 1 "${NSEG}"); do
    if [ "${SEG}" = "01" ]; then
        PREV="${PROJECT_ROOT}/simulations/${NAME}/equil_free/equil_free.rst7"
    else
        PREV_SEG=$(printf "%03d" $((10#${SEG} - 1)))
        PREV="${PROD}/prod_${PREV_SEG}.rst7"
    fi

    OUT="${PROD}/prod_${SEG}"

    # Skip if already done
    [ -f "${OUT}.rst7" ] && [ -f "${OUT}.nc" ] && { log "Seg ${SEG}: exists, skipping"; continue; }

    log "Seg ${SEG}: starting from $(basename ${PREV})"

    ${PMEMD_CUDA} -O \
        -i "${INP}" \
        -o "${OUT}.out" \
        -p "${PARM}" \
        -c "${PREV}" \
        -r "${OUT}.rst7" \
        -x "${OUT}.nc" \
        -e "${OUT}.mden" \
        -ref "${PREV}" \
        -inf "${OUT}.mdinfo"

    if [ $? -ne 0 ]; then
        log "[FAIL] Segment ${SEG}"
        exit 1
    fi

    NS_DAY=$(grep "ns/day" "${OUT}.out" | tail -1 | awk '{print $4}')
    log "[OK] Seg ${SEG}: ${NS_DAY} ns/day"
done

log "=== Production COMPLETE: $((NSEG * 10)) ns ==="
```

**CHECK:**

```bash
# Performance
grep "ns/day" simulations/mol01_pose01/prod/prod_01.out | tail -1
# Expected: 60-180 ns/day depending on GPU

# Trajectory integrity
cpptraj -p systems/mol01_pose01/solvated_hmr.parm7 << EOF
trajin simulations/mol01_pose01/prod/prod_*.nc
trajinfo
quit
EOF
```

---

## Phase 7: Analysis

### Step 7.1 — cpptraj structural analysis

See `analysis/cpptraj_scripts/` for the template scripts. Run:

```bash
#!/bin/bash
# scripts/analyze_all.sh
set -euo pipefail
source "$(dirname "$0")/../config.env"

for dir in "${PROJECT_ROOT}"/simulations/*/prod; do
    NAME=$(basename "$(dirname "$dir")")
    PARM="${PROJECT_ROOT}/systems/${NAME}/solvated_hmr.parm7"

    [ -f "${PARM}" ] || continue
    ls "${dir}"/prod_*.nc &>/dev/null || continue

    echo "Analyzing: ${NAME}"

    cpptraj -p "${PARM}" << EOF
trajin ${dir}/prod_*.nc
rms fit @CA,C,N first
rms RMSD @CA,C,N first out analysis/results/rmsd/${NAME}_rmsd.dat
atomicfluct RMSF @CA out analysis/results/rmsf/${NAME}_rmsf.dat byres
rms LIG :LIG&!@H= first out analysis/results/rmsd/${NAME}_lig.dat
distance NiNi :601@NI :602@NI out analysis/results/ni_dist/${NAME}_nini.dat
distance FLAP :322@CA :601@NI out analysis/results/rmsd/${NAME}_flap.dat
hbond HB :LIG avgout analysis/results/hbonds/${NAME}_hb_avg.dat
run
quit
EOF
done
```

### Step 7.2 — MM-GBSA binding free energy

```bash
#!/bin/bash
# scripts/run_mmpbsa.sh <system_name>
set -euo pipefail
source "$(dirname "$0")/../config.env"

NAME=${1:?"Usage: bash run_mmpbsa.sh <name>"}
PARM="${PROJECT_ROOT}/systems/${NAME}/solvated_hmr.parm7"
TRAJ="${PROJECT_ROOT}/simulations/${NAME}/prod"
OUT="${PROJECT_ROOT}/analysis/mmpbsa/results/${NAME}"

mkdir -p "${OUT}"

# Generate stripped topologies
ante-MMPBSA.py -p "${PARM}" -c "${OUT}/complex.parm7" -s ':WAT:Na+:Cl-'
ante-MMPBSA.py -p "${PARM}" -c "${OUT}/receptor.parm7" -s ':WAT:Na+:Cl-:LIG'
ante-MMPBSA.py -p "${PARM}" -c "${OUT}/ligand.parm7" -s ':WAT:Na+:Cl-' -l ':LIG'

# Run MMPBSA
MMPBSA.py -O \
    -i "${PROJECT_ROOT}/analysis/mmpbsa/mmpbsa.in" \
    -o "${OUT}/mmpbsa.dat" \
    -sp "${PARM}" \
    -cp "${OUT}/complex.parm7" \
    -rp "${OUT}/receptor.parm7" \
    -lp "${OUT}/ligand.parm7" \
    -y "${TRAJ}"/prod_*.nc

echo "[DONE] Results: ${OUT}/mmpbsa.dat"
grep "DELTA TOTAL" "${OUT}/mmpbsa.dat"
```

**MMPBSA input** (`analysis/mmpbsa/mmpbsa.in`):

```
&general
  startframe=1, endframe=500, interval=1,
  verbose=2,
/
&gb
  igb=5, saltcon=0.150,
/
```

**CHECK:**

```bash
# Binding energy should be negative for binders
grep "DELTA TOTAL" analysis/mmpbsa/results/mol01_pose01/mmpbsa.dat
# Expected: ΔG < 0 kcal/mol for active inhibitors
```

**WHY:** MM-GBSA provides an approximate binding free energy that correlates with experimental IC₅₀ (typically R² = 0.3–0.6 for congeneric series). The `igb=5` (OBC2) GB model with 150 mM salt is the most widely validated combination. The stripped topologies must keep the Ni center in the receptor (not stripped with water) because it's part of the binding site.
