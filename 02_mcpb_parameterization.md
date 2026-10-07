# Phase 2: Ni²⁺ Active Site Parameterization (MCPB.py + ORCA 6.1.1)

**Goal:** Generate bonded-model force field parameters for the dinuclear Ni²⁺ center.  
**Time estimate:** 3–7 days (ORCA frequency calculation dominates)  
**Do this ONCE — parameters apply to all 12 active sites in all 78 systems.**

**Current status (2026-03-22):**
- Steps 2.1–2.4: Done
- Step 2.5a (FM optimization): Running with CalcHess (ORCA, 28 cores)
- Step 2.5b–c (BS-DFT freq + ESP): Queued via watchdog auto-launch
- Step 2.7b (multi-frame RESP): GFN2-xTB MD running (8 cores, parallel)
- See `02b_multiframe_resp.md` for the multi-frame RESP extension

**Charge/multiplicity note:** MCPB.py auto-assigns charge=+2, mult=2 for the
71-atom model. This is incorrect for the dinuclear Ni²⁺ system. The correct
values are charge=+1, mult=5 (FM quintet) and charge=+1, mult=1 (BS singlet).
328 electrons (even) is consistent with odd multiplicities (1, 3, 5).

**Inputs:**
- `structures/prepared/6ZJA_prepared.pdb`

**Outputs:**
- `params/ni_center/final/urease.frcmod` — bond/angle force constants (Seminario)
- `params/ni_center/final/urease.lib` — AMBER library
- `params/ni_center/final/*.mol2` — residue mol2 files with multi-frame RESP charges

---

## Step 2.1 — Extract the QM model cluster

Cut a ~60–80 atom model from ONE active site containing only the first coordination shell.

**DO:**

```python
#!/usr/bin/env python3
"""
scripts/extract_qm_model.py
---------------------------------------------------------------------------
Extract a small QM model cluster from one active site of the dodecamer.

The model contains:
  - Both Ni²⁺ ions
  - KCX219 (carbamylated lysine side chain from Cγ onward)
  - Bridging OH⁻
  - His136, His138, His221, His248, His274 → truncated to methylimidazole
  - Asp362 → truncated to acetate
  - Capping H at truncation points (added later in Avogadro/GaussView)

This is deliberately smaller than the 257-atom DFT cluster. MCPB.py only
needs the immediate coordination sphere for the Seminario method.

Usage:
    python3 scripts/extract_qm_model.py \
        structures/prepared/6ZJA_prepared.pdb \
        params/ni_center/model_cluster/active_site_model.pdb \
        --chain B --ni1 601 --ni2 602

Inputs:
    Prepared PDB from Phase 1.

Outputs:
    active_site_model.pdb — ~60-80 atom model for MCPB.py
---------------------------------------------------------------------------
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# Atoms to KEEP from each coordinating residue (side-chain ring + Cβ)
HIS_KEEP: Set[str] = {
    "CB", "CG", "ND1", "CD2", "CE1", "NE2",
    "HD1", "HD2", "HE1", "HE2", "HB2", "HB3",
}

ASP_KEEP: Set[str] = {
    "CB", "CG", "OD1", "OD2", "HB2", "HB3",
}

KCX_SKIP: Set[str] = {
    "N", "H", "CA", "HA", "C", "O",  # Skip backbone; keep full side chain
}


def parse_pdb(pdb_path: Path) -> List[dict]:
    """Parse PDB into list of atom dicts."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atoms.append({
                "line": line,
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21].strip(),
                "resnum": int(line[22:26]),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "element": line[76:78].strip() if len(line) > 76 else "",
            })
    return atoms


def extract_model(
    atoms: List[dict],
    chain: str,
    ni1_res: int,
    ni2_res: int,
    his_residues: List[int],
    kcx_res: int,
    asp_res: int,
    oh_res: int,
) -> List[dict]:
    """Extract first-shell coordination sphere atoms."""
    target_residues = set([ni1_res, ni2_res, kcx_res, asp_res, oh_res] + his_residues)
    model = []

    for atom in atoms:
        if atom["chain"] != chain:
            continue
        if atom["resnum"] not in target_residues:
            continue

        resnum = atom["resnum"]
        name = atom["name"]

        # Filter by residue type
        if resnum in his_residues:
            if name not in HIS_KEEP:
                continue
        elif resnum == asp_res:
            if name not in ASP_KEEP:
                continue
        elif resnum == kcx_res:
            if name in KCX_SKIP:
                continue
        # Ni and OH: keep all atoms

        model.append(atom)

    return model


def write_pdb(atoms: List[dict], output_path: Path) -> None:
    """Write atom list as PDB with sequential numbering."""
    with open(output_path, "w") as f:
        for i, atom in enumerate(atoms, 1):
            line = atom["line"]
            # Renumber atom serial
            f.write(f"{line[:6]}{i:5d}{line[11:]}")
        f.write("END\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract QM model cluster")
    parser.add_argument("input_pdb", type=Path)
    parser.add_argument("output_pdb", type=Path)
    parser.add_argument("--chain", default="B")
    parser.add_argument("--ni1", type=int, default=601)
    parser.add_argument("--ni2", type=int, default=602)
    parser.add_argument("--his", type=int, nargs="+", default=[136, 138, 221, 248, 274])
    parser.add_argument("--kcx", type=int, default=219)
    parser.add_argument("--asp", type=int, default=362)
    parser.add_argument("--oh", type=int, default=700)
    args = parser.parse_args()

    log.info(f"Reading {args.input_pdb}")
    atoms = parse_pdb(args.input_pdb)
    log.info(f"  Total atoms in PDB: {len(atoms)}")

    model = extract_model(
        atoms, args.chain,
        args.ni1, args.ni2, args.his, args.kcx, args.asp, args.oh
    )
    log.info(f"  Model cluster: {len(model)} atoms")

    if len(model) < 30:
        log.error(f"Model too small ({len(model)} atoms). Check chain and residue numbers.")
        log.error("Did you update residue numbers after pdb4amber renumbering?")
        return 1

    if len(model) > 150:
        log.warning(f"Model larger than expected ({len(model)} atoms). Review atom selection.")

    args.output_pdb.parent.mkdir(parents=True, exist_ok=True)
    write_pdb(model, args.output_pdb)
    log.info(f"  Written to {args.output_pdb}")

    # Summary
    resnums = set(a["resnum"] for a in model)
    resnames = set(a["resname"] for a in model)
    ni_count = sum(1 for a in model if a["element"] == "NI" or a["name"] == "NI")
    log.info(f"  Residues: {sorted(resnums)}")
    log.info(f"  Residue types: {sorted(resnames)}")
    log.info(f"  Ni atoms: {ni_count}")

    log.info("")
    log.info("[NEXT STEPS]")
    log.info("  1. Open in Avogadro or PyMOL")
    log.info("  2. Add capping H at each truncation point (Cβ → CH₃)")
    log.info("  3. Verify Ni coordination geometry is intact")
    log.info("  4. Save final model → proceed to Step 2.2")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**RUN:**

```bash
source config.env

python3 scripts/extract_qm_model.py \
    "${PROJECT_ROOT}/structures/prepared/${PDB_ID}_prepared.pdb" \
    "${PROJECT_ROOT}/params/ni_center/model_cluster/active_site_model.pdb" \
    --chain B --ni1 601 --ni2 602
```

**CHECK:**

```bash
# 1. Atom count is in expected range (60-80 before capping H)
grep -c "^ATOM\|^HETATM" params/ni_center/model_cluster/active_site_model.pdb
# Expected: 50-80

# 2. Both Ni atoms present
grep "NI" params/ni_center/model_cluster/active_site_model.pdb | wc -l
# Expected: 2

# 3. Visual check (MANDATORY)
echo "pymol params/ni_center/model_cluster/active_site_model.pdb"
echo "Verify: 2 Ni, 5 His rings, 1 Asp carboxylate, 1 KCX, 1 OH"
```

**WHY:** The QM model cluster must be small enough for a DFT frequency calculation (hours, not days) but large enough to capture the full first-shell coordination. The Seminario method needs accurate Hessian elements for Ni–N and Ni–O bonds only — atoms beyond the first shell contribute noise, not signal.

---

## Step 2.2 — Prepare individual residue mol2 files

MCPB.py requires one mol2 file per non-standard residue and per metal ion.

**DO:**

```bash
cd "${PROJECT_ROOT}/params/ni_center/model_cluster/residue_mol2"
```

Create **NI1.mol2** (first nickel ion):
```
@<TRIPOS>MOLECULE
NI1
 1 0 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 NI        0.0000    0.0000    0.0000 Ni      1  NI1        2.0000
@<TRIPOS>BOND
```

Create **NI2.mol2** (second nickel — DIFFERENT name):
```
@<TRIPOS>MOLECULE
NI2
 1 0 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 NI        3.2000    0.0000    0.0000 Ni      1  NI2        2.0000
@<TRIPOS>BOND
```

**CRITICAL:** NI1 and NI2 MUST have different residue names. MCPB.py will crash with a `KeyError` if both metals share the same name in a dinuclear system. This was reported on the AMBER mailing list for Fe₂S₂ clusters in December 2023.

Generate KCX and OH mol2 files:

```bash
# Extract KCX atoms from model cluster, generate mol2
antechamber -i kcx_extracted.pdb -fi pdb -o KCX.mol2 -fo mol2 \
    -at amber -nc -1 -rn KCX -s 2 2>&1 | tee kcx_antechamber.log

# Extract bridging OH, generate mol2
antechamber -i oh_extracted.pdb -fi pdb -o OHB.mol2 -fo mol2 \
    -at amber -nc -1 -rn OHB -s 2 2>&1 | tee oh_antechamber.log
```

Also create placeholder frcmod files (MCPB.py will populate them):

```bash
for res in KCX OHB; do
    cat > "${res}.frcmod" << 'FRCEOF'
Remark placeholder — MCPB.py will generate actual parameters
MASS

BOND

ANGLE

DIHE

IMPROPER

NONBON

FRCEOF
done
```

**CHECK:**

```bash
# Verify all mol2 files have correct residue names
for f in NI1.mol2 NI2.mol2 KCX.mol2 OHB.mol2; do
    name=$(head -1 "$f" | tr -d '[:space:]')
    resname=$(grep "^@<TRIPOS>MOLECULE" -A1 "$f" | tail -1 | tr -d '[:space:]')
    echo "$f → residue name: $resname"
done
# Expected: NI1, NI2, KCX, OHB (all different)

# Verify frcmod placeholders exist
ls -la *.frcmod
```

**WHY:** MCPB.py builds the bonded model by connecting these mol2 residues. Each must have a unique identity so MCPB.py can track which metal binds to which ligand atoms. The frcmod placeholders are required by the MCPB.py input parser even though Step 2 will overwrite them.

---

## Step 2.3 — Write the MCPB.py input file

**DO:**

```bash
cd "${PROJECT_ROOT}/params/ni_center/mcpb_input"

cat > urease.in << 'MCPBEOF'
# -------------------------------------------------------
# MCPB.py input for dinuclear Ni²⁺ urease active site
# H. pylori urease (6ZJA) — one active site model
#
# NOTES:
#   - Two metal ions → list both in ion_ids
#   - force_field ff19SB → current AMBER standard
#   - software_version g16 → we convert to ORCA manually
#     (ORCA native interface is not fully validated)
# -------------------------------------------------------

original_pdb    ../model_cluster/active_site_model.pdb
group_name      urease
cut_off         2.8
ion_ids         601 602
ion_mol2files   ../model_cluster/residue_mol2/NI1.mol2 ../model_cluster/residue_mol2/NI2.mol2
naa_mol2files   ../model_cluster/residue_mol2/KCX.mol2 ../model_cluster/residue_mol2/OHB.mol2
frcmod_files    ../model_cluster/residue_mol2/KCX.frcmod ../model_cluster/residue_mol2/OHB.frcmod
large_opt       1
force_field     ff19SB
software_version g16
MCPBEOF
```

**CHECK:**

```bash
# Verify all referenced files exist
grep -E "original_pdb|mol2files|frcmod_files" urease.in | while read key val rest; do
    for f in $val $rest; do
        [ -f "$f" ] && echo "[OK] $f" || echo "[MISSING] $f"
    done
done
```

**WHY:**
- `cut_off 2.8` — Detects all Ni–ligand bonds (typical Ni–N = 2.05 Å, Ni–O = 2.05 Å) with sufficient margin for slightly elongated bonds but avoids picking up second-shell contacts.
- `large_opt 1` — Uses the larger model (more atoms around the metal center) for RESP charge fitting, giving a better electrostatic potential representation.
- `software_version g16` — Generates Gaussian-format QM input files. We manually convert these to ORCA format (Step 2.5) because the native ORCA interface in MCPB.py is incomplete for RESP charge extraction.

---

## Step 2.4 — Run MCPB.py Step 1 (generate QM templates)

**DO:**

```bash
cd "${PROJECT_ROOT}/params/ni_center/mcpb_input"

MCPB.py -i urease.in -s 1 2>&1 | tee "${LOG_DIR}/mcpb_step1.log"
```

**CHECK:**

```bash
# Three Gaussian input files should be generated
ls -la urease_small_opt.com urease_small_fc.com urease_large_mk.com

# Each should have coordinates and charge/multiplicity
head -20 urease_small_opt.com
# Expected line: "1 5" (charge=+1, multiplicity=5 for high-spin dinuclear Ni²⁺)
```

**WHY:** MCPB.py Step 1 analyses the model PDB to identify metal-ligand bonds (based on `cut_off`), determines the small and large models, and writes QM input files for: (a) geometry optimization of the small model, (b) frequency calculation on the small model (for Hessian/Seminario), (c) ESP calculation on the large model (for RESP charges). The charge of +1 and multiplicity of 5 reflect two ferromagnetically coupled Ni²⁺ (each d⁸, S=1) giving S_total=2, 2S+1=5.

---

## Step 2.5 — Convert to ORCA and run QM calculations

### 2.5a — Convert Gaussian input to ORCA XYZ

**DO:**

```python
#!/usr/bin/env python3
"""
scripts/gaussian_to_orca_xyz.py
---------------------------------------------------------------------------
Extract coordinates from MCPB.py-generated Gaussian .com file
and write as ORCA-compatible XYZ file.

Usage:
    python3 scripts/gaussian_to_orca_xyz.py urease_small_opt.com small_model.xyz
---------------------------------------------------------------------------
"""

import sys
import logging
from pathlib import Path
from typing import List, Tuple

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def parse_gaussian_com(com_file: Path) -> Tuple[int, int, List[str]]:
    """
    Parse Gaussian .com file.
    Returns: (charge, multiplicity, coordinate_lines)
    """
    with open(com_file) as f:
        lines = f.readlines()

    charge, mult = 0, 1
    coords: List[str] = []
    found_chg_mult = False

    for line in lines:
        stripped = line.strip()

        # Blank line after route → next non-blank with 2 integers is charge/mult
        if not found_chg_mult:
            parts = stripped.split()
            if len(parts) == 2:
                try:
                    charge = int(parts[0])
                    mult = int(parts[1])
                    found_chg_mult = True
                    continue
                except ValueError:
                    pass
        elif found_chg_mult:
            parts = stripped.split()
            if len(parts) >= 4:
                # Atom line: element x y z
                try:
                    float(parts[1])
                    float(parts[2])
                    float(parts[3])
                    coords.append(stripped)
                except ValueError:
                    continue
            elif stripped == "" and len(coords) > 0:
                break  # End of coordinate block

    return charge, mult, coords


def write_xyz(coords: List[str], output_file: Path, comment: str = "") -> None:
    """Write standard XYZ file."""
    with open(output_file, "w") as f:
        f.write(f"{len(coords)}\n")
        f.write(f"{comment}\n")
        for line in coords:
            f.write(line + "\n")


def main() -> int:
    if len(sys.argv) != 3:
        print(f"Usage: python3 {sys.argv[0]} <input.com> <output.xyz>")
        return 1

    com_file = Path(sys.argv[1])
    xyz_file = Path(sys.argv[2])

    charge, mult, coords = parse_gaussian_com(com_file)
    log.info(f"Parsed {com_file.name}: charge={charge}, mult={mult}, atoms={len(coords)}")

    write_xyz(coords, xyz_file, f"charge={charge} mult={mult} from {com_file.name}")
    log.info(f"Written {xyz_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**RUN:**

```bash
source config.env

# Convert coordinates for all three QM jobs
for job in small_opt small_fc large_mk; do
    python3 scripts/gaussian_to_orca_xyz.py \
        "params/ni_center/mcpb_input/urease_${job}.com" \
        "params/ni_center/orca_calcs/${job}/${job}.xyz"
done
```

### 2.5b — Write ORCA input files

**Small model optimization** (`params/ni_center/orca_calcs/small_opt/small_opt.inp`):

```
# -------------------------------------------------------
# ORCA 6.1.1 — Geometry optimization of small model
# Dinuclear Ni²⁺ urease active site
# Charge: +1  Multiplicity: 5 (high-spin quintet)
# -------------------------------------------------------

! B3LYP def2-SVP def2/J RIJCOSX TightSCF Opt D3BJ DEFGRID2
! KeepDens

%maxcore 2800

%pal
  nprocs 16
end

%scf
  MaxIter 300
end

%geom
  MaxIter 200
end

* xyzfile 1 5 small_opt.xyz
```

**Small model frequency** (`params/ni_center/orca_calcs/small_freq/small_freq.inp`):

```
# -------------------------------------------------------
# ORCA 6.1.1 — Analytical frequency calculation (Hessian)
# MUST use the OPTIMIZED geometry from small_opt
# The Hessian at a stationary point gives force constants
# for the Seminario method in MCPB.py Step 2
# -------------------------------------------------------

! B3LYP def2-SVP def2/J RIJCOSX TightSCF Freq D3BJ DEFGRID2
! KeepDens

%maxcore 2800

%pal
  nprocs 16
end

%scf
  MaxIter 300
end

* xyzfile 1 5 small_opt_optimized.xyz
```

**Large model ESP** (`params/ni_center/orca_calcs/large_mk/large_mk.inp`):

```
# -------------------------------------------------------
# ORCA 6.1.1 — Single-point ESP for RESP charge fitting
# Uses the larger model for better electrostatic potential
# CHELPG charges computed natively by ORCA
# Multiwfn used for two-stage RESP fitting (Step 2.7)
# -------------------------------------------------------

! B3LYP def2-SVP def2/J RIJCOSX TightSCF D3BJ DEFGRID2
! KeepDens CHELPG

%maxcore 2800

%pal
  nprocs 16
end

%scf
  MaxIter 300
end

%output
  Print[P_Hirshfeld] 1
end

* xyzfile 1 5 large_mk.xyz
```

### 2.5c — Run ORCA calculations sequentially

**DO:**

```bash
#!/bin/bash
# scripts/run_orca_mcpb.sh
# -------------------------------------------------------
# Run all three ORCA calculations for MCPB.py in sequence:
#   1. Small model optimization
#   2. Small model frequency (Hessian)
#   3. Large model ESP (for RESP charges)
#
# Usage:
#     bash scripts/run_orca_mcpb.sh
#
# Expected runtime: 4-12 hours total (dominated by Freq)
# -------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config.env"

ORCA="${ORCA_PATH}"
CALC_DIR="${PROJECT_ROOT}/params/ni_center/orca_calcs"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_DIR}/orca_mcpb.log"; }

# === Step 1: Optimization ===
log "=== STEP 1: Small model optimization ==="
cd "${CALC_DIR}/small_opt"

${ORCA} small_opt.inp > small_opt.out 2>&1

if grep -q "OPTIMIZATION HAS CONVERGED" small_opt.out; then
    log "[OK] Optimization converged"
    ENERGY=$(grep "FINAL SINGLE POINT ENERGY" small_opt.out | tail -1 | awk '{print $5}')
    log "  Final energy: ${ENERGY} Eh"
else
    log "[FAIL] Optimization did not converge. Check small_opt.out"
    exit 1
fi

# === Step 2: Frequency calculation (Hessian) ===
log "=== STEP 2: Frequency calculation ==="
cd "${CALC_DIR}/small_freq"

# Copy optimized geometry
cp "${CALC_DIR}/small_opt/small_opt.xyz" small_opt_optimized.xyz

${ORCA} small_freq.inp > small_freq.out 2>&1

if [ -f "small_freq.hess" ]; then
    log "[OK] Hessian file generated: small_freq.hess"
    # Check for imaginary frequencies
    N_IMAG=$(grep "***imaginary mode***" small_freq.out | wc -l)
    if [ "$N_IMAG" -eq 0 ]; then
        log "  No imaginary frequencies (clean minimum)"
    elif [ "$N_IMAG" -le 3 ]; then
        log "  WARNING: ${N_IMAG} small imaginary modes (usually OK for this purpose)"
    else
        log "  WARNING: ${N_IMAG} imaginary modes — geometry may not be at minimum"
    fi
else
    log "[FAIL] No .hess file produced. Check small_freq.out"
    exit 1
fi

# === Step 3: Large model ESP ===
log "=== STEP 3: Large model ESP ==="
cd "${CALC_DIR}/large_mk"

${ORCA} large_mk.inp > large_mk.out 2>&1

if grep -q "SCF CONVERGED" large_mk.out; then
    log "[OK] SCF converged"
    # Generate molden file for Multiwfn
    ${ORCA_2MKL} large_mk -molden 2>/dev/null
    if [ -f "large_mk.molden.input" ]; then
        log "[OK] Molden file generated for Multiwfn RESP fitting"
    else
        log "[WARNING] Molden file not generated — run orca_2mkl manually"
    fi
else
    log "[FAIL] SCF did not converge. Check large_mk.out"
    exit 1
fi

log ""
log "=== ALL ORCA CALCULATIONS COMPLETE ==="
log "Next: Run RESP fitting with Multiwfn (Step 2.7)"
log "      Convert Hessian to fchk format (Step 2.6)"
```

**CHECK:**

```bash
# After each ORCA job, verify:

# Optimization converged
grep "OPTIMIZATION HAS CONVERGED" params/ni_center/orca_calcs/small_opt/small_opt.out

# Hessian file exists and is non-empty
ls -la params/ni_center/orca_calcs/small_freq/small_freq.hess
# Expected: >10 KB for a 60-atom system

# ESP calculation converged
grep "SCF CONVERGED" params/ni_center/orca_calcs/large_mk/large_mk.out
```

**WHY:** These three calculations provide the raw QM data that MCPB.py needs:
1. **Optimization** → ensures the model is at an energy minimum so the Hessian is meaningful
2. **Frequency/Hessian** → Cartesian second derivatives for the Seminario method (extracts bond/angle force constants by projecting the Hessian onto internal coordinates)
3. **ESP** → electrostatic potential on the molecular surface for RESP charge fitting (partial charges that reproduce the QM electrostatics)

---

## Step 2.6 — Convert ORCA Hessian to Gaussian fchk format

MCPB.py Step 2 reads a Gaussian `.fchk` file. We convert the ORCA `.hess` file.

**DO:**

```python
#!/usr/bin/env python3
"""
scripts/orca_hess_to_fchk.py
---------------------------------------------------------------------------
Convert ORCA .hess file to Gaussian-style .fchk for MCPB.py Step 2.

MCPB.py's Seminario method reader expects:
  - Number of atoms
  - Atomic numbers
  - Cartesian coordinates (Bohr)
  - Cartesian force constants (lower-triangular Hessian in Hartree/Bohr²)

Usage:
    python3 scripts/orca_hess_to_fchk.py \
        params/ni_center/orca_calcs/small_freq/small_freq.hess \
        params/ni_center/mcpb_input/urease_small_fc.fchk

Inputs:
    ORCA .hess file from frequency calculation

Outputs:
    Gaussian-format .fchk with Hessian (read by MCPB.py Step 2)
---------------------------------------------------------------------------
"""

import sys
import logging
import numpy as np
from pathlib import Path
from typing import List, Tuple

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# Atomic number lookup
ELEMENT_TO_Z = {
    "H": 1, "He": 2, "C": 6, "N": 7, "O": 8, "F": 9, "S": 16,
    "Cl": 17, "Ni": 28, "Zn": 30, "Cu": 29, "Fe": 26, "Co": 27,
}


def parse_orca_hess(hess_file: Path) -> Tuple[int, List[int], np.ndarray, np.ndarray]:
    """
    Parse ORCA .hess file.
    Returns: (n_atoms, atomic_numbers, coords_bohr, hessian_matrix)
    """
    with open(hess_file) as f:
        content = f.read()

    # Split into sections
    sections = {}
    current_key = None
    current_lines: List[str] = []

    for line in content.split("\n"):
        if line.startswith("$"):
            if current_key:
                sections[current_key] = current_lines
            current_key = line.strip().lower()
            current_lines = []
        else:
            current_lines.append(line)
    if current_key:
        sections[current_key] = current_lines

    # Parse $atoms
    atom_lines = sections.get("$atoms", [])
    n_atoms = int(atom_lines[0].strip())
    atomic_numbers = []
    for i in range(1, n_atoms + 1):
        parts = atom_lines[i].split()
        element = parts[0]
        z = ELEMENT_TO_Z.get(element, int(float(parts[1])))
        atomic_numbers.append(z)

    # Parse $actual_temperature or $coordinates
    coords = np.zeros(3 * n_atoms)
    if "$coordinates" in sections:
        coord_lines = sections["$coordinates"]
        n_coords = int(coord_lines[0].strip())
        for i in range(1, n_coords + 1):
            coords[i - 1] = float(coord_lines[i].strip())

    # Parse $hessian
    hess_lines = sections.get("$hessian", [])
    n_dof = int(hess_lines[0].strip())
    hessian = np.zeros((n_dof, n_dof))

    i = 1
    col_offset = 0
    while i < len(hess_lines):
        line = hess_lines[i].strip()
        if not line:
            i += 1
            continue
        parts = line.split()
        # Column header line (all integers)
        try:
            if all(p.lstrip("-").replace(".", "").isdigit() and "." not in p for p in parts):
                col_offset = int(parts[0])
                i += 1
                continue
        except ValueError:
            pass
        # Data line
        try:
            row_idx = int(parts[0])
            for j, val in enumerate(parts[1:]):
                hessian[row_idx][col_offset + j] = float(val)
        except (ValueError, IndexError):
            pass
        i += 1

    log.info(f"Parsed: {n_atoms} atoms, {n_dof}×{n_dof} Hessian")
    return n_atoms, atomic_numbers, coords, hessian


def write_fchk(
    fchk_file: Path,
    n_atoms: int,
    atomic_numbers: List[int],
    coords_bohr: np.ndarray,
    hessian: np.ndarray
) -> None:
    """Write Gaussian-format fchk file."""
    n_dof = 3 * n_atoms
    n_hess = n_dof * (n_dof + 1) // 2

    with open(fchk_file, "w") as f:
        f.write("ORCA-to-fchk conversion for MCPB.py Seminario method\n")
        f.write("SP        RB3LYP                          def2-SVP\n")

        # Atom count
        f.write(f"Number of atoms                            I{n_atoms:>12d}\n")

        # Atomic numbers
        f.write(f"Atomic numbers                             I   N={n_atoms:>12d}\n")
        for i in range(0, n_atoms, 6):
            chunk = atomic_numbers[i:i + 6]
            f.write("".join(f"{z:12d}" for z in chunk) + "\n")

        # Coordinates (Bohr)
        f.write(f"Current cartesian coordinates              R   N={3 * n_atoms:>12d}\n")
        for i in range(0, len(coords_bohr), 5):
            chunk = coords_bohr[i:i + 5]
            f.write("".join(f"{x:16.8E}" for x in chunk) + "\n")

        # Hessian (lower triangular)
        f.write(f"Cartesian Force Constants                  R   N={n_hess:>12d}\n")
        hess_lt = []
        for i in range(n_dof):
            for j in range(i + 1):
                hess_lt.append(hessian[i][j])

        for i in range(0, len(hess_lt), 5):
            chunk = hess_lt[i:i + 5]
            f.write("".join(f"{x:16.8E}" for x in chunk) + "\n")

    log.info(f"Written {fchk_file}: {n_hess} Hessian elements")


def main() -> int:
    if len(sys.argv) != 3:
        print(f"Usage: python3 {sys.argv[0]} <input.hess> <output.fchk>")
        return 1

    hess_path = Path(sys.argv[1])
    fchk_path = Path(sys.argv[2])

    n_atoms, at_nums, coords, hessian = parse_orca_hess(hess_path)
    write_fchk(fchk_path, n_atoms, at_nums, coords, hessian)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**RUN:**

```bash
python3 scripts/orca_hess_to_fchk.py \
    params/ni_center/orca_calcs/small_freq/small_freq.hess \
    params/ni_center/mcpb_input/urease_small_fc.fchk
```

**CHECK:**

```bash
# fchk file exists and has Hessian data
grep "Cartesian Force Constants" params/ni_center/mcpb_input/urease_small_fc.fchk
# Should show "R   N=<number>" where number = 3N(3N+1)/2

# For 70 atoms: 3×70 = 210 DOF → 210×211/2 = 22155 elements
```

**WHY:** The Gaussian fchk format is the native input for MCPB.py's Seminario method implementation. The Hessian contains the same physics regardless of which QM program computed it (ORCA or Gaussian) — it's the second derivative matrix of the energy with respect to Cartesian displacements. The conversion is purely a format translation.

---

## Step 2.7 — RESP charge fitting

**UPDATED:** This step now uses multi-frame RESP fitting instead of single-frame.
See `02b_multiframe_resp.md` for the full multi-frame workflow.

**Summary of multi-frame approach:**

1. GFN2-xTB MD (ORCA native `! XTB2`) generates conformational ensemble
2. RMSD clustering selects 10 representative frames
3. ORCA CHELPG ESP at each frame (BS-DFT, same level as `large_mk.inp`)
4. Averaged CHELPG charges across all frames
5. Charges inserted into MCPB.py mol2 files before Step 2.8

**Single-frame fallback (Multiwfn):**

If multi-frame is not needed, use the single-frame Multiwfn approach:

```bash
cd "${PROJECT_ROOT}/params/ni_center/resp_charges"
cp ../orca_calcs/large_mk/large_mk.molden.input .

cat > multiwfn_commands.txt << 'EOF'
7
18
1
q
q
EOF

"${MULTIWFN_PATH}" large_mk.molden.input < multiwfn_commands.txt \
    > multiwfn_resp.log 2>&1
```

**CHECK:**

```bash
grep "Sum of charges" multiwfn_resp.log
# Expected: close to +1.000
```

**WHY:** Multi-frame RESP fitting produces charges that better represent the
electrostatic potential across the conformational landscape of the active site.
This is particularly important for flexible coordinating residues (His, Asp, KCX)
whose orientation changes during MD. Single-frame fitting can over-fit to one
conformation.

---

## Step 2.8 — Run MCPB.py Steps 2–4

**DO:**

```bash
cd "${PROJECT_ROOT}/params/ni_center/mcpb_input"

# Step 2: Extract force constants via Seminario method
MCPB.py -i urease.in -s 2 2>&1 | tee "${LOG_DIR}/mcpb_step2.log"

# Step 3: RESP charges (may need Multiwfn charges inserted first)
MCPB.py -i urease.in -s 3 2>&1 | tee "${LOG_DIR}/mcpb_step3.log"

# Step 4: Generate tleap input and library files
MCPB.py -i urease.in -s 4 2>&1 | tee "${LOG_DIR}/mcpb_step4.log"
```

**If Step 3 fails** (cannot read ORCA ESP data):
1. Manually update the mol2 files in `residue_mol2/` with the Multiwfn RESP charges from Step 2.7
2. Skip Step 3 and proceed to Step 4
3. MCPB.py Step 4 uses whatever charges are in the mol2 files

**CHECK:**

```bash
# Step 2 output: frcmod with force constants
echo "=== Force Constants ==="
grep "NI" urease.frcmod | head -20
# Expected: Ni-N bonds: k ≈ 50-120 kcal/mol/Å², r₀ ≈ 2.05 Å
# Expected: Ni-O bonds: k ≈ 60-130 kcal/mol/Å², r₀ ≈ 2.05 Å

# Step 4 output: library and tleap script
ls -la urease.lib urease_tleap.in
# Both should exist

# Verify charge sums
python3 << 'PYEOF'
"""Verify RESP charges sum to correct total for each residue."""
import re
from pathlib import Path

for mol2 in sorted(Path(".").glob("*.mol2")):
    charges = []
    in_atoms = False
    for line in mol2.read_text().splitlines():
        if "@<TRIPOS>ATOM" in line:
            in_atoms = True
            continue
        if in_atoms and line.startswith("@"):
            in_atoms = False
        if in_atoms and line.strip():
            parts = line.split()
            if len(parts) >= 9:
                charges.append(float(parts[8]))
    total = sum(charges)
    print(f"  {mol2.name:<15s} {len(charges):3d} atoms  charge = {total:+.4f}")
PYEOF
```

**Copy final parameters:**

```bash
cp urease.frcmod urease.lib urease_tleap.in *.mol2 \
   "${PROJECT_ROOT}/params/ni_center/final/"

echo "[DONE] Phase 2 complete."
echo "Parameters in: params/ni_center/final/"
echo "Next: Phase 3 (inhibitor parameterization)"
```

**WHY:** The final output is a complete set of files that define the Ni²⁺ bonded model for AMBER: force constants (frcmod), residue definitions (lib), and partial charges (mol2). These parameters maintain the correct octahedral/distorted-square-pyramidal Ni coordination geometry throughout MD by treating Ni–ligand interactions as explicit covalent bonds with DFT-derived spring constants.
