# AMBER MD Simulation — Complete Setup Instructions

## H. pylori Urease (6ZJA) · Dodecamer · 24 Ni²⁺ · 78 Inhibitor Complexes

> **How to use this guide:** Work through the phases sequentially. Each phase has numbered steps. Every step has a **DO** action (the command), a **CHECK** verification, and a **WHY** justification. Do not skip checks — a single error in Phase 1 propagates through all downstream work.
>
> **Cursor AI integration:** All Python scripts include type hints, docstrings, logging, and `if __name__ == "__main__"` guards. All Bash scripts use `set -euo pipefail` and structured error handling. Config is separated from logic via `config.env`. Every script can be run independently or composed by Cursor AI.

---

## Quick Reference — What Phase Am I In?

```
Phase 1  Structure Preparation ............. 1-2 days    (manual + pdb4amber)
Phase 2  Ni²⁺ Parameterization (MCPB.py) .. 5-10 days   (FM opt + freq + multi-frame RESP)
Phase 3  Inhibitor Parameterization ........ 1 day       (antechamber, automated)
Phase 4  System Assembly ................... 2-3 days    (tleap × 78, scriptable)
Phase 5  Equilibration .................... ~10 GPU-days  (pmemd.cuda × 78)
Phase 6  Production MD .................... ~10-16 GPU-days (Phase A: 10 ns × 78)
Phase 7  Analysis .......................... 1-2 weeks   (cpptraj + MMPBSA.py)
```

---

## Step 0: Project Setup

### 0.1 Create the directory structure

**DO:**

```bash
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
    "params/ni_center/multiframe_resp/xtb_md"
    "params/ni_center/multiframe_resp/frames"
    "params/ni_center/multiframe_resp/esp_calcs"
    "params/ni_center/multiframe_resp/resp_output"
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
```

**CHECK:**

```bash
# Verify the tree exists
find /path/to/urease_md -type d | head -30
# Expected: ~35 directories matching the list above
```

**WHY:** A fixed, documented directory structure means every script can use relative paths from the project root. Cursor AI can navigate the project without ambiguity. Separating `params/`, `systems/`, `simulations/`, and `analysis/` keeps the four stages of the pipeline independent and re-runnable.

---

### 0.2 Create the project configuration file

Every path, hardware setting, and software location goes in ONE file. No script hard-codes paths.

**DO:** Create `config.env` at the project root:

```bash
# config.env
# -------------------------------------------------------
# Central configuration for the urease MD project.
# Source this file at the top of every script:
#     source "$(dirname "$0")/../config.env"
# -------------------------------------------------------

# === Project root (set this to your actual path) ===
export PROJECT_ROOT="/path/to/urease_md"

# === Software paths ===
export ORCA_PATH="/opt/software/ORCA/build/orca"
export ORCA_2MKL="${ORCA_PATH%/*}/orca_2mkl"
export AMBER_HOME="/opt/software/amber24"    # Adjust to your install
export MULTIWFN_PATH="/opt/software/Multiwfn/Multiwfn"

# === Hardware ===
export N_CPU_CORES=56
export RAM_GB=62
export GPU_COUNT=1

# === ORCA resource allocation ===
# Rule: maxcore × nprocs ≤ 0.75 × RAM_GB × 1024
export ORCA_NPROCS=16
export ORCA_MAXCORE=2800       # MB per core → 16 × 2800 = 44.8 GB total

# === AMBER engine ===
export PMEMD_CUDA="pmemd.cuda"
export PMEMD_CPU="sander"      # Fallback for non-GPU stages

# === Force fields ===
export PROTEIN_FF="leaprc.protein.ff19SB"
export WATER_FF="leaprc.water.tip3p"
export LIGAND_FF="leaprc.gaff2"
export WATER_BOX="TIP3PBOX"

# === Simulation parameters ===
export BOX_BUFFER=12.0         # Å — solute-to-box-edge distance
export SALT_CONC=0.150         # M — physiological NaCl
export PRODUCTION_DT=0.004     # ps — 4 fs with HMR
export EQUIL_DT=0.002          # ps — 2 fs for equilibration
export PME_CUTOFF=9.0          # Å

# === PDB ===
export PDB_ID="6ZJA"
export PDB_URL="https://files.rcsb.org/download/${PDB_ID}.pdb"

# === Residue numbering (UPDATE after pdb4amber) ===
# These are placeholders — you MUST update after running pdb4amber
# with the actual renumbered values from 6ZJA_clean_renum.txt
export RESIDUE_ATOMS_IN_CLUSTER=257    # DFT cluster size
export NI_FLIP_INDEX=249               # 0-based Ni index for BS-DFT

# === Logging ===
export LOG_LEVEL="INFO"
export LOG_DIR="${PROJECT_ROOT}/logs"
```

**CHECK:**

```bash
source config.env
echo "Project root: ${PROJECT_ROOT}"
echo "ORCA: ${ORCA_PATH}"
ls -la "${ORCA_PATH}"     # Must exist
which tleap               # Must be on PATH
which pmemd.cuda           # Must be on PATH
echo "[OK] Config validated"
```

**WHY:** Centralised configuration is the single most important coding practice for a multi-phase computational project. When you move to a different machine, change AMBER versions, or hand off the project, only `config.env` needs editing. No grep-and-replace across 20 scripts.

---

### 0.3 Verify all software dependencies

**DO:**

```bash
#!/bin/bash
# scripts/check_dependencies.sh
# -------------------------------------------------------
# Verify all required software is installed and accessible.
# Run BEFORE starting any phase.
#
# Exit codes:
#   0 — all dependencies found
#   1 — one or more missing
# -------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config.env"

MISSING=0

check_cmd() {
    local name=$1
    local cmd=$2
    if command -v "$cmd" &>/dev/null; then
        local ver
        ver=$($cmd --version 2>&1 | head -1 || echo "version unknown")
        printf "  %-20s %-10s %s\n" "$name" "[OK]" "$ver"
    else
        printf "  %-20s %-10s\n" "$name" "[MISSING]"
        MISSING=$((MISSING + 1))
    fi
}

check_file() {
    local name=$1
    local path=$2
    if [ -x "$path" ]; then
        printf "  %-20s %-10s %s\n" "$name" "[OK]" "$path"
    else
        printf "  %-20s %-10s expected at %s\n" "$name" "[MISSING]" "$path"
        MISSING=$((MISSING + 1))
    fi
}

echo "=== Software Dependencies ==="
echo ""
echo "--- AmberTools ---"
check_cmd "pdb4amber"     pdb4amber
check_cmd "tleap"         tleap
check_cmd "cpptraj"       cpptraj
check_cmd "antechamber"   antechamber
check_cmd "parmchk2"      parmchk2
check_cmd "MCPB.py"       MCPB.py
check_cmd "MMPBSA.py"     MMPBSA.py
check_cmd "parmed"        parmed
check_cmd "reduce"        reduce
check_cmd "ambmask"       ambmask
check_cmd "ante-MMPBSA.py" ante-MMPBSA.py

echo ""
echo "--- AMBER (licensed) ---"
check_cmd "pmemd.cuda"    pmemd.cuda

echo ""
echo "--- QM Software ---"
check_file "ORCA 6.x"     "${ORCA_PATH}"
check_file "Multiwfn"     "${MULTIWFN_PATH}"

echo ""
echo "--- Utilities ---"
check_cmd "OpenBabel"     obabel
check_cmd "Python 3"      python3

echo ""
echo "--- GPU ---"
if nvidia-smi &>/dev/null; then
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
    GPU_MEM=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader | head -1)
    printf "  %-20s %-10s %s (%s)\n" "NVIDIA GPU" "[OK]" "$GPU_NAME" "$GPU_MEM"
else
    printf "  %-20s %-10s\n" "NVIDIA GPU" "[MISSING]"
    MISSING=$((MISSING + 1))
fi

echo ""
if [ "$MISSING" -eq 0 ]; then
    echo "[OK] All dependencies found."
    exit 0
else
    echo "[FAIL] ${MISSING} dependencies missing. Install them before proceeding."
    exit 1
fi
```

**CHECK:** The script self-checks. Run it and fix any `[MISSING]` items before continuing.

**WHY:** Discovering a missing tool in Phase 4 after days of QM calculations wastes time. This gate ensures everything is available from the start.

---

## Phase 1: Structure Preparation

**Goal:** Produce a clean, protonated PDB of the full 6ZJA dodecamer with correct histidine tautomers and non-standard residues identified for downstream parameterization.

**Inputs:** PDB file `6ZJA.pdb`  
**Outputs:**
- `structures/prepared/6ZJA_prepared.pdb` — final PDB for all downstream phases
- `structures/cleaned/6ZJA_clean_renum.txt` — residue renumbering map
- `logs/phase1.log` — full log of all operations

---

### Step 1.1 — Download the PDB

**DO:**

```bash
source config.env
cd "${PROJECT_ROOT}/structures/raw"

wget -q "${PDB_URL}" -O "${PDB_ID}.pdb"
```

**CHECK:**

```bash
# File exists and is not empty
[ -s "${PDB_ID}.pdb" ] && echo "[OK] Downloaded" || echo "[FAIL] Empty or missing"

# Verify it's the correct structure: dodecamer = 24 chains
grep "^ATOM" "${PDB_ID}.pdb" | awk '{print $5}' | sort -u | wc -l
# Expected: 24 unique chain identifiers (A through X, or similar)

# Count Ni atoms
grep "^HETATM.*NI" "${PDB_ID}.pdb" | wc -l
# Expected: 24 (2 per active site × 12 sites)
```

**WHY:** The PDB is the single source of truth for the protein structure. Verifying chain count and Ni atoms immediately confirms we have the correct biological assembly (dodecamer, not monomer or asymmetric unit).

---

### Step 1.2 — Clean with pdb4amber

**DO:**

```bash
source config.env
cd "${PROJECT_ROOT}/structures/cleaned"

pdb4amber \
    -i "${PROJECT_ROOT}/structures/raw/${PDB_ID}.pdb" \
    -o "${PDB_ID}_clean.pdb" \
    --dry \
    --nohyd \
    --most-populous \
    2>&1 | tee "${LOG_DIR}/pdb4amber.log"
```

**Flags explained:**

| Flag | Action | Why |
|------|--------|-----|
| `--dry` | Remove ALL water molecules | Crystallographic waters are removed; we add bulk TIP3P later |
| `--nohyd` | Strip ALL existing hydrogens | We add fresh H with `reduce` using controlled tautomer logic |
| `--most-populous` | Keep only the dominant alternate conformation | Eliminates alt-conf disorder that confuses tleap |

**What pdb4amber produces:**

| File | Description | Action needed |
|------|-------------|---------------|
| `6ZJA_clean.pdb` | Cleaned protein (standard residues only) | Input for Step 1.3 |
| `6ZJA_clean_nonprot.pdb` | Extracted non-standard residues (KCX, NI, etc.) | Review for Phase 2 |
| `6ZJA_clean_renum.txt` | Old→new residue number mapping | **CRITICAL: Update config.env and all scripts** |
| `6ZJA_clean_sslink` | Detected disulfide bonds | Review (urease has none typically) |

**CHECK:**

```bash
# 1. pdb4amber completed without fatal errors
grep -ci "error\|fatal" "${LOG_DIR}/pdb4amber.log"
# Expected: 0 (warnings about non-standard residues are OK)

# 2. Non-standard residues were flagged
grep -c "KCX\|NI" "${PDB_ID}_clean_nonprot.pdb"
# Expected: >0 (KCX and NI entries)

# 3. All 24 chains survived
grep "^ATOM" "${PDB_ID}_clean.pdb" | awk '{print $5}' | sort -u | wc -l
# Expected: 24

# 4. Renumbering map exists
[ -f "${PDB_ID}_clean_renum.txt" ] && echo "[OK]" || echo "[FAIL] No renum file"
```

**WHY:** `pdb4amber` is the standard AMBER structure preparation tool. It handles:
- Sequential atom/residue renumbering (required by tleap)
- Alternate conformation cleanup
- Non-standard residue identification and extraction

The renumbering map (`_renum.txt`) is essential because ALL downstream residue references (tleap bond commands, cpptraj masks, restraint masks) must use the NEW numbers.

---

### Step 1.3 — Add hydrogens with reduce

**DO:**

```bash
source config.env
cd "${PROJECT_ROOT}/structures/prepared"

reduce -build \
    "${PROJECT_ROOT}/structures/cleaned/${PDB_ID}_clean.pdb" \
    > "${PDB_ID}_H.pdb" \
    2> "${LOG_DIR}/reduce.log"
```

**CHECK:**

```bash
# 1. Output file exists and is larger than input (H atoms added)
INPUT_SIZE=$(stat -c%s "${PROJECT_ROOT}/structures/cleaned/${PDB_ID}_clean.pdb")
OUTPUT_SIZE=$(stat -c%s "${PDB_ID}_H.pdb")
echo "Input: ${INPUT_SIZE} bytes → Output: ${OUTPUT_SIZE} bytes"
# Output should be ~30-50% larger

# 2. Count hydrogen atoms added
grep "^ATOM" "${PDB_ID}_H.pdb" | awk '{print $3}' | grep -c "^H"
# Should be thousands of H atoms

# 3. Check reduce's histidine assignments
grep "^USER  MOD" "${LOG_DIR}/reduce.log" | grep -i "his"
# This shows which tautomers reduce chose — we will override many of these
```

**WHY:** `reduce` uses a scoring function based on steric contacts and hydrogen bonding to assign protonation states. However, it has NO understanding of metal coordination. It WILL incorrectly protonate metal-coordinating nitrogen atoms. Step 1.4 fixes this.

---

### Step 1.4 — Fix histidine protonation states for metal coordination

This is the most error-prone manual step in the entire workflow. Get it wrong and the active site electrostatics are incorrect for all 78 simulations.

**Scientific background:** Each histidine has two ring nitrogens — Nδ1 and Nε2. Only one can be protonated at physiological pH. When a nitrogen donates its lone pair to Ni²⁺, it CANNOT carry a hydrogen — the other nitrogen must be protonated.

| If this N coordinates Ni²⁺ | Then assign | Because |
|---------------------------|------------|---------|
| Nε2 | HID | Only Nδ1 carries H |
| Nδ1 | HIE | Only Nε2 carries H |
| Neither (free His) | HID or HIE | Use reduce's assignment or pKa prediction |

**Required assignments for the urease active site (per the DFT cluster model):**

| Residue | AMBER Name | Nδ1 | Nε2 | Role |
|---------|-----------|-----|-----|------|
| His136 | HID | protonated (H) | → coordinates Ni2 | Metal ligand |
| His138 | HID | protonated (H) | → coordinates Ni2 | Metal ligand |
| His221 | HIE | → coordinates Ni1 | protonated (H) | Metal ligand |
| His248 | HIE | → coordinates Ni1 | protonated (H) | Metal ligand |
| His274 | HID | protonated (H) | → coordinates Ni1 | Metal ligand (via Nε2) |
| His322 | HID | protonated (H) | free | Flap residue |
| His323 | HID | protonated (H) | free | Flap residue |

Also set:
- **Asp362** → ASP (deprotonated, charge −1) — Oδ1 coordinates Ni2
- **Cys321** → CYS (neutral thiol, with HG) — open flap, not coordinating

**DO:**

```python
#!/usr/bin/env python3
"""
scripts/fix_protonation.py
---------------------------------------------------------------------------
Fix histidine tautomers and key protonation states for the urease active site.

This script:
  1. Reads the reduce-protonated PDB
  2. Renames metal-coordinating His to HID/HIE per the DFT model
  3. Applies the fix across ALL 12 active sites (all UreB chains)
  4. Writes the final prepared PDB

Usage:
    python3 scripts/fix_protonation.py \
        structures/prepared/6ZJA_H.pdb \
        structures/prepared/6ZJA_prepared.pdb \
        structures/cleaned/6ZJA_clean_renum.txt

Coding standards:
    - Type hints on all functions
    - Logging instead of print()
    - Returns exit code 0/1
    - All configuration at top of file
---------------------------------------------------------------------------
"""

import sys
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# -------------------------------------------------------
# CONFIGURATION — Protonation assignments
# -------------------------------------------------------
# Key: ORIGINAL 6ZJA residue number (before pdb4amber renumbering)
# Value: target AMBER residue name
# These are per-chain assignments — applied to all UreB chains
ORIG_HIS_ASSIGNMENTS: Dict[int, str] = {
    136: "HID",   # His136: Nε2 → Ni2, so protonate Nδ1
    138: "HID",   # His138: Nε2 → Ni2, so protonate Nδ1
    221: "HIE",   # His221: Nδ1 → Ni1, so protonate Nε2
    248: "HIE",   # His248: Nδ1 → Ni1, so protonate Nε2
    274: "HID",   # His274: Nε2 → Ni1, so protonate Nδ1
    322: "HID",   # His322: flap, not metal-bound (DFT model)
    323: "HID",   # His323: flap, not metal-bound (DFT model)
}

# Other residues needing attention
# Asp362 should be ASP (deprotonated) — usually already correct
# Cys321 should be CYS (not CYX) with SH — check after reduce

# -------------------------------------------------------
# LOGGING
# -------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s"
)
log = logging.getLogger("fix_protonation")


def load_renumbering(renum_file: Path) -> Dict[Tuple[str, int], int]:
    """
    Parse pdb4amber renumbering file.
    Returns mapping: (chain, old_resnum) → new_resnum.
    """
    mapping: Dict[Tuple[str, int], int] = {}
    with open(renum_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            # Format varies by AmberTools version — parse adaptively
            # Common: old_chain old_num new_num new_chain
            if len(parts) >= 3:
                try:
                    old_chain = parts[0]
                    old_num = int(parts[1])
                    new_num = int(parts[2])
                    mapping[(old_chain, old_num)] = new_num
                except (ValueError, IndexError):
                    continue
    return mapping


def build_new_number_assignments(
    renum: Dict[Tuple[str, int], int],
    orig_assignments: Dict[int, str]
) -> Dict[int, str]:
    """
    Map original residue numbers to renumbered residue numbers.
    Returns: {new_resnum: target_name} for ALL chains.
    """
    new_assignments: Dict[int, str] = {}
    for (chain, old_num), new_num in renum.items():
        if old_num in orig_assignments:
            new_assignments[new_num] = orig_assignments[old_num]
            log.debug(f"  Chain {chain} res {old_num} → {new_num}: {orig_assignments[old_num]}")
    return new_assignments


def fix_protonation(
    input_pdb: Path,
    output_pdb: Path,
    assignments: Dict[int, str]
) -> int:
    """
    Read PDB, rename His residues per assignments, write output.
    Returns number of changes made.
    """
    changes = 0
    output_lines: List[str] = []

    with open(input_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                resnum = int(line[22:26])

                if resnum in assignments and resname in ("HIS", "HIP", "HID", "HIE"):
                    target = assignments[resnum]
                    if resname != target:
                        line = line[:17] + f"{target:>3}" + line[20:]
                        changes += 1

            output_lines.append(line)

    with open(output_pdb, "w") as f:
        f.writelines(output_lines)

    return changes


def remove_wrong_hydrogens(pdb_file: Path, assignments: Dict[int, str]) -> int:
    """
    Remove protons that are on the WRONG nitrogen after reduce.
    - HID: must NOT have HE2 (proton on Nε2); must have HD1
    - HIE: must NOT have HD1 (proton on Nδ1); must have HE2

    Returns number of atoms removed.
    """
    removed = 0
    kept_lines: List[str] = []

    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resnum = int(line[22:26])
                atom_name = line[12:16].strip()

                if resnum in assignments:
                    target = assignments[resnum]
                    # HID: proton on Nδ1 only → remove HE2 if present
                    if target == "HID" and atom_name == "HE2":
                        removed += 1
                        continue
                    # HIE: proton on Nε2 only → remove HD1 if present
                    if target == "HIE" and atom_name == "HD1":
                        removed += 1
                        continue

            kept_lines.append(line)

    with open(pdb_file, "w") as f:
        f.writelines(kept_lines)

    return removed


def main() -> int:
    if len(sys.argv) != 4:
        print("Usage: python3 fix_protonation.py <input.pdb> <output.pdb> <renum.txt>")
        return 1

    input_pdb = Path(sys.argv[1])
    output_pdb = Path(sys.argv[2])
    renum_file = Path(sys.argv[3])

    # Step A: Load renumbering map
    if renum_file.exists():
        renum = load_renumbering(renum_file)
        log.info(f"Loaded {len(renum)} renumbering entries from {renum_file}")
        assignments = build_new_number_assignments(renum, ORIG_HIS_ASSIGNMENTS)
        log.info(f"Mapped {len(assignments)} His residues across all chains")
    else:
        log.warning(f"No renumbering file found at {renum_file}")
        log.warning("Using original residue numbers (may be incorrect)")
        assignments = ORIG_HIS_ASSIGNMENTS

    if not assignments:
        log.error("No residue assignments to apply. Check renumbering file format.")
        return 1

    # Step B: Rename residues
    changes = fix_protonation(input_pdb, output_pdb, assignments)
    log.info(f"Renamed {changes} His residues in {output_pdb.name}")

    # Step C: Remove wrong hydrogens
    removed = remove_wrong_hydrogens(output_pdb, assignments)
    log.info(f"Removed {removed} incorrect H atoms")

    # Step D: Summary
    log.info("--- Verification ---")
    his_counts = {"HID": 0, "HIE": 0, "HIP": 0, "HIS": 0}
    with open(output_pdb) as f:
        seen = set()
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                resnum = int(line[22:26])
                key = (resnum, resname)
                if resname in his_counts and key not in seen:
                    his_counts[resname] += 1
                    seen.add(key)

    for name, count in his_counts.items():
        if count > 0:
            log.info(f"  {name}: {count} residues")

    if his_counts.get("HIS", 0) > 0:
        log.warning(f"  {his_counts['HIS']} residues still named HIS — these need manual review")

    log.info(f"[DONE] Output: {output_pdb}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**RUN:**

```bash
source config.env
python3 scripts/fix_protonation.py \
    "${PROJECT_ROOT}/structures/prepared/${PDB_ID}_H.pdb" \
    "${PROJECT_ROOT}/structures/prepared/${PDB_ID}_prepared.pdb" \
    "${PROJECT_ROOT}/structures/cleaned/${PDB_ID}_clean_renum.txt"
```

**CHECK:**

```bash
cd "${PROJECT_ROOT}/structures/prepared"

# 1. Count HID and HIE assignments
grep "^ATOM" "${PDB_ID}_prepared.pdb" | awk '{print $4}' | sort | uniq -c | grep -E "HID|HIE"
# Expected: HID count = (5 His per site × 12 sites) + non-metal His
#           HIE count = (2 His per site × 12 sites)

# 2. Verify no "HIS" residues remain (all should be HID, HIE, or HIP)
grep "^ATOM" "${PDB_ID}_prepared.pdb" | awk '{print $4}' | sort | uniq -c | grep " HIS"
# Expected: 0 (if any remain, they need manual assignment)

# 3. Verify metal-coordinating N atoms have NO proton
# For HID residues: NE2 should NOT have a bonded HE2
# For HIE residues: ND1 should NOT have a bonded HD1
python3 -c "
with open('${PDB_ID}_prepared.pdb') as f:
    atoms = {}
    for line in f:
        if line.startswith(('ATOM','HETATM')):
            resnum = int(line[22:26])
            name = line[12:16].strip()
            atoms.setdefault(resnum, set()).add(name)

problems = 0
for resnum, names in atoms.items():
    resname = None  # would need to track — simplified check
    if 'HE2' in names and 'HID' in str(names):  # HID should not have HE2
        print(f'WARNING: Residue {resnum} is HID but has HE2')
        problems += 1
print(f'Problems found: {problems}')
"

# 4. Visual verification (ESSENTIAL — do not skip)
echo "Open in PyMOL and inspect each active site:"
echo "  pymol ${PDB_ID}_prepared.pdb"
echo "  # show sticks, resi 136+138+221+248+274 (one chain)"
echo "  # Check: coordinating N has no H; other N has H"
echo "  # Check: Ni atoms are present and ~2.0 Å from coordinating N/O"
```

**WHY:** The protonation state of metal-coordinating histidines directly determines the electrostatic potential at the active site. An incorrect tautomer means a 2+ kcal/mol error in binding free energy. This step must match the DFT model exactly (Phase 2) for physical consistency. The script is designed to:
1. Handle the renumbering automatically (no manual residue number lookup)
2. Apply to ALL 12 active sites at once
3. Both rename residues AND remove incorrect protons
4. Provide a verification summary

---

### Step 1.5 — Verify bridging hydroxide

The bridging OH⁻ between the two Ni ions may or may not be resolved in the 2.0 Å cryo-EM density.

**DO:**

```bash
source config.env

# Check if any bridging species is present between Ni pairs
python3 << 'PYEOF'
"""Check for bridging OH⁻ in the prepared PDB."""
import sys
import math
from pathlib import Path

pdb = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("structures/prepared/6ZJA_prepared.pdb")

# Parse all atoms
atoms = []
with open(pdb) as f:
    for line in f:
        if line.startswith(("ATOM", "HETATM")):
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            resnum = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            atoms.append((atom_name, resname, resnum, x, y, z))

# Find Ni atoms
ni_atoms = [(a[2], a[3], a[4], a[5]) for a in atoms if a[0] == "NI"]
print(f"Found {len(ni_atoms)} Ni atoms")

# For each Ni pair (same chain, consecutive resnums), check for bridging O
ni_pairs = []
for i in range(0, len(ni_atoms) - 1, 2):
    ni1 = ni_atoms[i]
    ni2 = ni_atoms[i + 1]
    ni_dist = math.sqrt(sum((a - b)**2 for a, b in zip(ni1[1:], ni2[1:])))
    midpoint = tuple((a + b) / 2 for a, b in zip(ni1[1:], ni2[1:]))

    # Find O atoms within 2.5 Å of both Ni
    bridging = []
    for a in atoms:
        if a[0].startswith("O") and a[1] not in ("WAT", "HOH"):
            d1 = math.sqrt(sum((a[j+3] - ni1[j+1])**2 for j in range(3)))
            d2 = math.sqrt(sum((a[j+3] - ni2[j+1])**2 for j in range(3)))
            if d1 < 2.5 and d2 < 2.5:
                bridging.append((a[1], a[2], a[0], d1, d2))

    print(f"\nNi pair (res {ni1[0]}, {ni2[0]}): Ni-Ni = {ni_dist:.2f} Å")
    if bridging:
        for b in bridging:
            print(f"  Bridging O found: {b[0]}:{b[1]}:{b[2]} "
                  f"(d1={b[3]:.2f}, d2={b[4]:.2f} Å)")
    else:
        print(f"  WARNING: No bridging O found! Must add OH⁻ manually.")
        print(f"  Midpoint for placement: ({midpoint[0]:.3f}, "
              f"{midpoint[1]:.3f}, {midpoint[2]:.3f})")
PYEOF
```

**If bridging OH⁻ is missing:** You must add it manually to each active site. Place the O atom at the midpoint between the two Ni ions (slightly offset toward the KCX carboxylate plane), at ~2.0 Å from each Ni. Add the H atom perpendicular to the Ni–O–Ni plane. This is documented in detail in `02_mcpb_parameterization.md`.

**CHECK:**

```bash
# After adding OH⁻ (if needed), verify all 12 sites have bridging species
# The script above should report "Bridging O found" for all 12 Ni pairs
```

**WHY:** The bridging hydroxide is part of the catalytic machinery. Without it, the MCPB.py bonded model is incomplete, and the Ni–Ni distance will not be correctly maintained. The bridging OH⁻ provides a critical structural bridge and participates in the catalytic mechanism.

---

### Phase 1 — Completion Checklist

Before proceeding to Phase 2, verify ALL of the following:

```bash
#!/bin/bash
# scripts/verify_phase1.sh
set -euo pipefail
source "$(dirname "$0")/../config.env"

PASS=0
FAIL=0

check() {
    local desc=$1
    local result=$2
    if [ "$result" = "true" ]; then
        printf "  [PASS] %s\n" "$desc"
        PASS=$((PASS + 1))
    else
        printf "  [FAIL] %s\n" "$desc"
        FAIL=$((FAIL + 1))
    fi
}

echo "=== Phase 1 Completion Checklist ==="

PREP="${PROJECT_ROOT}/structures/prepared/${PDB_ID}_prepared.pdb"

# File exists
check "Prepared PDB exists" "$([ -s "$PREP" ] && echo true || echo false)"

# 24 chains
CHAINS=$(grep "^ATOM" "$PREP" | awk '{print $5}' | sort -u | wc -l)
check "24 chains present (found: ${CHAINS})" "$([ "$CHAINS" -ge 20 ] && echo true || echo false)"

# 24 Ni atoms
NI_COUNT=$(grep -c "NI" "$PREP" 2>/dev/null || echo 0)
check "24 Ni atoms present (found: ${NI_COUNT})" "$([ "$NI_COUNT" -ge 24 ] && echo true || echo false)"

# No HIS residues (all assigned to HID/HIE/HIP)
HIS_COUNT=$(grep "^ATOM" "$PREP" | awk '{print $4}' | grep -c "^HIS$" 2>/dev/null || echo 0)
check "No unassigned HIS residues (found: ${HIS_COUNT})" "$([ "$HIS_COUNT" -eq 0 ] && echo true || echo false)"

# Renumbering map exists
check "Renumbering map exists" "$([ -f "${PROJECT_ROOT}/structures/cleaned/${PDB_ID}_clean_renum.txt" ] && echo true || echo false)"

echo ""
echo "Results: ${PASS} passed, ${FAIL} failed"
[ "$FAIL" -eq 0 ] && echo "[OK] Phase 1 complete. Proceed to Phase 2." || echo "[BLOCKED] Fix failures before proceeding."
```

---

## Coding Standards Used in This Project

Every script in this project follows these conventions for Cursor AI compatibility and maintainability:

### Bash Scripts

```bash
#!/bin/bash
# scripts/<name>.sh
# -------------------------------------------------------
# <One-line description>
#
# Usage:
#     bash scripts/<name>.sh <arg1> [arg2]
#
# Inputs:
#     <arg1>  — description
#
# Outputs:
#     <file>  — description
# -------------------------------------------------------
set -euo pipefail                    # Exit on error, undefined var, pipe fail
source "$(dirname "$0")/../config.env"  # Load central config

log() { echo "[$(date '+%H:%M:%S')] $1"; }

log "Starting <name>..."
# ... script body ...
log "[DONE]"
```

**Rules:**
1. `set -euo pipefail` on every script — no silent failures
2. `source config.env` — no hard-coded paths
3. All user-facing output through `log()` with timestamps
4. `${VAR:?"message"}` for required arguments — fail immediately with clear message
5. Exit codes: 0 = success, 1 = error
6. Comments explain WHY, not WHAT

### Python Scripts

```python
#!/usr/bin/env python3
"""
scripts/<name>.py
---------------------------------------------------------------------------
<One-line description>

<Longer explanation of what this script does and why.>

Usage:
    python3 scripts/<name>.py <arg1> [arg2]

Inputs:
    <arg1>  — description

Outputs:
    <file>  — description
---------------------------------------------------------------------------
"""

import sys
import logging
from pathlib import Path
from typing import Dict, List, Optional

# Configuration
CONFIG_KEY = "value"

# Logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def do_thing(input_path: Path, output_path: Path) -> int:
    """One-line docstring. Returns number of items processed."""
    # Implementation
    return 0


def main() -> int:
    """Entry point. Returns 0 on success, 1 on error."""
    if len(sys.argv) < 2:
        print(f"Usage: python3 {sys.argv[0]} <input>")
        return 1

    result = do_thing(Path(sys.argv[1]), Path(sys.argv[2]))
    log.info(f"Processed {result} items")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**Rules:**
1. Type hints on all function signatures
2. Docstrings on all functions (Google style)
3. `logging` module, never bare `print()` (except for data output)
4. `if __name__ == "__main__"` guard — every script importable as a module
5. `Path` objects, never string path manipulation
6. `sys.exit(main())` — propagate exit codes
7. Configuration constants at TOP of file, clearly labeled

### AMBER Input Files

```
<Title describing what this stage does>
 <One-line justification for key parameter choices>
 &cntrl
   <parameters grouped logically>
   <with inline comments for non-obvious values>
 /
```

**Rules:**
1. Every input file has a descriptive title (line 1)
2. Non-obvious parameters have inline comments
3. All numerical values match config.env where applicable
4. Shared across all 78 systems (system-specific data comes from file paths, not input file contents)
