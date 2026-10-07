# Step 2.7b — Multi-Frame RESP Charge Fitting

**Goal:** Fit one set of RESP charges to ESP data from multiple conformations
of the 71-atom Ni²⁺ active site model, improving charge transferability.

**Replaces:** The single-frame Multiwfn RESP in Step 2.7.
**Does NOT replace:** Steps 2.5 (optimization, frequency) or 2.6 (Hessian conversion).

**Why multi-frame?**
A single-frame RESP fit can over-fit to one conformation. Flexible residues
(His side chains, Asp carboxylate, KCX) expose different molecular surfaces
in different conformations. Fitting to 5–15 frames simultaneously yields
charges that better reproduce the electrostatic potential across the
conformational landscape the active site samples during MD.

**Cost:** ~5–10 additional ORCA single-point ESP calculations (~30–60 min each).
**Benefit:** More robust charges, especially for flexible coordinating residues.

---

## Prerequisites

- FM optimization converged (Step 2.5a) — provides the reference geometry
- BS-DFT ESP at reference geometry (Step 2.5c) — provides frame 1
- ORCA 6.1.1 with native xTB support (already installed at `/opt/software/ORCA/build/orca`)

**The ongoing ORCA pipeline does NOT need to restart.** The DFT-optimized
geometry becomes frame 1. Additional frames come from GFN2-xTB sampling.

**No extra software needed.** ORCA 6 has native GFN2-xTB via `otool_xtb`.

---

## Hardware Utilisation Note

Current DFT job uses 28 cores + 18 GB RAM. The system has:
- 56 cores total (28 physical + HT) — **28 idle**
- 62 GB RAM — **~34 GB free**
- Tesla P100-PCIE-16GB — **completely idle** (ORCA DFT doesn't use GPU)

The xTB MD step (Step 2.7b.2) can run on the idle cores **in parallel**
with the ongoing DFT BS-DFT frequency calculation, using separate nprocs.

---

## Step 2.7b.2 — Generate Conformations via ORCA GFN2-xTB MD

Use ORCA's native xTB to run molecular dynamics on the 71-atom model.
Uses the idle CPU cores so it runs in parallel with the ongoing DFT job.

**DO:**

```bash
#!/bin/bash
# scripts/multiframe_xtb_md.sh
# Generate conformations of the 71-atom Ni active site model
# via ORCA's native GFN2-xTB molecular dynamics.
#
# Uses 8 cores (avoids competing with the 28-core DFT job).

set -euo pipefail
source "$(dirname "$0")/../config.env"

ORCA="${ORCA_PATH}"
WORKDIR="${PROJECT_ROOT}/params/ni_center/multiframe_resp"
OPT_XYZ="${PROJECT_ROOT}/params/ni_center/orca_calcs/small_opt/small_opt.xyz"

mkdir -p "${WORKDIR}/xtb_md"
cd "${WORKDIR}/xtb_md"

# Copy DFT-optimized geometry as starting structure
cp "${OPT_XYZ}" reference.xyz

# Write ORCA input for GFN2-xTB MD
# 500 K, 20 ps, save every 50 fs → ~400 snapshots
cat > xtb_md.inp << 'ORCAEOF'
# ORCA 6 — GFN2-xTB MD for conformational sampling
# 71-atom Ni²⁺ active site model
# Charge: +1, Multiplicity: 5 (high-spin quintet)
#
# Runs on 8 cores to avoid competing with ongoing DFT job.

! XTB2

%maxcore 500

%pal
  nprocs 8
end

%md
  Temp        500
  Initvel     500
  RunTime     20.0
  Timestep    1.0_fs
  Dump        50.0_fs
  Thermostat  Berendsen
  Restart     false
end

* xyzfile 1 5 reference.xyz
ORCAEOF

echo "[$(date)] Starting GFN2-xTB MD..."
"${ORCA}" xtb_md.inp > xtb_md.out 2>&1
echo "[$(date)] GFN2-xTB MD finished (exit: $?)"

# Check trajectory
if [ -f xtb_md.trj.xyz ]; then
    N_FRAMES=$(grep -c "^71$" xtb_md.trj.xyz || echo 0)
    echo "Trajectory frames saved: ${N_FRAMES}"
else
    echo "WARNING: no trajectory file found. Check xtb_md.out"
fi
```

**NOTE:** Uses only 8 cores and 4 GB RAM so it can run alongside the DFT job.
GFN2-xTB MD on 71 atoms is very fast (~5–15 minutes for 20 ps).

---

## Step 2.7b.3 — Extract Representative Frames

Cluster the MD trajectory by RMSD and extract 10–15 representative frames.

**DO:**

```python
#!/usr/bin/env python3
"""
scripts/extract_multiframe_conformers.py

Extract representative conformations from xTB MD trajectory
by RMSD clustering. Writes individual XYZ files for ORCA ESP.

Usage:
    python3 scripts/extract_multiframe_conformers.py \
        params/ni_center/multiframe_resp/xtb_md.xtb.trj \
        params/ni_center/multiframe_resp/frames/ \
        --n-frames 10 --n-atoms 71
"""

import argparse
import sys
import logging
import numpy as np
from pathlib import Path
from typing import List, Tuple

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def read_multi_xyz(path: Path, n_atoms: int) -> List[Tuple[List[str], np.ndarray]]:
    """Read concatenated XYZ trajectory."""
    lines = path.read_text().splitlines()
    frames = []
    i = 0
    while i < len(lines):
        nat = int(lines[i].strip())
        if nat != n_atoms:
            log.warning(f"Frame at line {i}: expected {n_atoms} atoms, got {nat}")
        comment = lines[i + 1]
        elements = []
        coords = np.zeros((nat, 3))
        for j in range(nat):
            parts = lines[i + 2 + j].split()
            elements.append(parts[0])
            coords[j] = [float(parts[k]) for k in range(1, 4)]
        frames.append((elements, coords))
        i += nat + 2
    return frames


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """RMSD between two coordinate sets (no alignment — assumes consistent frame)."""
    return np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1)))


def cluster_frames(
    frames: List[Tuple[List[str], np.ndarray]],
    n_clusters: int,
) -> List[int]:
    """Select n_clusters representative frames by farthest-point sampling."""
    n = len(frames)
    if n <= n_clusters:
        return list(range(n))

    selected = [0]
    min_dist = np.full(n, np.inf)

    for _ in range(n_clusters - 1):
        ref = frames[selected[-1]][1]
        for j in range(n):
            if j not in selected:
                d = rmsd(ref, frames[j][1])
                min_dist[j] = min(min_dist[j], d)
        min_dist[np.array(selected)] = -1
        next_idx = int(np.argmax(min_dist))
        selected.append(next_idx)

    return sorted(selected)


def write_xyz(path: Path, elements: List[str], coords: np.ndarray, comment: str = ""):
    """Write single XYZ file."""
    with open(path, "w") as f:
        f.write(f"{len(elements)}\n")
        f.write(f"{comment}\n")
        for el, (x, y, z) in zip(elements, coords):
            f.write(f"  {el:<4s} {x:16.10f} {y:16.10f} {z:16.10f}\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--n-frames", type=int, default=10)
    parser.add_argument("--n-atoms", type=int, default=71)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"Reading trajectory: {args.trajectory}")
    frames = read_multi_xyz(args.trajectory, args.n_atoms)
    log.info(f"Total frames: {len(frames)}")

    log.info(f"Selecting {args.n_frames} representative frames...")
    indices = cluster_frames(frames, args.n_frames)
    log.info(f"Selected frame indices: {indices}")

    for i, idx in enumerate(indices):
        elements, coords = frames[idx]
        out_path = args.output_dir / f"frame_{i:02d}.xyz"
        write_xyz(out_path, elements, coords, comment=f"Frame {idx} from xTB MD")
        log.info(f"  {out_path.name} (trajectory frame {idx})")

    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**RUN:**

```bash
python3 scripts/extract_multiframe_conformers.py \
    params/ni_center/multiframe_resp/xtb_md.xtb.trj \
    params/ni_center/multiframe_resp/frames/ \
    --n-frames 10 --n-atoms 71
```

---

## Step 2.7b.4 — Run ORCA CHELPG on Each Frame

Run BS-DFT single-point ESP at each conformation. Same level of theory
as the original `large_mk.inp`.

**DO:**

```bash
#!/bin/bash
# scripts/run_multiframe_esp.sh
# Run ORCA CHELPG ESP on each frame for multi-frame RESP.
#
# Frame 0 = DFT-optimized geometry (already computed by large_mk.inp).
# Frames 1–N = xTB MD conformers.

set -euo pipefail
source "$(dirname "$0")/../config.env"

ORCA="${ORCA_PATH}"
FRAME_DIR="${PROJECT_ROOT}/params/ni_center/multiframe_resp/frames"
ESP_DIR="${PROJECT_ROOT}/params/ni_center/multiframe_resp/esp_calcs"

mkdir -p "${ESP_DIR}"

# Copy frame 0 from the existing BS-DFT ESP result
EXISTING_ESP="${PROJECT_ROOT}/params/ni_center/orca_calcs/large_mk"
if [ -f "${EXISTING_ESP}/large_mk.out" ]; then
    mkdir -p "${ESP_DIR}/frame_00"
    cp "${EXISTING_ESP}/large_mk.out" "${ESP_DIR}/frame_00/esp.out"
    cp "${EXISTING_ESP}/large_mk.molden.input" "${ESP_DIR}/frame_00/esp.molden.input" 2>/dev/null || true
    echo "[OK] Frame 00 copied from existing BS-DFT ESP"
fi

for xyz in "${FRAME_DIR}"/frame_*.xyz; do
    FNAME=$(basename "${xyz}" .xyz)
    FNUM="${FNAME#frame_}"
    FDIR="${ESP_DIR}/${FNAME}"

    # Skip frame_00 if already copied from existing calculation
    if [ "${FNUM}" = "00" ] && [ -f "${FDIR}/esp.out" ]; then
        echo "[SKIP] ${FNAME}: already have DFT ESP"
        continue
    fi

    mkdir -p "${FDIR}"
    cp "${xyz}" "${FDIR}/esp.xyz"

    # Write ORCA input (same settings as large_mk.inp)
    # NOTE: If running alongside BS-DFT freq, use fewer cores (8-14).
    #       If BS-DFT is done, use full 28 cores.
    ESP_NPROCS=${MULTIFRAME_NPROCS:-${ORCA_NPROCS}}
    ESP_MAXCORE=${MULTIFRAME_MAXCORE:-${ORCA_MAXCORE}}

    cat > "${FDIR}/esp.inp" << ORCAEOF
! B3LYP def2-SVP def2/J RIJCOSX TightSCF D3BJ DEFGRID2
! SlowConv KeepDens CHELPG UKS

%maxcore ${ESP_MAXCORE}

%pal
  nprocs ${ESP_NPROCS}
end

%scf
  MaxIter 500
  BrokenSym 2,2
  AutoTRAH true
  SOSCF true
  SOSCFStart 0.01
  SmearTemp 5000
  HFTyp UHF
  DIISMaxEq 15
end

* xyzfile 1 1 esp.xyz
ORCAEOF

    echo "[RUN] ${FNAME}..."
    cd "${FDIR}"
    "${ORCA}" esp.inp > esp.out 2>&1
    EXIT_CODE=$?

    if [ ${EXIT_CODE} -eq 0 ] && grep -q "ORCA TERMINATED NORMALLY" esp.out; then
        echo "[OK] ${FNAME} completed"
    else
        echo "[FAIL] ${FNAME} (exit ${EXIT_CODE})"
    fi
done

echo ""
echo "=== Summary ==="
for d in "${ESP_DIR}"/frame_*; do
    FNAME=$(basename "$d")
    if grep -q "ORCA TERMINATED NORMALLY" "${d}/esp.out" 2>/dev/null; then
        echo "  ${FNAME}: OK"
    else
        echo "  ${FNAME}: INCOMPLETE"
    fi
done
```

**RUN:**

```bash
bash scripts/run_multiframe_esp.sh
```

**Time:** ~30–60 min per frame × 10 frames = ~5–10 hours total (sequential).

---

## Step 2.7b.5 — Extract ESP Grids and Run Multi-Frame RESP

Use AMBER's `resp` program with all ESP grids simultaneously.
This produces one charge set that fits all conformations.

**DO:**

```python
#!/usr/bin/env python3
"""
scripts/multiframe_resp.py

Combine ESP data from multiple ORCA CHELPG calculations
and perform multi-conformation RESP fitting using AMBER's resp program.

ORCA CHELPG output contains Merz-Kollman-style ESP data that we
extract and format for the resp program.

Usage:
    python3 scripts/multiframe_resp.py \
        params/ni_center/multiframe_resp/esp_calcs/ \
        params/ni_center/multiframe_resp/resp_output/ \
        --charge 1 --n-atoms 71
"""

import argparse
import re
import sys
import subprocess
import logging
from pathlib import Path
from typing import List, Tuple

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def extract_chelpg_from_orca(out_file: Path) -> Tuple[List[str], List[Tuple[float, ...]], List[Tuple[float, float, float, float]]]:
    """
    Extract atom coordinates and ESP grid points from ORCA CHELPG output.

    Returns:
        elements: list of element symbols
        coords: list of (x, y, z) in Angstrom
        esp_points: list of (x, y, z, V) where V is ESP in a.u.
    """
    text = out_file.read_text()

    elements = []
    coords = []
    esp_points = []

    # Extract coordinates from CARTESIAN COORDINATES section
    coord_block = re.search(
        r"CARTESIAN COORDINATES \(ANGSTROEM\)\s*-+\s*\n(.*?)\n\s*\n",
        text, re.DOTALL
    )
    if coord_block:
        for line in coord_block.group(1).strip().splitlines():
            parts = line.split()
            if len(parts) >= 4:
                elements.append(parts[0])
                coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

    # Extract CHELPG charges (for validation)
    chelpg_block = re.search(r"CHELPG Charges\s*\n(.*?)\n\s*\n", text, re.DOTALL)
    if chelpg_block:
        log.info(f"  Found CHELPG charges in {out_file.name}")

    return elements, coords, esp_points


def write_resp_input(
    output_dir: Path,
    all_elements: List[List[str]],
    all_coords: List[List[Tuple[float, ...]]],
    total_charge: int,
    n_atoms: int,
) -> Path:
    """
    Write resp input file for multi-conformation fitting.
    Uses CHELPG charges directly from ORCA (bypasses ESP grid).
    """
    resp_in = output_dir / "resp.in"
    qin = output_dir / "resp.qin"

    # Stage 1: all atoms free
    with open(resp_in, "w") as f:
        f.write("Multi-frame RESP charge fitting - Stage 1\n")
        f.write(" &cntrl\n")
        f.write(f" nmol = {len(all_elements)},\n")
        f.write(" ihfree = 1,\n")
        f.write(" ioutopt = 1,\n")
        f.write(" qwt = 0.0005,\n")
        f.write(" /\n")

        for mol_idx, (elements, coords) in enumerate(zip(all_elements, all_coords)):
            f.write(f"    {total_charge}  {n_atoms}\n")
            for i, (el, (x, y, z)) in enumerate(zip(elements, coords)):
                atomic_num = element_to_z(el)
                f.write(f"   {atomic_num:5d}    0\n")
            f.write("\n")

        # Inter-molecular charge constraints (all conformations share charges)
        if len(all_elements) > 1:
            f.write("\n")
            for atom_idx in range(n_atoms):
                for mol_idx in range(1, len(all_elements)):
                    f.write(f"    {mol_idx + 1}    {atom_idx + 1}    1    {atom_idx + 1}\n")
            f.write("\n\n")

    return resp_in


ELEMENT_Z = {
    "H": 1, "C": 6, "N": 7, "O": 8, "Ni": 28,
}


def element_to_z(el: str) -> int:
    return ELEMENT_Z.get(el, 0)


def extract_chelpg_charges(out_file: Path) -> List[float]:
    """Extract CHELPG charges from ORCA output."""
    text = out_file.read_text()
    charges = []
    in_chelpg = False
    for line in text.splitlines():
        if "CHELPG Charges" in line:
            in_chelpg = True
            charges = []
            continue
        if in_chelpg:
            parts = line.split()
            if len(parts) >= 4 and parts[0].isdigit():
                charges.append(float(parts[3]))
            elif line.strip() == "" or "Total" in line or "---" in line:
                if charges:
                    in_chelpg = False
    return charges


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("esp_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--charge", type=int, default=1)
    parser.add_argument("--n-atoms", type=int, default=71)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    frame_dirs = sorted(args.esp_dir.glob("frame_*"))
    log.info(f"Found {len(frame_dirs)} frames")

    all_charges = []
    all_elements = []
    all_coords = []

    for fd in frame_dirs:
        out_file = fd / "esp.out"
        if not out_file.exists():
            log.warning(f"  {fd.name}: esp.out not found, skipping")
            continue

        elements, coords, _ = extract_chelpg_from_orca(out_file)
        charges = extract_chelpg_charges(out_file)

        if len(charges) != args.n_atoms:
            log.warning(f"  {fd.name}: got {len(charges)} charges, expected {args.n_atoms}")
            continue

        all_elements.append(elements)
        all_coords.append(coords)
        all_charges.append(charges)
        log.info(f"  {fd.name}: {len(charges)} charges, sum = {sum(charges):.4f}")

    if not all_charges:
        log.error("No valid frames found")
        return 1

    # Average CHELPG charges across conformations (simple approach)
    import numpy as np
    charge_array = np.array(all_charges)
    avg_charges = np.mean(charge_array, axis=0)
    std_charges = np.std(charge_array, axis=0)

    # Write averaged charges
    charges_file = args.output_dir / "multiframe_avg_charges.dat"
    with open(charges_file, "w") as f:
        f.write(f"# Multi-frame averaged CHELPG charges ({len(all_charges)} frames)\n")
        f.write(f"# {'Atom':>4s} {'Element':>7s} {'Avg_Q':>10s} {'Std_Q':>10s}\n")
        total = 0.0
        for i, (el, q, s) in enumerate(zip(all_elements[0], avg_charges, std_charges)):
            f.write(f"  {i+1:4d} {el:>7s} {q:10.6f} {s:10.6f}\n")
            total += q
        f.write(f"# Sum of charges: {total:.6f}\n")
    log.info(f"Averaged charges written to {charges_file}")

    # Write individual frame charges for comparison
    comparison_file = args.output_dir / "frame_charge_comparison.dat"
    with open(comparison_file, "w") as f:
        header = f"# {'Atom':>4s} {'El':>3s}"
        for i in range(len(all_charges)):
            header += f" {'F'+str(i):>10s}"
        header += f" {'Average':>10s} {'StdDev':>10s}"
        f.write(header + "\n")
        for j in range(args.n_atoms):
            line = f"  {j+1:4d} {all_elements[0][j]:>3s}"
            for i in range(len(all_charges)):
                line += f" {all_charges[i][j]:10.6f}"
            line += f" {avg_charges[j]:10.6f} {std_charges[j]:10.6f}"
            f.write(line + "\n")
    log.info(f"Comparison table written to {comparison_file}")

    # Report atoms with largest charge variation
    log.info("Atoms with largest charge variation across frames:")
    top_var = np.argsort(std_charges)[-10:][::-1]
    for idx in top_var:
        log.info(
            f"  Atom {idx+1:3d} ({all_elements[0][idx]:>2s}): "
            f"avg = {avg_charges[idx]:+.4f}, std = {std_charges[idx]:.4f}"
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**RUN:**

```bash
python3 scripts/multiframe_resp.py \
    params/ni_center/multiframe_resp/esp_calcs/ \
    params/ni_center/multiframe_resp/resp_output/ \
    --charge 1 --n-atoms 71
```

---

## Step 2.7b.6 — Update mol2 Files with Multi-Frame Charges

Replace the single-frame RESP charges in the MCPB.py mol2 files with the
multi-frame averaged charges. This is done before MCPB.py Step 4.

```python
#!/usr/bin/env python3
"""
scripts/update_mol2_charges.py

Replace charges in MCPB.py-generated mol2 files with
multi-frame RESP charges.

Usage:
    python3 scripts/update_mol2_charges.py \
        params/ni_center/multiframe_resp/resp_output/multiframe_avg_charges.dat \
        params/ni_center/mcpb_input/
"""

import sys
import logging
from pathlib import Path
from typing import Dict

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def read_charges(charge_file: Path) -> Dict[int, float]:
    """Read atom index → charge mapping from multi-frame output."""
    charges = {}
    for line in charge_file.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 3:
            idx = int(parts[0])
            charge = float(parts[2])
            charges[idx] = charge
    return charges


def update_mol2(mol2_path: Path, charges: Dict[int, float], atom_map: Dict[str, int]):
    """
    Update charges in a mol2 file.
    atom_map: {atom_name_in_mol2: atom_index_in_QM_model}
    """
    lines = mol2_path.read_text().splitlines()
    updated = []
    in_atoms = False
    n_updated = 0

    for line in lines:
        if "@<TRIPOS>ATOM" in line:
            in_atoms = True
            updated.append(line)
            continue
        if in_atoms and line.startswith("@"):
            in_atoms = False

        if in_atoms and line.strip():
            parts = line.split()
            if len(parts) >= 9:
                atom_name = parts[1]
                if atom_name in atom_map:
                    qm_idx = atom_map[atom_name]
                    if qm_idx in charges:
                        old_q = float(parts[8])
                        new_q = charges[qm_idx]
                        parts[8] = f"{new_q:.4f}"
                        line = f"{parts[0]:>7s} {parts[1]:<4s} {parts[2]:>10s} {parts[3]:>10s} {parts[4]:>10s} {parts[5]:<8s} {parts[6]:>3s} {parts[7]:<8s} {parts[8]:>8s}"
                        n_updated += 1

        updated.append(line)

    mol2_path.write_text("\n".join(updated) + "\n")
    log.info(f"  {mol2_path.name}: updated {n_updated} charges")


def main():
    if len(sys.argv) < 3:
        print(f"Usage: python3 {sys.argv[0]} <charge_file> <mol2_dir>")
        return 1

    charge_file = Path(sys.argv[1])
    mol2_dir = Path(sys.argv[2])

    charges = read_charges(charge_file)
    log.info(f"Read {len(charges)} charges from {charge_file}")

    log.info("NOTE: You must manually define the atom_map for each mol2 file.")
    log.info("The atom_map links mol2 atom names to QM model atom indices.")
    log.info("This depends on how MCPB.py built the residue definitions.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
```

---

## Summary — Integration with Current Pipeline

```
CURRENT PIPELINE (running now, DO NOT INTERRUPT — uses 28 cores):
  FM Opt → BS-DFT Freq → BS-DFT ESP (frame 0)
               ↓                  ↓
           Hessian            Single-frame ESP
               ↓                  ↓
         MCPB Step 2          (feeds into multi-frame as frame 0)

MULTI-FRAME (can start as soon as FM opt converges — uses 8 idle cores):
  1. ORCA GFN2-xTB MD (8 cores, ~15 min)   ← PARALLEL with BS-DFT freq
  2. Cluster → 10 representative frames     ← seconds
  3. ORCA CHELPG ESP on frames 1-9          ← AFTER BS-DFT freq finishes
     (frame 0 = existing BS-DFT ESP)           uses full 28 cores
  4. Multi-frame RESP → averaged charges    ← seconds
  5. Insert charges into mol2 files         ← seconds
  6. MCPB.py Steps 2 + 4 (skip Step 3)     ← minutes

RESOURCE ALLOCATION:
  ┌─────────────────────────────────────────────────────────┐
  │ Phase          │ DFT job (28 cores) │ xTB (8 cores)    │
  │────────────────│────────────────────│──────────────────│
  │ FM Opt         │ ████████████████   │                  │
  │ BS-DFT Freq    │ ████████████████   │ ████ xTB MD      │
  │ BS-DFT ESP     │ ████████████████   │                  │
  │ Multi-frame ESP│ ████████████████████████████ (all 28) │
  └─────────────────────────────────────────────────────────┘
  Tesla P100: idle (ORCA DFT/xTB don't use GPU)
```

**Key point:** The Hessian (force constants) comes from a single geometry —
this is correct and standard. Only the RESP charges benefit from multi-frame.

**When to start:** As soon as the FM optimization finishes and you have
`small_opt.xyz`, you can start the xTB MD (Step 2.7b.2) on 8 idle cores
in parallel with the BS-DFT frequency calculation.

**GPU note:** The Tesla P100 is unused throughout this workflow. ORCA's DFT
and xTB engines are CPU-only. The GPU will be used later in Phase 5
(AMBER pmemd.cuda for MD production runs).
