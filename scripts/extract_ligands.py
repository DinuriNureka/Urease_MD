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

import sys
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
        log.warning(
            f"  {xyz_file.name}: {total} atoms, expected > {RESIDUE_ATOMS}"
        )
        return 0

    out = output_dir / f"{xyz_file.stem}_lig.xyz"
    with open(out, "w") as f:
        f.write(f"{n_lig}\n")
        f.write(f"Ligand from {xyz_file.name}\n")
        for i in range(2, 2 + n_lig):
            f.write(lines[i] + "\n")

    return n_lig


def main() -> int:
    """Entry point. Returns 0 on success, 1 on error."""
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
