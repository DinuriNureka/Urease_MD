#!/usr/bin/env python3
"""
scripts/orca_hess_to_fchk.py
---------------------------------------------------------------------------
Convert ORCA .hess file to Gaussian-style .fchk for MCPB.py Step 2.

MCPB.py's Seminario method reader expects:
  - Number of atoms
  - Atomic numbers
  - Cartesian coordinates (Bohr)
  - Cartesian force constants (lower-triangular Hessian in Hartree/Bohr^2)

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

ELEMENT_TO_Z = {
    "H": 1, "He": 2, "C": 6, "N": 7, "O": 8, "F": 9, "S": 16,
    "Cl": 17, "Ni": 28, "NI": 28, "Zn": 30, "Cu": 29, "Fe": 26, "Co": 27,
}


def parse_orca_hess(
    hess_file: Path,
) -> Tuple[int, List[int], np.ndarray, np.ndarray]:
    """
    Parse ORCA .hess file.
    Returns: (n_atoms, atomic_numbers, coords_bohr, hessian_matrix)
    """
    with open(hess_file) as f:
        content = f.read()

    sections: dict = {}
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

    atom_lines = sections.get("$atoms", [])
    n_atoms = int(atom_lines[0].strip())
    atomic_numbers = []
    for i in range(1, n_atoms + 1):
        parts = atom_lines[i].split()
        element = parts[0].strip()
        # ORCA may write "Ni" or "NI"; never fall back to mass (int) for Ni (58.7 -> 58)
        z = ELEMENT_TO_Z.get(element, ELEMENT_TO_Z.get(element.upper()))
        if z is None:
            z = int(float(parts[1]))
        atomic_numbers.append(z)

    # ORCA .hess stores Cartesian coordinates in the $atoms section:
    #   label  mass  x  y  z   (in bohr)
    coords = np.zeros(3 * n_atoms)
    for i in range(1, n_atoms + 1):
        parts = atom_lines[i].split()
        coords[3 * (i - 1) + 0] = float(parts[2])
        coords[3 * (i - 1) + 1] = float(parts[3])
        coords[3 * (i - 1) + 2] = float(parts[4])

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
        try:
            if all(
                p.lstrip("-").replace(".", "").isdigit() and "." not in p
                for p in parts
            ):
                col_offset = int(parts[0])
                i += 1
                continue
        except ValueError:
            pass
        try:
            row_idx = int(parts[0])
            for j, val in enumerate(parts[1:]):
                hessian[row_idx][col_offset + j] = float(val)
        except (ValueError, IndexError):
            pass
        i += 1

    log.info(f"Parsed: {n_atoms} atoms, {n_dof}x{n_dof} Hessian")
    return n_atoms, atomic_numbers, coords, hessian


def write_fchk(
    fchk_file: Path,
    n_atoms: int,
    atomic_numbers: List[int],
    coords_bohr: np.ndarray,
    hessian: np.ndarray,
) -> None:
    """Write Gaussian-format fchk file."""
    n_dof = 3 * n_atoms
    n_hess = n_dof * (n_dof + 1) // 2

    with open(fchk_file, "w") as f:
        f.write("ORCA-to-fchk conversion for MCPB.py Seminario method\n")
        f.write("SP        RB3LYP                          def2-SVP\n")

        f.write(
            f"Number of atoms                            I{n_atoms:>12d}\n"
        )

        f.write(
            f"Atomic numbers                             I   N="
            f"{n_atoms:>12d}\n"
        )
        for i in range(0, n_atoms, 6):
            chunk = atomic_numbers[i : i + 6]
            f.write("".join(f"{z:12d}" for z in chunk) + "\n")

        f.write(
            f"Current cartesian coordinates              R   N="
            f"{3 * n_atoms:>12d}\n"
        )
        for i in range(0, len(coords_bohr), 5):
            chunk = coords_bohr[i : i + 5]
            f.write("".join(f"{x:16.8E}" for x in chunk) + "\n")

        f.write(
            f"Cartesian Force Constants                  R   N="
            f"{n_hess:>12d}\n"
        )
        hess_lt = []
        for i in range(n_dof):
            for j in range(i + 1):
                hess_lt.append(hessian[i][j])

        for i in range(0, len(hess_lt), 5):
            chunk = hess_lt[i : i + 5]
            f.write("".join(f"{x:16.8E}" for x in chunk) + "\n")

    log.info(f"Written {fchk_file}: {n_hess} Hessian elements")


def main() -> int:
    """Entry point. Returns 0 on success, 1 on error."""
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
