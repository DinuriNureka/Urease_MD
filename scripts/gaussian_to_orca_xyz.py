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
                try:
                    float(parts[1])
                    float(parts[2])
                    float(parts[3])
                    coords.append(stripped)
                except ValueError:
                    continue
            elif stripped == "" and len(coords) > 0:
                break

    return charge, mult, coords


def write_xyz(coords: List[str], output_file: Path, comment: str = "") -> None:
    """Write standard XYZ file."""
    with open(output_file, "w") as f:
        f.write(f"{len(coords)}\n")
        f.write(f"{comment}\n")
        for line in coords:
            f.write(line + "\n")


def main() -> int:
    """Entry point. Returns 0 on success, 1 on error."""
    if len(sys.argv) != 3:
        print(f"Usage: python3 {sys.argv[0]} <input.com> <output.xyz>")
        return 1

    com_file = Path(sys.argv[1])
    xyz_file = Path(sys.argv[2])

    charge, mult, coords = parse_gaussian_com(com_file)
    log.info(
        f"Parsed {com_file.name}: charge={charge}, mult={mult}, "
        f"atoms={len(coords)}"
    )

    write_xyz(
        coords, xyz_file,
        f"charge={charge} mult={mult} from {com_file.name}"
    )
    log.info(f"Written {xyz_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
