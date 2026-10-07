#!/usr/bin/env python3
"""Print 1-based index of last residue before solvent (WAT/HOH/ions) for Amber prmtop."""
from __future__ import annotations

import sys


def last_solute_residue(prmtop_path: str) -> int:
    with open(prmtop_path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    i = 0
    resnum = 0
    while i < len(lines):
        line = lines[i]
        if line.strip() == "%FLAG RESIDUE_LABEL":
            i += 2  # %FORMAT
            while i < len(lines) and not lines[i].startswith("%"):
                row = lines[i].rstrip("\n")
                for k in range(0, len(row), 4):
                    chunk = row[k : k + 4].strip()
                    if not chunk:
                        continue
                    resnum += 1
                    if chunk in ("WAT", "HOH", "Na+", "Cl-", "K+", "Br-", "Cs+", "I-"):
                        return resnum - 1
                i += 1
            break
        i += 1

    return resnum


def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: _prmtop_last_solute_residue.py <prmtop>", file=sys.stderr)
        sys.exit(2)
    n = last_solute_residue(sys.argv[1])
    if n < 1:
        print("[FAIL] no residues parsed", file=sys.stderr)
        sys.exit(1)
    print(n)


if __name__ == "__main__":
    main()
