#!/usr/bin/env python3
"""
apply_pc_chelpg_to_mol2.py

Replace partial charges in an antechamber mol2 with values from ORCA resp.pc_chelpg
(point-charge lines: Q  charge  x  y  z), preserving atom order (must match ORCA).
"""

from __future__ import annotations

import argparse
from pathlib import Path


def read_pc_chelpg_charges(path: Path) -> list[float]:
    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError("empty pc_chelpg file")
    n = int(lines[0].strip())
    charges: list[float] = []
    for i in range(2, 2 + n):
        parts = lines[i].split()
        if len(parts) < 2:
            raise ValueError(f"bad pc_chelpg line: {lines[i]}")
        charges.append(float(parts[1]))
    if len(charges) != n:
        raise ValueError("charge count mismatch")
    return charges


def apply_charges(mol2_path: Path, charges: list[float], out_path: Path) -> None:
    lines = mol2_path.read_text().splitlines()
    out: list[str] = []
    i = 0
    in_atom = False
    ai = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom = True
            out.append(line)
            i += 1
            continue
        if in_atom:
            if line.strip() == "" or line.startswith("@"):
                in_atom = False
                out.append(line)
                i += 1
                continue
            parts = line.split()
            if len(parts) >= 9:
                if ai >= len(charges):
                    raise ValueError("more mol2 atoms than charges")
                parts[8] = f"{charges[ai]:.6f}"
                ai += 1
                out.append(
                    "{:>7} {:<4} {:>10.4f} {:>10.4f} {:>10.4f} {:<8} {:>3} {:>3} {:>10.6f}".format(
                        int(parts[0]),
                        parts[1],
                        float(parts[2]),
                        float(parts[3]),
                        float(parts[4]),
                        parts[5],
                        int(parts[6]),
                        parts[7],
                        float(parts[8]),
                    )
                )
            else:
                out.append(line)
            i += 1
            continue
        out.append(line)
        i += 1
    if ai != len(charges):
        raise ValueError(f"atom/charge count mismatch: mol2={ai} pc={len(charges)}")
    out_path.write_text("\n".join(out) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mol2", required=True)
    ap.add_argument("--pc", required=True, help="resp.pc_chelpg from ORCA")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    chg = read_pc_chelpg_charges(Path(args.pc))
    apply_charges(Path(args.mol2), chg, Path(args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
