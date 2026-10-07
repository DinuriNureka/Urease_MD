#!/usr/bin/env python3
"""Axis-aligned solute extent from a PDB; print tleap setBox centers {cx cy cz} for target box edges (Å).

Usage:
  compute_clearance_for_box.py <solute.pdb> <a> <b> <c>

Clearance per axis: (target_edge - extent) / 2. If any clearance is negative, the solute does not fit.
"""
from __future__ import annotations

import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 5:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    pdb = Path(sys.argv[1])
    ta, tb, tc = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    lines = pdb.read_text(errors="replace").splitlines()
    xs, ys, zs = [], [], []
    for ln in lines:
        if ln.startswith(("ATOM", "HETATM")) and len(ln) >= 38:
            try:
                xs.append(float(ln[30:38]))
                ys.append(float(ln[38:46]))
                zs.append(float(ln[46:54]))
            except ValueError:
                continue
    if not xs:
        print("No ATOM/HETATM coordinates parsed.", file=sys.stderr)
        return 1
    ex = max(xs) - min(xs)
    ey = max(ys) - min(ys)
    ez = max(zs) - min(zs)
    cx = (ta - ex) / 2.0
    cy = (tb - ey) / 2.0
    cz = (tc - ez) / 2.0
    if min(cx, cy, cz) < 0:
        print(
            f"Solute extent (Å): {ex:.4f} x {ey:.4f} x {ez:.4f}; "
            f"does not fit in {ta} x {tb} x {tc}.",
            file=sys.stderr,
        )
        return 1
    print(f"{{{cx:.7f} {cy:.7f} {cz:.7f}}}")
    print(
        f"# extent_xyz= {ex:.6f} {ey:.6f} {ez:.6f}  target= {ta} {tb} {tc}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
