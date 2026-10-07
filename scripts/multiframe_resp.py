#!/usr/bin/env python3
"""
Compute multi-frame averaged CHELPG charges from ORCA outputs.

This follows the current project's pragmatic route:
1) run CHELPG on multiple frames
2) average per-atom charges across frames
3) report per-atom standard deviation
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple

import numpy as np


def read_chelpg_table(out_text: str) -> List[Tuple[int, str, float]]:
    rows: List[Tuple[int, str, float]] = []
    in_block = False
    for line in out_text.splitlines():
        if "CHELPG Charges" in line:
            in_block = True
            continue
        if not in_block:
            continue
        s = line.strip()
        if not s:
            if rows:
                break
            continue
        parts = s.split()
        if len(parts) >= 4 and parts[0].isdigit():
            rows.append((int(parts[0]), parts[1], float(parts[3])))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("esp_dir", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--n-atoms", type=int, default=71)
    args = ap.parse_args()

    frame_dirs = sorted(args.esp_dir.glob("frame_*"))
    if not frame_dirs:
        raise RuntimeError(f"No frame_* directories in {args.esp_dir}")

    per_frame = []
    elements_ref: List[str] = []
    atom_ids_ref: List[int] = []

    for fd in frame_dirs:
        out = fd / "esp.out"
        if not out.exists():
            continue
        rows = read_chelpg_table(out.read_text())
        if len(rows) != args.n_atoms:
            continue
        atom_ids = [r[0] for r in rows]
        elements = [r[1] for r in rows]
        charges = [r[2] for r in rows]
        if not per_frame:
            atom_ids_ref = atom_ids
            elements_ref = elements
        else:
            if atom_ids != atom_ids_ref or elements != elements_ref:
                raise RuntimeError(f"Atom ordering mismatch in {out}")
        per_frame.append((fd.name, np.array(charges, dtype=float)))

    if not per_frame:
        raise RuntimeError("No valid esp.out frames with full CHELPG table found")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    names = [n for n, _ in per_frame]
    mat = np.vstack([c for _, c in per_frame])
    avg = mat.mean(axis=0)
    std = mat.std(axis=0)

    avg_path = args.out_dir / "multiframe_avg_charges.dat"
    with avg_path.open("w") as fh:
        fh.write(f"# averaged CHELPG charges across {len(names)} frames\n")
        fh.write("# atom  element     avg_charge      std\n")
        for atom_id, el, q, s in zip(atom_ids_ref, elements_ref, avg, std):
            fh.write(f"{atom_id:4d} {el:>3s} {q:14.8f} {s:12.8f}\n")
        fh.write(f"# total_charge {avg.sum():.8f}\n")

    cmp_path = args.out_dir / "frame_charge_comparison.dat"
    with cmp_path.open("w") as fh:
        header = "atom element " + " ".join(f"{n:>12s}" for n in names) + "    avg        std"
        fh.write(header + "\n")
        for i, (atom_id, el) in enumerate(zip(atom_ids_ref, elements_ref)):
            vals = " ".join(f"{mat[j, i]:12.6f}" for j in range(mat.shape[0]))
            fh.write(f"{atom_id:4d} {el:>3s} {vals} {avg[i]:10.6f} {std[i]:10.6f}\n")

    print(f"Used frames: {', '.join(names)}")
    print(f"Wrote {avg_path}")
    print(f"Wrote {cmp_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

