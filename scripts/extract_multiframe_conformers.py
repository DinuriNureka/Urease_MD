#!/usr/bin/env python3
"""
Extract representative conformations from a multi-XYZ trajectory.

This script reads ORCA xTB MD trajectory.xyz, downsamples in time, then
selects diverse frames via farthest-point sampling on heavy-atom RMSD.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Sequence, Tuple

import numpy as np


Frame = Tuple[List[str], np.ndarray, str]


def read_multi_xyz(path: Path, n_atoms: int) -> List[Frame]:
    lines = path.read_text().splitlines()
    frames: List[Frame] = []
    i = 0
    while i < len(lines):
        nat = int(lines[i].strip())
        if nat != n_atoms:
            raise ValueError(f"Expected {n_atoms} atoms, found {nat} at line {i + 1}")
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        elements: List[str] = []
        coords = np.zeros((nat, 3), dtype=float)
        for j in range(nat):
            parts = lines[i + 2 + j].split()
            elements.append(parts[0])
            coords[j] = [float(parts[1]), float(parts[2]), float(parts[3])]
        frames.append((elements, coords, comment))
        i += nat + 2
    return frames


def kabsch_align_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    # Align b onto a and compute RMSD.
    a_cent = a - a.mean(axis=0)
    b_cent = b - b.mean(axis=0)
    h = b_cent.T @ a_cent
    u, _, vt = np.linalg.svd(h)
    r = u @ vt
    if np.linalg.det(r) < 0:
        u[:, -1] *= -1
        r = u @ vt
    b_rot = b_cent @ r
    diff = a_cent - b_rot
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def representative_indices(
    frames: Sequence[Frame], n_pick: int, stride: int, heavy_only: bool = True
) -> List[int]:
    candidate = list(range(0, len(frames), max(1, stride)))
    if not candidate:
        return []
    if len(candidate) <= n_pick:
        return candidate

    # Build masks from first frame.
    elems0 = frames[candidate[0]][0]
    heavy_mask = np.array([el.upper() != "H" for el in elems0], dtype=bool)

    def coords(idx: int) -> np.ndarray:
        c = frames[idx][1]
        return c[heavy_mask] if heavy_only else c

    selected = [candidate[0]]
    min_dist = {idx: float("inf") for idx in candidate}

    while len(selected) < n_pick:
        ref = coords(selected[-1])
        for idx in candidate:
            if idx in selected:
                continue
            d = kabsch_align_rmsd(ref, coords(idx))
            if d < min_dist[idx]:
                min_dist[idx] = d
        next_idx = max((idx for idx in candidate if idx not in selected), key=min_dist.get)
        selected.append(next_idx)

    return sorted(selected)


def write_xyz(path: Path, elements: Sequence[str], xyz: np.ndarray, comment: str) -> None:
    with path.open("w") as fh:
        fh.write(f"{len(elements)}\n")
        fh.write(f"{comment}\n")
        for el, (x, y, z) in zip(elements, xyz):
            fh.write(f"{el:>2s} {x:16.10f} {y:16.10f} {z:16.10f}\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("trajectory", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--n-atoms", type=int, default=71)
    ap.add_argument("--n-frames", type=int, default=10)
    ap.add_argument("--stride", type=int, default=10)
    args = ap.parse_args()

    frames = read_multi_xyz(args.trajectory, args.n_atoms)
    if not frames:
        raise RuntimeError("No frames found in trajectory")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    picks = representative_indices(frames, args.n_frames, args.stride, heavy_only=True)
    for out_idx, src_idx in enumerate(picks):
        elems, xyz, comment = frames[src_idx]
        out = args.out_dir / f"frame_{out_idx:02d}.xyz"
        write_xyz(out, elems, xyz, f"source_frame={src_idx}; {comment}".strip())
    print(f"Wrote {len(picks)} representative frames to {args.out_dir}")
    print("Selected source indices:", ",".join(str(i) for i in picks))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

