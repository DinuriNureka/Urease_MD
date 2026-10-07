#!/usr/bin/env python3
"""
Create systems/<compound>_pose<pose>_full/<compound>_pose<pose>_gaff2.mol2 by taking
GAFF2 topology (atom names, types, bonds) from gaff2_params/<mol_id>/<mol_id>.mol2 and
replacing coordinates from a ligand .xyz (e.g. extracted from dft_batch/dft_jobs/*_pose*/cluster.xyz).

When the .xyz atom order matches the mol2 line order (same element sequence), coordinates
are copied directly. Otherwise a Kabsch-aligned Hungarian assignment (same element only) is used.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from build_pipeline_ligand_mol2 import docked_path, elem_from_name, parse_mol2  # noqa: E402


def read_ligand_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    lines = path.read_text().splitlines()
    n = int(lines[0].strip())
    elems: list[str] = []
    coords: list[list[float]] = []
    for line in lines[2 : 2 + n]:
        p = line.split()
        if len(p) < 4:
            continue
        elems.append(elem_from_name(p[0]))
        coords.append([float(p[1]), float(p[2]), float(p[3])])
    if len(elems) != n:
        raise ValueError(f"{path}: expected {n} atoms, parsed {len(elems)}")
    return elems, np.array(coords, dtype=float)


def kabsch_align_p_to_q(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return R, t such that R @ P + t best matches Q (row-wise atoms)."""
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)
    P0 = P - Pc
    Q0 = Q - Qc
    H = P0.T @ Q0
    U, _S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = Qc - R @ Pc
    return R, t, (R @ P.T).T + t


def hungarian_same_element(
    elems_gaff: list[str],
    coords_gaff: np.ndarray,
    elems_xyz: list[str],
    coords_xyz: np.ndarray,
) -> list[int]:
    """Return mapping dock_idx -> gaff_idx (permutation) minimizing squared distance after Kabsch."""
    n = len(elems_gaff)
    if len(elems_xyz) != n:
        raise ValueError(f"Atom count mismatch: gaff {n} vs xyz {len(elems_xyz)}")
    if sorted(elems_gaff) != sorted(elems_xyz):
        raise ValueError("Element multisets differ between mol2 and xyz")
    # direct order match
    if elems_gaff == elems_xyz:
        return list(range(n))
    _, _, G2 = kabsch_align_p_to_q(coords_gaff, coords_xyz)
    C = np.full((n, n), 1e18, dtype=float)
    for i in range(n):
        for j in range(n):
            if elems_gaff[i] == elems_xyz[j]:
                C[i, j] = np.sum((G2[i] - coords_xyz[j]) ** 2)
    row, col = linear_sum_assignment(C)
    if not np.array_equal(np.sort(row), np.arange(n)):
        raise RuntimeError("linear_sum_assignment row mismatch")
    # row[i] is gaff index, col[i] is xyz index — we want for each gaff i, which xyz j
    gaff_to_xyz = [0] * n
    for r, c in zip(row, col):
        gaff_to_xyz[r] = c
    # dock mol2 order follows gaff mol2 order; coords from xyz[gaff_to_xyz[i]] for gaff atom i
    return gaff_to_xyz


def write_mol2_new_xyz(
    template_lines: list[str],
    n_atoms: int,
    new_xyz: list[tuple[float, float, float]],
    out_path: Path,
) -> None:
    if len(new_xyz) != n_atoms:
        raise ValueError(f"new_xyz len {len(new_xyz)} vs atoms {n_atoms}")
    out_lines: list[str] = []
    atom_i = 0
    i = 0
    while i < len(template_lines):
        line = template_lines[i]
        if line.startswith("@<TRIPOS>ATOM"):
            out_lines.append(line)
            i += 1
            while i < len(template_lines) and not template_lines[i].startswith("@"):
                raw = template_lines[i]
                parts = raw.split()
                if len(parts) >= 9 and atom_i < n_atoms:
                    x, y, z = new_xyz[atom_i]
                    parts[2] = f"{x:.6f}"
                    parts[3] = f"{y:.6f}"
                    parts[4] = f"{z:.6f}"
                    out_lines.append(" ".join(parts))
                    atom_i += 1
                else:
                    out_lines.append(raw)
                i += 1
            continue
        out_lines.append(line)
        i += 1
    if atom_i != n_atoms:
        raise ValueError(f"Atom lines written {atom_i} vs expected {n_atoms}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(out_lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", type=Path, default=Path("/path/to/urease_md"))
    ap.add_argument("--mol-id", required=True, help="e.g. mol01")
    ap.add_argument("--compound", required=True, help="e.g. aha")
    ap.add_argument("--pose", required=True, help="pose number, e.g. 1")
    ap.add_argument("--lig-xyz", type=Path, required=True, help="ligand-only .xyz coordinates")
    args = ap.parse_args()
    root = args.project_root.resolve()
    mol_id = args.mol_id.strip()
    compound = args.compound.strip()
    pose = str(int(args.pose))

    gaff_mol2 = root / "gaff2_params" / mol_id / f"{mol_id}.mol2"
    if not gaff_mol2.is_file():
        print(f"ERROR missing {gaff_mol2}", file=sys.stderr)
        return 1
    ga, _gb, glines = parse_mol2(gaff_mol2)
    elems_xyz, coords_xyz = read_ligand_xyz(args.lig_xyz.resolve())

    coords_gaff = np.array([[a["x"], a["y"], a["z"]] for a in ga], dtype=float)
    elems_gaff = [a["elem"] for a in ga]

    if len(elems_gaff) != len(elems_xyz):
        print(
            f"ERROR atom count mol2 {len(elems_gaff)} vs xyz {len(elems_xyz)}",
            file=sys.stderr,
        )
        return 1

    gaff_to_xyz = hungarian_same_element(elems_gaff, coords_gaff, elems_xyz, coords_xyz)

    new_xyz: list[tuple[float, float, float]] = []
    for i in range(len(ga)):
        j = gaff_to_xyz[i]
        x, y, z = coords_xyz[j]
        new_xyz.append((float(x), float(y), float(z)))

    out = docked_path(root, compound, pose)
    write_mol2_new_xyz(glines, len(ga), new_xyz, out)
    print(f"OK wrote {out} ({len(ga)} atoms) from {gaff_mol2.name} + {args.lig_xyz}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
