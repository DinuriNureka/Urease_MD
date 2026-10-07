#!/usr/bin/env python3
"""
Build params/inhibitors/mol2/<compound>_gaff2.mol2 for each ligand in ligand_name_mapping.csv.

For each row: merge Multiwfn RESP charges from gaff2_params/<mol_id>/<mol_id>.mol2 onto
docked coordinates from systems/<compound>_pose<pose>_full/<compound>_pose<pose>_gaff2.mol2
when the docked file exists (graph isomorphism on element + connectivity).

If graphs differ but atom count and element multiset match, a Kabsch + Hungarian assignment
(minimize distance between same-element pairs) can transfer charges onto the docked pose
(see GEOMETRY_CHARGE_MAP_* warnings — RESP was computed on the gas-phase graph; use with care).

If no docked mol2 exists, copy gaff2_params source as a fallback (non-docked geometry) and record
a warning — see params/inhibitors/mol2/PIPELINE_LIGAND_BUILD_WARNINGS.txt
"""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
from collections import Counter
from pathlib import Path

import networkx as nx
import numpy as np
from networkx.algorithms import isomorphism
from scipy.optimize import linear_sum_assignment

B_TO_A = 0.529177210903  # unused; reserved


def elem_from_name(name: str) -> str:
    name = name.strip()
    m = re.match(r"^([A-Z][a-z]?)", name)
    if not m:
        return name[0]
    return m.group(1)


def parse_mol2(path: Path) -> tuple[list[dict], list[tuple[int, int]], list[str]]:
    """Return atoms (dicts), bonds (0-based), and raw lines for round-trip header preservation."""
    lines = path.read_text().splitlines()
    atoms: list[dict] = []
    bonds: list[tuple[int, int]] = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("@<TRIPOS>ATOM"):
            i += 1
            while i < len(lines) and not lines[i].startswith("@"):
                raw = lines[i]
                parts = raw.split()
                if len(parts) >= 9:
                    atoms.append(
                        {
                            "raw": raw,
                            "name": parts[1],
                            "x": float(parts[2]),
                            "y": float(parts[3]),
                            "z": float(parts[4]),
                            "charge": float(parts[8]),
                            "elem": elem_from_name(parts[1]),
                        }
                    )
                i += 1
            continue
        if lines[i].startswith("@<TRIPOS>BOND"):
            i += 1
            while i < len(lines) and not lines[i].startswith("@"):
                parts = lines[i].split()
                if len(parts) >= 4:
                    a, b = int(parts[1]) - 1, int(parts[2]) - 1
                    bonds.append((a, b))
                i += 1
            continue
        i += 1
    return atoms, bonds, lines


def build_graph(atoms: list[dict], bonds: list[tuple[int, int]]) -> nx.Graph:
    g = nx.Graph()
    for i, a in enumerate(atoms):
        g.add_node(i, elem=a["elem"])
    for a, b in bonds:
        g.add_edge(a, b)
    return g


def isomorphism_map(gaff_atoms: list[dict], gaff_bonds, dock_atoms: list[dict], dock_bonds) -> dict[int, int] | None:
    g1 = build_graph(gaff_atoms, gaff_bonds)
    g2 = build_graph(dock_atoms, dock_bonds)

    def nm(d1: dict, d2: dict) -> bool:
        return d1["elem"] == d2["elem"]

    gm = isomorphism.GraphMatcher(g1, g2, node_match=nm)
    if not gm.is_isomorphic():
        return None
    return next(gm.isomorphisms_iter())


def kabsch_align_p_to_q(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Return coordinates of P rigidly aligned onto Q (same shape, row = atom)."""
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
    return (R @ P.T).T + t


def geometry_charge_map(
    gaff_atoms: list[dict],
    dock_atoms: list[dict],
    mean_max_angstrom: tuple[float, float] = (2.6, 6.5),
) -> tuple[list[float] | None, float, float]:
    """
    Map RESP charges from gaff_atoms onto dock_atoms by same-element Hungarian assignment
    after aligning gaff coordinates to the docked frame. Returns (charges, mean_d, max_d) or (None, ...).
    """
    if len(gaff_atoms) != len(dock_atoms):
        return None, float("nan"), float("nan")
    elems_g = [a["elem"] for a in gaff_atoms]
    elems_d = [a["elem"] for a in dock_atoms]
    if Counter(elems_g) != Counter(elems_d):
        return None, float("nan"), float("nan")
    G = np.array([[a["x"], a["y"], a["z"]] for a in gaff_atoms], dtype=float)
    D = np.array([[a["x"], a["y"], a["z"]] for a in dock_atoms], dtype=float)
    n = len(gaff_atoms)
    G2 = kabsch_align_p_to_q(G, D)
    C = np.full((n, n), 1e18, dtype=float)
    for i in range(n):
        for j in range(n):
            if elems_g[i] == elems_d[j]:
                C[i, j] = float(np.sum((G2[i] - D[j]) ** 2))
    row, col = linear_sum_assignment(C)
    dists = np.sqrt(C[row, col])
    mean_d = float(dists.mean())
    max_d = float(dists.max())
    mean_lim, max_lim = mean_max_angstrom
    if mean_d > mean_lim or max_d > max_lim:
        return None, mean_d, max_d
    gaff_for_dock = [0] * n
    for k in range(n):
        gaff_for_dock[col[k]] = row[k]
    charges = [float(gaff_atoms[gaff_for_dock[j]]["charge"]) for j in range(n)]
    return charges, mean_d, max_d


def docked_path(project_root: Path, compound: str, pose: str) -> Path:
    return (
        project_root
        / "systems"
        / f"{compound}_pose{pose}_full"
        / f"{compound}_pose{pose}_gaff2.mol2"
    )


def write_merged_mol2(
    docked_lines: list[str],
    docked_atoms: list[dict],
    charges_by_dock_idx: list[float],
    out_path: Path,
) -> None:
    """Replace ATOM block charges; preserve line structure except last field."""
    out_lines: list[str] = []
    atom_i = 0
    i = 0
    while i < len(docked_lines):
        line = docked_lines[i]
        if line.startswith("@<TRIPOS>ATOM"):
            out_lines.append(line)
            i += 1
            while i < len(docked_lines) and not docked_lines[i].startswith("@"):
                raw = docked_lines[i]
                parts = raw.split()
                if len(parts) >= 9 and atom_i < len(charges_by_dock_idx):
                    parts[8] = f"{charges_by_dock_idx[atom_i]:.6f}"
                    out_lines.append(" ".join(parts))
                    atom_i += 1
                else:
                    out_lines.append(raw)
                i += 1
            continue
        out_lines.append(line)
        i += 1
    if atom_i != len(docked_atoms):
        raise ValueError(f"Atom count mismatch writing mol2: {atom_i} vs {len(docked_atoms)}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(out_lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", type=Path, default=Path("/path/to/urease_md"))
    ap.add_argument(
        "--mapping",
        type=Path,
        default=None,
        help="ligand_name_mapping.csv",
    )
    args = ap.parse_args()
    root = args.project_root.resolve()
    mapping = args.mapping or (root / "params/inhibitors/orca_opt/ligand_name_mapping.csv")
    gaff_root = root / "gaff2_params"
    out_dir = root / "params/inhibitors/mol2"
    warn_path = out_dir / "PIPELINE_LIGAND_BUILD_WARNINGS.txt"
    warnings: list[str] = []

    rows: list[dict] = []
    with mapping.open(newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append(row)

    override_csv = root / "params/inhibitors/docked_mol2_overrides.csv"
    dock_overrides: dict[str, Path] = {}
    if override_csv.is_file() and override_csv.stat().st_size > 20:
        with override_csv.open(newline="") as f:
            for row in csv.DictReader(f):
                c = row.get("compound", "").strip()
                p = row.get("docked_mol2_path", "").strip()
                if c and p:
                    dock_overrides[c] = Path(p).expanduser()

    for row in rows:
        mol_id = row["mol_id"].strip()
        compound = row["compound"].strip()
        pose = str(int(row["selected_pose"]))  # normalize 1 not 01
        net = int(row["net_charge"])

        gaff_mol2 = gaff_root / mol_id / f"{mol_id}.mol2"
        target = out_dir / f"{compound}_gaff2.mol2"

        if not gaff_mol2.is_file():
            print(f"ERROR missing gaff2 source {gaff_mol2}", file=sys.stderr)
            return 1

        dock = dock_overrides.get(compound) or docked_path(root, compound, pose)
        if dock.is_file():
            ga, gb, _ = parse_mol2(gaff_mol2)
            da, db, dlines = parse_mol2(dock)
            if len(ga) != len(da):
                print(f"ERROR atom count {mol_id} gaff {len(ga)} docked {len(da)}", file=sys.stderr)
                return 1
            m = isomorphism_map(ga, gb, da, db)
            if m is None:
                geom_charges, mean_d, max_d = geometry_charge_map(ga, da)
                if geom_charges is not None:
                    s = sum(geom_charges)
                    if abs(s - net) > 0.05:
                        warnings.append(
                            f"{compound} ({mol_id}): charge sum {s:.4f} vs CSV net_charge {net} (geometry map)"
                        )
                    msg = (
                        f"{compound} ({mol_id}): graph mismatch — used geometry-mapped RESP onto {dock.name} "
                        f"(mean d={mean_d:.3f} Å, max={max_d:.3f} Å; RESP from gas-phase {gaff_mol2.name})."
                    )
                    warnings.append(msg)
                    print(f"WARN {msg}")
                    write_merged_mol2(dlines, da, geom_charges, target)
                    print(f"OK merge (geometry map) {mol_id} -> {target.name} (docked {dock.name})")
                else:
                    shutil.copy2(gaff_mol2, target)
                    msg = (
                        f"{compound} ({mol_id}): graph mismatch gaff vs {dock.name} — copied {gaff_mol2.name} "
                        f"(gaff2_params geometry; docked topology differs from RESP source"
                        + (
                            f"; geometry map failed (mean={mean_d:.3f}, max={max_d:.3f} Å)."
                            if not (np.isnan(mean_d) or np.isnan(max_d))
                            else "."
                        )
                    )
                    warnings.append(msg)
                    print(f"WARN {msg}")
            else:
                inv: dict[int, int] = {}
                for g_idx, d_idx in m.items():
                    inv[d_idx] = g_idx
                charges = [ga[inv[i]]["charge"] for i in range(len(da))]
                s = sum(charges)
                if abs(s - net) > 0.05:
                    warnings.append(
                        f"{compound} ({mol_id}): charge sum {s:.4f} vs CSV net_charge {net} (after merge)"
                    )
                write_merged_mol2(dlines, da, charges, target)
                print(f"OK merge {mol_id} -> {target.name} (docked {dock.name})")
        else:
            shutil.copy2(gaff_mol2, target)
            msg = (
                f"{compound} ({mol_id}): no docked file {dock} — copied {gaff_mol2} "
                f"(ORCA/gaff2_params geometry, not pocket-docked)."
            )
            warnings.append(msg)
            print(f"WARN {msg}")

    if warnings:
        warn_path.write_text("\n".join(warnings) + "\n")
        print(f"Wrote {warn_path} ({len(warnings)} warning(s))")
    elif warn_path.is_file():
        warn_path.unlink()

    # GAFF2 frcmod for tleap: pipeline expects params/inhibitors/frcmod/<compound>.frcmod
    frc_dir = root / "params/inhibitors/frcmod"
    frc_dir.mkdir(parents=True, exist_ok=True)
    for row in rows:
        mol_id = row["mol_id"].strip()
        compound = row["compound"].strip()
        src = gaff_root / mol_id / f"{mol_id}.frcmod"
        dst = frc_dir / f"{compound}.frcmod"
        if not src.is_file():
            print(f"WARN missing frcmod source {src}", file=sys.stderr)
            continue
        shutil.copy2(src, dst)
        print(f"OK frcmod {mol_id} -> {dst.name}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
