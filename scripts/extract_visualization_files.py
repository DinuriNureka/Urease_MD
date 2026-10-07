#!/usr/bin/env python3
"""
extract_visualization_files.py

Create convenient PDBs for visualization from an AMBER solvated system:

1) Full system PDB (protein + ligand + Ni center + water + ions) via ambpdb.
2) "No water/ions" PDB: protein + ligand + Ni center only.
3) "Water shell" PDB: protein + ligand + Ni center + waters within R Å of ligand.

All three PDBs are copied to /mnt/user-data/outputs/ for download.

Notes:
  - Uses AmberTools' ambpdb (or bashambpdb if that's what exists on PATH).
  - Uses ParmEd + NumPy to compute the ligand-centered water shell robustly.
  - Default residue names assume this project conventions (LIG ligand, WAT water, Na+/Cl- ions, NI metals).
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np


DEFAULT_OUTPUTS_DIR = Path("/mnt/user-data/outputs")

DEFAULT_WATER_RESNAMES = ("WAT", "HOH", "TIP3", "TIP4", "SOL")
DEFAULT_ION_RESNAMES = ("Na+", "Cl-", "NA", "CL", "K+", "BR-", "I-", "CS+")
DEFAULT_METAL_RESNAMES = ("NI",)


def die(msg: str) -> "None":
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(2)


def which_ambpdb(explicit: str | None) -> str:
    if explicit:
        p = shutil.which(explicit) or explicit
        return p
    # User request mentions bashambpdb; AmberTools normally provides ambpdb.
    for cand in ("bashambpdb", "ambpdb"):
        p = shutil.which(cand)
        if p:
            return p
    die("Could not find ambpdb (or bashambpdb) on PATH. Set --ambpdb /path/to/ambpdb.")


def run_ambpdb_to_file(ambpdb: str, prmtop: Path, crd: Path, out_pdb: Path) -> None:
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    cmd = [ambpdb, "-p", str(prmtop), "-c", str(crd)]
    try:
        with out_pdb.open("w") as f:
            subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        die(f"ambpdb failed: {' '.join(cmd)}\n{(e.stderr or '')[:2000]}")


def pdb_filter_nowater_noions(full_pdb: Path, out_pdb: Path, drop_resnames: set[str]) -> None:
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    with full_pdb.open() as fin, out_pdb.open("w") as fout:
        for line in fin:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                # Keep header-ish records, TER, END, etc.
                fout.write(line)
                continue
            # PDB residue name columns 18-20 (1-indexed), i.e. [17:20] 0-indexed
            res = line[17:20].strip()
            if res in drop_resnames:
                continue
            fout.write(line)


def load_parmed(prmtop: Path, crd: Path):
    try:
        import parmed as pmd
    except ImportError:
        die("ParmEd is not importable (python -c 'import parmed' failed). Install AmberTools/ParmEd.")
    try:
        return pmd.load_file(str(prmtop), str(crd))
    except Exception as e:  # noqa: BLE001
        die(f"ParmEd failed to load prmtop+crd: {e}")


def coords_array(struct) -> np.ndarray:
    return np.array([[a.xx, a.xy, a.xz] for a in struct.atoms], dtype=np.float64)


def select_atom_indices_water_shell(
    struct,
    ligand_resname: str,
    water_resnames: set[str],
    ion_resnames: set[str],
    metal_resnames: set[str],
    radius_a: float,
) -> list[int]:
    """
    Return atom indices to keep for: (protein + ligand + metals) + (waters within radius of ligand).
    """
    crd = coords_array(struct)

    lig_idx = [i for i, a in enumerate(struct.atoms) if a.residue.name.strip() == ligand_resname]
    if not lig_idx:
        die(f"No ligand atoms found for residue name {ligand_resname!r}")

    # Always keep: everything that is NOT water or ion.
    keep = set()
    for i, a in enumerate(struct.atoms):
        r = a.residue.name.strip()
        if r in ion_resnames:
            continue
        if r in water_resnames:
            continue
        # keep proteins, ligand, NI/OHB/KCX, etc.
        keep.add(i)

    # Map water residue -> oxygen atom indices (fallback: all atoms if no oxygen).
    water_res_to_atoms: dict[int, list[int]] = {}
    water_res_to_o: dict[int, list[int]] = {}
    for i, a in enumerate(struct.atoms):
        r = a.residue.name.strip()
        if r not in water_resnames:
            continue
        resid = a.residue.idx  # 0-based residue index in ParmEd
        water_res_to_atoms.setdefault(resid, []).append(i)
        if a.name.strip().upper().startswith("O"):
            water_res_to_o.setdefault(resid, []).append(i)

    lig_xyz = crd[lig_idx]
    r2 = float(radius_a) ** 2

    for resid, atom_idx in water_res_to_atoms.items():
        probe = water_res_to_o.get(resid) or atom_idx
        w_xyz = crd[probe]
        # min distance from any water probe atom to any ligand atom
        d2min = np.min(np.sum((w_xyz[:, None, :] - lig_xyz[None, :, :]) ** 2, axis=2))
        if d2min <= r2:
            keep.update(atom_idx)

    return sorted(keep)


def write_pdb_from_parmed_subset(struct, atom_indices: list[int], out_pdb: Path) -> None:
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    # ParmEd supports slicing by atom indices list.
    sub = struct[atom_indices]
    sub.save(str(out_pdb), overwrite=True)


def copy_outputs(files: list[Path], outputs_dir: Path) -> bool:
    """Copy PDBs to outputs_dir. Returns False if the directory cannot be created (e.g. permission)."""
    try:
        outputs_dir.mkdir(parents=True, exist_ok=True)
    except (PermissionError, OSError) as e:
        print(f"WARNING: could not create outputs directory {outputs_dir}: {e}", file=sys.stderr)
        return False
    for p in files:
        if not p.is_file():
            die(f"Expected output file missing: {p}")
        shutil.copy2(p, outputs_dir / p.name)
    return True


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Create visualization PDBs (full / no-water / ligand watershell) and copy to outputs."
    )
    ap.add_argument("--prmtop", required=True, type=Path, help="AMBER topology (.prmtop)")
    ap.add_argument("--crd", required=True, type=Path, help="AMBER coordinates (.rst7/.inpcrd/.ncrst)")
    ap.add_argument("--out-prefix", default="system", help="Output prefix for PDBs (default: system)")
    ap.add_argument("--ligand-resname", default="LIG", help="Ligand residue name (default: LIG)")
    ap.add_argument("--radius-a", type=float, default=5.0, help="Water shell radius in Å (default: 5.0)")
    ap.add_argument(
        "--water-resnames",
        default=",".join(DEFAULT_WATER_RESNAMES),
        help=f"Comma-separated water residue names (default: {','.join(DEFAULT_WATER_RESNAMES)})",
    )
    ap.add_argument(
        "--ion-resnames",
        default=",".join(DEFAULT_ION_RESNAMES),
        help=f"Comma-separated ion residue names to drop (default: {','.join(DEFAULT_ION_RESNAMES)})",
    )
    ap.add_argument(
        "--metal-resnames",
        default=",".join(DEFAULT_METAL_RESNAMES),
        help=f"Comma-separated metal residue names (informational; default: {','.join(DEFAULT_METAL_RESNAMES)})",
    )
    ap.add_argument(
        "--ambpdb",
        default=None,
        help="Path or name of ambpdb executable (default: auto-detect bashambpdb/ambpdb on PATH)",
    )
    ap.add_argument(
        "--outputs-dir",
        type=Path,
        default=DEFAULT_OUTPUTS_DIR,
        help=f"Directory to copy outputs for download (default: {DEFAULT_OUTPUTS_DIR})",
    )
    ap.add_argument(
        "--workdir",
        type=Path,
        default=Path.cwd(),
        help="Directory to write intermediate PDBs (default: current directory)",
    )

    args = ap.parse_args()
    prmtop: Path = args.prmtop
    crd: Path = args.crd
    workdir: Path = args.workdir

    if not prmtop.is_file():
        die(f"Missing prmtop: {prmtop}")
    if not crd.is_file():
        die(f"Missing crd: {crd}")

    ambpdb = which_ambpdb(args.ambpdb)

    prefix = args.out_prefix
    full_pdb = workdir / f"{prefix}_full.pdb"
    nowater_pdb = workdir / f"{prefix}_nowater.pdb"
    shell_pdb = workdir / f"{prefix}_watershell_{args.radius_a:.1f}A.pdb"

    water_res = {x.strip() for x in str(args.water_resnames).split(",") if x.strip()}
    ion_res = {x.strip() for x in str(args.ion_resnames).split(",") if x.strip()}
    metal_res = {x.strip() for x in str(args.metal_resnames).split(",") if x.strip()}

    # 1) Full system PDB via ambpdb
    run_ambpdb_to_file(ambpdb, prmtop, crd, full_pdb)

    # 2) No-water/ions PDB (filtering the full PDB)
    drop = set(water_res) | set(ion_res) | {"HOH"}  # keep HOH for safety even if user omitted
    pdb_filter_nowater_noions(full_pdb, nowater_pdb, drop_resnames=drop)

    # 3) Watershell PDB via ParmEd selection
    struct = load_parmed(prmtop, crd)
    keep_idx = select_atom_indices_water_shell(
        struct=struct,
        ligand_resname=args.ligand_resname,
        water_resnames=water_res,
        ion_resnames=ion_res,
        metal_resnames=metal_res,
        radius_a=float(args.radius_a),
    )
    write_pdb_from_parmed_subset(struct, keep_idx, shell_pdb)

    # Copy for download
    copied = copy_outputs([full_pdb, nowater_pdb, shell_pdb], outputs_dir=args.outputs_dir)
    fallback = workdir / "visualization_download"
    if not copied and args.outputs_dir == DEFAULT_OUTPUTS_DIR:
        if copy_outputs([full_pdb, nowater_pdb, shell_pdb], outputs_dir=fallback):
            print(f"NOTE: copied to fallback directory (Cursor download path unavailable): {fallback.resolve()}", flush=True)

    print("Wrote:")
    print(f"  {full_pdb}")
    print(f"  {nowater_pdb}")
    print(f"  {shell_pdb}")
    if copied:
        print(f"Copied to: {args.outputs_dir.resolve()}")
    elif (fallback / full_pdb.name).is_file():
        pass  # already printed fallback note
    else:
        print("Copy to download dir skipped; PDBs are only under --workdir.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

