#!/usr/bin/env python3
"""
select_pose1.py

Purpose:
  Extract pose 1 (first record) from each GNINA multi-model SDF, generate clean XYZ
  files for ORCA, and write a charge/molecule manifest for batch RESP jobs.

Inputs:
  - Ligand SDF files under --ligand_dir (default: ligands/)

Outputs:
  - resp_inputs/molNN_pose1.xyz
  - resp_inputs/job_manifest.csv
  - resp_inputs/failed_parse.log (for molecules that fail parsing)

Example:
  python3 select_pose1.py --ligand_dir ligands/
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass
class JobRow:
    mol_id: str
    sdf_file: str
    xyz_file: str
    net_charge: int
    multiplicity: int
    n_atoms: int
    n_heavy_atoms: int


def detect_tools() -> tuple[bool, bool]:
    has_rdkit = False
    has_obabel = False
    try:
        import rdkit  # noqa: F401

        has_rdkit = True
    except Exception:
        has_rdkit = False

    try:
        subprocess.run(
            ["obabel", "-V"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        has_obabel = True
    except FileNotFoundError:
        has_obabel = False
    return has_rdkit, has_obabel


def first_sdf_record(sdf_text: str) -> str:
    parts = sdf_text.split("$$$$")
    if not parts:
        raise ValueError("SDF has no molecule records")
    first = parts[0].strip("\n\r\t ")
    if not first:
        raise ValueError("First SDF record is empty")
    return first + "\n$$$$\n"


def parse_charge_from_molblock_fallback(molblock: str) -> int:
    """
    Fallback parser for formal charges from MDL charge lines if RDKit is unavailable.
    Handles:
      - M  CHG lines
      - V2000 atom-block charge code column
    """
    lines = molblock.splitlines()
    if len(lines) < 4:
        return 0

    total_charge = 0
    # Prefer explicit M  CHG entries if present
    chg_lines = [ln for ln in lines if ln.startswith("M  CHG")]
    if chg_lines:
        for ln in chg_lines:
            toks = ln.split()
            if len(toks) < 4:
                continue
            try:
                n_pairs = int(toks[2])
            except ValueError:
                continue
            expected = 3 + 2 * n_pairs
            vals = toks[3:expected]
            for i in range(1, len(vals), 2):
                try:
                    total_charge += int(vals[i])
                except ValueError:
                    continue
        return total_charge

    # V2000 atom block charge code mapping
    try:
        nat = int(lines[3][0:3].strip())
    except Exception:
        return 0
    atom_start = 4
    atom_end = atom_start + nat
    code_to_charge = {
        0: 0,
        1: +3,
        2: +2,
        3: +1,
        4: 0,
        5: -1,
        6: -2,
        7: -3,
    }
    for ln in lines[atom_start:atom_end]:
        toks = ln.split()
        if len(toks) < 6:
            continue
        try:
            code = int(toks[5])
        except ValueError:
            continue
        total_charge += code_to_charge.get(code, 0)
    return total_charge


def rdkit_process_first_record(record: str) -> tuple[list[tuple[str, float, float, float]], int, int]:
    from rdkit import Chem

    mol = Chem.MolFromMolBlock(record, sanitize=True, removeHs=False, strictParsing=False)
    if mol is None:
        raise ValueError("RDKit failed to parse first SDF record")

    # Ensure explicit Hs are present while preserving 3D coordinates.
    mol = Chem.AddHs(mol, addCoords=True)

    conf = mol.GetConformer()
    atoms_xyz: list[tuple[str, float, float, float]] = []
    net_charge = 0
    n_heavy = 0
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        sym = a.GetSymbol()
        atoms_xyz.append((sym, float(p.x), float(p.y), float(p.z)))
        net_charge += int(a.GetFormalCharge())
        if sym.upper() != "H":
            n_heavy += 1
    return atoms_xyz, net_charge, n_heavy


def obabel_xyz_from_record(record: str) -> list[tuple[str, float, float, float]]:
    with tempfile.TemporaryDirectory() as td:
        tdir = Path(td)
        in_sdf = tdir / "pose1.sdf"
        out_xyz = tdir / "pose1.xyz"
        in_sdf.write_text(record)
        subprocess.run(
            ["obabel", "-isdf", str(in_sdf), "-oxyz", "-O", str(out_xyz), "-h"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        lines = out_xyz.read_text().splitlines()
        if len(lines) < 3:
            raise ValueError("OpenBabel generated malformed XYZ")
        atoms: list[tuple[str, float, float, float]] = []
        for ln in lines[2:]:
            toks = ln.split()
            if len(toks) < 4:
                continue
            atoms.append((toks[0], float(toks[1]), float(toks[2]), float(toks[3])))
        if not atoms:
            raise ValueError("OpenBabel generated empty XYZ atom list")
        return atoms


def write_xyz(path: Path, atoms: Iterable[tuple[str, float, float, float]], comment: str) -> tuple[int, int]:
    atom_list = list(atoms)
    n_heavy = sum(1 for sym, *_ in atom_list if sym.upper() != "H")
    with path.open("w") as fh:
        fh.write(f"{len(atom_list)}\n")
        fh.write(comment + "\n")
        for sym, x, y, z in atom_list:
            fh.write(f"{sym:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")
    return len(atom_list), n_heavy


def main() -> int:
    ap = argparse.ArgumentParser(description="Extract pose1 from SDFs and write RESP manifest.")
    ap.add_argument("--ligand_dir", default="ligands", help="Directory containing GNINA SDF files.")
    ap.add_argument("--out_dir", default="resp_inputs", help="Output directory for XYZ and manifest.")
    args = ap.parse_args()

    lig_dir = Path(args.ligand_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    failed_log = out_dir / "failed_parse.log"
    manifest = out_dir / "job_manifest.csv"

    sdf_files = sorted(lig_dir.glob("*.sdf"))
    if not sdf_files:
        raise SystemExit(f"No SDF files found in {lig_dir}")

    has_rdkit, has_obabel = detect_tools()
    if not has_rdkit and not has_obabel:
        raise RuntimeError(
            "Neither RDKit nor OpenBabel is available. Install RDKit (preferred) "
            "or OpenBabel ('obabel') to continue."
        )

    rows: list[JobRow] = []
    failed_entries: list[str] = []
    print(
        f"{'mol_id':<8} {'sdf_file':<24} {'net_q':>6} {'mult':>5} {'n_atoms':>8} {'n_heavy':>8} {'status':>10}"
    )
    print("-" * 78)

    for i, sdf_path in enumerate(sdf_files, start=1):
        mol_id = f"mol{i:02d}"
        xyz_name = f"{mol_id}_pose1.xyz"
        xyz_path = out_dir / xyz_name
        try:
            record = first_sdf_record(sdf_path.read_text())
            if has_rdkit:
                atoms_xyz, net_charge, n_heavy = rdkit_process_first_record(record)
            else:
                atoms_xyz = obabel_xyz_from_record(record)
                net_charge = parse_charge_from_molblock_fallback(record)
                n_heavy = sum(1 for sym, *_ in atoms_xyz if sym.upper() != "H")

            n_atoms, n_heavy2 = write_xyz(
                xyz_path,
                atoms_xyz,
                comment=f"Pose1 from {sdf_path.name}",
            )
            if n_heavy2 != n_heavy:
                n_heavy = n_heavy2

            row = JobRow(
                mol_id=mol_id,
                sdf_file=str(sdf_path),
                xyz_file=str(xyz_path),
                net_charge=int(net_charge),
                multiplicity=1,
                n_atoms=n_atoms,
                n_heavy_atoms=n_heavy,
            )
            rows.append(row)
            print(
                f"{mol_id:<8} {sdf_path.name:<24} {row.net_charge:>6d} {row.multiplicity:>5d} {row.n_atoms:>8d} {row.n_heavy_atoms:>8d} {'OK':>10}"
            )
        except Exception as e:
            msg = f"{mol_id},{sdf_path},{type(e).__name__}:{e}"
            failed_entries.append(msg)
            print(f"{mol_id:<8} {sdf_path.name:<24} {'-':>6} {'-':>5} {'-':>8} {'-':>8} {'FAILED':>10}")

    with manifest.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "mol_id",
                "sdf_file",
                "xyz_file",
                "net_charge",
                "multiplicity",
                "n_atoms",
                "n_heavy_atoms",
            ]
        )
        for r in rows:
            w.writerow(
                [r.mol_id, r.sdf_file, r.xyz_file, r.net_charge, r.multiplicity, r.n_atoms, r.n_heavy_atoms]
            )

    if failed_entries:
        failed_log.write_text("\n".join(failed_entries) + "\n")
    else:
        if failed_log.exists():
            failed_log.unlink()

    print("-" * 78)
    print(f"Total SDF: {len(sdf_files)} | Success: {len(rows)} | Failed: {len(failed_entries)}")
    print(f"Manifest: {manifest}")
    if failed_entries:
        print(f"Failed log: {failed_log}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

