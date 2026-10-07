#!/usr/bin/env python3
"""
prepare_resp_inputs.py

Purpose:
  Validate ligand-only XYZ files, perform charge sanity checks, copy validated
  XYZ files into a RESP working directory, and create a manifest CSV.

Inputs:
  --xyz_dir   Input directory containing *.xyz ligand files
  --out_dir   Output directory for validated XYZ + logs + manifest

Outputs:
  - <out_dir>/molNN.xyz
  - <out_dir>/job_manifest.csv
  - <out_dir>/prepare.log

Example:
  python3 prepare_resp_inputs.py --xyz_dir /path/to/data/free_ligands
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable

ALLOWED_ELEMENTS = {"C", "H", "N", "O", "S", "F", "CL", "BR", "I", "P"}


def norm_elem(elem: str) -> str:
    e = elem.strip()
    if not e:
        return e
    if len(e) == 1:
        return e.upper()
    return e[0].upper() + e[1:].lower()


def validate_xyz(path: Path) -> tuple[bool, str, list[tuple[str, float, float, float]]]:
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        return False, "Too few lines for XYZ", []
    try:
        n_atoms = int(lines[0].strip())
    except Exception:
        return False, "Line 1 is not a valid positive integer atom count", []
    if n_atoms <= 0:
        return False, "Atom count must be positive", []
    if len(lines) != n_atoms + 2:
        return False, f"Expected {n_atoms + 2} lines but found {len(lines)}", []
    if n_atoms < 5:
        return False, "N < 5 (likely corrupt/empty ligand)", []

    atoms: list[tuple[str, float, float, float]] = []
    for i, ln in enumerate(lines[2:], start=3):
        parts = ln.split()
        if len(parts) != 4:
            return False, f"Line {i}: expected 4 fields, got {len(parts)}", []
        elem = parts[0]
        if not elem.isalpha():
            return False, f"Line {i}: element token is not alphabetic ({elem})", []
        e = norm_elem(elem)
        if e.upper() not in ALLOWED_ELEMENTS:
            return False, f"Line {i}: disallowed element {e} (metal/non-organic found)", []
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except Exception:
            return False, f"Line {i}: coordinates are not valid floats", []
        atoms.append((e, x, y, z))
    return True, "OK", atoms


def detect_charge_rdkit(xyz_path: Path, default_charge: int = 0) -> tuple[int, str]:
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except Exception:
        return default_charge, "RDKit not found, defaulting charge=0"

    mol = Chem.MolFromXYZFile(str(xyz_path))
    if mol is None:
        return default_charge, "RDKit failed to parse XYZ, defaulting charge=0"
    try:
        rdDetermineBonds.DetermineBonds(mol, charge=0)
    except Exception as e:
        return default_charge, f"RDKit DetermineBonds failed ({e}); defaulting charge=0"
    detected = sum(int(a.GetFormalCharge()) for a in mol.GetAtoms())
    return detected, "RDKit charge detection OK"


def write_xyz(path: Path, atoms: Iterable[tuple[str, float, float, float]], comment: str) -> None:
    atom_list = list(atoms)
    with path.open("w") as fh:
        fh.write(f"{len(atom_list)}\n")
        fh.write(f"{comment}\n")
        for e, x, y, z in atom_list:
            fh.write(f"{e:2s} {x: .10f} {y: .10f} {z: .10f}\n")


def log(log_path: Path, msg: str) -> None:
    with log_path.open("a") as fh:
        fh.write(msg + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--xyz_dir",
        default="/path/to/data/free_ligands",
        help="Input directory with ligand-only XYZ files",
    )
    ap.add_argument("--out_dir", default="resp_inputs", help="Output directory")
    args = ap.parse_args()

    xyz_dir = Path(args.xyz_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "prepare.log"
    manifest_path = out_dir / "job_manifest.csv"
    if log_path.exists():
        log_path.unlink()

    all_xyz_files = sorted(xyz_dir.glob("*.xyz"))
    # Default behavior: ignore trajectory snapshot files such as *_trj.xyz.
    xyz_files = [p for p in all_xyz_files if not p.stem.lower().endswith("_trj")]
    if not xyz_files:
        raise SystemExit(f"No *.xyz files found in {xyz_dir}")

    rows: list[list[str]] = []
    fail_count = 0

    header = [
        "mol_id",
        "original_filename",
        "xyz_file",
        "net_charge",
        "multiplicity",
        "n_atoms",
        "validation_status",
    ]
    skipped = len(all_xyz_files) - len(xyz_files)
    if skipped > 0:
        log(log_path, f"Skipped {skipped} files matching *_trj.xyz")

    print(
        f"{'mol_id':<8} {'original_filename':<28} {'charge':>6} {'mult':>5} {'n_atoms':>8} {'status':>12}"
    )
    print("-" * 74)

    for i, xyz in enumerate(xyz_files, start=1):
        mol_id = f"mol{i:02d}"
        out_xyz = out_dir / f"{mol_id}.xyz"
        ok, msg, atoms = validate_xyz(xyz)
        n_atoms = len(atoms)
        if not ok:
            status = "FAIL_FORMAT"
            fail_count += 1
            net_charge = 0
            log(log_path, f"{mol_id} {xyz.name}: FAIL_FORMAT: {msg}")
        else:
            detected_charge, charge_msg = detect_charge_rdkit(xyz, default_charge=0)
            net_charge = int(detected_charge)
            status = "OK"
            if n_atoms > 120:
                status = "WARN_SIZE"
            if detected_charge != 0:
                status = "WARN_CHARGE"
                log(
                    log_path,
                    f"{mol_id} {xyz.name}: RDKit detected charge={detected_charge}, expected 0 (verify manually)",
                )
            if "defaulting charge=0" in charge_msg:
                log(log_path, f"{mol_id} {xyz.name}: {charge_msg}")
            write_xyz(out_xyz, atoms, f"validated copy of {xyz.name}")
            log(log_path, f"{mol_id} {xyz.name}: {status}, n_atoms={n_atoms}, net_charge={net_charge}")

        rows.append(
            [
                mol_id,
                xyz.name,
                str(out_xyz),
                str(net_charge),
                "1",
                str(n_atoms),
                status,
            ]
        )
        print(f"{mol_id:<8} {xyz.name:<28} {net_charge:>6} {1:>5} {n_atoms:>8} {status:>12}")

    with manifest_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)

    print("-" * 74)
    if skipped > 0:
        print(f"Skipped *_trj.xyz files: {skipped}")
    print(f"Manifest written: {manifest_path}")
    print(f"Log written: {log_path}")
    print(f"Total: {len(rows)} | FAIL_FORMAT: {fail_count}")
    return 1 if fail_count > 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())

