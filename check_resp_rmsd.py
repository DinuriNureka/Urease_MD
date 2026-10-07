#!/usr/bin/env python3
"""
check_resp_rmsd.py

Purpose:
  Validate geometry integrity for ORCA single-point jobs (RMSD should be near 0),
  extract SCF energies and HOMO/LUMO gap sanity metrics, and gate antechamber.

Inputs:
  --results_dir (default: resp_results)
  --inputs_dir  (default: resp_inputs)
  --manifest    (default: resp_inputs/job_manifest.csv)

Outputs:
  - resp_results/rmsd_check.csv
  - Colorized summary table on terminal (no color if non-TTY)
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path


def read_xyz(path: Path) -> list[tuple[str, float, float, float]]:
    lines = path.read_text().splitlines()
    n = int(lines[0].strip())
    atoms: list[tuple[str, float, float, float]] = []
    for ln in lines[2 : 2 + n]:
        t = ln.split()
        if len(t) < 4:
            continue
        atoms.append((t[0], float(t[1]), float(t[2]), float(t[3])))
    return atoms


def extract_orca_coords(resp_out: Path) -> list[tuple[str, float, float, float]]:
    lines = resp_out.read_text(errors="ignore").splitlines()
    starts = [i for i, ln in enumerate(lines) if "CARTESIAN COORDINATES (ANGSTROEM)" in ln]
    if not starts:
        return []
    i = starts[-1] + 2
    coords: list[tuple[str, float, float, float]] = []
    while i < len(lines):
        ln = lines[i].strip()
        if not ln:
            break
        t = lines[i].split()
        # Handle either "Elem idx x y z" or "Elem x y z"
        if len(t) >= 5 and t[1].isdigit():
            elem, x, y, z = t[0], float(t[2]), float(t[3]), float(t[4])
        elif len(t) >= 4:
            elem, x, y, z = t[0], float(t[1]), float(t[2]), float(t[3])
        else:
            i += 1
            continue
        coords.append((elem, x, y, z))
        i += 1
    return coords


def centroid(points: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    n = len(points)
    return (
        sum(p[0] for p in points) / n,
        sum(p[1] for p in points) / n,
        sum(p[2] for p in points) / n,
    )


def heavy_centroid_rmsd(
    a: list[tuple[str, float, float, float]],
    b: list[tuple[str, float, float, float]],
) -> tuple[float, int]:
    ah = [(e, x, y, z) for e, x, y, z in a if e.upper() != "H"]
    bh = [(e, x, y, z) for e, x, y, z in b if e.upper() != "H"]
    if len(ah) != len(bh) or not ah:
        return float("nan"), 0
    for i in range(len(ah)):
        if ah[i][0].upper() != bh[i][0].upper():
            return float("nan"), len(ah)
    ca = centroid([(x, y, z) for _, x, y, z in ah])
    cb = centroid([(x, y, z) for _, x, y, z in bh])
    s = 0.0
    for i in range(len(ah)):
        ax, ay, az = ah[i][1] - ca[0], ah[i][2] - ca[1], ah[i][3] - ca[2]
        bx, by, bz = bh[i][1] - cb[0], bh[i][2] - cb[1], bh[i][3] - cb[2]
        s += (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
    return math.sqrt(s / len(ah)), len(ah)


def extract_energy(resp_out: Path) -> float:
    energy = float("nan")
    for ln in resp_out.read_text(errors="ignore").splitlines():
        if "FINAL SINGLE POINT ENERGY" in ln:
            try:
                energy = float(ln.split()[-1])
            except Exception:
                pass
    return energy


def extract_homo_lumo(resp_out: Path) -> tuple[float, float, float]:
    txt = resp_out.read_text(errors="ignore").splitlines()
    idxs = [i for i, ln in enumerate(txt) if "ORBITAL ENERGIES" in ln]
    if not idxs:
        return float("nan"), float("nan"), float("nan")
    i = idxs[-1] + 1
    occ_energy: list[tuple[float, float]] = []
    num_re = re.compile(r"^\s*\d+\s+([-\d\.Ee+]+)\s+([-\d\.Ee+]+)")
    # Expect columns like: NO OCC E(Eh) E(eV)
    while i < len(txt):
        ln = txt[i]
        if not ln.strip():
            if occ_energy:
                break
            i += 1
            continue
        m = num_re.match(ln)
        if m:
            try:
                occ = float(m.group(1))
                e_eh = float(m.group(2))
                occ_energy.append((occ, e_eh))
            except Exception:
                pass
        i += 1
    if not occ_energy:
        return float("nan"), float("nan"), float("nan")
    homo = None
    lumo = None
    for occ, eh in occ_energy:
        if occ > 0.5:
            homo = eh
        elif lumo is None:
            lumo = eh
    if homo is None or lumo is None:
        return float("nan"), float("nan"), float("nan")
    homo_ev = homo * 27.2114
    lumo_ev = lumo * 27.2114
    return homo_ev, lumo_ev, lumo_ev - homo_ev


def colorize(s: str, status: str, use_color: bool) -> str:
    if not use_color:
        return s
    if status == "PASS":
        return f"\033[32m{s}\033[0m"
    if status in {"WARN"}:
        return f"\033[33m{s}\033[0m"
    return f"\033[31m{s}\033[0m"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results_dir", default="resp_results")
    ap.add_argument("--inputs_dir", default="resp_inputs")
    ap.add_argument("--manifest", default="resp_inputs/job_manifest.csv")
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    inputs_dir = Path(args.inputs_dir)
    manifest = Path(args.manifest)
    out_csv = results_dir / "rmsd_check.csv"
    results_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    with manifest.open() as fh:
        rdr = csv.DictReader(fh)
        for row in rdr:
            mol_id = row["mol_id"]
            orig = row["original_filename"]
            n_atoms = int(row["n_atoms"])
            n_heavy = 0
            inp_xyz = inputs_dir / f"{mol_id}.xyz"
            done = (results_dir / mol_id / "DONE").exists()
            resp_out = results_dir / mol_id / "resp.out"
            if not done:
                rows.append(
                    [mol_id, orig, n_atoms, "", "NaN", "ORCA_FAILED", "NaN", "NaN", "NaN", "NaN"]
                )
                continue
            a = read_xyz(inp_xyz)
            b = extract_orca_coords(resp_out)
            rmsd, n_heavy = heavy_centroid_rmsd(a, b)
            if math.isnan(rmsd):
                status = "FAIL"
            elif rmsd < 0.01:
                status = "PASS"
            elif rmsd < 0.1:
                status = "WARN"
            else:
                status = "FAIL"
            e = extract_energy(resp_out)
            h, l, g = extract_homo_lumo(resp_out)
            rows.append(
                [
                    mol_id,
                    orig,
                    n_atoms,
                    n_heavy,
                    f"{rmsd:.6f}" if not math.isnan(rmsd) else "NaN",
                    status,
                    f"{e:.10f}" if not math.isnan(e) else "NaN",
                    f"{h:.6f}" if not math.isnan(h) else "NaN",
                    f"{l:.6f}" if not math.isnan(l) else "NaN",
                    f"{g:.6f}" if not math.isnan(g) else "NaN",
                ]
            )

    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "mol_id",
                "original_filename",
                "n_atoms",
                "n_heavy_atoms",
                "rmsd_angstrom",
                "geometry_status",
                "scf_energy_Ha",
                "homo_eV",
                "lumo_eV",
                "gap_eV",
            ]
        )
        w.writerows(rows)

    use_color = sys.stdout.isatty()
    print(
        f"{'mol_id':<8} {'status':<12} {'RMSD(Å)':>10} {'E(Ha)':>18} {'gap(eV)':>10} {'n_heavy':>8}"
    )
    print("-" * 72)
    counts = {"PASS": 0, "WARN": 0, "FAIL": 0, "ORCA_FAILED": 0}
    ready = 0
    for r in rows:
        st = r[5]
        counts[st] = counts.get(st, 0) + 1
        if st in {"PASS", "WARN"}:
            ready += 1
        line = f"{r[0]:<8} {st:<12} {r[4]:>10} {r[6]:>18} {r[9]:>10} {str(r[3]):>8}"
        print(colorize(line, st, use_color))
    print("-" * 72)
    print(
        f"PASS: {counts['PASS']} | WARN: {counts['WARN']} | FAIL: {counts['FAIL']} | ORCA_FAILED: {counts['ORCA_FAILED']}"
    )
    print(f"Molecules ready for antechamber: {ready}")

    return 1 if (counts["FAIL"] > 0 or counts["ORCA_FAILED"] > 0) else 0


if __name__ == "__main__":
    raise SystemExit(main())

