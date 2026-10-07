#!/usr/bin/env python3
"""
calc_esp_rmse.py

Purpose:
  Validate RESP charge quality by comparing the quantum-mechanical molecular
  electrostatic potential (ESP) on a grid from ORCA with the classical ESP from
  RESP partial charges at the same grid points. Reports RMSE, max error, and R².

Inputs:
  - resp_results/molNN/resp.out (and optionally resp.vpot)
  - gaff2_params/molNN/molNN.mol2 (RESP charges + coordinates)
  - resp_inputs/job_manifest.csv (mol_id → original filename, atom counts)

Outputs:
  - resp_results/esp_rmse_report.csv
  - Summary table on stdout

Example:
  python3 calc_esp_rmse.py --results_dir resp_results --params_dir gaff2_params

For **Multiwfn RESP** fitting metrics (RMSE in e and RRMSE from the `RMSE: … RRMSE:` summary line), use
`scripts/collect_multiwfn_rmse.py` after `scripts/run_multiwfn_resp.sh`.

Theory:
  RMSE measures the root-mean-square deviation between V_QM and the potential
  from fitted point charges Σ q_j/|r−r_j| on the same surface/grid. Values
  below ~2 kcal/mol/e indicate an excellent fit; ~2–5 is good for routine MD;
  ~5–10 can be acceptable for large/flexible systems; above ~10 suggests the
  point-charge model is struggling and parameters should be reviewed.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Any

try:
    import numpy as np
except ImportError as e:
    raise RuntimeError("numpy is required: pip install numpy") from e

# Hartree (Eh) per electron → kcal/mol per unit charge (same numerical factor
# as Hartree → kcal/mol for energy).
HA_TO_KCAL_PER_E = 627.509

# 1 Å = 1.8897261254535... Bohr (CODATA); user specified 1.88973
ANG_TO_BOHR = 1.88973

# Minimum distance (Bohr) from any atom center to accept a grid point
MIN_DIST_BOHR = 0.1


def parse_vpot_file(path: Path) -> np.ndarray | None:
    """
    Parse ORCA companion file resp.vpot.

    ORCA 6.x (two integers on line 1):
      Line 1: NATOMS  NGRID
      Lines 2..NATOMS+1: atom X Y Z (Bohr)
      Next NGRID lines: V_ESP  X  Y  Z   (Hartree/e first, then Bohr)

    Legacy (single integer on line 1):
      Line 1: N (grid points only)
      Lines 2..N+1: X Y Z V (Bohr, Hartree/e)

    Returns:
      Array of shape (N, 4) with columns [x, y, z, v_qm] or None if invalid.
    """
    try:
        lines = path.read_text().splitlines()
    except OSError:
        return None
    if not lines:
        return None
    head = lines[0].split()
    try:
        if len(head) >= 2:
            nat = int(head[0])
            ngrid = int(head[1])
        else:
            nat = -1
            ngrid = int(head[0])
    except ValueError:
        return None

    rows: list[list[float]] = []

    if nat >= 0 and ngrid > 0:
        # ORCA 6: skip atom coordinate block, then read V x y z
        start = 1 + nat
        if len(lines) < start + ngrid:
            return None
        for i in range(start, start + ngrid):
            parts = lines[i].split()
            if len(parts) < 4:
                return None
            try:
                v, x, y, z = float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
            except ValueError:
                return None
            rows.append([x, y, z, v])
        return np.asarray(rows, dtype=np.float64)

    # Legacy: N grid lines of x y z v
    n = ngrid
    if n <= 0 or len(lines) < n + 1:
        return None
    for i in range(1, n + 1):
        parts = lines[i].split()
        if len(parts) < 4:
            return None
        try:
            x, y, z, v = float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            return None
        rows.append([x, y, z, v])
    return np.asarray(rows, dtype=np.float64)


def _find_esp_block_lines(lines: list[str]) -> tuple[int, int] | None:
    """
    Locate the ELECTROSTATIC POTENTIAL table in ORCA resp.out.

    Prefer a block whose preceding context mentions electric / ORCA property
    calculation; otherwise use the last matching header.
    Returns (start_data_idx, end_exclusive) for lines containing numeric data.
    """
    header_pat = re.compile(r"ELECTROSTATIC\s+POTENTIAL", re.IGNORECASE)
    candidates: list[int] = []
    for i, ln in enumerate(lines):
        if header_pat.search(ln.strip()):
            candidates.append(i)

    if not candidates:
        return None

    def score(idx: int) -> int:
        s = 0
        for j in range(max(0, idx - 12), idx):
            low = lines[j].upper()
            if "ORCA ELECTRIC PROPERTIES" in low or "ELECTRIC PROPERTIES" in low:
                s += 2
            if "PROPERTY" in low and "ELECTRIC" in low:
                s += 1
        return s

    best = max(candidates, key=score)
    if score(best) == 0:
        best = candidates[-1]

    i = best + 1
    # Skip separator dashes and column header
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i < len(lines) and re.match(r"^[-\s]+$", lines[i]):
        i += 1
    while i < len(lines):
        low = lines[i].upper()
        if "X(" in low and "BOHR" in low and "ESP" in low:
            i += 1
            break
        if re.match(r"^\s*\d", lines[i]) or (
            len(lines[i].split()) >= 4 and re.match(r"^[-\d.\sEe+]+$", lines[i].replace(" ", ""))
        ):
            break
        i += 1

    start = i
    end = start
    while end < len(lines):
        ln = lines[end].strip()
        if not ln:
            break
        parts = ln.split()
        if len(parts) >= 4:
            try:
                float(parts[0])
                float(parts[1])
                float(parts[2])
                float(parts[3])
            except ValueError:
                if end > start:
                    break
                end += 1
                continue
            end += 1
            continue
        if end > start:
            break
        end += 1
        break

    if end <= start:
        return None
    return start, end


def parse_esp_from_resp_out(text: str) -> np.ndarray | None:
    """
    Extract ESP grid from ORCA stdout.

    Looks for 'ELECTROSTATIC POTENTIAL' then rows: X Y Z ESP (Bohr, a.u.).
    """
    lines = text.splitlines()
    span = _find_esp_block_lines(lines)
    if span is None:
        return None
    start, end = span
    rows: list[list[float]] = []
    for k in range(start, end):
        parts = lines[k].split()
        if len(parts) < 4:
            continue
        try:
            x, y, z, v = float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
        rows.append([x, y, z, v])
    if not rows:
        return None
    return np.asarray(rows, dtype=np.float64)


def load_esp_grid(job_dir: Path) -> tuple[np.ndarray | None, str]:
    """
    Load QM ESP grid from resp.vpot (preferred) or embedded resp.out block.

    Returns:
      (array (N,4) or None, source label: 'vpot' | 'resp.out' | 'none')
    """
    vpot = job_dir / "resp.vpot"
    if vpot.is_file():
        arr = parse_vpot_file(vpot)
        if arr is not None:
            return arr, "vpot"
    out = job_dir / "resp.out"
    if out.is_file():
        arr = parse_esp_from_resp_out(out.read_text(errors="replace"))
        if arr is not None:
            return arr, "resp.out"
    return None, "none"


def parse_mol2_atoms(mol2_path: Path) -> tuple[np.ndarray, np.ndarray, int, float] | None:
    """
    Parse TRIPOS ATOM block: coordinates (Å) and partial charges (e).

    Returns:
      coords_ang (M,3), charges (M,), n_heavy, charge_sum
    """
    try:
        text = mol2_path.read_text(errors="replace")
    except OSError:
        return None
    if "@<TRIPOS>ATOM" not in text:
        return None
    parts = text.split("@<TRIPOS>ATOM")
    if len(parts) < 2:
        return None
    block = parts[1].split("@<TRIPOS>")[0].strip().splitlines()
    coords: list[list[float]] = []
    charges: list[float] = []
    n_heavy = 0
    for ln in block:
        if not ln.strip():
            continue
        tok = ln.split()
        if len(tok) < 9:
            continue
        try:
            x, y, z = float(tok[2]), float(tok[3]), float(tok[4])
            q = float(tok[8])
        except (ValueError, IndexError):
            continue
        coords.append([x, y, z])
        charges.append(q)
        el = tok[5].split(".")[0].split()[0].upper()
        if el != "H":
            n_heavy += 1
    if not coords:
        return None
    c = np.asarray(coords, dtype=np.float64)
    q = np.asarray(charges, dtype=np.float64)
    return c, q, n_heavy, float(np.sum(q))


def compute_metrics(
    grid_bohr_v: np.ndarray,
    coords_ang: np.ndarray,
    charges: np.ndarray,
) -> tuple[float, float, float, int] | None:
    """
    Vectorized ESP RMSE (Ha), max |delta| (Ha), R², and number of grid points used.
    """
    r_grid = grid_bohr_v[:, :3]
    v_qm = grid_bohr_v[:, 3]
    r_atoms = coords_ang * ANG_TO_BOHR

    # (N, M, 3)
    diff = r_grid[:, None, :] - r_atoms[None, :, :]
    dist = np.linalg.norm(diff, axis=2)
    dmin = np.min(dist, axis=1)
    mask = dmin >= MIN_DIST_BOHR
    if not np.any(mask):
        return None

    dist = dist[mask]
    v_qm = v_qm[mask]
    # Avoid div by zero (masked points already excluded)
    v_resp = np.sum(charges[None, :] / dist, axis=1)
    delta = v_qm - v_resp

    rmse = float(np.sqrt(np.mean(delta**2)))
    max_e = float(np.max(np.abs(delta)))
    mean_q = float(np.mean(v_qm))
    sst = float(np.sum((v_qm - mean_q) ** 2))
    if sst <= 0:
        r2 = float("nan")
    else:
        r2 = float(1.0 - np.sum(delta**2) / sst)
    npts = int(delta.shape[0])
    return rmse, max_e, r2, npts


def status_from_rmse(rmse_kcal: float) -> str:
    if rmse_kcal < 2.0:
        return "EXCELLENT"
    if rmse_kcal < 5.0:
        return "GOOD"
    if rmse_kcal < 10.0:
        return "ACCEPTABLE"
    return "POOR"


def color_line(line: str, status: str, use_color: bool) -> str:
    if not use_color:
        return line
    if status in {"EXCELLENT", "GOOD"}:
        return f"\033[32m{line}\033[0m"
    if status == "ACCEPTABLE":
        return f"\033[33m{line}\033[0m"
    return f"\033[31m{line}\033[0m"


def main() -> int:
    ap = argparse.ArgumentParser(description="ESP RMSE validation for RESP charges vs ORCA ESP grid.")
    ap.add_argument("--results_dir", default="resp_results")
    ap.add_argument("--params_dir", default="gaff2_params")
    ap.add_argument("--manifest", default="resp_inputs/job_manifest.csv")
    ap.add_argument("--out_csv", default="resp_results/esp_rmse_report.csv")
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    params_dir = Path(args.params_dir)
    manifest_path = Path(args.manifest)
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    with manifest_path.open(newline="") as fh:
        rdr = csv.DictReader(fh)
        manifest = sorted(list(rdr), key=lambda r: r["mol_id"])

    use_color = sys.stdout.isatty()

    print(
        f"{'mol_id':<8} {'orig':<22} {'RMSE':>8} {'max|d|':>8} {'R²':>7} {'pts':>6} {'status':<12}"
    )
    print("-" * 88)

    for row in manifest:
        mol_id = row["mol_id"]
        orig = row.get("original_filename", "")
        try:
            n_atoms = int(row.get("n_atoms", 0))
        except ValueError:
            n_atoms = 0

        mol2 = params_dir / mol_id / f"{mol_id}.mol2"
        job_dir = results_dir / mol_id

        if not mol2.is_file():
            rows.append(
                {
                    "mol_id": mol_id,
                    "original_filename": orig,
                    "n_atoms": n_atoms,
                    "n_heavy_atoms": "",
                    "n_esp_points": "",
                    "rmse_kcal_mol_per_e": "",
                    "max_error_kcal_mol_per_e": "",
                    "r_squared": "",
                    "charge_sum": "",
                    "status": "NO_MOL2",
                }
            )
            line = f"{mol_id:<8} {orig[:20]:<22} {'—':>8} {'—':>8} {'—':>7} {'—':>6} {'NO_MOL2':<12}"
            print(color_line(line, "NO_MOL2", use_color))
            continue

        parsed = parse_mol2_atoms(mol2)
        if parsed is None:
            rows.append(
                {
                    "mol_id": mol_id,
                    "original_filename": orig,
                    "n_atoms": n_atoms,
                    "n_heavy_atoms": "",
                    "n_esp_points": "",
                    "rmse_kcal_mol_per_e": "",
                    "max_error_kcal_mol_per_e": "",
                    "r_squared": "",
                    "charge_sum": "",
                    "status": "NO_MOL2",
                }
            )
            line = f"{mol_id:<8} {orig[:20]:<22} {'—':>8} {'—':>8} {'—':>7} {'—':>6} {'NO_MOL2':<12}"
            print(color_line(line, "NO_MOL2", use_color))
            continue

        coords_ang, charges, n_heavy, charge_sum = parsed

        grid, src = load_esp_grid(job_dir)
        if grid is None:
            rows.append(
                {
                    "mol_id": mol_id,
                    "original_filename": orig,
                    "n_atoms": n_atoms,
                    "n_heavy_atoms": n_heavy,
                    "n_esp_points": "",
                    "rmse_kcal_mol_per_e": "",
                    "max_error_kcal_mol_per_e": "",
                    "r_squared": "",
                    "charge_sum": f"{charge_sum:.6f}",
                    "status": "NO_ESP_DATA",
                }
            )
            line = f"{mol_id:<8} {orig[:20]:<22} {'—':>8} {'—':>8} {'—':>7} {'—':>6} {'NO_ESP_DATA':<12}"
            print(color_line(line, "NO_ESP_DATA", use_color))
            continue

        met = compute_metrics(grid, coords_ang, charges)
        if met is None:
            rows.append(
                {
                    "mol_id": mol_id,
                    "original_filename": orig,
                    "n_atoms": n_atoms,
                    "n_heavy_atoms": n_heavy,
                    "n_esp_points": 0,
                    "rmse_kcal_mol_per_e": "",
                    "max_error_kcal_mol_per_e": "",
                    "r_squared": "",
                    "charge_sum": f"{charge_sum:.6f}",
                    "status": "NO_ESP_DATA",
                }
            )
            line = f"{mol_id:<8} {orig[:20]:<22} {'—':>8} {'—':>8} {'—':>7} {'0':>6} {'NO_ESP_DATA':<12}"
            print(color_line(line, "NO_ESP_DATA", use_color))
            continue

        rmse_ha, max_ha, r2, npts = met
        rmse_k = rmse_ha * HA_TO_KCAL_PER_E
        max_k = max_ha * HA_TO_KCAL_PER_E
        st = status_from_rmse(rmse_k)

        rows.append(
            {
                "mol_id": mol_id,
                "original_filename": orig,
                "n_atoms": n_atoms,
                "n_heavy_atoms": n_heavy,
                "n_esp_points": npts,
                "rmse_kcal_mol_per_e": f"{rmse_k:.6f}",
                "max_error_kcal_mol_per_e": f"{max_k:.6f}",
                "r_squared": f"{r2:.6f}" if not np.isnan(r2) else "",
                "charge_sum": f"{charge_sum:.6f}",
                "status": st,
            }
        )
        r2s = f"{r2:.4f}" if not np.isnan(r2) else "nan"
        line = (
            f"{mol_id:<8} {orig[:20]:<22} {rmse_k:8.3f} {max_k:8.3f} {r2s:>7} {npts:6d} {st:<12}"
        )
        print(color_line(line, st, use_color))

    print("-" * 88)

    counts: dict[str, int] = {
        "EXCELLENT": 0,
        "GOOD": 0,
        "ACCEPTABLE": 0,
        "POOR": 0,
        "NO_ESP_DATA": 0,
        "NO_MOL2": 0,
    }
    rmses: list[float] = []
    for r in rows:
        st = r["status"]
        counts[st] = counts.get(st, 0) + 1
        if st in {"EXCELLENT", "GOOD", "ACCEPTABLE", "POOR"} and r.get("rmse_kcal_mol_per_e"):
            rmses.append(float(r["rmse_kcal_mol_per_e"]))

    n_total = len(rows)
    no_data = counts["NO_ESP_DATA"] + counts["NO_MOL2"]
    mean_rmse = float(np.mean(rmses)) if rmses else float("nan")
    ready = counts["EXCELLENT"] + counts["GOOD"] + counts["ACCEPTABLE"]

    print(
        f"EXCELLENT: {counts['EXCELLENT']} | GOOD: {counts['GOOD']} | "
        f"ACCEPTABLE: {counts['ACCEPTABLE']} | POOR: {counts['POOR']} | "
        f"NO_DATA: {counts['NO_ESP_DATA'] + counts['NO_MOL2']}"
    )
    if not np.isnan(mean_rmse):
        print(f"Mean RMSE: {mean_rmse:.2f} kcal/mol/e  (over molecules with data)")
    else:
        print("Mean RMSE: n/a  (no RMSE rows)")
    print(f"Molecules ready for MD: {ready}  (EXCELLENT + GOOD + ACCEPTABLE)")

    fieldnames = [
        "mol_id",
        "original_filename",
        "n_atoms",
        "n_heavy_atoms",
        "n_esp_points",
        "rmse_kcal_mol_per_e",
        "max_error_kcal_mol_per_e",
        "r_squared",
        "charge_sum",
        "status",
    ]
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"\nWrote: {out_csv}")

    poor = counts["POOR"] > 0
    frac_no_esp = counts["NO_ESP_DATA"] / n_total if n_total else 0.0
    any_no_mol2 = counts["NO_MOL2"] > 0
    # Spec: fail on any POOR, or if >20% of molecules lack ESP grid data, or if
    # RESP mol2 outputs are missing (cannot validate charges at all).
    exit1 = poor or (frac_no_esp > 0.20) or any_no_mol2
    return 1 if exit1 else 0


if __name__ == "__main__":
    raise SystemExit(main())


# INTERPRETATION GUIDE
# RMSE < 2   → excellent fit, charges reliably reproduce QM electrostatics
# RMSE 2-5   → good fit, suitable for standard MD
# RMSE 5-10  → acceptable for flexible/large molecules, check manually
# RMSE > 10  → poor fit, consider: wrong charge state, poor geometry,
#              or unusual electronic structure — do not use without review
