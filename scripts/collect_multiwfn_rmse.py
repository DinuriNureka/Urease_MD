#!/usr/bin/env python3
"""
collect_multiwfn_rmse.py

Aggregate resp_multiwfn_metrics.txt from each molNN dir into one CSV.

RMSE and RRMSE are read from Multiwfn's post-fit summary line
  RMSE: <x>   RRMSE: <y>
in resp_multiwfn_resp.log (also stored in metrics as rmse_charge_e and rrmse).
Units: elementary charge (e) for RMSE; RRMSE is dimensionless (Multiwfn convention).

Default pass rules (aligned with Multiwfn RESP examples):
  rmse_charge_e < 0.005 e
  rrmse < 0.3
"""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


def parse_metrics(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for ln in path.read_text(errors="replace").splitlines():
        if "=" not in ln:
            continue
        k, v = ln.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results_dir", default="resp_results")
    ap.add_argument("--out_csv", default="resp_results/multiwfn_resp_rmse_report.csv")
    ap.add_argument(
        "--rmse_max_e",
        type=float,
        default=0.005,
        help="Pass if rmse_charge_e is strictly below this (e).",
    )
    ap.add_argument(
        "--rrmse_max",
        type=float,
        default=0.3,
        help="Pass if RRMSE is strictly below this.",
    )
    args = ap.parse_args()

    root = Path(args.results_dir)
    rows: list[list[object]] = []
    for d in sorted(root.glob("mol*/")):
        mol_id = d.name
        mpath = d / "resp_multiwfn_metrics.txt"
        if not mpath.is_file():
            rows.append(
                [mol_id, "", "", "False", "False", "INCOMPLETE"],
            )
            continue
        m = parse_metrics(mpath)
        rmse_charge_s = (m.get("rmse_charge_e") or "").strip() or (m.get("rmse_au") or "").strip()
        rrmse_s = (m.get("rrmse") or "").strip()
        try:
            rmse_v = float(rmse_charge_s) if rmse_charge_s else float("nan")
        except ValueError:
            rmse_v = float("nan")
        try:
            rrmse_v = float(rrmse_s) if rrmse_s else float("nan")
        except ValueError:
            rrmse_v = float("nan")

        pass_rmse = (
            rmse_v == rmse_v and rmse_v < args.rmse_max_e
        ) if rmse_charge_s else False
        pass_rrmse = (
            rrmse_v == rrmse_v and rrmse_v < args.rrmse_max
        ) if rrmse_s else False
        status = "OK"
        if not rmse_charge_s or not rrmse_s:
            status = "INCOMPLETE"
        elif not pass_rmse or not pass_rrmse:
            status = "WARN"

        rows.append(
            [
                mol_id,
                rmse_charge_s,
                rrmse_s,
                f"{pass_rmse}",
                f"{pass_rrmse}",
                status,
            ]
        )

    out = Path(args.out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "mol_id",
                "rmse_charge_e",
                "rrmse",
                "pass_rmse",
                "pass_rrmse",
                "status",
            ]
        )
        w.writerows(rows)

    print(f"Wrote {out}")

    # Summary statistics (valid numeric rows only)
    rmse_vals: list[float] = []
    rrmse_vals: list[float] = []
    n_ok = 0
    n_warn = 0
    n_incomplete = 0
    n_pass_both = 0
    for r in rows:
        st = r[5]
        if st == "OK":
            n_ok += 1
        elif st == "WARN":
            n_warn += 1
        else:
            n_incomplete += 1
        if r[3] == "True" and r[4] == "True":
            n_pass_both += 1
        rs, rr = r[1], r[2]
        if rs:
            try:
                rmse_vals.append(float(rs))
            except ValueError:
                pass
        if rr:
            try:
                rrmse_vals.append(float(rr))
            except ValueError:
                pass

    n_total = len(rows)
    print()
    print("========== Multiwfn RESP RMSE report summary ==========")
    print(f"Total molecules processed: {n_total}")
    if rmse_vals:
        print(f"RMSE (e): mean = {statistics.mean(rmse_vals):.6f}")
        print(f"RMSE (e): min  = {min(rmse_vals):.6f}")
        print(f"RMSE (e): max  = {max(rmse_vals):.6f}")
    else:
        print("RMSE (e): no valid values")
    if rrmse_vals:
        print(f"RRMSE:    mean = {statistics.mean(rrmse_vals):.6f}")
    else:
        print("RRMSE:    no valid values")
    print(
        f"Pass both thresholds (RMSE < {args.rmse_max_e} e, "
        f"RRMSE < {args.rrmse_max}): {n_pass_both}"
    )
    print(f"Status OK:     {n_ok}")
    print(f"Status WARN:   {n_warn}")
    print(f"Status INCOMPLETE: {n_incomplete}")
    print("=======================================================")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
