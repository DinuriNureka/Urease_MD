#!/usr/bin/env python3
"""
Verify 36/36 mol2 updated, charges sum to net charge from manifest (default 0).

Optional soft check on RMSE in elementary charge (rmse_charge_e in report CSV):
values >= threshold print a warning (default 0.1 e) but do not fail the script.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def mol2_charge_sum(path: Path) -> float:
    text = path.read_text(errors="replace")
    if "@<TRIPOS>ATOM" not in text:
        return float("nan")
    block = text.split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>")[0]
    s = 0.0
    for ln in block.splitlines():
        tok = ln.split()
        if len(tok) >= 9:
            s += float(tok[8])
    return s


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gaff2_dir", default="gaff2_params")
    ap.add_argument("--manifest", default="resp_inputs/job_manifest.csv")
    ap.add_argument("--report", default="resp_results/multiwfn_resp_rmse_report.csv")
    ap.add_argument(
        "--rmse_charge_warn_max",
        type=float,
        default=0.1,
        help="Soft warning if rmse_charge_e >= this (e); does not change exit code.",
    )
    args = ap.parse_args()

    gaff2 = Path(args.gaff2_dir)
    manifest = Path(args.manifest)
    net_by_id: dict[str, int] = {}
    with manifest.open() as fh:
        for row in csv.DictReader(fh):
            net_by_id[row["mol_id"]] = int(row["net_charge"])

    n_ok = 0
    n_bad = 0
    for mol_id in sorted(net_by_id.keys()):
        p = gaff2 / mol_id / f"{mol_id}.mol2"
        if not p.is_file():
            print(f"MISSING {p}")
            n_bad += 1
            continue
        sm = mol2_charge_sum(p)
        expected = float(net_by_id[mol_id])
        if abs(sm - expected) < 1e-3:
            n_ok += 1
        else:
            print(f"CHARGE_SUM {mol_id}: got {sm} expected {expected}")
            n_bad += 1

    print(f"mol2 present + charge OK: {n_ok}/{len(net_by_id)}")
    rep = Path(args.report)
    if rep.is_file():
        print(f"RMSE report: {rep} ({rep.stat().st_size} bytes)")
        n_rmse_warn = 0
        with rep.open() as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames and "rmse_charge_e" in reader.fieldnames:
                for row in reader:
                    mol = row.get("mol_id", "")
                    raw = (row.get("rmse_charge_e") or "").strip()
                    if not raw:
                        continue
                    try:
                        v = float(raw)
                    except ValueError:
                        continue
                    if v >= args.rmse_charge_warn_max:
                        print(
                            f"WARNING rmse_charge_e: {mol} value {v} >= soft threshold "
                            f"{args.rmse_charge_warn_max} (e)"
                        )
                        n_rmse_warn += 1
                if n_rmse_warn:
                    print(
                        f"rmse_charge_e soft warnings: {n_rmse_warn} "
                        f"(threshold {args.rmse_charge_warn_max} e)"
                    )
            else:
                print(
                    "NOTE: report has no rmse_charge_e column; "
                    "re-run collect_multiwfn_rmse.py after refreshing metrics."
                )
    else:
        print(f"WARNING: missing {rep}")

    return 0 if n_bad == 0 and n_ok == len(net_by_id) else 1


if __name__ == "__main__":
    raise SystemExit(main())
