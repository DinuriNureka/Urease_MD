#!/usr/bin/env python3
"""
Parse Multiwfn stdout/stderr from RESP fitting and write:
  - resp_multiwfn.chg   (1-based index, charge per line after header)
  - resp_multiwfn_metrics.txt (rmse_au=..., rmse_charge_e=..., rrmse=...)

Multiwfn text varies by version; this uses flexible regex heuristics.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def parse_rmse_charge_e(text: str) -> float | None:
    """
    RESP summary line after the charge table (Multiwfn prints no ESP energy unit here):
        RMSE:    0.002217   RRMSE:    0.120512
    Sir/advisor interpretation: RMSE is in elementary charge (e), not Hartree or kcal/mol.
    """
    m = re.search(
        r"^\s*RMSE:\s+([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s+RRMSE:",
        text,
        re.MULTILINE | re.IGNORECASE,
    )
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def parse_rmse_rrmse(text: str) -> tuple[float | None, float | None]:
    rmse = None
    rrmse = None
    # Examples: "RMSE ... 0.001234 a.u."  "RRMSE ... 0.25"
    for pat in (
        r"RMSE[^\d\-]*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)",
        r"root\s+mean\s+square\s+error[^\d\-]*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)",
    ):
        m = re.search(pat, text, re.IGNORECASE)
        if m:
            try:
                rmse = float(m.group(1))
                break
            except ValueError:
                pass
    for pat in (
        r"RRMSE[^\d\-]*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)",
        r"relative\s+root\s+mean\s+square\s+error[^\d\-]*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)",
        r"relative\s+RMSE[^\d\-]*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)",
    ):
        m = re.search(pat, text, re.IGNORECASE)
        if m:
            try:
                rrmse = float(m.group(1))
                break
            except ValueError:
                pass
    return rmse, rrmse


def parse_atomic_charges(text: str) -> list[tuple[int, float]] | None:
    """
    Try to extract (1-based index, charge) in atom order.
    """
    lines = text.splitlines()
    charges: list[tuple[int, float]] = []

    # Strategy A: block after "Charge" header with two columns index charge
    in_block = False
    for ln in lines:
        low = ln.lower()
        if "charge" in low and ("atom" in low or "center" in low or "final" in low):
            in_block = True
            continue
        if in_block:
            if not ln.strip():
                if charges:
                    break
                continue
            # "   1(C )    0.123456"  or "  1     0.123"
            m = re.match(
                r"^\s*(\d+)\s*(?:\([^)]*\))?\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*$", ln
            )
            if m:
                charges.append((int(m.group(1)), float(m.group(2))))
                continue
            m2 = re.match(
                r"^\s*(\d+)\s+[A-Za-z]+\s+([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)", ln
            )
            if m2:
                charges.append((int(m2.group(1)), float(m2.group(2))))
                continue
            if re.match(r"^\s*-{3,}", ln) or ln.strip().startswith("==="):
                if charges:
                    break

    if len(charges) >= 3:
        charges.sort(key=lambda x: x[0])
        return charges

    # Strategy B: lines like " Atom 1 Charge : 0.123"
    charges = []
    for ln in lines:
        m = re.search(
            r"(?:Atom|Center)\s+(\d+).*?([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)", ln, re.I
        )
        if m:
            charges.append((int(m.group(1)), float(m.group(2))))
    if len(charges) >= 3:
        charges.sort(key=lambda x: x[0])
        return charges

    # Strategy C: last table of floats in RESP section
    return None


def parse_charges_sequential(text: str) -> list[float] | None:
    """If indices missing, collect floats that look like charges after RESP keyword."""
    if "RESP" not in text.upper():
        return None
    part = text.upper().split("RESP", 1)[-1]
    floats: list[float] = []
    for ln in part.splitlines():
        for tok in ln.split():
            try:
                v = float(tok)
            except ValueError:
                continue
            if -4.0 < v < 4.0 and (tok not in ("0", "1", "2", "3")):
                floats.append(v)
    if 5 <= len(floats) <= 500:
        return floats
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("log_file", help="Multiwfn RESP log (stdout capture)")
    ap.add_argument("--out-chg", default="resp_multiwfn.chg")
    ap.add_argument("--out-metrics", default="resp_multiwfn_metrics.txt")
    args = ap.parse_args()

    path = Path(args.log_file)
    text = path.read_text(errors="replace")

    rmse, rrmse = parse_rmse_rrmse(text)
    rmse_charge_e = parse_rmse_charge_e(text)
    if rmse_charge_e is None and rmse is not None:
        rmse_charge_e = rmse
    pairs = parse_atomic_charges(text)
    seq: list[float] | None = None

    out_chg = Path(args.out_chg)
    out_met = Path(args.out_metrics)

    if pairs:
        pairs.sort(key=lambda x: x[0])
        lines_out = ["# Multiwfn RESP charges (1-based index, charge)", str(len(pairs))]
        for idx, q in pairs:
            lines_out.append(f"{idx} {q:.8f}")
        out_chg.write_text("\n".join(lines_out) + "\n")
    else:
        seq = parse_charges_sequential(text)
        if seq:
            lines_out = ["# Multiwfn RESP charges (sequential order)", str(len(seq))]
            for i, q in enumerate(seq, start=1):
                lines_out.append(f"{i} {q:.8f}")
            out_chg.write_text("\n".join(lines_out) + "\n")
        else:
            out_chg.write_text(
                "# PARSE_FAILED — inspect Multiwfn log; edit charges manually\n"
            )

    parse_ok = bool(pairs) or bool(seq)
    with out_met.open("w") as fh:
        fh.write(f"rmse_au={rmse if rmse is not None else ''}\n")
        fh.write(
            f"rmse_charge_e={rmse_charge_e if rmse_charge_e is not None else ''}\n"
        )
        fh.write(f"rrmse={rrmse if rrmse is not None else ''}\n")
        fh.write(f"parse_charges={'OK' if parse_ok else 'FAIL'}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
