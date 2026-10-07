#!/usr/bin/env python3
"""
Scan Amber mdout (.out) files under simulations/*/ for GPU binary, NaN lines, and max clean NSTEP.

Assumed project root: parent of scripts/ (override with PROJECT_ROOT).
"""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

# Default: .../MD when script lives in .../MD/scripts/
ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parent.parent))
SIM_DIR = ROOT / "simulations"
LOG_DIR = ROOT / "logs"
OUT_LOG = LOG_DIR / "audit_pmemd_outs.tsv"

EXEC_RE = re.compile(r"Executable\s*(?:path)?\s*:\s*(\S+)", re.IGNORECASE)
NSTEP_RE = re.compile(r"NSTEP\s*=\s*(\d+)")


def classify_binary(path_str: str) -> str:
    """Return 'old', 'cuda128', or 'other' for summary counts."""
    if not path_str:
        return "other"
    p = path_str.replace("\\", "/")
    if "pmemd24_cuda128" in p:
        return "cuda128"
    if "/pmemd24/bin/" in p or p.endswith("pmemd24/bin/pmemd.cuda"):
        return "old"
    return "other"


def analyze_out(path: Path) -> tuple[str, int, int, bool]:
    """
    Returns (binary, nan_count, max_clean_nstep, completed).
    binary may be empty if not found in header.
    """
    binary = ""
    nan_count = 0
    max_clean = 0

    try:
        raw = path.read_bytes()
    except OSError:
        return "", 0, 0, False

    text = raw.decode("utf-8", errors="replace")
    lines = text.splitlines()

    for i, line in enumerate(lines[:50]):
        if "Executable" in line and ":" in line:
            m = EXEC_RE.search(line)
            if m:
                binary = m.group(1).strip()
                break

    for line in lines:
        if "NaN" in line:
            nan_count += 1

    for line in lines:
        if "NaN" in line:
            continue
        if "NSTEP" not in line or "=" not in line:
            continue
        m = NSTEP_RE.search(line)
        if m:
            max_clean = max(max_clean, int(m.group(1)))

    completed = max_clean >= 2_490_000
    return binary, nan_count, max_clean, completed


def main() -> int:
    if not SIM_DIR.is_dir():
        print(f"[FAIL] Missing simulations directory: {SIM_DIR}", file=sys.stderr)
        return 1

    rows: list[tuple[str, str, str, int, int, bool]] = []

    for out_path in sorted(SIM_DIR.rglob("*.out")):
        if not out_path.is_file():
            continue
        try:
            rel = out_path.relative_to(SIM_DIR)
        except ValueError:
            continue
        parts = rel.parts
        if len(parts) < 2:
            continue
        # system = top-level folder under simulations/ (e.g. mol1_pose2_monomer, md_output)
        # segment = immediate parent directory of the .out file (e.g. equil1, heat)
        system = parts[0]
        segment = parts[-2]

        binary, nan_count, max_nstep, completed = analyze_out(out_path)
        rows.append((system, segment, binary or "", nan_count, max_nstep, completed))

    rows.sort(key=lambda r: (r[0], r[1]))

    lines_out: list[str] = []
    header = "system\tsegment\tbinary\tnan_count\tmax_nstep\tcompleted"
    lines_out.append(header)

    old_ct = 0
    cuda128_ct = 0
    other_bin_ct = 0
    nan_files = 0

    for system, segment, binary, nan_count, max_nstep, completed in rows:
        bcls = classify_binary(binary)
        if bcls == "old":
            old_ct += 1
        elif bcls == "cuda128":
            cuda128_ct += 1
        else:
            other_bin_ct += 1
        if nan_count > 0:
            nan_files += 1

        # TSV: avoid tabs inside fields
        safe_bin = binary.replace("\t", " ")
        comp_str = "1" if completed else "0"
        lines_out.append(
            f"{system}\t{segment}\t{safe_bin}\t{nan_count}\t{max_nstep}\t{comp_str}"
        )

    table_text = "\n".join(lines_out) + "\n"
    sys.stdout.write(table_text)

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    OUT_LOG.write_text(table_text, encoding="utf-8")

    n = len(rows)
    sys.stdout.write(
        "\n--- Summary ---\n"
        f"Total .out files scanned: {n}\n"
        f"Files with nan_count > 0: {nan_files}\n"
        f"Files (header binary) — pmemd24/bin (old): {old_ct}\n"
        f"Files (header binary) — pmemd24_cuda128: {cuda128_ct}\n"
        f"Files (header binary) — other/unknown: {other_bin_ct}\n"
        f"Wrote: {OUT_LOG}\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
