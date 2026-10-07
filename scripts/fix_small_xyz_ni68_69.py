#!/usr/bin/env python3
"""Build urease_no_ohb_small.xyz with Ni at atoms 68-69 from MCPB Step 1 large-model xyz."""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "src",
        type=Path,
        help="Source xyz (e.g. urease_no_ohb_large.xyz from MCPB step 1)",
    )
    p.add_argument("dst", type=Path, help="Output xyz path")
    args = p.parse_args()

    lines = args.src.read_text().splitlines()
    natom = int(lines[0].strip())
    title = lines[1] if len(lines) > 1 else ""
    body = lines[2 : 2 + natom]
    if len(body) != natom:
        raise SystemExit(f"Expected {natom} atom lines, got {len(body)}")

    out_lines = [
        str(natom),
        "Small model — Ni at 68-69 (from MCPB step 1 geometry; NI1/NI2 sites)",
    ]
    for i, line in enumerate(body, start=1):
        parts = line.split()
        if not parts:
            continue
        if i in (68, 69):
            parts[0] = "Ni"
        out_lines.append(
            f"  {parts[0]:<3} {float(parts[1]):>18.15f} {float(parts[2]):>18.15f} {float(parts[3]):>18.15f}"
        )

    args.dst.parent.mkdir(parents=True, exist_ok=True)
    args.dst.write_text("\n".join(out_lines) + "\n")
    print(f"Wrote {args.dst} ({natom} atoms, Ni at 68-69)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
