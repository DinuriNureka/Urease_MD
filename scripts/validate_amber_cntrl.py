#!/usr/bin/env python3
"""Validate key=value pairs in the first &cntrl block of an AMBER mdin file.

AMBER 24 pmemd: ntt=11 is the Bussi (CSVR) thermostat (use tautp) — valid for equilibration
and production. ntt=3 is Langevin (requires gamma_ln). Production with Bussi must not set
gamma_ln.

Usage:
  validate_amber_cntrl.py <file.in> key=expected [key=expected ...]
  validate_amber_cntrl.py simulations/input_files/production.in
      (shorthand: runs default production checks; fails if gamma_ln is present)

Exits 0 if all match, 1 with stderr messages if any check fails.
Comments (!) are stripped before parsing. Float keys (dt) compared numerically.
"""
from __future__ import annotations

import os
import re
import sys


def strip_comments(text: str) -> str:
    out = []
    for line in text.splitlines():
        if "!" in line:
            line = line.split("!", 1)[0]
        out.append(line)
    return "\n".join(out)


def parse_first_cntrl(path: str) -> dict[str, str]:
    raw = open(path, encoding="utf-8", errors="replace").read()
    chunk = strip_comments(raw)
    m = re.search(r"&cntrl\s*(.*?)\s*/", chunk, re.S | re.I)
    if not m:
        raise ValueError(f"No &cntrl ... / block found in {path}")
    body = m.group(1).replace("\n", " ")
    pairs: dict[str, str] = {}
    # Split on commas; tokens look like "imin=0" or "nstlim=50000"
    for tok in body.split(","):
        tok = tok.strip()
        if not tok or "=" not in tok:
            continue
        k, v = tok.split("=", 1)
        k = k.strip().lower()
        v = v.strip()
        if k not in pairs:
            pairs[k] = v
    return pairs


def norm_pair(key: str, expected: str, got: str) -> bool:
    key = key.lower()
    if key == "dt":
        return abs(float(got) - float(expected)) < 1e-9
    # Thermostat/barostat and other floats commonly compared in gates
    float_keys = (
        "tautp",
        "taup",
        "gamma_ln",
        "pres0",
        "cut",
        "temp0",
        "tempi",
        "restraint_wt",
        "mcbarint",
    )
    if key in float_keys:
        try:
            return abs(float(got) - float(expected)) < 1e-5
        except ValueError:
            return False
    # integers / enums (ntt=11 Bussi, ntc=2, etc.)
    return got.strip() == expected.strip()


def _default_production_checks(path: str) -> tuple[list[tuple[str, str]], bool]:
    """Return (required key=expected list, forbid_gamma_ln)."""
    base = os.path.basename(path)
    if base != "production.in":
        return [], False
    checks = [
        ("dt", "0.002"),
        ("nstlim", "10000000"),
        ("ntc", "2"),
        ("ntf", "2"),
        ("barostat", "2"),
        ("mcbarint", "25"),
        ("ntt", "11"),
        ("tautp", "2.0"),
    ]
    return checks, True


def main() -> int:
    if len(sys.argv) < 2:
        print(
            "Usage: validate_amber_cntrl.py <file.in> [key=expected ...]\n"
            "       validate_amber_cntrl.py simulations/input_files/production.in",
            file=sys.stderr,
        )
        return 2

    path = sys.argv[1]
    forbid_gamma = False
    if len(sys.argv) == 2:
        checks, forbid_gamma = _default_production_checks(path)
        if not checks:
            print(
                "Usage: validate_amber_cntrl.py <file.in> key=expected ...",
                file=sys.stderr,
            )
            return 2
    else:
        checks = []
        for arg in sys.argv[2:]:
            if "=" not in arg:
                print(f"Bad argument (need key=value): {arg}", file=sys.stderr)
                return 2
            k, v = arg.split("=", 1)
            checks.append((k.strip(), v.strip()))

    try:
        pairs = parse_first_cntrl(path)
    except Exception as e:
        print(f"[GATE FAIL] {e}", file=sys.stderr)
        return 1

    errors: list[str] = []
    if forbid_gamma and "gamma_ln" in pairs:
        errors.append(
            f"{path}: gamma_ln must not appear in production.in with Bussi (ntt=11); "
            f"remove gamma_ln (Langevin-only)."
        )

    for key, expected in checks:
        kl = key.lower()
        got = pairs.get(kl)
        if got is None:
            errors.append(f"{path}: &cntrl missing '{kl}' (required {kl}={expected})")
            continue
        if not norm_pair(kl, expected, got):
            errors.append(f"{path}: {kl} expected {expected!r} got {got!r}")
    if errors:
        for e in errors:
            print(f"[GATE FAIL] {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
