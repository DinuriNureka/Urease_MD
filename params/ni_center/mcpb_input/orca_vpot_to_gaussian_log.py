#!/usr/bin/env python3
"""
Build a Gaussian-like ESP log from ORCA .vpot for MCPB.py Step 3.

MCPB.py's parser (pymsmt.mol.gauio.get_esp_from_gau) only requires:
1) A marker line: "Electrostatic Properties Using The SCF Density"
2) "Atomic Center" lines with XYZ at fixed columns 32:62 (Angstrom)
3) "ESP Fit Center" lines with XYZ at fixed columns 32:62 (Angstrom)
4) A marker line: "Electrostatic Properties (Atomic Units)"
5) At least 6 following lines, then "Atom" and "Fit" lines where the
   last token is an ESP value (a.u.); counts must match centers.

The ORCA .vpot format already contains:
  line 1: "<n_atoms> <n_fit_points>"
  next n_atoms lines: atom-center coordinates in Bohr
  next n_fit_points lines: "<esp> <x> <y> <z>" in Bohr
"""

from __future__ import annotations

from pathlib import Path
import sys

B_TO_A = 0.529177210903


def read_vpot(path: Path) -> tuple[list[tuple[float, float, float]], list[tuple[float, float, float, float]]]:
    lines = path.read_text().splitlines()
    header = lines[0].split()
    n_atoms = int(header[0])
    n_fit = int(header[1])

    atom_coords_bohr = []
    for i in range(1, 1 + n_atoms):
        x, y, z = (float(v) for v in lines[i].split())
        atom_coords_bohr.append((x, y, z))

    fit_data_bohr = []
    start = 1 + n_atoms
    end = start + n_fit
    for i in range(start, end):
        esp, x, y, z = (float(v) for v in lines[i].split())
        fit_data_bohr.append((esp, x, y, z))

    if len(atom_coords_bohr) != n_atoms or len(fit_data_bohr) != n_fit:
        raise ValueError("Failed parsing .vpot sizes.")
    return atom_coords_bohr, fit_data_bohr


def to_ang(coords_bohr: tuple[float, float, float]) -> tuple[float, float, float]:
    x, y, z = coords_bohr
    return x * B_TO_A, y * B_TO_A, z * B_TO_A


def fmt_center_line(prefix: str, idx: int, xyz_ang: tuple[float, float, float]) -> str:
    # MCPB parser reads fixed slices [32:42], [42:52], [52:62].
    # Build line with explicit placement so coordinates are always aligned.
    x, y, z = xyz_ang
    line = [" "] * 90
    tag = f"{prefix:>19s}{idx:5d}"
    for i, ch in enumerate(tag):
        line[i] = ch
    x_s = f"{x:10.6f}"
    y_s = f"{y:10.6f}"
    z_s = f"{z:10.6f}"
    for i, ch in enumerate(x_s):
        line[32 + i] = ch
    for i, ch in enumerate(y_s):
        line[42 + i] = ch
    for i, ch in enumerate(z_s):
        line[52 + i] = ch
    return "".join(line).rstrip()


def build_fake_log(vpot_path: Path, out_log_path: Path) -> None:
    atom_coords_bohr, fit_data_bohr = read_vpot(vpot_path)

    out: list[str] = []
    out.append(" Entering Gaussian System, Link 0=g16")
    out.append(" ------------------------------------------------------------")
    out.append(" Electrostatic Properties Using The SCF Density")
    out.append(" ------------------------------------------------------------")

    for i, crd_bohr in enumerate(atom_coords_bohr, start=1):
        out.append(fmt_center_line("Atomic Center", i, to_ang(crd_bohr)))

    for i, (_, x, y, z) in enumerate(fit_data_bohr, start=1):
        out.append(fmt_center_line("ESP Fit Center", i, to_ang((x, y, z))))

    out.append(" ")
    out.append(" Electrostatic Properties (Atomic Units)")
    # Parser starts reading at line +6 after this marker.
    out.append(" ------------------------------------------------------------")
    out.append(" Center      Value")
    out.append(" ------------------------------------------------------------")
    out.append(" NOTE: synthetic Gaussian-compatible log from ORCA .vpot")
    out.append(" NOTE: atom ESP values below are placeholders for parser size check")

    # Atom ESP entries: parser only uses count; values are not written to .esp.
    for i, _ in enumerate(atom_coords_bohr, start=1):
        out.append(f" Atom {i:6d}   0.00000000E+00")

    for i, (esp, _, _, _) in enumerate(fit_data_bohr, start=1):
        out.append(f" Fit  {i:6d}  {esp: .8E}")

    out.append(" Normal termination of Gaussian 16")
    out_log_path.write_text("\n".join(out) + "\n")


def main() -> int:
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.vpot> <output.log>")
        return 1
    in_vpot = Path(sys.argv[1])
    out_log = Path(sys.argv[2])
    build_fake_log(in_vpot, out_log)
    print(f"Wrote Gaussian-like ESP log: {out_log}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
