#!/usr/bin/env python3
"""
Batch builder for full dodecamer systems (molX_poseY_full).

For each ligand XYZ in params/inhibitors/structures/mol*_pose*_lig.xyz:
  - select the lowest pose index per molX (configurable)
  - write ligand PDB (HETATM) for antechamber
  - run antechamber (GAFF2 + AM1-BCC, nc=0) and parmchk2
  - create systems/molX_poseY_full/ with:
      - ligand mol2/frcmod copied in
      - reference protein/metal patch assets copied in (prepared_no_djm.pdb, set/bonds, frcmod patch)
      - tleap input generated (same as mol6_pose2 protocol) and executed
      - molecule-table repair executed (writes solvated_optB_fixedmol.prmtop)

This reuses the validated 12-site MCPB patching/bonding derived for the full dodecamer.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import time
from datetime import datetime
from contextlib import suppress


ROOT = Path("/path/to/urease_md")

# CPU throttling: keep room for concurrent MD equilibration stages.
# Default is 20 threads unless overridden.
MAX_CPU_THREADS = int(os.environ.get("MAX_CPU_THREADS", "20"))

# Force override (do not use setdefault), because the MD pipeline may set OMP_NUM_THREADS=56.
os.environ["OMP_NUM_THREADS"] = str(MAX_CPU_THREADS)
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

ANTECHAMBER = Path("/usr/local/bin/antechamber")
PARMCHK2 = Path("/usr/local/bin/parmchk2")
TLEAP = Path("/usr/local/bin/tleap")

TEMPLATE_SYS = ROOT / "systems" / "mol6_pose2_full"

XYZ_DIR = ROOT / "params" / "inhibitors" / "structures"

OUT_SYSTEMS_DIR = ROOT / "systems"

RESUME_FROM_MOL = os.environ.get("RESUME_FROM_MOL", "").strip()  # e.g. "mol17"


def is_nonempty(path: Path) -> bool:
    try:
        return path.exists() and path.stat().st_size > 0
    except OSError:
        return False


def mol2_has_DU(mol2_path: Path) -> bool:
    if not is_nonempty(mol2_path):
        return False
    # Antechamber can emit GAFF atom type "DU" when bond/valence perception fails.
    # parmchk2 cannot handle DU, so detect early.
    return " DU" in mol2_path.read_text(errors="ignore")


def write_mdl_with_bonds_from_xyz(xyz_path: Path, mdl_path: Path, *, charge: int = 0) -> None:
    """
    Use RDKit to infer bonds/bond orders from 3D XYZ coordinates and write an MDL molfile.
    This provides explicit connectivity for antechamber and avoids DU atom types in cases
    where PDB-based bond guessing fails.
    """
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import rdDetermineBonds  # type: ignore
    except Exception as e:
        raise RuntimeError(
            "RDKit is required for the DU-recovery path but could not be imported."
        ) from e

    xyz_block = xyz_path.read_text()
    mol = Chem.MolFromXYZBlock(xyz_block)
    if mol is None:
        raise RuntimeError(f"RDKit MolFromXYZBlock failed for {xyz_path}")

    # Infer bonds and bond orders. Hueckel helps for aromatic systems.
    with suppress(TypeError):
        rdDetermineBonds.DetermineBonds(mol, charge=charge, useHueckel=True)  # type: ignore[arg-type]
    if mol.GetNumBonds() == 0:
        rdDetermineBonds.DetermineBonds(mol, charge=charge)  # type: ignore[arg-type]

    mdl_path.write_text(Chem.MolToMolBlock(mol))


@dataclass(frozen=True)
class LigandSpec:
    mol: str  # e.g. "mol6"
    pose: int  # e.g. 2
    xyz: Path

    @property
    def tag(self) -> str:
        return f"{self.mol}_pose{self.pose}"

    @property
    def system_dirname(self) -> str:
        return f"{self.tag}_full"


def run(cmd: list[str], *, cwd: Path | None = None, log: Path | None = None) -> None:
    # Force the CPU caps even if some upstream environment variables were set later.
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(MAX_CPU_THREADS)
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"

    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("w") as fh:
            fh.write(" ".join(cmd) + "\n")
            fh.flush()
            subprocess.run(
                cmd,
                cwd=str(cwd) if cwd else None,
                stdout=fh,
                stderr=subprocess.STDOUT,
                check=True,
                env=env,
            )
    else:
        subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True, env=env)


def parse_xyz(xyz: Path) -> list[tuple[str, float, float, float]]:
    lines = xyz.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"XYZ too short: {xyz}")
    try:
        n = int(lines[0].strip())
    except Exception as e:
        raise ValueError(f"Bad XYZ atom count line in {xyz}: {lines[0]!r}") from e
    atoms: list[tuple[str, float, float, float]] = []
    for ln in lines[2 : 2 + n]:
        parts = ln.split()
        if len(parts) < 4:
            continue
        sym = parts[0]
        x, y, z = map(float, parts[1:4])
        atoms.append((sym, x, y, z))
    if len(atoms) != n:
        raise ValueError(f"Parsed {len(atoms)} atoms but expected {n} in {xyz}")
    return atoms


def xyz_to_pdb(atoms: list[tuple[str, float, float, float]], *, resname: str, chain: str, resid: int) -> str:
    # Minimal PDB with HETATM records. Atom name is element with index.
    out = []
    serial = 1
    for sym, x, y, z in atoms:
        elem = sym.strip().title()
        name = f"{elem}{serial}".ljust(4)[:4]
        out.append(
            f"HETATM{serial:5d} {name:>4s} {resname:>3s} {chain:1s}{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}"
        )
        serial += 1
    out.append("END")
    return "\n".join(out) + "\n"


def discover_ligands(xyz_dir: Path) -> list[LigandSpec]:
    pat = re.compile(r"^(mol\d+)_pose(\d+)_lig\.xyz$")
    specs: list[LigandSpec] = []
    for p in sorted(xyz_dir.glob("mol*_pose*_lig.xyz")):
        m = pat.match(p.name)
        if not m:
            continue
        mol, pose_s = m.group(1), m.group(2)
        specs.append(LigandSpec(mol=mol, pose=int(pose_s), xyz=p))
    return specs


def pick_lowest_pose_per_mol(specs: Iterable[LigandSpec]) -> list[LigandSpec]:
    best: dict[str, LigandSpec] = {}
    for s in specs:
        cur = best.get(s.mol)
        if cur is None or s.pose < cur.pose:
            best[s.mol] = s
    return [best[k] for k in sorted(best.keys(), key=lambda x: int(x.replace("mol", "")))]


def copy_if_missing(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        return
    try:
        if src.resolve() == dst.resolve():
            return
    except FileNotFoundError:
        # Fall back to the copy attempt below to raise a clear error.
        pass
    shutil.copy2(src, dst)


def write_tleap_in(system_dir: Path, lig_mol2: Path, lig_frcmod: Path) -> Path:
    # Build a tleap script identical in logic to systems/mol6_pose2_full/tleap_optB_geom.in
    # but parameterized for this system directory and ligand files.
    set_patch = (TEMPLATE_SYS / "set_patch_optB.txt").read_text()
    bonds = (TEMPLATE_SYS / "bonds_optB.txt").read_text()

    tleap_lines = [
        "source leaprc.protein.ff19SB",
        "source leaprc.gaff2",
        "source leaprc.water.tip3p",
        "addAtomTypes {",
        '        { "M1"  "Ni" "sp3" }',
        '        { "M2"  "Ni" "sp3" }',
        '        { "Y1"  "O" "sp3" }',
        '        { "Y2"  "N" "sp3" }',
        '        { "Y3"  "N" "sp3" }',
        '        { "Y4"  "O" "sp3" }',
        '        { "Y5"  "N" "sp3" }',
        '        { "Y6"  "N" "sp3" }',
        '        { "Y7"  "O" "sp3" }',
        '        { "Y8"  "O" "sp3" }',
        "}",
        "KCX = loadmol2 params/ni_center/mcpb_input/KCX_full_stopgap.mol2",
        "loadamberparams params/ni_center/mcpb_input/KCX.frcmod",
        "loadamberparams params/ni_center/mcpb_input/OHB.frcmod",
        "loadamberparams params/ni_center/mcpb_input/urease_mcpbpy.frcmod",
        f"loadamberparams {lig_frcmod.as_posix()}",
        f"loadamberparams {system_dir.as_posix()}/provisional_autopatch_optB.frcmod",
        f"protein = loadpdb {system_dir.as_posix()}/prepared_no_djm.pdb",
        f"LIG = loadmol2 {lig_mol2.as_posix()}",
        "mol = combine { protein LIG }",
        set_patch.rstrip(),
        bonds.rstrip(),
        "check mol",
        "solvateoct mol TIP3PBOX 12.0",
        "addionsrand mol Na+ 151",
        "addionsrand mol Na+ 140",
        "addionsrand mol Cl- 140",
        f"savepdb mol {system_dir.as_posix()}/solvated_optB.pdb",
        f"saveamberparm mol {system_dir.as_posix()}/solvated_optB.prmtop {system_dir.as_posix()}/solvated_optB.inpcrd",
        "quit",
        "",
    ]
    tleap_in = system_dir / "tleap_optB_geom.in"
    tleap_in.write_text("\n".join(tleap_lines))
    return tleap_in


def repair_molecule_table(system_dir: Path) -> None:
    # Reuse the validated repair script by copying it into the new system dir.
    # The script is now portable (uses __file__ parent).
    script_src = TEMPLATE_SYS / "repair_molecule_table.py"
    script_dst = system_dir / "repair_molecule_table.py"
    try:
        if script_src.resolve() != script_dst.resolve():
            shutil.copy2(script_src, script_dst)
    except FileNotFoundError:
        shutil.copy2(script_src, script_dst)
    run(["python3", str(script_dst)], cwd=ROOT, log=system_dir / "repair_molecule_table.log")


def build_one(spec: LigandSpec) -> None:
    system_dir = OUT_SYSTEMS_DIR / spec.system_dirname
    system_dir.mkdir(parents=True, exist_ok=True)

    # If the system was already fully built, skip it.
    if is_nonempty(system_dir / "solvated_optB_fixedmol.prmtop") and is_nonempty(system_dir / "solvated_optB.inpcrd"):
        print(f"[skip] {spec.system_dirname}: solvated_optB_fixedmol.prmtop already present")
        return

    # Copy protein + patch assets used by the working mol6 system
    copy_if_missing(TEMPLATE_SYS / "prepared_no_djm.pdb", system_dir / "prepared_no_djm.pdb")
    copy_if_missing(TEMPLATE_SYS / "provisional_autopatch_optB.frcmod", system_dir / "provisional_autopatch_optB.frcmod")
    copy_if_missing(TEMPLATE_SYS / "hmr.py", system_dir / "hmr.py")

    # Ligand: XYZ -> PDB for antechamber
    atoms = parse_xyz(spec.xyz)
    lig_pdb = system_dir / f"{spec.tag}_lig.pdb"
    lig_pdb.write_text(xyz_to_pdb(atoms, resname="LIG", chain="A", resid=9001))

    # antechamber + parmchk2 (write into the system dir to keep per-system isolation)
    lig_mol2 = system_dir / f"{spec.tag}_gaff2.mol2"
    lig_frcmod = system_dir / f"{spec.tag}.frcmod"

    if not is_nonempty(lig_mol2) or mol2_has_DU(lig_mol2):
        run(
            [
                str(ANTECHAMBER),
                "-i",
                str(lig_pdb),
                "-fi",
                "pdb",
                "-o",
                str(lig_mol2),
                "-fo",
                "mol2",
                "-c",
                "bcc",
                "-s",
                "2",
                "-nc",
                "0",
                "-at",
                "gaff2",
                "-pf",
                "yes",
            ],
            cwd=ROOT,
            log=system_dir / "antechamber.log",
        )

    if not is_nonempty(lig_frcmod):
        try:
            run(
                [str(PARMCHK2), "-i", str(lig_mol2), "-f", "mol2", "-o", str(lig_frcmod)],
                cwd=ROOT,
                log=system_dir / "parmchk2.log",
            )
        except subprocess.CalledProcessError:
            # If parmchk2 fails due to DU atom types, regenerate connectivity with RDKit
            # (XYZ -> MDL with explicit bonds) and rerun antechamber/parmchk2.
            if mol2_has_DU(lig_mol2):
                mdl = system_dir / f"{spec.tag}_rdkit.mdl"
                write_mdl_with_bonds_from_xyz(spec.xyz, mdl, charge=0)
                run(
                    [
                        str(ANTECHAMBER),
                        "-i",
                        str(mdl),
                        "-fi",
                        "mdl",
                        "-o",
                        str(lig_mol2),
                        "-fo",
                        "mol2",
                        "-c",
                        "bcc",
                        "-s",
                        "2",
                        "-nc",
                        "0",
                        "-at",
                        "gaff2",
                        "-pf",
                        "yes",
                    ],
                    cwd=ROOT,
                    log=system_dir / "antechamber_rdkit_mdl.log",
                )
                run(
                    [str(PARMCHK2), "-i", str(lig_mol2), "-f", "mol2", "-o", str(lig_frcmod)],
                    cwd=ROOT,
                    log=system_dir / "parmchk2_rdkit_mdl.log",
                )
            else:
                raise

    # tleap
    tleap_in = write_tleap_in(system_dir, lig_mol2=lig_mol2, lig_frcmod=lig_frcmod)
    run([str(TLEAP), "-f", str(tleap_in)], cwd=ROOT, log=system_dir / "tleap.log")

    # molecule table repair (needed for stable NPT later)
    repair_molecule_table(system_dir)


def main() -> int:
    specs = discover_ligands(XYZ_DIR)
    if not specs:
        raise SystemExit(f"No ligands found in {XYZ_DIR}")

    selected = pick_lowest_pose_per_mol(specs)
    total = len(selected)
    t0_all = time.time()

    (ROOT / "systems").mkdir(exist_ok=True)
    print(f"Found {len(specs)} xyz files; building {total} systems (lowest pose per mol).")

    started = False if RESUME_FROM_MOL else True
    for idx, s in enumerate(selected, start=1):
        if not started:
            if s.mol == RESUME_FROM_MOL:
                started = True
            else:
                continue
        sys_start = time.time()
        print(f"==> Building {s.system_dirname} ({s.xyz.name}) [{idx}/{total}]")
        build_one(s)
        sys_elapsed = time.time() - sys_start
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(
            f"[throughput] {ts} compound={s.system_dirname} elapsed_sec={sys_elapsed:.1f} "
            f"done={idx}/{total} total_elapsed_sec={time.time() - t0_all:.1f}"
        )

    print("DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

