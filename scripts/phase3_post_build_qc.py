#!/usr/bin/env python3
"""
Post–Phase 3 QC before GPU Phase 4. Reads built monomer systems under systems/<compound>_pose*_monomer/
and validates solvation, charges, HMR, geometry (Ni–Ni), ligand placement (centroid vs protein CA centroid,
protein–ligand contacts), tleap logs, and consistency with run_full_pipeline.sh MD templates.

Ligand placement (critical):
  - Centroid distance (ligand all-atom centroid vs protein CA centroid) must be ≤ 30 Å.
  - At least 5 distinct protein residues must contact the ligand (any ligand atom ↔ any protein atom
    within 5 Å via cpptraj nativecontacts); protein = not water/ions/ligand/NI.

Informational (not pass/fail):
  - Minimum distance from any ligand atom to any Ni (metal site sanity).

Exit 0 if all critical checks pass; exit 1 otherwise (stderr has details).
"""
from __future__ import annotations

import os
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", "/path/to/urease_md"))
MAPPING_CSV = PROJECT_ROOT / "params/inhibitors/orca_opt/ligand_name_mapping.csv"
PIPELINE_SH = PROJECT_ROOT / "scripts/run_full_pipeline.sh"
_ah = Path(os.environ.get("AMBER_HOME", "/opt/software/AMBER/pmemd24"))
CPPTRAJ = (
    os.environ.get("CPPTRAJ")
    or os.environ.get("CPPTRAJ_BIN")
    or shutil.which("cpptraj")
    or str(_ah / "bin/cpptraj")
)
BOX_BUFFER_EXPECT = float(os.environ.get("BOX_BUFFER", "12.0"))

WATER_RES = {"WAT", "HOH", "TIP3", "TIP4", "SOL"}
ION_RES = {"NA+", "CL-", "NA", "CL", "K+", "BR-", "I-", "CS+", "MG2", "CA2", "ZN2", "CU2", "CU1"}
STANDARD_AA = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "ASH",
    "GLH",
    "LYN",
    "CYX",
    "TYM",
}


def load_mapping_rows() -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    import csv

    with MAPPING_CSV.open(newline="") as f:
        for r in csv.DictReader(f):
            c = r["compound"].strip()
            p = str(int(r["selected_pose"]))
            rows.append((c, p))
    return rows


def try_parmed():
    try:
        import parmed as pmd  # noqa: F401

        return True
    except ImportError:
        return False


def infer_ligand_resname(p) -> str:
    """Return residue name used for the small-molecule ligand (LIG, MOL, etc.)."""
    import parmed

    assert isinstance(p, parmed.Structure)
    counts: dict[str, int] = {}
    for r in p.residues:
        n = r.name.strip()
        counts[n] = counts.get(n, 0) + 1
    for candidate in ("LIG", "MOL"):
        if candidate in counts:
            return candidate
    for r in p.residues:
        n = r.name.strip()
        if n in WATER_RES or n in ION_RES:
            continue
        if n in ("NI", "KCX", "OHB"):
            continue
        if n not in STANDARD_AA:
            return n
    return "LIG"


def prmtop_stats(p, prmtop_path: Path | None = None) -> dict:
    """Charge and atom counts from a loaded parmed Structure (topology)."""
    q_total = float(sum(p.parm_data["CHARGE"]))
    natoms = len(p.atoms)
    ion_res = {"Na+", "Cl-", "NA", "CL"}
    q_solute = 0.0
    n_na = n_cl = 0
    for i, atom in enumerate(p.atoms):
        res = atom.residue.name.strip()
        if res in ("WAT", "HOH", "TIP3", "TIP4"):
            continue
        if res in ion_res:
            if "Na" in res or res == "NA":
                n_na += 1
            elif "Cl" in res or res == "CL":
                n_cl += 1
            continue
        q_solute += float(p.parm_data["CHARGE"][i])
    size = prmtop_path.stat().st_size if prmtop_path is not None and prmtop_path.is_file() else 0
    return {
        "natoms": natoms,
        "q_total": q_total,
        "q_solute_before_ions": q_solute,
        "n_na": n_na,
        "n_cl": n_cl,
        "size": size,
    }


def hmr_sample_mean_H_mass(p) -> float:
    """Non-water hydrogens after HMR (atomic_number==1); expect mean ~3.024 Da."""
    masses = []
    for atom in p.atoms:
        if atom.residue.name.strip() in ("WAT", "HOH", "TIP3", "TIP4"):
            continue
        if getattr(atom, "atomic_number", None) != 1:
            continue
        masses.append(float(atom.mass))
    if not masses:
        return float("nan")
    return sum(masses) / len(masses)


def coords_array(p) -> np.ndarray:
    return np.array([[a.xx, a.xy, a.xz] for a in p.atoms], dtype=np.float64)


def ni_ni_distance(p) -> float | None:
    """Shortest Ni–Ni distance; prefer dinuclear pair in ~3.0–4.5 Å if present."""
    idx = [i for i, a in enumerate(p.atoms) if getattr(a, "atomic_number", None) == 28]
    if len(idx) < 2:
        return None
    crd = coords_array(p)
    in_range: list[float] = []
    best = None
    for i in range(len(idx)):
        for j in range(i + 1, len(idx)):
            d = float(np.linalg.norm(crd[idx[i]] - crd[idx[j]]))
            if 3.0 <= d <= 4.5:
                in_range.append(d)
            if best is None or d < best:
                best = d
    if len(in_range) == 1:
        return in_range[0]
    if len(in_range) > 1:
        return min(in_range, key=lambda x: abs(x - 3.5))
    return best


def ligand_min_distance_to_ni(p, lig_resname: str) -> float | None:
    """Minimum distance from any ligand atom to any Ni (informational)."""
    crd = coords_array(p)
    ni_idx = [i for i, a in enumerate(p.atoms) if getattr(a, "atomic_number", None) == 28]
    lig_idx = [i for i, a in enumerate(p.atoms) if a.residue.name.strip() == lig_resname]
    if not ni_idx or not lig_idx:
        return None
    lig_c = crd[lig_idx]
    ni_c = crd[ni_idx]
    dmin = np.min(np.linalg.norm(lig_c[:, None, :] - ni_c[None, :, :], axis=2))
    return float(dmin)


def cpptraj_resnums_for_ligand(p, lig_resname: str) -> set[int]:
    """cpptraj residue numbers are parmed Residue.idx + 1 for this project."""
    return {r.idx + 1 for r in p.residues if r.name.strip() == lig_resname}


def resnum_to_resname_map(p) -> dict[int, str]:
    return {r.idx + 1: r.name.strip() for r in p.residues}


def is_protein_contact_residue(resname: str, lig_resname: str) -> bool:
    n = resname.strip().upper()
    if n == lig_resname.strip().upper():
        return False
    if n in {x.upper() for x in WATER_RES}:
        return False
    if n in {x.upper() for x in ION_RES}:
        return False
    if n == "NI":
        return False
    return True


def centroid_ligand_vs_protein_ca(p, lig_resname: str) -> float | None:
    """Distance between ligand (all atoms) centroid and protein CA centroid."""
    crd = coords_array(p)
    lig_mask = np.array([a.residue.name.strip() == lig_resname for a in p.atoms])
    ca_mask = np.array(
        [
            (a.name.strip() == "CA")
            and is_protein_contact_residue(a.residue.name.strip(), lig_resname)
            for a in p.atoms
        ]
    )
    if not np.any(lig_mask) or not np.any(ca_mask):
        return None
    lig_c = crd[lig_mask].mean(axis=0)
    prot_c = crd[ca_mask].mean(axis=0)
    return float(np.linalg.norm(lig_c - prot_c))


_CONTACT_PAIR_RE = re.compile(r":(\d+)@(\S+)_:(\d+)@(\S+)")


def count_protein_contacts_nativecontacts(prmtop: Path, crd: Path, p, lig_resname: str) -> tuple[int | None, str | None]:
    """
    cpptraj nativecontacts; return count of unique protein residue numbers within 5 Å of ligand.
    """
    if not Path(CPPTRAJ).is_file():
        return None, "cpptraj binary missing"
    import tempfile

    prmtop = Path(prmtop).resolve()
    crd = Path(crd).resolve()
    lig_nums = cpptraj_resnums_for_ligand(p, lig_resname)
    if not lig_nums:
        return None, "no ligand residues"
    rmap = resnum_to_resname_map(p)
    tmp = tempfile.mkdtemp(prefix="qc_nc_")
    try:
        out = os.path.join(tmp, "nc.dat")
        cpp = os.path.join(tmp, "run.cpptraj")
        mask_lig = f":{lig_resname}"
        mask_prot = f"!(:WAT,HOH,TIP3,TIP4,Na+,Cl-,NA,CL,:{lig_resname})"
        Path(cpp).write_text(
            f"parm {prmtop}\ntrajin {crd}\n"
            f"nativecontacts {mask_lig} {mask_prot} distance 5.0 writecontacts {out}\nrun\n"
        )
        proc = subprocess.run([CPPTRAJ, "-i", cpp], capture_output=True, text=True, timeout=300)
        if proc.returncode != 0:
            return None, (proc.stderr or proc.stdout or "cpptraj failed")[:500]
        if not Path(out).is_file():
            return None, "writecontacts file missing"
        protein_res: set[int] = set()
        for ln in Path(out).read_text(errors="ignore").splitlines():
            if ln.strip().startswith("#") or not ln.strip():
                continue
            m = _CONTACT_PAIR_RE.search(ln)
            if not m:
                continue
            a, b = int(m.group(1)), int(m.group(3))
            if a in lig_nums and b not in lig_nums:
                pr = b
            elif b in lig_nums and a not in lig_nums:
                pr = a
            else:
                continue
            name = rmap.get(pr, "")
            if is_protein_contact_residue(name, lig_resname):
                protein_res.add(pr)
        return len(protein_res), None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def scan_tleap_warnings() -> list[str]:
    """Flag likely problems only (KCX connect0 / auto-patch remarks are common and usually benign)."""
    issues: list[str] = []
    critical = re.compile(
        r"close contact|missing parameter|Could not find (bond|angle|dihedral)|"
        r"unusual bond|FATAL|fatal error|Error\!",
        re.I,
    )
    for log in sorted(PROJECT_ROOT.glob("systems/*_monomer/tleap_pipeline_*.log")):
        text = log.read_text(errors="ignore")
        hit: list[str] = []
        for line in text.splitlines():
            if critical.search(line) and "KCX" not in line and "connect0" not in line.lower():
                hit.append(line.strip()[:220])
        if hit:
            issues.append(f"{log.relative_to(PROJECT_ROOT)}: " + "; ".join(hit[:4]))
    return issues


def verify_pipeline_templates() -> list[str]:
    """Check run_full_pipeline.sh for ntc/ntf/dt/ig in production-style blocks."""
    bad: list[str] = []
    text = PIPELINE_SH.read_text(errors="ignore")
    if "ntc=2" not in text or "ntf=2" not in text:
        bad.append("run_full_pipeline.sh: expected ntc=2 and ntf=2")
    if "dt=0.004" not in text:
        bad.append("run_full_pipeline.sh: expected dt=0.004 in MD sections")
    if "ig=-1" not in text:
        bad.append("run_full_pipeline.sh: expected ig=-1 in production-style inputs")
    return bad


def main() -> int:
    if not try_parmed():
        print("CRITICAL: parmed not importable", file=sys.stderr)
        return 1

    import parmed as pmd

    rows = load_mapping_rows()
    if len(rows) != 36:
        print(f"CRITICAL: expected 36 mapping rows, got {len(rows)}", file=sys.stderr)
        return 1

    built: list[tuple[str, dict]] = []
    missing: list[str] = []

    for compound, pose in rows:
        name = f"{compound}_pose{pose}_monomer"
        sdir = PROJECT_ROOT / "systems" / name
        ptop = sdir / "solvated_hmr.prmtop"
        crd = sdir / "solvated_hmr.rst7"
        if not ptop.is_file() or not crd.is_file():
            missing.append(name)
            continue
        p = pmd.load_file(str(ptop), str(crd))
        st = prmtop_stats(p, prmtop_path=ptop)
        st["name"] = name
        st["compound"] = compound
        st["lig_resname"] = infer_ligand_resname(p)
        st["hmr_mean_H"] = hmr_sample_mean_H_mass(p)
        st["ni_ni"] = ni_ni_distance(p)
        st["lig_min_d_ni"] = ligand_min_distance_to_ni(p, st["lig_resname"])
        st["centroid_ca_dist"] = centroid_ligand_vs_protein_ca(p, st["lig_resname"])
        nprot, err = count_protein_contacts_nativecontacts(ptop, crd, p, st["lig_resname"])
        st["n_prot_res_contacts_5A"] = nprot
        st["contacts_err"] = err
        tin = sorted(sdir.glob("tleap_pipeline_*.in"))
        st["solvateoct_ok"] = False
        if tin:
            ttxt = tin[-1].read_text()
            st["solvateoct_ok"] = (
                "solvateoct" in ttxt and "TIP3PBOX" in ttxt and str(BOX_BUFFER_EXPECT) in ttxt.replace(" ", "")
            )
        built.append((name, st))

    print("=== Phase 3 post-build QC ===", flush=True)
    print(f"Systems with solvated_hmr.prmtop + .rst7: {len(built)} / 36", flush=True)
    if missing:
        print("Missing solvated build:", ", ".join(missing), flush=True)

    if len(built) < 36:
        print("CRITICAL: not all 36 systems built.", file=sys.stderr, flush=True)
        return 1

    natoms_list = [b[1]["natoms"] for b in built]
    natoms_list.sort()
    med = natoms_list[len(natoms_list) // 2]
    low, hi = med * 0.95, med * 1.05
    outliers = [b[0] for b in built if not (low <= b[1]["natoms"] <= hi)]

    print(f"Median atom count (solvated HMR): {med}; 5% band [{low:.0f}, {hi:.0f}]")
    if outliers:
        print("CRITICAL: atom-count outliers (>5% from median):", outliers, file=sys.stderr)

    print(
        "\nPer-system summary (prmtop bytes | natoms | lig | q_total | q_solute* | mean H | "
        "Ni–Ni Å | min(lig–Ni) Å | CA-centroid dist Å | prot res ≤5Å | solvateoct)"
    )
    crit = []
    for name, st in sorted(built, key=lambda x: x[0]):
        qtot = st["q_total"]
        if abs(qtot) > 0.05:
            crit.append(f"{name}: total charge {qtot} not ~0")
        if abs(st["q_solute_before_ions"] - round(st["q_solute_before_ions"])) > 0.2:
            crit.append(f"{name}: solute charge {st['q_solute_before_ions']:.4f} not near integer")
        ni = st["ni_ni"]
        if ni is not None and not (3.0 <= ni <= 4.5):
            crit.append(f"{name}: Ni–Ni {ni:.3f} Å outside 3.0–4.5")
        cd = st["centroid_ca_dist"]
        if cd is None:
            crit.append(f"{name}: centroid ligand vs protein CA could not be computed")
        elif cd > 30.0:
            crit.append(f"{name}: ligand vs protein CA centroid {cd:.2f} Å > 30 (outside protein)")
        nc = st["n_prot_res_contacts_5A"]
        if nc is None:
            crit.append(f"{name}: nativecontacts failed ({st.get('contacts_err') or 'unknown'})")
        elif nc < 5:
            crit.append(f"{name}: only {nc} protein residues contact ligand within 5 Å (need ≥5)")
        hm = st["hmr_mean_H"]
        if hm == hm and not (2.9 <= hm <= 3.15):
            crit.append(f"{name}: mean H mass {hm:.4f} not ~3.024 (HMR)")
        if not st.get("solvateoct_ok"):
            crit.append(f"{name}: tleap .in missing solvateoct TIP3PBOX {BOX_BUFFER_EXPECT}")

        ld = st["lig_min_d_ni"]
        ld_s = f"{ld:.3f}" if ld is not None else "NA"
        cd_s = f"{cd:.3f}" if cd is not None else "NA"
        nc_s = str(nc) if nc is not None else "NA"
        print(
            f"  {name}: {st['size']} | {st['natoms']} | {st['lig_resname']} | {qtot:.4f} | "
            f"{st['q_solute_before_ions']:.4f} | {st['hmr_mean_H']:.4f} | "
            f"{ni if ni is not None else 'NA'} | {ld_s} | {cd_s} | {nc_s} | {st['solvateoct_ok']}"
        )

    tw = scan_tleap_warnings()
    if tw:
        print("\nTleap log hits (warning / close contact / missing):")
        for t in tw[:40]:
            print(" ", t)
        if len(tw) > 40:
            print(f"  ... and {len(tw) - 40} more")

    pt = verify_pipeline_templates()
    if pt:
        print("CRITICAL pipeline template issues:", pt, file=sys.stderr)
        crit.extend(pt)

    if outliers:
        crit.append("atom-count outliers")

    if crit:
        print("\nCRITICAL failures:", file=sys.stderr)
        for c in crit:
            print(" ", c, file=sys.stderr)
        return 1

    rng = random.Random(42)
    sample = rng.sample([b[0] for b in built], k=min(3, len(built)))
    print("\nRandom sample (seed=42) — Ni–Ni, centroid, contacts, LIG–Ni min:")
    for sname in sample:
        st = next(x[1] for x in built if x[0] == sname)
        print(
            f"  {sname}: Ni–Ni={st['ni_ni']}, centroid={st['centroid_ca_dist']}, "
            f"n_contacts={st['n_prot_res_contacts_5A']}, lig–Ni_min={st['lig_min_d_ni']}"
        )

    print("\nPASS: Phase 3 post-build QC (critical checks).")
    print("Note: min(lig–Ni) is informational only (not used for pass/fail).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
