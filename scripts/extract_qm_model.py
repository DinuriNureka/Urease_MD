#!/usr/bin/env python3
"""
scripts/extract_qm_model.py
---------------------------------------------------------------------------
Extract a small QM model cluster from one active site of the dodecamer.

The model contains:
  - Both Ni2+ ions
  - KCX219 (carbamylated lysine side chain from Cg onward)
  - Bridging OH-
  - His136, His138, His221, His248, His274 -> truncated to methylimidazole
  - Asp362 -> truncated to acetate
  - Capping H at truncation points (added later in Avogadro/GaussView)

This is deliberately smaller than the 257-atom DFT cluster. MCPB.py only
needs the immediate coordination sphere for the Seminario method.

Usage:
    python3 scripts/extract_qm_model.py \
        structures/prepared/6ZJA_prepared.pdb \
        params/ni_center/model_cluster/active_site_model.pdb \
        --chain B --ni1 601 --ni2 602

Inputs:
    Prepared PDB from Phase 1.

Outputs:
    active_site_model.pdb — ~60-80 atom model for MCPB.py
---------------------------------------------------------------------------
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# Atoms to KEEP from each coordinating residue (side-chain ring + Cb)
HIS_KEEP: Set[str] = {
    "CB", "CG", "ND1", "CD2", "CE1", "NE2",
    "HD1", "HD2", "HE1", "HE2", "HB2", "HB3",
}

ASP_KEEP: Set[str] = {
    "CB", "CG", "OD1", "OD2", "HB2", "HB3",
}

KCX_SKIP: Set[str] = {
    "N", "H", "CA", "HA", "C", "O",
}


def parse_pdb(pdb_path: Path) -> List[dict]:
    """Parse PDB into list of atom dicts."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atoms.append({
                "line": line,
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21].strip(),
                "resnum": int(line[22:26]),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "element": line[76:78].strip() if len(line) > 76 else "",
            })
    return atoms


def extract_model(
    atoms: List[dict],
    chain: str,
    ni1_res: int,
    ni2_res: int,
    his_residues: List[int],
    kcx_res: int,
    asp_res: int,
    oh_res: int,
) -> List[dict]:
    """Extract first-shell coordination sphere atoms."""
    target_residues = set(
        [ni1_res, ni2_res, kcx_res, asp_res, oh_res] + his_residues
    )
    model = []

    for atom in atoms:
        if atom["chain"] != chain:
            continue
        if atom["resnum"] not in target_residues:
            continue

        resnum = atom["resnum"]
        name = atom["name"]

        if resnum in his_residues:
            if name not in HIS_KEEP:
                continue
        elif resnum == asp_res:
            if name not in ASP_KEEP:
                continue
        elif resnum == kcx_res:
            if name in KCX_SKIP:
                continue

        model.append(atom)

    return model


def write_pdb(atoms: List[dict], output_path: Path) -> None:
    """Write atom list as PDB with sequential numbering."""
    with open(output_path, "w") as f:
        for i, atom in enumerate(atoms, 1):
            line = atom["line"]
            f.write(f"{line[:6]}{i:5d}{line[11:]}")
        f.write("END\n")


def main() -> int:
    """Entry point. Returns 0 on success, 1 on error."""
    parser = argparse.ArgumentParser(description="Extract QM model cluster")
    parser.add_argument("input_pdb", type=Path)
    parser.add_argument("output_pdb", type=Path)
    parser.add_argument("--chain", default="B")
    parser.add_argument("--ni1", type=int, default=601)
    parser.add_argument("--ni2", type=int, default=602)
    parser.add_argument(
        "--his", type=int, nargs="+", default=[136, 138, 221, 248, 274]
    )
    parser.add_argument("--kcx", type=int, default=219)
    parser.add_argument("--asp", type=int, default=362)
    parser.add_argument("--oh", type=int, default=700)
    args = parser.parse_args()

    log.info(f"Reading {args.input_pdb}")
    atoms = parse_pdb(args.input_pdb)
    log.info(f"  Total atoms in PDB: {len(atoms)}")

    model = extract_model(
        atoms, args.chain,
        args.ni1, args.ni2, args.his, args.kcx, args.asp, args.oh
    )
    log.info(f"  Model cluster: {len(model)} atoms")

    if len(model) < 30:
        log.error(
            f"Model too small ({len(model)} atoms). "
            "Check chain and residue numbers."
        )
        log.error("Did you update residue numbers after pdb4amber renumbering?")
        return 1

    if len(model) > 150:
        log.warning(
            f"Model larger than expected ({len(model)} atoms). "
            "Review atom selection."
        )

    args.output_pdb.parent.mkdir(parents=True, exist_ok=True)
    write_pdb(model, args.output_pdb)
    log.info(f"  Written to {args.output_pdb}")

    resnums = set(a["resnum"] for a in model)
    resnames = set(a["resname"] for a in model)
    ni_count = sum(
        1 for a in model if a["element"] == "NI" or a["name"] == "NI"
    )
    log.info(f"  Residues: {sorted(resnums)}")
    log.info(f"  Residue types: {sorted(resnames)}")
    log.info(f"  Ni atoms: {ni_count}")

    log.info("")
    log.info("[NEXT STEPS]")
    log.info("  1. Open in Avogadro or PyMOL")
    log.info("  2. Add capping H at each truncation point (Cb -> CH3)")
    log.info("  3. Verify Ni coordination geometry is intact")
    log.info("  4. Save final model -> proceed to Step 2.2")
    return 0


if __name__ == "__main__":
    sys.exit(main())
