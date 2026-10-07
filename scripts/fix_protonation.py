#!/usr/bin/env python3
"""
scripts/fix_protonation.py
---------------------------------------------------------------------------
Fix histidine tautomers and key protonation states for the urease active site.

This script:
  1. Reads the reduce-protonated PDB
  2. Renames metal-coordinating His to HID/HIE per the DFT model
  3. Applies the fix across ALL 12 active sites (all UreB chains)
  4. Writes the final prepared PDB

Usage:
    python3 scripts/fix_protonation.py \
        structures/prepared/6ZJA_H.pdb \
        structures/prepared/6ZJA_prepared.pdb \
        structures/cleaned/6ZJA_clean_renum.txt

Coding standards:
    - Type hints on all functions
    - Logging instead of print()
    - Returns exit code 0/1
    - All configuration at top of file
---------------------------------------------------------------------------
"""

import sys
import logging
from pathlib import Path
from typing import Dict, List, Tuple

# -------------------------------------------------------
# CONFIGURATION — Protonation assignments
# -------------------------------------------------------
# Key: ORIGINAL 6ZJA residue number (before pdb4amber renumbering)
# Value: target AMBER residue name
# These are per-chain assignments — applied to all UreB chains
ORIG_HIS_ASSIGNMENTS: Dict[int, str] = {
    136: "HID",   # His136: Ne2 -> Ni2, so protonate Nd1
    138: "HID",   # His138: Ne2 -> Ni2, so protonate Nd1
    221: "HIE",   # His221: Nd1 -> Ni1, so protonate Ne2
    248: "HIE",   # His248: Nd1 -> Ni1, so protonate Ne2
    274: "HID",   # His274: Ne2 -> Ni1, so protonate Nd1
    322: "HID",   # His322: flap, not metal-bound (DFT model)
    323: "HID",   # His323: flap, not metal-bound (DFT model)
}

# -------------------------------------------------------
# LOGGING
# -------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s"
)
log = logging.getLogger("fix_protonation")


def load_renumbering(renum_file: Path) -> Dict[Tuple[str, int], int]:
    """
    Parse pdb4amber renumbering file.
    Returns mapping: (chain, old_resnum) -> new_resnum.
    """
    mapping: Dict[Tuple[str, int], int] = {}
    with open(renum_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            # Format: OLD_RESNAME CHAIN OLD_RESNUM NEW_RESNAME NEW_RESNUM
            if len(parts) >= 5:
                try:
                    old_chain = parts[1]
                    old_num = int(parts[2])
                    new_num = int(parts[4])
                    mapping[(old_chain, old_num)] = new_num
                except (ValueError, IndexError):
                    continue
    return mapping


def build_new_number_assignments(
    renum: Dict[Tuple[str, int], int],
    orig_assignments: Dict[int, str]
) -> Dict[int, str]:
    """
    Map original residue numbers to renumbered residue numbers.
    Returns: {new_resnum: target_name} for ALL chains.
    """
    new_assignments: Dict[int, str] = {}
    for (chain, old_num), new_num in renum.items():
        if old_num in orig_assignments:
            new_assignments[new_num] = orig_assignments[old_num]
            log.debug(
                f"  Chain {chain} res {old_num} -> {new_num}: "
                f"{orig_assignments[old_num]}"
            )
    return new_assignments


def fix_protonation(
    input_pdb: Path,
    output_pdb: Path,
    assignments: Dict[int, str]
) -> int:
    """
    Read PDB, rename His residues per assignments, write output.
    Returns number of changes made.
    """
    changes = 0
    output_lines: List[str] = []

    with open(input_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                resnum = int(line[22:26])

                if resnum in assignments and resname in ("HIS", "HIP", "HID", "HIE"):
                    target = assignments[resnum]
                    if resname != target:
                        line = line[:17] + f"{target:>3}" + line[20:]
                        changes += 1

            output_lines.append(line)

    with open(output_pdb, "w") as f:
        f.writelines(output_lines)

    return changes


def remove_wrong_hydrogens(pdb_file: Path, assignments: Dict[int, str]) -> int:
    """
    Remove protons that are on the WRONG nitrogen after reduce.
    - HID: must NOT have HE2 (proton on Ne2); must have HD1
    - HIE: must NOT have HD1 (proton on Nd1); must have HE2

    Returns number of atoms removed.
    """
    removed = 0
    kept_lines: List[str] = []

    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resnum = int(line[22:26])
                atom_name = line[12:16].strip()

                if resnum in assignments:
                    target = assignments[resnum]
                    if target == "HID" and atom_name == "HE2":
                        removed += 1
                        continue
                    if target == "HIE" and atom_name == "HD1":
                        removed += 1
                        continue

            kept_lines.append(line)

    with open(pdb_file, "w") as f:
        f.writelines(kept_lines)

    return removed


def main() -> int:
    """Entry point. Returns 0 on success, 1 on error."""
    if len(sys.argv) != 4:
        print("Usage: python3 fix_protonation.py <input.pdb> <output.pdb> <renum.txt>")
        return 1

    input_pdb = Path(sys.argv[1])
    output_pdb = Path(sys.argv[2])
    renum_file = Path(sys.argv[3])

    if renum_file.exists():
        renum = load_renumbering(renum_file)
        log.info(f"Loaded {len(renum)} renumbering entries from {renum_file}")
        assignments = build_new_number_assignments(renum, ORIG_HIS_ASSIGNMENTS)
        log.info(f"Mapped {len(assignments)} His residues across all chains")
    else:
        log.warning(f"No renumbering file found at {renum_file}")
        log.warning("Using original residue numbers (may be incorrect)")
        assignments = ORIG_HIS_ASSIGNMENTS

    if not assignments:
        log.error("No residue assignments to apply. Check renumbering file format.")
        return 1

    changes = fix_protonation(input_pdb, output_pdb, assignments)
    log.info(f"Renamed {changes} His residues in {output_pdb.name}")

    removed = remove_wrong_hydrogens(output_pdb, assignments)
    log.info(f"Removed {removed} incorrect H atoms")

    log.info("--- Verification ---")
    his_counts: Dict[str, int] = {"HID": 0, "HIE": 0, "HIP": 0, "HIS": 0}
    with open(output_pdb) as f:
        seen: set = set()
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                resnum = int(line[22:26])
                key = (resnum, resname)
                if resname in his_counts and key not in seen:
                    his_counts[resname] += 1
                    seen.add(key)

    for name, count in his_counts.items():
        if count > 0:
            log.info(f"  {name}: {count} residues")

    if his_counts.get("HIS", 0) > 0:
        log.warning(
            f"  {his_counts['HIS']} residues still named HIS — "
            "these need manual review"
        )

    log.info(f"[DONE] Output: {output_pdb}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
