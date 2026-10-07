# Batch build: full dodecamer systems for all ligands

This workflow builds `systems/molX_poseY_full/` directories automatically from the ligand-only XYZs in:

- `params/inhibitors/structures/mol*_pose*_lig.xyz`

It follows the **same full-dodecamer MCPB patch + solvation protocol** validated for `mol6_pose2_full`:

- ff19SB + GAFF2 + TIP3P
- 12-site Ni center patching via `set_patch_optB.txt` + `bonds_optB.txt`
- truncated octahedral solvation with 12 Å buffer
- the same ionization commands as the reference build
- writes `solvated_optB.{prmtop,inpcrd,pdb}`
- immediately repairs the topology molecule table, writing `solvated_optB_fixedmol.prmtop`

## What it does per compound

For each `molX`, it selects the **lowest pose index** available and:

1. Converts XYZ → ligand PDB (HETATM records).
2. Runs `antechamber` to generate GAFF2 + AM1-BCC mol2 (assumes ligand charge 0).
3. Runs `parmchk2` to generate frcmod.
4. Generates a per-system `tleap_optB_geom.in` and runs `tleap`.
5. Runs a molecule-table repair (ParmEd) so NPT will not fail later.

Outputs go into:

- `systems/molX_poseY_full/`

## Run it

```bash
python3 scripts/build_full_systems_batch.py
```

## Notes / limitations

- **Charge**: the script assumes each ligand has **net charge 0** (`-nc 0`). If any ligand is charged, we need to extend the script to infer/override charge per mol.
- **MCPB**: this does **not** rerun QM/MCPB for each ligand; it reuses the validated full-dodecamer Ni-site patching derived for the protein (same as `mol6_pose2_full`).
- **Active site placement**: the ligand is loaded as a separate residue (`LIG`) and does not alter the fixed Ni-site patch/bond scheme.

