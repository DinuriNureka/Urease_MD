# Pipeline ligand mol2 (`params/inhibitors/mol2/*_gaff2.mol2`)

## Build script

`scripts/build_pipeline_ligand_mol2.py` reads `params/inhibitors/orca_opt/ligand_name_mapping.csv` and writes one file per compound: `params/inhibitors/mol2/<compound>_gaff2.mol2`.

- **RESP charges** come from `gaff2_params/<mol_id>/<mol_id>.mol2` (Multiwfn RESP workflow).
- **Coordinates** prefer the docked pose: `systems/<compound>_pose<pose>_full/<compound>_pose<pose>_gaff2.mol2`. Charges are mapped onto docked atoms via **graph isomorphism** (element + connectivity).

### Seeding missing `systems/..._gaff2.mol2` from DFT cluster ligands

If a docked mol2 is absent but you have a ligand-only `.xyz` (e.g. `params/inhibitors/structures/<compound>_pose<pose>_lig.xyz`, extracted from `/path/to/data/dft_batch/dft_jobs/<compound>_pose<pose>/<compound>_pose<pose>.xyz`), run:

```bash
python3 scripts/seed_docked_system_mol2_from_lig_xyz.py --mol-id mol01 --compound aha --pose 1 \
  --lig-xyz params/inhibitors/structures/aha_pose1_lig.xyz
```

The script keeps GAFF atom order and bond topology from `gaff2_params/<mol_id>/<mol_id>.mol2` and assigns coordinates via direct line order when element sequences match, or Kabsch + Hungarian assignment when they do not.

## Fallbacks (see `params/inhibitors/mol2/PIPELINE_LIGAND_BUILD_WARNINGS.txt`)

1. **Gaff vs docked graph mismatch** (mol21 / mol15): **geometry-mapped RESP** — same element multiset and atom count; charges from gas-phase `mol15.mol2` are assigned to docked atoms by minimizing distance after Kabsch alignment (thresholds in `geometry_charge_map`). This is weaker than isomorphism-based mapping; consider re-running RESP on the docked `mol21_pose2_gaff2.mol2` geometry for publication-quality charges.
2. **No docked mol2** and no manual seed: copy `gaff2_params` mol2 (ORCA-isolated geometry).

The same script copies GAFF2 parameter files to `params/inhibitors/frcmod/<compound>.frcmod` from `gaff2_params/<mol_id>/<mol_id>.frcmod` (names aligned with `ligand_name_mapping.csv`).

## Rebuild after RESP updates

```bash
python3 scripts/build_pipeline_ligand_mol2.py
```

## AMBER

`config.env` sets `AMBER_HOME=/opt/software/AMBER/pmemd24` so `${AMBER_HOME}/bin/pmemd.cuda` resolves correctly.
