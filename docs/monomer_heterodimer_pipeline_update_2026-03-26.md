# Monomer Heterodimer Pipeline Update (2026-03-26)

This update switches the simulation unit from the full dodecamer to a single monomer heterodimer (one UreA + one UreB subunit).

## Key configuration changes

`config.env` now sets:

- `SYSTEM_TYPE="monomer_heterodimer"`
- `N_SUBUNITS=1`
- `N_NI_TOTAL=2`
- `N_NI_PER_SITE=2`

## System building

- Added `scripts/build_monomer_system.sh`
  - Uses monomer Ni-center patch/bond snippets from `params/ni_center/monomer/`
  - Applies only one active site (Ni atoms at `mol.1046` and `mol.1047`, OHB at `mol.1048` in the un-solvated tleap unit)
  - Solvates with TIP3P truncated octahedron using the existing `BOX_BUFFER=12.0`

- Updated `scripts/build_system.sh` to delegate to `scripts/build_monomer_system.sh` when `SYSTEM_TYPE=monomer_heterodimer`.

## Equilibration / production topology handling (no deletion of working files)

Updated `scripts/run_equilibration.sh` to support monomer systems that start from:

- `systems/<name>/monomer_<inhib>.prmtop`
- `systems/<name>/monomer_<inhib>.inpcrd`

If an HMR topology is missing, the script generates:

- `systems/<name>/solvated_hmr.prmtop`

without modifying or deleting the original `monomer_<inhib>.prmtop/inpcrd`.

`scripts/run_production.sh` and the downstream analysis scripts were updated to prefer `solvated_hmr.prmtop` over `solvated_hmr.parm7`.

## Analysis mask updates

`scripts/analyze_all.sh` now uses monomer residue numbering when `SYSTEM_TYPE=monomer_heterodimer`:

- Ni–Ni distance: `:570@NI` to `:571@NI` (rebuilt monomer; verify with `cpptraj resinfo :NI`)
- FLAP distance: `:322@CA` to `:570@NI`

These Ni residue numbers were verified from `systems/mol6_pose2_monomer/monomer_mol6.prmtop` via `cpptraj resinfo`.

