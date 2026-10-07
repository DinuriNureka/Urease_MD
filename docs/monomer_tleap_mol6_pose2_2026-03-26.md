## Monomer `tleap` build script (mol6_pose2 test)

This step writes a monomer `tleap` input script for building a solvated monomer system
(chains **B+Q**) using the existing `mol6_pose2` ligand parameters as a test case.

### Files

- `systems/mol6_pose2_monomer/tleap_monomer.in`
- `systems/mol6_pose2_monomer/provisional_autopatch_optB.frcmod` (copied from `systems/mol6_pose2_full/`)

### Notes

- The script uses the monomer-specific Ni-site patch snippets:
  - `params/ni_center/monomer/set_patch_monomer.txt`
  - `params/ni_center/monomer/bonds_monomer.txt`
- It is written for review only (do not run yet).

