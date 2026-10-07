## Monomer Ni-site patch files (chain B site)

This step derives monomer-specific `tleap` patch snippets for the **chain B** dinuclear Ni site
from the validated dodecamer patch files.

### Source (dodecamer)

- `systems/mol6_pose2_full/bonds_optB.txt` (first 9 lines correspond to chain B site)
- `systems/mol6_pose2_full/set_patch_optB.txt` (first block corresponds to chain B site)

### Change for monomer apo template

In the extracted monomer PDB (`structures/monomer/monomer_BQ_apo.pdb`), the bridging hydroxide is:

- `OHB B9801`

The original dodecamer patch used residue `9709` for the bridging hydroxide, so we replaced:

- `9709` → `9801`

### Outputs

- `params/ni_center/monomer/bonds_monomer.txt`
- `params/ni_center/monomer/set_patch_monomer.txt`

