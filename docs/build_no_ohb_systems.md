## Build all 36 no-OHB solvated monomer systems

### What this does
Creates **36** solvated monomer systems (chains **B+Q**) in:

- `systems/<compound>_pose<pose>_monomer/`

using the **no-OHB** Ni-center MCPB outputs:

- `params/ni_center_no_ohb/mcpb_input/urease_no_ohb_mcpbpy.frcmod`
- `params/ni_center_no_ohb/mcpb_input/KCX.frcmod`
- `params/ni_center_no_ohb/mcpb_input/KCX_full_stopgap.mol2`
- `params/ni_center_no_ohb/monomer/bonds_monomer_no_ohb.txt`
- `params/ni_center_no_ohb/monomer/set_patch_monomer_no_ohb.txt`

and **does not load OHB** anywhere.

### Script
- `scripts/build_all_no_ohb_systems.sh`

### Input ligand list
- `params/inhibitors/orca_opt/ligand_name_mapping.csv`

### Key build settings
- **Force fields**: `ff19SB` (`leaprc.protein.ff19SB`), `gaff2`, `tip3p`
- **Solvation**: `solvateoct mol TIP3PBOX 12.0`
- **Neutralization**: `addionsrand mol Na+ 0` then `addionsrand mol Cl- 0`
- **Fix atom ordering**: runs `scripts/fix_atomorder_amber.sh` (cpptraj `fixatomorder`)
- **HMR**: via `parmed` → outputs `solvated_hmr.prmtop` and `solvated_hmr.rst7`

### ORCA pose gate (pre-build)
For each row, the build only proceeds if the selected ORCA pose is accepted as a minimum:

- expects `params/inhibitors/orca_opt/<mol_id>/<compound>_pose<pose>_lig.xyz`
- requires `params/inhibitors/orca_opt/<mol_id>/<mol_id>.log` contains `ORCA TERMINATED NORMALLY`
- rejects logs containing `***imaginary mode***`

### Logs
- Per system: `logs/build_no_ohb_<sys>.log`
- Summary table: `logs/build_no_ohb_summary.txt`

### Verification (post-build)
For each system, checks:
1. `solvated_hmr.prmtop` exists and is **>10 MB**
2. `solvated_hmr.rst7` exists and is **>1 MB**
3. Ni–Ni distance in `solvated_hmr.rst7` is **3.0–4.5 Å** (cpptraj `distance`)
4. No `OHB` residue is present (`cpptraj resinfo :OHB` should select **0 residues**)

### Run (after review)

```bash
cd /path/to/urease_md
bash scripts/build_all_no_ohb_systems.sh
```

