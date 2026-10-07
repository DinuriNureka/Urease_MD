# ORCA Opt + Frequency Update (36 ligands)

This step updated the ORCA keyword line used for ligand geometry optimizations to include a frequency calculation.

## Change

Replaced in:
- `params/inhibitors/orca_opt_template.inp`
- all ligand optimization inputs under `params/inhibitors/orca_opt/molNN/*.inp`

The ORCA method line:

- from: `! wB97X-D4 def2-SVP Opt`
- to:   `! wB97X-D4 def2-SVP Opt Freq`

## Verification

Checked that **all 72** generated input files (36 `molNN.inp` plus 36 `input.inp`) now have the `Opt Freq` keyword line on the first line.
