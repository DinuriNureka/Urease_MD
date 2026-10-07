# MCPB Continuation with Multi-frame Charges (2026-03-23)

## Goal

Continue from completed multi-frame RESP averaging to the next MCPB stage:
- inject averaged charges into MCPB residue mol2 files
- run MCPB Step 2 (Seminario bonded parameters)
- run MCPB Step 4 (tleap/lib generation)

## Actions Performed

1. Applied multi-frame averaged charges from:
   - `params/ni_center/multiframe_resp/resp_output/multiframe_avg_charges.dat`
2. Updated:
   - `KCX.mol2`
   - `OHB.mol2`
   - `NI1.mol2`
   - `NI2.mol2`
3. Created backups before editing:
   - `*.mol2.pre_multiframe.bak`
4. Restored `NI1.mol2` and `NI2.mol2` to formal +2 charge for MCPB Step 2 compatibility.
5. Set OHB oxygen atom type from `O.3` to `O` in `OHB.mol2` to avoid MCPB parser/type failures.
6. Added explicit water model line in `urease.in`:
   - `water_model     tip3p`
7. Copied ORCA-derived fchk file to MCPB expected names:
   - `urease_small_opt.fchk`
   - `urease_small_fc.fchk`
8. Ran MCPB Step 4 successfully (via runtime patch wrapper; see below), generating:
   - `urease_tleap.in`
   - `urease_mcpbpy.pdb`

## MCPB Runtime Compatibility Issue

This MCPB.py installation has two issues in this workflow:

1. **Ni element recognition with `NI1`/`NI2` residues**
   - PDB parser does not infer `Ni` from residue names `NI1` and `NI2` unless mapping is patched.
2. **Undefined exception symbol**
   - `gene_pre_frcmod_file.py` raises `pymsmtError`, but that symbol is undefined in the current install path.

Because `/usr/local/lib/python3.12/dist-packages` is not writable, Step 2/4 were launched through a Python wrapper that patches these at runtime:
- adds `('NI1','NI')` and `('NI2','NI')` to `METAL_PDB`
- maps `pymsmtError` to `RuntimeError`

## Current Status

- **Step 2:** completed successfully after fixing the ORCA->fchk coordinate parsing.
- **Step 4:** completed; leap input regenerated.
- **Outputs present:**
  - `urease_mcpbpy.frcmod`
  - `urease_mcpbpy_pre.frcmod`
  - `urease_tleap.in`
  - `urease_mcpbpy.pdb`

## Root Cause and Rigorous Fix

- Root cause of previous Step 2 failure: `scripts/orca_hess_to_fchk.py` read coordinates from a non-existent `$coordinates` block in ORCA `.hess`, writing all-zero coordinates to fchk.
- ORCA stores coordinates in `$atoms` (`label mass x y z`, bohr).
- Fix applied: parse coordinates from `$atoms`, regenerate `small_freq.fchk`, and recopy to:
  - `urease_small_opt.fchk`
  - `urease_small_fc.fchk`
- Verification: non-zero coordinates now present in `Current cartesian coordinates` section; MCPB Step 2 produced finite Ni-ligand distances and force constants.

## Remaining Technical Note

This MCPB.py installation still has environment-specific parser bugs for `NI1/NI2` recognition and internal exception naming in the pre-frcmod path. Since `/usr/local` is not writable, MCPB was executed through a runtime wrapper patch during these runs.

This does not change the QM theory; it only patches parser behavior to let MCPB read the existing Ni labels correctly in this environment.

## Post-Step2 Follow-up

- Attempted MCPB Step 3 to generate renamed residue mol2 files (`HD1.mol2`, `HD2.mol2`, `KX1.mol2`, etc.).
- Step 3 failed because MCPB expects Gaussian ESP output:
  - missing `urease_large_mk.log`
  - current workflow uses ORCA output + external multi-frame RESP.
- `tleap` run from `urease_tleap.in` confirms this dependency:
  - fails at first missing renamed mol2 (`HD1.mol2`).

### Implication

For a fully rigorous continuation, we must now provide MCPB-compatible renamed residue mol2 files with:
- custom atom types (`Y1..Y8`, `M1/M2`) from MCPB fingerprints
- multi-frame RESP charges (already available)

Then `urease_tleap.in` can be executed successfully with the bonded model.

## Finalization Performed

Implemented the above and completed tleap model build:

1. Generated renamed residue mol2 files:
   - `HD1.mol2`, `HD2.mol2`, `KX1.mol2`, `HE1.mol2`, `HD3.mol2`, `AP1.mol2`, `N11.mol2`, `N21.mol2`, `OB1.mol2`
2. Source data used for generation:
   - coordinates from `urease_mcpbpy.pdb`
   - atom-type reassignment from `urease_standard.fingerprint`
   - charges from `multiframe_avg_charges.dat`
3. Resolved missing tleap valence terms in `urease_mcpbpy.frcmod`:
   - added `Y4-H` bond (TIP3P OH analog)
   - added `N-C-Y1` and `N-C-Y7` angles (copied from `N-C-O` analog)
4. Ran:
   - `tleap -f urease_tleap.in`
   - Completed with **Errors = 0**

### Final Output Files

- Dry system:
  - `urease_dry.prmtop`
  - `urease_dry.inpcrd`
  - `urease_dry.pdb`
- Solvated system:
  - `urease_solv.prmtop`
  - `urease_solv.inpcrd`
  - `urease_solv.pdb`

