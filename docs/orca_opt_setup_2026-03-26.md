## ORCA ligand optimization setup (36 jobs)

Prepared ORCA geometry optimization inputs for all 36 ligands and created a parallel launcher script.

### ORCA executable

- From `config.env`: `/opt/software/ORCA/build/orca`
- `which orca` also resolves to: `/opt/software/ORCA/build/orca`

### Setup generated

- Root directory: `params/inhibitors/orca_opt/`
- Per-ligand working dirs: `params/inhibitors/orca_opt/mol01` ... `mol36`
- In each dir:
  - `molNN.xyz` (copied from selected structure pose)
  - `molNN.inp` (template below)

Template used in each `molNN.inp`:

```text
! wB97X-D4 def2-SVP Opt
%maxcore 3000
%pal
  nprocs 6
end
* xyzfile <net_charge> <multiplicity> <molNN>.xyz *
```

`net_charge` and `multiplicity` are from `resp_inputs/job_manifest.csv`.

### Name mapping and pose selection

- Mapping file: `params/inhibitors/orca_opt/ligand_name_mapping.csv`
- Because `structures/` contains multiple poses per ligand name, one structure was selected per manifest entry using:
  - `compound = original_filename` with `_free.xyz` removed
  - selected file = lowest available `pose` index for that compound

Examples:
- `mol32` in manifest (`mol6_free.xyz`) -> `mol6_pose2_lig.xyz`
- `mol35` in manifest (`mol9_free.xyz`) -> `mol9_pose2_lig.xyz`
- `mol05` in manifest (`mol12_free.xyz`) -> `mol12_pose3_lig.xyz`

### Batch runner

- Script: `scripts/run_orca_opt_all_ligands.sh`
- Launches up to 8 jobs in parallel using background jobs and `wait -n` throttling.
- Each job logs to `params/inhibitors/orca_opt/molNN/molNN.log`.
- On success, copies `molNN.gbw` to `resp_results/molNN/resp.gbw`.
- Final summary reports succeeded/failed counts.

### Status

- Inputs and script written.
- ORCA jobs not started yet.
