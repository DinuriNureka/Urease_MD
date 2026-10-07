## ORCA ligand run: `mol2` pose 5 (`mol24_pose5`)

### Goal
Run ORCA **geometry optimization + frequency** for **compound `mol2` pose 5** using the same settings as the latest ligand ORCA runs.

### Inputs
- **Structure**: `params/inhibitors/structures/mol2_pose5_lig.xyz`
- **Workdir**: `params/inhibitors/orca_opt/mol24_pose5/`

### ORCA settings (match latest ligand runs)
- **Method**: `wB97X-D4`
- **Basis**: `def2-SVP`
- **Job**: `Opt Freq`
- **Parallel**: `%pal nprocs 28 end`
- **Memory**: `%maxcore 2000` (MB per core)
- **ORCA binary**: `/opt/software/ORCA/build/orca`

### Files created
- `params/inhibitors/orca_opt/mol24_pose5/mol24_pose5.xyz` (copied from pose5 structure)
- `params/inhibitors/orca_opt/mol24_pose5/mol24_pose5.inp` (derived from `params/inhibitors/orca_opt/mol24/mol24.inp`)
- `params/inhibitors/orca_opt/mol24_pose5/mol24_pose5.log` (runtime output)

### Launch command
Run from project root:

```bash
cd /path/to/urease_md
ORCA=/opt/software/ORCA/build/orca

cd params/inhibitors/orca_opt/mol24_pose5
nohup env OMPI_MCA_rmaps_base_oversubscribe=1 \
          OMPI_MCA_hwloc_base_binding_policy=none \
          ${ORCA} mol24_pose5.inp > mol24_pose5.log 2>&1 &
echo "mol24_pose5 PID: $!"
```

### Completion criteria
Confirm `params/inhibitors/orca_opt/mol24_pose5/mol24_pose5.log` contains:
- `***********************HURRAY********************`
- `****ORCA TERMINATED NORMALLY****`

