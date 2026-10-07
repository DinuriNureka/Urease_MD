# mol02 ORCA reoptimization (tight settings)

**Purpose:** Reduce the soft-torsion imaginary frequency (~−18 cm⁻¹) on the sulfonamide ligand (dmj, mol02) by reoptimizing with tighter SCF, optimization, and DFT grid settings.

**Input:** `resp_results/mol02/mol02_reopt.inp`  
**Starting geometry:** `params/inhibitors/orca_opt/mol02/mol02.xyz` (copied to `resp_results/mol02/mol02.xyz` before the run).

**Theory line:** `! wB97X-D4 def2-SVP Opt Freq TightOpt TightSCF DEFGRID3` (same functional/basis as batch ligand jobs). ORCA 6.1 does not accept legacy `Grid5` / `FinalGrid6` on the simple input line; `DEFGRID3` is the supported high-density grid (see ORCA manual, numerical integration).

**Resources:** `%pal nprocs 6`, `%maxcore 3000` (matches other `orca_opt/molNN` inputs).

**Runner:**

```bash
bash scripts/reopt_mol02.sh
```

**Outputs (under `resp_results/mol02/`):** `mol02_reopt.log`, `mol02_reopt.xyz`, etc.

**Post-run:** If the log contains no imaginary modes, `scripts/reopt_mol02.sh` overwrites:

- `params/inhibitors/orca_opt/mol02/mol02.xyz`
- `resp_results/mol02/mol02.xyz`

with the reoptimized coordinates from `mol02_reopt.xyz`. No other molecules are modified.
