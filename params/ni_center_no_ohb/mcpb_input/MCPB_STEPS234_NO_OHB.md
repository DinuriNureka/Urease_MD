# MCPB.py steps 2–4 (no-OHB Ni center) — 2026-03-30

## Prerequisites satisfied

- `AMBERHOME` / `AMBER_HOME` must point at the AMBER/AmberTools tree that contains `dat/leap/parm/` (here: `/usr/local`, not `/usr/local/amber`).
- `MCPB.py` alone fails on this system for `NI1`/`NI2` residues (`METAL_PDB` lacks those names) and a missing `pymsmtError` import in `gene_pre_frcmod_file.py`. Runs used the project wrapper: `scripts/run_mcpb_patched.py`.
- Step 2 (Seminario) expects a Gaussian-format `.fchk` with Hessian. Generated from ORCA:  
  `scripts/orca_hess_to_fchk.py ../orca_calcs/small_freq/small_freq.hess urease_no_ohb_small_opt.fchk`
- Step 3 (`resp`) needs an ESP grid in Gaussian-log form. Built from ORCA:  
  `params/ni_center/mcpb_input/orca_vpot_to_gaussian_log.py ../orca_calcs/large_mk/large_mk.vpot urease_no_ohb_large_mk.log`

## Commands run (from this directory)

```bash
export AMBERHOME=/usr/local
export AMBER_HOME=/usr/local
python3 /path/to/urease_md/scripts/run_mcpb_patched.py -i urease_no_ohb.in -s 2 2>&1 | tee mcpb_step2.log
python3 /path/to/urease_md/scripts/run_mcpb_patched.py -i urease_no_ohb.in -s 3 2>&1 | tee mcpb_step3.log
python3 /path/to/urease_md/scripts/run_mcpb_patched.py -i urease_no_ohb.in -s 4 2>&1 | tee mcpb_step4.log
```

Exit codes: all **0**.

## Outputs

| File | Role |
|------|------|
| `urease_no_ohb_mcpbpy.frcmod` | Final bonded + NB parameters |
| `urease_no_ohb_mcpbpy_pre.frcmod` | Pre-Seminario frcmod |
| `HD1.mol2` … `N21.mol2` | Renamed metal-site residues after RESP |
| `urease_no_ohb_tleap.in` | tLeap script from step 4 |
| `mcpb_step2.log` … `mcpb_step4.log` | Full stdout/stderr per step |

## Note on “Y4”

`Y4` appears in the new frcmod as the **custom type for a histidine NE2** (HID residue 1) coordinated to Ni2 (`M2`). It is **not** the removed bridging hydroxide residue; the no-OHB model has no `OB1` mol2.

## QC warning

Seminario output lists some **non-physical bond lengths** (e.g. ~3–11 Å) and several **negative Hessian-derived** terms zeroed by MCPB. If those persist after checking the `.fchk` and small-model PDB, consider re-optimizing the small model or tightening the QM level.
