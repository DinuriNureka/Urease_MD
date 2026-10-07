# RESP Pipeline for 36 Urease Inhibitors

## Overview
This pipeline converts GFN2-xTB optimized free-ligand XYZ files into AMBER-ready GAFF2 parameters using ORCA single-point ESP and antechamber RESP charges. The workflow is: ligand XYZ -> ORCA `wB97X-D4/def2-SVP` single-point (`def2/J`, `TightSCF`) -> RESP charge fitting in antechamber -> `mol2` and `frcmod` outputs for tLeap. Single-point is used (no ORCA re-optimization) because the ligand geometries are already GFN2-xTB optimized and preserving them avoids potential gas-phase over-relaxation away from the docking-consistent conformer.

## Prerequisites
- ORCA 6.1.1: ensure `/opt/software/ORCA/build/orca` is executable.
- AmberTools 23+: provides `antechamber` and `parmchk2`.
- Python 3.10+: required for input prep and RMSD checks.
- RDKit (recommended): charge sanity checks in `prepare_resp_inputs.py` (`pip install rdkit` via conda/rdkit channel is typical).
- OpenBabel (optional): useful for manual charge sanity checks (`sudo apt install openbabel` or conda package).

## Directory layout

Before:

```text
/path/to/data/free_ligands/   # input XYZ files (read-only)
prepare_resp_inputs.py
run_resp_batch.sh
check_resp_rmsd.py
run_antechamber_resp.sh
RESP_PIPELINE.md
```

After:

```text
/path/to/data/free_ligands/   # input XYZ files (read-only)
resp_inputs/
  job_manifest.csv
  prepare.log
  mol01.xyz ... mol36.xyz
resp_results/
  resp_batch.log
  rmsd_check.csv
  mol01/ ... mol36/
    resp.inp
    resp.out
    DONE or FAILED
gaff2_params/
  antechamber.log
  mol01/ ... mol36/
    molNN.mol2
    molNN.frcmod
```

## Usage
Run these commands in order from project root:

```bash
python3 prepare_resp_inputs.py --xyz_dir /path/to/data/free_ligands
bash run_resp_batch.sh
python3 check_resp_rmsd.py
bash run_antechamber_resp.sh
```

## Interpreting RMSD results
Because ORCA is run as a single-point calculation (no `Opt`), geometry should remain effectively unchanged and RMSD should be near zero (typically floating-point noise level). `PASS` means expected geometry integrity. `WARN` indicates a small but unusual shift that should be spot-checked. `FAIL` indicates likely file mismatch/corruption or wrong coordinate block parsing and should be fixed before antechamber. `ORCA_FAILED` means SCF did not complete; resolve ORCA job first.

## Troubleshooting
- SCF convergence failure (`FAILED` in batch): add `SlowConv` or `KDIIS` to `resp.inp` and rerun with `--retry-failed`.
- `def2/J` not found: check ORCA auxiliary basis setup via `ORCA_AUXBASIS`, or use `AutoAux` as fallback.
- `antechamber` says `orcaout not recognized`: verify AmberTools version (`antechamber -h`); older versions may require an alternate RESP route.
- `WARN_CHARGE` in `prepare_resp_inputs.py`: manually inspect with OpenBabel, e.g. `obabel input.xyz -omol2 -O check.mol2`, then verify chemistry/valence.

## Output files reference

| Filename | Location | Used by | Description |
|---|---|---|---|
| `job_manifest.csv` | `resp_inputs/` | scripts 2, 3, 4 | Molecule index, file mapping, charge, multiplicity, validation status |
| `prepare.log` | `resp_inputs/` | user QC | Validation and RDKit charge-sanity log |
| `resp.inp` | `resp_results/molNN/` | ORCA | Single-point ESP input (`wB97X-D4/def2-SVP`, `def2/J`, `TightSCF`) |
| `resp.out` | `resp_results/molNN/` | scripts 3, 4 | ORCA output used for SCF status, geometry, and RESP charge fitting |
| `DONE` / `FAILED` | `resp_results/molNN/` | script 2 | Resume markers for crash-safe reruns |
| `resp_batch.log` | `resp_results/` | user QC | Timestamped ORCA batch execution log |
| `rmsd_check.csv` | `resp_results/` | script 4 | Geometry/energy/gap QC table and gating for parameterization |
| `molNN.mol2` | `gaff2_params/molNN/` | tLeap | GAFF2 atom types + CHELPG (ESP-fit) charges from ORCA |
| `resp.vpot` | `resp_results/molNN/` | `calc_esp_rmse.py` | ORCA 6 ESP grid (binary-style text; see script) |
| `resp.pc_chelpg` | `resp_results/molNN/` | `apply_pc_chelpg_to_mol2.py` | Point charges for CHELPG (Merz–Kollman-type) |
| `molNN.frcmod` | `gaff2_params/molNN/` | tLeap | Missing bonded parameter supplements for GAFF2 |
| `antechamber.log` | `gaff2_params/` | user QC | antechamber/parmchk2 full log |

### ORCA 6.1.1 input actually used (`run_resp_batch.sh`)

The batch script generates `resp.inp` with **`CHELPG`** on the main method line, **`%elprop`** (`Dipole` / `Quadrupole`), and **`%output`** with `PrintLevel 3` plus `Print[P_Hirshfeld]` / `Print[P_Mulliken]`. On this ORCA build, tokens **`GTTEN`** and **`Print[P_ESP]`** are **rejected** by the parser; they are omitted. With `%elprop` + `PrintLevel 3`, ORCA still writes **`resp.vpot`** (ESP grid) and **`resp.pc_chelpg`** (CHELPG charges).

### Antechamber note (AmberTools `orcout` + charges)

`antechamber -fi orcaout` is not valid here; use **`-fi orcout`**. The **`resp`** charge model (`-c resp`) requires Gaussian `gout`/`gesp`. This project therefore uses **`-c mul`** to obtain GAFF2 typing from ORCA, then **`apply_pc_chelpg_to_mol2.py`** injects the **CHELPG** values from **`resp.pc_chelpg`** into the mol2 (ESP-fitted, Merz–Kollman-type), which `parmchk2` then consumes.

## Step 3b: ESP RMSE Validation

After `run_antechamber_resp.sh` has produced `gaff2_params/molNN/molNN.mol2` (partial charges from ORCA CHELPG as above), you can quantify how well those point charges reproduce the quantum **electrostatic potential (ESP)** that ORCA computed on a grid. Script `calc_esp_rmse.py` reads the QM ESP grid from **`resp_results/molNN/resp.vpot`** (ORCA 6 two-line header) or, if absent, the `ELECTROSTATIC POTENTIAL` table in `resp.out`. It reads atom coordinates and charges from the mol2 and compares the QM potential to the classical potential Σ q/|r−r_j| at the same grid points. It reports RMSE and related statistics in **kcal/mol/e** (Hartree-based ESP × 627.509).

**Command (from project root):**

```bash
python3 calc_esp_rmse.py --results_dir resp_results --params_dir gaff2_params \
  --manifest resp_inputs/job_manifest.csv --out_csv resp_results/esp_rmse_report.csv
```

**Output:** `resp_results/esp_rmse_report.csv` with per-molecule RMSE, max error, R², and a status (`EXCELLENT` / `GOOD` / `ACCEPTABLE` / `POOR` / `NO_ESP_DATA` / `NO_MOL2`).

**Interpretation:** RMSE below **2** kcal/mol/e is an excellent fit; **2–5** is good for typical MD; **5–10** may be acceptable for large or flexible ligands but deserves a manual check; **10 or higher** (or status `POOR`) means the charge model is not reproducing the QM ESP well—investigate charge state, grid presence, or geometry before using the parameters.

**If `NO_ESP_DATA` for all or many jobs:** confirm `resp.vpot` exists and that `calc_esp_rmse.py` is up to date (ORCA 6 uses a **two-integer** first line in `resp.vpot`). Re-run `run_resp_batch.sh` after fixing `resp.inp` if needed.

**If `POOR`:** inspect the corresponding `resp.out`, `mol2`, and (if present) `resp.vpot`; compare charge sum to the intended net charge; consider re-running ORCA with tighter SCF or a different grid, or revisiting the ligand geometry.

**Exit codes:** the script exits **0** only when there are no missing mol2 files, `NO_ESP_DATA` is not more than 20% of molecules, and no molecule is classified as `POOR`. Otherwise it exits **1**.

