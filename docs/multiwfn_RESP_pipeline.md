# Multiwfn RESP pipeline (ORCA → molden → Multiwfn → mol2)

This replaces ORCA **CHELPG** charges in `gaff2_params/molNN/molNN.mol2` with **Multiwfn RESP** charges derived from the **SCF wavefunction** in `resp.gbw` (via `resp.molden`).

## Prerequisites

- **ORCA** `orca_2mkl` (default: `/opt/software/ORCA/build/orca_2mkl`)
- **Multiwfn** executable — set `MULTIWFN` to the full path if not in `PATH`
- Python 3.10+
- Existing `resp_results/molNN/resp.gbw` and GAFF2 `gaff2_params/molNN/molNN.mol2`

## ORCA molden note

`orca_2mkl resp -molden` writes **`resp.molden.input`**. The batch script renames it to **`resp.molden`** for Multiwfn.

## RESP fitting metrics (RMSE / RRMSE)

After a RESP or one-stage ESP fit, Multiwfn prints a **single summary line** (no unit labels on that line), for example:

```text
 RMSE:    0.002217   RRMSE:    0.120512
```

- **RMSE** — reported here as **`rmse_charge_e`** in project outputs. Multiwfn does not print a unit on this line; it is interpreted as **root mean square error in elementary charge units (e)** for the fit, consistent with Multiwfn RESP/ESP fitting examples.
- **RRMSE** — relative RMSE, **dimensionless**.

**Acceptable thresholds** (project defaults, aligned with Multiwfn RESP example workflows):

| Quantity        | Threshold   |
|----------------|------------|
| RMSE           | **&lt; 0.005 e** |
| RRMSE          | **&lt; 0.3**     |

`scripts/parse_multiwfn_resp.py` captures these from the log into `resp_multiwfn_metrics.txt`; `scripts/collect_multiwfn_rmse.py` writes `resp_results/multiwfn_resp_rmse_report.csv`.

## Steps

1. **GBW → molden**

   ```bash
   bash scripts/gbw_to_molden.sh
   ```

2. **Multiwfn RESP** (non-interactive menu: main menu 7 → 18 → … as in `scripts/run_multiwfn_resp.sh`)

   ```bash
   export MULTIWFN=/path/to/Multiwfn
   bash scripts/run_multiwfn_resp.sh
   ```

   Outputs per molecule: `resp_multiwfn_resp.log`, `resp_multiwfn.chg`, `resp_multiwfn_metrics.txt`.

3. **ESP map (optional)**

   ```bash
   bash scripts/run_multiwfn_esp_map.sh
   ```

   Uses `xvfb-run` if available. If mapping fails, check `esp_map_run.log` in each `molNN` directory.

4. **Inject charges into mol2**

   ```bash
   bash scripts/apply_multiwfn_charges_batch.sh
   ```

5. **RMSE / RRMSE report**

   ```bash
   python3 scripts/collect_multiwfn_rmse.py
   ```

   CSV columns: `mol_id`, `rmse_charge_e`, `rrmse`, `pass_rmse`, `pass_rrmse`, `status`. The script also prints a short summary (means, min/max RMSE, pass and WARN counts).

6. **Verify**

   ```bash
   python3 scripts/verify_multiwfn_pipeline.py
   ```

   Checks that every `gaff2_params/molNN/molNN.mol2` exists and that partial charges sum to **`net_charge`** from `resp_inputs/job_manifest.csv`. Optionally warns if `rmse_charge_e` exceeds a soft ceiling (default 0.1 e).

**Note (mol36 / thiourea):** Analytic harmonic frequencies for isolated thiourea can show a spurious **~−180 cm⁻¹** mode along a **strongly anharmonic** NH₂ coordinate; **RESP metrics for mol36 are acceptable** (e.g. **rmse_charge_e ≈ 0.0024 e**). See `docs/reopt_mol36.md` — **no extra RESP redo** is required for mol36 for the current MD workflow.

## One-command orchestration

```bash
export MULTIWFN=/path/to/Multiwfn
bash scripts/run_multiwfn_full_pipeline.sh
```

## AmberTools note

`antechamber -c resp -fi orcout` is not supported for Gaussian-style RESP on this installation; **Multiwfn** performs the RESP fit on the molden wavefunction. If `parse_multiwfn_resp.py` cannot parse your Multiwfn version, edit patterns in that script or paste charges manually into `resp_multiwfn.chg`.
