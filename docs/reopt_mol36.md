# mol36 (thiourea) ORCA reoptimization from planar NH₂ guess

## Files

| File | Role |
|------|------|
| `resp_results/mol36/mol36_planar.xyz` | Starting geometry: H–N–C–S dihedrals set to **0°** and **180°** per NH₂ (flattened from `params/inhibitors/orca_opt/mol36/mol36.xyz`). |
| `resp_results/mol36/mol36_reopt.inp` | `ωB97X-D4` / `def2-SVP`, `Opt Freq`, `TightOpt`, `TightSCF`, `DEFGRID3`, `nprocs 6`, `maxcore 3000`. |
| `scripts/reopt_mol36.sh` | Runs ORCA; if frequencies are clean (no imaginary modes), updates both `mol36.xyz` copies and writes `REDO_MULTIWFN_RESP.txt`. |

**Grid keywords:** ORCA 6.1 does not accept `Grid5` / `FinalGrid6` on the simple line; **`DEFGRID3`** is used (same convention as `mol02_reopt.inp`).

## Runner

```bash
bash scripts/reopt_mol36.sh
```

## Frequency check

If analytic frequencies still show an imaginary mode near **±180 cm⁻¹**, the stationary point can remain a **saddle** along the **NH₂ inversion / umbrella** coordinate despite a converged optimizer — a known difficulty for shallow, anharmonic PESs at harmonic frequency level. Options for follow-up: numerical frequencies, larger basis, or literature gas-phase minimum structure as a starting guess.

## Project decision (thiourea / mol36)

For **gas-phase thiourea** at ωB97X-D4/def2-SVP, a **~−180 cm⁻¹** harmonic mode is consistent with a **strongly anharmonic** NH₂ inversion / umbrella coordinate: the harmonic approximation is a poor model here, not necessarily a “bad” geometry for **RESP partial charges**.

**Accepted for production:** keep the **current Multiwfn RESP charges** for **mol36**. The batch report lists **rmse_charge_e ≈ 0.00242 e** and **RRMSE ≈ 0.095** (well within usual thresholds). **No further ORCA reoptimization or RESP redo is required** for mol36 unless the scientific question changes (e.g. explicit solvation, different level of theory).
