## Post–equil3 long NPT + production (36 systems)

### Protocol (matches validated 20 ns crash test)
Per-segment / per-replicate MDIN (`equil_seg.in` and `production.in` are identical):

| Parameter | Value |
|-----------|--------|
| `nstlim` | 5,000,000 |
| `dt` | 0.004 ps → **20 ns**/segment |
| `ntc` / `ntf` | 2 / 2 |
| `ntb` / `ntp` | 2 / 1 |
| `barostat` | 2 (Monte Carlo) |
| `mcbarint` | **100** |
| `pres0` | 1.0 atm |
| `ntt` / `tautp` / `temp0` | 11 / 2.0 ps / 300 K |
| `cut` | 9.0 Å |
| `ntpr` / `ntwx` / `ntwr` | 5000 / 5000 / 25000 |
| `ioutfm` / `ntxo` | 1 / 2 |
| `nmropt` | 1 |
| `ntr` | 0 |
| `ig` | -1 |
| `DISANG` | **absolute** path → short symlink `${PROJECT_ROOT}/nmr_disang/<SYS>.RST` (points at `simulations/<SYS>/restraints/ni_ni.RST`; avoids pmemd Fortran line-length truncation) |

**Executable:** `/opt/software/AMBER/pmemd24/bin/pmemd.cuda`  
**GPU:** `CUDA_VISIBLE_DEVICES=1`

### `ni_ni.RST`
Flat-bottom Ni–Ni restraint; `iat=` filled per system via `cpptraj` (`mask @NI`) on a matching `rst7`:

```
rk2=50.0, rk3=50.0
```

### Scripts
| Script | Role |
|--------|------|
| `scripts/post_equil_lib.sh` | Shared helpers (NI extract, mdin, pmemd wrapper, mdout checks) |
| `scripts/build_inputs_for_system.sh <SYS>` | Requires `equil3/equil3.rst7`; writes `restraints/ni_ni.RST` + `input_files/equil_seg.in` + `input_files/production.in` |
| `scripts/run_mol13_now.sh` | mol13 only: seg1 = existing `crashtest/patched_test/patched_test.rst7`; runs seg2–5 + rep1–3 |
| `scripts/run_system_full.sh <SYS>` | Generic: seg1←`equil3.rst7`, seg2–5 chain, rep1–3←`equil_seg5.rst7` |
| `scripts/queue_all_36.sh` | Serial queue over all `systems/*/solvated_hmr.prmtop` (**skips mol13**) |
| `scripts/status_all_36.sh` | E3 / S1–S5 / R1–R3 file checks |

### Launch (typical)
```bash
cd /path/to/urease_md
nohup bash scripts/run_mol13_now.sh > logs/mol13_nohup.log 2>&1 &
# After mol13 finishes:
nohup bash scripts/queue_all_36.sh > logs/queue_all_36_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```

### Completion checks (per segment)
- `pmemd` exit 0  
- `grep -c NaN` on `*.out` = 0  
- `grep -E 'NSTEP[[:space:]]*=[[:space:]]*5000000'`
