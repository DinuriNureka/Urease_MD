# Production protocol + queue launcher

## Production protocol (standard for all compounds)

- **Minimization**
  - `min1.in`: restrained solvent relaxation
  - `min2.in`: unrestrained minimization
- **Heating**
  - `heat_1ns.in`: **1 ns** NVT temperature ramp from 0→300 K
  - `nstlim=500000`, `dt=0.002 ps`
- **Equilibration**
  - `equil1.in`, `equil2.in`, `equil3.in`: **2 ns each** (NPT)
  - `nstlim=1000000`, `dt=0.002 ps`
  - total equil time: **6 ns**
- **Production**
  - `prod_10ns.in`: 10 ns (HMR topology, `dt=0.004 ps`)
  - `prod_extend_40ns.in`: +40 ns extension
  - total production time: **50 ns**
  - **Selected compounds** (mol1, mol2, mol3, mol4, mol5, mol6, mol7, mol8, mol15, mol16, mol17, mol18):
    - add `prod_extend_50ns.in`: +50 ns extension
    - total production time: **100 ns**

## Per-compound entrypoint

Each compound should have an executable:

- `simulations/<compound_run_dir>/run_production.sh`

The queue launcher auto-discovers these and runs them.

## Queue launcher

Script:
- `scripts/run_production_queue.sh`

Behavior:
- Runs **one job per GPU** (permitted GPUs) and queues the rest.
- Enforces CPU headroom: uses at most **80% of cores** (leaves ≥20% free) by setting `OMP_NUM_THREADS`.
- Enforces GPU headroom: only launches a job on a GPU if free VRAM ≥ `MIN_FREE_VRAM_MB`.
- Periodically re-discovers new jobs: rescans `simulations/*/run_production.sh` every `RE_DISCOVER_SECONDS` (default: **1800s = 30 min**) so newly built systems are picked up automatically.
- Automatically sets `PROD_TARGET_NS=100` for the selected compounds list above, otherwise `PROD_TARGET_NS=50`.

Single-GPU safety (to avoid OOM):
- If you run with a single allowed GPU (e.g. `GPU_ALLOWLIST="1"`), the queue will not start a new job while `pmemd.cuda` is currently running.

### Run

```bash
bash scripts/run_production_queue.sh
```

### Common tuning

- Run only on P100 (GPU 1):

```bash
GPU_ALLOWLIST="1" bash scripts/run_production_queue.sh
```

- Require more free VRAM (MB) before launching a new job:

```bash
MIN_FREE_VRAM_MB=9000 bash scripts/run_production_queue.sh
```

- Change polling interval:

```bash
POLL_SECONDS=60 bash scripts/run_production_queue.sh
```

