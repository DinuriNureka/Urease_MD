## Stopped runs summary (2026-03-26)

This note records what was running at the time of a manual stop request, and what was stopped.

### What was running

- **Production queue launcher**:
  - `scripts/run_production_queue.sh` was running and launching jobs on **GPU 1 (Tesla P100)**.

- **Active MD job at stop time**:
  - `simulations/mol1_pose2_full_2026-03-25/run_production.sh`
  - Current stage (active `pmemd.cuda`): `min1` (`-i simulations/mol1_pose2_full_2026-03-25/min1.in`)

- **mol6 status (from outputs)**
  - `simulations/mol6_pose2_full_2026-03-24/heat_1ns.out` shows completion:
    - `NSTEP = 500000`, `TIME(PS) = 1000.000` (heat done)

- **Batch system builder**
  - The resumed batch build (`RESUME_FROM_MOL=mol17`, `MAX_CPU_THREADS=20`) had reached **DONE**:
    - `scripts/batch_build_resume_mol17_20260326_034526.log` ends with `done=33/33` and `DONE`.

### Stop actions taken

Stopped all of the following:

- Queue launcher: `pkill -f "scripts/run_production_queue.sh"`
- Per-system run wrappers: `pkill -f "/path/to/urease_md/simulations/.*/run_production.sh"`
- All running AMBER GPU engines: `pkill -f "pmemd.cuda"`

### After-stop check

After the stop, `nvidia-smi pmon` no longer showed any `pmemd.cuda` compute processes running on the P100.

