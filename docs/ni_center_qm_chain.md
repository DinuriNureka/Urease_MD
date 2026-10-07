# Ni-center ORCA chain (small_opt → small_freq → large_mk)

## Orchestration

| Entry point | Behavior |
|-------------|----------|
| `scripts/run_all_qm.sh` | Runs the Ni-center track (and optionally the ligand batch). The Ni-center steps are implemented by calling `scripts/run_ni_center_chain.sh --from-scratch` so all three ORCA jobs run in order with shared success checks. |
| `scripts/run_ni_center_chain.sh` | **Default:** waits for an already-running `orca small_opt_wB97XD4.inp` (e.g. manual `nohup`), verifies convergence, then runs `small_freq` and `large_mk`. **`--from-scratch`:** runs `small_opt` first (do not use while another `small_opt` is already running). |

Persistent log for the dedicated chain script (and mirrored when invoked from `run_all_qm.sh`):

`logs/ni_center_chain.log`

## Success checks (chain script)

- **small_opt:** `HURRAY`, `THE OPTIMIZATION HAS CONVERGED`, or `OPTIMIZATION RUN DONE` in `small_opt_wB97XD4.out`; GBW present.
- **small_freq:** `.hess` output and `ORCA TERMINATED NORMALLY` or `VIBRATIONAL FREQUENCIES` in `small_freq_wB97XD4.out`.
- **large_mk:** `SCF CONVERGED` and `ORCA TERMINATED NORMALLY` in `large_mk_wB97XD4.out`.

## Environment

- `NI_CENTER_WAIT_PID` — optional PID of the main ORCA driver to wait on.
- `NI_STATUS_FILE` — set by `run_all_qm.sh` for per-step status lines (`small_opt:OK`, etc.).
- `RUN_ALL_QM_CHILD=1` — set internally when the chain is invoked from `run_all_qm.sh` (avoids double-tee of stdout).
