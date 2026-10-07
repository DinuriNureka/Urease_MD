# Resilient queue prerequisites (Ni–Ni DISANG)

Date: 2026-04-10 (segment protocol updated 2026-04-16)  
Working directory: `/path/to/urease_md`

## Segment layout (`run_resilient_queue.sh`)

Long NPT after standard equil is **100 ns** of `equil_free`–style dynamics split into **10 × 10 ns** segments (`equil_seg01` … `equil_seg10`, `nstlim=2,500,000`, `dt=0.004` ps). Production is **60 ns** per system: **3 replicates × 2 × 10 ns** (`prod_rep01_seg01` … `prod_rep03_seg02`). Each segment uses the same thermostat / barostat / nonbonded / NMR settings as before; only the **segment count and `nstlim`** changed from the older 5×20 ns equil + 3×20 ns prod layout.

| Phase | Segments | nstlim | dt (ps) | Sim time each | Total |
|-------|------------|--------|---------|-----------------|-------|
| equil_seg01–10 | 10 | 2,500,000 | 0.004 | 10 ns | 100 ns |
| prod_rep01_seg01–02 | 2 | 2,500,000 | 0.004 | 10 ns | 20 ns |
| prod_rep02_seg01–02 | 2 | 2,500,000 | 0.004 | 10 ns | 20 ns |
| prod_rep03_seg01–02 | 2 | 2,500,000 | 0.004 | 10 ns | 20 ns |
| **Total per system** | **16** | | | | **160 ns** |

Launch (when ready): `nohup bash scripts/run_resilient_queue.sh > logs/resilient_queue.log 2>&1 &`

**Migration:** Older runs used `equil_seg1` … `equil_seg5` (no zero padding, 20 ns per folder). The refactored queue only recognizes **`equil_seg01` … `equil_seg10`**. Existing `equil_seg1` directories are **not** read for skip/resume; the first new segment still starts from `equil3` / `gpu_bridge2` / `equil2` / `equil1` as before.

## Task 1 — `systems/<name>/ni_ni.RST`

For every directory under `systems/` that contains `solvated_hmr.prmtop`, `ni_ni.RST` was generated with an `&rst` block (flat-bottom Ni–Ni distance restraint).

**cpptraj logic**

1. Prefer Amber atom indices from `atominfo :570@NI` and `atominfo :571@NI` (matches `NI1_RES` / `NI2_RES` in `config.env` for standard monomer builds).
2. If either mask selects no atom (non-standard numbering or dodecamer topology), fall back to `atominfo @NI`:
   - exactly two Ni atoms → use those indices;
   - more than two → use the two **smallest** atom indices (first dinuclear pair in the current full-system topologies).

**Regenerate** (after `source config.env`):

```bash
cd /path/to/urease_md
source config.env
bash scripts/gen_ni_ni_rst_all.sh
```

## Task 2 — Priority queue readiness

The five systems in `scripts/run_resilient_queue.sh` `QUEUE=(...)` were checked for:

- `systems/<name>/solvated_hmr.prmtop`
- `systems/<name>/ni_ni.RST`
- Starting coordinates, in order: `equil3.rst7` → `gpu_bridge2.rst7` → `equil2.rst7` → `equil1.rst7`; if none, `min2/min2.rst7` (still in early equil).

## Task 3 — DISANG path (`mol13` and others)

`run_resilient_queue.sh` resolves the restraint file with `resolve_ni_rst`:

1. `${PROJECT_ROOT}/systems/${sys}/ni_ni.RST` (primary)
2. `${PROJECT_ROOT}/simulations/${sys}/restraints/ni_ni.RST` (fallback)

The absolute path is written into the generated mdin as `DISANG=...` (see `append_nmropt_disang`).

For `mol13_pose1_monomer`, DISANG therefore points at `systems/mol13_pose1_monomer/ni_ni.RST` once that file exists.

## Counts (2026-04-10)

- **38** system directories had `solvated_hmr.prmtop`; all received `ni_ni.RST`.
- **41** system directories had **no** `solvated_hmr.prmtop` (mostly `_full` builds, `_no_ohb` variants, or poses not yet built).

**Note:** `mol19_pose1_monomer` and `mol33_pose1_monomer` are named in the queue script but have **no** `systems/` or `simulations/` directories on disk; the closest built monomers are `mol19_pose2_monomer` and `mol33_pose3_monomer`. Update `QUEUE` in `run_resilient_queue.sh` when those systems exist or when you intend to run the alternate pose.

**Note:** `rst7_skip_ok` / `cpptraj_nini_ok` in the same script still measures Ni–Ni with masks `:570@NI` and `:571@NI` from `config.env`. That matches monomer prmtops used in the queue but would **not** match `mol21_pose2_monomer` (Ni on residues 808/809) if that system were ever added to the queue without changing the validation masks.
Queue launched successfully: Thu Apr 16 02:59:53 PM EDT 2026
mol13 equil_seg01 running — PID 24905, GPU1, ~160 ns/day
DISANG fix: rst/<sys>.RST real-file copies, path ≤75 chars
