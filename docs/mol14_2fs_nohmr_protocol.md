# mol14_pose2_monomer — 2 fs no-HMR Parallel Test

**Date created:** 2026-04-21  
**Purpose:** Supervisor-requested validation run at dt=0.002 (2 fs) without HMR, independent of the completed 4 fs HMR run.

---

## Overview

| Parameter | Value |
|---|---|
| System | mol14_pose2_monomer |
| Timestep | 0.002 ps (2 fs) |
| HMR | **No** — standard prmtop |
| Equilibration | 100 ns, single run (no segments) |
| Production | 3 replicates × 20 ns |
| Total MD time | 100 ns equil + 60 ns production |
| Output directory | `simulations/mol14_pose2_monomer_2fs/` |
| Topology | `systems/mol14_pose2_monomer/solvated.prmtop` |
| Starting coords | `systems/mol14_pose2_monomer/solvated.inpcrd` |

> The 4 fs HMR run lives in `simulations/mol14_pose2_monomer/` and is **not touched** by this test.

---

## Files Created

| File | Description |
|---|---|
| `scripts/build_nonhmr_prmtop.sh` | Builds `solvated.prmtop` + `solvated.inpcrd` (no HMR) from tLeap |
| `simulations/input_files/equil_free_2fs.in` | 100 ns NPT dt=0.002, no segments |
| `simulations/input_files/production_2fs_rep.in` | 20 ns NPT dt=0.002, single replicate |
| `scripts/run_mol14_2fs.sh` | Full pipeline orchestrator |

---

## Protocol Stages

```
solvated.inpcrd
    │
    ▼
[min1]          sander, 2000 steps  (solvent+ions restrained)
    │
    ▼
[min2]          sander, 2000 steps  (active site restrained)
    │
    ▼
[heat_nvt_cpu]  pmemd.cuda, 1 ps NVT  0→50 K  (active site restrained)
    │
    ▼
[heat]          pmemd.cuda, 100 ps NVT  50→300 K  (backbone + active site restrained)
    │
    ▼
[equil_cpu_npt] pmemd.cuda, 100 ps NPT  MC barostat  (active site restrained)
    │
    ▼
[equil_free_2fs] pmemd.cuda, 100 ns NPT  dt=0.002  (equil_free_2fs.in)
    │
    ├──► prod_rep1  20 ns  ig=-1
    ├──► prod_rep2  20 ns  ig=-1
    └──► prod_rep3  20 ns  ig=-1
```

All GPU stages: `CUDA_VISIBLE_DEVICES=1`, `$PMEMD_CUDA_BIN` (cuda128 build).

---

## Input File Comparison

### equil_free_2fs.in vs equil_free.in

| Parameter | `equil_free.in` (4 fs HMR) | `equil_free_2fs.in` (2 fs no-HMR) |
|---|---|---|
| `dt` | 0.004 | **0.002** |
| `nstlim` | 2,500,000 | **50,000,000** |
| duration | 10 ns (per segment) | **100 ns (single run)** |
| `mcbarint` | 100 | **25** |
| `iwrap` | — (not set) | **1** |
| `ig` | — | **-1** |
| `ntr` | — | **0** |

### production_2fs_rep.in vs production_2fs_nohmr.in

| Parameter | `production_2fs_nohmr.in` | `production_2fs_rep.in` |
|---|---|---|
| `nstlim` | 5,000,000 (10 ns) | **10,000,000 (20 ns)** |
| header | "10 ns … no segments" | **"20 ns … single run"** |

---

## Step-by-Step Execution

### Step 1 — Build the non-HMR topology (run once)

```bash
cd /path/to/urease_md
bash scripts/build_nonhmr_prmtop.sh mol14_pose2_monomer
```

This will:
- Re-run the same tLeap input (`tleap_pipeline_2026-03-28_122828.in`) with output redirected to `solvated_nonhmr_tmp.*`
- Run `cpptraj fixatomorder`
- Save `systems/mol14_pose2_monomer/solvated.prmtop` + `solvated.inpcrd`
- Create symlinks at `simulations/mol14_pose2_monomer/solvated.prmtop`
- Verify Ni–Ni distance (~3.5 Å)

> **Note:** tLeap's `addionsrand` places ions randomly, so the solvated box will differ from the 4 fs HMR run. This is expected and correct for an independent test.

### Step 2 — Launch the full 2 fs pipeline

```bash
nohup bash scripts/run_mol14_2fs.sh > logs/mol14_2fs.log 2>&1 &
```

The script is **resumable**: if `simulations/mol14_pose2_monomer_2fs/equil/equil.rst7` already exists, all equilibration stages are skipped and only production is run.

### Step 3 — Monitor progress

```bash
tail -f logs/mol14_2fs.log
```

### Step 4 — Check ns/day (after equil completes)

```bash
grep "ns/day" simulations/mol14_pose2_monomer_2fs/equil/equil.out | tail -3
```

---

## Safety Constraints

- `simulations/mol14_pose2_monomer/` is **never touched** by these scripts
- `production.in`, `equil_free.in`, and all HMR topology files are **not modified**
- `CUDA_VISIBLE_DEVICES=1` — uses GPU1 only
- Currently running mol1 queue on GPU1 is **not interrupted** (script queues after)
