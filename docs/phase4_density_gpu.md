# Phase 4d NPT density on GPU (pmemd.cuda)

## Symptom

`density.log` stops right after **`4. RESULTS`** with no `NSTEP` lines, `pmemd` exits, and no `density.rst7`. Stderr can show:

```text
ERROR: Calculation halted. Periodic box dimensions have changed too much from their initial values.
...
The GPU code does not automatically reorganize grid cells...
```

This is **not** OOM and not a bad HMR topology by default: the GPU PME grid is fixed at startup, and the jump from **NVT heat** to **NPT** (Berendsen) can change the box quickly enough that **4 fs** triggers the halt.

## Fix (in `run_full_pipeline.sh`)

Phase 4d **density** uses **`dt=0.002`** and **`nstlim=1000000`** so the simulated time stays **2 ns** (same as the old `500000 × 0.004`). Restraints, barostat, and masks are unchanged.

**Phase 4e release** also used **GPU halt at 4 fs** after density; it now uses **`dt=0.002`**, **`nstlim=1000000`** (2 ns). **Phase 4f test5ns** uses **`dt=0.002`**, **`nstlim=2500000`** (5 ns), **`ntwx=5000`** (same ~10 ps spacing as before).

### Truncated-octahedron barostat

**Monte Carlo barostat** (`ntp=2`, `barostat=2`) is **not allowed** for this non-orthorhombic box; pmemd errors with *“Nonisotropic scaling on nonorthorhombic unit cells is not permitted.”* **test5ns**, Phase 5 **eq100ns**, and **prod** inputs use **Berendsen** (`ntp=1`, `barostat=1`) instead.

## mdout parsing

Some pmemd versions append a timing/summary block with a spurious line `Density = 0.0009` (wrong volume). Do not use `grep Density | tail -1`; `run_full_pipeline.sh` uses `amber_last_density_g_cm3()` and keeps only values in **0.85–1.2 g/cm³**.

## Manual check

Use short symlinks for topology (`density.ptop`), not a bare `solvated_hmr.prmtop` name in the phase directory, if paths are long (Fortran `OPEN` limits).

```bash
cd .../phase4_.../
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES=1
/opt/software/AMBER/pmemd24/bin/pmemd.cuda -O -i density.in -p density.ptop -c heat.rst7 -r density.rst7 -o density.log -ref heat.rst7
```

## System logs

`dmesg` may require **root** (`sudo dmesg`) to read the kernel ring buffer; absence of OOM lines in unprivileged checks does not rule out other issues, but this failure mode is reported by **pmemd** on stderr.
