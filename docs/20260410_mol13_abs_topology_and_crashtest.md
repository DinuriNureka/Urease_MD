## Goal
Build a valid ABS topology for `mol13_pose1_monomer` and launch the **20 ns ABS crash test** (NPT, **MC barostat**) starting from `gpu_bridge2/gpu_bridge2.rst7` on **GPU 1**.

## Work directory
`/path/to/urease_md`

## Step 1 — Diagnose `tleap` failure (root cause)
`tleap` initially produced a **0-byte** `systems/mol13_pose1_monomer/monomer_mol13.prmtop`.

Root cause in `systems/mol13_pose1_monomer/tleap_abs_20260410.log`:
- `Error! ** No torsion terms for atom types: Y7-CO-2C-XC`
- `Parameter file was not saved.`

## Step 2 — Patch missing torsion term and rebuild
Added ABS-only frcmod patch:
- `systems/mol13_pose1_monomer/abs_missing_dihedrals.frcmod`
  - Defines `DIHE` term `Y7-CO-2C-XC` with **0.000** barrier (neutral “no-op” torsion; satisfies LEaP).

Loaded it in:
- `systems/mol13_pose1_monomer/tleap_no_ohb_mol13_pose1_monomer.in`

Re-ran:

```bash
cd systems/mol13_pose1_monomer
tleap -f tleap_no_ohb_mol13_pose1_monomer.in > tleap_abs_20260410.log 2>&1
```

Result:
- `monomer_mol13.prmtop`: **~19 MB**
- `monomer_mol13.inpcrd`: **~3.5 MB**
- `Exiting LEaP: Errors = 0`

## Step 3 — `fixatomorder` (write matching prmtop + rst7)
Used `fixatomorder parmout` to ensure the reordered coordinates and topology stay consistent:

```bash
cpptraj -p systems/mol13_pose1_monomer/monomer_mol13.prmtop <<'EOF'
trajin systems/mol13_pose1_monomer/monomer_mol13.inpcrd
fixatomorder parmout /tmp/mol13_abs_fixed.prmtop
trajout /tmp/mol13_abs_fixed.rst7 restart
run
quit
EOF
```

Ni atoms after fixatomorder are in residues **570** and **571**:
- `mask @NI` → `:570@NI` and `:571@NI`

Ni–Ni distance checks:
- Pre-fix (original topology, residues 808/809): **3.5105 Å**
- Post-fix (fixed topology, residues 570/571): **3.5105 Å**

## Step 4 — Apply HMR (correct prmtop/rst7 pairing)
HMR inputs:
- `/tmp/mol13_abs_fixed.prmtop`
- `/tmp/mol13_abs_fixed.rst7`

HMR outputs:
- `systems/mol13_pose1_monomer/solvated_abs_hmr.prmtop`
- `systems/mol13_pose1_monomer/solvated_abs_hmr.rst7`

Verification:
- Max H mass: **3.024**
- Ni masses: **[58.69, 58.69]**
- Ni–Ni in HMR topology (`:570@NI`–`:571@NI`): **3.5105 Å**

## Step 5 — Crash test input
Used existing:
- `simulations/input_files/crashtest_20ns.in`

Key settings:
- `nstlim=5000000, dt=0.004` (20 ns)
- `ntb=2, ntp=1`
- `ntt=11, tautp=2.0`
- `barostat=2, mcbarint=100, pres0=1.0`
- `DISANG=/path/to/urease_md/ni_ni.RST`

## Step 6 — Launch ABS crash test (GPU 1)
Starting restart (hard constraint):
- `simulations/mol13_pose1_monomer/gpu_bridge2/gpu_bridge2.rst7`

Launched:

```bash
CUDA_VISIBLE_DEVICES=1 nohup ${AMBER_HOME}/bin/pmemd.cuda -O \
  -i simulations/input_files/crashtest_20ns.in \
  -p systems/mol13_pose1_monomer/solvated_abs_hmr.prmtop \
  -c simulations/mol13_pose1_monomer/gpu_bridge2/gpu_bridge2.rst7 \
  -o simulations/mol13_pose1_monomer/crashtest/abs_test/abs_test.out \
  -r simulations/mol13_pose1_monomer/crashtest/abs_test/abs_test.rst7 \
  -x simulations/mol13_pose1_monomer/crashtest/abs_test/abs_test.nc \
  > logs/abs_crashtest_20260410_035900.log 2>&1 &
```

PID:
- **2989049** (also recorded in `simulations/mol13_pose1_monomer/crashtest/abs_test/abs_test.pid`)

## Step 7 — First energy / ANGLE sanity
First printed block (NSTEP=5000):
- **ANGLE = 6783.5139**
- NaN count: **0**

