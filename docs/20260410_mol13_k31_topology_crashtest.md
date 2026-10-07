## Purpose
Build **mol13** solvated topology with uniform **k = 31 kcal/mol/rad²** Ni angle terms (`urease_no_ohb_mcpbpy_k31.frcmod`), apply **fixatomorder** + **HMR**, and launch the **20 ns** `crashtest_20ns.in` run (NPT, MC barostat, `mcbarint=100`) on **GPU 1** from `gpu_bridge2/gpu_bridge2.rst7`.

## k31 frcmod repair
`params/ni_center_no_ohb/mcpb_input/urease_no_ohb_mcpbpy_k31.frcmod` had two **corrupted** `ANGL` header lines (atom types were merged with the force constant), which made `tleap` report dozens of “Could not find angle parameter” errors.

Corrected to match the working ABS frcmod atom labels:

- `C -Y1-M1        31.00     122.82`
- `C -Y6-M2        31.00     122.07`

## Tleap input
Added `systems/mol13_pose1_monomer/tleap_k31_20260410.in`:

- Loads `urease_no_ohb_mcpbpy_k31.frcmod`
- Keeps `abs_missing_dihedrals.frcmod` (same Y7-CO-2C-XC patch as ABS build)
- Writes **`monomer_mol13_k31.prmtop` / `monomer_mol13_k31.inpcrd`** so **`monomer_mol13.*` (ABS reference)** is not overwritten

Log: `systems/mol13_pose1_monomer/tleap_k31_20260410.log` → **Errors = 0**

## fixatomorder + HMR
Used the same safe pattern as the ABS workflow:

```bash
cpptraj -p systems/mol13_pose1_monomer/monomer_mol13_k31.prmtop <<'EOF'
trajin systems/mol13_pose1_monomer/monomer_mol13_k31.inpcrd
fixatomorder parmout /tmp/mol13_k31_fixed.prmtop
trajout /tmp/mol13_k31_fixed.rst7 restart
run
quit
EOF
```

Then `parmed` HMR on `/tmp/mol13_k31_fixed.{prmtop,rst7}` →

- `systems/mol13_pose1_monomer/solvated_k31_hmr.prmtop`
- `systems/mol13_pose1_monomer/solvated_k31_hmr.rst7`

## Geometry / HMR checks
- Ni–Ni (`:570@NI`–`:571@NI`, `noimage`): **3.5105 Å** (before and after HMR on the built-inpcrd frame)
- HMR: max H mass **3.024**; Ni masses **[58.69, 58.69]**
- Total charge (cpptraj): **-0.356352**

## Crash test launch
Stopped any existing `pmemd.cuda` tied to `mol13_pose1_monomer`, then:

```bash
CUDA_VISIBLE_DEVICES=1 nohup ${AMBER_HOME}/bin/pmemd.cuda -O \
  -i simulations/input_files/crashtest_20ns.in \
  -p systems/mol13_pose1_monomer/solvated_k31_hmr.prmtop \
  -c simulations/mol13_pose1_monomer/gpu_bridge2/gpu_bridge2.rst7 \
  -o simulations/mol13_pose1_monomer/crashtest/k31_test/k31_test.out \
  -r simulations/mol13_pose1_monomer/crashtest/k31_test/k31_test.rst7 \
  > logs/k31_crashtest_20260410_063427.log 2>&1 &
```

PID recorded in `simulations/mol13_pose1_monomer/crashtest/k31_test/k31_test.pid` (first line).

## First printed energy block (NSTEP = 5000)
- **ANGLE = 6621.8476**
- No **NaN** in early output
