## Patched (k=50/30) 20 ns crash test — launch

### Correct `pmemd.cuda` binary
After sourcing `config.env`, the intended GPU executable is:

- **`/opt/software/AMBER/pmemd24/bin/pmemd.cuda`**  
  (symlink → `pmemd.cuda_SPFP`)

This matches `AMBER_HOME` in `config.env` and `which pmemd.cuda` on this machine.

**Note:** Using `find … | head -1` to pick a binary is unsafe (order-dependent). Scripts should set:

```bash
PMEMD="${AMBER_HOME}/bin/pmemd.cuda"
```

### Inputs (unchanged)
- MDin: `simulations/input_files/crashtest_20ns.in` (`barostat=2`, `mcbarint=100`, `DISANG=…/ni_ni.RST`)
- Topology: `systems/mol13_pose1_monomer/solvated_hmr.prmtop` (mtime **2026-04-09 07:27**, ~19.4 MB)
- Coords: `simulations/mol13_pose1_monomer/gpu_bridge2/gpu_bridge2.rst7`

### Ni–Ni check (topology + bridge restart)
`cpptraj`, `:570@NI`–`:571@NI`, `noimage`:

- **3.5454 Å**

### Launch command
Run from project root with **`CUDA_VISIBLE_DEVICES=1`**:

```bash
nohup "${AMBER_HOME}/bin/pmemd.cuda" -O \
  -i simulations/input_files/crashtest_20ns.in \
  -p systems/mol13_pose1_monomer/solvated_hmr.prmtop \
  -c simulations/mol13_pose1_monomer/gpu_bridge2/gpu_bridge2.rst7 \
  -o simulations/mol13_pose1_monomer/crashtest/patched_test/patched_test.out \
  -r simulations/mol13_pose1_monomer/crashtest/patched_test/patched_test.rst7 \
  > logs/patched_crashtest_YYYYMMDD_HHMMSS.log 2>&1 &
```

PID / paths recorded in `simulations/mol13_pose1_monomer/crashtest/patched_test/patched_test.pid`.

### ~90 s sanity (this run)
- Process **running**
- First printed block **NSTEP = 5000**: **ANGLE = 6726.6500**
- **NaN** count: **0**

See also `logs/patched_crashtest_verify_20260410.txt`.
