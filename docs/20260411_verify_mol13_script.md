## `scripts/verify_mol13.sh`

Read-only diagnostics for `mol13_pose1_monomer` while **equil_seg2** (or other stages) run. Does **not** write under `simulations/mol13_pose1_monomer/`; temp files use `/tmp/mol13_*`.

### Run

```bash
cd /path/to/urease_md
bash scripts/verify_mol13.sh
```

Full log: `logs/mol13_verification_YYYYMMDD_HHMMSS.txt`

### cpptraj path

The script first tries `/opt/software/AMBER/pmemd24/bin/cpptraj` (per project convention). On this host that binary is often absent; it **falls back** to `$(command -v cpptraj)` (typically `/usr/local/bin/cpptraj`) and records a WARN at the top of the report.

### Checks (summary)

1. **Residue masks** — `resinfo` for key residues + empty `:OHB`
2. **Ni–Ni trajectory** — `patched_test.nc` or `equil_seg1.nc` if present; else WARN
3. **DISANG / NMR** — header + non-zero NMR bond energy in `equil_seg2.out`
4. **equil_cpu_npt** — `:572` identity + restraint mask notes (no auto-edit)
5. **HMR** — max solute H mass via `atominfo` (expect ~3.024 Da)
6. **Box** — orthorhombic vs octahedral from cpptraj + prmtop
7. **`nmr_disang/` symlink** — resolves to `restraints/ni_ni.RST`, `iat=`, DISANG line length
8. **Density / temperature** — from live `equil_seg2.out` (early-run density uses a slightly wider WARN band if &lt;500 samples)
9. **GPU** — `CUDA_VISIBLE_DEVICES` and `nvidia-smi` PID→UUID vs GPU index 1

Exit **0** if no FAIL (WARN allowed); **1** if any FAIL.
