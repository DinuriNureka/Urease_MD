# Troubleshooting Guide

Every issue follows: **SYMPTOM → CAUSE → FIX**.

---

## Phase 1

### T1.1 — pdb4amber warns about non-standard residues
**SYMPTOM:** Warning about KCX, NI. **CAUSE:** Expected — these are non-standard. **FIX:** No action. Handled by MCPB.py in Phase 2.

### T1.2 — reduce protonates metal-coordinating nitrogen
**SYMPTOM:** HIP (doubly protonated His) near Ni. **CAUSE:** reduce ignores metal coordination. **FIX:** Run `fix_protonation.py` (Step 1.4). Verify in PyMOL: coordinating N must have no H.

### T1.3 — Bridging OH⁻ absent
**SYMPTOM:** No O atom within 2.5 Å of both Ni. **CAUSE:** Low cryo-EM resolution. **FIX:** Add OH⁻ at midpoint, ~2.0 Å from each Ni. Use PyMOL Builder or Avogadro.

---

## Phase 2

### T2.1 — MCPB.py KeyError on second Ni
**SYMPTOM:** `KeyError: 'NI'`. **CAUSE:** Both metals share residue name. **FIX:** NI1 and NI2 must have DIFFERENT residue names in mol2 and PDB.

### T2.2 — Imaginary frequencies in ORCA output
**SYMPTOM:** `***imaginary mode***`. **CAUSE:** Not at true minimum. **FIX:** If <3 small modes (<50 cm⁻¹): acceptable. If large (>100 cm⁻¹): re-optimise with `VeryTightOpt`.

### T2.3 — ORCA out of memory
**SYMPTOM:** SHARK allocation error. **CAUSE:** maxcore × nprocs > RAM. Gradient steps use ~2× SCF memory. **FIX:** Reduce `%maxcore` to 1200 with `nprocs 28` (=33.6 GB, leaves ~30 GB for gradient/OS). We hit OOM at maxcore=1700 on our 62 GB system.

### T2.4 — MCPB.py Step 3 fails with ORCA
**SYMPTOM:** Zero charges or crash. **CAUSE:** Incomplete ORCA ESP reader. **FIX:** Use multi-frame RESP (Step 2.7b in `02b_multiframe_resp.md`) or single-frame Multiwfn. Update mol2 manually. Skip Step 3, run Step 4.

### T2.5 — Unreasonable force constants
**SYMPTOM:** k < 20 or > 300 kcal/mol/Å². **CAUSE:** Hessian conversion error. **FIX:** Re-check `orca_hess_to_fchk.py` atom count. Verify optimised geometry was used. Expected: Ni–N 50–120, Ni–O 60–130 kcal/mol/Å².

### T2.6 — MCPB.py assigns wrong charge/multiplicity
**SYMPTOM:** MCPB.py Gaussian templates show charge=+2, mult=2. **CAUSE:** MCPB.py's automatic charge assignment doesn't handle dinuclear metal sites correctly. For 71 atoms (Z=329, odd): charge=+2 gives 327 electrons (odd), requiring even multiplicities; but the correct FM quintet (mult=5) needs even electrons. **FIX:** Use charge=+1, mult=5 (FM) and charge=+1, mult=1 (BS singlet). These give 328 electrons (even), consistent with both multiplicities.

### T2.7 — FM geometry optimization oscillates / won't converge
**SYMPTOM:** MAX gradient stuck at ~0.01, energy oscillating ±0.0001 Eh for 50+ cycles. **CAUSE:** Approximate Hessian doesn't capture coupling in dinuclear Ni²⁺ center. The optimizer takes steps in wrong directions. **FIX:** Add `Calc_Hess true` and `Recalc_Hess 10` to `%geom` block. This computes the exact Hessian at the first step and recomputes every 10 cycles. Costs ~1–2 h upfront but converges in ~10–20 cycles after. Do NOT reduce MaxStep too aggressively (e.g. 0.003) — steps become too small to make progress.

### T2.8 — ORCA GBW corrupted after kill
**SYMPTOM:** SCF starts from scratch (100+ iterations) after restarting. **CAUSE:** Killing ORCA mid-cycle corrupts the GBW file. The next run cannot read converged orbitals. **FIX:** Avoid killing ORCA unless necessary. Each restart wastes ~15–20 min rebuilding orbitals. If you must restart, back up the output but understand the GBW may be unusable.

### T2.9 — ORCA MD `Dump` keyword error
**SYMPTOM:** `Unknown identifier: DUMP` or `Did not find a match for function "Dump"`. **CAUSE:** ORCA 6 MD module uses SANscript syntax, not keyword=value. **FIX:** Use `dump position stride 50 filename "trajectory.xyz"` and `run 20000` instead of `Dump 50.0_fs` and `RunTime 20.0`.

---

## Phase 3

### T3.1 — antechamber atom typing fails
**SYMPTOM:** "Cannot determine atom types". **CAUSE:** Missing H or bad bonds. **FIX:** Add H: `obabel -h`. Try SDF format. Check bond orders.

### T3.2 — sqm AM1-BCC fails to converge
**SYMPTOM:** SCF failure in sqm.out. **CAUSE:** AM1 struggles with certain groups. **FIX:** Re-optimise geometry first. Use ORCA HF/6-31G* RESP for IC₅₀ molecules.

### T3.3 — parmchk2 ATTN warnings
**SYMPTOM:** "ATTN: need revision" in frcmod. **CAUSE:** Parameters by analogy. **FIX:** Acceptable if far from binding site. For active-site-facing atoms: compute from QM.

---

## Phase 4

### T4.1 — Missing bond parameter in tleap
**SYMPTOM:** "Could not find bond parameter". **CAUSE:** Atom types don't match. **FIX:** Verify `loadoff urease.lib` and `loadamberparams urease.frcmod` in tleap script. Use `desc mol.<resnum>`.

### T4.2 — Non-zero total charge
**SYMPTOM:** |charge| > 0.01. **CAUSE:** Wrong mol2 charges or antechamber `-nc`. **FIX:** Check `charge mol` in tleap before ions. Verify all mol2 charge sums.

---

## Phase 5

### T5.1 — Crash during heating (SHAKE/NaN)
**SYMPTOM:** SHAKE error or NaN in output. **CAUSE:** Clashes, poor minimisation. **FIX:** More minimisation (maxcyc=20000). Slower heating (200 ps). Add intermediate MIN stage.

### T5.2 — Density doesn't converge
**SYMPTOM:** Density stuck at 0.85–0.95. **CAUSE:** Slow pressure coupling. **FIX:** Reduce taup to 0.5 ps. Extend EQUIL1 to 500 ps.

### T5.3 — Ni coordination distorts
**SYMPTOM:** Ni–N > 2.5 Å or Ni–Ni > 4.5 Å. **CAUSE:** Bonds not created in tleap. **FIX:** Check tleap.log for "Creating bond". Increase metal restraints to 50.0 during heating.

---

## Phase 6

### T6.1 — GPU memory error
**SYMPTOM:** Allocation failure. **CAUSE:** Unlikely for 200K atoms on ≥16 GB GPU. **FIX:** Reduce `cut=8.0`. Check `nvidia-smi` for contention.

### T6.2 — Slow performance
**SYMPTOM:** <20 ns/day. **CAUSE:** Wrong engine or GPU contention. **FIX:** Verify `pmemd.cuda` (not sander). Check `nvidia-smi dmon -s u`.

---

## Phase 7

### T7.1 — MMPBSA.py atom count mismatch
**SYMPTOM:** Crash at trajectory load. **CAUSE:** Topology mismatch. **FIX:** Regenerate stripped topologies from the EXACT production topology (solvated_hmr.parm7).

### T7.2 — cpptraj "Mask has no atoms"
**SYMPTOM:** Empty selection. **CAUSE:** Residue renumbering. **FIX:** Query `resinfo :NI*` in cpptraj. Update scripts with correct residue numbers.

---

## General Debugging

```bash
grep -iE "error|fatal|nan|shake|vlimit" *.out  # Scan for problems
watch -n 5 nvidia-smi                           # Monitor GPU
cpptraj -p sys.parm7 -y traj.nc <<< "trajinfo" # Check trajectory
process_mdout.perl *.out                        # Energy trends
```
