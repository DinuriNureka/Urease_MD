# Auto-start batch build after `min2` completes

This document records the automatic trigger:

- The script waits until `pmemd.cuda` is no longer running with input `min2.in` for:
  - `simulations/mol6_pose2_full_2026-03-24/min2.in`
- Then it starts:
  - `scripts/build_full_systems_batch.py`

The batch build is compute-heavy on the CPU (antechamber + tleap), so we delay it until the current GPU-driven equilibration/minimization stage is no longer using the node.

