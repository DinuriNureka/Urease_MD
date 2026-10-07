# Fixed rectangular box (mol6 monomer + LIG)

## Why not exactly 100 × 110 × 100 Å?

Axis-aligned bounding box of the **dry** solute (after `bonds_monomer.txt` / `set_patch_monomer.txt`) is about **68 × 80 × 103 Å**. A **100 Å** edge along **z** is **too small** (103 Å > 100 Å), so a 100 × 110 × 100 cell cannot enclose the solute without rotation.

Use a **minimum** third edge **≥ ~104 Å**, e.g. **100 × 110 × 108** Å.

## Build (automated)

```bash
bash /path/to/urease_md/scripts/build_rect_box_custom.sh 100 110 108
```

This:

1. Runs tleap to write a **dry** PDB (`mol6_dry_for_box_100x110x108.pdb`).
2. Computes `setBox mol centers {cx cy cz}` so periodic edges match the target (tleap: extent + 2×clearance per axis).
3. Runs `solvatebox mol TIP3PBOX 0.0`, neutralizing ions, then **`ambpdb`** for a full-system PDB.

Outputs (example for `100 110 108`):

- `params/box_comparison/rect_100x110x108.{prmtop,inpcrd}`
- `params/box_comparison/rect_100x110x108_full.pdb`
- `params/box_comparison/box_rect_100x110x108_report.txt`

**Note:** After `solvatebox`, tleap reports a **total VdW box** that is **larger** than the `setBox` target (e.g. target **100 × 110 × 108** Å can yield **~103 × 113 × 111** Å in `BOX_DIMENSIONS` and the last line of `inpcrd`). That is the **periodic** box used by the engine, not a bug. To land closer to a desired edge length, iterate the target passed to `build_rect_box_custom.sh` (smaller target → smaller final VdW box) or adjust the solvation strategy.

## Helpers

- `scripts/compute_clearance_for_box.py` — prints `{cx cy cz}` for `setBox` from a dry PDB and target **a b c**.
