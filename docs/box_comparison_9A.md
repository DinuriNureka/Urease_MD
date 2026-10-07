# 9 Å buffer: octahedral vs rectangular (mol6 monomer)

Same workflow as 10 Å / 12 Å box comparisons: **ff19SB + GAFF2 + TIP3P**, mol6 monomer recipe, `addionsrand` to neutralize.

## Files

| Input | Command | Outputs |
|-------|---------|---------|
| `params/box_comparison/tleap_oct_9A.in` | `solvateoct mol TIP3PBOX 9.0` | `oct_9A.prmtop`, `oct_9A.inpcrd` |
| `params/box_comparison/tleap_rect_9A.in` | `solvatebox mol TIP3PBOX 9.0` | `rect_9A.prmtop`, `rect_9A.inpcrd` |

## Build + report + PDBs

```bash
bash /path/to/urease_md/scripts/compare_box_9A.sh
```

Writes `box_comparison_9A_report.txt`, `oct_9A_full.pdb`, `rect_9A_full.pdb`, and tleap logs `tleap_oct_9A.log`, `tleap_rect_9A.log`.

Use `/usr/local/bin/tleap` if AMBER tools are not on `PATH`.
