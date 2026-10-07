# 10 Å buffer: octahedral vs rectangular (mol6 monomer)

This step builds two **test** solvated systems for comparison only (no changes to production `simulations/` files).

## Inputs

| Script | Solvation | Buffer | Outputs |
|--------|-----------|--------|---------|
| `params/box_comparison/tleap_oct_10A.in` | `solvateoct` | 10.0 Å (`TIP3PBOX 10.0`) | `oct_10A.prmtop`, `oct_10A.inpcrd` |
| `params/box_comparison/tleap_rect_10A.in` | `solvatebox` | 10.0 Å | `rect_10A.prmtop`, `rect_10A.inpcrd` |

Force fields and assembly match production monomer builds: **ff19SB**, **GAFF2** (ligand), **TIP3P**, MCPB Ni-center mol2/frcmod, `monomer_BQ_apo.pdb` + LIG, `bonds_monomer.txt` / `set_patch_monomer.txt`.

Ions: `addionsrand mol Na+ 0` and `addionsrand mol Cl- 0` to neutralize net charge (same pattern as `tleap_monomer.in`).

## Run comparison

```bash
/path/to/urease_md/scripts/compare_box_10A.sh
```

This runs both tleap jobs, then writes **`params/box_comparison/box_comparison_10A_report.txt`** with NATOM, WAT count (from tleap log), volume line, box dimensions (parmed), and percent atom-count difference (rect vs oct).

Logs: `tleap_oct_10A.log`, `tleap_rect_10A.log` in the same directory.

## tleap path

Use `/usr/local/bin/tleap` if your shell does not already put AMBER tools on `PATH`.
