# Box comparison (truncated octahedron vs rectangular)

This comparison builds **two** solvated versions of the **same** urease monomer + mol6 ligand system:

- **Truncated octahedron** (`solvateOct`)
- **Rectangular / orthorhombic** (`solvateBox`)

Both use the **same buffer distance** as the current mol6 monomer setup (**12.0 Å**, TIP3P).

Nothing under `simulations/` is modified; all outputs go to `params/box_comparison/`.

## Files

- `params/box_comparison/tleap_octahedral.in`
- `params/box_comparison/tleap_rectangular.in`
- `scripts/compare_box_types.sh` (runs both and writes a report)

## Run

```bash
cd /path/to/urease_md
bash scripts/compare_box_types.sh
```

## Outputs

- `params/box_comparison/oct/solvated.prmtop`
- `params/box_comparison/oct/solvated.inpcrd`
- `params/box_comparison/oct/solvated_full.pdb`
- `params/box_comparison/rect/solvated.prmtop`
- `params/box_comparison/rect/solvated.inpcrd`
- `params/box_comparison/rect/solvated_full.pdb`
- `params/box_comparison/box_comparison_report.txt`

