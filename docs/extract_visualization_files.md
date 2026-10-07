# Extract visualization PDBs

This project often needs three PDBs for clean visualization from an AMBER solvated system:

- **Full system**: protein + ligand + Ni center + water + ions
- **No-water**: protein + ligand + Ni center (no waters/ions)
- **Ligand watershell**: protein + ligand + Ni center + waters within **5 Å** of ligand

The script `scripts/extract_visualization_files.py` generates all three and copies them to
`/mnt/user-data/outputs/` for download when that path exists and is writable.

If that directory is unavailable (typical on a bare-metal workstation), the script copies the same PDBs to
`--workdir/visualization_download/` instead and prints a short note.

## Requirements

- AmberTools on `PATH` providing **`ambpdb`** (or `bashambpdb`)
- Python can import **`parmed`** and **`numpy`**

## Example (mol6 pose2 monomer)

```bash
cd /path/to/urease_md

python3 scripts/extract_visualization_files.py \
  --prmtop systems/mol6_pose2_monomer/solvated_hmr.prmtop \
  --crd systems/mol6_pose2_monomer/solvated_hmr.rst7 \
  --out-prefix mol6_pose2_monomer \
  --ligand-resname LIG \
  --radius-a 5.0 \
  --workdir resp_results
```

Outputs (written in `--workdir` and copied to `/mnt/user-data/outputs/`):

- `mol6_pose2_monomer_full.pdb`
- `mol6_pose2_monomer_nowater.pdb`
- `mol6_pose2_monomer_watershell_5.0A.pdb`

## Notes

- The no-water PDB is produced by filtering the `ambpdb` output, but it does **not** rely on `grep`.
- The watershell is computed on coordinates by finding water residues whose oxygen (or any atom if no oxygen)
  comes within \(R\) Å of any ligand atom.

