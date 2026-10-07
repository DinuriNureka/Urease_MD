# Phase 3 post-build QC

Script: `scripts/phase3_post_build_qc.py`

## Ligand placement (critical)

1. **Centroid test** — Distance between the ligand center of mass (all ligand atoms) and the center of mass of protein Cα atoms (excluding water, ions, and the ligand). **Fail** if distance &gt; 30 Å.
2. **Contacts** — `cpptraj` `nativecontacts :<LIGRES> !(:WAT,HOH,TIP3,TIP4,Na+,Cl-,NA,CL,:<LIGRES>) distance 5.0 writecontacts …`. Count unique contacting residues that are not water, ions, ligand, or NI. **Fail** if fewer than 5.

Ligand residue name is inferred from the topology (`LIG` or `MOL` preferred, else first nonstandard non-solvent residue).

## Informational (not pass/fail)

- **min(ligand–Ni)** — Minimum distance from any ligand atom to any nickel (atomic number 28), from loaded coordinates.

## Geometry

- **Ni–Ni** — From coordinates: among all Ni atoms, prefer a pair in 3.0–4.5 Å (dinuclear); else shortest Ni–Ni distance. **Fail** if outside 3.0–4.5 Å when a value is available.

## Other checks

Unchanged: total charge, solute integer charge, HMR mean hydrogen mass, `solvateoct` buffer in tleap input, atom-count median band, tleap log keyword scan, `run_full_pipeline.sh` template keywords.

## Run

```bash
cd /path/to/urease_md
python3 scripts/phase3_post_build_qc.py | tee logs/phase3_post_build_qc_latest.log
```

`CPPTRAJ` may be set if `cpptraj` is not on `PATH` (e.g. `export CPPTRAJ=/usr/local/bin/cpptraj`).
