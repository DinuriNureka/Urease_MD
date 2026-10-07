# mol6 monomer: tLeap rebuild (9 Å orthorhombic) + `fixatomorder` + install

**Date:** 2026-04-02

## What was done

1. **tLeap** rebuilt `mol` from `monomer_BQ_apo.pdb` + mol6 ligand, Ni patches, 9 Å `solvatebox`, neutralizing ions → `/tmp/rebuild_mol6.{prmtop,inpcrd}`.
2. **cpptraj** `fixatomorder` (non-contiguous molecule: segments 1–8610 and 12391–12394) → `/tmp/rebuild_fixed.rebuild_mol6.prmtop`.
3. **Reordered coordinates** written with **`trajout ... restart`** to `/tmp/rebuild_fixed.rebuild_mol6.inpcrd` (required so atom order matches the fixed topology).
4. **Installed** over `systems/mol6_pose2_monomer/monomer_mol6.{prmtop,inpcrd}`; backups: `*.BAD_atomorder`.
5. **Removed** old equilibration trees under `simulations/mol6_pose2_monomer/` (`min1`, `min2`, `heat`, `heat_nvt_cpu`, `equil_cpu`, `equil_cpu_npt`, `gpu_bridge`, `equil1`–`equil3`, `cpu_fix`, `npt_ramp`).

## Geometry check (installed files)

| Quantity | Å (frame 1) |
|----------|----------------|
| Ni–Ni (`:570@NI`–`:571@NI`) | 3.5105 |
| Ni1–O (`:570@NI`–`:572@O`)   | 2.0500 |
| Ni2–O (`:571@NI`–`:572@O`)   | 2.0503 |

Atom IDs on fixed topology: Ni **8611**, **8612**; OHB O **8613**, H **8614** (residues 570–572).

## Notes for copy-paste scripts

- **`fixatomorder` alone does not write a matching `.inpcrd`.** Add e.g.  
  `trajout /tmp/rebuild_fixed.rebuild_mol6.inpcrd restart` before `run`.
- This **cpptraj** build (v5.1.0) has no `printdata`; use **`writedata <file> <sets...>`** after `run`.
- **tLeap** reported non-integral / non-zero net charge (~0.017 e) before `saveamberparm`; same as prior builds unless you re-neutralize explicitly.

## Ready next

Restart equilibration from **min1** using the updated `monomer_mol6` files.
