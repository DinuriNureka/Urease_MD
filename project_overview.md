# Urease MD Simulation Project

## System

- **Protein:** H. pylori urease (PDB: 6ZJA), monomer heterodimer (one UreA + one UreB subunit)
- **Metal center:** 2 Ni atoms (dinuclear Ni2+ active site, 1 site)
- **Inhibitors:** 36 compounds in the screen — chemistry IDs are **mol1–mol33** plus **AHA**, **DMJ**, and **thiourea** (not a single 1–36 ligand index). Pipeline/ORCA folders use sequential batch IDs `mol01`–`mol36`; see `params/inhibitors/ligand_name_mapping.csv` (`compound` = canonical name). 78 pose combinations
- **Force fields:** ff19SB (protein), GAFF2 (inhibitors), MCPB.py bonded model (Ni center)
- **Water model:** TIP3P, truncated octahedral box, 12 A buffer
- **Salt:** 150 mM NaCl

## Hardware

- 56 CPU cores (28 physical, HT), 62 GB RAM
- GPUs: Quadro K2200 (4 GB, display), Tesla P100-PCIE-16GB (compute)
- ORCA 6.1.1 (CPU-only DFT + native GFN2-xTB), AMBER (pmemd.cuda), AmberTools

## Pipeline

| Phase | Description | Scripts | Status |
|-------|-------------|---------|--------|
| 0 | Project setup | `setup_project.sh`, `check_dependencies.sh` | Done |
| 1 | Structure preparation | `fix_protonation.py`, `verify_phase1.sh` | Done |
| 2 | Ni2+ parameterization | See sub-steps below | **In progress** |
| 3 | Inhibitor parameterization | `extract_ligands.py`, `batch_param_inhibitors.sh` | Not started |
| 4 | System assembly | `build_monomer_system.sh` | Not started |
| 5 | Equilibration | `run_equilibration.sh` + 7 input files | Not started |
| 6 | Production MD | `run_production.sh` + `production.in` | Not started |
| 7 | Analysis | `analyze_all.sh`, `run_mmpbsa.sh` | Not started |

### Phase 2 Sub-steps

| Step | Description | Status |
|------|-------------|--------|
| 2.1 | Extract QM model (71 atoms, chain B) | Done |
| 2.2 | Prepare NI1, NI2, KCX, OHB mol2 files | Done |
| 2.3 | Write MCPB.py input file | Done |
| 2.4 | Run MCPB.py Step 1 (generate QM templates) | Done |
| 2.5a | FM optimization (B3LYP/def2-SVP, CalcHess) | **Running** (cycle 51, CP-SCF iter 11) |
| 2.5b | BS-DFT frequency (Hessian for Seminario) | Queued (watchdog auto-launch) |
| 2.5c | BS-DFT ESP (CHELPG for RESP) | Queued (watchdog auto-launch) |
| 2.7b | Multi-frame RESP (GFN2-xTB MD + multi-ESP) | **xTB MD done** (20 ps, trajectory.xyz) |
| 2.6 | Convert ORCA Hessian to fchk | Pending |
| 2.7 | Multi-frame RESP charge fitting | Pending |
| 2.8 | Run MCPB.py Steps 2-4, copy final params | Pending |

### Resource Allocation

```
DFT FM optimization:     28 cores, ~18 GB RAM (running; xTB MD done)
Tesla P100:               idle (used in Phase 5 for pmemd.cuda)
```

## Key Files

```
config.env                              Central configuration
scripts/                                All executable scripts
params/ni_center/                       MCPB.py Ni2+ parameters
params/ni_center/orca_calcs/            ORCA DFT calculations
params/ni_center/orca_calcs/small_opt/  FM geometry optimization
params/ni_center/orca_calcs/small_freq/ BS-DFT frequency (queued)
params/ni_center/orca_calcs/large_mk/   BS-DFT ESP (queued)
params/ni_center/multiframe_resp/       Multi-frame RESP workflow
params/ni_center/multiframe_resp/xtb_md/  GFN2-xTB MD (done; trajectory.xyz)
params/inhibitors/                      GAFF2 inhibitor parameters
structures/                             PDB preparation pipeline
systems/                                Per-system assembled topologies
simulations/input_files/                AMBER input files (shared)
simulations/<name>/                     Per-system simulation output
analysis/                               cpptraj and MMPBSA results
```

## Automation

- `watch_and_launch.sh` — Watchdog: polls FM convergence every 5 min,
  auto-launches BS-DFT frequency then ESP sequentially.
- `launch_bs_jobs.sh` — Copies optimized geometry + orbitals, runs BS jobs.
- `run_multiframe_esp.sh` — Runs ORCA CHELPG on xTB MD frames.

## Execution Order

1. `bash scripts/setup_project.sh /path/to/urease_md`
2. `bash scripts/check_dependencies.sh`
3. Download PDB, run pdb4amber, reduce (manual steps in 00_MASTER_GUIDE.md)
4. `python3 scripts/fix_protonation.py ...` then `bash scripts/verify_phase1.sh`
5. `python3 scripts/extract_qm_model.py ...` then MCPB.py + ORCA workflow
6. `python3 scripts/extract_ligands.py ...` then `bash scripts/batch_param_inhibitors.sh`
7. `bash scripts/build_monomer_system.sh <name>` (x78)
8. `bash scripts/run_equilibration.sh <name>` (x78)
9. `bash scripts/run_production.sh <name> <segments>` (x78)
10. `bash scripts/analyze_all.sh` then `bash scripts/run_mmpbsa.sh <name>` (x78)

## Documentation

- `00_MASTER_GUIDE.md` — Phase 0–1 detailed instructions
- `02_mcpb_parameterization.md` — Phase 2 MCPB.py/ORCA workflow (single-frame)
- `02b_multiframe_resp.md` — Phase 2 multi-frame RESP extension
- `03_through_07_phases.md` — Phases 3–7 instructions
- `08_troubleshooting.md` — Symptom/cause/fix for all phases
- `params/ni_center/orca_calcs/README_bs_setup.md` — BS-DFT job setup notes
