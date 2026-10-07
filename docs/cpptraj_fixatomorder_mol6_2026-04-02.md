# cpptraj `fixatomorder` for mol6 test topology (`test_mol6`)

## Context

`test_mol6.prmtop` reported one molecule with **non-contiguous atom segments** (8610 atoms + 4 atoms at the end). `fixatomorder` reorders atoms so molecules are contiguous and writes a matching topology.

## Commands

**Step 1 — reorder topology and coordinates**

```bash
cpptraj << EOF
parm /tmp/test_mol6.prmtop
trajin /tmp/test_mol6.inpcrd
fixatomorder outprefix /tmp/fixed
run
EOF
```

- **`parm`**: loads the Amber topology.
- **`trajin`**: loads coordinates (restart).
- **`fixatomorder outprefix /tmp/fixed`**: writes reordered topology as `/tmp/fixed.test_mol6.prmtop` (prefix + original parm basename).

**Optional:** write a reordered restart for use with the fixed prmtop:

```bash
cpptraj << EOF
parm /tmp/test_mol6.prmtop
trajin /tmp/test_mol6.inpcrd
fixatomorder outprefix /tmp/fixed
trajout /tmp/fixed.rst7 restart
run
EOF
```

**Step 2 — Ni / OHB atom numbering on fixed topology**

```bash
cpptraj -p /tmp/fixed.test_mol6.prmtop << EOF
atominfo :NI
atominfo :OHB
EOF
```

## Results (2026-04-02 run)

| Role   | Atom # | Residue # | Residue name |
|--------|--------|-----------|--------------|
| Ni (M1)| 8611   | 570       | NI           |
| Ni (M2)| 8612   | 571       | NI           |
| O (bridge)| 8613 | 572    | OHB          |
| H      | 8614   | 572       | OHB          |

**Distances** (same frame, `fixatomorder` then `@8611`–`@8614`):

| Quantity | Å     |
|----------|-------|
| Ni–Ni    | 3.5105 |
| Ni1–O    | 2.0500 |
| Ni2–O    | 2.0503 |

These are chemically reasonable for a bridging hydroxide in a dinuclear Ni site.

## Files produced

- `/tmp/fixed.test_mol6.prmtop` — reordered topology
- `/tmp/fixed.rst7` — optional reordered restart (if `trajout` used)

Use the **fixed** prmtop with the **reordered** restart for any downstream tool that assumes contiguous molecules.
