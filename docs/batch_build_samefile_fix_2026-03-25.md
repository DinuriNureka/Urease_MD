## Batch build crash fix (SameFileError)

### Problem

`scripts/build_full_systems_batch.py` could crash with:

- `shutil.SameFileError: ... 'systems/mol6_pose2_full/repair_molecule_table.py' and ... 'systems/mol6_pose2_full/repair_molecule_table.py'`

This happens when the template system directory (`TEMPLATE_SYS`, currently `systems/mol6_pose2_full`) is also one of the systems being built, so the script tries to copy a file **onto itself**.

### Fix

Updated the script to **skip copying** when source and destination resolve to the same path:

- `copy_if_missing()`: now checks `src.resolve() == dst.resolve()` and returns early.
- `repair_molecule_table()`: now only copies `repair_molecule_table.py` if `script_src.resolve() != script_dst.resolve()`.

### Restart batch build (20 threads cap)

Run with:

```bash
MAX_CPU_THREADS=20 python3 scripts/build_full_systems_batch.py |& tee "scripts/batch_build_restart3_$(date +%Y%m%d_%H%M%S).log"
```

