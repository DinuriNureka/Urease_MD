#!/usr/bin/env bash
# Run MCPB.py steps 2-4 after ORCA completes
# Converts ORCA outputs to Gaussian format, runs Seminario + charge fitting
set -euo pipefail
cd /path/to/urease_md/params/ni_center_no_ohb/mcpb_input

MCPB="${MCPB:-/usr/local/bin/MCPB.py}"

echo "=== MCPB.py Steps 2-4 ==="
echo "Step 2 (Seminario force constants)..."
"${MCPB}" -i urease_no_ohb.in -s 2
echo "Step 3 (charge fitting)..."
"${MCPB}" -i urease_no_ohb.in -s 3
echo "Step 4 (frcmod generation)..."
"${MCPB}" -i urease_no_ohb.in -s 4
echo "=== Done: $(date) ==="
ls -la *.frcmod *.mol2 2>/dev/null || true
