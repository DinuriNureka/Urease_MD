#!/usr/bin/env bash
# Rebuild all 36 monomer systems WITHOUT OHB
# Uses new no-OHB frcmod + charges from MCPB.py
# Applies fixatomorder after each tLeap build
set -euo pipefail
cd /path/to/urease_md
# shellcheck source=/dev/null
source config.env

echo "=== Rebuild all 36 systems without OHB ==="
echo "Started: $(date)"

# TODO: Update paths after MCPB.py step 4 produces:
#   - new frcmod (no Y4 terms, proper Ni force constants)
#   - new mol2 files for NI residues
#   - new set_patch charges (no OHB entries)

# For each system:
# 1. Generate tleap input WITHOUT OHB mol2, bonds, or patches
# 2. Run tleap
# 3. Apply fixatomorder
# 4. Apply HMR via parmed
# 5. Verify Ni-Ni = 3.40-3.65 Å

echo "Script skeleton ready — fill in after MCPB.py completes"
