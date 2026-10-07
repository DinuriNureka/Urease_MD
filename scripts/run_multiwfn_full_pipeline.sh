#!/usr/bin/env bash
# Run Multiwfn RESP pipeline: molden → Multiwfn RESP → apply charges → collect RMSE.
# Optional: ESP maps (step 3) via run_multiwfn_esp_map.sh separately.
#
# Usage (from project root):
#   export MULTIWFN=/path/to/Multiwfn
#   bash scripts/run_multiwfn_full_pipeline.sh
#
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${ROOT}"

bash "${SCRIPT_DIR}/gbw_to_molden.sh" "${ROOT}/resp_results"
bash "${SCRIPT_DIR}/run_multiwfn_resp.sh" "${ROOT}/resp_results"
bash "${SCRIPT_DIR}/apply_multiwfn_charges_batch.sh"
python3 "${SCRIPT_DIR}/collect_multiwfn_rmse.py" --results_dir "${ROOT}/resp_results" \
  --out_csv "${ROOT}/resp_results/multiwfn_resp_rmse_report.csv"
python3 "${SCRIPT_DIR}/verify_multiwfn_pipeline.py" --gaff2_dir "${ROOT}/gaff2_params" \
  --manifest "${ROOT}/resp_inputs/job_manifest.csv" \
  --report "${ROOT}/resp_results/multiwfn_resp_rmse_report.csv"
