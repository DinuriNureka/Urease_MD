#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../config.env"

ORCA="${ORCA_PATH}"
PROJECT="${PROJECT_ROOT}"
FRAME_DIR="${PROJECT}/params/ni_center/multiframe_resp/frames"
ESP_DIR="${PROJECT}/params/ni_center/multiframe_resp/esp_calcs"
BASE_ESP="${PROJECT}/params/ni_center/orca_calcs/large_mk/large_mk.out"

mkdir -p "${ESP_DIR}"

# Frame 00 comes from the already completed BS-DFT ESP.
if [[ -f "${BASE_ESP}" ]]; then
  mkdir -p "${ESP_DIR}/frame_00"
  cp "${BASE_ESP}" "${ESP_DIR}/frame_00/esp.out"
  echo "[OK] frame_00 copied from completed large_mk.out"
fi

for xyz in "${FRAME_DIR}"/frame_*.xyz; do
  [[ -e "${xyz}" ]] || { echo "No frames found in ${FRAME_DIR}"; exit 1; }
  name="$(basename "${xyz}" .xyz)"
  idx="${name#frame_}"
  work="${ESP_DIR}/${name}"

  if [[ "${idx}" == "00" && -f "${work}/esp.out" ]]; then
    echo "[SKIP] ${name} already present"
    continue
  fi

  mkdir -p "${work}"
  cp "${xyz}" "${work}/esp.xyz"

  cat > "${work}/esp.inp" <<'EOF'
! B3LYP def2-SVP def2/J RIJCOSX TightSCF D3BJ DEFGRID2
! SlowConv KeepDens CHELPG UKS

%maxcore 1200

%pal
  nprocs 28
end

%scf
  MaxIter 500
  BrokenSym 2,2
  AutoTRAH true
  SOSCFStart 0.01
  SmearTemp 5000
  HFTyp UHF
  DIISMaxEq 15
end

* xyzfile 1 1 esp.xyz
EOF

  echo "[RUN] ${name}"
  (cd "${work}" && "${ORCA}" esp.inp > esp.out 2>&1)
  if grep -q "ORCA TERMINATED NORMALLY" "${work}/esp.out"; then
    echo "[OK] ${name} done"
  else
    echo "[FAIL] ${name} did not terminate normally"
    exit 2
  fi
done

echo "All frame ESP jobs completed."

