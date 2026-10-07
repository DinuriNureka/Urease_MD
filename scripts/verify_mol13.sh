#!/usr/bin/env bash
# scripts/verify_mol13.sh — read-only mol13 verification (safe while equil_seg2 runs)
# Does not write under simulations/mol13_pose1_monomer/; temp files use /tmp/mol13_*
set -uo pipefail

_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.env
source "${_SCRIPT_DIR}/../config.env"

SYS="mol13_pose1_monomer"
PRMTOP="${PROJECT_ROOT}/systems/${SYS}/solvated_hmr.prmtop"
EQUIL2_OUT="${PROJECT_ROOT}/simulations/${SYS}/equil_seg2/equil_seg2.out"
PAT_NC="${PROJECT_ROOT}/simulations/${SYS}/crashtest/patched_test/patched_test.nc"
SEG1_NC="${PROJECT_ROOT}/simulations/${SYS}/equil_seg1/equil_seg1.nc"
RST_BRIDGE="${PROJECT_ROOT}/simulations/${SYS}/gpu_bridge2/gpu_bridge2.rst7"
NI_RST="${PROJECT_ROOT}/simulations/${SYS}/restraints/ni_ni.RST"
EQUIL_SEG_IN="${PROJECT_ROOT}/simulations/${SYS}/input_files/equil_seg.in"
EQUIL_CPU_IN="${PROJECT_ROOT}/simulations/input_files/equil_cpu_npt.in"
DISANG_LINK="${PROJECT_ROOT}/nmr_disang/mol13_pose1_monomer.RST"

CPPTRAJ_REQ="/opt/software/AMBER/pmemd24/bin/cpptraj"
if [[ -x "${CPPTRAJ_REQ}" ]]; then
  CPPTRAJ="${CPPTRAJ_REQ}"
else
  CPPTRAJ="$(command -v cpptraj 2>/dev/null || echo /usr/local/bin/cpptraj)"
fi

LOG="${PROJECT_ROOT}/logs/mol13_verification_$(date +%Y%m%d_%H%M%S).txt"

# Summary lines (filled by checks)
S1=S2=S3=S4=S5=S6=S7=S8=S9=S10=""
M2_MEAN="" M8_MEAN="" M9_MEAN="" M5_HMAX="" M10_GPU=""

FAIL_COUNT=0
warn_count=0

inc_fail() { FAIL_COUNT=$((FAIL_COUNT + 1)); }
inc_warn() { warn_count=$((warn_count + 1)); }

run_cpptraj() {
  "${CPPTRAJ}" "$@" 2>&1
}

main() {
  echo "========================================"
  echo "mol13_pose1_monomer verification"
  echo "Started: $(date -Is)"
  echo "CPPTRAJ: ${CPPTRAJ}"
  if [[ ! -x "${CPPTRAJ}" ]]; then
    echo "CHECK cpptraj binary: FAIL (not executable: ${CPPTRAJ})"
    inc_fail
    return 1
  fi
  if [[ "${CPPTRAJ}" != "${CPPTRAJ_REQ}" ]]; then
    echo "WARN: requested ${CPPTRAJ_REQ} missing — using ${CPPTRAJ}"
    inc_warn
  fi
  echo "========================================"

  # ----- CHECK 1: residue identity -----
  echo ""
  echo "=== CHECK 1 — Residue identity (numeric masks) ==="
  local c1_out="/tmp/mol13_check1_resinfo.txt"
  run_cpptraj -p "${PRMTOP}" <<'EOF' > "${c1_out}"
resinfo :136
resinfo :138
resinfo :219
resinfo :248
resinfo :274
resinfo :362
resinfo :570
resinfo :571
resinfo :572
resinfo :573
resinfo :NI
resinfo :KCX
resinfo :OHB
quit
EOF
  cat "${c1_out}"
  local c1=PASS
  # Parse key residues from resinfo blocks
  local n136 n138 n219 n248 n274 n362 n570 n571 n572
  # Data line after "#Res  Name First ..." is: "<resnum> <ResName> <First> ..." → name is $2
  n136=$(awk '/resinfo :136/{f=1} f&&/#Res/{getline; print $2; exit}' "${c1_out}")
  n138=$(awk '/resinfo :138/{f=1} f&&/#Res/{getline; print $2; exit}' "${c1_out}")
  n219=$(awk '/resinfo :219/{f=1} f&&/#Res/{getline; print $2; exit}' "${c1_out}")
  n248=$(awk '/resinfo :248/{f=1} f&&/#Res/{getline; print $2; exit}' "${c1_out}")
  n274=$(awk '/resinfo :274/{f=1} f&&/#Res/{getline; print $2; exit}' "${c1_out}")
  n362=$(awk '/resinfo :362/{f=1} f&&/#Res/{getline; print $2; exit}' "${c1_out}")
  n570=$(awk '/resinfo :570/{f=1} f&&/#Res/{getline; print $2; exit}' "${c1_out}")
  n571=$(awk '/resinfo :571/{f=1} f&&/#Res/{getline; print $2; exit}' "${c1_out}")
  n572=$(awk '/resinfo :572/{f=1} f&&/#Res/{getline; print $2; exit}' "${c1_out}")

  if awk '/resinfo :OHB/{f=1} f&&/#Res/{print; exit}' "${c1_out}" | grep -q '#Res'; then
    echo "CHECK 1: FAIL — :OHB residue block present (OHB should be absent)"
    c1=FAIL
    inc_fail
  else
    echo "CHECK 1: :OHB selection empty — OK"
  fi
  [[ "${n570}" == "NI" && "${n571}" == "NI" ]] || {
    echo "CHECK 1: FAIL — NI residues at :570/:571 (got ${n570} ${n571})"
    c1=FAIL
    inc_fail
  }
  [[ "${n572}" != "OHB" && -n "${n572}" ]] || {
    echo "CHECK 1: FAIL — :572 should not be OHB (got '${n572}')"
    c1=FAIL
    inc_fail
  }
  # His/Asp coordination sanity
  local his_asp_warn=0
  [[ "${n136}" =~ ^HI[ED]$ || "${n136}" == "HID" || "${n136}" == "HIE" || "${n136}" == "HIP" ]] || his_asp_warn=1
  [[ "${n138}" =~ ^HI[ED]$ || "${n138}" == "HID" || "${n138}" == "HIE" || "${n138}" == "HIP" ]] || his_asp_warn=1
  [[ "${n248}" =~ ^HI[ED]$ || "${n248}" == "HID" || "${n248}" == "HIE" || "${n248}" == "HIP" ]] || his_asp_warn=1
  [[ "${n274}" =~ ^HI[ED]$ || "${n274}" == "HID" || "${n274}" == "HIE" || "${n274}" == "HIP" ]] || his_asp_warn=1
  [[ "${n362}" == "ASP" || "${n362}" == "ASH" ]] || his_asp_warn=1
  [[ "${n219}" == "KCX" ]] || his_asp_warn=1
  if [[ "${his_asp_warn}" -eq 1 ]]; then
    echo "CHECK 1: WARN — coordinating residue names unexpected (136=${n136} 138=${n138} 219=${n219} 248=${n248} 274=${n274} 362=${n362})"
    [[ "${c1}" == "PASS" ]] && c1=WARN
    inc_warn
  fi
  if [[ "${c1}" == "PASS" ]]; then
    echo "CHECK 1: PASS — coordinating residues + NI + no OHB"
  elif [[ "${c1}" == "WARN" ]]; then
    echo "CHECK 1: WARN — see coordinating-residue note above"
  else
    echo "CHECK 1: FAIL — see messages above"
  fi
  S1="${c1}"

  # ----- CHECK 2: Ni–Ni trajectory (seg1 / patched_test) -----
  echo ""
  echo "=== CHECK 2 — Ni–Ni distance (seg1 trajectory) ==="
  local traj="" c2=FAIL
  if [[ -f "${PAT_NC}" ]]; then
    traj="${PAT_NC}"
  elif [[ -f "${SEG1_NC}" ]]; then
    traj="${SEG1_NC}"
    echo "NOTE: patched_test.nc missing — using equil_seg1.nc for Ni–Ni stats"
  fi
  if [[ -z "${traj}" ]]; then
    echo "CHECK 2: WARN — no patched_test.nc or equil_seg1.nc (cannot reproduce seg1 trajectory stats)"
    c2=WARN
    M2_MEAN="n/a"
    inc_warn
  else
    local dfile="/tmp/mol13_nini_seg1.dat"
    run_cpptraj -p "${PRMTOP}" <<EOF > /tmp/mol13_check2_cpptraj.log
trajin ${traj}
distance nini @8611 @8612 out ${dfile} noimage
run
quit
EOF
    M2_MEAN=$(awk 'NR>1 && $2!~/nan/i && $2+0==$2 {
      v=$2+0
      if(!n){min=max=v}
      if(v<min)min=v; if(v>max)max=v; sum+=v; n++
    } END {
      if(n>0){printf "%.4f", sum/n} else {print "n/a"}
    }' "${dfile}" 2>/dev/null || echo "n/a")
    read -r min_v max_v nfr <<< "$(awk 'NR>1 && $2!~/nan/i && $2+0==$2 {
      v=$2+0
      if(!n){min=max=v}
      if(v<min)min=v; if(v>max)max=v; n++
    } END { if(n>0) print min, max, n; else print "nan nan 0" }' "${dfile}")"
    echo "Ni–Ni from ${traj}: frames=${nfr} min=${min_v} max=${max_v} mean=${M2_MEAN} Å"
    if [[ "${nfr}" == "0" || "${M2_MEAN}" == "n/a" ]]; then
      echo "CHECK 2: WARN — no valid distance frames"
      c2=WARN
      inc_warn
    else
      c2=$(awk -v mean="${M2_MEAN}" -v min="${min_v}" -v max="${max_v}" 'BEGIN{
        m=mean+0; lo=min+0; hi=max+0;
        if(lo<2.5 || hi>4.5) {print "FAIL"; exit 0}
        if(m>=3.3 && m<=3.9) {print "PASS"; exit 0}
        if(m>=3.0 && m<=4.2) {print "WARN"; exit 0}
        print "FAIL"; exit 0
      }')
      if [[ "${c2}" == "PASS" ]]; then echo "CHECK 2: PASS"
      elif [[ "${c2}" == "WARN" ]]; then echo "CHECK 2: WARN"; inc_warn
      else echo "CHECK 2: FAIL"; inc_fail
      fi
    fi
  fi
  echo "CHECK 2: ${c2}"
  S2="${c2}"

  # ----- CHECK 3: NMR restraint in equil_seg2 -----
  echo ""
  echo "=== CHECK 3 — NMR / DISANG in equil_seg2.out ==="
  local c3=FAIL
  if [[ ! -f "${EQUIL2_OUT}" ]]; then
    echo "CHECK 3: FAIL — missing ${EQUIL2_OUT}"
    c3=FAIL
    inc_fail
  else
    echo "--- grep RESTRAINTS / DISANG / NMR (head) ---"
    grep -E "RESTRAINTS|DISANG|NMR|Begin reading energy" "${EQUIL2_OUT}" 2>/dev/null | head -25 || true
    if ! grep -q "DISANG" "${EQUIL2_OUT}"; then
      echo "CHECK 3: FAIL — DISANG not found in mdout header"
      c3=FAIL
      inc_fail
    else
      local nonzero
      nonzero=$(awk 'BEGIN{c=0} /NMR restraints: Bond/ {
        for(i=1;i<=NF;i++) if($i=="Bond") {v=$(i+2)+0; if(v>0.001) c++; next}
      } END{print c+0}' "${EQUIL2_OUT}")
      if [[ "${nonzero}" -gt 0 ]]; then
        echo "CHECK 3: PASS — DISANG parsed and non-zero NMR Bond energy observed"
        c3=PASS
      else
        echo "CHECK 3: WARN — DISANG present but all sampled NMR Bond energies ~0 (flat region or inactive)"
        c3=WARN
        inc_warn
      fi
    fi
  fi
  S3="${c3}"

  # ----- CHECK 4: equil_cpu_npt mask vs residue 572 -----
  echo ""
  echo "=== CHECK 4 — equil_cpu_npt restraint mask & :572 identity ==="
  run_cpptraj -p "${PRMTOP}" <<'EOF' 2>&1 | tee /tmp/mol13_check4_res572.txt
resinfo :572
quit
EOF
  echo "--- equil_cpu_npt.in restraint lines ---"
  grep -E "restraintmask|restraint_wt|restraint_mask" "${EQUIL_CPU_IN}" 2>/dev/null || echo "(no restraint lines found)"
  local r572
  r572=$(awk '/resinfo :572/{f=1} f&&/#Res/{getline; print $2; exit}' /tmp/mol13_check4_res572.txt)
  echo "Residue :572 = ${r572}"
  local c4="INFO"
  if [[ "${r572}" == "OHB" ]]; then
    echo "CHECK 4: FAIL — :572 is still OHB"
    c4=FAIL
    inc_fail
  else
    echo "CHECK 4: INFO — :572 is ${r572} (not OHB); historical :572 OHB mask issue does not apply to this topology"
  fi
  local maskline
  maskline=$(grep -E "restraintmask|restraint_mask" "${EQUIL_CPU_IN}" 2>/dev/null | head -1 || true)
  if echo "${maskline}" | grep -q ':572'; then
    echo "CHECK 4: NOTE — restraint mask still references :572 (stale). Suggested fix (do not auto-apply):"
    echo "  sed -i \"s/:572//g\" simulations/input_files/equil_cpu_npt.in   # then review commas"
  elif echo "${maskline}" | grep -q ":NI,KCX'"; then
    echo "CHECK 4: NOTE — mask includes :NI,KCX pattern (plus other residues per line)"
  fi
  if ! echo "${maskline}" | grep -q "restraintmask=':NI,KCX'"; then
    echo "CHECK 4: NOTE — restraintmask is not exactly ':NI,KCX' only. Optional tighten:"
    echo "  sed -i \"s/restraintmask=.*$/restraintmask=':NI,KCX',/\" simulations/input_files/equil_cpu_npt.in"
  fi
  S4="${c4}"

  # ----- CHECK 5: HMR -----
  echo ""
  echo "=== CHECK 5 — HMR (H mass) ==="
  local c5=FAIL
  run_cpptraj -p "${PRMTOP}" <<'EOF' >/tmp/mol13_check5_atominfo.log 2>&1
atominfo @H= &!:WAT out /tmp/mol13_atominfo_h.txt
quit
EOF
  # atominfo columns: #Atom Name #Res Name ... Mass is typically field 8 on data lines
  M5_HMAX=$(awk '/^[[:space:]]*#Atom/ {next} /^[[:space:]]*[0-9]+/ && NF>=8 {
      v=$8+0; if(v>m)m=v
    } END{printf "%.4f", m+0}' /tmp/mol13_atominfo_h.txt)
  echo "Max protein/solute H mass (non-water): ${M5_HMAX} Da"
  awk -v m="${M5_HMAX}" 'BEGIN{if(m+0>2.0) exit 0; exit 1}' && { c5=PASS; echo "CHECK 5: PASS"; } || { c5=FAIL; echo "CHECK 5: FAIL — HMR not detected (expect max H > 2.0 Da)"; inc_fail; }
  S5="${c5}"

  # ----- CHECK 6: box type -----
  echo ""
  echo "=== CHECK 6 — Box type (rectangular vs octahedron) ==="
  local c6=FAIL
  local bfile="/tmp/mol13_box_cpptraj.log"
  if [[ -f "${RST_BRIDGE}" ]]; then
    run_cpptraj -p "${PRMTOP}" <<EOF > "${bfile}" 2>&1
trajin ${RST_BRIDGE}
box
run
quit
EOF
  else
    run_cpptraj -p "${PRMTOP}" <<'EOF' > "${bfile}" 2>&1
box
quit
EOF
  fi
  cat "${bfile}"
  if grep -qi 'truncat.*octa\|IFBOX.*2\|Octahedron' "${bfile}"; then
    echo "CHECK 6: FAIL — truncated octahedron / IFBOX 2 suspected"
    c6=FAIL
    inc_fail
  elif grep -qi 'Orthorhombic\|IFBOX.*1' "${bfile}" || grep -qi 'Orthorhombic' "${bfile}"; then
    echo "CHECK 6: PASS — orthorhombic / rectangular box"
    c6=PASS
  else
    echo "CHECK 6: WARN — could not classify box from cpptraj log"
    c6=WARN
    inc_warn
  fi
  S6="${c6}"

  # ----- CHECK 7: DISANG symlink + RST + mdin -----
  echo ""
  echo "=== CHECK 7 — DISANG symlink + ni_ni.RST + equil_seg.in ==="
  local c7=PASS
  ls -la "${DISANG_LINK}" 2>&1 || true
  local resolved
  resolved=$(readlink -f "${DISANG_LINK}" 2>/dev/null || true)
  echo "readlink -f: ${resolved}"
  echo "--- simulations/${SYS}/restraints/ni_ni.RST ---"
  cat "${NI_RST}" 2>/dev/null || echo "(missing ${NI_RST})"
  echo "--- DISANG in equil_seg.in ---"
  grep "DISANG" "${EQUIL_SEG_IN}" 2>/dev/null || echo "(missing equil_seg.in)"
  local disang_line disang_val
  disang_line=$(grep '^DISANG=' "${EQUIL_SEG_IN}" 2>/dev/null | head -1 || true)
  disang_val="${disang_line#DISANG=}"
  local len=${#disang_val}
  if [[ ! -L "${DISANG_LINK}" ]]; then
    echo "CHECK 7: FAIL — nmr_disang symlink missing"
    c7=FAIL
    inc_fail
  elif [[ "${resolved}" != "${NI_RST}" ]]; then
    echo "CHECK 7: FAIL — symlink does not resolve to ${NI_RST}"
    c7=FAIL
    inc_fail
  elif ! grep -qE 'iat=[0-9]+,[0-9]+,' "${NI_RST}"; then
    echo "CHECK 7: FAIL — ni_ni.RST iat= line missing or invalid"
    c7=FAIL
    inc_fail
  elif [[ "${len}" -ge 80 ]]; then
    echo "CHECK 7: FAIL — DISANG path length ${len} (>= 80 Fortran limit)"
    c7=FAIL
    inc_fail
  elif ! echo "${disang_val}" | grep -q 'nmr_disang/mol13_pose1_monomer.RST'; then
    echo "CHECK 7: WARN — DISANG= does not use nmr_disang/mol13_pose1_monomer.RST substring"
    c7=WARN
    inc_warn
  else
    echo "CHECK 7: PASS — symlink, iat=, path length OK"
  fi
  S7="${c7}"

  # ----- CHECK 8: density -----
  echo ""
  echo "=== CHECK 8 — Density (equil_seg2, live) ==="
  local c8=PASS
  if [[ ! -f "${EQUIL2_OUT}" ]]; then
    c8=FAIL
    M8_MEAN="n/a"
    echo "CHECK 8: FAIL — no equil_seg2.out"
    inc_fail
  else
    M8_MEAN=$(grep 'Density' "${EQUIL2_OUT}" | grep -v '^|' | awk '{print $NF}' | awk 'NF && $1+0==$1 {
      v=$1+0
      if(!n){min=max=v}
      if(v<min)min=v; if(v>max)max=v; sum+=v; n++
    } END {
      if(n>0) printf "%.4f", sum/n; else print "n/a"
    }')
    local n8
    n8=$(grep 'Density' "${EQUIL2_OUT}" | grep -v '^|' | awk '{print $NF}' | awk 'NF && $1+0==$1 {c++} END{print c+0}')
    echo "Density stats: lines=${n8} mean=${M8_MEAN} g/cm³"
    if [[ "${M8_MEAN}" == "n/a" || "${n8}" == "0" ]]; then
      echo "CHECK 8: WARN — no density lines yet (run just started?)"
      c8=WARN
      M8_MEAN="n/a"
      inc_warn
    else
      local e8
      # Few density lines = early live run — relax upper bound slightly (still flag severe drift)
      if [[ "${n8}" -lt 500 ]]; then
        e8=$(awk -v m="${M8_MEAN}" 'BEGIN{x=m+0; if(x>=0.990&&x<=1.010) print 0; else if(x>=0.970&&x<=1.030) print 1; else print 2}')
      else
        e8=$(awk -v m="${M8_MEAN}" 'BEGIN{x=m+0; if(x>=0.990&&x<=1.010) print 0; else if(x>=0.980&&x<=1.020) print 1; else print 2}')
      fi
      if [[ "${e8}" == "0" ]]; then echo "CHECK 8: PASS"
      elif [[ "${e8}" == "1" ]]; then
        echo "CHECK 8: WARN — density mean ${M8_MEAN} (relaxed/early-run band; target 0.990–1.010 long-run)"
        c8=WARN
        inc_warn
      else echo "CHECK 8: FAIL — density mean ${M8_MEAN} outside acceptable range"; c8=FAIL; inc_fail
      fi
    fi
  fi
  S8="${c8}"

  # ----- CHECK 9: temperature -----
  echo ""
  echo "=== CHECK 9 — Temperature (equil_seg2) ==="
  local c9=PASS
  if [[ ! -f "${EQUIL2_OUT}" ]]; then
    c9=FAIL
    M9_MEAN="n/a"
    inc_fail
  else
    M9_MEAN=$(grep -E 'TEMP\(K\)' "${EQUIL2_OUT}" | sed -n 's/.*TEMP(K)[[:space:]]*=[[:space:]]*\([0-9.]*\).*/\1/p' | awk 'NF && $1+0==$1 {
      v=$1+0
      if(!n){min=max=v}
      if(v<min)min=v; if(v>max)max=v; sum+=v; n++
    } END {
      if(n>0) printf "%.2f", sum/n; else print "n/a"
    }')
    read -r tmin tmax tn <<< "$(grep -E 'TEMP\(K\)' "${EQUIL2_OUT}" | sed -n 's/.*TEMP(K)[[:space:]]*=[[:space:]]*\([0-9.]*\).*/\1/p' | awk 'NF && $1+0==$1 {
      v=$1+0
      if(!n){min=max=v}
      if(v<min)min=v; if(v>max)max=v; n++
    } END { if(n>0) print min, max, n; else print "nan nan 0" }')"
    echo "Temperature: frames=${tn} min=${tmin} max=${tmax} mean=${M9_MEAN} K"
    if [[ "${tn}" == "0" || "${M9_MEAN}" == "n/a" ]]; then
      c9=WARN
      inc_warn
    else
      local rng
      rng=$(awk -v a="${tmin}" -v b="${tmax}" 'BEGIN{print b-a}')
      local e9
      e9=$(awk -v m="${M9_MEAN}" -v r="${rng}" 'BEGIN{
        mm=m+0; rr=r+0
        if(mm>=295&&mm<=305&&rr<20) print 0
        else if(mm>=290&&mm<=310) print 1
        else print 2
      }')
      if [[ "${e9}" == "0" ]]; then echo "CHECK 9: PASS"
      elif [[ "${e9}" == "1" ]]; then echo "CHECK 9: WARN"; c9=WARN; inc_warn
      else echo "CHECK 9: FAIL"; c9=FAIL; inc_fail
      fi
    fi
  fi
  S9="${c9}"

  # ----- CHECK 10: GPU -----
  echo ""
  echo "=== CHECK 10 — GPU assignment (P100 vs K2200) ==="
  local c10=FAIL
  PMEMD_PID=""
  for pid in $(pgrep -f 'pmemd\.cuda' 2>/dev/null || true); do
    if [[ -r "/proc/${pid}/cmdline" ]] && tr '\0' ' ' < "/proc/${pid}/cmdline" | grep -q 'mol13_pose1_monomer'; then
      PMEMD_PID="${pid}"
      break
    fi
  done
  if [[ -z "${PMEMD_PID}" ]]; then
    PMEMD_PID=$(pgrep -f 'pmemd\.cuda' 2>/dev/null | head -1 || true)
  fi
  echo "Selected pmemd PID: ${PMEMD_PID:-none}"
  if [[ -n "${PMEMD_PID}" ]] && [[ -r "/proc/${PMEMD_PID}/exe" ]]; then
    readlink "/proc/${PMEMD_PID}/exe" 2>/dev/null || true
    echo "--- CUDA_* from process environ ---"
    tr '\0' '\n' < "/proc/${PMEMD_PID}/environ" 2>/dev/null | grep -E '^CUDA' || true
  fi
  echo "--- nvidia-smi (GPUs) ---"
  nvidia-smi --query-gpu=index,name,uuid --format=csv,noheader 2>/dev/null || true
  echo "--- nvidia-smi (compute apps) ---"
  nvidia-smi --query-compute-apps=pid,gpu_uuid,used_memory --format=csv,noheader 2>/dev/null || true

  local cuda_vis=""
  if [[ -n "${PMEMD_PID}" ]]; then
    cuda_vis=$(tr '\0' '\n' < "/proc/${PMEMD_PID}/environ" 2>/dev/null | grep '^CUDA_VISIBLE_DEVICES=' | cut -d= -f2- || true)
  fi
  M10_GPU="${cuda_vis:-?}"
  # P100 UUID substring (common); K2200 often "Graphics"
  local p100_uuid k2200_uuid
  k2200_uuid=$(nvidia-smi --query-gpu=index,uuid --format=csv,noheader 2>/dev/null | awk -F',' '$1+0==0 {gsub(/^ +| +$/,"",$2); print $2; exit}')
  p100_uuid=$(nvidia-smi --query-gpu=index,uuid --format=csv,noheader 2>/dev/null | awk -F',' '$1+0==1 {gsub(/^ +| +$/,"",$2); print $2; exit}')
  local app_uuid=""
  if [[ -n "${PMEMD_PID}" ]]; then
    app_uuid=$(nvidia-smi --query-compute-apps=pid,gpu_uuid --format=csv,noheader 2>/dev/null | awk -F',' -v p="${PMEMD_PID}" '
      { gsub(/^ +| +$/,"",$1); if ($1 == p) { gsub(/^ +| +$/,"",$2); print $2; exit } }')
  fi
  if [[ -z "${PMEMD_PID}" ]]; then
    echo "CHECK 10: WARN — no pmemd.cuda process found (cannot verify GPU binding live)"
    c10=WARN
    inc_warn
  elif [[ "${cuda_vis}" == "1" ]] && [[ -n "${app_uuid}" ]] && [[ -n "${p100_uuid}" ]] && [[ "${app_uuid}" == "${p100_uuid}" ]]; then
    c10=PASS
    echo "CHECK 10: PASS — CUDA_VISIBLE_DEVICES=1 and process maps to GPU index 1 UUID"
  elif [[ "${cuda_vis}" == "1" ]]; then
    echo "CHECK 10: WARN — CUDA_VISIBLE_DEVICES=1 but PID→GPU UUID mismatch or parse issue (app_uuid=${app_uuid} gpu1=${p100_uuid})"
    c10=WARN
    inc_warn
  elif [[ "${cuda_vis}" == "0" ]]; then
    echo "CHECK 10: FAIL — CUDA_VISIBLE_DEVICES=0 (expected 1 for P100)"
    c10=FAIL
    inc_fail
  else
    echo "CHECK 10: FAIL — CUDA_VISIBLE_DEVICES is '${cuda_vis}' (expected 1)"
    c10=FAIL
    inc_fail
  fi
  S10="${c10}"

  # ----- FINAL SUMMARY -----
  echo ""
  echo "========================================"
  echo "mol13 verification summary"
  echo "========================================"
  printf "CHECK 1  Residue identity:        %s\n" "${S1}"
  printf "CHECK 2  Ni-Ni distance seg1:     %s  (mean=%s Å)\n" "${S2}" "${M2_MEAN:-n/a}"
  printf "CHECK 3  NMR restraint active:    %s\n" "${S3}"
  printf "CHECK 4  equil_cpu_npt mask:      %s (residue 572 = %s)\n" "${S4}" "${r572:-?}"
  printf "CHECK 5  HMR applied:             %s  (max H mass = %s Da)\n" "${S5}" "${M5_HMAX}"
  printf "CHECK 6  Rectangular box:         %s\n" "${S6}"
  printf "CHECK 7  DISANG symlink:          %s\n" "${S7}"
  printf "CHECK 8  Density:                 %s  (mean=%s g/cm³)\n" "${S8}" "${M8_MEAN:-n/a}"
  printf "CHECK 9  Temperature:             %s  (mean=%s K)\n" "${S9}" "${M9_MEAN:-n/a}"
  printf "CHECK 10 GPU assignment:          %s  (CUDA_VISIBLE_DEVICES=%s)\n" "${S10}" "${M10_GPU}"
  echo "========================================"
  if [[ "${FAIL_COUNT}" -eq 0 ]]; then
    echo "OVERALL: PASS (no FAIL checks; WARNs=${warn_count})"
  else
    echo "OVERALL: FAIL (${FAIL_COUNT} check(s) failed; WARNs=${warn_count})"
  fi
  echo "Full report: ${LOG}"
  echo "========================================"
  [[ "${FAIL_COUNT}" -eq 0 ]] && return 0
  return 1
}

set -o pipefail
main 2>&1 | tee "${LOG}"
exit "${PIPESTATUS[0]}"
