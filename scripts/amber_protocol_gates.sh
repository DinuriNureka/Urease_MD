#!/usr/bin/env bash
# amber_protocol_gates.sh — enforce MD protocol before launch.
# Source from run_equilibration.sh / run_production.sh after config.env.
# GPU dynamics binary: set PMEMD_CUDA_BIN in config.env (canonical CUDA 12.8 build); do not rely on $AMBER_HOME/bin/pmemd.cuda.
#
# Protocol (post-Apr 8 2026 upgrade):
# - Bussi thermostat (ntt=11, tautp=2.0) throughout all dynamics stages
# - MC barostat (barostat=2, mcbarint=25) for short NPT stages (equil1–3)
# - MC barostat (barostat=2) for equil_free and production (resilient queue uses mcbarint=100 there)
# - Ni-Ni flat-bottom NMR restraint (nmropt=1) in all dynamics stages
# - HMR topology (solvated_hmr.prmtop) for all stages (CPU + GPU)
# - dt=0.002 for equilibration, dt=0.004 for equil_free and production
# - Production segment: nstlim=2,500,000 * dt=0.004 ps = 10 ns per segment (2 segments × 3 reps = 60 ns)
#
# Exits the parent shell with status 1 on failure (callers should use set -e).

SCRIPT_DIR_AMBER_GATE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VALIDATE_CNTRL="${SCRIPT_DIR_AMBER_GATE}/validate_amber_cntrl.py"

_require_validate_cntrl() {
  [[ -x "${VALIDATE_CNTRL}" ]] || chmod +x "${VALIDATE_CNTRL}" 2>/dev/null || true
  [[ -f "${VALIDATE_CNTRL}" ]] || {
    echo "[GATE FAIL] Missing ${VALIDATE_CNTRL}" >&2
    exit 1
  }
}

# Check that a dynamics input file has nmropt=1 and DISANG line so the Ni-Ni
# NMR restraint is actually loaded. Without this check an accidental removal
# would silently produce unrestrained Ni-Ni dynamics (the mol13 Berendsen bug).
_require_nmropt() {
  local p="$1"
  if ! grep -qE '^\s*nmropt\s*=\s*1\b' "${p}"; then
    echo "[GATE FAIL] ${p}: missing nmropt=1 — Ni-Ni restraint will not load." >&2
    exit 1
  fi
  if ! grep -q "DISANG=" "${p}"; then
    echo "[GATE FAIL] ${p}: missing DISANG= line for Ni-Ni restraint file." >&2
    exit 1
  fi
}

# Short CPU stages: Bussi thermostat (ntt=11), no barostat (NVT heat_nvt_cpu)
# or MC barostat (equil_cpu_npt). Both use dt=0.002 and require nmropt=1.
gate_cpu_md_inputs() {
  local inp_dir="$1"
  local f p
  _require_validate_cntrl

  # heat_nvt_cpu: NVT, no barostat
  p="${inp_dir}/heat_nvt_cpu.in"
  [[ -f "${p}" ]] || {
    echo "[GATE FAIL] Missing CPU staging input: ${p}" >&2
    exit 1
  }
  python3 "${VALIDATE_CNTRL}" "${p}" dt=0.002 ntc=2 ntf=2 ntt=11 tautp=2.0 || exit 1
  _require_nmropt "${p}"

  # equil_cpu_npt: NPT with MC barostat
  p="${inp_dir}/equil_cpu_npt.in"
  [[ -f "${p}" ]] || {
    echo "[GATE FAIL] Missing CPU staging input: ${p}" >&2
    exit 1
  }
  python3 "${VALIDATE_CNTRL}" "${p}" dt=0.002 ntc=2 ntf=2 ntt=11 tautp=2.0 barostat=2 mcbarint=25 || exit 1
  _require_nmropt "${p}"
}

# Main equilibration stages: Bussi (ntt=11) throughout. MC barostat in NPT stages.
gate_equilibration_md_inputs() {
  local inp_dir="$1"
  local f p
  _require_validate_cntrl
  gate_cpu_md_inputs "${inp_dir}"

  # heat and gpu_bridge: NVT, dt=0.002, Bussi, no barostat
  for f in heat.in gpu_bridge.in; do
    p="${inp_dir}/${f}"
    [[ -f "${p}" ]] || {
      echo "[GATE FAIL] Missing equilibration input: ${p}" >&2
      exit 1
    }
    python3 "${VALIDATE_CNTRL}" "${p}" dt=0.002 ntc=2 ntf=2 ntt=11 tautp=2.0 || exit 1
    _require_nmropt "${p}"
  done

  # equil1-3: NPT, dt=0.002, Bussi, MC barostat
  for f in equil1.in equil2.in equil3.in; do
    p="${inp_dir}/${f}"
    [[ -f "${p}" ]] || {
      echo "[GATE FAIL] Missing equilibration input: ${p}" >&2
      exit 1
    }
    python3 "${VALIDATE_CNTRL}" "${p}" dt=0.002 ntc=2 ntf=2 ntt=11 tautp=2.0 barostat=2 mcbarint=25 || exit 1
    _require_nmropt "${p}"
  done

  # equil_free: NPT, dt=0.004 (HMR), Bussi, MC barostat (mcbarint=100 matches resilient queue / equil_free.in)
  p="${inp_dir}/equil_free.in"
  [[ -f "${p}" ]] || { echo "[GATE FAIL] Missing equilibration input: ${p}" >&2; exit 1; }
  python3 "${VALIDATE_CNTRL}" "${p}" dt=0.004 ntc=2 ntf=2 ntt=11 tautp=2.0 barostat=2 mcbarint=100 || exit 1
  _require_nmropt "${p}"
}

# Deprecated: kept as no-op for compatibility with any callers that still reference it.
# Pre-Apr-8-2026 protocol required non-HMR prmtop for equilibration; current protocol
# uses solvated_hmr.prmtop for all stages (HMR is safe at dt=0.002 for equilibration).
gate_equil_topology_monomer() {
  return 0
}

# Production: HMR topology + Bussi + MC barostat. 10 ns per segment (nstlim=2500000 * dt=0.004 = 10 ns).
# Optional LJ smoothing: fswitch=8.0 with cut=9.0 is recommended with barostat=2 + PME (not enforced here).
# Checkpoint interval: prefer ntwr=5000 (20 ps checkpoints at dt=0.004), not validated by this gate.
gate_production_protocol() {
  local name="$1"
  local sdir="$2"
  local prod_in="$3"
  local parm_path="$4"

  _require_validate_cntrl

  [[ -f "${prod_in}" ]] || {
    echo "[GATE FAIL] [${name}] Missing production input: ${prod_in}" >&2
    exit 1
  }
  python3 "${VALIDATE_CNTRL}" "${prod_in}" \
    dt=0.004 nstlim=2500000 ntc=2 ntf=2 barostat=2 mcbarint=25 ntt=11 tautp=2.0 || {
    echo "[GATE FAIL] [${name}] production.in must have dt=0.004, nstlim=2500000, ntc=2, ntf=2, barostat=2, mcbarint=25, ntt=11 (Bussi), tautp=2.0." >&2
    exit 1
  }
  _require_nmropt "${prod_in}"

  local hmr_top=""
  if [[ -f "${sdir}/solvated_hmr.prmtop" ]]; then
    hmr_top="${sdir}/solvated_hmr.prmtop"
  elif [[ -f "${sdir}/solvated_hmr.parm7" ]]; then
    hmr_top="${sdir}/solvated_hmr.parm7"
  fi

  if [[ "${SYSTEM_TYPE:-}" == "monomer_heterodimer" ]]; then
    [[ -n "${hmr_top}" ]] || {
      echo "[GATE FAIL] [${name}] monomer_heterodimer production requires solvated_hmr.prmtop or solvated_hmr.parm7 under ${sdir}." >&2
      echo "  Build or run equilibration first so HMR topology exists." >&2
      exit 1
    }
    [[ -n "${parm_path}" && -f "${parm_path}" ]] || {
      echo "[GATE FAIL] [${name}] Invalid production topology (-p): ${parm_path:-missing}" >&2
      exit 1
    }
    case "${parm_path}" in
      *solvated_hmr.prmtop* | *solvated_hmr.parm7*) ;;
      *)
        echo "[GATE FAIL] [${name}] Production -p must be solvated_hmr topology (got ${parm_path})." >&2
        exit 1
        ;;
    esac
    [[ "$(readlink -f "${parm_path}")" == "$(readlink -f "${hmr_top}")" ]] || {
      echo "[GATE FAIL] [${name}] -p topology does not resolve to ${hmr_top}" >&2
      exit 1
    }
  fi
}

# After equil_free: density / RMSD / Ni–Ni convergence (Python + cpptraj).
# Usage: gate_equil_convergence <systems_dir> <inhib_tag> <prmtop> <equil_stage_dir> <system_name>
#   <equil_stage_dir> = .../simulations/<name>/equil_free (contains equil_free.out, equil_free.nc)
# NOTE: EQUIL_OHB defaults to 0 (no OHB in current no-OHB topology). The Python checker
# should skip the OHB analysis when ohb=0.
gate_equil_convergence() {
  local SYSDIR="$1"
  local _INHIB="$2"
  local PRMTOP="$3"
  local EQUIL_DIR="$4"
  local SNAME="$5"
  local CHECK="${SCRIPT_DIR_AMBER_GATE}/check_equil_convergence.py"
  [[ -f "${CHECK}" ]] || {
    echo "[GATE FAIL] Missing ${CHECK}" >&2
    return 1
  }
  chmod +x "${CHECK}" 2>/dev/null || true
  local MDOUT="${EQUIL_DIR}/equil_free.out"
  [[ -f "${MDOUT}" ]] || MDOUT="${EQUIL_DIR}/equil_free.mdout"
  local TRAJ="${EQUIL_DIR}/equil_free.nc"
  if [[ ! -f "${MDOUT}" ]]; then
    echo "[GATE FAIL] equil_free mdout not found (expected ${EQUIL_DIR}/equil_free.out or equil_free.mdout)" >&2
    return 1
  fi
  if [[ ! -f "${TRAJ}" ]]; then
    echo "[GATE FAIL] equil_free.nc not found: ${TRAJ}" >&2
    return 1
  fi
  if [[ ! -f "${PRMTOP}" ]]; then
    echo "[GATE FAIL] prmtop not found: ${PRMTOP}" >&2
    return 1
  fi
  export EQUIL_SYSTEM_NAME="${SNAME}"
  local REPORT="${SYSDIR}/equil_convergence_report.txt"
  local EQIN="${PROJECT_ROOT}/simulations/input_files/equil_free.in"
  local NI1="${EQUIL_NI1:-570}"
  local NI2="${EQUIL_NI2:-571}"
  local OHB="${EQUIL_OHB:-0}"
  if ! python3 "${CHECK}" \
    --mdout "${MDOUT}" \
    --prmtop "${PRMTOP}" \
    --traj "${TRAJ}" \
    --report "${REPORT}" \
    --equil-mdin "${EQIN}" \
    --ni1 "${NI1}" \
    --ni2 "${NI2}" \
    --ohb "${OHB}"; then
    echo "[GATE FAIL] Equilibration convergence check failed." >&2
    echo "  Review: ${REPORT}" >&2
    return 1
  fi
  echo "[GATE PASS] Equilibration convergence confirmed."
  return 0
}
