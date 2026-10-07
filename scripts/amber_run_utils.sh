#!/usr/bin/env bash
set -euo pipefail

# Utilities for robust Amber pmemd.cuda execution.
#
# Key feature: automatic restart/continue loop when pmemd.cuda halts with:
#   "Periodic box dimensions have changed too much from their initial values."
#
# This is a known early-density-equilibration issue on GPU when grid cell
# decomposition needs rebuilding; restarting from the last restart regenerates
# grid cells and allows the simulation to continue.

pmemd_boxchange_restart_loop() {
  # Arguments:
  #   1: stage name (string)
  #   2: pmemd executable (path)
  #   3: mdin file (path)
  #   4: prmtop (path)
  #   5: initial coord/restart input (path)  (-c)
  #   6: restart output base (path)          (-r)
  #   7: mdout base (path)                   (-o)
  #   8: netcdf traj base (path or empty)    (-x) optional
  #   9: reference coord (path or empty)     (-ref) optional
  #
  # Behavior:
  # - Runs pmemd.cuda normally.
  # - If it exits non-zero AND mdout contains the box-change error, it will:
  #     - save the last restart as <rst>.boxretryN
  #     - rerun using that restart as the new -c
  #   repeated until success or MAX_BOX_RESTARTS is reached.

  local stage="$1"; shift
  local PMEMD="$1"; shift
  local mdin="$1"; shift
  local prmtop="$1"; shift
  local c_in="$1"; shift
  local rst_out="$1"; shift
  local out_base="$1"; shift
  local nc_out="${1:-}"; shift || true
  local ref_in="${1:-}"; shift || true

  local MAX_BOX_RESTARTS="${MAX_BOX_RESTARTS:-20}"
  # Match a stable substring; exact wording can wrap across lines.
  local BOX_ERR="Periodic box dimensions have changed too much"
  local GRID_HINT="generate new grid cells"
  local DENS_HINT="density has likely changed by a large amount"

  local attempt=0
  local cur_c="${c_in}"
  local saw_boxchange=0

  while true; do
    attempt=$((attempt+1))
    local out="${out_base}"
    if [[ "${attempt}" -gt 1 ]]; then
      out="${out_base%.out}.boxretry${attempt}.out"
    fi

    echo "[$(date '+%F %T')] START ${stage} (attempt ${attempt}, -c ${cur_c})"

    # Build pmemd command.
    local cmd=("${PMEMD}" -O -i "${mdin}" -p "${prmtop}" -c "${cur_c}" -r "${rst_out}" -o "${out}")
    if [[ -n "${nc_out}" ]]; then
      local nc="${nc_out}"
      if [[ "${attempt}" -gt 1 ]]; then
        nc="${nc_out%.nc}.boxretry${attempt}.nc"
      fi
      cmd+=(-x "${nc}")
    fi
    if [[ -n "${ref_in}" ]]; then
      cmd+=(-ref "${ref_in}")
    fi

    # pmemd.cuda prints some fatal GPU messages to stderr (not always captured in -o mdout).
    # Append stdout/stderr to the mdout file so we can reliably detect the box-change condition.
    set +e
    "${cmd[@]}" >>"${out}" 2>&1
    local rc=$?
    set -e

    if [[ "${rc}" -eq 0 ]]; then
      echo "[$(date '+%F %T')] DONE  ${stage} (attempt ${attempt})"
      # For convenience, keep a canonical out file if we had retries.
      if [[ "${attempt}" -gt 1 ]]; then
        cp -f "${out}" "${out_base}"
      fi
      return 0
    fi

    # Auto-retry the specific box-change condition (sometimes only partially written).
    local is_boxchange=0
    if [[ -f "${out}" ]] && ( grep -Fq "${BOX_ERR}" "${out}" || grep -Fq "${GRID_HINT}" "${out}" || grep -Fq "${DENS_HINT}" "${out}" ); then
      is_boxchange=1
    elif [[ "${saw_boxchange}" -eq 1 && "${rc}" -eq 255 ]]; then
      # After an initial confirmed box-change halt, subsequent retries can still exit 255
      # even if the message is not fully flushed into mdout before termination.
      is_boxchange=1
    elif [[ "${rc}" -eq 255 && -s "${rst_out}" ]]; then
      # Last-resort heuristic: CUDA NPT box-change halts often exit 255. If we at least have
      # a usable restart and there is no obvious different CUDA error, try a restart.
      if [[ -f "${out}" ]] && ! grep -Fqi "illegal memory access" "${out}" && ! grep -Fqi "Bad topology file" "${out}"; then
        is_boxchange=1
      fi
    fi

    if [[ "${is_boxchange}" -eq 1 ]]; then
      saw_boxchange=1
      if [[ "${attempt}" -ge "${MAX_BOX_RESTARTS}" ]]; then
        echo "ERROR: ${stage} hit box-change restart limit (${MAX_BOX_RESTARTS})." >&2
        return 10
      fi
      if [[ ! -s "${rst_out}" ]]; then
        echo "ERROR: ${stage} box-change detected but restart file missing/empty: ${rst_out}" >&2
        return 11
      fi
      # Save the last restart and continue from it.
      cp -f "${rst_out}" "${rst_out}.boxretry${attempt}"
      cur_c="${rst_out}"
      echo "[$(date '+%F %T')] RETRY ${stage}: box changed too much; continuing from ${cur_c}"
      continue
    fi

    echo "ERROR: ${stage} failed (rc=${rc}); not a recognized box-change case. See ${out}" >&2
    return "${rc}"
  done
}

