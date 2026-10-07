#!/usr/bin/env bash
# scripts/post_equil_lib.sh — shared helpers for post-equil3 / long-NPT pipeline
# shellcheck disable=SC2034  # sourced; not all vars used in every caller
# Source from project scripts after config.env (PROJECT_ROOT set).
#
# Fortran/pmemd truncates path buffers (~80 chars). Short absolute symlinks live
# under ${PROJECT_ROOT}/{m,p,r}/ (mdin, prmtop, input coordinates).

post_equil__path_len() {
  echo -n "$1" | wc -c
}

post_equil__assert_path_short() {
  local label=$1 path=$2
  local n
  n=$(post_equil__path_len "${path}")
  if [[ "${n}" -ge 80 ]]; then
    echo "[FAIL] ${label} path length ${n} (>=80): ${path}" >&2
    return 1
  fi
  return 0
}

post_equil__cpptraj_bin() {
  command -v cpptraj 2>/dev/null || echo "/usr/local/bin/cpptraj"
}

# Extract first two amber atom numbers for atoms named NI (@NI) using one trajectory frame.
# Args: prmtop rst7 mask_out_file
# Sets globals: NI1_ATOM NI2_ATOM
post_equil_extract_ni_atoms() {
  local prmtop=$1 rst7=$2 mask_out=$3
  local cpptraj
  cpptraj="$(post_equil__cpptraj_bin)"
  [[ -x "${cpptraj}" || -f "${cpptraj}" ]] || {
    echo "[FAIL] cpptraj not found (tried PATH and /usr/local/bin/cpptraj)" >&2
    return 1
  }
  [[ -f "${prmtop}" && -f "${rst7}" ]] || {
    echo "[FAIL] Missing prmtop or rst7: ${prmtop} / ${rst7}" >&2
    return 1
  }
  mkdir -p "$(dirname "${mask_out}")"
  "${cpptraj}" -p "${prmtop}" <<EOF >/dev/null 2>&1
trajin ${rst7}
mask @NI out ${mask_out}
run
quit
EOF
  mapfile -t _ni_nums < <(awk '$3=="NI" {print $2}' "${mask_out}" | head -2)
  if [[ ${#_ni_nums[@]} -ne 2 ]]; then
    echo "[FAIL] Expected 2 @NI atoms, got ${#_ni_nums[@]} (see ${mask_out})" >&2
    return 1
  fi
  NI1_ATOM="${_ni_nums[0]}"
  NI2_ATOM="${_ni_nums[1]}"
  return 0
}

# Short absolute DISANG path for pmemd Fortran line-length limits.
post_equil_disang_short_abs() {
  local sys=$1 real_rst=$2
  local d="${PROJECT_ROOT}/nmr_disang"
  mkdir -p "${d}"
  local link="${d}/${sys}.RST"
  ln -sf "$(readlink -f "${real_rst}" 2>/dev/null || echo "${real_rst}")" "${link}"
  printf '%s' "${link}"
}

# Write flat-bottom Ni–Ni restraint (rk2=rk3=50 kcal/mol/Å²).
post_equil_write_ni_ni_rst() {
  local path=$1 ni1=$2 ni2=$3
  mkdir -p "$(dirname "${path}")"
  cat > "${path}" <<EOF
# Ni-Ni flat-bottom harmonic distance restraint (post_equil pipeline)
# iat = amber atom indices for this topology
&rst
  iat=${ni1},${ni2},
  r1=2.5, r2=3.3, r3=3.7, r4=4.5,
  rk2=50.0, rk3=50.0,
/
EOF
}

post_equil__write_one_mdin() {
  local out=$1 nstlim=$2 disang_abs=$3
  mkdir -p "$(dirname "${out}")"
  local mdin
  mdin=$'Long NPT segment / production — post_equil protocol\n'
  mdin+=$' &cntrl\n'
  mdin+="  imin=0, nstlim=${nstlim}, dt=0.004,"$'\n'
  mdin+=$'  ntx=5, irest=1,\n'
  mdin+=$'  ntb=2, ntp=1,\n'
  mdin+=$'  ntt=11, tautp=2.0, temp0=300.0,\n'
  mdin+=$'  ntc=2, ntf=2,\n'
  mdin+=$'  cut=9.0,\n'
  mdin+=$'  ntpr=5000, ntwx=5000, ntwr=25000,\n'
  mdin+=$'  ioutfm=1, ntxo=2,\n'
  mdin+=$'  barostat=2, mcbarint=100, pres0=1.0,\n'
  mdin+=$'  nmropt=1,\n'
  mdin+=$'  ntr=0,\n'
  mdin+=$'  ig=-1,\n'
  mdin+=$' /\n'
  mdin+=$' &wt type=\'END\' /\n'
  mdin+="DISANG=${disang_abs}"$'\n'
  printf '%s' "${mdin}" > "${out}"
}

# equil_seg.in: nstlim=2500000; production.in: nstlim=5000000 (same MC / thermostat / DISANG).
post_equil_write_mdin_pair() {
  local disang_abs=$1 equil_path=$2 prod_path=$3
  post_equil__write_one_mdin "${equil_path}" 2500000 "${disang_abs}"
  post_equil__write_one_mdin "${prod_path}" 5000000 "${disang_abs}"
}

# Short symlinks for mdin + prmtop (<80 char absolute paths). Call after per-system mdin files exist.
# Optional second arg: real coordinate file (rst7) — updates r/${sys}_in.rst7 before each pmemd -c.
post_equil_ensure_short_links() {
  local sys=$1
  local real_crd=${2:-}
  local equil_real="${PROJECT_ROOT}/simulations/${sys}/input_files/equil_seg.in"
  local prod_real="${PROJECT_ROOT}/simulations/${sys}/input_files/production.in"
  local prm_real="${PROJECT_ROOT}/systems/${sys}/solvated_hmr.prmtop"
  local md_dir="${PROJECT_ROOT}/m"
  local parm_dir="${PROJECT_ROOT}/p"
  mkdir -p "${md_dir}" "${parm_dir}"

  [[ -f "${equil_real}" ]] || {
    echo "[FAIL] Missing equil mdin: ${equil_real}" >&2
    return 1
  }
  [[ -f "${prod_real}" ]] || {
    echo "[FAIL] Missing production mdin: ${prod_real}" >&2
    return 1
  }
  [[ -f "${prm_real}" ]] || {
    echo "[FAIL] Missing topology: ${prm_real}" >&2
    return 1
  }

  local e_short="${md_dir}/${sys}_eq.in"
  local p_short="${md_dir}/${sys}_pr.in"
  local pr_short="${parm_dir}/${sys}.prmtop"

  ln -sf "${equil_real}" "${e_short}"
  ln -sf "${prod_real}" "${p_short}"
  ln -sf "${prm_real}" "${pr_short}"

  local el pl prl
  el=$(post_equil__path_len "${e_short}")
  pl=$(post_equil__path_len "${p_short}")
  prl=$(post_equil__path_len "${pr_short}")
  if [[ "${el}" -ge 80 ]] || [[ "${pl}" -ge 80 ]] || [[ "${prl}" -ge 80 ]]; then
    echo "[FAIL] Short path >=80 chars: equil_mdin=${el} prod_mdin=${pl} prmtop=${prl}" >&2
    echo "       Shorten PROJECT_ROOT or use a shorter symlink directory name." >&2
    return 1
  fi
  echo "[OK] Short paths (${el}/${pl}/${prl} chars):"
  echo "     -i ${e_short}"
  echo "     -p ${pr_short}"
  if [[ -n "${real_crd}" ]]; then
    post_equil_short_crd_path "${sys}" "${real_crd}" >/dev/null || return 1
    echo "     -c ${PROJECT_ROOT}/r/${sys}_in.rst7 -> $(readlink -f "${real_crd}" 2>/dev/null || echo "${real_crd}")"
  fi
  return 0
}

# Per-run coordinate symlink (<80 chars): overwrites r/${sys}_in.rst7
post_equil_short_crd_path() {
  local sys=$1 real_crd=$2
  local r_dir="${PROJECT_ROOT}/r"
  mkdir -p "${r_dir}"
  local tgt
  tgt=$(readlink -f "${real_crd}" 2>/dev/null || echo "${real_crd}")
  local l="${r_dir}/${sys}_in.rst7"
  ln -sf "${tgt}" "${l}"
  post_equil__assert_path_short "crd -c" "${l}" || return 1
  printf '%s' "${l}"
}

# Args: run_dir SYS equil|prod real_crd prefix pmemd_bin
post_equil_run_pmemd_segment() {
  local run_dir=$1 sys=$2 mkind=$3 real_crd=$4 prefix=$5 pmemd=$6
  local mdin prmtop crd

  if ! post_equil_ensure_short_links "${sys}" "${real_crd}"; then
    return 1
  fi

  prmtop="${PROJECT_ROOT}/p/${sys}.prmtop"
  crd="${PROJECT_ROOT}/r/${sys}_in.rst7"
  [[ -e "${crd}" ]] || {
    echo "[FAIL] Short crd symlink missing: ${crd}" >&2
    return 1
  }
  post_equil__assert_path_short "crd -c" "${crd}" || return 1

  case "${mkind}" in
    equil) mdin="${PROJECT_ROOT}/m/${sys}_eq.in" ;;
    prod) mdin="${PROJECT_ROOT}/m/${sys}_pr.in" ;;
    *)
      echo "[FAIL] mkind must be equil or prod (got ${mkind})" >&2
      return 1
      ;;
  esac

  if [[ ! -e "${mdin}" ]]; then
    echo "[FAIL] Short mdin missing: ${mdin}" >&2
    echo "       Run: bash scripts/build_inputs_for_system.sh ${sys}" >&2
    return 1
  fi

  post_equil__assert_path_short "mdin -i" "${mdin}" || return 1
  post_equil__assert_path_short "prmtop -p" "${prmtop}" || return 1

  mkdir -p "${run_dir}"
  (
    cd "${run_dir}" || exit 1
    "${pmemd}" -O \
      -i "${mdin}" \
      -p "${prmtop}" \
      -c "${crd}" \
      -o "${prefix}.out" \
      -r "${prefix}.rst7" \
      -x "${prefix}.nc"
    rc=$?
    if [[ -f mdinfo ]]; then
      mv -f mdinfo "${prefix}.mdinfo"
    fi
    exit "${rc}"
  )
}

# Args: mdout expected_final_nstep (e.g. 2500000 or 5000000)
post_equil_verify_mdout() {
  local out=$1 expected=${2:-2500000}
  [[ -f "${out}" ]] || {
    echo "[FAIL] Missing mdout: ${out}" >&2
    return 1
  }
  local nan_c
  nan_c=$(grep -c 'NaN' "${out}" 2>/dev/null || true)
  if [[ "${nan_c}" != "0" ]]; then
    echo "[FAIL] NaN count=${nan_c} in ${out}" >&2
    return 1
  fi
  if ! grep -qE "NSTEP[[:space:]]*=[[:space:]]*${expected}" "${out}"; then
    echo "[FAIL] No final NSTEP=${expected} line in ${out}" >&2
    return 1
  fi
  return 0
}
