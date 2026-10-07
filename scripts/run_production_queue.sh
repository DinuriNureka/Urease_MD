#!/usr/bin/env bash
set -euo pipefail

# Queue launcher for production MD runs across available GPUs with headroom rules.
#
# Policy:
# - Use at most 80% of CPU cores (leave >=20% free).
# - Only start a new GPU job if a GPU has enough free VRAM.
# - Auto-discover jobs by finding executable run_production.sh scripts under simulations/*/
# - Start next job when a running job finishes.
#
# Usage:
#   bash scripts/run_production_queue.sh
#
# Optional env vars:
#   GPU_ALLOWLIST="1"            # e.g. "1" or "0,1"
#   MIN_FREE_VRAM_MB=6000        # require at least this much free VRAM to launch a job on a GPU
#   MAX_CPU_FRACTION=0.80        # default 0.80
#   POLL_SECONDS=20              # queue polling interval
#   MAX_BOX_RESTARTS=20          # max automatic box-change restarts per stage

ROOT=/path/to/urease_md
cd "${ROOT}"

# Make CUDA device numbering match nvidia-smi (stable mapping by PCI bus ID).
export CUDA_DEVICE_ORDER=PCI_BUS_ID

GPU_ALLOWLIST=${GPU_ALLOWLIST:-""}
MIN_FREE_VRAM_MB=${MIN_FREE_VRAM_MB:-6000}
MAX_CPU_FRACTION=${MAX_CPU_FRACTION:-0.80}
POLL_SECONDS=${POLL_SECONDS:-20}
RE_DISCOVER_SECONDS=${RE_DISCOVER_SECONDS:-1800}

SINGLE_GPU_MODE=false
if [[ -n "${GPU_ALLOWLIST}" ]]; then
  IFS=',' read -r -a _allow <<<"${GPU_ALLOWLIST//[[:space:]]/}"
  if [[ ${#_allow[@]} -eq 1 ]]; then
    SINGLE_GPU_MODE=true
  fi
fi

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing required command: $1" >&2; exit 1; }
}
need_cmd nvidia-smi
need_cmd awk
need_cmd sort
need_cmd sed
need_cmd bash
need_cmd python3

n_cores=$(nproc)
max_cores=$(python3 - <<PY
import math
n=${n_cores}
f=float("${MAX_CPU_FRACTION}")
print(max(1, int(math.floor(n*f))))
PY
)
leave_free=$(( n_cores - max_cores ))

echo "[queue] cores: ${n_cores}; using <= ${max_cores} (leave >= ${leave_free} free)"
echo "[queue] MIN_FREE_VRAM_MB=${MIN_FREE_VRAM_MB}; POLL_SECONDS=${POLL_SECONDS}"

is_50ns_compound() {
  # Return 0 (true) if compound should run 50 ns production total.
  # All other compounds run 10 ns total.
  # 50 ns list: mol1, mol2, mol3, mol4, mol5, mol6, mol7, mol8, mol15, mol16, mol17, mol18
  local mol="$1"
  case "$mol" in
    mol1|mol2|mol3|mol4|mol5|mol6|mol7|mol8|mol15|mol16|mol17|mol18) return 0 ;;
    *) return 1 ;;
  esac
}

compound_from_job() {
  # Extract "molN" from the job directory name, if present.
  local job="$1"
  local d
  d=$(basename "$(dirname "$job")")
  if [[ "$d" =~ (mol[0-9]+) ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo ""
  fi
}

ensure_prod_extend_10ns_in() {
  local job="$1"
  local d
  d=$(dirname "$job")
  local f="${d}/prod_extend_10ns.in"
  if [[ -f "$f" ]]; then
    return 0
  fi
  cat >"$f" <<'EOF'
10 ns extension from 10 ns restart (HMR, 4 fs)
&cntrl
 imin=0,irest=1,ntx=5,
 nstlim=2500000,dt=0.004,
 ntc=1,ntf=2,
 ntb=2,ntp=1,pres0=1.0,taup=1.0,
 barostat=2,
 ntt=3,gamma_ln=2.0,temp0=300.0,
 ntr=0,
 cut=9.0,
 ntpr=5000,ntwx=5000,ntwr=50000,
 ioutfm=1,ntxo=2,
 ig=-1,
/
EOF
}

discover_jobs() {
  # Print absolute paths to run_production.sh scripts (sorted) that are executable.
  # We avoid 'find' to keep it simple/portable here by using bash globbing.
  shopt -s nullglob
  local jobs=()
  for f in simulations/*/run_production.sh; do
    if [[ -x "$f" ]]; then
      jobs+=("$ROOT/$f")
    fi
  done
  printf "%s\n" "${jobs[@]}" | sort
}

now_epoch() { date +%s; }

LAST_DISCOVER=0
declare -a JOBS=()

rediscover_jobs() {
  local discovered=()
  mapfile -t discovered < <(discover_jobs)
  if [[ ${#discovered[@]} -eq 0 ]]; then
    echo "[queue] No jobs found. Expected executable simulations/*/run_production.sh" >&2
    return 1
  fi
  JOBS=("${discovered[@]}")
  echo "[queue] discovered ${#JOBS[@]} jobs"
  return 0
}

rediscover_jobs || exit 1

gpu_list() {
  # Print GPU indices, one per line.
  nvidia-smi --query-gpu=index --format=csv,noheader,nounits | sed 's/[[:space:]]//g'
}

gpu_name() {
  local idx="$1"
  nvidia-smi --query-gpu=name --format=csv,noheader,nounits -i "$idx" | head -n 1
}

gpu_free_mb() {
  local idx="$1"
  # total-used
  local line
  line=$(nvidia-smi --query-gpu=memory.total,memory.used --format=csv,noheader,nounits -i "$idx" | head -n 1)
  local total used
  total=$(awk -F',' '{gsub(/ /,"",$1); print $1}' <<<"$line")
  used=$(awk -F',' '{gsub(/ /,"",$2); print $2}' <<<"$line")
  echo $(( total - used ))
}

allow_gpu() {
  local idx="$1"
  if [[ -z "${GPU_ALLOWLIST}" ]]; then
    return 0
  fi
  # GPU_ALLOWLIST supports "0,1" or "1"
  IFS=',' read -r -a allow <<<"${GPU_ALLOWLIST}"
  for a in "${allow[@]}"; do
    a="${a//[[:space:]]/}"
    if [[ "$a" == "$idx" ]]; then
      return 0
    fi
  done
  return 1
}

job_is_running() {
  local job="$1"
  # If any process cmdline contains the job script path, treat as running.
  # (This prevents duplicate launches when restarting the queue mid-run.)
  pgrep -af "${job}" >/dev/null 2>&1
}

any_pmemd_cuda_running() {
  pgrep -af "pmemd\\.cuda" >/dev/null 2>&1
}

declare -A RUNNING_PID_BY_GPU=()
declare -A RUNNING_JOB_BY_GPU=()

active_job_count() {
  local c=0
  for g in "${!RUNNING_PID_BY_GPU[@]}"; do
    if kill -0 "${RUNNING_PID_BY_GPU[$g]}" 2>/dev/null; then
      c=$((c+1))
    fi
  done
  echo "$c"
}

cleanup_finished() {
  for g in "${!RUNNING_PID_BY_GPU[@]}"; do
    local pid="${RUNNING_PID_BY_GPU[$g]}"
    if ! kill -0 "$pid" 2>/dev/null; then
      echo "[queue] GPU${g} finished job: ${RUNNING_JOB_BY_GPU[$g]} (pid ${pid})"
      unset "RUNNING_PID_BY_GPU[$g]"
      unset "RUNNING_JOB_BY_GPU[$g]"
    fi
  done
}

launch_job_on_gpu() {
  local job="$1"
  local gpu="$2"
  local mol
  mol=$(compound_from_job "$job")
  local target_ns=10
  if [[ -n "$mol" ]] && is_50ns_compound "$mol"; then
    target_ns=50
  else
    # Keep extension script available for later manual extension if desired.
    ensure_prod_extend_10ns_in "$job"
  fi

  local running
  running=$(active_job_count)
  local new_running=$(( running + 1 ))
  local omp=$(( max_cores / new_running ))
  if [[ "$omp" -lt 1 ]]; then omp=1; fi

  local log_dir
  log_dir=$(dirname "$job")
  local ts
  ts=$(date '+%Y%m%d_%H%M%S')
  local log="${log_dir}/queue_${ts}_gpu${gpu}.log"

  echo "[queue] launch on GPU${gpu} (free=$(gpu_free_mb "$gpu")MB, omp=${omp}): ${job}"
  (
    export CUDA_VISIBLE_DEVICES="${gpu}"
    export OMP_NUM_THREADS="${omp}"
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export MAX_BOX_RESTARTS="${MAX_BOX_RESTARTS:-20}"
    export PROD_TARGET_NS="${target_ns}"
    cd "${ROOT}"
    bash "${job}"
  ) >"${log}" 2>&1 &

  local pid=$!
  RUNNING_PID_BY_GPU["$gpu"]="$pid"
  RUNNING_JOB_BY_GPU["$gpu"]="$job"
  echo "[queue] started pid=${pid}; log=${log}"
}

job_done_marker() {
  local job="$1"
  # Treat a job as "done" if the appropriate final restart exists in its folder.
  local d
  d=$(dirname "$job")
  local mol
  mol=$(compound_from_job "$job")
  if [[ -n "$mol" ]] && is_50ns_compound "$mol"; then
    [[ -f "${d}/prod_50ns.rst7" ]]
    return
  fi
  [[ -f "${d}/prod_10ns.rst7" ]]
}

JOB_IDX=0

while true; do
  cleanup_finished

  # Periodically re-discover jobs so newly built systems get picked up automatically.
  # This is especially useful since systems can be generated while the queue is running.
  now=$(now_epoch)
  if (( now - LAST_DISCOVER >= RE_DISCOVER_SECONDS )); then
    rediscover_jobs || true
    LAST_DISCOVER="$now"
  fi

  # If all jobs already done, finish.
  local_all_done=true
  for j in "${JOBS[@]}"; do
    if ! job_done_marker "$j"; then
      local_all_done=false
      break
    fi
  done
  if [[ "${local_all_done}" == "true" ]]; then
    echo "[queue] all jobs complete (found prod_50ns.rst7 markers)"
    exit 0
  fi

  # Try to launch jobs on any free GPU.
  # In single-GPU mode, be conservative: if pmemd.cuda is already running, do not start a second job.
  # This avoids OOM crashes caused by concurrent dodecamer-scale topologies.
  if [[ "${SINGLE_GPU_MODE}" == "true" ]]; then
    if any_pmemd_cuda_running; then
      sleep "${POLL_SECONDS}"
      continue
    fi
  fi
  for g in $(gpu_list); do
    if ! allow_gpu "$g"; then
      continue
    fi
    # skip if GPU already has a running job
    if [[ -n "${RUNNING_PID_BY_GPU[$g]:-}" ]] && kill -0 "${RUNNING_PID_BY_GPU[$g]}" 2>/dev/null; then
      continue
    fi
    free=$(gpu_free_mb "$g")
    if [[ "$free" -lt "$MIN_FREE_VRAM_MB" ]]; then
      continue
    fi

    # Find next not-yet-done job.
    launched=false
    for j in "${JOBS[@]}"; do
      if job_done_marker "$j"; then
        continue
      fi
      # Avoid launching same job twice
      already=false
      for gg in "${!RUNNING_JOB_BY_GPU[@]}"; do
        if [[ "${RUNNING_JOB_BY_GPU[$gg]}" == "$j" ]]; then
          already=true
          break
        fi
      done
      if [[ "${already}" == "true" ]]; then
        continue
      fi

      # Skip jobs already running elsewhere (e.g., when restarting the queue).
      if job_is_running "$j"; then
        continue
      fi

      launch_job_on_gpu "$j" "$g"
      launched=true
      break
    done
  done

  sleep "${POLL_SECONDS}"
done

