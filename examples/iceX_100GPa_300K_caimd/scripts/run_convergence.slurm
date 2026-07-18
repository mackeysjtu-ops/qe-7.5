#!/bin/bash
#SBATCH --account=ICT26_CMSP_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --job-name=icex_matrix
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

ROOT="/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919"
EXAMPLE_ROOT="$ROOT/example_iceX_100GPa_300K_caimd"
PW="$ROOT/install-current-gpu/bin/pw.x"
EXPECTED_PW_SHA256="c9588df640013f098a33e33e1eec6e1f0820ed53f7488efffb0a83ec7ed25b48"
CASES_FILE="${CASES_FILE:-convergence/cases.txt}"

cd "$EXAMPLE_ROOT"

if [[ ! -x "$PW" ]]; then
  echo "ERROR: pw.x not executable at $PW" >&2
  exit 10
fi
actual_sha="$(sha256sum "$PW" | awk '{print $1}')"
if [[ "$actual_sha" != "$EXPECTED_PW_SHA256" ]]; then
  echo "ERROR: pw.x checksum mismatch: got $actual_sha expected $EXPECTED_PW_SHA256" >&2
  exit 11
fi
source "$ROOT/scripts/load_build_environment.sh"
if ldd "$PW" | grep -qi "not found"; then
  echo "ERROR: unresolved pw.x shared-library dependency" >&2
  ldd "$PW" >&2
  exit 12
fi
if [[ ! -f "$CASES_FILE" ]]; then
  echo "ERROR: missing cases file $CASES_FILE" >&2
  exit 13
fi
if [[ ! -r pseudo/H.upf || ! -r pseudo/O.upf ]]; then
  echo "ERROR: remote pseudo symlinks pseudo/H.upf and pseudo/O.upf are required" >&2
  exit 14
fi

mkdir -p evidence
summary="evidence/matrix-${SLURM_JOB_ID:-manual}.tsv"
printf "case\tstatus\texit_code\tjob_done\tcaimd_rows\tstart\tend\n" > "$summary"

monitor_log="evidence/nvidia-smi-matrix-${SLURM_JOB_ID:-manual}.log"
(
  while true; do
    date -Is
    nvidia-smi --query-gpu=timestamp,index,name,utilization.gpu,utilization.memory,memory.used,power.draw --format=csv,noheader,nounits || true
    sleep 30
  done
) > "$monitor_log" 2>&1 &
monitor_pid=$!
cleanup() {
  rc=$?
  kill "$monitor_pid" 2>/dev/null || true
  wait "$monitor_pid" 2>/dev/null || true
  exit "$rc"
}
trap cleanup EXIT

failures=0
echo "JOB_CONTEXT host=$(hostname) pwd=$(pwd) job=${SLURM_JOB_ID:-manual} cases=$CASES_FILE pw=$PW"
echo "PW_SHA256 $actual_sha"
nvidia-smi || true

while IFS= read -r case_dir; do
  [[ -z "$case_dir" || "$case_dir" =~ ^# ]] && continue
  if [[ ! -f "$case_dir/qe.in" ]]; then
    printf "%s\tmissing_qe_input\t99\t0\t0\t%s\t%s\n" "$case_dir" "$(date -Is)" "$(date -Is)" >> "$summary"
    failures=$((failures + 1))
    continue
  fi
  start="$(date -Is)"
  echo "=== CASE START $case_dir $start ==="
  mkdir -p "$case_dir/tmp" "$case_dir/evidence"
  set +e
  (cd "$case_dir" && srun --mpi=pmix_v3 "$PW" -in qe.in < /dev/null > qe.out 2>&1)
  rc=$?
  set -e
  end="$(date -Is)"
  job_done=0
  caimd_rows=0
  if grep -q "JOB DONE." "$case_dir/qe.out" 2>/dev/null; then
    job_done=1
  fi
  prefix="$(awk -F"'" '/prefix[[:space:]]*=/{print $2; exit}' "$case_dir/qe.in")"
  if [[ -n "$prefix" && -f "$case_dir/tmp/${prefix}.caimd" ]]; then
    cp "$case_dir/tmp/${prefix}.caimd" "$case_dir/${prefix}.caimd"
    caimd_rows="$(awk 'BEGIN{data=0} /^[[:space:]]*#/ {next} NF>=26 {data++} END{print data+0}' "$case_dir/tmp/${prefix}.caimd")"
  fi
  status="ok"
  if [[ "$rc" -ne 0 || "$job_done" -ne 1 ]]; then
    status="failed"
    failures=$((failures + 1))
  fi
  if grep -Eiq "caimd[[:space:]]*=[[:space:]]*\\.true\\." "$case_dir/qe.in" && grep -Eiq "caimd_mask\\([0-9]+\\).*\\.true\\." "$case_dir/qe.in"; then
    if [[ "$caimd_rows" -lt 1 ]]; then
      status="failed_caimd_missing"
      failures=$((failures + 1))
    fi
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$case_dir" "$status" "$rc" "$job_done" "$caimd_rows" "$start" "$end" >> "$summary"
  echo "=== CASE END $case_dir $status rc=$rc job_done=$job_done caimd_rows=$caimd_rows $end ==="
done < "$CASES_FILE"

if [[ "$failures" -ne 0 ]]; then
  echo "ERROR: $failures case(s) failed; see $summary" >&2
  exit 30
fi

echo "JOB DONE: Ice X bounded matrix completed"
