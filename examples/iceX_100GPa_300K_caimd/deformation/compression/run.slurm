#!/bin/bash
#SBATCH --account=ICT26_CMSP_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --job-name=icex_compression_xx
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

ROOT="/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919"
PW="$ROOT/install-current-gpu/bin/pw.x"
EXPECTED_PW_SHA256="7cacabbef8d32aa83fa561cd124e319ced62714e0a66e2a85f2b1d99def4a6e8"

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

if [[ ! -f qe.in ]]; then
  echo "ERROR: no qe.in in $(pwd)" >&2
  exit 13
fi

mkdir -p tmp evidence
if [[ ! -r ../../pseudo/H.upf && ! -r ../../../pseudo/H.upf ]]; then
  echo "ERROR: expected remote pseudo/H.upf symlink is missing from this case" >&2
  exit 14
fi
if [[ ! -r ../../pseudo/O.upf && ! -r ../../../pseudo/O.upf ]]; then
  echo "ERROR: expected remote pseudo/O.upf symlink is missing from this case" >&2
  exit 15
fi

monitor_log="nvidia-smi-${SLURM_JOB_ID:-manual}.log"
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

echo "JOB_CONTEXT host=$(hostname) pwd=$(pwd) job=${SLURM_JOB_ID:-manual} pw=$PW"
echo "PW_SHA256 $actual_sha"
nvidia-smi || true

srun --mpi=pmix_v3 "$PW" -in qe.in > qe.out 2>&1

if ! grep -q "JOB DONE." qe.out; then
  echo "ERROR: qe.out does not contain JOB DONE." >&2
  exit 20
fi

prefix="$(awk -F"'" '/prefix[[:space:]]*=/{print $2; exit}' qe.in)"
if [[ -n "$prefix" && -f "tmp/${prefix}.caimd" ]]; then
  cp "tmp/${prefix}.caimd" "${prefix}.caimd"
fi

if grep -Eiq "caimd[[:space:]]*=[[:space:]]*\\.true\\." qe.in && grep -Eiq "caimd_mask\\([0-9]+\\).*\\.true\\." qe.in; then
  caimd_file="tmp/${prefix}.caimd"
  if [[ ! -s "$caimd_file" ]]; then
    echo "ERROR: active CAIMD input did not produce $caimd_file" >&2
    exit 21
  fi
  awk 'BEGIN{data=0} /^[[:space:]]*#/ {next} NF>=26 {data++} END{exit(data>0?0:1)}' "$caimd_file" || {
    echo "ERROR: active CAIMD file has no data rows" >&2
    exit 22
  }
fi

echo "JOB DONE: icex_compression_xx completed"
