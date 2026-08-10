#!/usr/bin/env bash
#SBATCH --job-name=d8_vcr74_ibrav
#SBATCH --account=ICT26_CMSP_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --time=02:00:00
#SBATCH --output=qe74/slurm-%j.out
#SBATCH --error=qe74/slurm-%j.err

# JOB74: Diamond8 vc-relax with the official Leonardo QE 7.4.1 module.
# Clean module environment. Independent job (no dependency).
set -uo pipefail

ROOT="/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_vcrelax_qe75_vs_qe74_ibrav_20260728T205518Z"
CASE="qe74"
QE74_MODULE="quantum-espresso/7.4.1--hpcx-mpi--2.19--nvhpc--24.5-openblas"
EXPECT_SHA="beb956cce572605941b791c37bc4dd30ac6cee714e18be7b93fd4ff0e47f4aca"
INPUT="$ROOT/inputs/diamond8_vcrelax_${CASE}.in"

cd "$ROOT"
mkdir -p "$CASE/tmp" "$CASE/analysis" evidence
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

# --- clean module environment: official Leonardo QE 7.4.1 ---
module purge
module load profile/chem-phys
module load "$QE74_MODULE"
PW_X="$(command -v pw.x)"

PY=python3
[ -x "$HOME/venvs/mace_py311/bin/python" ] && PY="$HOME/venvs/mace_py311/bin/python"

if [ -z "$PW_X" ]; then
  echo "BLOCKED_EXECUTABLE_MISMATCH case=$CASE reason=pw.x_not_resolved_after_module_load" \
      > "$ROOT/$CASE/EXECUTABLE_MISMATCH.txt"
  exit 3
fi
ACTUAL_SHA="$(sha256sum "$PW_X" | awk '{print $1}')"
INPUT_SHA="$(sha256sum "$INPUT" | awk '{print $1}')"
{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "hostname=$(hostname)"
  echo "case=$CASE"
  echo "module=$QE74_MODULE"
  echo "executable=$PW_X"
  echo "executable_sha256=$ACTUAL_SHA"
  echo "executable_sha256_expected=$EXPECT_SHA"
  echo "input=$INPUT"
  echo "input_sha256=$INPUT_SHA"
  echo "srun_command=srun --mpi=pmix_v3 -n 1 $PW_X -in $INPUT"
  echo "--- module list ---"
  module list
  echo "--- pw.x --version ---"
  timeout 30 "$PW_X" --version < /dev/null || echo "(version probe produced no output)"
  echo "--- ldd ---"
  ldd "$PW_X" || true
} > "$ROOT/$CASE/RUN_ENVIRONMENT.txt" 2>&1

if [ "$ACTUAL_SHA" != "$EXPECT_SHA" ]; then
  echo "BLOCKED_EXECUTABLE_MISMATCH case=$CASE actual=$ACTUAL_SHA expected=$EXPECT_SHA" \
      > "$ROOT/$CASE/EXECUTABLE_MISMATCH.txt"
  exit 3
fi

# --- run vc-relax ---
cd "$ROOT/$CASE"
srun --mpi=pmix_v3 -n 1 "$PW_X" -in "$INPUT" > qe.out 2> qe.err
RC=$?
echo "$RC" > "$ROOT/$CASE/QE_RETURN_CODE.txt"
cd "$ROOT"

# --- analysis + scientific gate ---
"$PY" scripts/analyze_vcrelax.py \
   --qe-out "$ROOT/$CASE/qe.out" \
   --qe-rc "$RC" \
   --case "$CASE" \
   --version-label "QE 7.4.1 official (Leonardo module $QE74_MODULE)" \
   --forc-conv-thr 1.0e-4 \
   --outdir "$ROOT/$CASE/analysis" \
   --report "$ROOT/$CASE/FINAL_VCRELAX_REPORT.md"
GATE=$?
{
  echo "RC_${CASE}=$RC"
  echo "GATE_${CASE}=$GATE"
} >> "$ROOT/evidence/gate_return_codes.txt"
exit $GATE
