#!/usr/bin/env bash
#SBATCH --job-name=d8_vcr75_ibrav
#SBATCH --account=ICT26_CMSP_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --time=02:00:00
#SBATCH --output=qe75/slurm-%j.out
#SBATCH --error=qe75/slurm-%j.err

# JOB75: Diamond8 vc-relax with the custom CAIMD-capable QE 7.5 build af7ebfc.
# CAIMD is absent from the input. Independent job (no dependency).
set -uo pipefail

ROOT="/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_vcrelax_qe75_vs_qe74_ibrav_20260728T205518Z"
CASE="qe75"
QE75_ROOT="/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc"
PW_X="$QE75_ROOT/install-current-gpu/bin/pw.x"
EXPECT_SHA="89d4fda827d2cc22a3b4c1670e90dc40ce629af42096c785b55c70c355a10e7d"
INPUT="$ROOT/inputs/diamond8_vcrelax_${CASE}.in"

cd "$ROOT"
mkdir -p "$CASE/tmp" "$CASE/analysis" evidence
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

# --- environment: reuse the proven custom QE 7.5 build environment ---
module purge
if [ -f "$QE75_ROOT/scripts/load_build_environment.sh" ]; then
  # shellcheck source=/dev/null
  source "$QE75_ROOT/scripts/load_build_environment.sh"
else
  module load profile/chem-phys
  module load nvhpc/24.5 hpcx-mpi/2.19
  module load fftw/3.3.10--hpcx-mpi--2.19--nvhpc--24.5 libxc/7.0.0--nvhpc--24.5
  module load openblas/0.3.26--nvhpc--24.5 netlib-scalapack/2.2.0--hpcx-mpi--2.19--nvhpc--24.5
fi
export LD_LIBRARY_PATH="$QE75_ROOT/local-libs/lib:${LD_LIBRARY_PATH:-}"

PY=python3
[ -x "$HOME/venvs/mace_py311/bin/python" ] && PY="$HOME/venvs/mace_py311/bin/python"

ACTUAL_SHA="$(sha256sum "$PW_X" | awk '{print $1}')"
INPUT_SHA="$(sha256sum "$INPUT" | awk '{print $1}')"
{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "hostname=$(hostname)"
  echo "case=$CASE"
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
   --version-label "QE 7.5 custom af7ebfc (CAIMD disabled/absent)" \
   --forc-conv-thr 1.0e-4 \
   --outdir "$ROOT/$CASE/analysis" \
   --report "$ROOT/$CASE/FINAL_VCRELAX_REPORT.md"
GATE=$?
{
  echo "RC_${CASE}=$RC"
  echo "GATE_${CASE}=$GATE"
} >> "$ROOT/evidence/gate_return_codes.txt"
exit $GATE
