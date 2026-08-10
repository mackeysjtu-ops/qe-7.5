#!/usr/bin/env bash
#SBATCH --job-name=dia8_init_vs_nose
#SBATCH --account=ICT26_CMSP_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# Single-factor thermostat isolation, one executable (QE 7.5 af7ebfc).
#   CONTROL : ion_temperature='initial'  (native 300 K init, no continuing thermostat)
#   TEST    : ion_temperature='nose'     (fnosep=50, nhpcl=4)
# Both: no explicit velocities, caimd=.false., identical structure/electronics/PR.
#
# NOTE: errexit is deliberately NOT enabled, so the Nose run always executes even
# if the control returns non-zero, and vice versa.

set -uo pipefail

WORKDIR="/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_af7ebfc_no_velocity_initial_vs_nose_300K_20260724T120726Z"
QE_ROOT="/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc"
PW_X="$QE_ROOT/install-current-gpu/bin/pw.x"

cd "$WORKDIR"
mkdir -p evidence analysis control_initial/tmp test_nose/tmp

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

{
  echo "job_id=${SLURM_JOB_ID:-unknown}"
  echo "job_name=${SLURM_JOB_NAME:-unknown}"
  echo "nodelist=${SLURM_JOB_NODELIST:-unknown}"
  echo "submit_dir=${SLURM_SUBMIT_DIR:-unknown}"
  echo "workdir=$WORKDIR"
  echo "omp_num_threads=$OMP_NUM_THREADS"
  hostname
  date -u
  env | grep -E '^(SLURM|CUDA|NVHPC|OMPI|OMP)_' | sort || true
} > evidence/slurm_job_environment.txt 2>&1

sha256sum inputs/control_initial_no_velocity.in \
          inputs/test_nose_no_velocity.in \
          pseudo/C.UPF > evidence/input_checksums.txt 2>&1

# ------------------------------------------------------------ environment ---
module purge
if [ -f "$QE_ROOT/scripts/load_build_environment.sh" ]; then
  # shellcheck source=/dev/null
  source "$QE_ROOT/scripts/load_build_environment.sh"
else
  module load profile/chem-phys
  module load nvhpc/24.5
  module load hpcx-mpi/2.19
  module load fftw/3.3.10--hpcx-mpi--2.19--nvhpc--24.5
  module load libxc/7.0.0--nvhpc--24.5
  module load openblas/0.3.26--nvhpc--24.5
  module load netlib-scalapack/2.2.0--hpcx-mpi--2.19--nvhpc--24.5
fi
export LD_LIBRARY_PATH="$QE_ROOT/local-libs/lib:${LD_LIBRARY_PATH:-}"

{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "hostname=$(hostname)"
  echo "executable=$PW_X"
  echo "loader=$QE_ROOT/scripts/load_build_environment.sh"
  echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
  echo "--- slurm vars ---"
  env | grep -E '^SLURM_' | sort || true
  echo "--- module list ---"
  module list
  echo "--- sha256 ---"
  sha256sum "$PW_X"
  echo "--- ldd ---"
  ldd "$PW_X" || true
} > evidence/executable_environment.txt 2>&1

{
  echo "pw_x=$PW_X"
  sha256sum "$PW_X"
} > evidence/executable_checksums.txt

# --------------------------------------------------------------- CONTROL ----
cd "$WORKDIR/control_initial"
set +e
srun --mpi=pmix_v3 -n 1 "$PW_X" -in "$WORKDIR/inputs/control_initial_no_velocity.in" \
     > qe.out 2> qe.err
RC_CONTROL=$?
cd "$WORKDIR"
echo "RC_CONTROL=$RC_CONTROL" > evidence/run_return_codes.txt

# ------------------------------------------------------------------ TEST ----
cd "$WORKDIR/test_nose"
srun --mpi=pmix_v3 -n 1 "$PW_X" -in "$WORKDIR/inputs/test_nose_no_velocity.in" \
     > qe.out 2> qe.err
RC_NOSE=$?
cd "$WORKDIR"

{
  echo "RC_CONTROL=$RC_CONTROL"
  echo "RC_NOSE=$RC_NOSE"
} > evidence/run_return_codes.txt

# -------------------------------------------------------------- analysis ----
PY=python3
if [ -x "$HOME/venvs/mace_py311/bin/python" ]; then
  PY="$HOME/venvs/mace_py311/bin/python"
fi
echo "analysis_python=$PY" >> evidence/run_return_codes.txt
"$PY" scripts/analyze_initial_vs_nose.py \
      > evidence/analysis_stdout.txt 2> evidence/analysis_stderr.txt
RC_ANALYSIS=$?
echo "RC_ANALYSIS=$RC_ANALYSIS" >> evidence/run_return_codes.txt

exit 0
