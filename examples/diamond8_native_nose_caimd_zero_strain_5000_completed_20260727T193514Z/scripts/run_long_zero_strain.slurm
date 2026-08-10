#!/usr/bin/env bash
#SBATCH --job-name=dia8_caimd_zero_5k
#SBATCH --account=ICT26_CMSP_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# Three 5000-step (10 ps) QE 7.5 af7ebfc trajectories, one executable, native init + Nose.
#   CASE A: caimd=.false.                    (ordinary control)
#   CASE B: caimd=.true., all masks .false.  (no-op boundary)
#   CASE C: caimd=.true., all masks .true.,  step_strain=0  (zero-strain hold)
# No explicit velocities. errexit is deliberately NOT enabled so every case runs.

set -uo pipefail

WORKDIR="/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_native_nose_caimd_zero_strain_5000_20260727T151734Z"
QE_ROOT="/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc"
PW_X="$QE_ROOT/install-current-gpu/bin/pw.x"

cd "$WORKDIR"
mkdir -p evidence analysis case_A_caimd_false/tmp case_B_caimd_noop/tmp case_C_caimd_zero_strain/tmp

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

{
  echo "job_id=${SLURM_JOB_ID:-unknown}"
  echo "job_name=${SLURM_JOB_NAME:-unknown}"
  echo "nodelist=${SLURM_JOB_NODELIST:-unknown}"
  echo "workdir=$WORKDIR"
  echo "omp_num_threads=$OMP_NUM_THREADS"
  hostname
  date -u
  env | grep -E '^(SLURM|CUDA|NVHPC|OMPI|OMP)_' | sort || true
} > evidence/slurm_job_environment.txt 2>&1

sha256sum inputs/case_A_caimd_false_nose_5000.in \
          inputs/case_B_caimd_true_all_false_nose_5000.in \
          inputs/case_C_caimd_true_all_fixed_zero_strain_nose_5000.in \
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
  echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
  echo "--- module list ---"; module list
  echo "--- sha256 ---"; sha256sum "$PW_X"
  echo "--- ldd ---"; ldd "$PW_X" || true
} > evidence/executable_environment.txt 2>&1

{
  echo "pw_x=$PW_X"
  sha256sum "$PW_X"
} > evidence/executable_checksums.txt

run_case () {   # $1=case dir  $2=input basename
  cd "$WORKDIR/$1"
  srun --mpi=pmix_v3 -n 1 "$PW_X" -in "$WORKDIR/inputs/$2" > qe.out 2> qe.err
  local rc=$?
  cd "$WORKDIR"
  return $rc
}

set +e
run_case case_A_caimd_false        case_A_caimd_false_nose_5000.in
RC_A=$?
run_case case_B_caimd_noop         case_B_caimd_true_all_false_nose_5000.in
RC_B=$?
run_case case_C_caimd_zero_strain  case_C_caimd_true_all_fixed_zero_strain_nose_5000.in
RC_C=$?
set +e

{
  echo "RC_A=$RC_A"
  echo "RC_B=$RC_B"
  echo "RC_C=$RC_C"
} > evidence/run_return_codes.txt

# -------------------------------------------------------------- analysis ----
PY=python3
if [ -x "$HOME/venvs/mace_py311/bin/python" ]; then
  PY="$HOME/venvs/mace_py311/bin/python"
fi
echo "analysis_python=$PY" >> evidence/run_return_codes.txt
"$PY" scripts/analyze_long_zero_strain.py > evidence/analysis_stdout.txt 2> evidence/analysis_stderr.txt
RC_ANALYSIS=$?
echo "RC_ANALYSIS=$RC_ANALYSIS" >> evidence/run_return_codes.txt

exit 0
