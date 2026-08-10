#!/usr/bin/env bash
#SBATCH --job-name=dia8_qe74_qe75_novel
#SBATCH --account=ICT26_CMSP_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# Diamond8 native from-scratch 300 K comparison, no explicit velocities.
#   run 1: official unmodified QE 7.4.1 (Leonardo module)
#   run 2: QE 7.5 af7ebfc with caimd=.false.
# Both use ion_temperature='initial', tempw=300.0 (the only 300 K
# initialization mode accepted by both versions for PW vc-md).

set -uo pipefail

WORKDIR="/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_no_velocity_qe74_vs_qe75_caimd_false_300K_20260724T092828Z"
QE74_MODULE="quantum-espresso/7.4.1--hpcx-mpi--2.19--nvhpc--24.5-openblas"
QE75_ROOT="/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc"
QE75_EXE="$QE75_ROOT/install-current-gpu/bin/pw.x"

cd "$WORKDIR"
mkdir -p evidence analysis qe74_official/tmp qe75_af7ebfc_caimd_false/tmp

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

sha256sum inputs/qe74_no_velocity.in \
          inputs/qe75_caimd_false_no_velocity.in \
          pseudo/C.UPF > evidence/input_checksums.txt 2>&1

# ---------------------------------------------------------------- QE 7.4.1 --
module purge
module load profile/chem-phys
module load "$QE74_MODULE"
QE74_EXE="$(command -v pw.x)"
{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "hostname=$(hostname)"
  echo "module=$QE74_MODULE"
  echo "executable=$QE74_EXE"
  echo "srun_command=srun --mpi=pmix_v3 -n 1 $QE74_EXE -in $WORKDIR/inputs/qe74_no_velocity.in"
  echo "--- slurm vars ---"
  env | grep -E '^SLURM_' | sort || true
  echo "--- module list ---"
  module list
  echo "--- pw.x --version ---"
  timeout 30 "$QE74_EXE" --version < /dev/null || echo "(version probe produced no output / non-zero)"
  echo "--- sha256 ---"
  sha256sum "$QE74_EXE"
  echo "--- ldd ---"
  ldd "$QE74_EXE" || true
} > evidence/qe74_environment.txt 2>&1

{
  echo "qe74_official_executable=$QE74_EXE"
  sha256sum "$QE74_EXE"
} > evidence/executable_checksums.txt

cd "$WORKDIR/qe74_official"
set +e
srun --mpi=pmix_v3 -n 1 "$QE74_EXE" -in "$WORKDIR/inputs/qe74_no_velocity.in" \
     > qe.out 2> qe.err
RC_QE74=$?
cd "$WORKDIR"
echo "RC_QE74=$RC_QE74" > evidence/run_return_codes.txt

# ------------------------------------------------------- QE 7.5 af7ebfc -----
module purge
if [ -f "$QE75_ROOT/scripts/load_build_environment.sh" ]; then
  # shellcheck source=/dev/null
  source "$QE75_ROOT/scripts/load_build_environment.sh"
else
  module load profile/chem-phys
  module load nvhpc/24.5
  module load hpcx-mpi/2.19
  module load fftw/3.3.10--hpcx-mpi--2.19--nvhpc--24.5
  module load libxc/7.0.0--nvhpc--24.5
  module load openblas/0.3.26--nvhpc--24.5
  module load netlib-scalapack/2.2.0--hpcx-mpi--2.19--nvhpc--24.5
fi
export LD_LIBRARY_PATH="$QE75_ROOT/local-libs/lib:${LD_LIBRARY_PATH:-}"

{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "hostname=$(hostname)"
  echo "executable=$QE75_EXE"
  echo "loader=$QE75_ROOT/scripts/load_build_environment.sh"
  echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
  echo "srun_command=srun --mpi=pmix_v3 -n 1 $QE75_EXE -in $WORKDIR/inputs/qe75_caimd_false_no_velocity.in"
  echo "--- slurm vars ---"
  env | grep -E '^SLURM_' | sort || true
  echo "--- module list ---"
  module list
  echo "--- pw.x --version ---"
  timeout 30 "$QE75_EXE" --version < /dev/null || echo "(version probe produced no output / non-zero)"
  echo "--- sha256 ---"
  sha256sum "$QE75_EXE"
  echo "--- ldd ---"
  ldd "$QE75_EXE" || true
} > evidence/qe75_environment.txt 2>&1

{
  echo "qe75_af7ebfc_executable=$QE75_EXE"
  sha256sum "$QE75_EXE"
} >> evidence/executable_checksums.txt

cd "$WORKDIR/qe75_af7ebfc_caimd_false"
set +e
srun --mpi=pmix_v3 -n 1 "$QE75_EXE" -in "$WORKDIR/inputs/qe75_caimd_false_no_velocity.in" \
     > qe.out 2> qe.err
RC_QE75=$?

cd "$WORKDIR"
{
  echo "RC_QE74=$RC_QE74"
  echo "RC_QE75_CAIMD_FALSE=$RC_QE75"
} > evidence/run_return_codes.txt

# ------------------------------------------------------------- analysis -----
PY=python3
if [ -x "$HOME/venvs/mace_py311/bin/python" ]; then
  PY="$HOME/venvs/mace_py311/bin/python"
fi
echo "analysis_python=$PY" >> evidence/run_return_codes.txt
set +e
"$PY" scripts/analyze_comparison.py > evidence/analysis_stdout.txt 2> evidence/analysis_stderr.txt
RC_ANALYSIS=$?
echo "RC_ANALYSIS=$RC_ANALYSIS" >> evidence/run_return_codes.txt

exit 0
