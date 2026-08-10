#!/usr/bin/env bash
#SBATCH --job-name=diamond8_nose_qe74_qe75
#SBATCH --account=ICT26_CMSP_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

WORKDIR="/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_ordinary_nose_qe74_vs_qe75_af7ebfc_20260723T175625Z"
QE74_MODULE="quantum-espresso/7.4.1--hpcx-mpi--2.19--nvhpc--24.5-openblas"
QE75_ROOT="/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc"
QE75_EXE="/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc/install-current-gpu/bin/pw.x"

cd "$WORKDIR"
mkdir -p evidence analysis qe74_official/tmp qe75_af7ebfc_no_caimd/tmp

{
  echo "job_id=${SLURM_JOB_ID:-unknown}"
  echo "job_name=${SLURM_JOB_NAME:-unknown}"
  echo "nodelist=${SLURM_JOB_NODELIST:-unknown}"
  echo "submit_dir=${SLURM_SUBMIT_DIR:-unknown}"
  echo "workdir=$WORKDIR"
  hostname
  date -u
  env | grep -E '^(SLURM|CUDA|NVHPC|OMPI|OMP)_' | sort || true
} > evidence/slurm_job_environment.txt 2>&1

sha256sum inputs/qe74.in inputs/qe75_af7ebfc.in pseudo/C.UPF > evidence/input_checksums.txt

module purge
module load profile/chem-phys
module load "$QE74_MODULE"
QE74_EXE="$(command -v pw.x)"
{
  echo "module=$QE74_MODULE"
  echo "executable=$QE74_EXE"
  echo "srun_command=srun --mpi=pmix_v3 -n 1 $QE74_EXE -in $WORKDIR/inputs/qe74.in"
  module list
  timeout 30 "$QE74_EXE" --version || true
  sha256sum "$QE74_EXE"
  ldd "$QE74_EXE" || true
} > evidence/qe74_environment.txt 2>&1

{
  echo "qe74_executable=$QE74_EXE"
  sha256sum "$QE74_EXE"
} > evidence/executable_checksums.txt

cd "$WORKDIR/qe74_official"
set +e
srun --mpi=pmix_v3 -n 1 "$QE74_EXE" -in "$WORKDIR/inputs/qe74.in" > qe.out 2> qe.err
RC_QE74=$?
set -e

cd "$WORKDIR"
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
  echo "executable=$QE75_EXE"
  echo "loader=$QE75_ROOT/scripts/load_build_environment.sh"
  echo "srun_command=srun --mpi=pmix_v3 -n 1 $QE75_EXE -in $WORKDIR/inputs/qe75_af7ebfc.in"
  module list
  timeout 30 "$QE75_EXE" --version || true
  sha256sum "$QE75_EXE"
  ldd "$QE75_EXE" || true
} > evidence/qe75_environment.txt 2>&1
{
  echo "qe75_af7ebfc_executable=$QE75_EXE"
  sha256sum "$QE75_EXE"
} >> evidence/executable_checksums.txt

cd "$WORKDIR/qe75_af7ebfc_no_caimd"
set +e
srun --mpi=pmix_v3 -n 1 "$QE75_EXE" -in "$WORKDIR/inputs/qe75_af7ebfc.in" > qe.out 2> qe.err
RC_QE75_AF7EBFC=$?
set -e

cd "$WORKDIR"
{
  echo "RC_QE74=$RC_QE74"
  echo "RC_QE75_AF7EBFC=$RC_QE75_AF7EBFC"
} > evidence/run_return_codes.txt

python3 scripts/analyze_comparison.py > evidence/analysis_stdout.txt 2> evidence/analysis_stderr.txt

exit 0
