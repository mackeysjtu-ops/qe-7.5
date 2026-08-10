#!/usr/bin/env bash
#SBATCH --job-name=d8_vcrelax
#SBATCH --account=ICT26_CMSP_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --time=02:00:00
#SBATCH --output=01_vc_relax/slurm-%j.out
#SBATCH --error=01_vc_relax/slurm-%j.err

set -uo pipefail
WORKDIR="/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_vcrelax_nvt_npt_caimd_chain_20260727T203049Z"
QE_ROOT="/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc"
PW_X="$QE_ROOT/install-current-gpu/bin/pw.x"
cd "$WORKDIR"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
module purge
if [ -f "$QE_ROOT/scripts/load_build_environment.sh" ]; then source "$QE_ROOT/scripts/load_build_environment.sh"; fi
export LD_LIBRARY_PATH="$QE_ROOT/local-libs/lib:${LD_LIBRARY_PATH:-}"
PY=python3; [ -x "$HOME/venvs/mace_py311/bin/python" ] && PY="$HOME/venvs/mace_py311/bin/python"
export PYTHONPATH="$WORKDIR/scripts:${PYTHONPATH:-}"
run_qe () { srun --mpi=pmix_v3 -n 1 "$PW_X" -in "$1" > "$2" 2> "${2%.out}.err"; return $?; }


# ---------------- JOB 1: vc-relax ----------------
cd "$WORKDIR/01_vc_relax"
cp "$WORKDIR/inputs/generated/01_vc_relax.in" ./qe.in
mkdir -p tmp
run_qe "$WORKDIR/inputs/generated/01_vc_relax.in" qe.out
RC=$?
cd "$WORKDIR"
echo "RC_VCRELAX=$RC" > evidence/rc_01.txt
"$PY" scripts/extract_final_structure.py --mode vcrelax_gate \
   --job1-out "$WORKDIR/01_vc_relax/qe.out" \
   --out-struct "$WORKDIR/01_vc_relax/final_structure.json" \
   --out-fragment "$WORKDIR/01_vc_relax/final_structure.in.fragment" \
   --out-report "$WORKDIR/01_vc_relax/FINAL_VCRELAX_REPORT.md" \
   --out-pressure "$WORKDIR/01_vc_relax/pressure_history.csv" \
   --out-volume "$WORKDIR/01_vc_relax/volume_history.csv"
GATE=$?
echo "GATE_VCRELAX=$GATE" >> evidence/rc_01.txt
exit $GATE
