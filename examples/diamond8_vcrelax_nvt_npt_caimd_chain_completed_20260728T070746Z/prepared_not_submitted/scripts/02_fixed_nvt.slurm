#!/usr/bin/env bash
#SBATCH --job-name=d8_fixed_nvt
#SBATCH --account=ICT26_CMSP_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --time=04:00:00
#SBATCH --output=02_fixed_nvt/slurm-%j.out
#SBATCH --error=02_fixed_nvt/slurm-%j.err

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


# ---------------- JOB 2: fixed-cell NVT ----------------
cd "$WORKDIR"
"$PY" scripts/extract_final_structure.py --mode gen_nvt \
   --struct-json "$WORKDIR/01_vc_relax/final_structure.json" \
   --out-input "$WORKDIR/inputs/generated/02_fixed_nvt.in" \
   --out-handoff "$WORKDIR/evidence/handoff_01_to_02.json" || { echo "handoff 01->02 failed"; exit 1; }
cd "$WORKDIR/02_fixed_nvt"; mkdir -p tmp
run_qe "$WORKDIR/inputs/generated/02_fixed_nvt.in" qe.out
RC=$?
cd "$WORKDIR"; echo "RC_NVT=$RC" > evidence/rc_02.txt
"$PY" scripts/analyze_stage.py --stage nvt --dir "$WORKDIR/02_fixed_nvt" \
   --prefix diamond8_nvt --out "$WORKDIR/02_fixed_nvt" --dt-fs 1.0 --rc $RC
GATE=$?
cp "$WORKDIR/02_fixed_nvt/analysis_summary.json" "$WORKDIR/02_fixed_nvt/analysis_summary.json" 2>/dev/null || true
# NVT report
"$PY" - "$WORKDIR/02_fixed_nvt" <<'PYEOF'
import json,sys,os
d=json.load(open(os.path.join(sys.argv[1],"analysis_summary.json")))
t=d.get("temperature",{})
open(os.path.join(sys.argv[1],"FINAL_NVT_REPORT.md"),"w").write(
 f"# Fixed-cell NVT report\n\nGate: {'PASS' if d.get('gate_pass') else 'FAIL'}\n\n"
 f"- steps {d.get('completed_steps')}, JOB DONE {d.get('job_done')}, ionic Nose {d.get('nose_ion_active')}\n"
 f"- final-2ps T {t.get('final_2ps_mean')} K; block {t.get('block_uncertainty',{}).get('block_mean')} +/- {t.get('block_uncertainty',{}).get('sem')}\n"
 f"- max T {t.get('max')} K; final-quarter slope {t.get('slope_final_quarter')} K/ps\n"
 f"- pressure mean {d.get('pressure',{}).get('mean')} kbar (diagnostic)\n"
 f"- gate failures: {d.get('gate_failures') or 'none'}\n")
# also emit final_structure.json for downstream (cell from JOB1, positions from NVT final .pos)
PYEOF
# structure passthrough for NPT uses JOB1 cell + NVT final positions (done at NPT gen time)
cp "$WORKDIR/01_vc_relax/final_structure.json" "$WORKDIR/02_fixed_nvt/final_structure.json"
exit $GATE
