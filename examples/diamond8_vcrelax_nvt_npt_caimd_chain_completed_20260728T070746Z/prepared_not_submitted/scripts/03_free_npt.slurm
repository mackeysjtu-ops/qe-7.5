#!/usr/bin/env bash
#SBATCH --job-name=d8_free_npt
#SBATCH --account=ICT26_CMSP_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --time=04:00:00
#SBATCH --output=03_free_npt/slurm-%j.out
#SBATCH --error=03_free_npt/slurm-%j.err

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


# ---------------- JOB 3: free-cell NPT (CAIMD disabled, ion+cell Nose) ----------------
cd "$WORKDIR"
REFVOL=$("$PY" -c "import json;print(json.load(open('$WORKDIR/01_vc_relax/final_structure.json'))['determinant_bohr3'])")
"$PY" scripts/extract_final_structure.py --mode gen_npt \
   --struct-json "$WORKDIR/01_vc_relax/final_structure.json" \
   --prev-cel "$WORKDIR/02_fixed_nvt/tmp/diamond8_nvt.cel" \
   --prev-pos "$WORKDIR/02_fixed_nvt/tmp/diamond8_nvt.pos" \
   --out-input "$WORKDIR/inputs/generated/03_free_npt.in" \
   --out-handoff "$WORKDIR/evidence/handoff_02_to_03.json" || { echo "handoff 02->03 failed"; exit 1; }
cd "$WORKDIR/03_free_npt"; mkdir -p tmp
run_qe "$WORKDIR/inputs/generated/03_free_npt.in" qe.out
RC=$?
cd "$WORKDIR"; echo "RC_NPT=$RC" > evidence/rc_03.txt
# verify cell Nose active in stdout (temph/fnoseh accepted, tnoseh path)
grep -qi "nose" 03_free_npt/qe.out && echo "cell/ion nose strings present" >> evidence/rc_03.txt
"$PY" scripts/analyze_stage.py --stage npt --dir "$WORKDIR/03_free_npt" \
   --prefix diamond8_npt --out "$WORKDIR/03_free_npt" --dt-fs 1.0 --ref-volume $REFVOL --rc $RC
GATE=$?
"$PY" - "$WORKDIR/03_free_npt" <<'PYEOF'
import json,sys,os
d=json.load(open(os.path.join(sys.argv[1],"analysis_summary.json")))
t=d.get("temperature",{});p=d.get("pressure",{});c=d.get("cell",{})
open(os.path.join(sys.argv[1],"FINAL_NPT_REPORT.md"),"w").write(
 f"# Free-cell NPT report (CAIMD disabled, ion+cell Nose)\n\nGate: {'PASS' if d.get('gate_pass') else 'FAIL'}\n\n"
 f"- steps {d.get('completed_steps')}, JOB DONE {d.get('job_done')}, CAIMD disabled {d.get('caimd_disabled')}, ionic Nose {d.get('nose_ion_active')}\n"
 f"- final-2ps T {t.get('final_2ps_mean')} K; block {t.get('block_uncertainty',{}).get('block_mean')} +/- {t.get('block_uncertainty',{}).get('sem')}\n"
 f"- pressure block mean {p.get('final_block_mean_kbar')} +/- {p.get('block_sem_kbar')} kbar; mean {p.get('mean')} kbar\n"
 f"- volume final {c.get('volume_final_bohr3')} bohr^3; max dev vs refvol {c.get('max_abs_volume_dev_vs_refvol_pct')}%; cond max {c.get('condition_number_max')}\n"
 f"- gate failures: {d.get('gate_failures') or 'none'}\n")
PYEOF
exit $GATE
