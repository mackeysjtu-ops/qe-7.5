#!/usr/bin/env bash
#SBATCH --job-name=d8_caimd_zn
#SBATCH --account=ICT26_CMSP_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --time=04:00:00
#SBATCH --output=04_caimd/slurm-%j.out
#SBATCH --error=04_caimd/slurm-%j.err

set -uo pipefail
WORKDIR="/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_vcrelax_nvt_npt_caimd_chain_20260728T120614Z"
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


# ---------------- JOB 4: CAIMD zero-strain (4A) then small nonzero (4B) ----------------
cd "$WORKDIR"
# Stage 4A input from JOB3 final cell/pos
"$PY" scripts/extract_final_structure.py --mode gen_caimd_zero \
   --prev-cel "$WORKDIR/03_free_npt/tmp/diamond8_npt.cel" \
   --prev-pos "$WORKDIR/03_free_npt/tmp/diamond8_npt.pos" \
   --out-input "$WORKDIR/inputs/generated/04A_caimd_zero.in" \
   --out-handoff "$WORKDIR/evidence/handoff_03_to_04.json" || { echo "handoff 03->04 failed"; exit 1; }
cd "$WORKDIR/04_caimd/zero_strain"; mkdir -p tmp
run_qe "$WORKDIR/inputs/generated/04A_caimd_zero.in" qe.out
RC=$?
cd "$WORKDIR"; echo "RC_CAIMD_ZERO=$RC" > evidence/rc_04.txt
"$PY" scripts/analyze_stage.py --stage caimd_zero --dir "$WORKDIR/04_caimd/zero_strain" \
   --prefix diamond8_caimd_zero --out "$WORKDIR/04_caimd/zero_strain" --dt-fs 1.0 --rc $RC
GATE_A=$?
"$PY" - "$WORKDIR/04_caimd/zero_strain" <<'PYEOF'
import json,sys,os
d=json.load(open(os.path.join(sys.argv[1],"analysis_summary.json")))
t=d.get("temperature",{});c=d.get("caimd",{})
open(os.path.join(sys.argv[1],"FINAL_ZERO_STRAIN_REPORT.md"),"w").write(
 f"# CAIMD zero-strain (Stage 4A) report\n\nGate: {'PASS' if d.get('gate_pass') else 'FAIL'}\n\n"
 f"- steps {d.get('completed_steps')}, JOB DONE {d.get('job_done')}, events {c.get('n_event_flags_true')}\n"
 f"- max strain err {c.get('max_abs_nominal_minus_actual')}; max vjump {c.get('max_velocity_jump_max')}; max skew {c.get('max_skew_norm')}\n"
 f"- det range {c.get('determinant_range')}; final-1.5ps T {t.get('final_1p5ps_mean')} K\n"
 f"- gate failures: {d.get('gate_failures') or 'none'}\n")
PYEOF
if [ $GATE_A -ne 0 ]; then
  echo "GATE_4A=FAIL - skipping 4B" >> evidence/rc_04.txt
  "$PY" scripts/build_overall_report.py "$WORKDIR" || true
  exit 1
fi
echo "GATE_4A=PASS" >> evidence/rc_04.txt
# Stage 4B from 4A final cell/pos
"$PY" scripts/extract_final_structure.py --mode gen_caimd_nonzero \
   --prev-cel "$WORKDIR/04_caimd/zero_strain/tmp/diamond8_caimd_zero.cel" \
   --prev-pos "$WORKDIR/04_caimd/zero_strain/tmp/diamond8_caimd_zero.pos" \
   --out-input "$WORKDIR/inputs/generated/04B_caimd_nonzero.in" \
   --out-handoff "$WORKDIR/evidence/handoff_zero_to_nonzero.json" || { echo "handoff zero->nonzero failed"; exit 1; }
cd "$WORKDIR/04_caimd/small_nonzero_strain"; mkdir -p tmp
run_qe "$WORKDIR/inputs/generated/04B_caimd_nonzero.in" qe.out
RC=$?
cd "$WORKDIR"; echo "RC_CAIMD_NONZERO=$RC" >> evidence/rc_04.txt
"$PY" scripts/analyze_stage.py --stage caimd_nonzero --dir "$WORKDIR/04_caimd/small_nonzero_strain" \
   --prefix diamond8_caimd_nonzero --out "$WORKDIR/04_caimd/small_nonzero_strain" --dt-fs 1.0 --rc $RC
GATE_B=$?
"$PY" - "$WORKDIR/04_caimd/small_nonzero_strain" <<'PYEOF'
import json,sys,os
d=json.load(open(os.path.join(sys.argv[1],"analysis_summary.json")))
t=d.get("temperature",{});c=d.get("caimd",{})
open(os.path.join(sys.argv[1],"FINAL_NONZERO_STRAIN_REPORT.md"),"w").write(
 f"# CAIMD small nonzero-strain (Stage 4B) report\n\nGate: {'PASS' if d.get('gate_pass') else 'FAIL'}\n\n"
 f"- steps {d.get('completed_steps')}, JOB DONE {d.get('job_done')}, events {c.get('stage4B_event_count')} (expected {c.get('stage4B_expected_events')})\n"
 f"- final nominal xx {c.get('stage4B_final_nominal_xx')}; final actual xx {c.get('stage4B_final_actual_xx')}; expected cumulative {c.get('stage4B_expected_cumulative_xx')}\n"
 f"- xx tracking err {c.get('stage4B_xx_tracking_error')}; max vjump {c.get('max_velocity_jump_max')}; max skew {c.get('max_skew_norm')}\n"
 f"- final-1.5ps T {t.get('final_1p5ps_mean')} K; note: {c.get('prompt_estimate_note')}\n"
 f"- gate failures: {d.get('gate_failures') or 'none'}\n")
PYEOF
"$PY" scripts/build_overall_report.py "$WORKDIR" || true
exit 0
