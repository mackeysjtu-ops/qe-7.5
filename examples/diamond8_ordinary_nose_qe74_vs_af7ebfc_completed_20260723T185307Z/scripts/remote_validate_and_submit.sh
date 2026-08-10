#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_ordinary_nose_qe74_vs_qe75_af7ebfc_20260723T175625Z"
QE74_MODULE="quantum-espresso/7.4.1--hpcx-mpi--2.19--nvhpc--24.5-openblas"
QE75_EXE="/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc/install-current-gpu/bin/pw.x"
QE75_ROOT="/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc"
JOB_SCRIPT="$WORKDIR/scripts/run_comparison.slurm"
SUBMISSION_JSON="$WORKDIR/evidence/SUBMISSION_RECORD.json"

cd "$WORKDIR"
mkdir -p evidence

squeue -u "$USER" -h -o "%i|%j|%T|%Z" > evidence/pre_submission_squeue.txt || true
if grep -E '\|diamond8_nose_qe74_qe75\|' evidence/pre_submission_squeue.txt >/dev/null 2>&1; then
  echo "BLOCKED_DUPLICATE_JOB"
  exit 45
fi

module purge
module load profile/chem-phys
if ! module load "$QE74_MODULE"; then
  module avail quantum-espresso > evidence/module_avail_quantum_espresso.txt 2>&1 || true
  echo "BLOCKED_MODULE"
  exit 42
fi
QE74_EXE="$(command -v pw.x)"
if [ -z "$QE74_EXE" ] || [ ! -x "$QE74_EXE" ]; then
  echo "BLOCKED_MODULE"
  exit 42
fi
timeout 30 "$QE74_EXE" --version > evidence/qe74_version_pre_submit.txt 2>&1 || true
sha256sum "$QE74_EXE" > evidence/qe74_executable_pre_submit.sha256

if [ ! -x "$QE75_EXE" ]; then
  echo "BLOCKED_EXECUTABLE"
  exit 43
fi
sha256sum "$QE75_EXE" > evidence/qe75_executable_pre_submit.sha256
module purge
if [ -f "$QE75_ROOT/scripts/load_build_environment.sh" ]; then
  echo "$QE75_ROOT/scripts/load_build_environment.sh" > evidence/qe75_loader_path.txt
fi
module load profile/chem-phys
module load nvhpc/24.5
module load hpcx-mpi/2.19
module load fftw/3.3.10--hpcx-mpi--2.19--nvhpc--24.5
module load libxc/7.0.0--nvhpc--24.5
module load openblas/0.3.26--nvhpc--24.5
module load netlib-scalapack/2.2.0--hpcx-mpi--2.19--nvhpc--24.5
export LD_LIBRARY_PATH="$QE75_ROOT/local-libs/lib:${LD_LIBRARY_PATH:-}"
timeout 30 "$QE75_EXE" --version > evidence/qe75_version_pre_submit.txt 2>&1 || true
if ! timeout 60 ldd "$QE75_EXE" > evidence/qe75_ldd_pre_submit.txt 2>&1; then
  echo "BLOCKED_EXECUTABLE"
  exit 43
fi
if grep -i "not found" evidence/qe75_ldd_pre_submit.txt >/dev/null 2>&1; then
  echo "BLOCKED_EXECUTABLE"
  exit 43
fi

python3 scripts/pre_submit_validate.py > evidence/pre_submit_validate_stdout.json
python3 -m py_compile scripts/analyze_comparison.py scripts/pre_submit_validate.py
bash -n "$JOB_SCRIPT"
sha256sum inputs/qe74.in inputs/qe75_af7ebfc.in pseudo/C.UPF > evidence/input_checksums.txt

JOB_ID="$(sbatch --parsable "$JOB_SCRIPT")"
scontrol show job "$JOB_ID" > evidence/submission_scontrol_show_job.txt || true

python3 - "$JOB_ID" "$QE74_EXE" <<'PY'
import json
import pathlib
import sys

job_id = sys.argv[1]
qe74_exe = sys.argv[2]
root = pathlib.Path("/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_ordinary_nose_qe74_vs_qe75_af7ebfc_20260723T175625Z")
record = {
    "status": "SUBMITTED",
    "job_id": job_id,
    "workdir": str(root),
    "slurm_script": str(root / "scripts" / "run_comparison.slurm"),
    "qe74_module": "quantum-espresso/7.4.1--hpcx-mpi--2.19--nvhpc--24.5-openblas",
    "qe74_executable": qe74_exe,
    "qe74_executable_sha256": (root / "evidence" / "qe74_executable_pre_submit.sha256").read_text().split()[0],
    "qe75_af7ebfc_executable": "/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc/install-current-gpu/bin/pw.x",
    "qe75_af7ebfc_executable_sha256": (root / "evidence" / "qe75_executable_pre_submit.sha256").read_text().split()[0],
    "input_checksums": (root / "evidence" / "input_checksums.txt").read_text().splitlines(),
    "pseudo_sha256": [line.split()[0] for line in (root / "evidence" / "input_checksums.txt").read_text().splitlines() if line.endswith("pseudo/C.UPF")][0],
    "expected_final_report": str(root / "FINAL_COMPARISON_REPORT.md"),
    "expected_summary": str(root / "analysis" / "comparison_summary.json"),
    "future_inspection_commands": [
        f"sacct -j {job_id} --format=JobID,JobName%32,State,ExitCode,Elapsed,NodeList%30",
        f"cat {root}/FINAL_COMPARISON_REPORT.md",
        f"cat {root}/analysis/comparison_summary.json",
    ],
    "asynchronous_rule": "sbatch accepted only; no result inspection performed",
}
(root / "evidence" / "SUBMISSION_RECORD.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
print("SUBMISSION_RECORD_BEGIN")
print(json.dumps(record, indent=2, sort_keys=True))
print("SUBMISSION_RECORD_END")
PY
