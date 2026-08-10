#!/usr/bin/env python3
"""Create the Leonardo home/workspace Diamond64 linear-defgrad campaign."""

from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
from pathlib import Path

TS = "20260801T192312Z"
HOME = Path("/leonardo/home/userexternal/mcheng00")
SOFTWARE_ROOT = HOME / "software" / f"qe-7.5-caimd-linear-defgrad-af7ebfc-{TS}"
CAMPAIGN_ROOT = HOME / "workspace" / "Diamond" / f"diamond64_3000K_111_linear_defgrad_reduced_cost_{TS}"
OLD_WORK_ROOT = Path("/leonardo_work/ICT26_CMSP_0/DIAMOND/runs") / f"diamond64_3000K_111_linear_defgrad_reduced_cost_{TS}"
PSEUDO_SOURCE = HOME / "workspace/Diamond/diamond64_caimd_zero_pmxx_pmxz_1pct_20260731T183444Z/pseudo/C.UPF"
STRUCT_SOURCE = HOME / "workspace/Diamond/diamond8_vcrelax_qe75_vs_qe74_ibrav_20260728T205518Z/qe75/analysis/final_structure.json"
OLD_PW = HOME / "software/qe-7.5-caimd-thermostat-af7ebfc/install-current-gpu/bin/pw.x"
NEW_PW = SOFTWARE_ROOT / "install-current-gpu/bin/pw.x"
TOOL_SOURCE = Path(__file__).with_name("leonardo_campaign_tools.py")

EXPECTED = {
    "pw_sha": "19524a9fb817f0b723be2cf4ab0989fa3734c42559f150ced036b9a482c7ac1a",
    "old_pw_sha": "89d4fda827d2cc22a3b4c1670e90dc40ce629af42096c785b55c70c355a10e7d",
    "pseudo_sha": "9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639",
    "structure_sha": "f9e251cd3a6ccaa955fa374b41f0129c4db6183332fc3302528c9e9dcaa9b199",
    "source_patch_sha": "ac5006766f0196d8750363bd6d802ec74f29b5ecc161ea143f8951172777fd9b",
    "source_tar_sha": "dfd8166addae10b48b88e9ca04fadc6e5e7cab00954a66aefdcf52be5a7f1eb3",
}


def sha256(path: Path | str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def det3(m):
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def matvec_row(fr, h):
    return [sum(fr[i] * h[i][j] for i in range(3)) for j in range(3)]


def distance_spectrum(cell, frac):
    values = []
    for i in range(len(frac)):
        for j in range(i + 1, len(frac)):
            fi = frac[i][1]
            fj = frac[j][1]
            best = 1.0e99
            for sx in (-1, 0, 1):
                for sy in (-1, 0, 1):
                    for sz in (-1, 0, 1):
                        df = [fj[0] - fi[0] + sx, fj[1] - fi[1] + sy, fj[2] - fi[2] + sz]
                        r = math.sqrt(sum(x * x for x in matvec_row(df, cell)))
                        if 1.0e-8 < r < best:
                            best = r
            values.append(best)
    values.sort()
    return values


def write_executable(path: Path, text: str):
    path.write_text(text)
    os.chmod(path, 0o755)


def slurm(name, job, walltime, body, ntasks=4, cpus=8, gpus=4):
    text = f"""#!/usr/bin/env bash
#SBATCH --job-name={job}
#SBATCH --account=ICT26_CMSP_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks={ntasks}
#SBATCH --cpus-per-task={cpus}
#SBATCH --gpus-per-node={gpus}
#SBATCH --time={walltime}
#SBATCH --output={CAMPAIGN_ROOT}/logs/%x-%j.out
#SBATCH --error={CAMPAIGN_ROOT}/logs/%x-%j.err

set -euo pipefail
export CAMPAIGN_ROOT="{CAMPAIGN_ROOT}"
export SOFTWARE_ROOT="{SOFTWARE_ROOT}"
cd "$CAMPAIGN_ROOT"
source "$SOFTWARE_ROOT/scripts/load_build_environment.sh"
export OMP_NUM_THREADS="${{SLURM_CPUS_PER_TASK:-{cpus}}}"
export OMPI_MCA_coll_hcoll_enable=0
PY="$HOME/venvs/mace_py311/bin/python"
if [ ! -x "$PY" ]; then PY=python3; fi
{body}
"""
    write_executable(CAMPAIGN_ROOT / "scripts" / name, text)


def make_structure():
    src = json.loads(STRUCT_SOURCE.read_text())
    a = src["celldm1_bohr"]
    cubic = [[2 * a, 0.0, 0.0], [0.0, 2 * a, 0.0], [0.0, 0.0, 2 * a]]
    x = [1 / math.sqrt(6), 1 / math.sqrt(6), -2 / math.sqrt(6)]
    y = [-1 / math.sqrt(2), 1 / math.sqrt(2), 0.0]
    z = [1 / math.sqrt(3), 1 / math.sqrt(3), 1 / math.sqrt(3)]
    basis = [[x[0], y[0], z[0]], [x[1], y[1], z[1]], [x[2], y[2], z[2]]]
    oriented = [[sum(cubic[i][k] * basis[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    frac = []
    for ix in range(2):
        for iy in range(2):
            for iz in range(2):
                for sym, pos in src["fractional_positions"]:
                    frac.append([sym, [(pos[0] + ix) / 2.0, (pos[1] + iy) / 2.0, (pos[2] + iz) / 2.0]])
    spec_c = distance_spectrum(cubic, frac)
    spec_o = distance_spectrum(oriented, frac)
    report = {
        "status": "ORIENTATION_PASS",
        "basis_original_cubic_columns": {"local_x_11m2": x, "local_y_m110": y, "local_z_111": z},
        "dot_products": {"x_y": dot(x, y), "x_z": dot(x, z), "y_z": dot(y, z)},
        "cross_x_y": cross(x, y),
        "det_basis": det3(basis),
        "cell_cubic_bohr_rows": cubic,
        "cell_oriented_bohr_rows": oriented,
        "volume_cubic_bohr3": det3(cubic),
        "volume_oriented_bohr3": det3(oriented),
        "volume_abs_difference": abs(det3(cubic) - det3(oriented)),
        "distance_spectrum_first20_bohr": spec_o[:20],
        "distance_spectrum_max_abs_delta_bohr": max(abs(a - b) for a, b in zip(spec_c, spec_o)),
        "miller_axis_checks": {
            "local_z_parallel_111_local_components": [dot([1, 1, 1], x), dot([1, 1, 1], y), dot([1, 1, 1], z)],
            "local_x_parallel_11m2_local_components": [dot([1, 1, -2], x), dot([1, 1, -2], y), dot([1, 1, -2], z)],
            "local_y_parallel_m110_local_components": [dot([-1, 1, 0], x), dot([-1, 1, 0], y), dot([-1, 1, 0], z)],
        },
        "n_atoms": len(frac),
        "source_structure_sha256": EXPECTED["structure_sha"],
    }
    if abs(report["det_basis"] - 1.0) > 1.0e-12 or report["distance_spectrum_max_abs_delta_bohr"] > 1.0e-9 or len(frac) != 64:
        report["status"] = "ORIENTATION_FAIL"
        raise SystemExit(json.dumps(report, indent=2))
    structure = {
        "cell_bohr_rows": oriented,
        "fractional_positions": frac,
        "n_atoms": 64,
        "source": "diamond8 relaxed 2x2x2, oriented local z=[111], local x=[11-2], local y=[-110]",
        "determinant_bohr3": det3(oriented),
        "orientation_report": "orientation_report.json",
    }
    (CAMPAIGN_ROOT / "source_structure/diamond64_oriented_111.json").write_text(json.dumps(structure, indent=2) + "\n")
    (CAMPAIGN_ROOT / "source_structure/orientation_report.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = ["CELL_PARAMETERS bohr"]
    for row in oriented:
        lines.append("  " + "  ".join(f"{v:.12f}" for v in row))
    lines.append("ATOMIC_POSITIONS crystal")
    for sym, pos in frac:
        lines.append(f"{sym}  {pos[0]:.12f}  {pos[1]:.12f}  {pos[2]:.12f}")
    (CAMPAIGN_ROOT / "source_structure/diamond64_oriented_111.in.fragment").write_text("\n".join(lines) + "\n")


def main():
    required = {
        "new_pw": (NEW_PW, EXPECTED["pw_sha"]),
        "old_pw": (OLD_PW, EXPECTED["old_pw_sha"]),
        "pseudo": (PSEUDO_SOURCE, EXPECTED["pseudo_sha"]),
        "structure": (STRUCT_SOURCE, EXPECTED["structure_sha"]),
        "source_patch": (SOFTWARE_ROOT / "src-current/caimd_linear_defgrad_phase1/evidence/linear_defgrad_phase1.patch", EXPECTED["source_patch_sha"]),
    }
    for label, (path, expected) in required.items():
        if not path.exists():
            raise SystemExit(f"missing {label}: {path}")
        actual = sha256(path)
        if actual != expected:
            raise SystemExit(f"{label} SHA mismatch: {actual} != {expected}")

    if CAMPAIGN_ROOT.exists() and any(CAMPAIGN_ROOT.iterdir()):
        raise SystemExit(f"campaign root exists and is non-empty: {CAMPAIGN_ROOT}")
    for directory in [
        CAMPAIGN_ROOT,
        CAMPAIGN_ROOT / "scripts",
        CAMPAIGN_ROOT / "inputs/generated",
        CAMPAIGN_ROOT / "pseudo",
        CAMPAIGN_ROOT / "source_structure",
        CAMPAIGN_ROOT / "validation/runs",
        CAMPAIGN_ROOT / "validation/benchmarks",
        CAMPAIGN_ROOT / "validation/analysis",
        CAMPAIGN_ROOT / "equilibration",
        CAMPAIGN_ROOT / "prepare_loading_start",
        CAMPAIGN_ROOT / "compression_111",
        CAMPAIGN_ROOT / "shear_111_11m2",
        CAMPAIGN_ROOT / "logs",
        CAMPAIGN_ROOT / "evidence",
    ]:
        directory.mkdir(parents=True, exist_ok=True)

    shutil.copy2(PSEUDO_SOURCE, CAMPAIGN_ROOT / "pseudo/C.UPF")
    shutil.copy2(STRUCT_SOURCE, CAMPAIGN_ROOT / "source_structure/diamond8_relaxed_source.json")
    make_structure()

    tool_text = TOOL_SOURCE.read_text()
    tool_text = tool_text.replace("__CAMPAIGN_ROOT__", str(CAMPAIGN_ROOT))
    tool_text = tool_text.replace("__SOFTWARE_ROOT__", str(SOFTWARE_ROOT))
    tool_text = tool_text.replace("__PW_SHA__", EXPECTED["pw_sha"])
    write_executable(CAMPAIGN_ROOT / "scripts/campaign_tools.py", tool_text)

    definition = {
        "task": "VALIDATE_LINEAR_DEFGRAD_BUILD_LEONARDO_AND_SUBMIT_GATED_3000K_111_LOADING_CAMPAIGN",
        "timestamp_utc": TS,
        "label": "REDUCED_COST_GAMMA_ONLY_PILOT",
        "software_root": str(SOFTWARE_ROOT),
        "campaign_root": str(CAMPAIGN_ROOT),
        "superseded_initial_work_root": str(OLD_WORK_ROOT),
        "pw_exe": str(NEW_PW),
        "pw_sha256": EXPECTED["pw_sha"],
        "old_pw_exe": str(OLD_PW),
        "old_pw_sha256": EXPECTED["old_pw_sha"],
        "pseudo": str(CAMPAIGN_ROOT / "pseudo/C.UPF"),
        "pseudo_sha256": EXPECTED["pseudo_sha"],
        "source_json": str(STRUCT_SOURCE),
        "source_json_sha256": EXPECTED["structure_sha"],
        "source_patch_sha256": EXPECTED["source_patch_sha"],
        "dt_fs": 2.0,
        "dt_qe": 41.3413733366,
        "caimd_stride": 10,
        "compression": {"D33_per_step": -2.0e-6, "nstep": 75000, "events": 7500, "final_F33": 0.85, "Ezz": -0.13875},
        "shear": {"D13_per_step": 2.0e-6, "nstep": 100000, "events": 10000, "final_F13": 0.20, "detF": 1.0, "Exz": 0.10, "Ezx": 0.10, "Ezz": 0.02},
        "home_workspace_layout_requested_by_user": True,
    }
    (CAMPAIGN_ROOT / "CAMPAIGN_DEFINITION.json").write_text(json.dumps(definition, indent=2) + "\n")
    (CAMPAIGN_ROOT / "README.md").write_text(
        "# Diamond64 3000 K [111] linear_defgrad CAIMD campaign\n\n"
        f"Software root: `{SOFTWARE_ROOT}`\n\n"
        f"Campaign root: `{CAMPAIGN_ROOT}`\n\n"
        "This uses the home/software plus home/workspace layout requested on 2026-08-01. "
        "The earlier `/leonardo_work/ICT26_CMSP_0` placeholder is not used for outputs.\n"
    )
    if OLD_WORK_ROOT.exists():
        (OLD_WORK_ROOT / "MOVED_TO_HOME_WORKSPACE.txt").write_text(
            f"Current campaign artifacts are under {CAMPAIGN_ROOT}\n"
            "No scientific results were moved or deleted from /leonardo_work by this setup step.\n"
        )

    slurm("validate_and_benchmark.slurm", "d64_ldg_validate", "06:00:00", '"$PY" scripts/campaign_tools.py validate')
    slurm("run_eq_0300.slurm", "d64_eq0300", "04:00:00", '"$PY" scripts/campaign_tools.py run-eq T0300 300.0')
    slurm("run_eq_1000.slurm", "d64_eq1000", "04:00:00", '"$PY" scripts/campaign_tools.py run-eq T1000 1000.0 --prev T0300')
    slurm("run_eq_2000.slurm", "d64_eq2000", "04:00:00", '"$PY" scripts/campaign_tools.py run-eq T2000 2000.0 --prev T1000')
    slurm("run_eq_3000.slurm", "d64_eq3000", "04:00:00", '"$PY" scripts/campaign_tools.py run-eq T3000 3000.0 --prev T2000')
    slurm("prepare_loading_start.slurm", "d64_prepare_load", "00:30:00", '"$PY" scripts/campaign_tools.py prepare-loading-start', ntasks=1, cpus=4, gpus=1)
    slurm("run_compression_111.slurm", "d64_comp111", "24:00:00", '"$PY" scripts/campaign_tools.py run-loading compression')
    slurm("run_shear_111_11m2.slurm", "d64_shear111", "24:00:00", '"$PY" scripts/campaign_tools.py run-loading shear')

    submit = f"""#!/usr/bin/env bash
set -euo pipefail
CAMPAIGN_ROOT="{CAMPAIGN_ROOT}"
cd "$CAMPAIGN_ROOT"
PY="$HOME/venvs/mace_py311/bin/python"
if [ ! -x "$PY" ]; then PY=python3; fi
for script in scripts/*.slurm scripts/submit_campaign.sh; do bash -n "$script"; done
"$PY" -m py_compile scripts/campaign_tools.py
"$PY" - <<'PRE'
import hashlib, json, pathlib
root = pathlib.Path("{CAMPAIGN_ROOT}")
def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()
checks = {{
    "pw_sha": sha("{NEW_PW}"),
    "pseudo_sha": sha(root / "pseudo/C.UPF"),
    "orientation_status": json.loads((root / "source_structure/orientation_report.json").read_text())["status"],
}}
if checks["pw_sha"] != "{EXPECTED['pw_sha']}":
    raise SystemExit("new pw.x SHA mismatch")
if checks["pseudo_sha"] != "{EXPECTED['pseudo_sha']}":
    raise SystemExit("C.UPF SHA mismatch")
if checks["orientation_status"] != "ORIENTATION_PASS":
    raise SystemExit("orientation validation failed")
definition = json.loads((root / "CAMPAIGN_DEFINITION.json").read_text())
assert definition["compression"]["events"] == 7500
assert abs(definition["compression"]["final_F33"] - 0.85) < 1.0e-15
assert definition["shear"]["events"] == 10000
assert abs(definition["shear"]["final_F13"] - 0.20) < 1.0e-15
(root / "PRE_SUBMISSION_VALIDATION.json").write_text(json.dumps({{"status": "PRE_SUBMISSION_PASS_DEFER_RUNTIME_TO_VALIDATION_JOB", "checks": checks}}, indent=2) + "\\n")
PRE
mkdir -p evidence logs
J_VALIDATION=$(sbatch --parsable scripts/validate_and_benchmark.slurm)
J_EQ300=$(sbatch --parsable --dependency=afterok:$J_VALIDATION scripts/run_eq_0300.slurm)
J_EQ1000=$(sbatch --parsable --dependency=afterok:$J_EQ300 scripts/run_eq_1000.slurm)
J_EQ2000=$(sbatch --parsable --dependency=afterok:$J_EQ1000 scripts/run_eq_2000.slurm)
J_EQ3000=$(sbatch --parsable --dependency=afterok:$J_EQ2000 scripts/run_eq_3000.slurm)
J_PREPARE=$(sbatch --parsable --dependency=afterok:$J_EQ3000 scripts/prepare_loading_start.slurm)
J_COMPRESSION=$(sbatch --parsable --dependency=afterok:$J_PREPARE scripts/run_compression_111.slurm)
J_SHEAR=$(sbatch --parsable --dependency=afterok:$J_PREPARE scripts/run_shear_111_11m2.slurm)
cat > evidence/SUBMISSION_JOB_IDS.env <<IDS
J_VALIDATION=$J_VALIDATION
J_EQ300=$J_EQ300
J_EQ1000=$J_EQ1000
J_EQ2000=$J_EQ2000
J_EQ3000=$J_EQ3000
J_PREPARE=$J_PREPARE
J_COMPRESSION=$J_COMPRESSION
J_SHEAR=$J_SHEAR
IDS
"$PY" - <<POST
import json, pathlib, time
root = pathlib.Path("$CAMPAIGN_ROOT")
ids = dict(line.strip().split("=", 1) for line in (root / "evidence/SUBMISSION_JOB_IDS.env").read_text().splitlines() if line.strip())
manifest = {{
    "status": "SUBMITTED",
    "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "campaign_root": str(root),
    "software_root": "{SOFTWARE_ROOT}",
    "job_ids": ids,
    "dependencies": {{
        "J_EQ300": "afterok:J_VALIDATION",
        "J_EQ1000": "afterok:J_EQ300",
        "J_EQ2000": "afterok:J_EQ1000",
        "J_EQ3000": "afterok:J_EQ2000",
        "J_PREPARE": "afterok:J_EQ3000",
        "J_COMPRESSION": "afterok:J_PREPARE",
        "J_SHEAR": "afterok:J_PREPARE"
    }},
    "validation_artifact": "validation/VALIDATION_PASS",
    "selected_settings": "validation/selected_settings.env"
}}
(root / "evidence/SUBMISSION_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\\n")
print(json.dumps(manifest, indent=2))
POST
    """
    write_executable(CAMPAIGN_ROOT / "scripts/submit_campaign.sh", submit)

    manifest = []
    for path in sorted(CAMPAIGN_ROOT.rglob("*")):
        if path.is_file():
            manifest.append(f"{sha256(path)}  {path.relative_to(CAMPAIGN_ROOT)}")
    (CAMPAIGN_ROOT / "evidence/initial_file_manifest.sha256").write_text("\n".join(manifest) + "\n")
    print(json.dumps({"status": "CAMPAIGN_ARTIFACTS_CREATED", "campaign_root": str(CAMPAIGN_ROOT), "software_root": str(SOFTWARE_ROOT), "files": len(manifest)}, indent=2))


if __name__ == "__main__":
    main()
