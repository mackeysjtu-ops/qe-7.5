#!/usr/bin/env python3
"""Run static validation checks for the Ice X CAIMD example tree."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


EXAMPLE_ROOT = Path(__file__).resolve().parents[1]
EVIDENCE_DIR = EXAMPLE_ROOT / "evidence"
VELOCITY_EVIDENCE = EVIDENCE_DIR / "velocity_initialization_validation.json"
REQUIRED_FILES = [
    "README.md",
    "PARAMETER_SELECTION.md",
    "ICE_X_CAIMD_EXAMPLE_REPORT.md",
    "structure/iceX_original.data",
    "structure/iceX_qe_structure.inc",
    "structure/structure_validation.json",
    "structure/structure_validation.md",
    "structure/iceX_velocities_300K.inc",
    "structure/iceX_velocities_300K.json",
    "pseudo/README.md",
    "pseudo/pseudopotential_checksums.txt",
    "scripts/run_qe_gpu.slurm",
    "scripts/run_convergence.slurm",
    "scripts/generate_300K_velocities.py",
    "scripts/analyze_convergence.py",
    "scripts/analyze_caimd_stress_strain.py",
    "scripts/validate_example.py",
    "evidence/velocity_initialization_validation.json",
    "deformation/tension/qe.in",
    "deformation/tension/run.slurm",
    "deformation/tension/README.md",
    "deformation/tension/expected_files.txt",
    "deformation/tension/analyze.sh",
    "deformation/compression/qe.in",
    "deformation/compression/run.slurm",
    "deformation/compression/README.md",
    "deformation/compression/expected_files.txt",
    "deformation/compression/analyze.sh",
    "deformation/shear/qe.in",
    "deformation/shear/run.slurm",
    "deformation/shear/README.md",
    "deformation/shear/expected_files.txt",
    "deformation/shear/analyze.sh",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fail(checks: list[dict[str, object]], name: str, detail: str) -> None:
    checks.append({"name": name, "ok": False, "detail": detail})


def ok(checks: list[dict[str, object]], name: str, detail: str) -> None:
    checks.append({"name": name, "ok": True, "detail": detail})


def check_required_files(checks: list[dict[str, object]]) -> None:
    missing = [path for path in REQUIRED_FILES if not (EXAMPLE_ROOT / path).exists()]
    if missing:
        fail(checks, "required_files", f"missing files: {missing}")
    else:
        ok(checks, "required_files", f"{len(REQUIRED_FILES)} required files present")


def check_structure(checks: list[dict[str, object]]) -> None:
    validation_path = EXAMPLE_ROOT / "structure" / "structure_validation.json"
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    problems: list[str] = []
    if validation["nat"] != 144:
        problems.append(f"nat={validation['nat']}")
    if validation["atom_counts"] != {"H": 96, "O": 48}:
        problems.append(f"counts={validation['atom_counts']}")
    if validation["source_sha256"] != "103ccf1fe76e3b3a737914a4eeddb4016d6fb8b308129516541ddc5de08aa6c6":
        problems.append("unexpected source checksum")
    if validation["qe_structure_sha256"] != sha256(EXAMPLE_ROOT / "structure" / "iceX_qe_structure.inc"):
        problems.append("structure include checksum mismatch")
    if validation["max_roundtrip_error_angstrom"] > 1.0e-10:
        problems.append(f"roundtrip error={validation['max_roundtrip_error_angstrom']}")
    if problems:
        fail(checks, "structure", "; ".join(problems))
    else:
        ok(checks, "structure", "144 atoms, 48 O, 96 H, source/structure hashes and roundtrip validated")


def check_velocities(checks: list[dict[str, object]]) -> None:
    metadata_path = EXAMPLE_ROOT / "structure" / "iceX_velocities_300K.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    problems: list[str] = []
    if abs(metadata["final_temperature_K"] - 300.0) > 1.0e-9:
        problems.append(f"final_temperature_K={metadata['final_temperature_K']}")
    if metadata["ndof"] != 429:
        problems.append(f"ndof={metadata['ndof']}")
    if metadata["velocity_include_sha256"] != sha256(EXAMPLE_ROOT / "structure" / "iceX_velocities_300K.inc"):
        problems.append("velocity include checksum mismatch")
    if metadata["max_abs_com_velocity_after_rescale"] > 1.0e-16:
        problems.append(f"COM residual={metadata['max_abs_com_velocity_after_rescale']}")
    if problems:
        fail(checks, "velocities", "; ".join(problems))
    else:
        ok(checks, "velocities", "deterministic 300 K velocities, ndof=429, COM removed")


def active_caimd_input(path: Path) -> bool:
    text = path.read_text(encoding="utf-8")
    return bool(re.search(r"caimd\s*=\s*\.true\.", text, re.IGNORECASE) and re.search(r"caimd_mask\([0-9]+\).*\.true\.", text, re.IGNORECASE))


def check_qe_inputs(checks: list[dict[str, object]]) -> None:
    problems: list[str] = []
    for input_path in sorted(EXAMPLE_ROOT.glob("**/qe.in")):
        text = input_path.read_text(encoding="utf-8")
        rel = input_path.relative_to(EXAMPLE_ROOT)
        variable_cell = bool(re.search(r"calculation\s*=\s*'vc-(md|relax)'", text, re.IGNORECASE))
        if variable_cell and "press = 1000.0" not in text:
            problems.append(f"{rel}: missing press=1000.0 kbar")
        if active_caimd_input(input_path):
            required = [
                r"calculation\s*=\s*'vc-md'",
                r"restart_mode\s*=\s*'from_scratch'",
                r"ion_temperature\s*=\s*'not_controlled'",
                r"cell_temperature\s*=\s*'not_controlled'",
                r"cell_dynamics\s*=\s*'pr'",
                r"cell_dofree\s*=\s*'all'",
            ]
            for pattern in required:
                if not re.search(pattern, text, re.IGNORECASE):
                    problems.append(f"{rel}: active CAIMD missing {pattern}")
            if re.search(r"^\s*[A-Za-z]{1,2}\s+[-+0-9.]+\s+[-+0-9.]+\s+[-+0-9.]+\s+0\s+0\s+0", text, re.MULTILINE):
                problems.append(f"{rel}: fixed-coordinate flags found in active CAIMD")
            if "ATOMIC_VELOCITIES { a.u }" not in text:
                problems.append(f"{rel}: active CAIMD missing explicit velocities")
    if problems:
        fail(checks, "qe_inputs", "; ".join(problems))
    else:
        ok(checks, "qe_inputs", "all QE inputs pass static CAIMD boundary and pressure checks")


def number(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def require_number(
    problems: list[str],
    payload: dict[str, object],
    key: str,
    label: str,
) -> float | None:
    value = number(payload.get(key))
    if value is None:
        problems.append(f"{label}: missing numeric {key}")
    return value


def check_velocity_initialization_evidence(checks: list[dict[str, object]]) -> None:
    if not VELOCITY_EVIDENCE.exists():
        ok(checks, "velocity_initialization_evidence", "no runtime velocity-initialization evidence present")
        return

    evidence = json.loads(VELOCITY_EVIDENCE.read_text(encoding="utf-8"))
    generic = evidence.get("generic_current_harness")
    problems: list[str] = []
    if not isinstance(generic, dict):
        problems.append("missing generic_current_harness")
    else:
        current = generic.get("current_vcmd_explicit")
        if not isinstance(current, dict):
            problems.append("missing current_vcmd_explicit")
        else:
            input_speed = require_number(problems, current, "input_velocity_max_abs", "current_vcmd_explicit")
            ratd = require_number(problems, current, "ratd_max_abs", "current_vcmd_explicit")
            error = require_number(problems, current, "input_to_initialized_max_abs", "current_vcmd_explicit")
            if input_speed is not None and input_speed <= 1.0e-12:
                problems.append("current_vcmd_explicit: input velocity is effectively zero")
            if ratd is not None and ratd <= 1.0e-12:
                problems.append("current_vcmd_explicit: initialized PR ratd is near zero")
            if error is not None and error > 1.0e-10:
                problems.append(f"current_vcmd_explicit: input/init mismatch {error}")

        nocard = generic.get("current_vcmd_nocard")
        if not isinstance(nocard, dict):
            problems.append("missing current_vcmd_nocard")
        else:
            ratd = require_number(problems, nocard, "ratd_max_abs", "current_vcmd_nocard")
            if ratd is not None and ratd > 1.0e-15:
                problems.append(f"current_vcmd_nocard: expected zero ratd, got {ratd}")

        baseline = generic.get("baseline_vcmd_explicit")
        if not isinstance(baseline, dict):
            problems.append("missing baseline_vcmd_explicit")
        else:
            input_speed = require_number(problems, baseline, "input_velocity_max_abs", "baseline_vcmd_explicit")
            ratd = require_number(problems, baseline, "ratd_max_abs", "baseline_vcmd_explicit")
            if input_speed is not None and input_speed <= 1.0e-12:
                problems.append("baseline_vcmd_explicit: input velocity is effectively zero")
            if ratd is not None and ratd > 1.0e-15:
                problems.append(f"baseline_vcmd_explicit: expected untouched PR ratd zero, got {ratd}")

    for case in evidence.get("iceX_runtime_cases", []):
        if not isinstance(case, dict):
            problems.append("iceX_runtime_cases: non-object case")
            continue
        name = str(case.get("name", "<unnamed>"))
        target = number(case.get("target_temperature_K")) or 300.0
        initialized = require_number(problems, case, "initialized_temperature_K", name)
        if initialized is not None and not (0.8 * target <= initialized <= 1.2 * target):
            problems.append(f"{name}: initialized_temperature_K={initialized} outside 20% of {target}")

    if problems:
        fail(checks, "velocity_initialization_evidence", "; ".join(problems))
    else:
        ok(checks, "velocity_initialization_evidence", "runtime evidence confirms explicit PR velocity initialization")


def check_slurm(checks: list[dict[str, object]]) -> None:
    problems: list[str] = []
    for rel in ("deformation/tension/run.slurm", "deformation/compression/run.slurm", "deformation/shear/run.slurm"):
        text = (EXAMPLE_ROOT / rel).read_text(encoding="utf-8")
        for snippet in (
            "#SBATCH --account=ICT26_CMSP_0",
            "#SBATCH --partition=boost_usr_prod",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=1",
            "#SBATCH --cpus-per-task=8",
            "#SBATCH --gres=gpu:1",
            "#SBATCH --time=24:00:00",
            "srun --mpi=pmix_v3 \"$PW\" -in qe.in",
            "nvidia-smi --query-gpu",
            "EXPECTED_PW_SHA256=\"c9588df640013f098a33e33e1eec6e1f0820ed53f7488efffb0a83ec7ed25b48\"",
        ):
            if snippet not in text:
                problems.append(f"{rel}: missing {snippet}")
    if problems:
        fail(checks, "slurm", "; ".join(problems))
    else:
        ok(checks, "slurm", "manual deformation Slurm scripts have required Leonardo resources and checks")


def check_pseudo(checks: list[dict[str, object]]) -> None:
    regular_upfs = sorted(
        path.name for path in (EXAMPLE_ROOT / "pseudo").glob("*.upf") if path.exists() and not path.is_symlink()
    )
    symlink_upfs = sorted(path.name for path in (EXAMPLE_ROOT / "pseudo").glob("*.upf") if path.is_symlink())
    if regular_upfs:
        fail(checks, "pseudo", f"regular pseudopotential binaries should not be in Git example: {regular_upfs}")
    else:
        ok(checks, "pseudo", f"no regular UPF binaries present; symlinks allowed/documented: {symlink_upfs}")


def write_report(checks: list[dict[str, object]]) -> int:
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    failed = [check for check in checks if not check["ok"]]
    payload = {"ok": not failed, "checks": checks}
    (EVIDENCE_DIR / "local_static_validation.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with (EVIDENCE_DIR / "local_static_validation.md").open("w", encoding="utf-8") as handle:
        handle.write("# Local static validation\n\n")
        handle.write(f"Overall: {'PASS' if not failed else 'FAIL'}\n\n")
        for check in checks:
            handle.write(f"- {'PASS' if check['ok'] else 'FAIL'} `{check['name']}`: {check['detail']}\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if not failed else 1


def main() -> int:
    checks: list[dict[str, object]] = []
    check_required_files(checks)
    check_structure(checks)
    check_velocities(checks)
    check_qe_inputs(checks)
    check_velocity_initialization_evidence(checks)
    check_slurm(checks)
    check_pseudo(checks)
    return write_report(checks)


if __name__ == "__main__":
    raise SystemExit(main())
