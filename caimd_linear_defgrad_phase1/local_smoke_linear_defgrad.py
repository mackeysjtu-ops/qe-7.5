#!/usr/bin/env python3
"""Local non-scientific smoke tests for CAIMD linear_defgrad input/runtime paths."""

from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "test-suite" / "caimd"))

from caimd_validation import qe_input  # noqa: E402


def floats(parts: list[str]) -> list[float]:
    return [float(part.replace("D", "E").replace("d", "e")) for part in parts]


def parse_defgrad(path: Path) -> list[dict[str, object]]:
    rows = []
    if not path.exists():
        return rows
    for line in path.read_text(errors="ignore").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        values = floats(parts[2:])
        rows.append(
            {
                "ionic_step": int(parts[0]),
                "event_index": int(parts[1]),
                "det_f_before": values[0],
                "det_f_after": values[1],
                "det_a_event": values[2],
                "f_before": values[3:12],
                "f_after": values[12:21],
                "a_event": values[21:30],
                "green_lagrange": values[30:39],
                "velocity_jump_max": values[39],
                "velocity_jump_rms": values[40],
            }
        )
    return rows


def assert_close(name: str, value: float, expected: float, tol: float = 1.0e-12) -> None:
    if abs(value - expected) > tol:
        raise AssertionError(f"{name}: got {value:.16e}, expected {expected:.16e}")


def assert_identity(name: str, values: list[float], tol: float = 1.0e-12) -> None:
    expected = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    for idx, (value, ref) in enumerate(zip(values, expected), 1):
        assert_close(f"{name}[{idx}]", value, ref, tol)


def run_case(pw: Path, work_root: Path, pseudo_dir: Path, name: str, caimd: str, nstep: int,
             expect_success: bool = True) -> dict[str, object]:
    run_dir = work_root / name
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    input_path = run_dir / f"{name}.in"
    input_path.write_text(
        qe_input(
            name,
            "vc-md",
            pseudo_dir,
            nstep=nstep,
            caimd=caimd,
            velocities=True,
            extra_system="nosym=.true., noinv=.true.,",
            dt=5.0,
            wmass=200.0,
        )
    )
    proc = subprocess.run(
        [str(pw), "-in", str(input_path)],
        cwd=run_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=180,
    )
    out_path = run_dir / f"{name}.out"
    out_path.write_text(proc.stdout)
    ok = proc.returncode == 0 and "JOB DONE" in proc.stdout
    if expect_success and not ok:
        raise AssertionError(f"{name} failed unexpectedly; return code {proc.returncode}")
    if not expect_success and ok:
        raise AssertionError(f"{name} succeeded unexpectedly")
    return {
        "name": name,
        "returncode": proc.returncode,
        "job_done": "JOB DONE" in proc.stdout,
        "run_dir": str(run_dir),
        "stdout": str(out_path),
        "defgrad_rows": parse_defgrad(run_dir / "tmp" / f"{name}.caimd_defgrad"),
    }


def require_matrix_echo(result: dict[str, object], row1: str, row2: str, row3: str) -> None:
    stdout = Path(result["stdout"]).read_text(errors="ignore")
    if "caimd_step_defgrad matrix, per MD step" not in stdout:
        raise AssertionError(f"{result['name']}: missing initialization matrix echo")
    for row in (row1, row2, row3):
        if row not in stdout:
            raise AssertionError(f"{result['name']}: missing matrix row {row!r}")


def main() -> None:
    pw = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/private/tmp/qe_linear_defgrad_build/bin/pw.x")
    pseudo_dir = ROOT / "pseudo"
    work_root = Path("/private/tmp/qe_linear_defgrad_local_smoke")
    work_root.mkdir(parents=True, exist_ok=True)
    cases = []

    cases.append(
        run_case(
            pw,
            work_root,
            pseudo_dir,
            "parser_d13",
            "caimd=.true., caimd_loading_mode='linear_defgrad', caimd_stride=1, caimd_step_defgrad(1,3)=2.0d-6,",
            1,
        )
    )
    require_matrix_echo(
        cases[-1],
        "0.000000E+00   0.000000E+00   2.000000E-06",
        "0.000000E+00   0.000000E+00   0.000000E+00",
        "0.000000E+00   0.000000E+00   0.000000E+00",
    )

    cases.append(
        run_case(
            pw,
            work_root,
            pseudo_dir,
            "parser_d33",
            "caimd=.true., caimd_loading_mode='linear_defgrad', caimd_stride=1, caimd_step_defgrad(3,3)=-2.0d-6,",
            1,
        )
    )
    require_matrix_echo(
        cases[-1],
        "0.000000E+00   0.000000E+00   0.000000E+00",
        "0.000000E+00   0.000000E+00   0.000000E+00",
        "0.000000E+00   0.000000E+00  -2.000000E-06",
    )

    cases.append(
        run_case(
            pw,
            work_root,
            pseudo_dir,
            "zero_identity",
            "caimd=.true., caimd_loading_mode='linear_defgrad', caimd_stride=1,",
            2,
        )
    )
    if len(cases[-1]["defgrad_rows"]) != 2:
        raise AssertionError("zero_identity: expected 2 identity events")
    for row in cases[-1]["defgrad_rows"]:
        assert_identity("zero F_before", row["f_before"])
        assert_identity("zero F_after", row["f_after"])
        assert_identity("zero A_event", row["a_event"])

    cases.append(
        run_case(
            pw,
            work_root,
            pseudo_dir,
            "short_compression",
            "caimd=.true., caimd_loading_mode='linear_defgrad', caimd_stride=10, caimd_step_defgrad(3,3)=-2.0d-6,",
            12,
        )
    )
    rows = cases[-1]["defgrad_rows"]
    if [row["ionic_step"] for row in rows] != [1, 11]:
        raise AssertionError(f"short_compression: wrong event steps {rows}")
    assert_close("short compression final F33", rows[-1]["f_after"][8], 0.99996)
    expected_a33 = (1.0 - 2 * 10 * 2.0e-6) / (1.0 - 1 * 10 * 2.0e-6)
    assert_close("short compression A33", rows[-1]["a_event"][8], expected_a33)

    cases.append(
        run_case(
            pw,
            work_root,
            pseudo_dir,
            "short_shear",
            "caimd=.true., caimd_loading_mode='linear_defgrad', caimd_stride=10, caimd_step_defgrad(1,3)=2.0d-6,",
            12,
        )
    )
    rows = cases[-1]["defgrad_rows"]
    if [row["ionic_step"] for row in rows] != [1, 11]:
        raise AssertionError(f"short_shear: wrong event steps {rows}")
    assert_close("short shear final F13", rows[-1]["f_after"][2], 0.00004)
    assert_close("short shear det", rows[-1]["det_f_after"], 1.0)

    invalid = run_case(
        pw,
        work_root,
        pseudo_dir,
        "invalid_determinant",
        "caimd=.true., caimd_loading_mode='linear_defgrad', caimd_stride=1, caimd_step_defgrad(3,3)=-2.0d0,",
        1,
        expect_success=False,
    )
    invalid_stdout = Path(invalid["stdout"]).read_text(errors="ignore")
    if "F_load after event has non-positive determinant" not in invalid_stdout:
        raise AssertionError("invalid_determinant: missing determinant rejection")
    cases.append(invalid)

    summary = {
        "status": "LOCAL_LINEAR_DEFGRAD_SMOKE_PASS",
        "pw": str(pw),
        "work_root": str(work_root),
        "cases": cases,
    }
    out = ROOT / "caimd_linear_defgrad_phase1" / "evidence" / "local_smoke_summary.json"
    out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
