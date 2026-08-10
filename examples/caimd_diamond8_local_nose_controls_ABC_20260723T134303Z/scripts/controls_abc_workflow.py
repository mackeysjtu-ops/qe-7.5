#!/usr/bin/env python3
"""Local Diamond-8 QE vc-md Nose controls A/B/C workflow.

This workflow is intentionally local-only. It never uses networking, SSH,
Slurm, remote mounts, package downloads, git push, or QE source edits. All
generated files are written under the workflow directory containing this script.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parents[1]
PREV_FOUR_DOF = REPO_ROOT / "examples" / "caimd_diamond8_local_nose_four_dof_tests_0p03_20260723T115714Z"
PREV_0P30 = REPO_ROOT / "examples" / "caimd_diamond8_local_nose_continuation_0p30_20260723T102331Z"
COMMON_SRC = PREV_FOUR_DOF / "01_common_equilibration"
PW_EXE = Path("/private/tmp/qe_diamond8_af7ebfc_build/bin/pw.x")
PSEUDO = REPO_ROOT / "pseudo" / "C.UPF"
PSEUDO_NAME = "C.UPF"

EXPECTED_PW_SHA = "7f7f5f815f53983c3b41dfab1570d96e8f8816d56cdcc574a3a2e33fcf492303"
EXPECTED_PSEUDO_SHA = "9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639"

NAT = 8
MASS_C_AMU = 12.0107
ND_I = 3 * NAT - 3
DT_FS = 2.0
DT_QE = 41.341381006
NSTEP = 5000
TEMP_TARGET_K = 300.0
FNOSEP_THZ = 50.0
NHPCL = 4
WMASS = 2.0e6
ECUTWFC = 25.0
ECUTRHO = 120.0
CONV_THR = 1.0e-4
ELECTRON_MAXSTEP = 40
MIXING_BETA = 0.3
CAIMD_STRIDE = 10

RYDBERG_SI = 2.1798723611035e-18
K_BOLTZMANN_SI = 1.380649e-23
AMU_SI = 1.66053906660e-27
ELECTRONMASS_SI = 9.1093837015e-31
AMU_RY = (AMU_SI / ELECTRONMASS_SI) / 2.0
RY_TO_KELVIN = RYDBERG_SI / K_BOLTZMANN_SI

TEMP_STOP_K = 2000.0
TEMP_STOP_CONSECUTIVE = 20
CELL_COND_STOP = 1.0e8
MIN_CC_STOP_BOHR = 0.7
MAX_CC_STOP_BOHR = 20.0

COMPONENTS = ("xx", "yy", "zz", "yz", "xz", "xy")
CONTROL_ORDER = ("A", "B", "C")
CONTROLS: dict[str, dict[str, Any]] = {
    "A": {
        "dirname": "02_control_A_ordinary_vcmd_nose",
        "prefix": "diamond8_control_A_ordinary_vcmd_nose",
        "title": "ordinary PW vc-md Nose, CAIMD disabled",
        "caimd": False,
        "mask": (False, False, False, False, False, False),
        "step_strain": (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    },
    "B": {
        "dirname": "03_control_B_caimd_all_false",
        "prefix": "diamond8_control_B_caimd_all_false",
        "title": "CAIMD enabled with all six masks false and zero increments",
        "caimd": True,
        "mask": (False, False, False, False, False, False),
        "step_strain": (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    },
    "C": {
        "dirname": "04_control_C_caimd_all_true_zero_load",
        "prefix": "diamond8_control_C_caimd_all_true_zero_load",
        "title": "active all-true zero-strain CAIMD projection",
        "caimd": True,
        "mask": (True, True, True, True, True, True),
        "step_strain": (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    },
}

FLOAT_RE = re.compile(
    r"[-+]?(?:(?:\d+\.\d*|\.\d+)(?:[EeDd][-+]?\d+)?|\d+(?:[EeDd][-+]?\d+)?)"
)

CAIMD_COLUMNS = (
    "ionic_step",
    "time_ps",
    "imposed_events_done",
    "following_move_is_event",
    "nominal_xx",
    "nominal_yy",
    "nominal_zz",
    "nominal_yz",
    "nominal_xz",
    "nominal_xy",
    "actual_xx",
    "actual_yy",
    "actual_zz",
    "actual_yz",
    "actual_xz",
    "actual_xy",
    "stress_xx_GPa",
    "stress_yy_GPa",
    "stress_zz_GPa",
    "stress_yz_GPa",
    "stress_xz_GPa",
    "stress_xy_GPa",
    "temperature_K",
    "volume_bohr3",
    "velocity_jump_max",
    "velocity_jump_rms",
)


def utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run_text(cmd: list[str], cwd: Path = REPO_ROOT, timeout: int | None = None) -> str:
    return subprocess.check_output(cmd, cwd=str(cwd), text=True, stderr=subprocess.STDOUT, timeout=timeout)


def run_capture(cmd: list[str], cwd: Path = REPO_ROOT, timeout: int | None = None) -> dict[str, Any]:
    start = time.time()
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
        )
        return {
            "cmd": cmd,
            "cwd": str(cwd),
            "returncode": proc.returncode,
            "stdout": proc.stdout,
            "timeout": False,
            "elapsed_s": time.time() - start,
        }
    except subprocess.TimeoutExpired as exc:
        out = exc.stdout or ""
        if isinstance(out, bytes):
            out = out.decode(errors="replace")
        return {
            "cmd": cmd,
            "cwd": str(cwd),
            "returncode": 124,
            "stdout": out + "\nTIMEOUT\n",
            "timeout": True,
            "elapsed_s": time.time() - start,
        }


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def dfloat(x: float) -> str:
    return f"{x:.16e}".replace("e", "d")


def extract_floats(text: str) -> list[float]:
    vals: list[float] = []
    for raw in FLOAT_RE.findall(text):
        try:
            vals.append(float(raw.replace("D", "E").replace("d", "e")))
        except ValueError:
            pass
    return vals


def parse_float_file(path: Path) -> list[float]:
    if not path.exists():
        return []
    return extract_floats(path.read_text(errors="ignore"))


def transpose(a: list[list[float]]) -> list[list[float]]:
    return [list(row) for row in zip(*a)]


def matvec(a: list[list[float]], x: list[float]) -> list[float]:
    return [sum(a[i][j] * x[j] for j in range(3)) for i in range(3)]


def matmul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def inv3(a: list[list[float]]) -> list[list[float]]:
    return np.linalg.inv(np.array(a, dtype=float)).tolist()


def det3(a: list[list[float]]) -> float:
    return float(np.linalg.det(np.array(a, dtype=float)))


def condition_number(a: list[list[float]]) -> float:
    return float(np.linalg.cond(np.array(a, dtype=float)))


def singular_values(a: list[list[float]]) -> list[float]:
    return [float(v) for v in np.linalg.svd(np.array(a, dtype=float), compute_uv=False)]


def norm(x: list[float]) -> float:
    return math.sqrt(sum(v * v for v in x))


def wrap_frac(frac: list[float]) -> list[float]:
    return [x - math.floor(x) for x in frac]


def cart_to_frac(avec: list[list[float]], cart: list[float]) -> list[float]:
    return matvec(inv3(avec), cart)


def frac_to_cart(avec: list[list[float]], frac: list[float]) -> list[float]:
    return matvec(avec, frac)


def strain_components(ref: list[list[float]], cur: list[list[float]]) -> list[float]:
    f = matmul(cur, inv3(ref))
    l = [[f[i][j] - (1.0 if i == j else 0.0) for j in range(3)] for i in range(3)]
    e = [[0.5 * (l[i][j] + l[j][i]) for j in range(3)] for i in range(3)]
    return [e[0][0], e[1][1], e[2][2], e[1][2], e[0][2], e[0][1]]


def min_max_periodic_distances(avec: list[list[float]], frac: list[list[float]]) -> tuple[float | None, float | None]:
    if len(frac) < 2:
        return None, None
    best_min = float("inf")
    best_max = 0.0
    for i in range(len(frac)):
        for j in range(i + 1, len(frac)):
            pair_best = float("inf")
            for sx in (-1, 0, 1):
                for sy in (-1, 0, 1):
                    for sz in (-1, 0, 1):
                        df = [
                            frac[j][0] - frac[i][0] + sx,
                            frac[j][1] - frac[i][1] + sy,
                            frac[j][2] - frac[i][2] + sz,
                        ]
                        pair_best = min(pair_best, norm(frac_to_cart(avec, df)))
            best_min = min(best_min, pair_best)
            best_max = max(best_max, pair_best)
    return best_min, best_max


def kinetic_temperature(velocities: list[list[float]], remove_com: bool = True) -> dict[str, Any]:
    vel = [row[:] for row in velocities]
    com = [sum(row[i] for row in vel) / len(vel) for i in range(3)]
    if remove_com:
        vel = [[row[i] - com[i] for i in range(3)] for row in vel]
    mass_ry = MASS_C_AMU * AMU_RY
    ekin = 0.5 * mass_ry * sum(sum(v * v for v in row) for row in vel)
    return {
        "ekin_Ry": ekin,
        "temperature_K": 2.0 * ekin / ND_I * RY_TO_KELVIN,
        "ndof": ND_I,
        "center_of_mass_velocity_au": com,
        "center_of_mass_max_abs_au": max(abs(x) for x in com),
        "remove_com": remove_com,
    }


def rescale_velocities_to_temperature(velocities: list[list[float]], target_k: float) -> tuple[list[list[float]], dict[str, Any]]:
    com = [sum(row[i] for row in velocities) / len(velocities) for i in range(3)]
    centered = [[row[i] - com[i] for i in range(3)] for row in velocities]
    before = kinetic_temperature(centered, remove_com=False)
    scale = math.sqrt(target_k / before["temperature_K"])
    scaled = [[scale * x for x in row] for row in centered]
    after = kinetic_temperature(scaled, remove_com=True)
    return scaled, {
        "original_temperature_K": before["temperature_K"],
        "removed_center_of_mass_velocity_au": com,
        "scaling_factor": scale,
        "final_temperature_K": after["temperature_K"],
        "final_center_of_mass_velocity_au": after["center_of_mass_velocity_au"],
        "ndof": ND_I,
    }


def parse_blocks(path: Path, rows_per_block: int) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    lines = [line.rstrip() for line in path.read_text(errors="ignore").splitlines() if line.strip()]
    out: list[dict[str, Any]] = []
    i = 0
    while i < len(lines):
        head = extract_floats(lines[i])
        if len(head) < 2:
            i += 1
            continue
        rows: list[list[float]] = []
        ok = True
        for j in range(rows_per_block):
            if i + 1 + j >= len(lines):
                ok = False
                break
            vals = extract_floats(lines[i + 1 + j])
            if len(vals) < 3:
                ok = False
                break
            rows.append(vals[:3])
        if not ok:
            break
        out.append({"step": int(round(head[0])), "time_raw": head[1], "rows": rows})
        i += 1 + rows_per_block
    return out


def artifact(run_dir: Path, prefix: str, ext: str) -> Path:
    direct = run_dir / f"{prefix}.{ext}"
    if direct.exists():
        return direct
    return run_dir / "tmp" / f"{prefix}.{ext}"


def parse_md_history(path: Path) -> dict[str, Any]:
    vals = parse_float_file(path)
    first_record_floats = 5 * 3 * NAT
    second_record_floats = 4 * 9
    if len(vals) < first_record_floats + second_record_floats:
        raise RuntimeError(f"malformed PR .md history: {path} has only {len(vals)} floats")
    p = 0
    rati = [vals[p + 3 * i : p + 3 * i + 3] for i in range(NAT)]
    p += 3 * NAT
    ratd = [vals[p + 3 * i : p + 3 * i + 3] for i in range(NAT)]
    p += 3 * NAT
    rat2d = [vals[p + 3 * i : p + 3 * i + 3] for i in range(NAT)]
    p += 3 * NAT
    rat2di = [vals[p + 3 * i : p + 3 * i + 3] for i in range(NAT)]
    p += 3 * NAT
    tauold = [vals[p + 3 * i : p + 3 * i + 3] for i in range(NAT)]
    p = first_record_floats
    aveci = [vals[p + 3 * i : p + 3 * i + 3] for i in range(3)]
    p += 9
    avecd = [vals[p + 3 * i : p + 3 * i + 3] for i in range(3)]
    p += 9
    avec2d = [vals[p + 3 * i : p + 3 * i + 3] for i in range(3)]
    p += 9
    avec2di = [vals[p + 3 * i : p + 3 * i + 3] for i in range(3)]
    return {
        "rati": rati,
        "ratd": ratd,
        "rat2d": rat2d,
        "rat2di": rat2di,
        "tauold": tauold,
        "aveci": aveci,
        "avecd": avecd,
        "avec2d": avec2d,
        "avec2di": avec2di,
        "float_count": len(vals),
    }


def parse_eal(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not path.exists():
        return rows
    for line in path.read_text(errors="ignore").splitlines():
        vals = extract_floats(line)
        if len(vals) >= 7:
            uta, eka, eta, utl, ekla, etl, nst = vals[:7]
            rows.append(
                {
                    "nst": int(round(nst)),
                    "step": int(round(nst)) + 1,
                    "time_ps": (int(round(nst)) + 1) * DT_FS / 1000.0,
                    "uta_Ry": uta,
                    "eka_Ry": eka,
                    "eta_Ry": eta,
                    "utl_Ry": utl,
                    "ekla_Ry": ekla,
                    "etl_Ry": etl,
                    "ionic_temperature_K": 2.0 * eka / ND_I * RY_TO_KELVIN,
                    "extended_total_like_Ry": uta + eka + eta + utl + ekla + etl,
                }
            )
    return rows


def parse_stdout_temperatures(path: Path) -> list[float]:
    if not path.exists():
        return []
    out: list[float] = []
    for line in path.read_text(errors="ignore").splitlines():
        m = re.search(r"\bT\s*=\s*([-+0-9.EeDd]+)\s*K", line)
        if m:
            out.append(float(m.group(1).replace("D", "E").replace("d", "e")))
    return out


def parse_caimd(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not path.exists():
        return rows
    for line in path.read_text(errors="ignore").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 26:
            continue
        row: dict[str, Any] = {
            "ionic_step": int(parts[0]),
            "time_ps": float(parts[1].replace("D", "E").replace("d", "e")),
            "imposed_events_done": int(parts[2]),
            "following_move_is_event": parts[3].upper().startswith("T"),
        }
        for name, val in zip(CAIMD_COLUMNS[4:], parts[4:]):
            row[name] = float(val.replace("D", "E").replace("d", "e"))
        rows.append(row)
    return rows


def parse_stress_blocks(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    lines = path.read_text(errors="ignore").splitlines()
    rows: list[dict[str, Any]] = []
    i = 0
    while i < len(lines):
        if "total   stress" not in lines[i]:
            i += 1
            continue
        head = extract_floats(lines[i])
        pressure_kbar = head[-1] if head else None
        mat: list[list[float]] = []
        for j in range(1, 4):
            if i + j >= len(lines):
                break
            vals = extract_floats(lines[i + j])
            if len(vals) >= 6:
                mat.append(vals[3:6])
        if len(mat) == 3:
            comp_kbar = [mat[0][0], mat[1][1], mat[2][2], mat[1][2], mat[0][2], mat[0][1]]
            row = {"pressure_kbar": pressure_kbar}
            row.update({f"stress_{name}_GPa": comp_kbar[idx] / 10.0 for idx, name in enumerate(COMPONENTS)})
            rows.append(row)
        i += 4
    return rows


def cell_series(run_dir: Path, prefix: str, reference_cell: list[list[float]] | None = None) -> list[dict[str, Any]]:
    cels = parse_blocks(artifact(run_dir, prefix, "cel"), 3)
    poss = parse_blocks(artifact(run_dir, prefix, "pos"), NAT)
    rows: list[dict[str, Any]] = []
    for idx, cel in enumerate(cels):
        avec = cel["rows"]
        frac: list[list[float]] = []
        if idx < len(poss):
            try:
                frac = [wrap_frac(cart_to_frac(avec, xyz)) for xyz in poss[idx]["rows"]]
            except Exception:
                frac = []
        min_cc, max_cc = min_max_periodic_distances(avec, frac) if frac else (None, None)
        svals = singular_values(avec)
        row: dict[str, Any] = {
            "step": cel["step"],
            "time_ps": cel["step"] * DT_FS / 1000.0,
            "determinant_bohr3": det3(avec),
            "volume_bohr3": det3(avec),
            "condition_number": condition_number(avec),
            "singular_min_bohr": min(svals),
            "singular_max_bohr": max(svals),
            "minimum_C_C_distance_bohr": min_cc,
            "maximum_C_C_distance_bohr": max_cc,
        }
        if reference_cell is not None:
            eps = strain_components(reference_cell, avec)
            row.update({f"strain_{name}": eps[i] for i, name in enumerate(COMPONENTS)})
            row["strain_max_abs"] = max(abs(x) for x in eps)
        rows.append(row)
    return rows


def summary_stats(values: list[float], times: list[float] | None = None) -> dict[str, Any]:
    vals = [v for v in values if math.isfinite(v)]
    out: dict[str, Any] = {
        "count": len(values),
        "finite": len(vals) == len(values),
        "mean": sum(vals) / len(vals) if vals else None,
        "std": statistics.pstdev(vals) if len(vals) > 1 else (0.0 if len(vals) == 1 else None),
        "min": min(vals) if vals else None,
        "max": max(vals) if vals else None,
        "first_500_mean": sum(vals[: min(500, len(vals))]) / min(500, len(vals)) if vals else None,
        "final_500_mean": sum(vals[-min(500, len(vals)) :]) / min(500, len(vals)) if vals else None,
        "final": vals[-1] if vals else None,
    }
    if times is not None and len(times) == len(values) and len(vals) >= 2:
        pairs = [(t, v) for t, v in zip(times, values) if math.isfinite(t) and math.isfinite(v)]
        if len(pairs) >= 2:
            mt = sum(t for t, _ in pairs) / len(pairs)
            mv = sum(v for _, v in pairs) / len(pairs)
            denom = sum((t - mt) ** 2 for t, _ in pairs)
            out["drift_per_ps"] = sum((t - mt) * (v - mv) for t, v in pairs) / denom if denom else None
        else:
            out["drift_per_ps"] = None
    else:
        out["drift_per_ps"] = None
    return out


def rolling(values: list[float], window: int) -> list[float]:
    out: list[float] = []
    acc = 0.0
    q: list[float] = []
    for v in values:
        q.append(v)
        acc += v
        if len(q) > window:
            acc -= q.pop(0)
        out.append(acc / len(q))
    return out


def block_means(values: list[float], block: int = 500) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for start in range(0, len(values), block):
        chunk = values[start : start + block]
        if not chunk:
            continue
        rows.append(
            {
                "block": len(rows) + 1,
                "start_step": start + 1,
                "end_step": start + len(chunk),
                "count": len(chunk),
                "mean": sum(chunk) / len(chunk),
                "min": min(chunk),
                "max": max(chunk),
            }
        )
    return rows


def longest_run(values: list[float], predicate) -> int:
    best = 0
    cur = 0
    for v in values:
        if predicate(v):
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


def first_index(values: list[float], predicate) -> int | None:
    for idx, value in enumerate(values, start=1):
        if predicate(value):
            return idx
    return None


def temperature_stats(values: list[float], times: list[float]) -> dict[str, Any]:
    s = summary_stats(values, times)
    b = block_means(values, 500)
    r25 = rolling(values, 25)
    r100 = rolling(values, 100)
    s.update(
        {
            "rolling25_min": min(r25) if r25 else None,
            "rolling25_max": max(r25) if r25 else None,
            "rolling100_min": min(r100) if r100 else None,
            "rolling100_max": max(r100) if r100 else None,
            "block500_min_mean": min((row["mean"] for row in b), default=None),
            "block500_max_mean": max((row["mean"] for row in b), default=None),
            "fraction_below_150K": sum(1 for v in values if v < 150.0) / len(values) if values else None,
            "fraction_above_600K": sum(1 for v in values if v > 600.0) / len(values) if values else None,
            "longest_below_150K_steps": longest_run(values, lambda v: v < 150.0),
            "longest_above_600K_steps": longest_run(values, lambda v: v > 600.0),
            "longest_above_2000K_steps": longest_run(values, lambda v: v > 2000.0),
            "first_step_above_600K": first_index(values, lambda v: v > 600.0),
            "first_step_above_2000K": first_index(values, lambda v: v > 2000.0),
        }
    )
    return s


def classify_control(job_done: bool, completed_steps: int, temps: list[float], temp_s: dict[str, Any]) -> str:
    if not temps or not temp_s.get("finite", False):
        return "FAIL"
    if temp_s.get("longest_above_2000K_steps", 0) >= TEMP_STOP_CONSECUTIVE:
        return "FAIL"
    if temp_s.get("longest_above_600K_steps", 0) >= 25:
        return "FAIL"
    if temp_s.get("longest_below_150K_steps", 0) >= 25:
        return "FAIL"
    if not job_done or completed_steps < NSTEP:
        return "PARTIAL"
    final500 = temp_s.get("final_500_mean")
    if final500 is not None and 200.0 <= final500 <= 450.0:
        return "PASS"
    return "WARN"


def source_ref(path: Path, start: int, end: int | None = None) -> str:
    full = path if path.is_absolute() else REPO_ROOT / path
    rel = full.relative_to(REPO_ROOT)
    if end is None or end == start:
        return f"{rel}:{start}"
    return f"{rel}:{start}-{end}"


def write_implementation_confirmation() -> None:
    text = f"""# Implementation Confirmation

Generated: {utc_stamp()}

This file records a read-only source inspection of the current QE tree. No
implementation file was modified.

## Ordinary `vc-md` Nose Call Path

- `PW/src/input.f90` maps `calculation='vc-md'` plus `cell_dynamics='pr'` to
  `calc='cd'` and enables MD, forces, and variable-cell motion
  ({source_ref(Path('PW/src/input.f90'), 621, 637)}).
- `PW/src/run_pwscf.f90` computes forces and stress, writes the CAIMD pre-move
  diagnostic if enabled, then calls `move_ions`
  ({source_ref(Path('PW/src/run_pwscf.f90'), 232, 257)}).
- `PW/src/move_ions.f90` dispatches variable-cell shape MD to `vcsmd`
  ({source_ref(Path('PW/src/move_ions.f90'), 386, 393)}).
- `PW/src/vcsmd.f90` calls `ions_nosevel` before `vcmove` and `ions_noseupd`
  after `vcmove` when `tnosep` is active
  ({source_ref(Path('PW/src/vcsmd.f90'), 334, 351)}).
- `PW/src/vcsubs.f90` applies the ionic Nose friction term to `rat2d` before
  the PR/Beeman update
  ({source_ref(Path('PW/src/vcsubs.f90'), 696, 712)}).

## CAIMD Allowlist And All-False No-Op

- Active CAIMD is defined as `caimd=.true.` with at least one mask component
  true. Only active CAIMD is restricted to `ion_temperature='not_controlled'`
  or `'nose'`, and it rejects cell thermostats
  ({source_ref(Path('PW/src/input.f90'), 1739, 1756)}).
- `caimd_controls_cell()` returns `caimd_enabled_ .AND. ANY(caimd_mask_)`
  ({source_ref(Path('PW/src/caimd_module.f90'), 113, 115)}). Therefore the
  all-false Control B is expected to write CAIMD diagnostics because
  `caimd_enabled_` is true, but it must not enter the projection/tangent path.

## Projection Location And State Representation

- Active final projection is called near the end of `vcmove` after the candidate
  PR cell update and before the final metric/energy refresh
  ({source_ref(Path('PW/src/vcsubs.f90'), 1120, 1150)}).
- `ratd` is the current fractional/lattice velocity history, `rat2d` is the
  current fractional/lattice acceleration, and `rat2di` is the previous
  acceleration history. During initial velocity input, Cartesian
  `ATOMIC_VELOCITIES` are temporarily held in `rat2di`, then converted to
  `ratd` through `sigma/vcell`
  ({source_ref(Path('PW/src/vcsubs.f90'), 185, 265)}).
- `caimd_transform_atomic_history` transforms `ratd`, `rat2d`, and `rat2di`
  from the PR candidate basis into the projected basis while preserving
  physical Cartesian velocities to roundoff
  ({source_ref(Path('PW/src/caimd_module.f90'), 321, 349)}).

## Energy And Nose Synchronization

- After active projection, `vc_refresh_projected_energies` recomputes ionic
  kinetic energy `eka`, cell kinetic energy `ekla`, pressure, and ionic
  temperature from the projected `ratd` and metric `g`
  ({source_ref(Path('PW/src/vcsubs.f90'), 1159, 1241)}).
- When Nose is active, `vc_refresh_projected_nose_ekin` refreshes `ekin2nhp`
  from the same projected `ratd`/`g` before `ions_noseupd`
  ({source_ref(Path('PW/src/vcsubs.f90'), 1245, 1278)}).

## `.caimd`, `.eal`, And DOF

- `.caimd` is written by `caimd_write_premove` from `run_pwscf` before
  `move_ions`, so each row describes the state before the following ionic move.
  Its temperature is `caimd_temperature_` if available, otherwise the input
  temperature fallback at the first row
  ({source_ref(Path('PW/src/run_pwscf.f90'), 241, 257)};
  {source_ref(Path('PW/src/caimd_module.f90'), 215, 267)}).
- `caimd_note_ionic_temperature` defines `.caimd temperature_K` as
  `2*eka/(max(3*(nat-1),1)*k_B)` with `eka` from `vcsmd`
  ({source_ref(Path('PW/src/caimd_module.f90'), 201, 213)}).
- `.eal` is written in `vcsmd` as
  `uta, eka, eta, utl, ekla, etl, nst` with format `6D12.5,I6`
  ({source_ref(Path('PW/src/vcsmd.f90'), 397, 411)};
  {source_ref(Path('PW/src/vcsmd.f90'), 498, 506)}).
- For the unconstrained 8-atom system the ionic DOF convention is
  `3*nat - 3 = 21`, matching `dynamics_module::get_ndof`
  ({source_ref(Path('PW/src/dynamics_module.f90'), 2503, 2527)}).

## Velocity History Caveat

The PR `.md` file is the authoritative terminal history for `ratd`; the `.vel`
trajectory blocks are also inspected in this workflow, but they are reported as
diagnostic only unless they agree with `.eal`/terminal `.md` kinetic energy.
"""
    (ROOT / "00_metadata" / "IMPLEMENTATION_CONFIRMATION.md").write_text(text)


def extract_previous_evidence() -> dict[str, Any]:
    common_eal = parse_eal(COMMON_SRC / "diamond8_common_equil.eal")
    eal_t = [row["ionic_temperature_K"] for row in common_eal]
    terminal = extract_terminal_raw_state()
    prev_test_rows: dict[str, Any] = {}
    for key, dirname in (("test1", "02_test1_zero_load_full_pr"), ("test2", "03_test2_load_all_components_constrained")):
        run = PREV_FOUR_DOF / dirname
        caimd_rows = parse_caimd(run / "tmp" / f"diamond8_{'test1_zero_load_full_pr' if key == 'test1' else 'test2_all_components_constrained'}.caimd")
        safety = run / "SAFETY_STOP_REPORT.md"
        temps = [row["temperature_K"] for row in caimd_rows]
        nominal_actual_errors: dict[str, float | None] = {}
        for comp in COMPONENTS:
            nominal_actual_errors[comp] = max(
                (abs(row[f"nominal_{comp}"] - row[f"actual_{comp}"]) for row in caimd_rows),
                default=None,
            )
        prev_test_rows[key] = {
            "safety_report_exists": safety.exists(),
            "safety_report_excerpt": safety.read_text(errors="ignore")[:1800] if safety.exists() else None,
            "caimd_rows": len(caimd_rows),
            "temperature_min_K": min(temps) if temps else None,
            "temperature_max_K": max(temps) if temps else None,
            "temperature_final_K": temps[-1] if temps else None,
            "first_step_above_600K": first_index(temps, lambda v: v > 600.0),
            "first_step_above_2000K": first_index(temps, lambda v: v > 2000.0),
            "nominal_actual_abs_error_max_by_component": nominal_actual_errors,
            "velocity_jump_max": max((row["velocity_jump_max"] for row in caimd_rows), default=None),
            "velocity_jump_rms_max": max((row["velocity_jump_rms"] for row in caimd_rows), default=None),
        }
    old_0p30 = PREV_0P30 / "SAFETY_STOP_REPORT.md"
    return {
        "previous_four_dof_dir": str(PREV_FOUR_DOF),
        "previous_0p30_dir": str(PREV_0P30),
        "four_dof_final_report_exists": (PREV_FOUR_DOF / "FINAL_VALIDATION_REPORT.md").exists(),
        "continuation_final_report_exists": (PREV_0P30 / "FINAL_VALIDATION_REPORT.md").exists(),
        "common_raw_terminal_step": terminal["source_step"],
        "common_terminal_velocity_temperature_K_raw": terminal["temperature_from_raw_velocities"]["temperature_K"],
        "common_eal_rows": len(common_eal),
        "common_eal_final_500_mean_K": sum(eal_t[-500:]) / 500 if len(eal_t) >= 500 else None,
        "common_eal_mean_K": sum(eal_t) / len(eal_t) if eal_t else None,
        "previous_tests": prev_test_rows,
        "continuation_safety_report_exists": old_0p30.exists(),
        "continuation_safety_report_excerpt": old_0p30.read_text(errors="ignore")[:2400] if old_0p30.exists() else None,
        "verified_utc": utc_stamp(),
    }


def extract_terminal_raw_state() -> dict[str, Any]:
    cel_path = COMMON_SRC / "tmp" / "diamond8_common_equil.cel"
    pos_path = COMMON_SRC / "tmp" / "diamond8_common_equil.pos"
    md_path = COMMON_SRC / "tmp" / "diamond8_common_equil.md"
    eal_path = COMMON_SRC / "diamond8_common_equil.eal"
    cel = parse_blocks(cel_path, 3)
    pos = parse_blocks(pos_path, NAT)
    if not cel or not pos:
        raise RuntimeError("previous common equilibration lacks synchronized cel/pos blocks")
    if cel[-1]["step"] != pos[-1]["step"]:
        raise RuntimeError(f"terminal cel/pos mismatch: {cel[-1]['step']} != {pos[-1]['step']}")
    md = parse_md_history(md_path)
    avec = cel[-1]["rows"]
    frac = [wrap_frac(cart_to_frac(avec, xyz)) for xyz in pos[-1]["rows"]]
    velocities = [matvec(avec, ratd) for ratd in md["ratd"]]
    min_cc, max_cc = min_max_periodic_distances(avec, frac)
    return {
        "source_run_dir": str(COMMON_SRC),
        "source_prefix": "diamond8_common_equil",
        "source_step": cel[-1]["step"],
        "source_files": {
            "cel": str(cel_path),
            "pos": str(pos_path),
            "md": str(md_path),
            "eal": str(eal_path),
        },
        "record_indices": {
            "cel_block_index_1based": len(cel),
            "pos_block_index_1based": len(pos),
            "md_record": "last written PR history",
        },
        "atom_labels": ["C"] * NAT,
        "cell_matrix_bohr_avec_rows": avec,
        "cell_parameters_bohr_rows_for_qe": transpose(avec),
        "cartesian_positions_bohr": pos[-1]["rows"],
        "fractional_positions": frac,
        "terminal_ratd_lattice_velocity": md["ratd"],
        "cartesian_velocities_au_raw": velocities,
        "temperature_from_raw_velocities": kinetic_temperature(velocities, remove_com=True),
        "volume_bohr3": det3(avec),
        "cell_condition_number": condition_number(avec),
        "cell_singular_values_bohr": singular_values(avec),
        "minimum_periodic_C_C_distance_bohr": min_cc,
        "maximum_periodic_C_C_distance_bohr": max_cc,
        "md_float_count": md["float_count"],
        "extracted_utc": utc_stamp(),
    }


def write_state_fragment(path: Path, state: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("CELL_PARAMETERS bohr\n")
    for row in transpose(state["cell_matrix_bohr_avec_rows"]):
        lines.append("  " + "  ".join(f"{x:.14f}" for x in row) + "\n")
    lines.append("ATOMIC_POSITIONS crystal\n")
    for label, xyz in zip(state["atom_labels"], state["fractional_positions"]):
        lines.append(f"{label}  " + "  ".join(f"{x:.14f}" for x in xyz) + "\n")
    lines.append("ATOMIC_VELOCITIES { a.u }\n")
    for label, xyz in zip(state["atom_labels"], state["cartesian_velocities_au"]):
        lines.append(f"{label}  " + "  ".join(f"{x:.14e}" for x in xyz) + "\n")
    path.write_text("".join(lines))


def caimd_cell_lines(control: dict[str, Any]) -> str:
    if not control["caimd"]:
        return "  caimd = .false.,\n"
    lines = ["  caimd = .true.,", f"  caimd_stride = {CAIMD_STRIDE},"]
    for i, val in enumerate(control["step_strain"], start=1):
        lines.append(f"  caimd_step_strain({i}) = {dfloat(val)},")
    for i, val in enumerate(control["mask"], start=1):
        lines.append(f"  caimd_mask({i}) = {'.true.' if val else '.false.'},")
    return "\n".join(lines) + "\n"


def write_qe_input(control_key: str, state: dict[str, Any]) -> None:
    control = CONTROLS[control_key]
    run_dir = ROOT / control["dirname"]
    prefix = control["prefix"]
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "tmp").mkdir(exist_ok=True)
    text = f"""&CONTROL
  calculation = 'vc-md',
  restart_mode = 'from_scratch',
  prefix = '{prefix}',
  pseudo_dir = '{REPO_ROOT / "pseudo"}',
  outdir = '{run_dir / "tmp"}',
  nstep = {NSTEP},
  dt = {dfloat(DT_QE)},
  iprint = 1,
  tstress = .true.,
  tprnfor = .true.,
  disk_io = 'low',
  verbosity = 'high',
/
&SYSTEM
  ibrav = 0,
  nat = {NAT},
  ntyp = 1,
  ecutwfc = {ECUTWFC:.1f},
  ecutrho = {ECUTRHO:.1f},
  occupations = 'smearing',
  smearing = 'gaussian',
  degauss = 0.01,
  nosym = .true.,
/
&ELECTRONS
  conv_thr = {dfloat(CONV_THR)},
  electron_maxstep = {ELECTRON_MAXSTEP},
  mixing_beta = {MIXING_BETA:.3f},
/
&IONS
  ion_dynamics = 'beeman',
  ion_temperature = 'nose',
  ion_velocities = 'from_input',
  tempw = {TEMP_TARGET_K:.1f},
  fnosep(1) = {FNOSEP_THZ:.8f},
  nhpcl = {NHPCL},
  ndega = 0,
/
&CELL
  cell_dynamics = 'pr',
  cell_dofree = 'all',
  press = 0.0,
  wmass = {WMASS:.8f},
  cell_temperature = 'not_controlled',
{caimd_cell_lines(control).rstrip()}
/
ATOMIC_SPECIES
C {MASS_C_AMU:.8f} {PSEUDO_NAME}
CELL_PARAMETERS bohr
"""
    for row in transpose(state["cell_matrix_bohr_avec_rows"]):
        text += "  " + "  ".join(f"{x:.14f}" for x in row) + "\n"
    text += "ATOMIC_POSITIONS crystal\n"
    for label, xyz in zip(state["atom_labels"], state["fractional_positions"]):
        text += f"{label}  " + "  ".join(f"{x:.14f}" for x in xyz) + "\n"
    text += "K_POINTS gamma\n"
    text += "ATOMIC_VELOCITIES { a.u }\n"
    for label, xyz in zip(state["atom_labels"], state["cartesian_velocities_au"]):
        text += f"{label}  " + "  ".join(f"{x:.14e}" for x in xyz) + "\n"
    (run_dir / "qe.in").write_text(text)
    write_json(
        run_dir / "INPUT_METADATA.json",
        {
            "control": control_key,
            "title": control["title"],
            "prefix": prefix,
            "input": str(run_dir / "qe.in"),
            "caimd": control["caimd"],
            "caimd_mask": control["mask"],
            "caimd_step_strain": control["step_strain"],
            "nstep": NSTEP,
            "dt_fs": DT_FS,
            "dt_qe": DT_QE,
            "tempw_K": TEMP_TARGET_K,
            "fnosep_THz": FNOSEP_THZ,
            "nhpcl": NHPCL,
            "wmass": WMASS,
        },
    )


def section_from_input(text: str, start: str, end_markers: tuple[str, ...]) -> str:
    lines = text.splitlines()
    out: list[str] = []
    grab = False
    for line in lines:
        stripped = line.strip()
        if stripped.upper().startswith(start.upper()):
            grab = True
        if grab:
            if out and any(stripped.upper().startswith(marker.upper()) for marker in end_markers):
                break
            out.append(line)
    return "\n".join(out).strip()


def compare_inputs() -> dict[str, Any]:
    sections = ("CELL_PARAMETERS", "ATOMIC_POSITIONS", "ATOMIC_VELOCITIES")
    ref_text = (ROOT / CONTROLS["A"]["dirname"] / "qe.in").read_text()
    out: dict[str, Any] = {}

    def normalized_namelists(text: str) -> list[str]:
        normalized: list[str] = []
        for raw in text.splitlines():
            line = raw.strip()
            low = line.lower()
            if low.startswith("prefix"):
                continue
            if low.startswith("outdir"):
                continue
            if low.startswith("caimd"):
                continue
            if low.startswith("atomic_species"):
                break
            normalized.append(line)
        return normalized

    ref_normalized = normalized_namelists(ref_text)
    for key in CONTROL_ORDER:
        text = (ROOT / CONTROLS[key]["dirname"] / "qe.in").read_text()
        sec_ok = {}
        for sec in sections:
            end = {
                "CELL_PARAMETERS": ("ATOMIC_POSITIONS",),
                "ATOMIC_POSITIONS": ("K_POINTS",),
                "ATOMIC_VELOCITIES": ("",),
            }[sec]
            if sec == "ATOMIC_VELOCITIES":
                a_sec = "\n".join(ref_text.split("ATOMIC_VELOCITIES", 1)[1].splitlines()).strip()
                b_sec = "\n".join(text.split("ATOMIC_VELOCITIES", 1)[1].splitlines()).strip()
            else:
                a_sec = section_from_input(ref_text, sec, end)
                b_sec = section_from_input(text, sec, end)
            sec_ok[sec] = a_sec == b_sec
        sec_ok["non_caimd_namelist_settings"] = normalized_namelists(text) == ref_normalized
        out[key] = sec_ok
    write_json(ROOT / "01_common_state" / "INPUT_IDENTITY_CHECK.json", out)
    return out


def write_readme_and_wrapper() -> None:
    readme = f"""# Diamond8 CAIMD Nose Controls A/B/C

Generated: {utc_stamp()}

Strictly local macOS diagnostic workflow for QE 7.5 CAIMD Nose behavior.

Controls:

- A: ordinary PW `vc-md`, Nose, CAIMD disabled.
- B: `caimd=.true.` all six masks false, zero increments, intended no-op.
- C: active all-true zero-strain CAIMD projection.

Run with:

```bash
python3 scripts/controls_abc_workflow.py prepare
python3 scripts/controls_abc_workflow.py run-all
python3 scripts/controls_abc_workflow.py analyze
python3 scripts/controls_abc_workflow.py regression
```

All generated files are inside this directory. The QE source tree is not
modified by this workflow.
"""
    (ROOT / "README.md").write_text(readme)
    wrapper = """#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/controls_abc_workflow.py prepare
python3 scripts/controls_abc_workflow.py run-all
python3 scripts/controls_abc_workflow.py analyze
python3 scripts/controls_abc_workflow.py regression
"""
    path = ROOT / "scripts" / "run_controls_abc.sh"
    path.write_text(wrapper)
    path.chmod(0o755)


def prepare() -> None:
    for path in (PW_EXE, PSEUDO):
        if not path.exists():
            raise SystemExit(f"required local file missing: {path}")
    metadata = {
        "created_utc": utc_stamp(),
        "repo_root": str(REPO_ROOT),
        "workflow_root": str(ROOT),
        "platform": platform.platform(),
        "branch": run_text(["git", "branch", "--show-current"]).strip(),
        "head_short": run_text(["git", "rev-parse", "--short", "HEAD"]).strip(),
        "head_full": run_text(["git", "rev-parse", "HEAD"]).strip(),
        "recent_commits": run_text(["git", "log", "--oneline", "-5"]).splitlines(),
        "git_status_short_initial": run_text(["git", "status", "--short", "--branch"]).splitlines(),
        "pw_executable": str(PW_EXE),
        "pw_sha256": sha256(PW_EXE),
        "pw_sha256_expected": EXPECTED_PW_SHA,
        "pseudo": str(PSEUDO),
        "pseudo_sha256": sha256(PSEUDO),
        "pseudo_sha256_expected": EXPECTED_PSEUDO_SHA,
        "local_only_restrictions": {
            "remote": False,
            "ssh": False,
            "slurm": False,
            "network": False,
            "git_push": False,
            "source_modification": False,
            "commit": False,
        },
    }
    write_json(ROOT / "00_metadata" / "WORKFLOW_METADATA.json", metadata)
    write_implementation_confirmation()
    previous = extract_previous_evidence()
    write_json(ROOT / "00_metadata" / "PREVIOUS_EVIDENCE_VERIFICATION.json", previous)

    raw = extract_terminal_raw_state()
    scaled_vel, scale_info = rescale_velocities_to_temperature(raw["cartesian_velocities_au_raw"], TEMP_TARGET_K)
    common = dict(raw)
    common["cartesian_velocities_au"] = scaled_vel
    common["temperature_scaling"] = scale_info
    common["temperature_from_scaled_velocities"] = kinetic_temperature(scaled_vel, remove_com=True)
    common["units"] = {
        "cell": "bohr, QE internal avec rows with lattice vectors in columns",
        "cell_parameters_for_qe": "bohr, rows are transpose of QE internal avec",
        "positions": "crystal fractional reconstructed from synchronized .cel/.pos",
        "velocities": "a.u. Cartesian written to ATOMIC_VELOCITIES",
        "temperature": "K, 3N-3 ionic DOF",
    }
    common["atom_order_verified"] = common["atom_labels"] == ["C"] * NAT and len(common["fractional_positions"]) == NAT
    common["input_identity_requirement"] = "A/B/C share identical cell, positions, velocities, electronic, Nose, and PR settings; CAIMD keywords differ."
    write_json(ROOT / "01_common_state" / "COMMON_STATE.json", common)
    write_state_fragment(ROOT / "01_common_state" / "COMMON_STATE_INPUT_FRAGMENT.in", common)
    for key in CONTROL_ORDER:
        write_qe_input(key, common)
    compare_inputs()
    write_readme_and_wrapper()


def control_dir(key: str) -> Path:
    return ROOT / CONTROLS[key]["dirname"]


def is_complete(key: str) -> bool:
    return (control_dir(key) / "CONTROL_COMPLETE.json").exists()


def next_incomplete_control() -> str | None:
    for key in CONTROL_ORDER:
        if not is_complete(key) and not (control_dir(key) / "SAFETY_STOP_CONTEXT.json").exists():
            return key
    return None


def latest_eal_temperature(run_dir: Path, prefix: str) -> tuple[int, float] | None:
    rows = parse_eal(artifact(run_dir, prefix, "eal"))
    if not rows:
        return None
    return rows[-1]["step"], rows[-1]["ionic_temperature_K"]


def inspect_runtime_safety(run_dir: Path, prefix: str, consecutive_state: dict[str, Any]) -> dict[str, Any] | None:
    out_text = (run_dir / "qe.out").read_text(errors="ignore") if (run_dir / "qe.out").exists() else ""
    bad_markers = ("NaN", "nan", "Infinity", "infinity", "Error in routine", "SCF correction compared to forces is large")
    if any(marker in out_text for marker in bad_markers):
        return {"reason": "catastrophic stdout marker", "marker_detected": True}
    rows = parse_eal(artifact(run_dir, prefix, "eal"))
    if rows:
        temp = rows[-1]["ionic_temperature_K"]
        if not math.isfinite(temp):
            return {"reason": "nonfinite eal ionic temperature", "step": rows[-1]["step"], "temperature_K": temp}
        if temp > TEMP_STOP_K:
            consecutive_state["hot"] = consecutive_state.get("hot", 0) + 1
            consecutive_state.setdefault("first_hot_step", rows[-1]["step"])
            consecutive_state.setdefault("first_hot_temperature_K", temp)
        else:
            consecutive_state["hot"] = 0
            consecutive_state.pop("first_hot_step", None)
            consecutive_state.pop("first_hot_temperature_K", None)
        if consecutive_state.get("hot", 0) >= TEMP_STOP_CONSECUTIVE:
            return {
                "reason": "temperature above 2000 K for sustained records",
                "first_threshold_step": consecutive_state.get("first_hot_step"),
                "first_threshold_temperature_K": consecutive_state.get("first_hot_temperature_K"),
                "latest_step": rows[-1]["step"],
                "latest_temperature_K": temp,
                "consecutive_records": consecutive_state["hot"],
            }
    cells = cell_series(run_dir, prefix)
    if cells:
        last = cells[-1]
        if last["determinant_bohr3"] <= 0.0 or not math.isfinite(last["determinant_bohr3"]):
            return {"reason": "nonpositive/nonfinite cell determinant", "latest_cell": last}
        if last["condition_number"] > CELL_COND_STOP:
            return {"reason": "extreme cell condition number", "latest_cell": last}
        dmin = last.get("minimum_C_C_distance_bohr")
        dmax = last.get("maximum_C_C_distance_bohr")
        if dmin is not None and dmin < MIN_CC_STOP_BOHR:
            return {"reason": "pathological short C-C distance", "latest_cell": last}
        if dmax is not None and dmax > MAX_CC_STOP_BOHR:
            return {"reason": "pathological long C-C distance", "latest_cell": last}
    return None


def terminate_process(proc: subprocess.Popen[str]) -> None:
    proc.terminate()
    try:
        proc.wait(timeout=20)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait(timeout=10)


def run_control(key: str) -> dict[str, Any]:
    run_dir = control_dir(key)
    control = CONTROLS[key]
    prefix = control["prefix"]
    complete = run_dir / "CONTROL_COMPLETE.json"
    if complete.exists():
        return read_json(complete)
    if (run_dir / "SAFETY_STOP_CONTEXT.json").exists():
        return read_json(run_dir / "SAFETY_STOP_CONTEXT.json")
    out_path = run_dir / "qe.out"
    err_path = run_dir / "qe.err"
    start = time.time()
    with out_path.open("w") as out, err_path.open("w") as err:
        proc = subprocess.Popen(
            [str(PW_EXE), "-in", "qe.in"],
            cwd=str(run_dir),
            stdout=out,
            stderr=err,
            text=True,
        )
        monitor_state: dict[str, Any] = {}
        safety: dict[str, Any] | None = None
        last_report_step = -1
        while proc.poll() is None:
            time.sleep(2.0)
            current = latest_eal_temperature(run_dir, prefix)
            if current and current[0] != last_report_step and current[0] % 100 == 0:
                last_report_step = current[0]
                print(f"[{utc_stamp()}] control {key}: step {current[0]} eal_T={current[1]:.3f} K", flush=True)
            safety = inspect_runtime_safety(run_dir, prefix, monitor_state)
            if safety is not None:
                safety["control"] = key
                safety["prefix"] = prefix
                safety["detected_utc"] = utc_stamp()
                terminate_process(proc)
                break
        returncode = proc.returncode
    elapsed = time.time() - start
    if safety is not None:
        context = {
            "status": "SAFETY_STOP",
            "control": key,
            "returncode": returncode,
            "elapsed_s": elapsed,
            "safety": safety,
            "preserved_outputs": preserved_outputs(run_dir, prefix),
        }
        write_json(run_dir / "SAFETY_STOP_CONTEXT.json", context)
        validate_control(key)
        return context
    result = {
        "status": "PROCESS_FINISHED",
        "control": key,
        "returncode": returncode,
        "elapsed_s": elapsed,
        "finished_utc": utc_stamp(),
    }
    write_json(run_dir / "PROCESS_FINISHED.json", result)
    validation = validate_control(key)
    if validation.get("job_done") and validation.get("completed_steps", 0) >= NSTEP and validation.get("control_classification") in ("PASS", "WARN"):
        write_json(complete, {**result, "validation": validation})
    return result


def preserved_outputs(run_dir: Path, prefix: str) -> list[str]:
    paths: list[str] = []
    for ext in ("out", "err", "e", "eal", "ave", "avec", "p", "tv"):
        p = run_dir / f"{prefix}.{ext}" if ext != "out" else run_dir / "qe.out"
        if ext == "err":
            p = run_dir / "qe.err"
        if p.exists():
            paths.append(str(p))
    for ext in ("caimd", "md", "cel", "pos", "vel", "for", "xml"):
        p = run_dir / "tmp" / f"{prefix}.{ext}"
        if p.exists():
            paths.append(str(p))
    return paths


def terminal_direct_temperature(run_dir: Path, prefix: str) -> dict[str, Any] | None:
    md_path = artifact(run_dir, prefix, "md")
    cel = parse_blocks(artifact(run_dir, prefix, "cel"), 3)
    if not md_path.exists() or not cel:
        return None
    md = parse_md_history(md_path)
    avec = cel[-1]["rows"]
    velocities = [matvec(avec, ratd) for ratd in md["ratd"]]
    return {
        "step": cel[-1]["step"],
        "source": ".md terminal ratd transformed by terminal .cel",
        "temperature": kinetic_temperature(velocities, remove_com=True),
        "velocity_norms_au": [norm(v) for v in velocities],
        "center_of_mass_removed": True,
    }


def vel_file_diagnostic(run_dir: Path, prefix: str) -> dict[str, Any]:
    blocks = parse_blocks(artifact(run_dir, prefix, "vel"), NAT)
    if not blocks:
        return {"exists": False, "usable_as_physical_direct_temperature": False}
    temps = [kinetic_temperature(block["rows"], remove_com=True)["temperature_K"] for block in blocks]
    first_rows = blocks[0]["rows"]
    constant = all(max(abs(blocks[i]["rows"][a][b] - first_rows[a][b]) for a in range(NAT) for b in range(3)) < 1.0e-14 for i in range(len(blocks)))
    return {
        "exists": True,
        "blocks": len(blocks),
        "temperature_stats_if_interpreted_as_cartesian_au": temperature_stats(temps, [b["step"] * DT_FS / 1000.0 for b in blocks]),
        "constant_over_time": constant,
        "usable_as_physical_direct_temperature": False if constant else None,
        "caveat": "QE PR .vel blocks were constant in previous workflows; terminal .md ratd is used for authoritative direct cross-check.",
    }


def validate_control(key: str) -> dict[str, Any]:
    control = CONTROLS[key]
    run_dir = control_dir(key)
    prefix = control["prefix"]
    out_text = (run_dir / "qe.out").read_text(errors="ignore") if (run_dir / "qe.out").exists() else ""
    eal = parse_eal(artifact(run_dir, prefix, "eal"))
    eal_t = [row["ionic_temperature_K"] for row in eal]
    eal_times = [row["time_ps"] for row in eal]
    stdout_t = parse_stdout_temperatures(run_dir / "qe.out")
    cells = cell_series(run_dir, prefix, read_json(ROOT / "01_common_state" / "COMMON_STATE.json")["cell_matrix_bohr_avec_rows"])
    caimd_rows = parse_caimd(artifact(run_dir, prefix, "caimd"))
    terminal = terminal_direct_temperature(run_dir, prefix)
    vel_diag = vel_file_diagnostic(run_dir, prefix)
    completed_steps = len(parse_blocks(artifact(run_dir, prefix, "cel"), 3))
    job_done = "JOB DONE" in out_text
    eal_s = temperature_stats(eal_t, eal_times) if eal_t else {}
    cclass = classify_control(job_done, completed_steps, eal_t, eal_s) if eal_t else "FAIL"
    dets = [row["determinant_bohr3"] for row in cells]
    conds = [row["condition_number"] for row in cells]
    strains = [row["strain_max_abs"] for row in cells if "strain_max_abs" in row]
    caimd_temp = [row["temperature_K"] for row in caimd_rows]
    validation: dict[str, Any] = {
        "control": key,
        "title": control["title"],
        "prefix": prefix,
        "validated_utc": utc_stamp(),
        "job_done": job_done,
        "returncode": read_json(run_dir / "PROCESS_FINISHED.json").get("returncode") if (run_dir / "PROCESS_FINISHED.json").exists() else None,
        "completed_steps": completed_steps,
        "simulated_time_ps": completed_steps * DT_FS / 1000.0,
        "control_classification": cclass,
        "primary_temperature_source": ".eal ionic kinetic energy",
        "eal_temperature_K": eal_s,
        "stdout_temperature_K": temperature_stats(stdout_t, [(i + 1) * DT_FS / 1000.0 for i in range(len(stdout_t))]) if stdout_t else None,
        "caimd_rows": len(caimd_rows),
        "caimd_temperature_K": temperature_stats(caimd_temp, [row["time_ps"] for row in caimd_rows]) if caimd_temp else None,
        "terminal_direct_temperature": terminal,
        "vel_file_diagnostic": vel_diag,
        "cell_determinant_bohr3": summary_stats(dets, [row["time_ps"] for row in cells]) if dets else {},
        "cell_condition_number": summary_stats(conds, [row["time_ps"] for row in cells]) if conds else {},
        "cell_strain_max_abs": summary_stats(strains, [row["time_ps"] for row in cells]) if strains else {},
        "volume_range_bohr3": [min(dets), max(dets)] if dets else None,
        "minimum_C_C_distance_bohr": summary_stats([row["minimum_C_C_distance_bohr"] for row in cells if row["minimum_C_C_distance_bohr"] is not None]) if cells else {},
        "maximum_C_C_distance_bohr": summary_stats([row["maximum_C_C_distance_bohr"] for row in cells if row["maximum_C_C_distance_bohr"] is not None]) if cells else {},
        "qe_error_lines": [line.strip() for line in out_text.splitlines() if "Error in routine" in line or "STOP" in line],
        "preserved_outputs": preserved_outputs(run_dir, prefix),
    }
    if caimd_rows:
        validation["caimd_event_count_last"] = caimd_rows[-1]["imposed_events_done"]
        validation["caimd_event_count_expected_for_completed_steps"] = completed_steps // CAIMD_STRIDE
        validation["caimd_velocity_jump_max"] = max((row["velocity_jump_max"] for row in caimd_rows), default=None)
        validation["caimd_velocity_jump_rms_max"] = max((row["velocity_jump_rms"] for row in caimd_rows), default=None)
        for comp in COMPONENTS:
            actual = [row[f"actual_{comp}"] for row in caimd_rows]
            nominal = [row[f"nominal_{comp}"] for row in caimd_rows]
            validation[f"caimd_{comp}_actual_abs_max"] = max((abs(v) for v in actual), default=None)
            validation[f"caimd_{comp}_nominal_actual_abs_error_max"] = max((abs(a - n) for a, n in zip(actual, nominal)), default=None)
    write_json(run_dir / "VALIDATION.json", validation)
    if cclass in ("FAIL", "PARTIAL") or (run_dir / "SAFETY_STOP_CONTEXT.json").exists():
        write_safety_report(run_dir, validation)
    return validation


def write_safety_report(run_dir: Path, validation: dict[str, Any]) -> None:
    context = read_json(run_dir / "SAFETY_STOP_CONTEXT.json") if (run_dir / "SAFETY_STOP_CONTEXT.json").exists() else None
    text = f"""# Safety Stop Report

Generated: {utc_stamp()}

Status: {validation.get('control_classification')}

Control: {validation.get('control')}

Run directory: `{run_dir}`

Prefix: `{validation.get('prefix')}`

Completed trajectory blocks: {validation.get('completed_steps')}

JOB DONE: {validation.get('job_done')}

Primary temperature source: {validation.get('primary_temperature_source')}

EAL temperature stats:

```json
{json.dumps(validation.get('eal_temperature_K'), indent=2)}
```

Safety context:

```json
{json.dumps(context, indent=2)}
```

QE error lines:

```json
{json.dumps(validation.get('qe_error_lines'), indent=2)}
```

The partial outputs are preserved in this directory and are not treated as a
successful 5000-step control.
"""
    (run_dir / "SAFETY_STOP_REPORT.md").write_text(text)


def run_all() -> None:
    for key in CONTROL_ORDER:
        if key != "A" and (control_dir("A") / "SAFETY_STOP_CONTEXT.json").exists():
            ctx = read_json(control_dir("A") / "SAFETY_STOP_CONTEXT.json")
            reason = ctx.get("safety", {}).get("reason", "")
            if "temperature above 2000" in reason or "nonfinite" in reason or "determinant" in reason:
                print(f"Stopping before control {key}: Control A had an unambiguous safety stop.", flush=True)
                break
        if is_complete(key):
            print(f"control {key} already complete", flush=True)
            continue
        if (control_dir(key) / "SAFETY_STOP_CONTEXT.json").exists():
            print(f"control {key} already safety-stopped", flush=True)
            continue
        print(f"[{utc_stamp()}] launching control {key}: {CONTROLS[key]['title']}", flush=True)
        run_control(key)


def compare_numeric_files(a: Path, b: Path, tolerance_abs: float = 1.0e-10, tolerance_rel: float = 1.0e-10) -> dict[str, Any]:
    av = parse_float_file(a)
    bv = parse_float_file(b)
    n = min(len(av), len(bv))
    if n == 0:
        return {
            "path_a": str(a),
            "path_b": str(b),
            "count_a": len(av),
            "count_b": len(bv),
            "max_abs": None,
            "max_rel": None,
            "rms": None,
            "first_tolerance_exceeded_index": None,
        }
    diffs = [av[i] - bv[i] for i in range(n)]
    max_abs = max(abs(x) for x in diffs)
    rels = [abs(diffs[i]) / max(abs(av[i]), abs(bv[i]), 1.0) for i in range(n)]
    rms = math.sqrt(sum(x * x for x in diffs) / n)
    first = None
    for i, (da, rel) in enumerate(zip(diffs, rels), start=1):
        if abs(da) > tolerance_abs and rel > tolerance_rel:
            first = i
            break
    return {
        "path_a": str(a),
        "path_b": str(b),
        "count_a": len(av),
        "count_b": len(bv),
        "compared_count": n,
        "same_count": len(av) == len(bv),
        "max_abs": max_abs,
        "max_rel": max(rels),
        "rms": rms,
        "first_tolerance_exceeded_index": first,
        "tolerance_abs": tolerance_abs,
        "tolerance_rel": tolerance_rel,
    }


def compare_a_b() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    a_dir = control_dir("A")
    b_dir = control_dir("B")
    a_prefix = CONTROLS["A"]["prefix"]
    b_prefix = CONTROLS["B"]["prefix"]
    files = {
        "e": ("e", False),
        "eal": ("eal", False),
        "ave": ("ave", False),
        "avec": ("avec", False),
        "p": ("p", False),
        "tv": ("tv", False),
        "cel": ("cel", True),
        "pos": ("pos", True),
        "vel": ("vel", True),
        "for": ("for", True),
        "md": ("md", True),
    }
    for name, (ext, tmp) in files.items():
        ap = (a_dir / "tmp" / f"{a_prefix}.{ext}") if tmp else (a_dir / f"{a_prefix}.{ext}")
        bp = (b_dir / "tmp" / f"{b_prefix}.{ext}") if tmp else (b_dir / f"{b_prefix}.{ext}")
        if ap.exists() and bp.exists():
            rows.append({"observable": name, **compare_numeric_files(ap, bp)})
    write_csv(ROOT / "analysis" / "control_A_vs_B_differences.csv", rows)
    return rows


def classify_workflow(validations: dict[str, dict[str, Any]], ab_rows: list[dict[str, Any]]) -> dict[str, Any]:
    a = validations.get("A")
    b = validations.get("B")
    c = validations.get("C")
    ab_ok = bool(ab_rows) and all((row.get("max_abs") is not None and row["max_abs"] <= 1.0e-8) for row in ab_rows)
    status = "PARTIAL"
    conclusion = "Insufficient completed controls for a full decision."
    earliest = None
    if a and a.get("control_classification") == "FAIL":
        status = "BASELINE_NOSE_FAIL"
        conclusion = "Control A already exhibits physical instability in ordinary PW vc-md + Nose; this does not isolate a CAIMD-specific bug."
        earliest = {"control": "A", "observable": "eal ionic temperature", "step": a.get("eal_temperature_K", {}).get("first_step_above_600K")}
    elif a and a.get("control_classification") in ("PASS", "WARN") and b and b.get("control_classification") == "FAIL":
        status = "CAIMD_NOOP_FAIL"
        conclusion = "Control A is bounded but Control B fails; the all-false CAIMD no-op path is suspect."
        earliest = {"control": "B", "observable": "eal ionic temperature", "step": b.get("eal_temperature_K", {}).get("first_step_above_600K")}
    elif a and b and a.get("control_classification") in ("PASS", "WARN") and b.get("control_classification") in ("PASS", "WARN") and not ab_ok:
        status = "CAIMD_NOOP_FAIL"
        conclusion = "Control B did not reproduce Control A within strict numeric tolerance."
        for row in ab_rows:
            if row.get("first_tolerance_exceeded_index"):
                earliest = {"control": "A_vs_B", "observable": row["observable"], "float_index": row["first_tolerance_exceeded_index"]}
                break
    elif a and b and c and a.get("control_classification") in ("PASS", "WARN") and b.get("control_classification") in ("PASS", "WARN") and ab_ok and c.get("control_classification") == "FAIL":
        status = "ACTIVE_PROJECTION_FAIL"
        conclusion = "Controls A/B are bounded and no-op matches, but Control C fails; failure is isolated to active zero-strain projection or Nose state synchronization."
        earliest = {"control": "C", "observable": "eal ionic temperature", "step": c.get("eal_temperature_K", {}).get("first_step_above_600K")}
    elif a and b and c and all(x.get("control_classification") in ("PASS", "WARN") for x in (a, b, c)) and ab_ok:
        c_strain = c.get("cell_strain_max_abs", {}).get("max")
        c_ok = c_strain is not None and c_strain <= 1.0e-8
        if c_ok:
            status = "SUCCESS_ALL_STABLE"
            conclusion = "A, B, and C completed with bounded temperatures; B matches A and C maintains zero strain."
        else:
            status = "ACTIVE_PROJECTION_FAIL"
            conclusion = "Control C completed but did not hold all zero-strain components to numerical precision."
            earliest = {"control": "C", "observable": "cell strain residual", "max": c_strain}
    diag = {
        "status": status,
        "conclusion": conclusion,
        "earliest_diverging_observable": earliest,
        "a_vs_b_noop_ok": ab_ok,
        "generated_utc": utc_stamp(),
    }
    write_json(ROOT / "analysis" / "DIAGNOSTIC_CLASSIFICATION.json", diag)
    return diag


def plot_series(path: Path, title: str, series: list[tuple[str, list[float], list[float]]], ylabel: str = "") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width, height = 1000, 650
    margin_l, margin_r, margin_t, margin_b = 90, 30, 70, 80
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf"]
    font = ImageFont.load_default()
    draw.text((margin_l, 20), title, fill="black", font=font)
    valid_points = [(x, y) for _, xs, ys in series for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if not valid_points:
        draw.text((margin_l, margin_t), "No data available", fill="black", font=font)
        img.save(path)
        return
    xs_all = [p[0] for p in valid_points]
    ys_all = [p[1] for p in valid_points]
    xmin, xmax = min(xs_all), max(xs_all)
    ymin, ymax = min(ys_all), max(ys_all)
    if xmin == xmax:
        xmax = xmin + 1.0
    if ymin == ymax:
        pad = max(abs(ymin) * 0.1, 1.0)
        ymin -= pad
        ymax += pad
    else:
        pad = 0.08 * (ymax - ymin)
        ymin -= pad
        ymax += pad
    x0, y0 = margin_l, height - margin_b
    x1, y1 = width - margin_r, margin_t
    draw.rectangle((x0, y1, x1, y0), outline="black")
    for i in range(6):
        x = x0 + i * (x1 - x0) / 5
        y = y0 - i * (y0 - y1) / 5
        draw.line((x, y0, x, y0 + 5), fill="black")
        draw.line((x0 - 5, y, x0, y), fill="black")
        draw.text((x - 20, y0 + 10), f"{xmin + i*(xmax-xmin)/5:.2g}", fill="black", font=font)
        draw.text((8, y - 6), f"{ymin + i*(ymax-ymin)/5:.3g}", fill="black", font=font)
    draw.text((width // 2 - 30, height - 35), "time / ps", fill="black", font=font)
    draw.text((10, 20), ylabel, fill="black", font=font)
    for idx, (label, xs, ys) in enumerate(series):
        color = colors[idx % len(colors)]
        pts = []
        for x, y in zip(xs, ys):
            if math.isfinite(x) and math.isfinite(y):
                px = x0 + (x - xmin) / (xmax - xmin) * (x1 - x0)
                py = y0 - (y - ymin) / (ymax - ymin) * (y0 - y1)
                pts.append((px, py))
        if len(pts) >= 2:
            draw.line(pts, fill=color, width=2)
        elif len(pts) == 1:
            px, py = pts[0]
            draw.ellipse((px - 3, py - 3, px + 3, py + 3), fill=color)
        draw.rectangle((margin_l + idx * 190, height - 58, margin_l + idx * 190 + 12, height - 46), fill=color)
        draw.text((margin_l + idx * 190 + 18, height - 60), label, fill="black", font=font)
    img.save(path)


def write_analysis_outputs(validations: dict[str, dict[str, Any]], ab_rows: list[dict[str, Any]], diag: dict[str, Any]) -> None:
    temp_rows: list[dict[str, Any]] = []
    cross_rows: list[dict[str, Any]] = []
    cell_rows: list[dict[str, Any]] = []
    nose_rows: list[dict[str, Any]] = []
    c_projection_rows: list[dict[str, Any]] = []
    eal_series: list[tuple[str, list[float], list[float]]] = []
    direct_series: list[tuple[str, list[float], list[float]]] = []
    rolling_series: list[tuple[str, list[float], list[float]]] = []
    volume_series: list[tuple[str, list[float], list[float]]] = []
    det_series: list[tuple[str, list[float], list[float]]] = []
    cond_series: list[tuple[str, list[float], list[float]]] = []
    kinetic_series: list[tuple[str, list[float], list[float]]] = []
    ext_series: list[tuple[str, list[float], list[float]]] = []
    velnorm_series: list[tuple[str, list[float], list[float]]] = []
    strain_series: list[tuple[str, list[float], list[float]]] = []
    common = read_json(ROOT / "01_common_state" / "COMMON_STATE.json")
    initial_temp = common["temperature_from_scaled_velocities"]["temperature_K"]
    for key in CONTROL_ORDER:
        control = CONTROLS[key]
        run_dir = control_dir(key)
        prefix = control["prefix"]
        val = validations.get(key, {})
        eal = parse_eal(artifact(run_dir, prefix, "eal"))
        eal_t = [row["ionic_temperature_K"] for row in eal]
        times = [row["time_ps"] for row in eal]
        cells = cell_series(run_dir, prefix, common["cell_matrix_bohr_avec_rows"])
        caimd_rows = parse_caimd(artifact(run_dir, prefix, "caimd"))
        tstats = val.get("eal_temperature_K", {})
        temp_rows.append(
            {
                "control": key,
                "completion_status": val.get("control_classification"),
                "completed_steps": val.get("completed_steps"),
                "simulated_time_ps": val.get("simulated_time_ps"),
                "initial_direct_temperature_K": initial_temp,
                "mean_direct_temperature_K": val.get("terminal_direct_temperature", {}).get("temperature", {}).get("temperature_K"),
                "mean_eal_temperature_K": tstats.get("mean"),
                "std_eal_temperature_K": tstats.get("std"),
                "first_500_mean_K": tstats.get("first_500_mean"),
                "final_500_mean_K": tstats.get("final_500_mean"),
                "min_K": tstats.get("min"),
                "max_K": tstats.get("max"),
                "drift_K_per_ps": tstats.get("drift_per_ps"),
                "fraction_below_150K": tstats.get("fraction_below_150K"),
                "fraction_above_600K": tstats.get("fraction_above_600K"),
                "first_step_above_600K": tstats.get("first_step_above_600K"),
                "first_step_above_2000K": tstats.get("first_step_above_2000K"),
                "determinant_min": val.get("cell_determinant_bohr3", {}).get("min"),
                "condition_number_max": val.get("cell_condition_number", {}).get("max"),
                "volume_min": val.get("volume_range_bohr3", [None, None])[0] if val.get("volume_range_bohr3") else None,
                "volume_max": val.get("volume_range_bohr3", [None, None])[1] if val.get("volume_range_bohr3") else None,
                "maximum_cell_strain": val.get("cell_strain_max_abs", {}).get("max"),
                "classification": val.get("control_classification"),
            }
        )
        terminal_t = val.get("terminal_direct_temperature", {})
        cross_rows.append(
            {
                "control": key,
                "source": "input_direct_scaled_velocities",
                "step": 0,
                "time_ps": 0.0,
                "temperature_K": initial_temp,
                "notes": "COMMON_STATE velocities after COM removal and uniform rescale",
            }
        )
        if terminal_t:
            cross_rows.append(
                {
                    "control": key,
                    "source": "terminal_md_direct",
                    "step": terminal_t.get("step"),
                    "time_ps": (terminal_t.get("step") or 0) * DT_FS / 1000.0,
                    "temperature_K": terminal_t.get("temperature", {}).get("temperature_K"),
                    "notes": terminal_t.get("source"),
                }
            )
        if eal:
            cross_rows.append(
                {
                    "control": key,
                    "source": "eal_terminal_ionic_ekin",
                    "step": eal[-1]["step"],
                    "time_ps": eal[-1]["time_ps"],
                    "temperature_K": eal[-1]["ionic_temperature_K"],
                    "notes": ".eal column 2 eka converted with 3N-3 DOF",
                }
            )
        if caimd_rows:
            best_same = None
            best_shift = None
            ctemps = [r["temperature_K"] for r in caimd_rows]
            for shift in (0, 1):
                pairs = []
                for idx, ct in enumerate(ctemps):
                    j = idx - shift
                    if 0 <= j < len(eal_t):
                        pairs.append(abs(ct - eal_t[j]))
                if pairs:
                    mean_abs = sum(pairs) / len(pairs)
                    if best_same is None or mean_abs < best_same:
                        best_same = mean_abs
                        best_shift = shift
            cross_rows.append(
                {
                    "control": key,
                    "source": "caimd_alignment",
                    "step": None,
                    "time_ps": None,
                    "temperature_K": None,
                    "notes": f"best caimd-vs-eal offset={best_shift}, mean_abs_diff_K={best_same}",
                }
            )
        for row in cells:
            cell_rows.append({"control": key, **row})
        for row in eal:
            nose_rows.append(
                {
                    "control": key,
                    "step": row["step"],
                    "time_ps": row["time_ps"],
                    "uta_Ry": row["uta_Ry"],
                    "eka_Ry": row["eka_Ry"],
                    "eta_Ry": row["eta_Ry"],
                    "utl_Ry": row["utl_Ry"],
                    "ekla_Ry": row["ekla_Ry"],
                    "etl_Ry": row["etl_Ry"],
                    "extended_total_like_Ry": row["extended_total_like_Ry"],
                    "ionic_temperature_K": row["ionic_temperature_K"],
                }
            )
        if key == "C":
            for row in caimd_rows:
                c_projection_rows.append(
                    {
                        "control": key,
                        "ionic_step": row["ionic_step"],
                        "time_ps": row["time_ps"],
                        "event_count": row["imposed_events_done"],
                        "following_move_is_event": row["following_move_is_event"],
                        "max_abs_actual_strain": max(abs(row[f"actual_{c}"]) for c in COMPONENTS),
                        "max_abs_nominal_strain": max(abs(row[f"nominal_{c}"]) for c in COMPONENTS),
                        "velocity_jump_max": row["velocity_jump_max"],
                        "velocity_jump_rms": row["velocity_jump_rms"],
                    }
                )
        eal_series.append((key, times, eal_t))
        rolling_series.append((key, times, rolling(eal_t, 100)))
        kinetic_series.append((key, times, [row["eka_Ry"] for row in eal]))
        ext0 = eal[0]["extended_total_like_Ry"] if eal else 0.0
        ext_series.append((key, times, [row["extended_total_like_Ry"] - ext0 for row in eal]))
        volume_series.append((key, [r["time_ps"] for r in cells], [r["volume_bohr3"] for r in cells]))
        det_series.append((key, [r["time_ps"] for r in cells], [r["determinant_bohr3"] for r in cells]))
        cond_series.append((key, [r["time_ps"] for r in cells], [r["condition_number"] for r in cells]))
        strain_series.append((key, [r["time_ps"] for r in cells], [r.get("strain_max_abs", 0.0) for r in cells]))
        if terminal_t:
            direct_series.append(
                (
                    key,
                    [0.0, (terminal_t.get("step") or 0) * DT_FS / 1000.0],
                    [initial_temp, terminal_t.get("temperature", {}).get("temperature_K")],
                )
            )
        vel_diag = val.get("vel_file_diagnostic", {})
        if vel_diag.get("exists"):
            vel_blocks = parse_blocks(artifact(run_dir, prefix, "vel"), NAT)
            velnorm_series.append(
                (
                    key,
                    [b["step"] * DT_FS / 1000.0 for b in vel_blocks],
                    [statistics.mean(norm(v) for v in b["rows"]) for b in vel_blocks],
                )
            )
    write_csv(ROOT / "analysis" / "temperature_statistics.csv", temp_rows)
    write_csv(ROOT / "analysis" / "temperature_crosscheck.csv", cross_rows)
    write_csv(ROOT / "analysis" / "cell_statistics.csv", cell_rows)
    nose_stat_rows: list[dict[str, Any]] = []
    for key in CONTROL_ORDER:
        key_rows = [row for row in nose_rows if row["control"] == key]
        for field in ("uta_Ry", "eka_Ry", "eta_Ry", "utl_Ry", "ekla_Ry", "etl_Ry", "extended_total_like_Ry", "ionic_temperature_K"):
            vals = [row[field] for row in key_rows]
            times = [row["time_ps"] for row in key_rows]
            s = summary_stats(vals, times)
            nose_stat_rows.append(
                {
                    "control": key,
                    "component": field,
                    "count": s.get("count"),
                    "mean": s.get("mean"),
                    "std": s.get("std"),
                    "min": s.get("min"),
                    "max": s.get("max"),
                    "initial": vals[0] if vals else None,
                    "final": vals[-1] if vals else None,
                    "final_minus_initial": vals[-1] - vals[0] if len(vals) >= 2 else None,
                    "drift_per_ps": s.get("drift_per_ps"),
                }
            )
    write_csv(ROOT / "analysis" / "nose_extended_energy_statistics.csv", nose_stat_rows)
    write_csv(ROOT / "analysis" / "nose_extended_energy_timeseries.csv", nose_rows)
    write_csv(ROOT / "analysis" / "control_C_projection_statistics.csv", c_projection_rows)
    plot_series(ROOT / "analysis" / "temperature_eal_ABC_vs_time.png", "EAL Ionic Temperature A/B/C", eal_series, "K")
    plot_series(ROOT / "analysis" / "temperature_direct_ABC_vs_time.png", "Direct Velocity Temperature Cross-Checks A/B/C", direct_series, "K")
    plot_series(ROOT / "analysis" / "temperature_rolling100_ABC.png", "100-Step Rolling EAL Temperature", rolling_series, "K")
    block_series = []
    for key in CONTROL_ORDER:
        eal = parse_eal(artifact(control_dir(key), CONTROLS[key]["prefix"], "eal"))
        vals = [r["ionic_temperature_K"] for r in eal]
        b = block_means(vals, 500)
        block_series.append((key, [row["end_step"] * DT_FS / 1000.0 for row in b], [row["mean"] for row in b]))
    plot_series(ROOT / "analysis" / "temperature_block500_ABC.png", "500-Step Block Mean Temperature", block_series, "K")
    a_eal = parse_eal(artifact(control_dir("A"), CONTROLS["A"]["prefix"], "eal"))
    b_eal = parse_eal(artifact(control_dir("B"), CONTROLS["B"]["prefix"], "eal"))
    n = min(len(a_eal), len(b_eal))
    plot_series(
        ROOT / "analysis" / "temperature_difference_A_minus_B.png",
        "Temperature Difference A-B",
        [("A-B", [a_eal[i]["time_ps"] for i in range(n)], [a_eal[i]["ionic_temperature_K"] - b_eal[i]["ionic_temperature_K"] for i in range(n)])],
        "K",
    )
    plot_series(ROOT / "analysis" / "cell_strains_ABC.png", "Max Absolute Cell Strain", strain_series, "strain")
    plot_series(ROOT / "analysis" / "cell_determinant_ABC.png", "Cell Determinant", det_series, "bohr^3")
    plot_series(ROOT / "analysis" / "cell_condition_number_ABC.png", "Cell Condition Number", cond_series, "cond")
    plot_series(ROOT / "analysis" / "volume_ABC.png", "Volume", volume_series, "bohr^3")
    plot_series(ROOT / "analysis" / "ionic_kinetic_energy_ABC.png", "Ionic Kinetic Energy EKA", kinetic_series, "Ry")
    plot_series(ROOT / "analysis" / "nose_extended_terms_ABC.png", "Extended Total-Like EAL Sum", ext_series, "Ry delta")
    plot_series(ROOT / "analysis" / "extended_energy_drift_ABC.png", "Extended Energy Drift", ext_series, "Ry delta")
    plot_series(ROOT / "analysis" / "velocity_norm_ABC.png", "QE .vel Mean Norm Diagnostic", velnorm_series, "a.u.")
    c_proj_series = (
        [("C velocity_jump_max", [r["time_ps"] for r in c_projection_rows], [r["velocity_jump_max"] for r in c_projection_rows])]
        if c_projection_rows
        else []
    )
    plot_series(ROOT / "analysis" / "velocity_jump_C.png", "Control C CAIMD Velocity Jump", c_proj_series, "a.u.")
    ab_plot_series = [("max_abs", list(range(1, len(ab_rows) + 1)), [row["max_abs"] or 0.0 for row in ab_rows])]
    plot_series(ROOT / "analysis" / "A_B_noop_differences.png", "A vs B Numeric File Differences", ab_plot_series, "max abs")
    c_rows = parse_caimd(artifact(control_dir("C"), CONTROLS["C"]["prefix"], "caimd"))
    if c_rows:
        c_eal = parse_eal(artifact(control_dir("C"), CONTROLS["C"]["prefix"], "eal"))
        m = min(len(c_rows), len(c_eal))
        plot_series(
            ROOT / "analysis" / "temperature_crosscheck_C.png",
            "Control C CAIMD vs EAL Temperature",
            [
                ("C caimd", [c_rows[i]["time_ps"] for i in range(m)], [c_rows[i]["temperature_K"] for i in range(m)]),
                ("C eal", [c_eal[i]["time_ps"] for i in range(m)], [c_eal[i]["ionic_temperature_K"] for i in range(m)]),
            ],
            "K",
        )
    else:
        plot_series(ROOT / "analysis" / "temperature_crosscheck_C.png", "Control C CAIMD vs EAL Temperature", [], "K")
    write_final_report(validations, diag, ab_rows)


def write_final_report(validations: dict[str, dict[str, Any]], diag: dict[str, Any], ab_rows: list[dict[str, Any]]) -> None:
    meta = read_json(ROOT / "00_metadata" / "WORKFLOW_METADATA.json")
    common = read_json(ROOT / "01_common_state" / "COMMON_STATE.json")
    prev = read_json(ROOT / "00_metadata" / "PREVIOUS_EVIDENCE_VERIFICATION.json")
    regression_path = ROOT / "evidence" / "REGRESSION_RESULTS.json"
    regression_result = read_json(regression_path) if regression_path.exists() else None
    max_ab = max((row.get("max_abs") or 0.0 for row in ab_rows), default=None)
    max_ab_row = max(ab_rows, key=lambda row: row.get("max_abs") or -1.0) if ab_rows else None
    c_val = validations.get("C", {})
    c_strain = c_val.get("cell_strain_max_abs", {}).get("max")
    c_jump = c_val.get("caimd_velocity_jump_max")
    lines = [
        "# Final Validation Report",
        "",
        f"Generated: {utc_stamp()}",
        "",
        f"Status: `{diag['status']}`",
        "",
        "## Outcome",
        "",
        diag["conclusion"],
        "",
        "## Repository And Local Inputs",
        "",
        f"- Repository root: `{meta['repo_root']}`",
        f"- Branch/HEAD: `{meta['branch']}` / `{meta['head_short']}`",
        f"- Executable: `{meta['pw_executable']}`",
        f"- Executable SHA256: `{meta['pw_sha256']}`",
        f"- Pseudopotential: `{meta['pseudo']}`",
        f"- Pseudopotential SHA256: `{meta['pseudo_sha256']}`",
        "- QE source modified: no",
        "- Git commit created: no",
        "- Leonardo/SSH/Slurm/network/git push: not used",
        "",
        "## Prior Evidence Rechecked",
        "",
        f"- Four-DOF final report exists: {prev['four_dof_final_report_exists']}",
        f"- 0.30 continuation final report exists: {prev['continuation_final_report_exists']}",
        f"- Common terminal velocity temperature from raw `.md`: {prev['common_terminal_velocity_temperature_K_raw']:.12f} K",
        f"- Common `.eal` final-500 mean: {prev['common_eal_final_500_mean_K']:.12f} K",
        f"- Previous Test1 max CAIMD temperature: {prev['previous_tests']['test1']['temperature_max_K']} K",
        f"- Previous Test2 max CAIMD temperature: {prev['previous_tests']['test2']['temperature_max_K']} K",
        "",
        "## Common State",
        "",
        f"- Source step: {common['source_step']}",
        f"- Original velocity temperature: {common['temperature_scaling']['original_temperature_K']:.12f} K",
        f"- Removed COM velocity au: {common['temperature_scaling']['removed_center_of_mass_velocity_au']}",
        f"- Velocity scaling factor: {common['temperature_scaling']['scaling_factor']:.16e}",
        f"- Final scaled temperature: {common['temperature_scaling']['final_temperature_K']:.12f} K",
        f"- Ionic DOF: {common['temperature_scaling']['ndof']}",
        f"- Minimum periodic C-C distance: {common['minimum_periodic_C_C_distance_bohr']:.12f} bohr",
        "",
        "## Parameters",
        "",
        f"- `calculation='vc-md'`, `cell_dynamics='pr'`, `cell_dofree='all'`, `press=0.0`",
        f"- `ion_temperature='nose'`, `tempw=300`, `fnosep(1)=50`, `nhpcl=4`, `wmass=2.0e6`",
        f"- timestep: {DT_FS} fs; requested steps: {NSTEP}; requested time: {NSTEP * DT_FS / 1000.0} ps",
        f"- cutoffs: ecutwfc={ECUTWFC} Ry, ecutrho={ECUTRHO} Ry; Gamma point",
        "",
        "## Control Results",
        "",
    ]
    for key in CONTROL_ORDER:
        val = validations.get(key)
        if not val:
            lines.append(f"- Control {key}: not run or no validation available.")
            continue
        ts = val.get("eal_temperature_K", {})
        terminal = val.get("terminal_direct_temperature", {}).get("temperature", {})
        caimd_ts = val.get("caimd_temperature_K")
        lines.extend(
            [
                f"- Control {key}: {val.get('control_classification')}; steps={val.get('completed_steps')}; JOB DONE={val.get('job_done')}",
                f"  EAL mean/final500/max K: {ts.get('mean')} / {ts.get('final_500_mean')} / {ts.get('max')}",
                f"  Terminal direct `.md` T K: {terminal.get('temperature_K')}",
                f"  CAIMD T mean K: {caimd_ts.get('mean') if caimd_ts else None}",
            ]
        )
    lines.extend(
        [
            "",
            "## No-Op And Projection Diagnostics",
            "",
            f"- A-vs-B maximum numeric difference: {max_ab} ({max_ab_row['observable'] if max_ab_row else None})",
            f"- A-vs-B no-op classified OK: {diag.get('a_vs_b_noop_ok')}",
            f"- Control C maximum cell-strain residual: {c_strain}",
            f"- Control C velocity jump max: {c_jump}",
            f"- Earliest diverging observable: `{diag.get('earliest_diverging_observable')}`",
            "",
            "## Temperature Alignment",
            "",
            "- `.caimd` rows are pre-move diagnostics written before `move_ions`.",
            "- `.eal` rows are written in `vcsmd` after the PR/Beeman update with `nst=istep-1`.",
            "- `analysis/temperature_crosscheck.csv` records zero/one-step offset checks where `.caimd` exists.",
            "",
            "## Required Questions",
            "",
            f"1. Ordinary PW `vc-md` + Nose stable? {answer_stability(validations.get('A'))}",
            f"2. All-false CAIMD valid no-op under Nose? {answer_noop(diag)}",
            f"3. Active all-true zero-strain CAIMD stable? {answer_stability(validations.get('C'))}",
            f"4. High temperature physically present in velocities? See terminal `.md` direct temperatures above and `temperature_crosscheck.csv`.",
            f"5. First divergent control: {diag.get('earliest_diverging_observable')}",
            f"6. Earliest differing physical variable: {diag.get('earliest_diverging_observable')}",
            f"7. Likely issue class: {diag['status']}",
            "",
            "## Local Checks",
            "",
            f"- Python syntax check: {regression_result.get('py_compile', {}).get('returncode') if regression_result else 'not run'}",
            f"- Shell syntax check: {regression_result.get('bash_n', {}).get('returncode') if regression_result else 'not run'}",
            f"- `git diff --check`: {regression_result.get('git_diff_check', {}).get('returncode') if regression_result else 'not run'}",
            f"- Tracked source diff empty: {regression_result.get('tracked_source_diff_empty') if regression_result else 'not checked'}",
            f"- Targeted local CAIMD+Nose harness return code: {regression_result.get('targeted_harness', {}).get('returncode') if regression_result else 'not run'}",
            f"- Full local CAIMD harness: {regression_result.get('full_harness', {}).get('status') if regression_result else 'not run'}",
            "",
            "## Output Files",
            "",
            "- `analysis/temperature_statistics.csv`",
            "- `analysis/temperature_crosscheck.csv`",
            "- `analysis/cell_statistics.csv`",
            "- `analysis/nose_extended_energy_statistics.csv`",
            "- `analysis/control_A_vs_B_differences.csv`",
            "- `analysis/control_C_projection_statistics.csv`",
            "- `analysis/DIAGNOSTIC_CLASSIFICATION.json`",
        ]
    )
    (ROOT / "FINAL_VALIDATION_REPORT.md").write_text("\n".join(lines) + "\n")


def answer_stability(val: dict[str, Any] | None) -> str:
    if not val:
        return "not run"
    cls = val.get("control_classification")
    if cls in ("PASS", "WARN"):
        return f"bounded/completed classification {cls}"
    return f"not stable/classification {cls}"


def answer_noop(diag: dict[str, Any]) -> str:
    if diag.get("a_vs_b_noop_ok"):
        return "yes within strict numeric comparison"
    if diag.get("status") == "BASELINE_NOSE_FAIL":
        return "not determined because Control A failed first"
    return "no or not yet determined"


def analyze() -> None:
    validations: dict[str, dict[str, Any]] = {}
    for key in CONTROL_ORDER:
        run_dir = control_dir(key)
        if any((run_dir / name).exists() for name in ("PROCESS_FINISHED.json", "SAFETY_STOP_CONTEXT.json", "CONTROL_COMPLETE.json")):
            validations[key] = validate_control(key)
    ab_rows = compare_a_b()
    diag = classify_workflow(validations, ab_rows)
    write_analysis_outputs(validations, ab_rows, diag)


def regression() -> None:
    evidence = ROOT / "evidence"
    evidence.mkdir(exist_ok=True)
    results: dict[str, Any] = {"started_utc": utc_stamp()}
    py_compile = run_capture([sys.executable, "-m", "py_compile", str(ROOT / "scripts" / "controls_abc_workflow.py")])
    bash_n = run_capture(["bash", "-n", str(ROOT / "scripts" / "run_controls_abc.sh")])
    diff_check = run_capture(["git", "diff", "--check"])
    source_diff = run_capture(["git", "diff", "--", "PW/src", "Modules", "test-suite"])
    targeted_cases = "ordinary_vcmd_nose,vcmd_allfalse_nose,active_nose_zero,active_nose_stride2_xx"
    harness_cmd = [
        sys.executable,
        str(REPO_ROOT / "test-suite" / "caimd" / "caimd_validation.py"),
        "--mode",
        "local",
        "--current-exe",
        str(PW_EXE),
        "--pseudo",
        str(PSEUDO),
        "--work-root",
        str(evidence),
        "--run-id",
        "caimd_harness_targeted_controls_ABC",
        "--cases",
        targeted_cases,
        "--timeout",
        "300",
    ]
    targeted = run_capture(harness_cmd, cwd=REPO_ROOT, timeout=1200)
    (evidence / "regression_targeted.out").write_text(targeted["stdout"])
    full_note = "Full harness not run automatically after long ABC controls; targeted Nose/CAIMD cases were run."
    results.update(
        {
            "py_compile": {k: v for k, v in py_compile.items() if k != "stdout"},
            "py_compile_stdout": py_compile["stdout"],
            "bash_n": {k: v for k, v in bash_n.items() if k != "stdout"},
            "bash_n_stdout": bash_n["stdout"],
            "git_diff_check": {k: v for k, v in diff_check.items() if k != "stdout"},
            "git_diff_check_stdout": diff_check["stdout"],
            "tracked_source_diff_empty": source_diff["stdout"] == "",
            "tracked_source_diff_stdout": source_diff["stdout"],
            "targeted_harness": {k: v for k, v in targeted.items() if k != "stdout"},
            "targeted_harness_stdout_path": str(evidence / "regression_targeted.out"),
            "full_harness": {"status": "not_run", "reason": full_note},
            "finished_utc": utc_stamp(),
        }
    )
    write_json(evidence / "REGRESSION_RESULTS.json", results)


def status() -> None:
    rows = []
    for key in CONTROL_ORDER:
        run_dir = control_dir(key)
        state = "pending"
        if (run_dir / "CONTROL_COMPLETE.json").exists():
            state = "complete"
        elif (run_dir / "SAFETY_STOP_CONTEXT.json").exists():
            state = "safety_stop"
        elif (run_dir / "PROCESS_FINISHED.json").exists():
            state = "process_finished"
        latest = latest_eal_temperature(run_dir, CONTROLS[key]["prefix"])
        rows.append({"control": key, "state": state, "latest": latest})
    print(json.dumps(rows, indent=2))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare", "run-next", "run-all", "analyze", "regression", "status", "validate"))
    args = parser.parse_args()
    if args.command == "prepare":
        prepare()
    elif args.command == "run-next":
        key = next_incomplete_control()
        if key is None:
            print("No incomplete control remains.")
        else:
            run_control(key)
    elif args.command == "run-all":
        run_all()
    elif args.command == "analyze":
        analyze()
    elif args.command == "regression":
        regression()
    elif args.command == "status":
        status()
    elif args.command == "validate":
        for key in CONTROL_ORDER:
            if (control_dir(key) / "qe.out").exists():
                validate_control(key)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
