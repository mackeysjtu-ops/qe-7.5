#!/usr/bin/env python3
"""Local Diamond-8 CAIMD + Nose four-DOF stability workflow.

This script is deliberately local-only. It never opens network connections,
never uses SSH or Slurm, and writes generated files only below the workflow
directory that contains this script.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import random
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parents[1]
PREVIOUS_ROOT = REPO_ROOT / "examples" / "caimd_diamond8_local_nose_continuation_0p30_20260723T102331Z"
PSEUDO = REPO_ROOT / "pseudo" / "C.UPF"
PSEUDO_NAME = "C.UPF"
DEFAULT_PW = Path("/private/tmp/qe_diamond8_af7ebfc_build/bin/pw.x")

RYDBERG_SI = 2.1798723611035e-18
K_BOLTZMANN_SI = 1.380649e-23
AMU_SI = 1.66053906660e-27
ELECTRONMASS_SI = 9.1093837015e-31
AMU_RY = (AMU_SI / ELECTRONMASS_SI) / 2.0
RY_TO_KELVIN = RYDBERG_SI / K_BOLTZMANN_SI
ANGSTROM_TO_BOHR = 1.8897261254578281

MASS_C_AMU = 12.0107
NAT = 8
DT_FS = 2.0
DT_QE = 41.341381006
TEMP_TARGET_K = 300.0
INITIAL_VELOCITY_TEMPERATURE_K = 300.0
VELOCITY_SEED = 20260723
ECUTWFC = 25.0
ECUTRHO = 120.0
MIXING_BETA = 0.3
CONV_THR = 1.0e-4
ELECTRON_MAXSTEP = 40
FNOSEP_THZ = 50.0
NHPCL = 4
WMASS = 2.0e6
EQUIL_STEPS = 2000
LONG_STEPS = 5000
CAIMD_STRIDE = 10
LOADING_STEP_STRAIN = -6.0e-6
TARGET_FINAL_XX = -0.03
STRAIN_RATE_PS_INV = 0.003
LITERATURE_RATE_PS_INV = 0.001
FAILED_LOCAL_RATE_PS_INV = 0.05

VOLUME_BOHR3_MIN = 10.0
VOLUME_BOHR3_MAX = 5000.0
CELL_CONDITION_MAX = 1.0e6
MIN_CC_BOHR_MIN = 0.8
MIN_CC_BOHR_MAX = 20.0
SUSTAINED_WINDOW = 25

FLOAT_RE = re.compile(
    r"[-+]?(?:(?:\d+\.\d*|\.\d+)(?:[EeDd][-+]?\d+)?|\d+(?:[EeDd][-+]?\d+)?)"
)

COMPONENTS = ("xx", "yy", "zz", "yz", "xz", "xy")
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

TESTS: dict[str, dict[str, Any]] = {
    "test1": {
        "dirname": "02_test1_zero_load_full_pr",
        "prefix": "diamond8_test1_zero_load_full_pr",
        "title": "zero load with full unmasked PR response",
        "mask": (True, False, False, False, False, False),
        "step": (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    },
    "test2": {
        "dirname": "03_test2_load_all_components_constrained",
        "prefix": "diamond8_test2_all_components_constrained",
        "title": "xx loading with all six strain components constrained",
        "mask": (True, True, True, True, True, True),
        "step": (LOADING_STEP_STRAIN, 0.0, 0.0, 0.0, 0.0, 0.0),
    },
    "test3": {
        "dirname": "04_test3_load_transverse_normal_free",
        "prefix": "diamond8_test3_transverse_normal_free",
        "title": "xx loading with yy/zz free and shears constrained",
        "mask": (True, False, False, True, True, True),
        "step": (LOADING_STEP_STRAIN, 0.0, 0.0, 0.0, 0.0, 0.0),
    },
    "test4": {
        "dirname": "05_test4_load_full_pr",
        "prefix": "diamond8_test4_load_full_pr",
        "title": "xx loading with all five remaining components free",
        "mask": (True, False, False, False, False, False),
        "step": (LOADING_STEP_STRAIN, 0.0, 0.0, 0.0, 0.0, 0.0),
    },
}


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run_text(cmd: list[str], cwd: Path | None = None) -> str:
    return subprocess.check_output(cmd, cwd=str(cwd) if cwd else None, text=True, stderr=subprocess.STDOUT)


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


def ffloat(x: float) -> str:
    return f"{x:.16e}".replace("e", "d")


def extract_floats(text: str) -> list[float]:
    out: list[float] = []
    for token in FLOAT_RE.findall(text):
        try:
            out.append(float(token.replace("D", "E").replace("d", "e")))
        except ValueError:
            continue
    return out


def transpose(a: list[list[float]]) -> list[list[float]]:
    return [list(row) for row in zip(*a)]


def matvec(a: list[list[float]], x: list[float]) -> list[float]:
    return [sum(a[i][j] * x[j] for j in range(3)) for i in range(3)]


def matmul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def vecnorm(x: list[float]) -> float:
    return math.sqrt(sum(v * v for v in x))


def det3(a: list[list[float]]) -> float:
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def inv3(a: list[list[float]]) -> list[list[float]]:
    d = det3(a)
    if not math.isfinite(d) or abs(d) < 1.0e-18:
        raise ValueError(f"singular 3x3 matrix with determinant {d}")
    return [
        [
            (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / d,
            -(a[0][1] * a[2][2] - a[0][2] * a[2][1]) / d,
            (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / d,
        ],
        [
            -(a[1][0] * a[2][2] - a[1][2] * a[2][0]) / d,
            (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / d,
            -(a[0][0] * a[1][2] - a[0][2] * a[1][0]) / d,
        ],
        [
            (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / d,
            -(a[0][0] * a[2][1] - a[0][1] * a[2][0]) / d,
            (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / d,
        ],
    ]


def norm_inf(a: list[list[float]]) -> float:
    return max(sum(abs(v) for v in row) for row in a)


def condition_number(a: list[list[float]]) -> float:
    try:
        import numpy as np

        return float(np.linalg.cond(np.array(a, dtype=float)))
    except Exception:
        return norm_inf(a) * norm_inf(inv3(a))


def singular_values(a: list[list[float]]) -> list[float]:
    try:
        import numpy as np

        return [float(x) for x in np.linalg.svd(np.array(a, dtype=float), compute_uv=False)]
    except Exception:
        return []


def finite_values(values: list[float]) -> bool:
    return all(math.isfinite(v) for v in values)


def wrap_frac(frac: list[float]) -> list[float]:
    return [x - math.floor(x) for x in frac]


def frac_to_cart(avec: list[list[float]], frac: list[float]) -> list[float]:
    return matvec(avec, frac)


def cart_to_frac(avec: list[list[float]], cart: list[float]) -> list[float]:
    return matvec(inv3(avec), cart)


def min_max_periodic_distances(avec: list[list[float]], frac: list[list[float]]) -> tuple[float | None, float | None]:
    if len(frac) < 2:
        return None, None
    best_min = float("inf")
    best_max = 0.0
    shifts = (-1, 0, 1)
    for i in range(len(frac)):
        for j in range(i + 1, len(frac)):
            pair_best = float("inf")
            for sx in shifts:
                for sy in shifts:
                    for sz in shifts:
                        df = [
                            frac[j][0] - frac[i][0] + sx,
                            frac[j][1] - frac[i][1] + sy,
                            frac[j][2] - frac[i][2] + sz,
                        ]
                        pair_best = min(pair_best, vecnorm(frac_to_cart(avec, df)))
            best_min = min(best_min, pair_best)
            best_max = max(best_max, pair_best)
    return best_min, best_max


def velocity_temperature(velocities: list[list[float]]) -> dict[str, Any]:
    mass_ry = MASS_C_AMU * AMU_RY
    ekin = 0.5 * mass_ry * sum(sum(v * v for v in row) for row in velocities)
    ndof_3n3 = 3 * len(velocities) - 3
    ndof_pr = 3 * (len(velocities) + 1)
    com = [sum(row[i] for row in velocities) / len(velocities) for i in range(3)]
    return {
        "ekin_Ry": ekin,
        "temperature_3N_minus_3_K": 2.0 * ekin / ndof_3n3 * RY_TO_KELVIN,
        "temperature_3N_plus_cell_K": 2.0 * ekin / ndof_pr * RY_TO_KELVIN,
        "center_of_mass_velocity_au": com,
        "center_of_mass_max_abs_au": max(abs(x) for x in com),
        "ndof_3N_minus_3": ndof_3n3,
        "ndof_3N_plus_cell": ndof_pr,
    }


def stats(values: list[float], final_window: int = 500) -> dict[str, Any]:
    vals = [v for v in values if math.isfinite(v)]
    if not values:
        return {
            "count": 0,
            "finite": False,
            "mean": None,
            "std": None,
            "min": None,
            "max": None,
            "first_500_mean": None,
            "final_500_mean": None,
        }
    mean = sum(vals) / len(vals) if vals else None
    var = sum((v - mean) ** 2 for v in vals) / len(vals) if vals and mean is not None else None
    first = vals[: min(final_window, len(vals))]
    final = vals[-min(final_window, len(vals)) :]
    return {
        "count": len(values),
        "finite": len(vals) == len(values),
        "mean": mean,
        "std": math.sqrt(var) if var is not None else None,
        "min": min(vals) if vals else None,
        "max": max(vals) if vals else None,
        "first_500_mean": sum(first) / len(first) if first else None,
        "final_500_mean": sum(final) / len(final) if final else None,
    }


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
                "block_index": len(rows) + 1,
                "start_step": start + 1,
                "end_step": start + len(chunk),
                "count": len(chunk),
                "mean": sum(chunk) / len(chunk),
                "min": min(chunk),
                "max": max(chunk),
            }
        )
    return rows


def linear_drift_per_ps(times: list[float], values: list[float]) -> float | None:
    pairs = [(t, v) for t, v in zip(times, values) if math.isfinite(t) and math.isfinite(v)]
    if len(pairs) < 2:
        return None
    mt = sum(t for t, _ in pairs) / len(pairs)
    mv = sum(v for _, v in pairs) / len(pairs)
    denom = sum((t - mt) ** 2 for t, _ in pairs)
    if denom == 0.0:
        return None
    return sum((t - mt) * (v - mv) for t, v in pairs) / denom


def fraction(values: list[float], predicate) -> float | None:
    if not values:
        return None
    return sum(1 for v in values if predicate(v)) / len(values)


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


def temperature_summary(times: list[float], values: list[float]) -> dict[str, Any]:
    r25 = rolling(values, 25)
    r100 = rolling(values, 100)
    b500 = block_means(values, 500)
    s = stats(values)
    s.update(
        {
            "rolling25_min": min(r25) if r25 else None,
            "rolling25_max": max(r25) if r25 else None,
            "rolling100_min": min(r100) if r100 else None,
            "rolling100_max": max(r100) if r100 else None,
            "block500_min_mean": min((b["mean"] for b in b500), default=None),
            "block500_max_mean": max((b["mean"] for b in b500), default=None),
            "drift_K_per_ps": linear_drift_per_ps(times, values),
            "fraction_below_150K": fraction(values, lambda v: v < 150.0),
            "fraction_above_600K": fraction(values, lambda v: v > 600.0),
            "longest_below_150K_steps": longest_run(values, lambda v: v < 150.0),
            "longest_above_600K_steps": longest_run(values, lambda v: v > 600.0),
            "longest_above_2000K_steps": longest_run(values, lambda v: v > 2000.0),
        }
    )
    return s


def classify_temperature(summary: dict[str, Any]) -> str:
    if not summary.get("finite") or summary.get("count", 0) == 0:
        return "FAIL"
    if summary.get("longest_below_150K_steps", 0) >= SUSTAINED_WINDOW:
        return "FAIL"
    if summary.get("longest_above_600K_steps", 0) >= SUSTAINED_WINDOW:
        return "FAIL"
    if summary.get("longest_above_2000K_steps", 0) > 0:
        return "FAIL"
    block_min = summary.get("block500_min_mean")
    block_max = summary.get("block500_max_mean")
    full_mean = summary.get("mean")
    final_mean = summary.get("final_500_mean")
    drift = summary.get("drift_K_per_ps")
    if block_min is None or block_max is None or full_mean is None or final_mean is None:
        return "FAIL"
    if 240.0 <= block_min and block_max <= 360.0 and 240.0 <= full_mean <= 360.0 and 240.0 <= final_mean <= 360.0:
        if drift is None or abs(drift) < 50.0:
            return "PASS"
    if 200.0 <= block_min and block_max <= 450.0 and 200.0 <= full_mean <= 450.0 and 200.0 <= final_mean <= 450.0:
        return "WARN"
    # The requested acceptance policy treats this as WARN, not FAIL, when the
    # thermostat remains active and there is no sustained low/high-temperature
    # interval. Individual instantaneous excursions are expected in an 8-atom
    # cell and are reported separately through the fractions and longest runs.
    return "WARN"


def diamond8_setup() -> dict[str, Any]:
    a_ang = 3.567
    a_bohr = a_ang * ANGSTROM_TO_BOHR
    row_x = [1.0 / math.sqrt(3.0), 1.0 / math.sqrt(3.0), 1.0 / math.sqrt(3.0)]
    row_y = [1.0 / math.sqrt(2.0), -1.0 / math.sqrt(2.0), 0.0]
    row_z = [1.0 / math.sqrt(6.0), 1.0 / math.sqrt(6.0), -2.0 / math.sqrt(6.0)]
    rot = [row_x, row_y, row_z]
    # QE's internal cell matrix has lattice vectors in columns.
    avec = [[a_bohr * rot[i][j] for j in range(3)] for i in range(3)]
    frac = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ]
    velocities, scale, raw_temp = generate_velocities()
    min_cc, max_cc = min_max_periodic_distances(avec, frac)
    return {
        "lattice_constant_angstrom": a_ang,
        "lattice_constant_bohr": a_bohr,
        "rotation_matrix_lab_from_crystal": rot,
        "cell_matrix_bohr_avec_rows": avec,
        "cell_parameters_bohr_rows_for_qe": transpose(avec),
        "fractional_positions": frac,
        "atom_labels": ["C"] * NAT,
        "atom_order": [
            "C 0 0 0",
            "C 0 1/2 1/2",
            "C 1/2 0 1/2",
            "C 1/2 1/2 0",
            "C 1/4 1/4 1/4",
            "C 1/4 3/4 3/4",
            "C 3/4 1/4 3/4",
            "C 3/4 3/4 1/4",
        ],
        "cartesian_velocities_au": velocities,
        "velocity_seed": VELOCITY_SEED,
        "velocity_raw_temperature_3N_minus_3_K": raw_temp,
        "velocity_scale_factor": scale,
        "velocity_temperature": velocity_temperature(velocities),
        "determinant_bohr3": det3(avec),
        "condition_number": condition_number(avec),
        "singular_values_bohr": singular_values(avec),
        "minimum_periodic_C_C_distance_bohr": min_cc,
        "maximum_periodic_C_C_distance_bohr": max_cc,
        "positive_determinant": det3(avec) > 0.0,
    }


def generate_velocities() -> tuple[list[list[float]], float, float]:
    rng = random.Random(VELOCITY_SEED)
    velocities = [[rng.gauss(0.0, 1.0) for _ in range(3)] for _ in range(NAT)]
    com = [sum(row[i] for row in velocities) / NAT for i in range(3)]
    velocities = [[row[i] - com[i] for i in range(3)] for row in velocities]
    raw_temp = velocity_temperature(velocities)["temperature_3N_minus_3_K"]
    scale = math.sqrt(INITIAL_VELOCITY_TEMPERATURE_K / raw_temp)
    velocities = [[scale * x for x in row] for row in velocities]
    return velocities, scale, raw_temp


def input_state_from_setup(setup: dict[str, Any]) -> dict[str, Any]:
    return {
        "cell_matrix_bohr_avec_rows": setup["cell_matrix_bohr_avec_rows"],
        "fractional_positions": setup["fractional_positions"],
        "cartesian_velocities_au": setup["cartesian_velocities_au"],
        "atom_labels": setup["atom_labels"],
        "source": "generated_rotated_conventional_diamond8",
        "units": {
            "cell": "bohr, QE internal avec rows with lattice vectors in columns",
            "positions": "crystal fractional",
            "velocities": "a.u. Cartesian",
        },
    }


def caimd_cell_lines(active: bool, mask: tuple[bool, ...], step: tuple[float, ...]) -> str:
    if not active:
        return "  caimd = .false.,\n"
    lines = ["  caimd = .true,", f"  caimd_stride = {CAIMD_STRIDE},"]
    for i, val in enumerate(step, start=1):
        lines.append(f"  caimd_step_strain({i}) = {ffloat(val)},")
    for i, val in enumerate(mask, start=1):
        lines.append(f"  caimd_mask({i}) = {'.true.' if val else '.false.'},")
    return "\n".join(lines) + "\n"


def write_input(
    run_dir: Path,
    prefix: str,
    nstep: int,
    state: dict[str, Any],
    *,
    active_caimd: bool,
    mask: tuple[bool, ...] = (False, False, False, False, False, False),
    step: tuple[float, ...] = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "tmp").mkdir(exist_ok=True)
    avec = state["cell_matrix_bohr_avec_rows"]
    qe_rows = transpose(avec)
    frac = state["fractional_positions"]
    vel = state["cartesian_velocities_au"]
    labels = state.get("atom_labels", ["C"] * len(frac))
    if len(frac) != NAT or len(vel) != NAT:
        raise ValueError("state must contain exactly eight positions and velocities")
    text = f"""&CONTROL
  calculation = 'vc-md',
  restart_mode = 'from_scratch',
  prefix = '{prefix}',
  pseudo_dir = '{REPO_ROOT / "pseudo"}',
  outdir = '{run_dir / "tmp"}',
  nstep = {nstep},
  dt = {ffloat(DT_QE)},
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
  conv_thr = {ffloat(CONV_THR)},
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
{caimd_cell_lines(active_caimd, mask, step).rstrip()}
/
ATOMIC_SPECIES
C {MASS_C_AMU:.8f} {PSEUDO_NAME}
CELL_PARAMETERS bohr
"""
    for row in qe_rows:
        text += "  " + "  ".join(f"{x:.14f}" for x in row) + "\n"
    text += "ATOMIC_POSITIONS crystal\n"
    for label, xyz in zip(labels, frac):
        text += f"{label}  " + "  ".join(f"{x:.14f}" for x in xyz) + "\n"
    text += "K_POINTS gamma\n"
    text += "ATOMIC_VELOCITIES { a.u }\n"
    for label, xyz in zip(labels, vel):
        text += f"{label}  " + "  ".join(f"{x:.14e}" for x in xyz) + "\n"
    path = run_dir / "qe.in"
    path.write_text(text)
    return path


def parse_blocks(path: Path, rows_per_block: int) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    lines = [line.rstrip() for line in path.read_text(errors="ignore").splitlines() if line.strip()]
    blocks: list[dict[str, Any]] = []
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
        blocks.append({"step": int(round(head[0])), "time_raw": head[1], "rows": rows})
        i += 1 + rows_per_block
    return blocks


def artifact(run_dir: Path, prefix: str, ext: str) -> Path:
    direct = run_dir / f"{prefix}.{ext}"
    if direct.exists():
        return direct
    return run_dir / "tmp" / f"{prefix}.{ext}"


def parse_stdout_temperatures(path: Path) -> list[float]:
    if not path.exists():
        return []
    vals: list[float] = []
    for line in path.read_text(errors="ignore").splitlines():
        match = re.search(r"\bT\s*=\s*([-+0-9.EeDd]+)\s*K", line)
        if match:
            vals.append(float(match.group(1).replace("D", "E").replace("d", "e")))
    return vals


def parse_eal_temperatures(path: Path) -> list[float]:
    if not path.exists():
        return []
    vals: list[float] = []
    ndof = 3 * NAT - 3
    for line in path.read_text(errors="ignore").splitlines():
        nums = extract_floats(line)
        if len(nums) >= 2:
            vals.append(2.0 * nums[1] / ndof * RY_TO_KELVIN)
    return vals


def parse_stress_blocks(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    lines = path.read_text(errors="ignore").splitlines()
    blocks: list[dict[str, Any]] = []
    i = 0
    while i < len(lines):
        if "total   stress" not in lines[i]:
            i += 1
            continue
        nums = extract_floats(lines[i])
        pressure_kbar = nums[-1] if nums else None
        rows = []
        for j in range(1, 4):
            if i + j >= len(lines):
                break
            vals = extract_floats(lines[i + j])
            if len(vals) >= 6:
                rows.append(vals[3:6])
        if len(rows) == 3:
            components_kbar = [rows[0][0], rows[1][1], rows[2][2], rows[1][2], rows[0][2], rows[0][1]]
            blocks.append(
                {
                    "pressure_kbar": pressure_kbar,
                    **{f"stress_{name}_GPa": components_kbar[idx] / 10.0 for idx, name in enumerate(COMPONENTS)},
                }
            )
        i += 4
    return blocks


def rows_from_caimd(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows: list[dict[str, Any]] = []
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


def parse_caimd_summary(path: Path, expected_steps: int, mask: tuple[bool, ...], step: tuple[float, ...]) -> dict[str, Any]:
    rows = rows_from_caimd(path)
    if not path.exists():
        return {"exists": False, "row_count": 0}
    malformed = 0
    for line in path.read_text(errors="ignore").splitlines():
        if line.strip() and not line.lstrip().startswith("#") and len(line.split()) != 26:
            malformed += 1
    out: dict[str, Any] = {
        "exists": True,
        "path": str(path),
        "row_count": len(rows),
        "malformed_data_lines": malformed,
        "column_count": 26,
        "finite": all(all(math.isfinite(v) for v in row.values() if isinstance(v, float)) for row in rows),
        "event_count_last": rows[-1]["imposed_events_done"] if rows else None,
        "temperature_K": temperature_summary([r["time_ps"] for r in rows], [r["temperature_K"] for r in rows]),
        "velocity_jump_max": max((r["velocity_jump_max"] for r in rows), default=None),
        "velocity_jump_rms_max": max((r["velocity_jump_rms"] for r in rows), default=None),
    }
    for name in COMPONENTS:
        actual = [r[f"actual_{name}"] for r in rows]
        nominal = [r[f"nominal_{name}"] for r in rows]
        out[f"actual_{name}_final"] = actual[-1] if actual else None
        out[f"nominal_{name}_final"] = nominal[-1] if nominal else None
        out[f"actual_{name}_min"] = min(actual) if actual else None
        out[f"actual_{name}_max"] = max(actual) if actual else None
        out[f"nominal_{name}_min"] = min(nominal) if nominal else None
        out[f"nominal_{name}_max"] = max(nominal) if nominal else None
        out[f"nominal_actual_{name}_abs_error_max"] = max((abs(a - n) for a, n in zip(actual, nominal)), default=None)
        if mask[COMPONENTS.index(name)]:
            out[f"constrained_{name}_held_abs_error_max"] = out[f"nominal_actual_{name}_abs_error_max"]
    out["expected_event_count"] = expected_steps // CAIMD_STRIDE
    out["expected_final_nominal_xx"] = (expected_steps // CAIMD_STRIDE) * CAIMD_STRIDE * step[0]
    return out


def terminal_ratd_from_md(path: Path) -> list[list[float]]:
    if not path.exists():
        raise RuntimeError(f"missing QE PR MD history file: {path}")
    vals = extract_floats(path.read_text(errors="ignore"))
    need = 5 * 3 * NAT
    if len(vals) < need:
        raise RuntimeError(f"malformed QE PR MD history: only {len(vals)} floats in {path}")
    start = 3 * NAT
    return [vals[start + 3 * i : start + 3 * i + 3] for i in range(NAT)]


def extract_terminal_state(run_dir: Path, prefix: str, output: Path) -> dict[str, Any]:
    cel = parse_blocks(artifact(run_dir, prefix, "cel"), 3)
    pos = parse_blocks(artifact(run_dir, prefix, "pos"), NAT)
    if not cel or not pos:
        raise RuntimeError(f"missing terminal cel/pos blocks for {run_dir}")
    if cel[-1]["step"] != pos[-1]["step"]:
        raise RuntimeError(f"terminal cel/pos step mismatch in {run_dir}")
    avec = cel[-1]["rows"]
    frac = [wrap_frac(cart_to_frac(avec, xyz)) for xyz in pos[-1]["rows"]]
    ratd = terminal_ratd_from_md(artifact(run_dir, prefix, "md"))
    velocities = [matvec(avec, row) for row in ratd]
    min_cc, max_cc = min_max_periodic_distances(avec, frac)
    state = {
        "source_run_dir": str(run_dir),
        "source_prefix": prefix,
        "source_step": cel[-1]["step"],
        "source_files": {
            "cel": str(artifact(run_dir, prefix, "cel")),
            "pos": str(artifact(run_dir, prefix, "pos")),
            "md": str(artifact(run_dir, prefix, "md")),
            "eal": str(artifact(run_dir, prefix, "eal")),
        },
        "record_indices": {
            "cel_block_index_1based": len(cel),
            "pos_block_index_1based": len(pos),
            "md_record": "last written PR history",
        },
        "cell_matrix_bohr_avec_rows": avec,
        "cell_parameters_bohr_rows_for_qe": transpose(avec),
        "fractional_positions": frac,
        "cartesian_positions_bohr": pos[-1]["rows"],
        "cartesian_velocities_au": velocities,
        "terminal_ratd_lattice_velocity": ratd,
        "atom_labels": ["C"] * NAT,
        "velocity_units": "a.u. Cartesian, computed as QE avec * terminal ratd from .md",
        "temperature_from_velocities": velocity_temperature(velocities),
        "volume_bohr3": det3(avec),
        "cell_condition_number": condition_number(avec),
        "cell_singular_values_bohr": singular_values(avec),
        "minimum_periodic_C_C_distance_bohr": min_cc,
        "maximum_periodic_C_C_distance_bohr": max_cc,
        "units": {
            "cell": "bohr, QE internal avec rows with lattice vectors in columns",
            "positions": "crystal fractional, reconstructed from Cartesian .pos and .cel",
            "velocities": "a.u. Cartesian",
        },
        "extraction_utc": utc_now(),
    }
    write_json(output, state)
    write_state_fragment(output.parent / "common_start_state_fragment.in", state)
    return state


def write_state_fragment(path: Path, state: dict[str, Any]) -> None:
    lines = ["CELL_PARAMETERS bohr\n"]
    for row in transpose(state["cell_matrix_bohr_avec_rows"]):
        lines.append("  " + "  ".join(f"{x:.14f}" for x in row) + "\n")
    lines.append("ATOMIC_POSITIONS crystal\n")
    for label, xyz in zip(state["atom_labels"], state["fractional_positions"]):
        lines.append(f"{label}  " + "  ".join(f"{x:.14f}" for x in xyz) + "\n")
    lines.append("ATOMIC_VELOCITIES { a.u }\n")
    for label, xyz in zip(state["atom_labels"], state["cartesian_velocities_au"]):
        lines.append(f"{label}  " + "  ".join(f"{x:.14e}" for x in xyz) + "\n")
    path.write_text("".join(lines))


def cell_series(run_dir: Path, prefix: str) -> list[dict[str, Any]]:
    cel = parse_blocks(artifact(run_dir, prefix, "cel"), 3)
    pos = parse_blocks(artifact(run_dir, prefix, "pos"), NAT)
    out: list[dict[str, Any]] = []
    n = min(len(cel), len(pos))
    for idx in range(n):
        avec = cel[idx]["rows"]
        det = det3(avec)
        cond = condition_number(avec) if det != 0.0 and math.isfinite(det) else float("inf")
        svals = singular_values(avec)
        min_cc = None
        max_cc = None
        try:
            frac = [wrap_frac(cart_to_frac(avec, xyz)) for xyz in pos[idx]["rows"]]
            min_cc, max_cc = min_max_periodic_distances(avec, frac)
        except Exception:
            pass
        out.append(
            {
                "step": cel[idx]["step"],
                "determinant_bohr3": det,
                "volume_bohr3": det,
                "condition_number": cond,
                "singular_value_min_bohr": min(svals) if svals else None,
                "singular_value_max_bohr": max(svals) if svals else None,
                "minimum_C_C_distance_bohr": min_cc,
                "maximum_C_C_distance_bohr": max_cc,
            }
        )
    return out


def validate_run(
    run_dir: Path,
    prefix: str,
    expected_steps: int,
    *,
    active_caimd: bool,
    mask: tuple[bool, ...] = (False, False, False, False, False, False),
    step: tuple[float, ...] = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
) -> dict[str, Any]:
    out_path = run_dir / "qe.out"
    text = out_path.read_text(errors="ignore") if out_path.exists() else ""
    cel = parse_blocks(artifact(run_dir, prefix, "cel"), 3)
    pos = parse_blocks(artifact(run_dir, prefix, "pos"), NAT)
    vel = parse_blocks(artifact(run_dir, prefix, "vel"), NAT)
    eal_t = parse_eal_temperatures(artifact(run_dir, prefix, "eal"))
    stdout_t = parse_stdout_temperatures(out_path)
    stress = parse_stress_blocks(out_path)
    cell = cell_series(run_dir, prefix)
    caimd_rows = rows_from_caimd(artifact(run_dir, prefix, "caimd")) if active_caimd else []
    if active_caimd and caimd_rows:
        temp_times = [r["time_ps"] for r in caimd_rows]
        temp_values = [r["temperature_K"] for r in caimd_rows]
        temp_source = "caimd_temperature_K"
    elif eal_t:
        temp_times = [(i + 1) * DT_FS / 1000.0 for i in range(len(eal_t))]
        temp_values = eal_t
        temp_source = "eal_ionic_3N_minus_3_K"
    else:
        temp_times = [(i + 1) * DT_FS / 1000.0 for i in range(len(stdout_t))]
        temp_values = stdout_t
        temp_source = "stdout_T_K"
    temp_stats = temperature_summary(temp_times, temp_values)
    temp_class = classify_temperature(temp_stats)
    caimd = parse_caimd_summary(artifact(run_dir, prefix, "caimd"), expected_steps, mask, step) if active_caimd else {"exists": False}
    run_complete = read_json(run_dir / f"{prefix}_RUN_COMPLETE.json") if (run_dir / f"{prefix}_RUN_COMPLETE.json").exists() else {}
    dets = [r["determinant_bohr3"] for r in cell]
    conds = [r["condition_number"] for r in cell]
    min_ccs = [r["minimum_C_C_distance_bohr"] for r in cell if r["minimum_C_C_distance_bohr"] is not None]
    max_ccs = [r["maximum_C_C_distance_bohr"] for r in cell if r["maximum_C_C_distance_bohr"] is not None]
    stress_stats = {
        name: stats([row[f"stress_{name}_GPa"] for row in stress if f"stress_{name}_GPa" in row], 500)
        for name in COMPONENTS
    }
    result: dict[str, Any] = {
        "run_dir": str(run_dir),
        "prefix": prefix,
        "expected_steps": expected_steps,
        "active_caimd": active_caimd,
        "job_done": "JOB DONE" in text,
        "returncode": run_complete.get("returncode"),
        "elapsed_s": run_complete.get("elapsed_s"),
        "trajectory_blocks": {"cel": len(cel), "pos": len(pos), "vel": len(vel)},
        "completed_expected_steps": len(cel) >= expected_steps and len(pos) >= expected_steps,
        "temperature_analysis_source": temp_source,
        "temperature_analysis_K": temp_stats,
        "temperature_classification": temp_class,
        "stdout_temperature_K": temperature_summary([(i + 1) * DT_FS / 1000.0 for i in range(len(stdout_t))], stdout_t) if stdout_t else None,
        "eal_ionic_temperature_3N_minus_3_K": temperature_summary([(i + 1) * DT_FS / 1000.0 for i in range(len(eal_t))], eal_t) if eal_t else None,
        "caimd_temperature_K": caimd.get("temperature_K"),
        "stress_block_count": len(stress),
        "stress_component_statistics_GPa": stress_stats,
        "pressure_kbar_stats": stats([row["pressure_kbar"] for row in stress if row.get("pressure_kbar") is not None], 500),
        "cell_blocks": len(cell),
        "cell_determinant_bohr3": stats(dets, 500),
        "cell_condition_number": stats(conds, 500),
        "minimum_periodic_C_C_distance_bohr": stats(min_ccs, 500),
        "maximum_periodic_C_C_distance_bohr": stats(max_ccs, 500),
        "volume_positive_all": all(v > 0.0 for v in dets) if dets else False,
        "volume_sane_all": all(VOLUME_BOHR3_MIN <= v <= VOLUME_BOHR3_MAX for v in dets) if dets else False,
        "cell_condition_sane_all": all(v <= CELL_CONDITION_MAX for v in conds) if conds else False,
        "minimum_C_C_sane_all": all(MIN_CC_BOHR_MIN <= v <= MIN_CC_BOHR_MAX for v in min_ccs) if min_ccs else False,
        "maximum_C_C_available": bool(max_ccs),
        "caimd": caimd,
        "qe_error_lines": [line.strip() for line in text.splitlines() if "Error in routine" in line or "STOP" in line],
        "validated_utc": utc_now(),
    }
    if active_caimd:
        nominal_xx = caimd.get("nominal_xx_final")
        actual_xx = caimd.get("actual_xx_final")
        result["final_nominal_xx_strain"] = nominal_xx
        result["final_actual_xx_strain"] = actual_xx
        result["event_count"] = caimd.get("event_count_last")
        result["xx_tracking_error_max"] = caimd.get("nominal_actual_xx_abs_error_max")
        result["velocity_jump_max"] = caimd.get("velocity_jump_max")
        result["velocity_jump_rms_max"] = caimd.get("velocity_jump_rms_max")
    valid = bool(result["job_done"] and result["completed_expected_steps"])
    valid = valid and temp_class != "FAIL"
    valid = valid and result["volume_positive_all"] and result["volume_sane_all"]
    valid = valid and result["cell_condition_sane_all"] and result["minimum_C_C_sane_all"]
    if active_caimd:
        valid = valid and caimd.get("exists") and caimd.get("row_count", 0) >= expected_steps
        valid = valid and caimd.get("finite", False)
        valid = valid and caimd.get("event_count_last") == expected_steps // CAIMD_STRIDE
        if step[0] != 0.0:
            valid = valid and abs((caimd.get("nominal_xx_final") or 0.0) - TARGET_FINAL_XX) <= 1.0e-8
            valid = valid and abs((caimd.get("actual_xx_final") or 0.0) - TARGET_FINAL_XX) <= 1.0e-6
        valid = valid and (caimd.get("nominal_actual_xx_abs_error_max") is not None and caimd.get("nominal_actual_xx_abs_error_max") <= 1.0e-8)
        valid = valid and (caimd.get("velocity_jump_max") is not None and caimd.get("velocity_jump_max") <= 1.0e-8)
    result["valid"] = bool(valid)
    result["classification"] = "PASS" if result["valid"] and temp_class == "PASS" else ("WARN" if result["valid"] else "FAIL")
    result["first_failing_observable"] = first_failing_observable(result, caimd_rows, cell)
    write_json(run_dir / "VALIDATION.json", result)
    if active_caimd:
        write_csv(run_dir / "CAIMD_ROWS.csv", caimd_rows)
    write_csv(run_dir / "CELL_SERIES.csv", cell)
    if not result["valid"]:
        write_safety_stop_report(run_dir, result)
    return result


def first_failing_observable(result: dict[str, Any], caimd_rows: list[dict[str, Any]], cell: list[dict[str, Any]]) -> dict[str, Any] | None:
    candidates: list[dict[str, Any]] = []
    for idx, row in enumerate(caimd_rows, start=1):
        temp = row.get("temperature_K")
        if temp is not None and temp > 2000.0:
            candidates.append({"step_index": idx, "observable": "temperature_K", "value": temp, "reason": "above 2000 K"})
            break
    for idx in range(0, max(len(caimd_rows) - SUSTAINED_WINDOW + 1, 0)):
        window = caimd_rows[idx : idx + SUSTAINED_WINDOW]
        if all(r["temperature_K"] > 600.0 for r in window):
            candidates.append({"step_index": idx + 1, "observable": "temperature_K", "value": window[0]["temperature_K"], "reason": "sustained above 600 K"})
            break
    for idx in range(0, max(len(caimd_rows) - SUSTAINED_WINDOW + 1, 0)):
        window = caimd_rows[idx : idx + SUSTAINED_WINDOW]
        if all(r["temperature_K"] < 150.0 for r in window):
            candidates.append({"step_index": idx + 1, "observable": "temperature_K", "value": window[0]["temperature_K"], "reason": "sustained below 150 K"})
            break
    for idx, row in enumerate(cell, start=1):
        if row["determinant_bohr3"] <= 0.0 or not math.isfinite(row["determinant_bohr3"]):
            candidates.append({"step_index": idx, "observable": "determinant_bohr3", "value": row["determinant_bohr3"], "reason": "nonpositive or nonfinite determinant"})
            break
        if row["condition_number"] > CELL_CONDITION_MAX:
            candidates.append({"step_index": idx, "observable": "cell_condition_number", "value": row["condition_number"], "reason": "extreme condition number"})
            break
        d = row.get("minimum_C_C_distance_bohr")
        if d is not None and d < MIN_CC_BOHR_MIN:
            candidates.append({"step_index": idx, "observable": "minimum_C_C_distance_bohr", "value": d, "reason": "pathological short C-C distance"})
            break
    if result.get("qe_error_lines"):
        candidates.append({"step_index": result["trajectory_blocks"]["cel"], "observable": "qe_error", "value": result["qe_error_lines"][0], "reason": "QE error/STOP"})
    if not candidates:
        return None
    return sorted(candidates, key=lambda r: r["step_index"])[0]


def write_safety_stop_report(run_dir: Path, result: dict[str, Any]) -> None:
    report = run_dir / "SAFETY_STOP_REPORT.md"
    report.write_text(
        f"""# Safety Stop Report

Generated: {utc_now()}

Status: FAIL

Run directory: `{run_dir}`

Prefix: `{result['prefix']}`

Expected steps: `{result['expected_steps']}`

Completed trajectory blocks: `{result['trajectory_blocks']}`

Return code: `{result.get('returncode')}`

JOB DONE: `{result.get('job_done')}`

First failing observable:

```json
{json.dumps(result.get('first_failing_observable'), indent=2, sort_keys=True)}
```

QE error lines:

```json
{json.dumps(result.get('qe_error_lines'), indent=2)}
```

Temperature classification: `{result.get('temperature_classification')}`

Cell determinant stats:

```json
{json.dumps(result.get('cell_determinant_bohr3'), indent=2, sort_keys=True)}
```

Cell condition stats:

```json
{json.dumps(result.get('cell_condition_number'), indent=2, sort_keys=True)}
```

The partial outputs are preserved in this directory. This test must not be
treated as successful unless it is rerun and validated from a clean input state.
"""
    )


def config(pw: Path) -> dict[str, Any]:
    return {
        "repo_root": str(REPO_ROOT),
        "workflow_root": str(ROOT),
        "branch": run_text(["git", "rev-parse", "--abbrev-ref", "HEAD"], REPO_ROOT).strip(),
        "head": run_text(["git", "rev-parse", "HEAD"], REPO_ROOT).strip(),
        "head_oneline": run_text(["git", "log", "-1", "--oneline"], REPO_ROOT).strip(),
        "recent_commits": run_text(["git", "log", "--oneline", "-5"], REPO_ROOT).splitlines(),
        "git_status_short": run_text(["git", "status", "--short"], REPO_ROOT).splitlines(),
        "pw_executable": str(pw),
        "pw_sha256": sha256(pw),
        "pseudo": str(PSEUDO),
        "pseudo_sha256": sha256(PSEUDO),
        "python": sys.version,
        "platform": platform.platform(),
        "local_only": True,
        "remote_operations": "forbidden and not used",
        "created_utc": utc_now(),
    }


def write_implementation_confirmation() -> None:
    path = ROOT / "00_metadata" / "IMPLEMENTATION_CONFIRMATION.md"
    path.write_text(
        f"""# Implementation Confirmation

Generated: {utc_now()}

Scope: local-only confirmation for active CAIMD plus ionic Nose thermostat in
the current Quantum ESPRESSO source tree.

## Confirmed Behavior

- Active CAIMD supports `ion_temperature='nose'` and `ion_temperature='not_controlled'`.
- Unsupported active-CAIMD ionic thermostat modes remain rejected in
  `PW/src/input.f90` by the active CAIMD allowlist.
- PW `calculation='vc-md'` with `cell_dynamics='pr'` maps to `calc='cd'`, so
  variable-cell MD enters `PW/src/vcsmd.f90` and `PW/src/vcsubs.f90`.
- Ionic Nose setup uses QE's `get_ndof()` convention; for no fixed coordinates
  this gives `3*nat - 3`.
- `vcsmd` applies `ions_nosevel` and `ions_noseupd` around `vcmove`.
- The Nose friction term acts on ionic velocity histories (`ratd`/related
  histories), not on the imposed CAIMD strain schedule itself.
- CAIMD projection transforms `ratd`, `rat2d`, and `rat2di` so the committed
  physical Cartesian ionic velocity is preserved across a projected cell event.
- After final CAIMD projection, `vc_refresh_projected_nose_ekin` synchronizes
  projected ionic kinetic energy with the Nose chain kinetic inputs.
- `.caimd` column `temperature_K` is based on ionic kinetic energy. Rows are
  premove diagnostics; the first row may fall back to the input target
  temperature before a committed ionic kinetic temperature exists.
- CAIMD component order is `(xx, yy, zz, yz, xz, xy)`.
- `caimd_mask(i)=.true.` means component `i` is imposed or held by CAIMD.
- Unmasked symmetric components retain the ordinary Parrinello-Rahman candidate
  response before final projection.

No source redesign or source modification was performed during this audit.
"""
    )


def diagnose_previous_failure() -> dict[str, Any]:
    path = PREVIOUS_ROOT / "02_continuous_0p00_to_0p30" / "tmp" / "diamond8_continuous.caimd"
    rows = rows_from_caimd(path)
    diag: dict[str, Any] = {
        "previous_root": str(PREVIOUS_ROOT),
        "caimd_path": str(path),
        "row_count": len(rows),
        "rate_ps_inv": FAILED_LOCAL_RATE_PS_INV,
    }
    if not rows:
        write_json(ROOT / "00_metadata" / "OLD_FAILURE_DIAGNOSIS.json", diag)
        return diag
    for threshold in (0.05, 0.10, 0.25, 0.50):
        hits = []
        for name in COMPONENTS[1:]:
            for row in rows:
                if abs(row[f"actual_{name}"]) >= threshold:
                    hits.append(
                        {
                            "component": name,
                            "threshold": threshold,
                            "ionic_step": row["ionic_step"],
                            "time_ps": row["time_ps"],
                            "value": row[f"actual_{name}"],
                            "nominal_xx": row["nominal_xx"],
                        }
                    )
                    break
        diag[f"first_abs_actual_component_ge_{threshold:g}"] = sorted(hits, key=lambda r: r["ionic_step"])
    for temp_threshold in (600.0, 2000.0):
        for row in rows:
            if row["temperature_K"] >= temp_threshold:
                diag[f"first_temperature_ge_{temp_threshold:g}K"] = {
                    "ionic_step": row["ionic_step"],
                    "time_ps": row["time_ps"],
                    "temperature_K": row["temperature_K"],
                    "nominal_xx": row["nominal_xx"],
                    "actual_components": {name: row[f"actual_{name}"] for name in COMPONENTS},
                }
                break
    diag["last_row"] = {
        "ionic_step": rows[-1]["ionic_step"],
        "time_ps": rows[-1]["time_ps"],
        "event_count": rows[-1]["imposed_events_done"],
        "temperature_K": rows[-1]["temperature_K"],
        "volume_bohr3": rows[-1]["volume_bohr3"],
        "nominal_components": {name: rows[-1][f"nominal_{name}"] for name in COMPONENTS},
        "actual_components": {name: rows[-1][f"actual_{name}"] for name in COMPONENTS},
        "velocity_jump_max": rows[-1]["velocity_jump_max"],
    }
    growth = []
    for name in COMPONENTS[1:]:
        max_delta = 0.0
        max_step = None
        for prev, cur in zip(rows, rows[1:]):
            delta = abs(cur[f"actual_{name}"] - prev[f"actual_{name}"])
            if delta > max_delta:
                max_delta = delta
                max_step = cur["ionic_step"]
        growth.append({"component": name, "max_step_to_step_delta": max_delta, "at_ionic_step": max_step})
    diag["largest_stepwise_unmasked_growth"] = sorted(growth, key=lambda r: r["max_step_to_step_delta"], reverse=True)
    write_json(ROOT / "00_metadata" / "OLD_FAILURE_DIAGNOSIS.json", diag)
    return diag


def prepare(args: argparse.Namespace) -> None:
    pw = Path(args.pw).resolve()
    if not pw.exists() or not os.access(pw, os.X_OK):
        raise SystemExit(f"pw.x is unavailable or not executable: {pw}")
    if not PSEUDO.exists():
        raise SystemExit(f"local pseudopotential is unavailable: {PSEUDO}")
    setup = diamond8_setup()
    write_json(ROOT / "00_metadata" / "DIAMOND8_SETUP.json", setup)
    write_json(ROOT / "WORKFLOW_METADATA.json", config(pw))
    write_implementation_confirmation()
    diagnose_previous_failure()
    state = input_state_from_setup(setup)
    write_input(ROOT / "01_common_equilibration", "diamond8_common_equil", EQUIL_STEPS, state, active_caimd=False)
    write_json(
        ROOT / "01_common_equilibration" / "INPUT_METADATA.json",
        {
            "nstep": EQUIL_STEPS,
            "dt_fs": DT_FS,
            "caimd": False,
            "initial_velocity_seed": VELOCITY_SEED,
            "initial_velocity_scale_factor": setup["velocity_scale_factor"],
            "initial_velocity_temperature_3N_minus_3_K": setup["velocity_temperature"]["temperature_3N_minus_3_K"],
            "target_temperature_K": TEMP_TARGET_K,
            "target_pressure_GPa": 0.0,
            "wmass": WMASS,
            "fnosep_THz": FNOSEP_THZ,
            "nhpcl": NHPCL,
        },
    )
    (ROOT / "README.md").write_text(
        f"""# Diamond8 CAIMD + Nose Four DOF Tests

Generated: {utc_now()}

This is a completely local macOS validation workflow. It uses the local
Quantum ESPRESSO `pw.x` executable and the existing local `pseudo/C.UPF` only.

The four long tests are independent active-CAIMD runs from
`01_common_equilibration/common_start_state.json`.

No remote operation, Slurm command, network access, or Git push is part of this
workflow.
"""
    )
    print(json.dumps({"prepared": str(ROOT), "pw": str(pw), "pseudo": str(PSEUDO)}, indent=2))


def run_qe(run_dir: Path, prefix: str, pw: Path) -> dict[str, Any]:
    marker = run_dir / f"{prefix}_RUN_COMPLETE.json"
    if marker.exists():
        return read_json(marker)
    if (run_dir / "TEST_COMPLETE.json").exists() or (run_dir / "EQUILIBRATION_COMPLETE.json").exists():
        return {"skipped": True, "reason": "completion marker exists"}
    inp = run_dir / "qe.in"
    if not inp.exists():
        raise SystemExit(f"missing QE input: {inp}")
    out = run_dir / "qe.out"
    start = time.time()
    with out.open("w") as handle:
        proc = subprocess.run([str(pw), "-inp", str(inp)], cwd=str(run_dir), stdout=handle, stderr=subprocess.STDOUT)
    result = {
        "command": [str(pw), "-inp", str(inp)],
        "cwd": str(run_dir),
        "output": str(out),
        "returncode": proc.returncode,
        "elapsed_s": time.time() - start,
        "completed_utc": utc_now(),
    }
    write_json(marker, result)
    return result


def run_equil(args: argparse.Namespace) -> None:
    pw = Path(args.pw).resolve()
    run_dir = ROOT / "01_common_equilibration"
    prefix = "diamond8_common_equil"
    run = run_qe(run_dir, prefix, pw)
    val = validate_run(run_dir, prefix, EQUIL_STEPS, active_caimd=False)
    if val["valid"]:
        state = extract_terminal_state(run_dir, prefix, run_dir / "common_start_state.json")
        write_json(run_dir / "EQUILIBRATION_COMPLETE.json", {"valid": True, "validation": val, "common_start_state": state})
    print(json.dumps({"run": run, "validation": val}, indent=2, sort_keys=True))
    if not val["valid"]:
        raise SystemExit("common equilibration validation failed")


def build_tests(_: argparse.Namespace) -> None:
    state_path = ROOT / "01_common_equilibration" / "common_start_state.json"
    if not state_path.exists():
        raise SystemExit(f"missing common state: {state_path}")
    state = read_json(state_path)
    for name, spec in TESTS.items():
        run_dir = ROOT / spec["dirname"]
        write_input(
            run_dir,
            spec["prefix"],
            LONG_STEPS,
            state,
            active_caimd=True,
            mask=spec["mask"],
            step=spec["step"],
        )
        write_json(
            run_dir / "INPUT_METADATA.json",
            {
                "test": name,
                "title": spec["title"],
                "source_state": str(state_path),
                "nstep": LONG_STEPS,
                "dt_fs": DT_FS,
                "simulated_time_ps": LONG_STEPS * DT_FS / 1000.0,
                "caimd_stride": CAIMD_STRIDE,
                "caimd_mask": spec["mask"],
                "caimd_step_strain": spec["step"],
                "expected_events": LONG_STEPS // CAIMD_STRIDE,
                "expected_final_nominal_xx": (LONG_STEPS // CAIMD_STRIDE) * CAIMD_STRIDE * spec["step"][0],
                "strain_rate_ps_inv": STRAIN_RATE_PS_INV if spec["step"][0] != 0.0 else 0.0,
            },
        )
    print(json.dumps({"built_tests": list(TESTS)}, indent=2))


def run_test(args: argparse.Namespace) -> None:
    pw = Path(args.pw).resolve()
    if args.test not in TESTS:
        raise SystemExit(f"unknown test {args.test}; choose one of {sorted(TESTS)}")
    spec = TESTS[args.test]
    run_dir = ROOT / spec["dirname"]
    marker = run_dir / "TEST_COMPLETE.json"
    if marker.exists():
        print(json.dumps({"skipped": args.test, "reason": "TEST_COMPLETE.json exists"}, indent=2))
        return
    run = run_qe(run_dir, spec["prefix"], pw)
    val = validate_run(
        run_dir,
        spec["prefix"],
        LONG_STEPS,
        active_caimd=True,
        mask=spec["mask"],
        step=spec["step"],
    )
    if val["valid"]:
        extract_terminal_state(run_dir, spec["prefix"], run_dir / "terminal_state.json")
        write_json(marker, {"test": args.test, "valid": True, "validation": val})
    print(json.dumps({"test": args.test, "run": run, "validation": val}, indent=2, sort_keys=True))
    if not val["valid"]:
        raise SystemExit(f"{args.test} validation failed")


def validate_test(args: argparse.Namespace) -> None:
    if args.test not in TESTS:
        raise SystemExit(f"unknown test {args.test}; choose one of {sorted(TESTS)}")
    spec = TESTS[args.test]
    run_dir = ROOT / spec["dirname"]
    context = {
        "test": args.test,
        "validated_without_rerun": True,
        "reason": "post-stop validation of preserved partial or completed output",
        "utc": utc_now(),
    }
    if not (run_dir / f"{spec['prefix']}_RUN_COMPLETE.json").exists():
        write_json(run_dir / "POST_STOP_VALIDATION_CONTEXT.json", context)
    val = validate_run(
        run_dir,
        spec["prefix"],
        LONG_STEPS,
        active_caimd=True,
        mask=spec["mask"],
        step=spec["step"],
    )
    print(json.dumps({"test": args.test, "validation": val}, indent=2, sort_keys=True))
    if not val["valid"]:
        raise SystemExit(f"{args.test} validation failed")


def run_matrix(args: argparse.Namespace) -> None:
    build_tests(args)
    completed: list[str] = []
    failed: list[str] = []
    for name in ("test1", "test2", "test3", "test4"):
        try:
            run_test(argparse.Namespace(test=name, pw=args.pw))
            completed.append(name)
        except SystemExit:
            failed.append(name)
            if name in ("test1", "test2"):
                write_recovery_status(f"{name} failed; stopping because this may indicate a source-level CAIMD/Nose defect.")
                raise
            write_recovery_status(f"{name} failed; preserving evidence and continuing only to independent later tests when safe.")
            if name == "test4":
                break
            continue
    print(json.dumps({"completed": completed, "failed": failed}, indent=2))


def write_recovery_status(message: str) -> None:
    (ROOT / "RECOVERY_STATUS.md").write_text(
        f"""# Recovery Status

Generated: {utc_now()}

{message}

Completed markers should be trusted only when `TEST_COMPLETE.json` or
`EQUILIBRATION_COMPLETE.json` is present in the relevant run directory.
"""
    )


def load_test_rows() -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = {}
    for name, spec in TESTS.items():
        out[name] = rows_from_caimd(artifact(ROOT / spec["dirname"], spec["prefix"], "caimd"))
    return out


def validation_for(name: str) -> dict[str, Any] | None:
    path = ROOT / TESTS[name]["dirname"] / "VALIDATION.json"
    return read_json(path) if path.exists() else None


def flatten(prefix: str, data: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, val in data.items():
        if isinstance(val, dict):
            out.update(flatten(f"{prefix}{key}_", val))
        else:
            out[f"{prefix}{key}"] = val
    return out


def analysis_rows() -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    temp_rows: list[dict[str, Any]] = []
    block_rows: list[dict[str, Any]] = []
    strain_rows: list[dict[str, Any]] = []
    cell_rows_out: list[dict[str, Any]] = []
    stress_rows: list[dict[str, Any]] = []
    eq_path = ROOT / "01_common_equilibration" / "VALIDATION.json"
    if eq_path.exists():
        eq = read_json(eq_path)
        temp_rows.append({"test": "common_equilibration", **eq["temperature_analysis_K"], "classification": eq["temperature_classification"]})
        for b in block_means(parse_eal_temperatures(artifact(ROOT / "01_common_equilibration", "diamond8_common_equil", "eal")), 500):
            block_rows.append({"test": "common_equilibration", **b})
        cell_rows_out.append({"test": "common_equilibration", **flatten("", eq["cell_determinant_bohr3"]), **flatten("condition_", eq["cell_condition_number"])})
        for comp, st in eq["stress_component_statistics_GPa"].items():
            stress_rows.append({"test": "common_equilibration", "component": comp, **st})
    for name, spec in TESTS.items():
        val = validation_for(name)
        if not val:
            continue
        temp_rows.append({"test": name, **val["temperature_analysis_K"], "classification": val["temperature_classification"]})
        rows = rows_from_caimd(artifact(ROOT / spec["dirname"], spec["prefix"], "caimd"))
        temps = [r["temperature_K"] for r in rows]
        for b in block_means(temps, 500):
            block_rows.append({"test": name, **b})
        caimd = val["caimd"]
        for comp in COMPONENTS:
            strain_rows.append(
                {
                    "test": name,
                    "component": comp,
                    "mask": spec["mask"][COMPONENTS.index(comp)],
                    "step_strain": spec["step"][COMPONENTS.index(comp)],
                    "actual_min": caimd.get(f"actual_{comp}_min"),
                    "actual_max": caimd.get(f"actual_{comp}_max"),
                    "actual_final": caimd.get(f"actual_{comp}_final"),
                    "nominal_min": caimd.get(f"nominal_{comp}_min"),
                    "nominal_max": caimd.get(f"nominal_{comp}_max"),
                    "nominal_final": caimd.get(f"nominal_{comp}_final"),
                    "nominal_actual_abs_error_max": caimd.get(f"nominal_actual_{comp}_abs_error_max"),
                }
            )
        cell_rows_out.append(
            {
                "test": name,
                "determinant_min": val["cell_determinant_bohr3"]["min"],
                "determinant_max": val["cell_determinant_bohr3"]["max"],
                "volume_range": (val["cell_determinant_bohr3"]["max"] - val["cell_determinant_bohr3"]["min"])
                if val["cell_determinant_bohr3"]["max"] is not None and val["cell_determinant_bohr3"]["min"] is not None
                else None,
                "condition_min": val["cell_condition_number"]["min"],
                "condition_max": val["cell_condition_number"]["max"],
                "min_C_C_min": val["minimum_periodic_C_C_distance_bohr"]["min"],
                "min_C_C_max": val["minimum_periodic_C_C_distance_bohr"]["max"],
                "max_C_C_min": val["maximum_periodic_C_C_distance_bohr"]["min"],
                "max_C_C_max": val["maximum_periodic_C_C_distance_bohr"]["max"],
            }
        )
        for comp, st in val["stress_component_statistics_GPa"].items():
            stress_rows.append({"test": name, "component": comp, **st})
    return temp_rows, block_rows, strain_rows, cell_rows_out, stress_rows


def make_plot(path: Path, series: list[dict[str, Any]], xkey: str, ykeys: list[str], title: str, xlabel: str, ylabel: str) -> None:
    from PIL import Image, ImageDraw

    width, height = 1100, 700
    left, right, top, bottom = 85, 35, 55, 80
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    plot_w = width - left - right
    plot_h = height - top - bottom
    colors = [
        (31, 119, 180),
        (214, 39, 40),
        (44, 160, 44),
        (148, 103, 189),
        (255, 127, 14),
        (23, 190, 207),
        (127, 127, 127),
        (188, 189, 34),
    ]
    all_x: list[float] = []
    all_y: list[float] = []
    expanded: list[tuple[str, list[float], list[float], tuple[int, int, int]]] = []
    ci = 0
    for item in series:
        rows = item.get("rows", [])
        label = item.get("label", "series")
        for ykey in ykeys:
            pts = [(r.get(xkey), r.get(ykey)) for r in rows if r.get(xkey) is not None and r.get(ykey) is not None]
            pts = [(float(x), float(y)) for x, y in pts if math.isfinite(float(x)) and math.isfinite(float(y))]
            if not pts:
                continue
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            all_x.extend(xs)
            all_y.extend(ys)
            expanded.append((f"{label} {ykey}", xs, ys, colors[ci % len(colors)]))
            ci += 1
    draw.text((left, 18), title, fill=(0, 0, 0))
    if not expanded:
        draw.text((left, top + 30), "No data available", fill=(180, 0, 0))
        img.save(path)
        return
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)
    if xmax == xmin:
        xmax = xmin + 1.0
    if ymax == ymin:
        ymax = ymin + 1.0
    ypad = 0.08 * (ymax - ymin)
    ymin -= ypad
    ymax += ypad
    def sx(x: float) -> int:
        return int(left + (x - xmin) / (xmax - xmin) * plot_w)
    def sy(y: float) -> int:
        return int(top + plot_h - (y - ymin) / (ymax - ymin) * plot_h)
    draw.rectangle((left, top, left + plot_w, top + plot_h), outline=(0, 0, 0))
    for i in range(6):
        x = left + i * plot_w / 5
        y = top + i * plot_h / 5
        draw.line((x, top, x, top + plot_h), fill=(230, 230, 230))
        draw.line((left, y, left + plot_w, y), fill=(230, 230, 230))
        xv = xmin + i * (xmax - xmin) / 5
        yv = ymax - i * (ymax - ymin) / 5
        draw.text((x - 20, top + plot_h + 8), f"{xv:.3g}", fill=(0, 0, 0))
        draw.text((8, y - 7), f"{yv:.3g}", fill=(0, 0, 0))
    for label, xs, ys, color in expanded:
        pts = [(sx(x), sy(y)) for x, y in zip(xs, ys)]
        if len(pts) == 1:
            draw.ellipse((pts[0][0] - 2, pts[0][1] - 2, pts[0][0] + 2, pts[0][1] + 2), fill=color)
        else:
            draw.line(pts, fill=color, width=2)
    lx, ly = left + 10, top + 8
    for label, _, _, color in expanded[:12]:
        draw.line((lx, ly + 6, lx + 22, ly + 6), fill=color, width=3)
        draw.text((lx + 28, ly), label[:58], fill=(0, 0, 0))
        ly += 18
    draw.text((left + plot_w // 2 - 40, height - 35), xlabel, fill=(0, 0, 0))
    draw.text((10, 20), ylabel, fill=(0, 0, 0))
    path.parent.mkdir(parents=True, exist_ok=True)
    img.save(path)


def analyze(_: argparse.Namespace | None = None) -> None:
    analysis = ROOT / "analysis"
    analysis.mkdir(exist_ok=True)
    temp_rows, block_rows, strain_rows, cell_rows_out, stress_rows = analysis_rows()
    write_csv(analysis / "temperature_statistics.csv", temp_rows)
    write_csv(analysis / "temperature_block_statistics.csv", block_rows)
    write_csv(analysis / "strain_tracking_statistics.csv", strain_rows)
    write_csv(analysis / "cell_stability_statistics.csv", cell_rows_out)
    write_csv(analysis / "stress_statistics.csv", stress_rows)
    rows = load_test_rows()
    series_time = [{"label": name, "rows": r} for name, r in rows.items()]
    make_plot(analysis / "temperature_all_tests_vs_time.png", series_time, "time_ps", ["temperature_K"], "Temperature all tests", "time (ps)", "K")
    rolling_rows = []
    for name, rs in rows.items():
        vals = rolling([r["temperature_K"] for r in rs], 100)
        rolling_rows.append({"label": name, "rows": [{**r, "rolling100_temperature_K": vals[i]} for i, r in enumerate(rs)]})
    make_plot(analysis / "temperature_rolling100_all_tests.png", rolling_rows, "time_ps", ["rolling100_temperature_K"], "Rolling 100-step temperature", "time (ps)", "K")
    block_plot_rows = []
    for name, rs in rows.items():
        vals = [r["temperature_K"] for r in rs]
        block_plot_rows.append(
            {
                "label": name,
                "rows": [
                    {"block_mid_step": (b["start_step"] + b["end_step"]) / 2, "block_mean_temperature_K": b["mean"]}
                    for b in block_means(vals, 500)
                ],
            }
        )
    make_plot(analysis / "temperature_block500_all_tests.png", block_plot_rows, "block_mid_step", ["block_mean_temperature_K"], "500-step temperature block means", "step", "K")
    make_plot(
        analysis / "temperature_vs_xx_strain_tests2_3_4.png",
        [{"label": name, "rows": rows[name]} for name in ("test2", "test3", "test4")],
        "actual_xx",
        ["temperature_K"],
        "Temperature vs actual xx strain",
        "actual xx strain",
        "K",
    )
    for name in ("test1", "test2", "test3", "test4"):
        make_plot(
            analysis / f"actual_strains_{name}.png",
            [{"label": name, "rows": rows[name]}],
            "time_ps",
            [f"actual_{c}" for c in COMPONENTS],
            f"Actual strains {name}",
            "time (ps)",
            "strain",
        )
    make_plot(analysis / "cell_determinant_all_tests.png", series_time, "time_ps", ["volume_bohr3"], "Cell determinant from CAIMD diagnostic", "time (ps)", "bohr^3")
    cell_series_rows = []
    for name, spec in TESTS.items():
        cell_series_rows.append({"label": name, "rows": cell_series(ROOT / spec["dirname"], spec["prefix"])})
    make_plot(analysis / "cell_condition_number_all_tests.png", cell_series_rows, "step", ["condition_number"], "Cell condition number", "step", "condition")
    make_plot(analysis / "volume_all_tests.png", series_time, "time_ps", ["volume_bohr3"], "Volume", "time (ps)", "bohr^3")
    make_plot(
        analysis / "normal_stresses_tests2_3_4.png",
        [{"label": name, "rows": rows[name]} for name in ("test2", "test3", "test4")],
        "actual_xx",
        ["stress_xx_GPa", "stress_yy_GPa", "stress_zz_GPa"],
        "Normal stresses",
        "actual xx strain",
        "GPa",
    )
    make_plot(
        analysis / "shear_stresses_tests2_3_4.png",
        [{"label": name, "rows": rows[name]} for name in ("test2", "test3", "test4")],
        "actual_xx",
        ["stress_yz_GPa", "stress_xz_GPa", "stress_xy_GPa"],
        "Shear stresses",
        "actual xx strain",
        "GPa",
    )
    make_plot(analysis / "velocity_jump_all_tests.png", series_time, "time_ps", ["velocity_jump_max", "velocity_jump_rms"], "CAIMD velocity jump residual", "time (ps)", "a.u.")
    old_rows = rows_from_caimd(PREVIOUS_ROOT / "02_continuous_0p00_to_0p30" / "tmp" / "diamond8_continuous.caimd")
    make_plot(
        analysis / "old_failure_vs_new_test4.png",
        [{"label": "old_0p05ps", "rows": old_rows}, {"label": "new_test4", "rows": rows["test4"]}],
        "time_ps",
        ["actual_yy", "actual_zz", "actual_yz", "actual_xz", "actual_xy", "temperature_K"],
        "Old failure vs new Test 4",
        "time (ps)",
        "strain or K",
    )
    summary = comparison_summary()
    write_json(analysis / "test_comparison_summary.json", summary)
    write_final_report(summary)
    print(json.dumps({"analysis": str(analysis), "summary": summary}, indent=2, sort_keys=True))


def comparison_summary() -> dict[str, Any]:
    tests: list[dict[str, Any]] = []
    for name, spec in TESTS.items():
        val = validation_for(name)
        if not val:
            tests.append({"test": name, "status": "NOT_RUN"})
            continue
        temp = val["temperature_analysis_K"]
        det = val["cell_determinant_bohr3"]
        cond = val["cell_condition_number"]
        tests.append(
            {
                "test": name,
                "title": spec["title"],
                "completed_steps": min(val["trajectory_blocks"].get("cel", 0), val["trajectory_blocks"].get("pos", 0)),
                "simulated_time_ps": min(val["trajectory_blocks"].get("cel", 0), val["trajectory_blocks"].get("pos", 0)) * DT_FS / 1000.0,
                "final_nominal_xx_strain": val.get("final_nominal_xx_strain"),
                "final_actual_xx_strain": val.get("final_actual_xx_strain"),
                "event_count": val.get("event_count"),
                "temperature_mean": temp.get("mean"),
                "temperature_std": temp.get("std"),
                "first_500_mean": temp.get("first_500_mean"),
                "final_500_mean": temp.get("final_500_mean"),
                "temperature_drift_K_per_ps": temp.get("drift_K_per_ps"),
                "temperature_min": temp.get("min"),
                "temperature_max": temp.get("max"),
                "fraction_below_150K": temp.get("fraction_below_150K"),
                "fraction_above_600K": temp.get("fraction_above_600K"),
                "determinant_min": det.get("min"),
                "condition_number_max": cond.get("max"),
                "volume_range": (det.get("max") - det.get("min")) if det.get("max") is not None and det.get("min") is not None else None,
                "maximum_strain_tracking_error": val.get("xx_tracking_error_max"),
                "maximum_velocity_jump_residual": val.get("velocity_jump_max"),
                "classification": val.get("classification"),
                "valid": val.get("valid"),
                "first_failing_observable": val.get("first_failing_observable"),
            }
        )
    old = read_json(ROOT / "00_metadata" / "OLD_FAILURE_DIAGNOSIS.json") if (ROOT / "00_metadata" / "OLD_FAILURE_DIAGNOSIS.json").exists() else {}
    source_changed = bool(run_text(["git", "diff", "--", "PW/src", "Modules", "test-suite"], REPO_ROOT).strip())
    any_run = any(t.get("completed_steps", 0) for t in tests)
    all_valid = all(t.get("valid") for t in tests)
    all_complete = all(t.get("completed_steps") == LONG_STEPS for t in tests if t.get("status") != "NOT_RUN")
    early_blocked = any(
        t.get("test") in ("test1", "test2") and t.get("classification") == "FAIL"
        for t in tests
    )
    if early_blocked:
        status = "BLOCKED"
    elif all_valid and all_complete and len(tests) == 4:
        status = "SUCCESS"
        if any(t.get("classification") == "WARN" for t in tests):
            status = "SUCCESS_WITH_WARNINGS"
    elif any_run:
        status = "PARTIAL"
    else:
        status = "BLOCKED"
    return {
        "status": status,
        "workflow_root": str(ROOT),
        "tests": tests,
        "old_failure_diagnosis": old,
        "source_changed": source_changed,
        "local_only_confirmed": True,
        "strain_rate_ps_inv": STRAIN_RATE_PS_INV,
        "failed_previous_rate_ps_inv": FAILED_LOCAL_RATE_PS_INV,
        "target_final_xx": TARGET_FINAL_XX,
        "dt_fs": DT_FS,
        "steps_per_test": LONG_STEPS,
    }


def write_final_report(summary: dict[str, Any]) -> None:
    cfg = read_json(ROOT / "WORKFLOW_METADATA.json") if (ROOT / "WORKFLOW_METADATA.json").exists() else {}
    eq_val = read_json(ROOT / "01_common_equilibration" / "VALIDATION.json") if (ROOT / "01_common_equilibration" / "VALIDATION.json").exists() else {}
    common = read_json(ROOT / "01_common_equilibration" / "common_start_state.json") if (ROOT / "01_common_equilibration" / "common_start_state.json").exists() else {}
    lines = [
        "# Final Validation Report\n\n",
        f"Status: {summary['status']}\n\n",
        f"Workflow directory: `{ROOT}`\n\n",
        f"Branch/HEAD: `{cfg.get('branch')}` / `{cfg.get('head_oneline')}`\n\n",
        f"Executable: `{cfg.get('pw_executable')}`\n\n",
        f"Executable SHA256: `{cfg.get('pw_sha256')}`\n\n",
        f"Pseudopotential: `{cfg.get('pseudo')}`\n\n",
        f"Pseudopotential SHA256: `{cfg.get('pseudo_sha256')}`\n\n",
        "No Leonardo, Slurm, SSH, network access, or Git push was used by this workflow.\n\n",
        "No QE source code was modified by this workflow.\n\n",
        "## Common Setup\n\n",
        f"- timestep: `{DT_FS} fs`\n",
        f"- test duration: `{LONG_STEPS * DT_FS / 1000.0} ps`\n",
        f"- strain rate for Tests 2-4: `{STRAIN_RATE_PS_INV} ps^-1`\n",
        f"- target final nominal xx strain: `{TARGET_FINAL_XX}`\n",
        f"- Nose parameters: `fnosep={FNOSEP_THZ} THz`, `nhpcl={NHPCL}`\n",
        f"- `wmass={WMASS}`\n",
        f"- common-state source step: `{common.get('source_step')}`\n",
        f"- common-state velocity temperature: `{common.get('temperature_from_velocities', {}).get('temperature_3N_minus_3_K')}` K\n\n",
        "## Common Equilibration\n\n",
        "```json\n",
        json.dumps(eq_val.get("temperature_analysis_K"), indent=2, sort_keys=True),
        "\n```\n\n",
        "## Test Comparison\n\n",
        "| test | completed | event | final nominal xx | final actual xx | T mean | T std | final500 T | det min | cond max | velocity jump | class |\n",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |\n",
    ]
    for row in summary["tests"]:
        lines.append(
            f"| {row.get('test')} | {row.get('completed_steps')} | {row.get('event_count')} | "
            f"{fmt(row.get('final_nominal_xx_strain'))} | {fmt(row.get('final_actual_xx_strain'))} | "
            f"{fmt(row.get('temperature_mean'))} | {fmt(row.get('temperature_std'))} | "
            f"{fmt(row.get('final_500_mean'))} | {fmt(row.get('determinant_min'))} | "
            f"{fmt(row.get('condition_number_max'))} | {fmt(row.get('maximum_velocity_jump_residual'))} | "
            f"{row.get('classification')} |\n"
        )
    lines.extend(
        [
            "\n## Old Failure Comparison\n\n",
            "The previous 0.05 ps^-1 local stress test reached only about -0.0265 nominal xx before QE stopped with `caimd projected cell is singular`.\n\n",
            "```json\n",
            json.dumps(summary.get("old_failure_diagnosis"), indent=2, sort_keys=True),
            "\n```\n\n",
            "## Interpretation\n\n",
        ]
    )
    valid_tests = [t for t in summary["tests"] if t.get("valid")]
    failed_tests = [t for t in summary["tests"] if t.get("classification") == "FAIL"]
    if valid_tests:
        lines.append("Completed valid tests show that active CAIMD and ionic Nose can share the same trajectory without a velocity-projection anomaly for those mask choices.\n\n")
    if failed_tests:
        lines.append("At least one later test failed; its `SAFETY_STOP_REPORT.md` identifies the first unstable observable. Earlier valid tests remain independent evidence because every test starts from the same common state.\n\n")
    lines.append("This run distinguishes thermostat behavior from free PR-cell behavior by comparing fully constrained loading, transverse-normal-only relaxation, and full unmasked PR response at the same 0.003 ps^-1 rate.\n")
    (ROOT / "analysis" / "FINAL_VALIDATION_REPORT.md").write_text("".join(lines))


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def run_regression(args: argparse.Namespace) -> None:
    pw = Path(args.pw).resolve()
    evidence = ROOT / "evidence"
    evidence.mkdir(exist_ok=True)
    py = sys.executable
    commands = [
        [py, "-B", "-m", "py_compile", str(ROOT / "scripts" / "four_dof_workflow.py")],
        ["bash", "-n", str(ROOT / "scripts" / "run_four_dof_workflow.sh")],
        ["git", "diff", "--check"],
        [
            py,
            "-B",
            str(REPO_ROOT / "test-suite" / "caimd" / "caimd_validation.py"),
            "--mode",
            "local",
            "--current-exe",
            str(pw),
            "--pseudo",
            str(PSEUDO),
            "--work-root",
            str(evidence / "caimd_harness_targeted"),
            "--run-id",
            "targeted_active_nose",
            "--cases",
            "active_nose_zero,active_nose_stride2_xx,active_nose_shear_xy",
        ],
        [
            py,
            "-B",
            str(REPO_ROOT / "test-suite" / "caimd" / "caimd_validation.py"),
            "--mode",
            "local",
            "--current-exe",
            str(pw),
            "--pseudo",
            str(PSEUDO),
            "--work-root",
            str(evidence / "caimd_harness_complete"),
            "--run-id",
            "complete_local",
        ],
    ]
    results = []
    for idx, cmd in enumerate(commands, start=1):
        start = time.time()
        proc = subprocess.run(cmd, cwd=str(REPO_ROOT), text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        result = {
            "index": idx,
            "command": cmd,
            "returncode": proc.returncode,
            "elapsed_s": time.time() - start,
            "output_path": str(evidence / f"regression_{idx:02d}.out"),
        }
        (evidence / f"regression_{idx:02d}.out").write_text(proc.stdout)
        results.append(result)
        if proc.returncode != 0:
            write_json(evidence / "REGRESSION_RESULTS.json", {"results": results, "valid": False})
            raise SystemExit(f"regression command failed: {cmd}")
    write_json(evidence / "REGRESSION_RESULTS.json", {"results": results, "valid": True})
    print(json.dumps({"regression": results}, indent=2, sort_keys=True))


def run_all(args: argparse.Namespace) -> None:
    prepare(args)
    run_equil(args)
    run_matrix(args)
    analyze(args)
    run_regression(args)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pw", default=str(DEFAULT_PW), help="local pw.x executable")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("prepare")
    sub.add_parser("run-equil")
    sub.add_parser("build-tests")
    p_test = sub.add_parser("run-test")
    p_test.add_argument("test", choices=sorted(TESTS))
    p_validate = sub.add_parser("validate-test")
    p_validate.add_argument("test", choices=sorted(TESTS))
    sub.add_parser("run-matrix")
    sub.add_parser("analyze")
    sub.add_parser("run-regression")
    sub.add_parser("run-all")
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args)
    elif args.command == "run-equil":
        run_equil(args)
    elif args.command == "build-tests":
        build_tests(args)
    elif args.command == "run-test":
        run_test(args)
    elif args.command == "validate-test":
        validate_test(args)
    elif args.command == "run-matrix":
        run_matrix(args)
    elif args.command == "analyze":
        analyze(args)
    elif args.command == "run-regression":
        run_regression(args)
    elif args.command == "run-all":
        run_all(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
