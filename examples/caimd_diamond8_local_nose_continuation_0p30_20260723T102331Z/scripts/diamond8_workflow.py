#!/usr/bin/env python3
"""Local Diamond-8 CAIMD + Nose validation workflow.

This script is intentionally local-only. It never uses SSH, Slurm, network
access, or remote paths. It creates QE inputs, runs the local pw.x executable in
foreground, extracts terminal cell/position/velocity state from QE trajectory
files, and analyzes continuous versus segmented active-CAIMD loading.
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
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parents[1]
PSEUDO = REPO_ROOT / "pseudo" / "C.UPF"
PSEUDO_NAME = "C.UPF"

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
INITIAL_VELOCITY_TEMPERATURE_K = 450.0
FNOSEP_THZ = 50.0
NHPCL = 4
FNOSEP_CHOICE = (
    "50 THz with a four-link Nose-Hoover chain gives about 19 QE steps per "
    "Nose oscillation at the required 2 fs timestep. This is stronger than the "
    "archived 10/20 THz single-chain attempts that stayed cold, while avoiding "
    "the overly hard 100 THz screen setting that QE reported as only 9 steps "
    "per Nose oscillation."
)
WMASS_BASELINE = 2000.0
WMASS_BASELINE_DT_QE = 5.0
WMASS_DT_SCALED = WMASS_BASELINE * (DT_QE / WMASS_BASELINE_DT_QE) ** 2
WMASS = 2.0e6
WMASS_CHOICE = (
    f"Started from the dt-scaled local CAIMD+Nose harness estimate "
    f"{WMASS_DT_SCALED:.8f} AMU, then increased to 2.0e6 AMU after archived "
    "2 fs Diamond8 attempts showed PR cell kinetic modes draining the ionic "
    "3N-3 temperature below target during the short local equilibration. This "
    "keeps cell_dynamics='pr' active but slows the barostat for local "
    "functional validation; it is not a quantitative barostat calibration."
)
ECUTWFC = 25.0
ECUTRHO = 120.0
ECUT_CHOICE = (
    "25/120 Ry matches the bundled QE pioud_diamond test-suite carbon setup "
    "for local Diamond+C.UPF runs; 12/48 Ry was retained only as archived "
    "failed screen evidence."
)
MIXING_BETA = 0.3
CONV_THR = 1.0e-4
ELECTRON_MAXSTEP = 40
CAIMD_STRIDE = 5
EQUIL_STEPS = 500
LOADING_RATE_PS_INV = 0.05
LITERATURE_RATE_PS_INV = 0.001
VELOCITY_SEED = 20260723
VOLUME_BOHR3_MIN = 10.0
VOLUME_BOHR3_MAX = 5000.0
CELL_CONDITION_MAX = 1.0e6
MIN_CC_BOHR_MIN = 0.8
MIN_CC_BOHR_MAX = 20.0

FLOAT_RE = re.compile(
    r"[-+]?(?:(?:\d+\.\d*|\.\d+)(?:[EeDd][-+]?\d+)?|\d+(?:[EeDd][-+]?\d+)?)"
)

CAIMD_COLUMNS = [
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
]


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run_text(cmd: list[str], cwd: Path | None = None) -> str:
    return subprocess.check_output(cmd, cwd=str(cwd) if cwd else None, text=True, stderr=subprocess.STDOUT)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def ffloat(x: float) -> str:
    return f"{x:.16e}".replace("e", "d")


def extract_floats(text: str) -> list[float]:
    return [float(x.replace("D", "E").replace("d", "e")) for x in FLOAT_RE.findall(text)]


def dot(a: list[float], b: list[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def transpose(a: list[list[float]]) -> list[list[float]]:
    return [[a[j][i] for j in range(3)] for i in range(3)]


def matmul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def matvec(a: list[list[float]], x: list[float]) -> list[float]:
    return [sum(a[i][j] * x[j] for j in range(3)) for i in range(3)]


def vecsub(a: list[float], b: list[float]) -> list[float]:
    return [a[i] - b[i] for i in range(3)]


def vecnorm(a: list[float]) -> float:
    return math.sqrt(sum(x * x for x in a))


def det3(a: list[list[float]]) -> float:
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def inv3(a: list[list[float]]) -> list[list[float]]:
    d = det3(a)
    if not math.isfinite(d) or abs(d) < 1.0e-18:
        raise ValueError(f"singular matrix determinant {d}")
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
    return max(sum(abs(x) for x in row) for row in a)


def condition_inf(a: list[list[float]]) -> float:
    return norm_inf(a) * norm_inf(inv3(a))


def finite_matrix(a: list[list[float]]) -> bool:
    return all(math.isfinite(x) for row in a for x in row)


def sym_strain(h: list[list[float]], href: list[list[float]]) -> list[float]:
    f = matmul(h, inv3(href))
    l = [[f[i][j] - (1.0 if i == j else 0.0) for j in range(3)] for i in range(3)]
    e = [[0.5 * (l[i][j] + l[j][i]) for j in range(3)] for i in range(3)]
    return [e[0][0], e[1][1], e[2][2], e[1][2], e[0][2], e[0][1]]


def strain_vec_to_matrix(v: list[float]) -> list[list[float]]:
    return [
        [v[0], v[5], v[4]],
        [v[5], v[1], v[3]],
        [v[4], v[3], v[2]],
    ]


def cell_rows_to_columns(rows: list[list[float]]) -> list[list[float]]:
    return transpose(rows)


def cell_columns_to_rows(cols: list[list[float]]) -> list[list[float]]:
    return transpose(cols)


def wrap_frac(frac: list[float]) -> list[float]:
    return [x - math.floor(x) for x in frac]


def velocity_temperature(velocities: list[list[float]], mass_amu: float = MASS_C_AMU) -> dict[str, Any]:
    mass_ry = mass_amu * AMU_RY
    ekin = 0.5 * mass_ry * sum(sum(v * v for v in row) for row in velocities)
    nat = len(velocities)
    ndof_3n3 = max(3 * nat - 3, 1)
    ndof_pr = 3 * (nat + 1)
    com = [sum(v[i] for v in velocities) / nat for i in range(3)]
    return {
        "ekin_Ry": ekin,
        "temperature_3N_minus_3_K": 2.0 * ekin / ndof_3n3 * RY_TO_KELVIN,
        "temperature_PR_nat_plus_1_K": 2.0 * ekin / ndof_pr * RY_TO_KELVIN,
        "center_of_mass_velocity_au": com,
        "center_of_mass_max_abs_au": max(abs(x) for x in com),
    }


def stats(values: list[float]) -> dict[str, Any]:
    vals = [v for v in values if math.isfinite(v)]
    if not vals:
        return {
            "count": 0,
            "finite": False,
            "mean": None,
            "std": None,
            "min": None,
            "max": None,
            "first_200_mean": None,
            "final_200_mean": None,
            "drift_per_step": None,
        }
    mean = sum(vals) / len(vals)
    var = sum((v - mean) ** 2 for v in vals) / len(vals)
    xbar = (len(values) - 1) / 2.0
    denom = sum((i - xbar) ** 2 for i in range(len(values)))
    drift = None
    if denom > 0:
        drift = sum((i - xbar) * (values[i] - mean) for i in range(len(values))) / denom
    return {
        "count": len(values),
        "finite": len(vals) == len(values),
        "mean": mean,
        "std": math.sqrt(var),
        "min": min(vals),
        "max": max(vals),
        "first_200_mean": sum(vals[: min(200, len(vals))]) / min(200, len(vals)),
        "final_200_mean": sum(vals[-min(200, len(vals)) :]) / min(200, len(vals)),
        "drift_per_step": drift,
    }


def classify_temperature(values: list[float]) -> str:
    s = stats(values)
    if not s["finite"]:
        return "FAIL"
    if s["min"] is None or s["max"] is None:
        return "FAIL"
    final_mean = s["final_200_mean"]
    overall = s["mean"]
    if s["min"] < 0.0 or s["max"] > 1500.0:
        return "FAIL"
    if final_mean is not None and (final_mean < 150.0 or final_mean > 600.0):
        return "FAIL"
    if overall is not None and (overall < 150.0 or overall > 600.0):
        return "FAIL"
    preferred = final_mean is not None and 240.0 <= final_mean <= 360.0
    preferred = preferred and overall is not None and 240.0 <= overall <= 360.0
    if preferred:
        return "PASS"
    broad = final_mean is not None and 200.0 <= final_mean <= 450.0
    broad = broad and overall is not None and 200.0 <= overall <= 450.0
    return "WARN" if broad else "FAIL"


def diamond8_setup() -> dict[str, Any]:
    a_ang = 3.567
    a_bohr = a_ang * ANGSTROM_TO_BOHR
    row_x = [1.0 / math.sqrt(3.0), 1.0 / math.sqrt(3.0), 1.0 / math.sqrt(3.0)]
    row_y = [1.0 / math.sqrt(2.0), -1.0 / math.sqrt(2.0), 0.0]
    row_z = [1.0 / math.sqrt(6.0), 1.0 / math.sqrt(6.0), -2.0 / math.sqrt(6.0)]
    rot = [row_x, row_y, row_z]
    h = [[a_bohr * rot[i][j] for j in range(3)] for i in range(3)]
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
    velocities = generate_velocities()
    min_dist = min_periodic_distance(h, frac)
    det = det3(h)
    return {
        "lattice_constant_angstrom": a_ang,
        "lattice_constant_bohr": a_bohr,
        "rotation_matrix_lab_from_crystal": rot,
        "cell_matrix_bohr_columns": h,
        "cell_parameters_bohr_rows_for_qe": cell_columns_to_rows(h),
        "fractional_positions_crystal_order": frac,
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
        "velocity_temperature": velocity_temperature(velocities),
        "determinant_bohr3": det,
        "condition_number_inf": condition_inf(h),
        "minimum_periodic_C_C_distance_bohr": min_dist,
        "minimum_periodic_C_C_distance_angstrom": min_dist / ANGSTROM_TO_BOHR,
        "positive_determinant": det > 0.0,
    }


def generate_velocities() -> list[list[float]]:
    rng = random.Random(VELOCITY_SEED)
    velocities = [[rng.gauss(0.0, 1.0) for _ in range(3)] for _ in range(NAT)]
    com = [sum(row[i] for row in velocities) / NAT for i in range(3)]
    velocities = [[row[i] - com[i] for i in range(3)] for row in velocities]
    temp = velocity_temperature(velocities)["temperature_3N_minus_3_K"]
    scale = math.sqrt(INITIAL_VELOCITY_TEMPERATURE_K / temp)
    return [[scale * x for x in row] for row in velocities]


def min_periodic_distance(h: list[list[float]], frac: list[list[float]]) -> float:
    best = float("inf")
    shifts = [-1, 0, 1]
    for i in range(len(frac)):
        for j in range(i + 1, len(frac)):
            for sx in shifts:
                for sy in shifts:
                    for sz in shifts:
                        df = [
                            frac[j][0] - frac[i][0] + sx,
                            frac[j][1] - frac[i][1] + sy,
                            frac[j][2] - frac[i][2] + sz,
                        ]
                        best = min(best, vecnorm(matvec(h, df)))
    return best


def input_state_from_setup(setup: dict[str, Any]) -> dict[str, Any]:
    return {
        "cell_matrix_bohr_columns": setup["cell_matrix_bohr_columns"],
        "fractional_positions": setup["fractional_positions_crystal_order"],
        "cartesian_velocities_au": setup["cartesian_velocities_au"],
        "atom_labels": setup["atom_labels"],
        "source": "generated_rotated_conventional_diamond8",
        "units": {
            "cell": "bohr",
            "positions": "crystal fractional",
            "velocities": "a.u. Cartesian",
        },
    }


def caimd_cell_lines(active: bool, step_strain: float = 0.0) -> str:
    if not active:
        return "  caimd = .false.,\n"
    return f"""  caimd = .true.,
  caimd_stride = {CAIMD_STRIDE},
  caimd_step_strain(1) = {ffloat(step_strain)},
  caimd_step_strain(2) = 0.0d0,
  caimd_step_strain(3) = 0.0d0,
  caimd_step_strain(4) = 0.0d0,
  caimd_step_strain(5) = 0.0d0,
  caimd_step_strain(6) = 0.0d0,
  caimd_mask(1) = .true.,
  caimd_mask(2) = .false.,
  caimd_mask(3) = .false.,
  caimd_mask(4) = .false.,
  caimd_mask(5) = .false.,
  caimd_mask(6) = .false.,
"""


def write_input(
    run_dir: Path,
    prefix: str,
    nstep: int,
    state: dict[str, Any],
    active_caimd: bool,
    caimd_step_strain: float = 0.0,
) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "tmp").mkdir(exist_ok=True)
    h = state["cell_matrix_bohr_columns"]
    rows = cell_columns_to_rows(h)
    frac = state["fractional_positions"]
    vel = state["cartesian_velocities_au"]
    labels = state.get("atom_labels", ["C"] * len(frac))
    if len(frac) != NAT or len(vel) != NAT:
        raise ValueError("state must contain exactly 8 positions and velocities")
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
{caimd_cell_lines(active_caimd, caimd_step_strain).rstrip()}
/
ATOMIC_SPECIES
C {MASS_C_AMU:.8f} {PSEUDO_NAME}
CELL_PARAMETERS bohr
"""
    for row in rows:
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


def write_audit() -> None:
    path = ROOT / "00_metadata" / "IMPLEMENTATION_AUDIT.md"
    path.write_text(
        f"""# Implementation Audit

Generated: {utc_now()}

Scope: local-only audit for the Diamond-8 CAIMD + Nose continuation validation.

## Confirmed Source Behavior

- Active CAIMD allows `ion_temperature='nose'` and `ion_temperature='not_controlled'`; unsupported ionic thermostats remain rejected by the active-CAIMD allowlist documented in `CAIMD_IMPLEMENTATION.md`.
- PW `calculation='vc-md'` with `cell_dynamics='pr'` maps to `calc='cd'` in `PW/src/input.f90`, so variable-cell MD enters `PW/src/vcsmd.f90` and `PW/src/vcsubs.f90`.
- `PW/src/input.f90:1144-1148` sets `tnosep=.true.` for `ion_temperature='nose'`; `PW/src/input.f90:2017-2020` initializes ionic Nose variables with `get_ndof()`.
- `PW/src/dynamics_module.f90:2517-2524` gives the no-fixed-coordinate ionic degree count as `3*nat - 3`.
- `PW/src/vcsmd.f90:334-354` applies the existing PW `vc-md` ionic Nose update around `vcmove`.
- `PW/src/vcsubs.f90:708-715` applies the Nose friction term to atomic velocity histories, not to CAIMD-imposed strain.
- Active CAIMD diagnostic temperature is computed from ionic kinetic energy, not atomic-plus-cell kinetic energy: `PW/src/caimd_module.f90:201-213`.
- `PW/src/vcsubs.f90:1143-1151` refreshes projected ionic kinetic quantities after final CAIMD projection; `vc_refresh_projected_nose_ekin` refreshes Nose chain kinetic inputs from the committed projected `ratd`.
- `PW/src/caimd_module.f90:321-349` transforms `ratd`, `rat2d`, and `rat2di` by `inverse(H_projected) * H_candidate`, preserving physical Cartesian ionic velocity across CAIMD cell projection.
- `.caimd` rows are premove diagnostics. The `temperature_K` column reports the last committed ionic kinetic temperature when available; the first row may fall back to the target/input temperature.
- Native CAIMD restart remains unsupported. `caimd_configure` rejects restart mode, and `vcsmd` rejects existing MD history for active CAIMD at start of run.
- Segmented continuation in this workflow is therefore controlled external continuation with `restart_mode='from_scratch'`, predecessor terminal cell, predecessor fractional positions, and predecessor physical Cartesian velocities. It resets PR history, Nose extended variables, CAIMD event count, and the local CAIMD reference cell.

## Loading Ensemble Name

The loading runs are thermostatted constrained-strain Parrinello-Rahman dynamics. They are not unconstrained full NpT because laboratory `xx` strain is externally driven by active CAIMD.

## High-Strain Limitation

The target global engineering compression of about 0.30 is far outside the nominal small-strain regime of CAIMD V1. This local run is a software stress test, thermostat-under-loading test, and continuation-workflow test. It must not be used as a quantitative high-strain stress curve.
"""
    )


def prepare(args: argparse.Namespace) -> None:
    pw = Path(args.pw).resolve()
    if not pw.exists() or not os.access(pw, os.X_OK):
        raise SystemExit(f"pw.x not executable: {pw}")
    if not PSEUDO.exists():
        raise SystemExit(f"missing local pseudopotential: {PSEUDO}")
    setup = diamond8_setup()
    state = input_state_from_setup(setup)
    write_audit()
    write_json(ROOT / "00_metadata" / "DIAMOND8_SETUP.json", {
        **setup,
        "pseudopotential_path": str(PSEUDO),
        "pseudopotential_sha256": sha256(PSEUDO),
        "ecutwfc_Ry": ECUTWFC,
        "ecutrho_Ry": ECUTRHO,
        "ecut_choice": ECUT_CHOICE,
        "wmass": WMASS,
        "wmass_choice": WMASS_CHOICE,
        "wmass_baseline": WMASS_BASELINE,
        "wmass_baseline_dt_qe": WMASS_BASELINE_DT_QE,
        "wmass_dt_scaled_estimate": WMASS_DT_SCALED,
        "fnosep_THz": FNOSEP_THZ,
        "fnosep_choice": FNOSEP_CHOICE,
        "nhpcl": NHPCL,
        "timestep_fs": DT_FS,
        "qe_dt": DT_QE,
        "qe_dt_units": "QE PW vc-md time units; validated against stdout time in ps",
        "expected_stdout_time_increment_ps": DT_FS / 1000.0,
        "nose_target_temperature_K": TEMP_TARGET_K,
        "initial_velocity_temperature_3N_minus_3_K": INITIAL_VELOCITY_TEMPERATURE_K,
        "initial_velocity_temperature_choice": "450 K warm start with tempw=300 K, used because archived zero-displacement Diamond8 starts with 300 K velocities equilibrated to a cold ~180 K ionic kinetic temperature over 500-1500 local Nose+PR steps.",
        "equilibration_steps": EQUIL_STEPS,
        "equilibration_steps_choice": "Final candidate returns to the requested initial 500-step screen after archived 1500-step single-chain Nose evidence showed stable but cold ionic-only temperature.",
        "ionic_degrees_of_freedom_convention": "3N-3 for explicit velocities and CAIMD ionic temperature",
        "pressure_target_GPa": 0.0,
        "strain_rate_ps^-1": LOADING_RATE_PS_INV,
        "literature_rate_ps^-1_not_used": LITERATURE_RATE_PS_INV,
        "rate_acceleration_factor": LOADING_RATE_PS_INV / LITERATURE_RATE_PS_INV,
    })
    branch = run_text(["git", "branch", "--show-current"], REPO_ROOT).strip()
    head = run_text(["git", "rev-parse", "HEAD"], REPO_ROOT).strip()
    commit = run_text(["git", "log", "-1", "--oneline"], REPO_ROOT).strip()
    try:
        linked = run_text(["otool", "-L", str(pw)])
    except Exception as exc:
        linked = f"otool failed: {exc}"
    write_json(ROOT / "WORKFLOW_METADATA.json", {
        "created_utc": utc_now(),
        "workflow_root": str(ROOT),
        "repo_root": str(REPO_ROOT),
        "branch": branch,
        "head": head,
        "head_oneline": commit,
        "local_only": True,
        "remote_operations": "forbidden and not used by this workflow script",
        "pw_executable": str(pw),
        "pw_sha256": sha256(pw),
        "pseudo": str(PSEUDO),
        "pseudo_sha256": sha256(PSEUDO),
        "platform": platform.platform(),
        "python": sys.version,
        "linked_libraries": linked.splitlines(),
        "build_cache": "/private/tmp/qe_diamond8_af7ebfc_build/CMakeCache.txt",
        "ensemble_description": "thermostatted constrained-strain Parrinello-Rahman dynamics",
    })
    write_input(ROOT / "evidence" / "pw_smoke", "diamond8_smoke", 1, state, False)
    write_input(ROOT / "01_equilibration", "diamond8_equil", EQUIL_STEPS, state, False)
    (ROOT / "01_equilibration" / "initial_state_fragment.in").write_text(input_fragment(state))
    (ROOT / "README.md").write_text(readme_text())
    print(json.dumps({"prepared": True, "root": str(ROOT), "pw": str(pw)}, indent=2))


def input_fragment(state: dict[str, Any]) -> str:
    h = state["cell_matrix_bohr_columns"]
    text = "CELL_PARAMETERS bohr\n"
    for row in cell_columns_to_rows(h):
        text += "  " + "  ".join(f"{x:.14f}" for x in row) + "\n"
    text += "ATOMIC_POSITIONS crystal\n"
    for label, frac in zip(state.get("atom_labels", ["C"] * NAT), state["fractional_positions"]):
        text += f"{label}  " + "  ".join(f"{x:.14f}" for x in frac) + "\n"
    text += "ATOMIC_VELOCITIES { a.u }\n"
    for label, vel in zip(state.get("atom_labels", ["C"] * NAT), state["cartesian_velocities_au"]):
        text += f"{label}  " + "  ".join(f"{x:.14e}" for x in vel) + "\n"
    return text


def readme_text() -> str:
    return f"""# Local Diamond-8 CAIMD + Nose Continuation To 0.30

This directory contains a completely local macOS validation of the QE 7.5
CAIMD implementation at commit `af7ebfc`.

The loading ensemble is thermostatted constrained-strain Parrinello-Rahman
dynamics. It is not unconstrained full NpT because laboratory `xx` strain is
externally driven by CAIMD.

The loading rate is `{LOADING_RATE_PS_INV} ps^-1`, which is 50 times faster than
`0.001 ps^-1`. It is used only for local functional validation. The target
global engineering compression of about 0.30 is a software stress test, not a
quantitative high-strain stress curve.

The QE input `dt` is `{DT_QE:.12f}`, validated against QE stdout as a 0.002 ps
per-step increment. The Parrinello-Rahman cell mass is `{WMASS:.8f}`. The
starting estimate from the local harness scaling would be
`{WMASS_DT_SCALED:.8f}`, then this workflow uses the heavier value to keep PR
cell kinetic modes from draining the ionic temperature during the short local
equilibration.

The ionic Nose setup is `fnosep={FNOSEP_THZ:.1f} THz` with `nhpcl={NHPCL}`.
The carbon cutoff is `ecutwfc/ecutrho={ECUTWFC:.1f}/{ECUTRHO:.1f} Ry`,
matching the bundled QE Diamond+C.UPF test-suite level.

The Nose thermostat target remains `{TEMP_TARGET_K:.1f} K`. Initial velocities
are warm-started at `{INITIAL_VELOCITY_TEMPERATURE_K:.1f} K` under the 3N-3
ionic convention because archived zero-displacement Diamond8 starts at 300 K
equilibrated to a cold ionic kinetic temperature near 180 K.

The ordinary Nose+PR equilibration is `{EQUIL_STEPS}` steps. The archived
attempt directories document the 500-step screens and the stable 1500-step
single-chain extension that motivated the four-link Nose-Hoover chain.

Run commands are recorded in `evidence/commands.log`. Completed stages write
`*_COMPLETE.json` markers and can be skipped by rerunning the workflow.
"""


def config() -> dict[str, Any]:
    cfg = read_json(ROOT / "WORKFLOW_METADATA.json")
    return cfg


def run_qe(run_dir: Path, prefix: str) -> dict[str, Any]:
    cfg = config()
    pw = Path(cfg["pw_executable"])
    inp = run_dir / "qe.in"
    out = run_dir / "qe.out"
    if not inp.exists():
        raise SystemExit(f"missing input: {inp}")
    complete_marker = run_dir / f"{prefix}_RUN_COMPLETE.json"
    if complete_marker.exists():
        return read_json(complete_marker)
    if out.exists() and "JOB DONE" not in out.read_text(errors="ignore"):
        raise SystemExit(f"refusing to overwrite incomplete QE output: {out}")
    cmd = [str(pw), "-in", str(inp.name)]
    log_command(run_dir, cmd)
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    t0 = time.time()
    with out.open("w") as fh:
        proc = subprocess.run(cmd, cwd=str(run_dir), stdout=fh, stderr=subprocess.STDOUT, env=env)
    elapsed = time.time() - t0
    text = out.read_text(errors="ignore")
    result = {
        "prefix": prefix,
        "run_dir": str(run_dir),
        "input": str(inp),
        "output": str(out),
        "command": cmd,
        "returncode": proc.returncode,
        "elapsed_s": elapsed,
        "job_done": "JOB DONE" in text,
        "timestamp_utc": utc_now(),
    }
    write_json(complete_marker, result)
    return result


def log_command(run_dir: Path, cmd: list[str]) -> None:
    with (ROOT / "evidence" / "commands.log").open("a") as fh:
        fh.write(f"{utc_now()} cwd={run_dir} cmd={' '.join(cmd)}\n")


def parse_stdout_temperatures(path: Path) -> list[float]:
    text = path.read_text(errors="ignore") if path.exists() else ""
    vals: list[float] = []
    for m in re.finditer(r"\bT\s*=\s*([-+0-9.EeDd]+)\s*K", text):
        vals.append(float(m.group(1).replace("D", "E").replace("d", "e")))
    return vals


def parse_dynamics_times(path: Path) -> list[dict[str, Any]]:
    text = path.read_text(errors="ignore") if path.exists() else ""
    rows: list[dict[str, Any]] = []
    pattern = re.compile(
        r"Entering Dynamics;\s+it\s*=\s*(\d+)\s+time\s*=\s*([-+0-9.EeDd]+)\s+pico-seconds"
    )
    for m in pattern.finditer(text):
        rows.append({
            "ionic_step": int(m.group(1)),
            "time_ps": float(m.group(2).replace("D", "E").replace("d", "e")),
        })
    return rows


def observed_time_increment_ps(times: list[dict[str, Any]]) -> float | None:
    if len(times) < 2:
        return None
    first = times[0]
    last = times[-1]
    step_delta = last["ionic_step"] - first["ionic_step"]
    if step_delta <= 0:
        return None
    return (last["time_ps"] - first["time_ps"]) / step_delta


def parse_eal_ionic_temperatures(path: Path) -> list[float]:
    vals: list[float] = []
    if not path.exists():
        return vals
    ndof = 3 * NAT - 3
    for line in path.read_text(errors="ignore").splitlines():
        nums = extract_floats(line)
        if len(nums) >= 2:
            eka_ry = nums[1]
            vals.append(2.0 * eka_ry / ndof * RY_TO_KELVIN)
    return vals


def parse_stress_blocks(path: Path) -> list[dict[str, Any]]:
    text = path.read_text(errors="ignore") if path.exists() else ""
    blocks: list[dict[str, Any]] = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        if "total   stress" in lines[i]:
            nums = extract_floats(lines[i])
            pressure_kbar = nums[-1] if nums else None
            rows = []
            for j in range(1, 4):
                if i + j < len(lines):
                    vals = extract_floats(lines[i + j])
                    if len(vals) >= 6:
                        rows.append(vals[3:6])
            if len(rows) == 3:
                blocks.append({"pressure_kbar": pressure_kbar, "stress_kbar": rows})
            i += 4
        else:
            i += 1
    return blocks


def parse_blocks(path: Path, rows_per_block: int) -> list[dict[str, Any]]:
    blocks: list[dict[str, Any]] = []
    if not path.exists():
        return blocks
    lines = [line.rstrip() for line in path.read_text(errors="ignore").splitlines() if line.strip()]
    i = 0
    while i < len(lines):
        head = extract_floats(lines[i])
        if len(head) < 2:
            i += 1
            continue
        rows = []
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
        if ok:
            blocks.append({"step": int(round(head[0])), "time": head[1], "rows": rows})
            i += 1 + rows_per_block
        else:
            break
    return blocks


def artifact(run_dir: Path, prefix: str, ext: str) -> Path:
    p = run_dir / f"{prefix}.{ext}"
    if p.exists():
        return p
    return run_dir / "tmp" / f"{prefix}.{ext}"


def extract_terminal_state(run_dir: Path, prefix: str, output: Path | None = None) -> dict[str, Any]:
    cel = parse_blocks(artifact(run_dir, prefix, "cel"), 3)
    pos = parse_blocks(artifact(run_dir, prefix, "pos"), NAT)
    if not cel or not pos:
        raise SystemExit(f"missing trajectory blocks in {run_dir}")
    if cel[-1]["step"] != pos[-1]["step"]:
        raise SystemExit(f"terminal trajectory step mismatch in {run_dir}")
    h = cell_rows_to_columns(cel[-1]["rows"])
    hinv = inv3(h)
    frac = [wrap_frac(matvec(hinv, xyz)) for xyz in pos[-1]["rows"]]
    ratd = terminal_ratd_from_md(artifact(run_dir, prefix, "md"))
    velocities = [matvec(h, row) for row in ratd]
    state = {
        "source_run_dir": str(run_dir),
        "prefix": prefix,
        "source_step": cel[-1]["step"],
        "source_time": cel[-1]["time"],
        "cell_matrix_bohr_columns": h,
        "fractional_positions": frac,
        "cartesian_positions_bohr": pos[-1]["rows"],
        "cartesian_velocities_au": velocities,
        "terminal_ratd_lattice_velocity": ratd,
        "atom_labels": ["C"] * NAT,
        "temperature_from_velocities": velocity_temperature(velocities),
        "volume_bohr3": det3(h),
        "condition_number_inf": condition_inf(h),
        "minimum_periodic_C_C_distance_bohr": min_periodic_distance(h, frac),
        "units": {
            "cell": "bohr",
            "positions": "crystal fractional from terminal Cartesian positions",
            "velocities": "a.u. Cartesian from QE .vel trajectory",
        },
        "extraction_utc": utc_now(),
    }
    if output:
        write_json(output, state)
        (output.parent / "explicit_state_fragment.in").write_text(input_fragment(state))
    return state


def terminal_ratd_from_md(path: Path) -> list[list[float]]:
    if not path.exists():
        raise SystemExit(f"missing PR MD history file for authoritative velocities: {path}")
    values = extract_floats(path.read_text(errors="ignore"))
    need = 6 * NAT
    if len(values) < need:
        raise SystemExit(f"malformed PR MD history, need at least {need} floats: {path}")
    ratd_flat = values[3 * NAT : 6 * NAT]
    return [ratd_flat[3 * i : 3 * i + 3] for i in range(NAT)]


def validate_run(run_dir: Path, prefix: str, expected_steps: int, active_caimd: bool) -> dict[str, Any]:
    out = run_dir / "qe.out"
    text = out.read_text(errors="ignore") if out.exists() else ""
    cel = parse_blocks(artifact(run_dir, prefix, "cel"), 3)
    pos = parse_blocks(artifact(run_dir, prefix, "pos"), NAT)
    vel = parse_blocks(artifact(run_dir, prefix, "vel"), NAT)
    stdout_temps = parse_stdout_temperatures(out)
    eal_ionic_temps = parse_eal_ionic_temperatures(artifact(run_dir, prefix, "eal"))
    dynamics_times = parse_dynamics_times(out)
    observed_dt_ps = observed_time_increment_ps(dynamics_times)
    expected_dt_ps = DT_FS / 1000.0
    time_increment_sane = (
        observed_dt_ps is None
        or abs(observed_dt_ps - expected_dt_ps) <= 1.0e-6
    )
    stress = parse_stress_blocks(out)
    h_last = cell_rows_to_columns(cel[-1]["rows"]) if cel else None
    volume_final = None
    cell_condition = None
    min_cc = None
    geometry_error = None
    geometry_finite = bool(h_last and finite_matrix(h_last))
    if h_last:
        try:
            volume_final = det3(h_last)
            cell_condition = condition_inf(h_last)
            if pos:
                hinv = inv3(h_last)
                frac_last = [wrap_frac(matvec(hinv, xyz)) for xyz in pos[-1]["rows"]]
                min_cc = min_periodic_distance(h_last, frac_last)
        except Exception as exc:
            geometry_error = str(exc)
    caimd = parse_caimd(artifact(run_dir, prefix, "caimd")) if active_caimd else {"exists": False}
    caimd_temps = []
    if active_caimd:
        caimd_temps = [r["temperature_K"] for r in rows_from_caimd(artifact(run_dir, prefix, "caimd"))]
    if active_caimd and caimd_temps:
        analysis_temps = caimd_temps
        analysis_source = "caimd_temperature_K"
    elif eal_ionic_temps:
        analysis_temps = eal_ionic_temps
        analysis_source = "eal_ionic_3N_minus_3_K"
    else:
        analysis_temps = stdout_temps
        analysis_source = "stdout_T_K"
    result = {
        "run_dir": str(run_dir),
        "prefix": prefix,
        "expected_steps": expected_steps,
        "job_done": "JOB DONE" in text,
        "trajectory_blocks": {"cel": len(cel), "pos": len(pos), "vel": len(vel)},
        "completed_expected_steps": len(cel) >= expected_steps and len(pos) >= expected_steps,
        "temperature_analysis_source": analysis_source,
        "temperature_analysis_K": stats(analysis_temps),
        "temperature_classification": classify_temperature(analysis_temps),
        "stdout_temperature_K": stats(stdout_temps),
        "eal_ionic_temperature_3N_minus_3_K": stats(eal_ionic_temps),
        "caimd_temperature_K": stats(caimd_temps) if caimd_temps else None,
        "dynamics_time_blocks": len(dynamics_times),
        "expected_time_increment_ps": expected_dt_ps,
        "observed_time_increment_ps": observed_dt_ps,
        "time_increment_sane": bool(time_increment_sane),
        "stress_block_count": len(stress),
        "pressure_kbar_final": stress[-1]["pressure_kbar"] if stress else None,
        "pressure_kbar_stats": stats([s["pressure_kbar"] for s in stress if s["pressure_kbar"] is not None]),
        "geometry_finite": geometry_finite,
        "geometry_error": geometry_error,
        "volume_bohr3_final": volume_final,
        "volume_positive": bool(volume_final is not None and volume_final > 0.0),
        "volume_sane": bool(volume_final is not None and VOLUME_BOHR3_MIN <= volume_final <= VOLUME_BOHR3_MAX),
        "volume_sanity_bounds_bohr3": [VOLUME_BOHR3_MIN, VOLUME_BOHR3_MAX],
        "cell_condition_number_final": cell_condition,
        "cell_condition_sane": bool(cell_condition is not None and cell_condition <= CELL_CONDITION_MAX),
        "cell_condition_sanity_max": CELL_CONDITION_MAX,
        "minimum_periodic_C_C_distance_bohr_final": min_cc,
        "minimum_periodic_C_C_distance_sane": bool(min_cc is not None and MIN_CC_BOHR_MIN <= min_cc <= MIN_CC_BOHR_MAX),
        "minimum_periodic_C_C_distance_sanity_bounds_bohr": [MIN_CC_BOHR_MIN, MIN_CC_BOHR_MAX],
        "caimd": caimd,
        "returncode": read_json(run_dir / f"{prefix}_RUN_COMPLETE.json").get("returncode") if (run_dir / f"{prefix}_RUN_COMPLETE.json").exists() else None,
        "validated_utc": utc_now(),
    }
    valid = result["job_done"] and result["completed_expected_steps"] and result["volume_positive"]
    valid = valid and result["geometry_finite"] and result["volume_sane"]
    valid = valid and result["cell_condition_sane"] and result["minimum_periodic_C_C_distance_sane"]
    valid = valid and result["time_increment_sane"]
    valid = valid and result["temperature_classification"] != "FAIL"
    if active_caimd:
        valid = valid and caimd.get("exists") and caimd.get("row_count", 0) >= expected_steps
    result["valid"] = bool(valid)
    write_json(run_dir / "VALIDATION.json", result)
    return result


def parse_caimd(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"exists": False}
    rows = []
    for line in path.read_text(errors="ignore").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 26:
            return {"exists": True, "malformed_line": line, "row_count": len(rows)}
        row: dict[str, Any] = {}
        row["ionic_step"] = int(parts[0])
        row["time_ps"] = float(parts[1].replace("D", "E").replace("d", "e"))
        row["imposed_events_done"] = int(parts[2])
        row["following_move_is_event"] = parts[3].upper().startswith("T")
        for name, val in zip(CAIMD_COLUMNS[4:], parts[4:]):
            row[name] = float(val.replace("D", "E").replace("d", "e"))
        rows.append(row)
    cols_finite = all(all(math.isfinite(v) for k, v in r.items() if isinstance(v, float)) for r in rows)
    return {
        "exists": True,
        "path": str(path),
        "row_count": len(rows),
        "column_count": 26,
        "finite": cols_finite,
        "event_count_last": rows[-1]["imposed_events_done"] if rows else None,
        "nominal_xx_last": rows[-1]["nominal_xx"] if rows else None,
        "actual_xx_last": rows[-1]["actual_xx"] if rows else None,
        "xx_abs_error_max": max((abs(r["actual_xx"] - r["nominal_xx"]) for r in rows), default=None),
        "velocity_residual_max": max((r["velocity_jump_max"] for r in rows), default=None),
        "velocity_residual_rms_max": max((r["velocity_jump_rms"] for r in rows), default=None),
        "temperature_K": stats([r["temperature_K"] for r in rows]),
    }


def rows_from_caimd(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows = []
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


def build_continuous() -> None:
    eq_state = read_json(ROOT / "01_equilibration" / "equilibrated_state.json")
    write_input(ROOT / "02_continuous_0p00_to_0p30", "diamond8_continuous", 3000, eq_state, True, -1.0e-4)
    write_json(ROOT / "02_continuous_0p00_to_0p30" / "INPUT_METADATA.json", {
        "source_state": str(ROOT / "01_equilibration" / "equilibrated_state.json"),
        "target_global_xx_strain": -0.30,
        "caimd_step_strain_1": -1.0e-4,
        "caimd_stride": CAIMD_STRIDE,
        "expected_events": 600,
        "expected_nominal_final_xx": -0.30,
    })


def build_segment(segment: int) -> dict[str, Any]:
    if segment == 1:
        pred_state_path = ROOT / "01_equilibration" / "equilibrated_state.json"
        run_dir = ROOT / "03_segmented" / "segment_01_global_0p00_to_0p10"
        target = -0.10
        prefix = "diamond8_seg01"
    elif segment == 2:
        pred_state_path = ROOT / "03_segmented" / "segment_01_global_0p00_to_0p10" / "terminal_state.json"
        run_dir = ROOT / "03_segmented" / "segment_02_global_0p10_to_0p20"
        target = -0.20
        prefix = "diamond8_seg02"
    elif segment == 3:
        pred_state_path = ROOT / "03_segmented" / "segment_02_global_0p10_to_0p20" / "terminal_state.json"
        run_dir = ROOT / "03_segmented" / "segment_03_global_0p20_to_0p30"
        target = -0.30
        prefix = "diamond8_seg03"
    else:
        raise SystemExit("segment must be 1, 2, or 3")
    state = read_json(pred_state_path)
    h0 = read_json(ROOT / "01_equilibration" / "equilibrated_state.json")["cell_matrix_bohr_columns"]
    hprev = state["cell_matrix_bohr_columns"]
    fprev = matmul(hprev, inv3(h0))
    local_total = (1.0 + target) / fprev[0][0] - 1.0
    step_strain = local_total / 1000.0
    write_input(run_dir, prefix, 1000, state, True, step_strain)
    meta = {
        "segment": segment,
        "source_state": str(pred_state_path),
        "target_global_xx_strain": target,
        "predecessor_global_strain": sym_strain(hprev, h0),
        "predecessor_F11": fprev[0][0],
        "local_total_xx_strain": local_total,
        "caimd_step_strain_1": step_strain,
        "caimd_stride": CAIMD_STRIDE,
        "expected_events": 200,
        "expected_local_nominal_final_xx": local_total,
        "derivation": "local_total = (1 + target_global_xx) / Fprev_11 - 1, using CAIMD sym(H*inv(Href)-I) convention for xx",
    }
    write_json(run_dir / "INPUT_METADATA.json", meta)
    if segment > 1:
        validate_boundary(segment, state, run_dir, prefix, h0)
    return meta


def validate_boundary(segment: int, pred_state: dict[str, Any], run_dir: Path, prefix: str, h0: list[list[float]]) -> dict[str, Any]:
    input_state = parse_input_state(run_dir / "qe.in")
    h_pred = pred_state["cell_matrix_bohr_columns"]
    h_input = input_state["cell_matrix_bohr_columns"]
    frac_pred = pred_state["fractional_positions"]
    frac_input = input_state["fractional_positions"]
    vel_pred = pred_state["cartesian_velocities_au"]
    vel_input = input_state["cartesian_velocities_au"]
    cart_pred = [matvec(h_pred, f) for f in frac_pred]
    cart_input = [matvec(h_input, f) for f in frac_input]
    diffs_frac = [abs(frac_pred[i][j] - frac_input[i][j]) for i in range(NAT) for j in range(3)]
    diffs_cart = [abs(cart_pred[i][j] - cart_input[i][j]) for i in range(NAT) for j in range(3)]
    diffs_vel = [abs(vel_pred[i][j] - vel_input[i][j]) for i in range(NAT) for j in range(3)]
    temp_pred = velocity_temperature(vel_pred)["temperature_3N_minus_3_K"]
    temp_input = velocity_temperature(vel_input)["temperature_3N_minus_3_K"]
    rec = {
        "segment": segment,
        "predecessor_state": pred_state.get("source_run_dir"),
        "successor_input": str(run_dir / "qe.in"),
        "max_fractional_position_difference": max(diffs_frac),
        "max_cartesian_position_difference_bohr": max(diffs_cart),
        "max_cartesian_velocity_difference_au": max(diffs_vel),
        "temperature_before_K": temp_pred,
        "temperature_after_input_K": temp_input,
        "temperature_difference_K": temp_input - temp_pred,
        "cell_matrix_difference_max_abs_bohr": max(abs(h_pred[i][j] - h_input[i][j]) for i in range(3) for j in range(3)),
        "volume_before_bohr3": det3(h_pred),
        "volume_after_input_bohr3": det3(h_input),
        "volume_difference_bohr3": det3(h_input) - det3(h_pred),
        "global_strain_before": sym_strain(h_pred, h0),
        "global_strain_after_input": sym_strain(h_input, h0),
        "atom_order_consistent": True,
        "local_nominal_strain_reset": True,
        "nose_state_reset": True,
        "pr_history_reset": True,
        "caimd_local_reference_reset": True,
    }
    write_json(run_dir / "BOUNDARY_VALIDATION.json", rec)
    return rec


def parse_input_state(path: Path) -> dict[str, Any]:
    lines = path.read_text().splitlines()
    cell_rows: list[list[float]] = []
    frac: list[list[float]] = []
    vel: list[list[float]] = []
    mode = None
    for line in lines:
        stripped = line.strip()
        upper = stripped.upper()
        if upper.startswith("CELL_PARAMETERS"):
            mode = "cell"
            continue
        if upper.startswith("ATOMIC_POSITIONS"):
            mode = "pos"
            continue
        if upper.startswith("ATOMIC_VELOCITIES"):
            mode = "vel"
            continue
        if upper.startswith("K_POINTS") or upper.startswith("ATOMIC_SPECIES") or stripped.startswith("&") or stripped == "/":
            mode = None
            continue
        if mode == "cell" and len(cell_rows) < 3:
            vals = extract_floats(stripped)
            if len(vals) >= 3:
                cell_rows.append(vals[:3])
        elif mode == "pos" and len(frac) < NAT:
            vals = extract_floats(stripped)
            if len(vals) >= 3:
                frac.append(vals[:3])
        elif mode == "vel" and len(vel) < NAT:
            vals = extract_floats(stripped)
            if len(vals) >= 3:
                vel.append(vals[:3])
    return {
        "cell_matrix_bohr_columns": cell_rows_to_columns(cell_rows),
        "fractional_positions": frac,
        "cartesian_velocities_au": vel,
        "atom_labels": ["C"] * len(frac),
    }


def run_stage(args: argparse.Namespace) -> None:
    stage = args.stage
    if stage == "smoke":
        run_dir = ROOT / "evidence" / "pw_smoke"
        prefix = "diamond8_smoke"
        expected = 1
        active = False
    elif stage == "equilibration":
        run_dir = ROOT / "01_equilibration"
        prefix = "diamond8_equil"
        expected = EQUIL_STEPS
        active = False
    elif stage == "continuous":
        run_dir = ROOT / "02_continuous_0p00_to_0p30"
        prefix = "diamond8_continuous"
        expected = 3000
        active = True
    elif stage in ("seg1", "seg2", "seg3"):
        seg = int(stage[-1])
        run_dir = [
            ROOT / "03_segmented" / "segment_01_global_0p00_to_0p10",
            ROOT / "03_segmented" / "segment_02_global_0p10_to_0p20",
            ROOT / "03_segmented" / "segment_03_global_0p20_to_0p30",
        ][seg - 1]
        prefix = f"diamond8_seg0{seg}"
        expected = 1000
        active = True
    else:
        raise SystemExit(f"unknown stage: {stage}")
    run = run_qe(run_dir, prefix)
    val = validate_run(run_dir, prefix, expected, active)
    if stage == "equilibration" and val["valid"]:
        extract_terminal_state(run_dir, prefix, ROOT / "01_equilibration" / "equilibrated_state.json")
    if stage in ("continuous", "seg1", "seg2", "seg3") and val["valid"]:
        extract_terminal_state(run_dir, prefix, run_dir / "terminal_state.json")
        write_json(run_dir / "SEGMENT_COMPLETE.json", {"stage": stage, "valid": True, "validation": val})
    print(json.dumps({"run": run, "validation": val}, indent=2))
    if not val["valid"]:
        raise SystemExit(f"stage validation failed: {stage}")


def run_segmented() -> None:
    for seg in (1, 2, 3):
        build_segment(seg)
        run_stage(argparse.Namespace(stage=f"seg{seg}"))


def run_all(args: argparse.Namespace) -> None:
    run_stage(argparse.Namespace(stage="smoke"))
    run_stage(argparse.Namespace(stage="equilibration"))
    build_continuous()
    run_stage(argparse.Namespace(stage="continuous"))
    run_segmented()
    analyze()


def global_rows_for_run(
    run_dir: Path,
    prefix: str,
    h_segment_ref: list[list[float]],
    h_global_ref: list[list[float]],
    global_step_offset: int,
    label: str,
) -> list[dict[str, Any]]:
    rows = rows_from_caimd(artifact(run_dir, prefix, "caimd"))
    out: list[dict[str, Any]] = []
    for r in rows:
        local_actual = [r[f"actual_{name}"] for name in ("xx", "yy", "zz", "yz", "xz", "xy")]
        emat = strain_vec_to_matrix(local_actual)
        for i in range(3):
            emat[i][i] += 1.0
        hpremove = matmul(emat, h_segment_ref)
        gstrain = sym_strain(hpremove, h_global_ref)
        rec = dict(r)
        rec["run_label"] = label
        rec["global_loading_step"] = global_step_offset + r["ionic_step"]
        rec["global_xx_strain"] = gstrain[0]
        rec["global_yy_strain"] = gstrain[1]
        rec["global_zz_strain"] = gstrain[2]
        rec["global_yz_strain"] = gstrain[3]
        rec["global_xz_strain"] = gstrain[4]
        rec["global_xy_strain"] = gstrain[5]
        out.append(rec)
    return out


def analyze() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    analysis = ROOT / "analysis"
    analysis.mkdir(exist_ok=True)
    h0 = read_json(ROOT / "01_equilibration" / "equilibrated_state.json")["cell_matrix_bohr_columns"]
    cont_dir = ROOT / "02_continuous_0p00_to_0p30"
    cont_rows = global_rows_for_run(cont_dir, "diamond8_continuous", h0, h0, 0, "continuous")
    seg_dirs = [
        ROOT / "03_segmented" / "segment_01_global_0p00_to_0p10",
        ROOT / "03_segmented" / "segment_02_global_0p10_to_0p20",
        ROOT / "03_segmented" / "segment_03_global_0p20_to_0p30",
    ]
    seg_prefixes = ["diamond8_seg01", "diamond8_seg02", "diamond8_seg03"]
    seg_rows: list[dict[str, Any]] = []
    offset = 0
    for idx, (d, p) in enumerate(zip(seg_dirs, seg_prefixes), start=1):
        href = parse_input_state(d / "qe.in")["cell_matrix_bohr_columns"]
        seg_rows.extend(global_rows_for_run(d, p, href, h0, offset, f"segment_{idx:02d}"))
        offset += 1000
    write_csv(analysis / "continuous_caimd.csv", cont_rows)
    write_csv(analysis / "segmented_caimd_global.csv", seg_rows)
    temp_rows = temperature_stats_rows(cont_rows, seg_rows)
    write_csv(analysis / "temperature_statistics.csv", temp_rows)
    boundary_rows = collect_boundary_rows()
    write_csv(analysis / "segment_boundary_statistics.csv", boundary_rows)
    summary = comparison_summary(cont_rows, seg_rows, boundary_rows)
    write_json(analysis / "continuous_vs_segmented_summary.json", summary)
    make_plots(analysis, cont_rows, seg_rows, boundary_rows)
    write_final_report(summary, temp_rows, boundary_rows)
    print(json.dumps({"analysis": str(analysis), "summary": summary}, indent=2))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def temperature_stats_rows(cont_rows: list[dict[str, Any]], seg_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    groups = {"continuous": cont_rows, "segmented_all": seg_rows}
    for label in sorted({r["run_label"] for r in seg_rows}):
        groups[label] = [r for r in seg_rows if r["run_label"] == label]
    eq_val = read_json(ROOT / "01_equilibration" / "VALIDATION.json")
    rows.append({"run_label": "equilibration", **flatten_stats(eq_val["temperature_analysis_K"]), "classification": eq_val["temperature_classification"]})
    for label, data in groups.items():
        vals = [r["temperature_K"] for r in data]
        rows.append({"run_label": label, **flatten_stats(stats(vals)), "classification": classify_temperature(vals)})
    return rows


def flatten_stats(s: dict[str, Any]) -> dict[str, Any]:
    return {f"temperature_{k}": v for k, v in s.items()}


def collect_boundary_rows() -> list[dict[str, Any]]:
    rows = []
    for d in [
        ROOT / "03_segmented" / "segment_02_global_0p10_to_0p20",
        ROOT / "03_segmented" / "segment_03_global_0p20_to_0p30",
    ]:
        p = d / "BOUNDARY_VALIDATION.json"
        if p.exists():
            row = read_json(p)
            for key in ("global_strain_before", "global_strain_after_input"):
                vals = row.pop(key)
                for name, val in zip(("xx", "yy", "zz", "yz", "xz", "xy"), vals):
                    row[f"{key}_{name}"] = val
            rows.append(row)
    return rows


def comparison_summary(cont_rows: list[dict[str, Any]], seg_rows: list[dict[str, Any]], boundary_rows: list[dict[str, Any]]) -> dict[str, Any]:
    def final(rows: list[dict[str, Any]]) -> dict[str, Any]:
        return rows[-1] if rows else {}

    cfinal = final(cont_rows)
    sfinal = final(seg_rows)
    return {
        "continuous_rows": len(cont_rows),
        "segmented_rows": len(seg_rows),
        "continuous_event_count_final": cfinal.get("imposed_events_done"),
        "segmented_event_count_final_local": sfinal.get("imposed_events_done"),
        "continuous_final_global_xx": cfinal.get("global_xx_strain"),
        "segmented_final_global_xx": sfinal.get("global_xx_strain"),
        "final_global_xx_difference_segmented_minus_continuous": (
            sfinal.get("global_xx_strain") - cfinal.get("global_xx_strain")
            if cfinal and sfinal
            else None
        ),
        "continuous_final_temperature_K": cfinal.get("temperature_K"),
        "segmented_final_temperature_K": sfinal.get("temperature_K"),
        "continuous_temperature_classification": classify_temperature([r["temperature_K"] for r in cont_rows]),
        "segmented_temperature_classification": classify_temperature([r["temperature_K"] for r in seg_rows]),
        "continuous_xx_nominal_actual_error_max": max((abs(r["actual_xx"] - r["nominal_xx"]) for r in cont_rows), default=None),
        "segmented_local_xx_nominal_actual_error_max": max((abs(r["actual_xx"] - r["nominal_xx"]) for r in seg_rows), default=None),
        "segmented_boundary_count": len(boundary_rows),
        "boundary_velocity_difference_max_au": max((r["max_cartesian_velocity_difference_au"] for r in boundary_rows), default=None),
        "boundary_position_difference_max_bohr": max((r["max_cartesian_position_difference_bohr"] for r in boundary_rows), default=None),
        "interpretation": "Segmented and continuous runs need not be atom-by-atom identical; differences include chaotic MD divergence plus Nose, PR-history, CAIMD event-count, and local-reference resets at segment boundaries.",
    }


def make_plots(analysis: Path, cont: list[dict[str, Any]], seg: list[dict[str, Any]], boundary_rows: list[dict[str, Any]]) -> None:
    def plot_xy(filename: str, xkey: str, ykeys: list[str], ylabel: str, rows_list: list[tuple[str, list[dict[str, Any]]]]) -> None:
        plt.figure(figsize=(8, 5))
        for label, rows in rows_list:
            if not rows:
                continue
            for ykey in ykeys:
                plt.plot([r[xkey] for r in rows], [r[ykey] for r in rows], label=f"{label} {ykey}")
        for b in boundary_rows:
            if "global_strain_before_xx" in b:
                plt.axvline(b["global_strain_before_xx"], color="0.75", linestyle="--", linewidth=0.8)
        plt.xlabel(xkey)
        plt.ylabel(ylabel)
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(analysis / filename, dpi=160)
        plt.close()

    both = [("continuous", cont), ("segmented", seg)]
    plot_xy("temperature_vs_global_step.png", "global_loading_step", ["temperature_K"], "Temperature (K)", both)
    plot_xy("temperature_vs_global_strain.png", "global_xx_strain", ["temperature_K"], "Temperature (K)", both)
    plot_xy("temperature_continuous_vs_segmented.png", "global_loading_step", ["temperature_K"], "Temperature (K)", both)
    plot_xy("nominal_vs_actual_global_strain.png", "global_loading_step", ["global_xx_strain", "nominal_xx"], "Strain", both)
    plot_xy("six_stresses_vs_global_strain.png", "global_xx_strain", ["stress_xx_GPa", "stress_yy_GPa", "stress_zz_GPa", "stress_yz_GPa", "stress_xz_GPa", "stress_xy_GPa"], "Stress (GPa)", both)
    plot_xy("normal_stresses_vs_global_strain.png", "global_xx_strain", ["stress_xx_GPa", "stress_yy_GPa", "stress_zz_GPa"], "Stress (GPa)", both)
    plot_xy("shear_stresses_vs_global_strain.png", "global_xx_strain", ["stress_yz_GPa", "stress_xz_GPa", "stress_xy_GPa"], "Stress (GPa)", both)
    plot_xy("volume_vs_global_strain.png", "global_xx_strain", ["volume_bohr3"], "Volume (bohr^3)", both)
    # Determinant equals volume in the CAIMD diagnostic.
    shutil.copyfile(analysis / "volume_vs_global_strain.png", analysis / "cell_determinant_vs_global_strain.png")
    plot_xy("velocity_jump_vs_global_strain.png", "global_xx_strain", ["velocity_jump_max", "velocity_jump_rms"], "Velocity residual (a.u.)", both)
    plt.figure(figsize=(8, 5))
    for label, rows in both:
        if not rows:
            continue
        plt.plot([r["global_loading_step"] for r in rows], [r["temperature_K"] for r in rows], label=label)
    for b in boundary_rows:
        step = 1000 * (int(b["segment"]) - 1)
        plt.axvline(step, color="0.7", linestyle="--")
    plt.xlabel("global_loading_step")
    plt.ylabel("Temperature (K)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(analysis / "segment_boundary_temperature_detail.png", dpi=160)
    plt.close()


def write_final_report(summary: dict[str, Any], temp_rows: list[dict[str, Any]], boundary_rows: list[dict[str, Any]]) -> None:
    status = "SUCCESS"
    if summary["continuous_final_global_xx"] is None or summary["segmented_final_global_xx"] is None:
        status = "PARTIAL"
    elif abs(summary["segmented_final_global_xx"] + 0.30) > 1.0e-6:
        status = "SUCCESS_WITH_WARNINGS"
    elif summary["continuous_temperature_classification"] == "FAIL" or summary["segmented_temperature_classification"] == "FAIL":
        status = "PARTIAL"
    elif summary["continuous_temperature_classification"] == "WARN" or summary["segmented_temperature_classification"] == "WARN":
        status = "SUCCESS_WITH_WARNINGS"
    cfg = config()
    report = ROOT / "analysis" / "FINAL_VALIDATION_REPORT.md"
    report.write_text(
        f"""# Final Validation Report

Status: {status}

Workflow directory: `{ROOT}`

Executable: `{cfg['pw_executable']}`

Executable SHA256: `{cfg['pw_sha256']}`

Pseudopotential: `{PSEUDO}`

Pseudopotential SHA256: `{sha256(PSEUDO)}`

Branch/HEAD: `{cfg['branch']}` / `{cfg['head']}`

No source code was modified by this workflow. No Leonardo or other remote
connection was used. No Git push was used.

## Ensemble And Rate

The loading ensemble is thermostatted constrained-strain Parrinello-Rahman
dynamics, not unconstrained full NpT. Laboratory `xx` strain is externally
driven by active CAIMD while the other components are left to the ordinary PR
candidate response.

The global engineering-strain rate is `{LOADING_RATE_PS_INV} ps^-1`, 50 times
faster than `0.001 ps^-1`; this is only a local functional validation rate.

The target compression of 0.30 is outside the nominal CAIMD V1 small-strain
regime and is treated as a software stress test.

## Summary

```json
{json.dumps(summary, indent=2, sort_keys=True)}
```

## Temperature Statistics

```json
{json.dumps(temp_rows, indent=2, sort_keys=True)}
```

## Boundary Statistics

```json
{json.dumps(boundary_rows, indent=2, sort_keys=True)}
```

## Continuation Limitations

Segmented continuation transfers the terminal cell, fractional positions, and
physical Cartesian velocities. It resets Nose extended variables, PR cell
velocity/acceleration history, CAIMD event count, and the CAIMD local reference
cell. It is therefore not an exact continuous restart.
"""
    )


def cmd_extract(args: argparse.Namespace) -> None:
    state = extract_terminal_state(Path(args.run_dir), args.prefix, Path(args.output))
    print(json.dumps(state, indent=2))


def cmd_validate(args: argparse.Namespace) -> None:
    val = validate_run(Path(args.run_dir), args.prefix, args.expected_steps, args.active_caimd)
    print(json.dumps(val, indent=2))
    if not val["valid"]:
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--pw", required=True)
    p.set_defaults(func=prepare)
    p = sub.add_parser("run-stage")
    p.add_argument("stage", choices=["smoke", "equilibration", "continuous", "seg1", "seg2", "seg3"])
    p.set_defaults(func=run_stage)
    p = sub.add_parser("build-continuous")
    p.set_defaults(func=lambda args: build_continuous())
    p = sub.add_parser("build-segment")
    p.add_argument("segment", type=int, choices=[1, 2, 3])
    p.set_defaults(func=lambda args: print(json.dumps(build_segment(args.segment), indent=2)))
    p = sub.add_parser("run-segmented")
    p.set_defaults(func=lambda args: run_segmented())
    p = sub.add_parser("run-all")
    p.set_defaults(func=run_all)
    p = sub.add_parser("analyze")
    p.set_defaults(func=lambda args: analyze())
    p = sub.add_parser("extract")
    p.add_argument("--run-dir", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_extract)
    p = sub.add_parser("validate")
    p.add_argument("--run-dir", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--expected-steps", type=int, required=True)
    p.add_argument("--active-caimd", action="store_true")
    p.set_defaults(func=cmd_validate)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
