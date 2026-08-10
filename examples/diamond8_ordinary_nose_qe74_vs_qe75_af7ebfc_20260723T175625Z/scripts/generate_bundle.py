#!/usr/bin/env python3
"""Generate the Leonardo QE 7.4.1 vs QE 7.5 ordinary vc-md Nose bundle."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parents[1]
CONTROL_ROOT = REPO_ROOT / "examples" / "caimd_diamond8_local_nose_controls_ABC_20260723T134303Z"
COMMON_STATE = CONTROL_ROOT / "01_common_state" / "COMMON_STATE.json"
PSEUDO_SOURCE = REPO_ROOT / "pseudo" / "C.UPF"
LOCAL_PW = Path("/private/tmp/qe_diamond8_af7ebfc_build/bin/pw.x")

REMOTE_WORKDIR = "/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_ordinary_nose_qe74_vs_qe75_af7ebfc_20260723T175625Z"
QE74_MODULE = "quantum-espresso/7.4.1--hpcx-mpi--2.19--nvhpc--24.5-openblas"
QE75_EXE = "/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc/install-current-gpu/bin/pw.x"
QE75_ROOT = "/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc"

NAT = 8
MASS_C_AMU = 12.0107
ND_I = 21
DT_QE = 41.34137333665
DT_FS = 2.0
NSTEP = 1000
TEMPW = 300.0
FNOSEP = 50.0
NHPCL = 4
WMASS = 2000000.0
ECUTWFC = 25.0
ECUTRHO = 120.0
CONV_THR = 1.0e-4
ELECTRON_MAXSTEP = 40
MIXING_BETA = 0.3

RYDBERG_SI = 2.1798723611035e-18
K_BOLTZMANN_SI = 1.380649e-23
AMU_SI = 1.66053906660e-27
ELECTRONMASS_SI = 9.1093837015e-31
AMU_RY = (AMU_SI / ELECTRONMASS_SI) / 2.0
RY_TO_KELVIN = RYDBERG_SI / K_BOLTZMANN_SI


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run_text(cmd: list[str], cwd: Path = REPO_ROOT) -> str:
    return subprocess.check_output(cmd, cwd=str(cwd), text=True, stderr=subprocess.STDOUT)


def write_json(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def dfloat(value: float) -> str:
    return f"{value:.16e}".replace("e", "d")


def kinetic_temperature(velocities: list[list[float]]) -> dict[str, object]:
    com = [sum(row[i] for row in velocities) / len(velocities) for i in range(3)]
    centered = [[row[i] - com[i] for i in range(3)] for row in velocities]
    mass_ry = MASS_C_AMU * AMU_RY
    ekin = 0.5 * mass_ry * sum(sum(v * v for v in row) for row in centered)
    return {
        "temperature_K": 2.0 * ekin / ND_I * RY_TO_KELVIN,
        "ekin_Ry": ekin,
        "center_of_mass_velocity_au": com,
        "center_of_mass_max_abs_au": max(abs(x) for x in com),
        "ndof": ND_I,
    }


def input_text(prefix: str, state: dict[str, object]) -> str:
    cell_rows = state["cell_parameters_bohr_rows_for_qe"]
    frac = state["fractional_positions"]
    velocities = state["cartesian_velocities_au"]
    lines = [
        "&CONTROL",
        "  calculation = 'vc-md',",
        "  restart_mode = 'from_scratch',",
        f"  prefix = '{prefix}',",
        "  pseudo_dir = '../pseudo',",
        "  outdir = './tmp',",
        f"  nstep = {NSTEP},",
        f"  dt = {dfloat(DT_QE)},",
        "  iprint = 1,",
        "  tstress = .true.,",
        "  tprnfor = .true.,",
        "  disk_io = 'low',",
        "  verbosity = 'high',",
        "/",
        "&SYSTEM",
        "  ibrav = 0,",
        f"  nat = {NAT},",
        "  ntyp = 1,",
        f"  ecutwfc = {ECUTWFC:.1f},",
        f"  ecutrho = {ECUTRHO:.1f},",
        "  occupations = 'smearing',",
        "  smearing = 'gaussian',",
        "  degauss = 0.01,",
        "  nosym = .true.,",
        "  noinv = .true.,",
        "/",
        "&ELECTRONS",
        f"  conv_thr = {dfloat(CONV_THR)},",
        f"  electron_maxstep = {ELECTRON_MAXSTEP},",
        f"  mixing_beta = {MIXING_BETA:.3f},",
        "/",
        "&IONS",
        "  ion_dynamics = 'beeman',",
        "  ion_temperature = 'nose',",
        "  ion_velocities = 'from_input',",
        f"  tempw = {TEMPW:.1f},",
        f"  fnosep(1) = {FNOSEP:.8f},",
        f"  nhpcl = {NHPCL},",
        "  ndega = 0,",
        "/",
        "&CELL",
        "  cell_dynamics = 'pr',",
        "  cell_dofree = 'all',",
        "  press = 0.0,",
        f"  wmass = {WMASS:.8f},",
        "  cell_temperature = 'not_controlled',",
        "/",
        "ATOMIC_SPECIES",
        f"C {MASS_C_AMU:.8f} C.UPF",
        "CELL_PARAMETERS bohr",
    ]
    for row in cell_rows:
        lines.append("  " + "  ".join(f"{x:.14f}" for x in row))
    lines.append("ATOMIC_POSITIONS crystal")
    for xyz in frac:
        lines.append("C  " + "  ".join(f"{x:.14f}" for x in xyz))
    lines.append("K_POINTS gamma")
    lines.append("ATOMIC_VELOCITIES { a.u }")
    for xyz in velocities:
        lines.append("C  " + "  ".join(f"{x:.14e}" for x in xyz))
    text = "\n".join(lines) + "\n"
    if "caimd" in text.lower():
        raise RuntimeError("generated input unexpectedly contains a forbidden keyword")
    return text


VALIDATE_SCRIPT = r'''#!/usr/bin/env python3

import hashlib
import json
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RYDBERG_SI = 2.1798723611035e-18
K_BOLTZMANN_SI = 1.380649e-23
AMU_SI = 1.66053906660e-27
ELECTRONMASS_SI = 9.1093837015e-31
AMU_RY = (AMU_SI / ELECTRONMASS_SI) / 2.0
RY_TO_KELVIN = RYDBERG_SI / K_BOLTZMANN_SI
MASS_C_AMU = 12.0107
NAT = 8
ND_I = 21
FLOAT_RE = re.compile(r"[-+]?(?:(?:\d+\.\d*|\.\d+)(?:[EeDd][-+]?\d+)?|\d+(?:[EeDd][-+]?\d+)?)")

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def floats(line):
    return [float(x.replace("D", "E").replace("d", "e")) for x in FLOAT_RE.findall(line)]

def parse_input(path):
    text = path.read_text()
    low = text.lower()
    forbidden = ["caimd", "caimd_stride", "caimd_mask", "caimd_step_strain"]
    hits = [word for word in forbidden if word in low]
    lines = text.splitlines()
    cell = []
    pos = []
    vel = []
    mode = None
    settings = {}
    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        up = line.upper()
        if up.startswith("CELL_PARAMETERS"):
            mode = "cell"
            continue
        if up.startswith("ATOMIC_POSITIONS"):
            mode = "pos"
            continue
        if up.startswith("ATOMIC_VELOCITIES"):
            mode = "vel"
            continue
        if up.startswith("K_POINTS"):
            mode = None
        if mode == "cell" and len(cell) < 3:
            vals = floats(line)
            if len(vals) >= 3:
                cell.append(vals[:3])
            continue
        if mode == "pos" and len(pos) < NAT:
            vals = floats(line)
            if len(vals) >= 3:
                pos.append(vals[-3:])
            continue
        if mode == "vel" and len(vel) < NAT:
            vals = floats(line)
            if len(vals) >= 3:
                vel.append(vals[-3:])
            continue
        if "=" in line and not line.startswith("!"):
            key = line.split("=", 1)[0].strip().lower()
            val = line.split("=", 1)[1].strip().rstrip(",").strip()
            if key not in ("prefix", "outdir"):
                settings[key] = val
    mass_ry = MASS_C_AMU * AMU_RY
    com = [sum(row[i] for row in vel) / len(vel) for i in range(3)] if vel else [None, None, None]
    centered = [[row[i] - com[i] for i in range(3)] for row in vel] if vel else []
    ekin = 0.5 * mass_ry * sum(sum(v * v for v in row) for row in centered)
    temp = 2.0 * ekin / ND_I * RY_TO_KELVIN if vel else None
    return {
        "path": str(path),
        "sha256": sha256(path),
        "forbidden_hits": hits,
        "cell": cell,
        "positions": pos,
        "velocities": vel,
        "settings": settings,
        "nat_from_positions": len(pos),
        "nat_from_velocities": len(vel),
        "temperature_K": temp,
        "ekin_Ry": ekin,
        "center_of_mass_velocity_au": com,
    }

def max_abs_diff(a, b):
    diffs = []
    for x, y in zip(a, b):
        if isinstance(x, list):
            diffs.append(max_abs_diff(x, y))
        else:
            diffs.append(abs(float(x) - float(y)))
    return max(diffs) if diffs else 0.0

def normalized_diff(a, b):
    sa = a["settings"]
    sb = b["settings"]
    keys = sorted(set(sa) | set(sb))
    return {key: [sa.get(key), sb.get(key)] for key in keys if sa.get(key) != sb.get(key)}

def main():
    qe74 = parse_input(ROOT / "inputs" / "qe74.in")
    qe75 = parse_input(ROOT / "inputs" / "qe75_af7ebfc.in")
    pseudo = ROOT / "pseudo" / "C.UPF"
    result = {
        "qe74": {k: v for k, v in qe74.items() if k not in ("cell", "positions", "velocities", "settings")},
        "qe75_af7ebfc": {k: v for k, v in qe75.items() if k not in ("cell", "positions", "velocities", "settings")},
        "normalized_diff": normalized_diff(qe74, qe75),
        "cell_difference_max_abs": max_abs_diff(qe74["cell"], qe75["cell"]),
        "position_difference_max_abs": max_abs_diff(qe74["positions"], qe75["positions"]),
        "velocity_difference_max_abs": max_abs_diff(qe74["velocities"], qe75["velocities"]),
        "thermostat_parameter_difference": {k: v for k, v in normalized_diff(qe74, qe75).items() if k in ("ion_temperature", "tempw", "fnosep(1)", "nhpcl", "ndega")},
        "pr_parameter_difference": {k: v for k, v in normalized_diff(qe74, qe75).items() if k in ("cell_dynamics", "cell_dofree", "press", "wmass", "cell_temperature", "dt", "nstep")},
        "electronic_parameter_difference": {k: v for k, v in normalized_diff(qe74, qe75).items() if k in ("ecutwfc", "ecutrho", "conv_thr", "electron_maxstep", "mixing_beta", "occupations", "smearing", "degauss", "nosym", "noinv")},
        "pseudo_sha256": sha256(pseudo),
        "all_physical_differences_zero": False,
    }
    result["all_physical_differences_zero"] = (
        result["cell_difference_max_abs"] == 0.0
        and result["position_difference_max_abs"] == 0.0
        and result["velocity_difference_max_abs"] == 0.0
        and not result["thermostat_parameter_difference"]
        and not result["pr_parameter_difference"]
        and not result["electronic_parameter_difference"]
        and not qe74["forbidden_hits"]
        and not qe75["forbidden_hits"]
        and qe74["nat_from_positions"] == 8
        and qe75["nat_from_positions"] == 8
        and abs(float(qe74["temperature_K"]) - 300.0) < 1.0e-9
        and abs(float(qe75["temperature_K"]) - 300.0) < 1.0e-9
    )
    out = ROOT / "evidence" / "input_identity_check.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if not result["all_physical_differences_zero"]:
        print(json.dumps(result, indent=2, sort_keys=True))
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
'''


ANALYZE_SCRIPT = r'''#!/usr/bin/env python3

import json
import math
import re
import statistics
import struct
import sys
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NAT = 8
ND_I = 21
MASS_C_AMU = 12.0107
DT_FS = 2.0
RYDBERG_SI = 2.1798723611035e-18
K_BOLTZMANN_SI = 1.380649e-23
AMU_SI = 1.66053906660e-27
ELECTRONMASS_SI = 9.1093837015e-31
AMU_RY = (AMU_SI / ELECTRONMASS_SI) / 2.0
RY_TO_KELVIN = RYDBERG_SI / K_BOLTZMANN_SI
FLOAT_RE = re.compile(r"[-+]?(?:(?:\d+\.\d*|\.\d+)(?:[EeDd][-+]?\d+)?|\d+(?:[EeDd][-+]?\d+)?)")

RUNS = {
    "qe74": {"dir": "qe74_official", "prefix": "diamond8_qe74_official", "label": "QE 7.4.1 official"},
    "qe75_af7ebfc": {"dir": "qe75_af7ebfc_no_caimd", "prefix": "diamond8_qe75_af7ebfc", "label": "QE 7.5 af7ebfc ordinary"},
}

def floats(text):
    return [float(x.replace("D", "E").replace("d", "e")) for x in FLOAT_RE.findall(text)]

def det3(a):
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )

def inv3(a):
    d = det3(a)
    if abs(d) < 1.0e-20 or not math.isfinite(d):
        raise ValueError("singular cell")
    return [
        [(a[1][1]*a[2][2]-a[1][2]*a[2][1])/d, -(a[0][1]*a[2][2]-a[0][2]*a[2][1])/d, (a[0][1]*a[1][2]-a[0][2]*a[1][1])/d],
        [-(a[1][0]*a[2][2]-a[1][2]*a[2][0])/d, (a[0][0]*a[2][2]-a[0][2]*a[2][0])/d, -(a[0][0]*a[1][2]-a[0][2]*a[1][0])/d],
        [(a[1][0]*a[2][1]-a[1][1]*a[2][0])/d, -(a[0][0]*a[2][1]-a[0][1]*a[2][0])/d, (a[0][0]*a[1][1]-a[0][1]*a[1][0])/d],
    ]

def matvec(a, x):
    return [sum(a[i][j] * x[j] for j in range(3)) for i in range(3)]

def norm_inf(a):
    return max(sum(abs(v) for v in row) for row in a)

def condition(a):
    return norm_inf(a) * norm_inf(inv3(a))

def parse_blocks(path, rows_per_block):
    if not path.exists():
        return []
    lines = [line.rstrip() for line in path.read_text(errors="ignore").splitlines() if line.strip()]
    blocks = []
    i = 0
    while i < len(lines):
        head = floats(lines[i])
        if len(head) < 2:
            i += 1
            continue
        rows = []
        ok = True
        for j in range(rows_per_block):
            if i + 1 + j >= len(lines):
                ok = False
                break
            vals = floats(lines[i + 1 + j])
            if len(vals) < 3:
                ok = False
                break
            rows.append(vals[:3])
        if not ok:
            break
        blocks.append({"step": int(round(head[0])), "time_raw": head[1], "rows": rows})
        i += 1 + rows_per_block
    return blocks

def artifact(run_dir, prefix, ext):
    direct = run_dir / f"{prefix}.{ext}"
    if direct.exists():
        return direct
    return run_dir / "tmp" / f"{prefix}.{ext}"

def parse_eal(path: Path):
    rows = []
    if not path.exists():
        return rows
    for line in path.read_text(errors="ignore").splitlines():
        vals = floats(line)
        if len(vals) >= 7:
            uta, eka, eta, utl, ekla, etl, nst = vals[:7]
            step = int(round(nst)) + 1
            rows.append({
                "step": step,
                "time_ps": step * DT_FS / 1000.0,
                "uta_Ry": uta,
                "eka_Ry": eka,
                "eta_Ry": eta,
                "utl_Ry": utl,
                "ekla_Ry": ekla,
                "etl_Ry": etl,
                "temperature_K": 2.0 * eka / ND_I * RY_TO_KELVIN,
            })
    return rows

def parse_md(path: Path):
    if not path.exists():
        return None
    vals = floats(path.read_text(errors="ignore"))
    need = 5 * 3 * NAT
    if len(vals) < need:
        return None
    p = 3 * NAT
    ratd = [vals[p + 3*i:p + 3*i + 3] for i in range(NAT)]
    return {"ratd": ratd, "float_count": len(vals), "record_count": 1}

def kinetic_temperature(velocities):
    com = [sum(row[i] for row in velocities) / len(velocities) for i in range(3)]
    centered = [[row[i] - com[i] for i in range(3)] for row in velocities]
    mass_ry = MASS_C_AMU * AMU_RY
    ekin = 0.5 * mass_ry * sum(sum(v * v for v in row) for row in centered)
    return {"temperature_K": 2.0 * ekin / ND_I * RY_TO_KELVIN, "ekin_Ry": ekin, "center_of_mass_velocity_au": com}

def stats(values, times=None):
    vals = [v for v in values if math.isfinite(v)]
    out = {
        "count": len(values),
        "finite": len(vals) == len(values),
        "mean": sum(vals) / len(vals) if vals else None,
        "std": statistics.pstdev(vals) if len(vals) > 1 else (0.0 if len(vals) == 1 else None),
        "min": min(vals) if vals else None,
        "max": max(vals) if vals else None,
        "initial": vals[0] if vals else None,
        "final": vals[-1] if vals else None,
    }
    if times and len(times) == len(values) and len(vals) > 1:
        pairs = [(t, v) for t, v in zip(times, values) if math.isfinite(t) and math.isfinite(v)]
        mt = sum(t for t, _ in pairs) / len(pairs)
        mv = sum(v for _, v in pairs) / len(pairs)
        denom = sum((t - mt) ** 2 for t, _ in pairs)
        out["drift_per_ps"] = sum((t - mt) * (v - mv) for t, v in pairs) / denom if denom else None
    else:
        out["drift_per_ps"] = None
    return out

def first_step(values, threshold):
    for i, value in enumerate(values, start=1):
        if value > threshold:
            return i
    return None

def longest_run(values, predicate):
    best = 0
    cur = 0
    for value in values:
        if predicate(value):
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best

def return_codes():
    path = ROOT / "evidence" / "run_return_codes.txt"
    out = {}
    if not path.exists():
        return out
    for line in path.read_text(errors="ignore").splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            try:
                out[k.strip()] = int(v.strip())
            except ValueError:
                out[k.strip()] = v.strip()
    return out

def classify(result):
    temps = result["temperatures_K"]
    terminal_direct = result.get("terminal_direct_temperature_K")
    terminal_eal = temps[-1] if temps else None
    if not temps:
        return "INCOMPLETE"
    if not all(math.isfinite(t) for t in temps):
        return "NUMERICAL_FAILURE"
    if terminal_direct is not None and terminal_eal is not None:
        if abs(terminal_direct - terminal_eal) > max(100.0, 0.1 * max(abs(terminal_eal), 1.0)):
            return "DIAGNOSTIC_MISMATCH"
    if max(temps) > 2000.0 and terminal_direct is not None and terminal_direct > 2000.0:
        return "PHYSICAL_HEATING"
    if result.get("return_code") not in (0, None) and len(temps) < 20:
        return "NUMERICAL_FAILURE"
    if longest_run(temps, lambda t: t > 600.0) >= 25:
        return "PHYSICAL_HEATING"
    if longest_run(temps, lambda t: t < 150.0) >= 25:
        return "NUMERICAL_FAILURE"
    final = stats(temps).get("final")
    mean = stats(temps).get("mean")
    if final is not None and mean is not None and 200.0 <= final <= 450.0 and 200.0 <= mean <= 450.0:
        return "STABLE"
    return "INCOMPLETE"

def analyze_run(key, cfg, rcs):
    run_dir = ROOT / cfg["dir"]
    prefix = cfg["prefix"]
    out_text = (run_dir / "qe.out").read_text(errors="ignore") if (run_dir / "qe.out").exists() else ""
    eal = parse_eal(artifact(run_dir, prefix, "eal"))
    temps = [row["temperature_K"] for row in eal]
    times = [row["time_ps"] for row in eal]
    cel = parse_blocks(artifact(run_dir, prefix, "cel"), 3)
    md = parse_md(artifact(run_dir, prefix, "md"))
    direct = None
    if md and cel:
        velocities = [matvec(cel[-1]["rows"], row) for row in md["ratd"]]
        direct = kinetic_temperature(velocities)
    dets = []
    conds = []
    for block in cel:
        try:
            dets.append(det3(block["rows"]))
            conds.append(condition(block["rows"]))
        except Exception:
            dets.append(float("nan"))
            conds.append(float("inf"))
    result = {
        "key": key,
        "label": cfg["label"],
        "run_dir": str(run_dir),
        "return_code": rcs.get(f"RC_{key.upper()}"),
        "job_done": "JOB DONE" in out_text,
        "trajectory_records": len(cel),
        "simulated_time_ps": len(cel) * DT_FS / 1000.0,
        "eal_exists": artifact(run_dir, prefix, "eal").exists(),
        "eal_rows": len(eal),
        "md_exists": artifact(run_dir, prefix, "md").exists(),
        "md_record_count": md["record_count"] if md else 0,
        "temperatures_K": temps,
        "temperature_stats": stats(temps, times),
        "initial_temperature_K": temps[0] if temps else None,
        "mean_temperature_K": sum(temps) / len(temps) if temps else None,
        "final_temperature_K": temps[-1] if temps else None,
        "maximum_temperature_K": max(temps) if temps else None,
        "first_step_above_600K": first_step(temps, 600.0),
        "first_step_above_2000K": first_step(temps, 2000.0),
        "temperature_drift_K_per_ps": stats(temps, times).get("drift_per_ps"),
        "ionic_ekin_initial_Ry": eal[0]["eka_Ry"] if eal else None,
        "ionic_ekin_final_Ry": eal[-1]["eka_Ry"] if eal else None,
        "terminal_direct_temperature_K": direct["temperature_K"] if direct else None,
        "terminal_direct_ekin_Ry": direct["ekin_Ry"] if direct else None,
        "terminal_eal_direct_abs_diff_K": abs((temps[-1] if temps else float("nan")) - direct["temperature_K"]) if direct and temps else None,
        "cell_determinant_min": min(dets) if dets else None,
        "cell_condition_number_max": max(conds) if conds else None,
        "volume_range": [min(dets), max(dets)] if dets else None,
        "qe_error_lines": [line.strip() for line in out_text.splitlines() if "Error in routine" in line or "STOP" in line or "convergence NOT achieved" in line],
        "first_clear_numerical_instability": None,
    }
    if result["first_step_above_2000K"]:
        result["first_clear_numerical_instability"] = {"observable": "temperature_K", "step": result["first_step_above_2000K"], "threshold": 2000.0}
    result["classification"] = classify(result)
    return result, eal, dets

def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w") as handle:
        handle.write(",".join(fields) + "\n")
        for row in rows:
            handle.write(",".join(str(row.get(field, "")) for field in fields) + "\n")

def png_chunk(kind, data):
    return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xffffffff)

def save_png(path, width, height, pixels):
    raw = b"".join(b"\x00" + bytes(row) for row in pixels)
    data = b"\x89PNG\r\n\x1a\n"
    data += png_chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
    data += png_chunk(b"IDAT", zlib.compress(raw, 9))
    data += png_chunk(b"IEND", b"")
    path.write_bytes(data)

def draw_line(pixels, x0, y0, x1, y1, color):
    w = len(pixels[0]) // 3
    h = len(pixels)
    dx = abs(x1 - x0)
    dy = -abs(y1 - y0)
    sx = 1 if x0 < x1 else -1
    sy = 1 if y0 < y1 else -1
    err = dx + dy
    while True:
        if 0 <= x0 < w and 0 <= y0 < h:
            p = x0 * 3
            pixels[y0][p:p+3] = color
        if x0 == x1 and y0 == y1:
            break
        e2 = 2 * err
        if e2 >= dy:
            err += dy
            x0 += sx
        if e2 <= dx:
            err += dx
            y0 += sy

def plot(path, title, series):
    width, height = 1000, 650
    pixels = [bytearray([255, 255, 255] * width) for _ in range(height)]
    ml, mr, mt, mb = 80, 30, 40, 70
    x0, x1 = ml, width - mr
    y0, y1 = height - mb, mt
    for x in range(x0, x1 + 1):
        pixels[y0][3*x:3*x+3] = b"\x00\x00\x00"
        pixels[y1][3*x:3*x+3] = b"\x00\x00\x00"
    for y in range(y1, y0 + 1):
        pixels[y][3*x0:3*x0+3] = b"\x00\x00\x00"
        pixels[y][3*x1:3*x1+3] = b"\x00\x00\x00"
    points = [(x, y) for _, xs, ys in series for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if not points:
        save_png(path, width, height, pixels)
        return
    xmin, xmax = min(x for x, _ in points), max(x for x, _ in points)
    ymin, ymax = min(y for _, y in points), max(y for _, y in points)
    if xmax == xmin:
        xmax = xmin + 1.0
    if ymax == ymin:
        ymax = ymin + 1.0
    pad = 0.08 * (ymax - ymin)
    ymin -= pad
    ymax += pad
    colors = [(31,119,180), (214,39,40), (44,160,44)]
    for idx, (_, xs, ys) in enumerate(series):
        pts = []
        for x, y in zip(xs, ys):
            if math.isfinite(x) and math.isfinite(y):
                px = int(x0 + (x - xmin) / (xmax - xmin) * (x1 - x0))
                py = int(y0 - (y - ymin) / (ymax - ymin) * (y0 - y1))
                pts.append((px, py))
        for a, b in zip(pts, pts[1:]):
            draw_line(pixels, a[0], a[1], b[0], b[1], colors[idx % len(colors)])
    save_png(path, width, height, pixels)

def overall(results):
    q74 = results["qe74"]["classification"]
    q75 = results["qe75_af7ebfc"]["classification"]
    if q74 == "PHYSICAL_HEATING" and q75 == "PHYSICAL_HEATING":
        return "Both executables physically heat; CAIMD activation is not required to reproduce the problem."
    if q74 == "STABLE" and q75 == "PHYSICAL_HEATING":
        return "Official QE 7.4.1 is stable while af7ebfc heats; next test should be pristine unmodified QE 7.5."
    if q74 == "PHYSICAL_HEATING" and q75 == "STABLE":
        return "Official QE 7.4.1 heats while af7ebfc is stable; inspect version/runtime differences."
    if q74 == "STABLE" and q75 == "STABLE":
        return "Both are stable; compare Leonardo and local inputs/build/runtime details."
    if "DIAGNOSTIC_MISMATCH" in (q74, q75):
        return "A diagnostic mismatch exists; resolve temperature-output semantics first."
    return "At least one run is incomplete or numerically failed; inspect output before drawing a thermostat conclusion."

def main():
    (ROOT / "analysis").mkdir(exist_ok=True)
    rcs = return_codes()
    results = {}
    eal_by_key = {}
    det_by_key = {}
    for key, cfg in RUNS.items():
        result, eal, dets = analyze_run(key, cfg, rcs)
        results[key] = {k: v for k, v in result.items() if k != "temperatures_K"}
        eal_by_key[key] = eal
        det_by_key[key] = dets
        write_csv(ROOT / "analysis" / ("qe74_temperature.csv" if key == "qe74" else "qe75_af7ebfc_temperature.csv"), eal)
    summary = {
        "runs": results,
        "overall_interpretation": overall(results),
        "input_identity": json.loads((ROOT / "evidence" / "input_identity_check.json").read_text()) if (ROOT / "evidence" / "input_identity_check.json").exists() else None,
    }
    (ROOT / "analysis" / "comparison_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    plot(ROOT / "analysis" / "temperature_vs_time.png", "temperature", [(k, [r["time_ps"] for r in eal], [r["temperature_K"] for r in eal]) for k, eal in eal_by_key.items()])
    plot(ROOT / "analysis" / "ionic_kinetic_energy_vs_time.png", "eka", [(k, [r["time_ps"] for r in eal], [r["eka_Ry"] for r in eal]) for k, eal in eal_by_key.items()])
    plot(ROOT / "analysis" / "cell_determinant_vs_time.png", "det", [(k, [(i+1)*DT_FS/1000.0 for i in range(len(dets))], dets) for k, dets in det_by_key.items()])
    plot(ROOT / "analysis" / "volume_vs_time.png", "volume", [(k, [(i+1)*DT_FS/1000.0 for i in range(len(dets))], dets) for k, dets in det_by_key.items()])
    q74 = results["qe74"]
    q75 = results["qe75_af7ebfc"]
    report = [
        "# Final Comparison Report",
        "",
        f"Generated by Slurm job analysis under `{ROOT}`.",
        "",
        f"Overall interpretation: {summary['overall_interpretation']}",
        "",
        "## Run Classifications",
        "",
        f"- QE 7.4.1 official: {q74['classification']}; steps={q74['trajectory_records']}; max T={q74['maximum_temperature_K']}; terminal direct T={q74['terminal_direct_temperature_K']}",
        f"- QE 7.5 af7ebfc ordinary: {q75['classification']}; steps={q75['trajectory_records']}; max T={q75['maximum_temperature_K']}; terminal direct T={q75['terminal_direct_temperature_K']}",
        "",
        "## Required Answers",
        "",
        f"1. Did official QE 7.4.1 ordinary `vc-md + Nose` physically heat? {q74['classification'] == 'PHYSICAL_HEATING'}",
        f"2. Did af7ebfc ordinary `vc-md + Nose` physically heat with the feature absent from input? {q75['classification'] == 'PHYSICAL_HEATING'}",
        f"3. Did `.eal` and direct terminal velocities agree? QE74 diff K={q74['terminal_eal_direct_abs_diff_K']}; QE75 diff K={q75['terminal_eal_direct_abs_diff_K']}",
        f"4. Which executable crossed 600 K and 2000 K first? QE74 600/2000={q74['first_step_above_600K']}/{q74['first_step_above_2000K']}; QE75 600/2000={q75['first_step_above_600K']}/{q75['first_step_above_2000K']}",
        f"5. Were the input states physically identical? {summary['input_identity']['all_physical_differences_zero'] if summary['input_identity'] else None}",
        f"6. Is activation of the new feature required to reproduce the problem? {'No' if q74['classification'] == 'PHYSICAL_HEATING' and q75['classification'] == 'PHYSICAL_HEATING' else 'Not determined'}",
        "7. Minimum scientifically justified next test: pristine unmodified QE 7.5 with the same input if QE74 and af7ebfc differ; otherwise investigate ordinary Nose/PR parameters and runtime.",
    ]
    (ROOT / "FINAL_COMPARISON_REPORT.md").write_text("\n".join(report) + "\n")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
'''


SLURM_SCRIPT = r'''#!/usr/bin/env bash
#SBATCH --job-name=diamond8_nose_qe74_qe75
#SBATCH --account=ICT26_CMSP_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

WORKDIR="@WORKDIR@"
QE74_MODULE="@QE74_MODULE@"
QE75_ROOT="@QE75_ROOT@"
QE75_EXE="@QE75_EXE@"

cd "$WORKDIR"
mkdir -p evidence analysis qe74_official/tmp qe75_af7ebfc_no_caimd/tmp

{
  echo "job_id=${SLURM_JOB_ID:-unknown}"
  echo "job_name=${SLURM_JOB_NAME:-unknown}"
  echo "nodelist=${SLURM_JOB_NODELIST:-unknown}"
  echo "submit_dir=${SLURM_SUBMIT_DIR:-unknown}"
  echo "workdir=$WORKDIR"
  hostname
  date -u
  env | grep -E '^(SLURM|CUDA|NVHPC|OMPI|OMP)_' | sort || true
} > evidence/slurm_job_environment.txt 2>&1

sha256sum inputs/qe74.in inputs/qe75_af7ebfc.in pseudo/C.UPF > evidence/input_checksums.txt

module purge
module load profile/chem-phys
module load "$QE74_MODULE"
QE74_EXE="$(command -v pw.x)"
{
  echo "module=$QE74_MODULE"
  echo "executable=$QE74_EXE"
  echo "srun_command=srun --mpi=pmix_v3 -n 1 $QE74_EXE -in $WORKDIR/inputs/qe74.in"
  module list
  timeout 30 "$QE74_EXE" --version || true
  sha256sum "$QE74_EXE"
  ldd "$QE74_EXE" || true
} > evidence/qe74_environment.txt 2>&1

{
  echo "qe74_executable=$QE74_EXE"
  sha256sum "$QE74_EXE"
} > evidence/executable_checksums.txt

cd "$WORKDIR/qe74_official"
set +e
srun --mpi=pmix_v3 -n 1 "$QE74_EXE" -in "$WORKDIR/inputs/qe74.in" > qe.out 2> qe.err
RC_QE74=$?
set -e

cd "$WORKDIR"
module purge
if [ -f "$QE75_ROOT/scripts/load_build_environment.sh" ]; then
  # shellcheck source=/dev/null
  source "$QE75_ROOT/scripts/load_build_environment.sh"
else
  module load profile/chem-phys
  module load nvhpc/24.5
  module load hpcx-mpi/2.19
  module load fftw/3.3.10--hpcx-mpi--2.19--nvhpc--24.5
  module load libxc/7.0.0--nvhpc--24.5
  module load openblas/0.3.26--nvhpc--24.5
  module load netlib-scalapack/2.2.0--hpcx-mpi--2.19--nvhpc--24.5
fi
export LD_LIBRARY_PATH="$QE75_ROOT/local-libs/lib:${LD_LIBRARY_PATH:-}"
{
  echo "executable=$QE75_EXE"
  echo "loader=$QE75_ROOT/scripts/load_build_environment.sh"
  echo "srun_command=srun --mpi=pmix_v3 -n 1 $QE75_EXE -in $WORKDIR/inputs/qe75_af7ebfc.in"
  module list
  timeout 30 "$QE75_EXE" --version || true
  sha256sum "$QE75_EXE"
  ldd "$QE75_EXE" || true
} > evidence/qe75_environment.txt 2>&1
{
  echo "qe75_af7ebfc_executable=$QE75_EXE"
  sha256sum "$QE75_EXE"
} >> evidence/executable_checksums.txt

cd "$WORKDIR/qe75_af7ebfc_no_caimd"
set +e
srun --mpi=pmix_v3 -n 1 "$QE75_EXE" -in "$WORKDIR/inputs/qe75_af7ebfc.in" > qe.out 2> qe.err
RC_QE75_AF7EBFC=$?
set -e

cd "$WORKDIR"
{
  echo "RC_QE74=$RC_QE74"
  echo "RC_QE75_AF7EBFC=$RC_QE75_AF7EBFC"
} > evidence/run_return_codes.txt

python3 scripts/analyze_comparison.py > evidence/analysis_stdout.txt 2> evidence/analysis_stderr.txt

exit 0
'''


REMOTE_SUBMIT_SCRIPT = r'''#!/usr/bin/env bash
set -euo pipefail

WORKDIR="@WORKDIR@"
QE74_MODULE="@QE74_MODULE@"
QE75_EXE="@QE75_EXE@"
QE75_ROOT="@QE75_ROOT@"
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
root = pathlib.Path("@WORKDIR@")
record = {
    "status": "SUBMITTED",
    "job_id": job_id,
    "workdir": str(root),
    "slurm_script": str(root / "scripts" / "run_comparison.slurm"),
    "qe74_module": "@QE74_MODULE@",
    "qe74_executable": qe74_exe,
    "qe74_executable_sha256": (root / "evidence" / "qe74_executable_pre_submit.sha256").read_text().split()[0],
    "qe75_af7ebfc_executable": "@QE75_EXE@",
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
'''


def main() -> int:
    state = json.loads(COMMON_STATE.read_text())
    for sub in ("inputs", "pseudo", "qe74_official", "qe75_af7ebfc_no_caimd", "scripts", "analysis", "evidence"):
        (ROOT / sub).mkdir(parents=True, exist_ok=True)
    shutil.copy2(PSEUDO_SOURCE, ROOT / "pseudo" / "C.UPF")

    qe74_input = input_text("diamond8_qe74_official", state)
    qe75_input = input_text("diamond8_qe75_af7ebfc", state)
    (ROOT / "inputs" / "qe74.in").write_text(qe74_input)
    (ROOT / "inputs" / "qe75_af7ebfc.in").write_text(qe75_input)

    (ROOT / "scripts" / "pre_submit_validate.py").write_text(VALIDATE_SCRIPT)
    (ROOT / "scripts" / "analyze_comparison.py").write_text(ANALYZE_SCRIPT)
    (ROOT / "scripts" / "run_comparison.slurm").write_text(
        SLURM_SCRIPT
        .replace("@WORKDIR@", REMOTE_WORKDIR)
        .replace("@QE74_MODULE@", QE74_MODULE)
        .replace("@QE75_EXE@", QE75_EXE)
        .replace("@QE75_ROOT@", QE75_ROOT)
    )
    (ROOT / "scripts" / "remote_validate_and_submit.sh").write_text(
        REMOTE_SUBMIT_SCRIPT
        .replace("@WORKDIR@", REMOTE_WORKDIR)
        .replace("@QE74_MODULE@", QE74_MODULE)
        .replace("@QE75_EXE@", QE75_EXE)
        .replace("@QE75_ROOT@", QE75_ROOT)
    )
    for path in ROOT.glob("scripts/*"):
        if path.suffix in (".py", ".sh", ".slurm") or path.name.endswith(".slurm"):
            path.chmod(0o755)

    metadata = {
        "created_utc": utc_now(),
        "local_repo_root": str(REPO_ROOT),
        "local_branch": run_text(["git", "branch", "--show-current"]).strip(),
        "local_head": run_text(["git", "rev-parse", "--short", "HEAD"]).strip(),
        "local_head_full": run_text(["git", "rev-parse", "HEAD"]).strip(),
        "local_git_status": run_text(["git", "status", "--short", "--branch"]).splitlines(),
        "local_control_workflow": str(CONTROL_ROOT),
        "common_state_source": {
            "path": str(COMMON_STATE),
            "source_step": state["source_step"],
            "source_files": state["source_files"],
            "velocity_scaling_factor": state["temperature_scaling"]["scaling_factor"],
            "starting_temperature_K": state["temperature_from_scaled_velocities"]["temperature_K"],
            "ionic_degrees_of_freedom": ND_I,
            "units": state["units"],
        },
        "local_pw_provenance": {"path": str(LOCAL_PW), "sha256": sha256(LOCAL_PW)},
        "pseudo": {"source": str(PSEUDO_SOURCE), "bundle_path": str(ROOT / "pseudo" / "C.UPF"), "sha256": sha256(ROOT / "pseudo" / "C.UPF")},
        "remote_workdir": REMOTE_WORKDIR,
        "qe74_module": QE74_MODULE,
        "qe75_af7ebfc_executable": QE75_EXE,
        "parameters": {
            "nstep": NSTEP,
            "dt_qe": DT_QE,
            "dt_fs": DT_FS,
            "simulated_time_ps": NSTEP * DT_FS / 1000.0,
            "tempw_K": TEMPW,
            "fnosep": FNOSEP,
            "nhpcl": NHPCL,
            "wmass": WMASS,
            "ecutwfc": ECUTWFC,
            "ecutrho": ECUTRHO,
            "conv_thr": CONV_THR,
            "electron_maxstep": ELECTRON_MAXSTEP,
            "mixing_beta": MIXING_BETA,
        },
    }
    write_json(ROOT / "WORKFLOW_METADATA.json", metadata)
    (ROOT / "FINAL_COMPARISON_REPORT.md").write_text("PENDING: generated by the Slurm job after both calculations finish.\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
