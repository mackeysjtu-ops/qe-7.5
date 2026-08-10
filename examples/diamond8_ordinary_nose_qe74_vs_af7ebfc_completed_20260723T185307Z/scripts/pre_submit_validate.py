#!/usr/bin/env python3

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
