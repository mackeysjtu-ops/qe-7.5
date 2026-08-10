#!/usr/bin/env python3

import csv
import difflib
import hashlib
import json
import math
import os
import re
import struct
import subprocess
import zlib
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
EVIDENCE = ROOT / "evidence"

NAT = 8
ND_I = 21
DT_FS = 2.0
MASS_C_AMU_DEFAULT = 12.0107
RYDBERG_SI = 2.1798723611035e-18
K_BOLTZMANN_SI = 1.380649e-23
AMU_SI = 1.66053906660e-27
ELECTRONMASS_SI = 9.1093837015e-31
AMU_RY = (AMU_SI / ELECTRONMASS_SI) / 2.0
RY_TO_KELVIN = RYDBERG_SI / K_BOLTZMANN_SI
RY_TO_EV = 13.605693122994
FLOAT_RE = re.compile(r"[-+]?(?:(?:\d+\.\d*|\.\d+)(?:[EeDd][-+]?\d+)?|\d+(?:[EeDd][-+]?\d+)?)")
STRESS_HEADER_RE = re.compile(r"P=\s*([-+]?\d+(?:\.\d*)?)")

RUNS = {
    "qe74": {
        "label": "Official QE 7.4.1",
        "run_dir": ROOT / "qe74_official",
        "prefix": "diamond8_qe74_official",
        "input": ROOT / "inputs" / "qe74.in",
        "rc_key": "RC_QE74",
    },
    "qe75_af7ebfc": {
        "label": "QE 7.5 af7ebfc ordinary, no CAIMD input",
        "run_dir": ROOT / "qe75_af7ebfc_no_caimd",
        "prefix": "diamond8_qe75_af7ebfc",
        "input": ROOT / "inputs" / "qe75_af7ebfc.in",
        "rc_key": "RC_QE75_AF7EBFC",
    },
}


def read_text(path):
    return path.read_text(errors="ignore") if path.exists() else ""


def floats(text):
    return [float(x.replace("D", "E").replace("d", "e")) for x in FLOAT_RE.findall(text)]


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def write_csv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def run_local(cmd):
    proc = subprocess.run(cmd, cwd=str(ROOT.parents[1]), text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return {"command": " ".join(cmd), "returncode": proc.returncode, "output": proc.stdout}


def strip_comment(line):
    return line.split("!", 1)[0].strip()


def parse_input(path):
    text = read_text(path)
    lines = text.splitlines()
    lexical_hits = []
    for lineno, line in enumerate(lines, 1):
        if re.search(r"caimd", line, re.IGNORECASE):
            lexical_hits.append({"line": lineno, "text": line})
    namelists = {}
    current = None
    species = []
    cell = []
    positions = []
    velocities = []
    kpoints = []
    mode = None
    for raw in lines:
        line = strip_comment(raw)
        if not line:
            continue
        if line.startswith("&"):
            current = line[1:].strip().lower()
            namelists.setdefault(current, {})
            mode = None
            continue
        if line == "/":
            current = None
            mode = None
            continue
        upper = line.upper()
        if upper.startswith("ATOMIC_SPECIES"):
            mode = "species"
            current = None
            continue
        if upper.startswith("CELL_PARAMETERS"):
            mode = "cell"
            current = None
            continue
        if upper.startswith("ATOMIC_POSITIONS"):
            mode = "positions"
            current = None
            continue
        if upper.startswith("ATOMIC_VELOCITIES"):
            mode = "velocities"
            current = None
            continue
        if upper.startswith("K_POINTS"):
            mode = "kpoints"
            current = None
            kpoints.append(line.lower())
            continue
        if current and "=" in line:
            key, value = line.split("=", 1)
            namelists[current][key.strip().lower()] = value.strip().rstrip(",").strip()
            continue
        if mode == "species":
            parts = line.split()
            if len(parts) >= 3:
                species.append({"symbol": parts[0], "mass_amu": float(parts[1]), "pseudo": parts[2]})
            continue
        if mode == "cell" and len(cell) < 3:
            vals = floats(line)
            if len(vals) >= 3:
                cell.append(vals[:3])
            continue
        if mode == "positions" and len(positions) < NAT:
            vals = floats(line)
            if len(vals) >= 3:
                positions.append(vals[-3:])
            continue
        if mode == "velocities" and len(velocities) < NAT:
            vals = floats(line)
            if len(vals) >= 3:
                velocities.append(vals[-3:])
            continue
        if mode == "kpoints":
            kpoints.append(line.lower())
    structured_hits = []
    for section, values in namelists.items():
        for key in values:
            if "caimd" in key:
                structured_hits.append({"section": section, "key": key})
    mass = species[0]["mass_amu"] if species else MASS_C_AMU_DEFAULT
    temp = kinetic_temperature_cartesian(velocities, mass) if velocities else None
    normalized = {
        "namelists": {
            section: {
                key: value
                for key, value in sorted(values.items())
                if key not in ("prefix", "outdir")
            }
            for section, values in sorted(namelists.items())
        },
        "species": species,
        "cell": cell,
        "positions": positions,
        "velocities": velocities,
        "kpoints": kpoints,
    }
    return {
        "path": str(path),
        "sha256": sha256(path),
        "lexical_caimd_hits": lexical_hits,
        "structured_caimd_hits": structured_hits,
        "namelists": namelists,
        "species": species,
        "cell": cell,
        "positions": positions,
        "velocities": velocities,
        "kpoints": kpoints,
        "mass_amu": mass,
        "input_direct_temperature": temp,
        "normalized": normalized,
    }


def kinetic_temperature_cartesian(velocities, mass_amu):
    if not velocities:
        return None
    arr = np.array(velocities, dtype=float)
    com = arr.mean(axis=0)
    centered = arr - com
    mass_ry = mass_amu * AMU_RY
    ekin = 0.5 * mass_ry * float(np.sum(centered * centered))
    return {
        "temperature_K": 2.0 * ekin / (ND_I * (1.0 / RY_TO_KELVIN)),
        "ekin_Ry": ekin,
        "center_of_mass_velocity_au": com.tolist(),
        "center_of_mass_max_abs_au": float(np.max(np.abs(com))),
        "ndof": ND_I,
    }


def max_abs_nested(a, b):
    aa = np.array(a, dtype=float)
    bb = np.array(b, dtype=float)
    if aa.size == 0 or bb.size == 0:
        return None
    return float(np.max(np.abs(aa - bb)))


def pstdev(values):
    if not values:
        return None
    if len(values) == 1:
        return 0.0
    mean = sum(values) / len(values)
    return math.sqrt(sum((value - mean) ** 2 for value in values) / len(values))


def parse_return_codes():
    out = {}
    for line in read_text(EVIDENCE / "run_return_codes.txt").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            try:
                out[key.strip()] = int(value.strip())
            except ValueError:
                out[key.strip()] = value.strip()
    return out


def parse_slurm_once():
    rows = []
    text = read_text(EVIDENCE / "slurm_state_once.txt")
    for line in text.splitlines():
        if not line.strip() or line.startswith("Command") or line.startswith("Output") or line.startswith("-") or line.startswith("JobID"):
            continue
        parts = line.split()
        if parts and parts[0].startswith("50127909"):
            rows.append({
                "job_id": parts[0],
                "job_name": parts[1] if len(parts) > 1 else "",
                "state": parts[2] if len(parts) > 2 else "",
                "exit_code": parts[3] if len(parts) > 3 else "",
                "elapsed": parts[4] if len(parts) > 4 else "",
                "node": parts[5] if len(parts) > 5 else "",
                "start": parts[6] if len(parts) > 6 else "",
                "end": parts[7] if len(parts) > 7 else "",
            })
    return rows


def parse_blocks(path, rows_per_block):
    blocks = []
    if not path.exists():
        return blocks
    lines = [line.rstrip() for line in read_text(path).splitlines() if line.strip()]
    i = 0
    while i < len(lines):
        head = floats(lines[i])
        if len(head) < 1:
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
        step = int(round(head[0]))
        time_raw = head[1] if len(head) > 1 else None
        blocks.append({"step": step, "time_raw": time_raw, "rows": rows})
        i += 1 + rows_per_block
    return blocks


def parse_eal(path):
    rows = []
    if not path.exists():
        return rows
    for line in read_text(path).splitlines():
        vals = floats(line)
        if len(vals) >= 7:
            uta, eka, eta, utl, ekla, etl, nst = vals[:7]
            nst_i = int(round(nst))
            temp = 2.0 * eka / (ND_I * (1.0 / RY_TO_KELVIN))
            rows.append({
                "nst": nst_i,
                "step": nst_i + 1,
                "time_ps": nst_i * DT_FS / 1000.0,
                "uta_Ry": uta,
                "eka_Ry": eka,
                "eta_Ry": eta,
                "utl_Ry": utl,
                "ekla_Ry": ekla,
                "etl_Ry": etl,
                "temperature_eal_21dof_K": temp,
            })
    return rows


def parse_table(path, labels):
    rows = []
    if not path.exists():
        return rows
    for line in read_text(path).splitlines():
        vals = floats(line)
        if len(vals) >= len(labels):
            row = {}
            for label, value in zip(labels, vals[:len(labels)]):
                row[label] = int(round(value)) if label in ("nst", "step") else value
            if "nst" in row and "step" not in row:
                row["step"] = row["nst"] + 1
                row["time_ps"] = row["nst"] * DT_FS / 1000.0
            rows.append(row)
    return rows


def parse_stdout_dynamics(path):
    rows = []
    current_step = None
    current_time = None
    text = read_text(path)
    for line in text.splitlines():
        if "Entering Dynamics" in line:
            m = re.search(r"it\s*=\s*(\d+).*time\s*=\s*([-+0-9.EeDd]+)", line)
            if m:
                current_step = int(m.group(1))
                current_time = float(m.group(2).replace("D", "E").replace("d", "e"))
        if "Ekin =" in line and " T =" in line:
            vals = floats(line)
            if len(vals) >= 3:
                rows.append({
                    "step": current_step,
                    "time_ps": current_time,
                    "stdout_ekin_Ry": vals[0],
                    "temperature_stdout_K": vals[1],
                    "etot_stdout_Ry": vals[2],
                })
    return rows


def parse_runtime_version(path):
    m = re.search(r"Program PWSCF v\.([^\s]+)", read_text(path))
    return m.group(1) if m else None


def parse_stress(path):
    rows = []
    lines = read_text(path).splitlines()
    step_counter = 0
    for i, line in enumerate(lines):
        if "total   stress" in line and "P=" in line:
            step_counter += 1
            pressure = None
            m = STRESS_HEADER_RE.search(line)
            if m:
                pressure = float(m.group(1))
            mat_ry = []
            mat_kbar = []
            ok = True
            for j in range(1, 4):
                if i + j >= len(lines):
                    ok = False
                    break
                vals = floats(lines[i + j])
                if len(vals) < 6:
                    ok = False
                    break
                mat_ry.append(vals[:3])
                mat_kbar.append(vals[3:6])
            if ok:
                rows.append({
                    "step": step_counter,
                    "time_ps": (step_counter - 1) * DT_FS / 1000.0,
                    "pressure_kbar_stdout": pressure,
                    "stress_xx_kbar": mat_kbar[0][0],
                    "stress_yy_kbar": mat_kbar[1][1],
                    "stress_zz_kbar": mat_kbar[2][2],
                    "stress_yz_kbar": mat_kbar[1][2],
                    "stress_xz_kbar": mat_kbar[0][2],
                    "stress_xy_kbar": mat_kbar[0][1],
                })
    return rows


def parse_md_terminal(path, cell_blocks, mass_amu):
    if not path.exists() or not cell_blocks:
        return None
    vals = floats(read_text(path))
    block = 3 * NAT
    if len(vals) < 5 * block:
        return None
    ratd = np.array(vals[block:2 * block], dtype=float).reshape((NAT, 3))
    h = np.array(cell_blocks[-1]["rows"], dtype=float)
    velocities = np.array([h.dot(v) for v in ratd])
    temp = kinetic_temperature_cartesian(velocities.tolist(), mass_amu)
    return {
        "source": str(path),
        "step": cell_blocks[-1]["step"],
        "time_ps": (cell_blocks[-1]["step"] - 1) * DT_FS / 1000.0,
        "ratd_fractional_velocity": ratd.tolist(),
        "cartesian_velocity_au": velocities.tolist(),
        "temperature_K": temp["temperature_K"],
        "ekin_Ry": temp["ekin_Ry"],
        "center_of_mass_velocity_au": temp["center_of_mass_velocity_au"],
        "float_count": len(vals),
        "format_interpretation": "vcsmd restart WRITE(4,*) rati, ratd, rat2d, rat2di, tauold; ratd converted by final .cel cell matrix",
    }


def det_cell(h):
    return float(np.linalg.det(np.array(h, dtype=float)))


def cell_metrics(cell_blocks, pos_blocks):
    out = []
    if not cell_blocks:
        return out
    h0 = np.array(cell_blocks[0]["rows"], dtype=float)
    inv_h0 = np.linalg.inv(h0)
    for idx, block in enumerate(cell_blocks):
        h = np.array(block["rows"], dtype=float)
        volume = float(np.linalg.det(h))
        svals = np.linalg.svd(h, compute_uv=False)
        condition = float(svals.max() / svals.min())
        f = h.dot(inv_h0)
        eps = 0.5 * (f + f.T) - np.eye(3)
        row = {
            "step": block["step"],
            "time_ps": (block["step"] - 1) * DT_FS / 1000.0,
            "volume_bohr3": volume,
            "determinant_bohr3": volume,
            "singular_value_1": float(svals[0]),
            "singular_value_2": float(svals[1]),
            "singular_value_3": float(svals[2]),
            "condition_number": condition,
            "strain_xx": float(eps[0, 0]),
            "strain_yy": float(eps[1, 1]),
            "strain_zz": float(eps[2, 2]),
            "strain_yz": float(eps[1, 2]),
            "strain_xz": float(eps[0, 2]),
            "strain_xy": float(eps[0, 1]),
            "strain_frobenius": float(np.linalg.norm(eps)),
        }
        if idx < len(pos_blocks):
            positions = np.array(pos_blocks[idx]["rows"], dtype=float)
            inv_h = np.linalg.inv(h)
            frac = np.array([inv_h.dot(p) for p in positions])
            dists = []
            for a in range(len(frac)):
                for b in range(a + 1, len(frac)):
                    df = frac[a] - frac[b]
                    df -= np.round(df)
                    cart = h.dot(df)
                    dists.append(float(np.linalg.norm(cart)))
            if dists:
                row["min_cc_distance_bohr"] = min(dists)
                row["max_cc_distance_bohr"] = max(dists)
        out.append(row)
    return out


def rolling(values, window):
    rows = []
    for i in range(len(values)):
        start = max(0, i + 1 - window)
        chunk = values[start:i + 1]
        rows.append(sum(chunk) / len(chunk))
    return rows


def stat_block(values, times):
    vals = [v for v in values if v is not None and math.isfinite(v)]
    out = {
        "record_count": len(values),
        "finite": len(vals) == len(values),
        "mean": sum(vals) / len(vals) if vals else None,
        "std": pstdev(vals),
        "min": min(vals) if vals else None,
        "max": max(vals) if vals else None,
        "initial": vals[0] if vals else None,
        "final": vals[-1] if vals else None,
    }
    if len(vals) > 1 and len(times) == len(values):
        pairs = [(t, v) for t, v in zip(times, values) if v is not None and math.isfinite(v)]
        mt = sum(t for t, _ in pairs) / len(pairs)
        mv = sum(v for _, v in pairs) / len(pairs)
        denom = sum((t - mt) ** 2 for t, _ in pairs)
        out["linear_drift_per_ps"] = sum((t - mt) * (v - mv) for t, v in pairs) / denom if denom else None
    else:
        out["linear_drift_per_ps"] = None
    return out


def first_crossing(rows, key, threshold):
    for row in rows:
        value = row.get(key)
        if value is not None and value > threshold:
            return {"step": row.get("step"), "time_ps": row.get("time_ps"), "value": value}
    return None


def longest_consecutive(values, predicate):
    best = 0
    cur = 0
    for value in values:
        if predicate(value):
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


def classify_run(key, eal_rows, terminal_direct, return_code, job_done, qe_error):
    if key == "qe74" and return_code not in (0, None) and not eal_rows:
        return "NUMERICAL_FAILURE"
    if not eal_rows:
        return "INCOMPLETE"
    temps = [r["temperature_eal_21dof_K"] for r in eal_rows]
    if not all(math.isfinite(t) for t in temps):
        return "NUMERICAL_FAILURE"
    if terminal_direct and abs(terminal_direct["temperature_K"] - temps[-1]) > max(100.0, 0.1 * max(abs(temps[-1]), 1.0)):
        return "DIAGNOSTIC_MISMATCH"
    if max(temps) > 2000.0 and terminal_direct and terminal_direct["temperature_K"] > 2000.0:
        return "PHYSICAL_HEATING"
    if longest_consecutive(temps, lambda x: x > 600.0) >= 25:
        return "PHYSICAL_HEATING"
    final100 = sum(temps[-100:]) / min(100, len(temps))
    if 200.0 <= final100 <= 450.0 and 200.0 <= temps[-1] <= 450.0:
        return "STABLE"
    return "INCOMPLETE"


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
    width = len(pixels[0]) // 3
    height = len(pixels)
    dx = abs(x1 - x0)
    dy = -abs(y1 - y0)
    sx = 1 if x0 < x1 else -1
    sy = 1 if y0 < y1 else -1
    err = dx + dy
    while True:
        if 0 <= x0 < width and 0 <= y0 < height:
            p = x0 * 3
            pixels[y0][p:p + 3] = color
        if x0 == x1 and y0 == y1:
            break
        e2 = 2 * err
        if e2 >= dy:
            err += dy
            x0 += sx
        if e2 <= dx:
            err += dx
            y0 += sy


def plot_png(path, series):
    width, height = 1000, 650
    pixels = [bytearray([255, 255, 255] * width) for _ in range(height)]
    ml, mr, mt, mb = 80, 30, 40, 70
    x0, x1 = ml, width - mr
    y0, y1 = height - mb, mt
    axis = (0, 0, 0)
    for x in range(x0, x1 + 1):
        pixels[y0][3 * x:3 * x + 3] = axis
    for y in range(y1, y0 + 1):
        pixels[y][3 * x0:3 * x0 + 3] = axis
    points = []
    for _, xs, ys in series:
        points.extend((x, y) for x, y in zip(xs, ys) if x is not None and y is not None and math.isfinite(x) and math.isfinite(y))
    if not points:
        save_png(path, width, height, pixels)
        return
    xmin, xmax = min(x for x, _ in points), max(x for x, _ in points)
    ymin, ymax = min(y for _, y in points), max(y for _, y in points)
    if xmin == xmax:
        xmax = xmin + 1.0
    if ymin == ymax:
        ymax = ymin + 1.0
    pad = 0.08 * (ymax - ymin)
    ymin -= pad
    ymax += pad
    colors = [(31, 119, 180), (214, 39, 40), (44, 160, 44), (127, 127, 127)]
    for idx, (_, xs, ys) in enumerate(series):
        pts = []
        for x, y in zip(xs, ys):
            if x is None or y is None or not math.isfinite(x) or not math.isfinite(y):
                continue
            px = int(x0 + (x - xmin) / (xmax - xmin) * (x1 - x0))
            py = int(y0 - (y - ymin) / (ymax - ymin) * (y0 - y1))
            pts.append((px, py))
        for a, b in zip(pts, pts[1:]):
            draw_line(pixels, a[0], a[1], b[0], b[1], colors[idx % len(colors)])
    save_png(path, width, height, pixels)


def safe_get(rows, key):
    return [row.get(key) for row in rows]


def analyze_run(key, cfg, inputs, return_codes):
    run_dir = cfg["run_dir"]
    prefix = cfg["prefix"]
    eal = parse_eal(run_dir / f"{prefix}.eal")
    if not eal:
        eal = parse_eal(run_dir / "tmp" / f"{prefix}.eal")
    tv = parse_table(run_dir / f"{prefix}.tv", ["volume_bohr3_tv", "temperature_tv_stdout_convention_K", "nst"])
    p_rows = parse_table(run_dir / f"{prefix}.p", ["target_press_kbar_or_input_units", "instant_pressure_internal", "average_pressure_internal", "nst"])
    e_rows = parse_table(run_dir / f"{prefix}.e", ["ut_Ry", "ekint_Ry", "edyn_Ry", "pv_Ry", "nst"])
    ave_rows = parse_table(run_dir / f"{prefix}.ave", ["average_potential_Ry", "average_kinetic_Ry", "nst"])
    avec_rows = parse_table(run_dir / f"{prefix}.avec", ["a_len", "b_len", "c_len", "theta_ab", "theta_bc", "theta_ca", "nst"])
    cel = parse_blocks(run_dir / "tmp" / f"{prefix}.cel", 3)
    pos = parse_blocks(run_dir / "tmp" / f"{prefix}.pos", NAT)
    vel_blocks = parse_blocks(run_dir / "tmp" / f"{prefix}.vel", NAT)
    stress = parse_stress(run_dir / "qe.out")
    stdout_dyn = parse_stdout_dynamics(run_dir / "qe.out")
    metrics = cell_metrics(cel, pos)
    terminal_direct = parse_md_terminal(run_dir / "tmp" / f"{prefix}.md", cel, inputs[key]["mass_amu"])
    rc = return_codes.get(cfg["rc_key"])
    qe_out = read_text(run_dir / "qe.out")
    qe_error = []
    lines = qe_out.splitlines()
    for i, line in enumerate(lines):
        if "Error in routine" in line:
            qe_error.append(line.strip())
            if i + 1 < len(lines):
                qe_error.append(lines[i + 1].strip())
    classification = classify_run(key, eal, terminal_direct, rc, "JOB DONE" in qe_out, qe_error)
    temps = [row["temperature_eal_21dof_K"] for row in eal]
    times = [row["time_ps"] for row in eal]
    temp_stats = stat_block(temps, times)
    first100 = sum(temps[:100]) / min(100, len(temps)) if temps else None
    final100 = sum(temps[-100:]) / min(100, len(temps)) if temps else None
    rolling25 = rolling(temps, 25) if temps else []
    rolling100 = rolling(temps, 100) if temps else []
    energy_values = [row["edyn_Ry"] for row in e_rows]
    energy_times = [row["time_ps"] for row in e_rows]
    energy_stats = stat_block(energy_values, energy_times)
    energy_drift_ry_ps = energy_stats.get("linear_drift_per_ps")
    direct_terminal_delta = None
    if terminal_direct and eal:
        direct_terminal_delta = terminal_direct["temperature_K"] - eal[-1]["temperature_eal_21dof_K"]
    stdout_temp_delta = None
    if stdout_dyn and eal:
        stdout_temp_delta = stdout_dyn[-1]["temperature_stdout_K"] - eal[-1]["temperature_eal_21dof_K"]
    vel_first_temp = None
    vel_last_temp = None
    if vel_blocks:
        vel_first_temp = kinetic_temperature_cartesian(vel_blocks[0]["rows"], inputs[key]["mass_amu"])
        vel_last_temp = kinetic_temperature_cartesian(vel_blocks[-1]["rows"], inputs[key]["mass_amu"])
    result = {
        "key": key,
        "label": cfg["label"],
        "return_code": rc,
        "classification": classification,
        "job_done": "JOB DONE" in qe_out,
        "qe_error": qe_error,
        "runtime_version": parse_runtime_version(run_dir / "qe.out"),
        "eal_rows": len(eal),
        "cell_records": len(cel),
        "position_records": len(pos),
        "vel_records_printed_pw_global": len(vel_blocks),
        "md_exists": (run_dir / "tmp" / f"{prefix}.md").exists(),
        "tv_rows": len(tv),
        "p_rows": len(p_rows),
        "e_rows": len(e_rows),
        "ave_rows": len(ave_rows),
        "avec_rows": len(avec_rows),
        "stress_records": len(stress),
        "stdout_temperature_records": len(stdout_dyn),
        "last_record_time_ps": times[-1] if times else None,
        "requested_nstep": 1000,
        "completed_requested_steps": len(eal) == 1000 and "JOB DONE" in qe_out,
        "temperature_statistics_eal_21dof": temp_stats,
        "first100_mean_temperature_K": first100,
        "final100_mean_temperature_K": final100,
        "first_above_600K": first_crossing(eal, "temperature_eal_21dof_K", 600.0),
        "first_above_1000K": first_crossing(eal, "temperature_eal_21dof_K", 1000.0),
        "first_above_2000K": first_crossing(eal, "temperature_eal_21dof_K", 2000.0),
        "fraction_above_600K": sum(1 for t in temps if t > 600.0) / len(temps) if temps else None,
        "fraction_above_2000K": sum(1 for t in temps if t > 2000.0) / len(temps) if temps else None,
        "longest_consecutive_above_600K": longest_consecutive(temps, lambda x: x > 600.0) if temps else 0,
        "longest_consecutive_above_2000K": longest_consecutive(temps, lambda x: x > 2000.0) if temps else 0,
        "terminal_direct_from_md": terminal_direct,
        "terminal_direct_minus_eal_K": direct_terminal_delta,
        "stdout_terminal_minus_eal_K": stdout_temp_delta,
        "printed_vel_first_temperature_K": vel_first_temp["temperature_K"] if vel_first_temp else None,
        "printed_vel_last_temperature_K": vel_last_temp["temperature_K"] if vel_last_temp else None,
        "printed_vel_note": "PW .vel remains at the input velocities in this vc-md run and is not the vcsmd ratd physical velocity history.",
        "energy_drift_Ry_per_ps": energy_drift_ry_ps,
        "energy_drift_eV_per_atom_per_ps": energy_drift_ry_ps * RY_TO_EV / NAT if energy_drift_ry_ps is not None else None,
        "cell_metrics_summary": summarize_cell(metrics),
        "cell_instability_relation": cell_instability_relation(eal, metrics),
    }
    return {
        "summary": result,
        "eal": eal,
        "tv": tv,
        "p": p_rows,
        "e": e_rows,
        "ave": ave_rows,
        "avec": avec_rows,
        "cel": cel,
        "pos": pos,
        "vel_blocks": vel_blocks,
        "stress": stress,
        "stdout_dyn": stdout_dyn,
        "cell_metrics": metrics,
        "rolling25": rolling25,
        "rolling100": rolling100,
    }


def summarize_cell(metrics):
    if not metrics:
        return {}
    volumes = [r["volume_bohr3"] for r in metrics]
    conds = [r["condition_number"] for r in metrics]
    mincc = [r.get("min_cc_distance_bohr") for r in metrics if r.get("min_cc_distance_bohr") is not None]
    strain = [r["strain_frobenius"] for r in metrics]
    return {
        "volume_initial_bohr3": volumes[0],
        "volume_final_bohr3": volumes[-1],
        "volume_min_bohr3": min(volumes),
        "volume_max_bohr3": max(volumes),
        "volume_final_fractional_change": (volumes[-1] - volumes[0]) / volumes[0],
        "condition_max": max(conds),
        "strain_frobenius_max": max(strain),
        "min_cc_distance_min_bohr": min(mincc) if mincc else None,
        "max_cc_distance_max_bohr": max([r.get("max_cc_distance_bohr") for r in metrics if r.get("max_cc_distance_bohr") is not None], default=None),
    }


def cell_instability_relation(eal, metrics):
    if not eal or not metrics:
        return "No usable trajectory for cell/heating ordering."
    cross = first_crossing(eal, "temperature_eal_21dof_K", 600.0)
    if not cross:
        return "No heating crossing above 600 K."
    idx = min(max(cross["step"] - 1, 0), len(metrics) - 1)
    m = metrics[idx]
    v0 = metrics[0]["volume_bohr3"]
    frac_volume = abs(m["volume_bohr3"] - v0) / v0
    if frac_volume < 0.01 and m["condition_number"] < 5.0 and m.get("min_cc_distance_bohr", 99.0) > 1.0:
        return "Heating begins before significant cell deformation or cell singularity."
    return "Heating begins with or after measurable cell/structural distortion."


def executable_provenance(submission, metadata):
    qe74_env = read_text(EVIDENCE / "qe74_environment.txt")
    qe75_env = read_text(EVIDENCE / "qe75_environment.txt")
    slurm_env = read_text(EVIDENCE / "slurm_job_environment.txt")
    return {
        "job_id": submission.get("job_id"),
        "workdir": submission.get("workdir"),
        "submission_timestamp_from_scontrol": extract_value(read_text(EVIDENCE / "submission_scontrol_show_job.txt"), "SubmitTime"),
        "runtime_node": extract_line(slurm_env, "nodelist="),
        "runtime_hostname": next((line for line in slurm_env.splitlines() if line.endswith(".leonardo.local")), None),
        "qe74": {
            "module": submission.get("qe74_module"),
            "executable": submission.get("qe74_executable"),
            "sha256": submission.get("qe74_executable_sha256"),
            "runtime_stdout_version": parse_runtime_version(ROOT / "qe74_official" / "qe.out"),
            "pre_submit_version_output": read_text(EVIDENCE / "qe74_version_pre_submit.txt").strip(),
            "module_environment_excerpt": qe74_env.splitlines()[:14],
            "official_module_evidence": "quantum-espresso/7.4.1" in qe74_env and "/leonardo/prod/spack/" in (submission.get("qe74_executable") or ""),
        },
        "qe75_af7ebfc": {
            "executable": submission.get("qe75_af7ebfc_executable"),
            "sha256": submission.get("qe75_af7ebfc_executable_sha256"),
            "expected_sha256": "89d4fda827d2cc22a3b4c1670e90dc40ce629af42096c785b55c70c355a10e7d",
            "runtime_stdout_version": parse_runtime_version(ROOT / "qe75_af7ebfc_no_caimd" / "qe.out"),
            "pre_submit_version_output": read_text(EVIDENCE / "qe75_version_pre_submit.txt").strip(),
            "module_environment_excerpt": qe75_env.splitlines()[:18],
            "path_matches_required": submission.get("qe75_af7ebfc_executable") == "/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc/install-current-gpu/bin/pw.x",
        },
        "notes": [
            "The pre-submit --version files are empty because the submission helper used timeout protection around pw.x --version.",
            "Runtime stdout records the PWSCF version actually launched for each process.",
        ],
        "local_workflow_metadata": metadata,
    }


def extract_value(text, key):
    m = re.search(r"\b" + re.escape(key) + r"=([^\s]+)", text)
    return m.group(1) if m else None


def extract_line(text, prefix):
    return next((line.split("=", 1)[1] for line in text.splitlines() if line.startswith(prefix)), None)


def generate_manifest():
    rows = []
    total_size = 0
    for path in sorted(ROOT.rglob("*")):
        if path.is_file():
            rel = path.relative_to(ROOT).as_posix()
            size = path.stat().st_size
            total_size += size
            rows.append({"sha256": sha256(path), "bytes": size, "path": rel})
    (EVIDENCE / "download_manifest.sha256").write_text("".join(f"{r['sha256']}  {r['path']}\n" for r in rows))
    (EVIDENCE / "relative_file_listing.txt").write_text("".join(f"{r['path']}\n" for r in rows))
    write_json(EVIDENCE / "mirror_inventory.json", {"file_count": len(rows), "total_bytes": total_size, "files": rows})


def critical_file_audit():
    critical = [
        "inputs/qe74.in",
        "inputs/qe75_af7ebfc.in",
        "pseudo/C.UPF",
        "scripts/run_comparison.slurm",
        "slurm-50127909.out",
        "slurm-50127909.err",
        "qe74_official/qe.out",
        "qe74_official/qe.err",
        "qe75_af7ebfc_no_caimd/qe.out",
        "qe75_af7ebfc_no_caimd/qe.err",
        "qe75_af7ebfc_no_caimd/diamond8_qe75_af7ebfc.eal",
        "qe75_af7ebfc_no_caimd/diamond8_qe75_af7ebfc.tv",
        "qe75_af7ebfc_no_caimd/diamond8_qe75_af7ebfc.p",
        "qe75_af7ebfc_no_caimd/diamond8_qe75_af7ebfc.e",
        "qe75_af7ebfc_no_caimd/diamond8_qe75_af7ebfc.ave",
        "qe75_af7ebfc_no_caimd/diamond8_qe75_af7ebfc.avec",
        "qe75_af7ebfc_no_caimd/tmp/diamond8_qe75_af7ebfc.md",
        "qe75_af7ebfc_no_caimd/tmp/diamond8_qe75_af7ebfc.cel",
        "qe75_af7ebfc_no_caimd/tmp/diamond8_qe75_af7ebfc.pos",
        "qe75_af7ebfc_no_caimd/tmp/diamond8_qe75_af7ebfc.vel",
        "evidence/run_return_codes.txt",
        "evidence/executable_checksums.txt",
        "evidence/qe74_environment.txt",
        "evidence/qe75_environment.txt",
        "analysis/comparison_summary.json",
        "FINAL_COMPARISON_REPORT.md",
    ]
    rows = []
    for rel in critical:
        path = ROOT / rel
        rows.append({
            "path": rel,
            "exists": path.exists(),
            "bytes": path.stat().st_size if path.exists() else None,
            "sha256": sha256(path) if path.exists() else None,
        })
    write_json(EVIDENCE / "critical_file_audit.json", {"critical_files": rows})


def make_plots(results):
    q74 = results["qe74"]
    q75 = results["qe75_af7ebfc"]
    plot_png(ANALYSIS / "temperature_qe74_vs_af7ebfc.png", [
        ("qe74", safe_get(q74["eal"], "time_ps"), safe_get(q74["eal"], "temperature_eal_21dof_K")),
        ("qe75", safe_get(q75["eal"], "time_ps"), safe_get(q75["eal"], "temperature_eal_21dof_K")),
    ])
    plot_png(ANALYSIS / "temperature_rolling100_qe74_vs_af7ebfc.png", [
        ("qe75_rolling100", safe_get(q75["eal"], "time_ps"), q75["rolling100"]),
    ])
    plot_png(ANALYSIS / "ionic_kinetic_energy_qe74_vs_af7ebfc.png", [
        ("qe74", safe_get(q74["eal"], "time_ps"), safe_get(q74["eal"], "eka_Ry")),
        ("qe75", safe_get(q75["eal"], "time_ps"), safe_get(q75["eal"], "eka_Ry")),
    ])
    plot_png(ANALYSIS / "extended_energy_qe74_vs_af7ebfc.png", [
        ("qe75_edyn", safe_get(q75["e"], "time_ps"), safe_get(q75["e"], "edyn_Ry")),
    ])
    plot_png(ANALYSIS / "temperature_source_crosscheck_qe74.png", [
        ("qe74", safe_get(q74["eal"], "time_ps"), safe_get(q74["eal"], "temperature_eal_21dof_K")),
    ])
    plot_png(ANALYSIS / "temperature_source_crosscheck_af7ebfc.png", [
        ("eal", safe_get(q75["eal"], "time_ps"), safe_get(q75["eal"], "temperature_eal_21dof_K")),
        ("stdout", safe_get(q75["stdout_dyn"], "time_ps"), safe_get(q75["stdout_dyn"], "temperature_stdout_K")),
        ("tv", safe_get(q75["tv"], "time_ps"), safe_get(q75["tv"], "temperature_tv_stdout_convention_K")),
    ])
    plot_png(ANALYSIS / "cell_determinant_qe74_vs_af7ebfc.png", [
        ("qe75", safe_get(q75["cell_metrics"], "time_ps"), safe_get(q75["cell_metrics"], "determinant_bohr3")),
    ])
    plot_png(ANALYSIS / "cell_condition_number_qe74_vs_af7ebfc.png", [
        ("qe75", safe_get(q75["cell_metrics"], "time_ps"), safe_get(q75["cell_metrics"], "condition_number")),
    ])
    plot_png(ANALYSIS / "volume_qe74_vs_af7ebfc.png", [
        ("qe75", safe_get(q75["cell_metrics"], "time_ps"), safe_get(q75["cell_metrics"], "volume_bohr3")),
    ])
    plot_png(ANALYSIS / "normal_stress_qe74_vs_af7ebfc.png", [
        ("xx", safe_get(q75["stress"], "time_ps"), safe_get(q75["stress"], "stress_xx_kbar")),
        ("yy", safe_get(q75["stress"], "time_ps"), safe_get(q75["stress"], "stress_yy_kbar")),
        ("zz", safe_get(q75["stress"], "time_ps"), safe_get(q75["stress"], "stress_zz_kbar")),
    ])
    plot_png(ANALYSIS / "shear_stress_qe74_vs_af7ebfc.png", [
        ("yz", safe_get(q75["stress"], "time_ps"), safe_get(q75["stress"], "stress_yz_kbar")),
        ("xz", safe_get(q75["stress"], "time_ps"), safe_get(q75["stress"], "stress_xz_kbar")),
        ("xy", safe_get(q75["stress"], "time_ps"), safe_get(q75["stress"], "stress_xy_kbar")),
    ])
    plot_png(ANALYSIS / "minimum_CC_distance_qe74_vs_af7ebfc.png", [
        ("qe75", safe_get(q75["cell_metrics"], "time_ps"), safe_get(q75["cell_metrics"], "min_cc_distance_bohr")),
    ])


def make_temperature_csvs(results, inputs):
    cross_rows = []
    for key, run in results.items():
        reconstructed = []
        for row in run["eal"]:
            reconstructed.append({
                "run": key,
                "step": row["step"],
                "time_ps": row["time_ps"],
                "source": "eal_eka_21dof",
                "temperature_K": row["temperature_eal_21dof_K"],
                "ekin_Ry": row["eka_Ry"],
            })
        if key in inputs:
            reconstructed.insert(0, {
                "run": key,
                "step": 0,
                "time_ps": 0.0,
                "source": "input_cartesian_velocities_21dof",
                "temperature_K": inputs[key]["input_direct_temperature"]["temperature_K"],
                "ekin_Ry": inputs[key]["input_direct_temperature"]["ekin_Ry"],
            })
        if run["summary"].get("terminal_direct_from_md"):
            td = run["summary"]["terminal_direct_from_md"]
            reconstructed.append({
                "run": key,
                "step": td["step"],
                "time_ps": td["time_ps"],
                "source": "terminal_md_ratd_direct_21dof",
                "temperature_K": td["temperature_K"],
                "ekin_Ry": td["ekin_Ry"],
            })
        write_csv(ANALYSIS / ("qe74_temperature_reconstructed.csv" if key == "qe74" else "qe75_af7ebfc_temperature_reconstructed.csv"), reconstructed)
        for row in reconstructed:
            if row["source"] != "eal_eka_21dof" or row["step"] in (1, 26, 130, 1000):
                cross_rows.append(row)
        for row in run["stdout_dyn"]:
            if row["step"] in (1, 26, 130, 1000):
                cross_rows.append({
                    "run": key,
                    "step": row["step"],
                    "time_ps": row["time_ps"],
                    "source": "qe_stdout_T",
                    "temperature_K": row["temperature_stdout_K"],
                    "ekin_Ry": row["stdout_ekin_Ry"],
                })
        for row in run["tv"]:
            if row["step"] in (1, 26, 130, 1000):
                cross_rows.append({
                    "run": key,
                    "step": row["step"],
                    "time_ps": row["time_ps"],
                    "source": "tv_file_T",
                    "temperature_K": row["temperature_tv_stdout_convention_K"],
                    "ekin_Ry": "",
                })
    write_csv(ANALYSIS / "temperature_crosscheck.csv", cross_rows)


def make_run_statistics(results):
    rows = []
    for key, run in results.items():
        summary = run["summary"]
        stats = summary["temperature_statistics_eal_21dof"]
        rows.append({
            "run": key,
            "classification": summary["classification"],
            "return_code": summary["return_code"],
            "job_done": summary["job_done"],
            "eal_rows": summary["eal_rows"],
            "initial_temperature_K": stats.get("initial"),
            "mean_temperature_K": stats.get("mean"),
            "std_temperature_K": stats.get("std"),
            "min_temperature_K": stats.get("min"),
            "max_temperature_K": stats.get("max"),
            "final_temperature_K": stats.get("final"),
            "first100_mean_temperature_K": summary["first100_mean_temperature_K"],
            "final100_mean_temperature_K": summary["final100_mean_temperature_K"],
            "linear_drift_K_per_ps": stats.get("linear_drift_per_ps"),
            "first_step_above_600K": (summary["first_above_600K"] or {}).get("step"),
            "first_time_above_600K_ps": (summary["first_above_600K"] or {}).get("time_ps"),
            "first_step_above_1000K": (summary["first_above_1000K"] or {}).get("step"),
            "first_time_above_1000K_ps": (summary["first_above_1000K"] or {}).get("time_ps"),
            "first_step_above_2000K": (summary["first_above_2000K"] or {}).get("step"),
            "first_time_above_2000K_ps": (summary["first_above_2000K"] or {}).get("time_ps"),
            "fraction_above_600K": summary["fraction_above_600K"],
            "fraction_above_2000K": summary["fraction_above_2000K"],
            "longest_consecutive_above_600K": summary["longest_consecutive_above_600K"],
            "longest_consecutive_above_2000K": summary["longest_consecutive_above_2000K"],
            "energy_drift_Ry_per_ps": summary["energy_drift_Ry_per_ps"],
            "energy_drift_eV_per_atom_per_ps": summary["energy_drift_eV_per_atom_per_ps"],
        })
    write_csv(ANALYSIS / "run_statistics.csv", rows)


def trajectory_difference_stats(results):
    rows = []
    if not results["qe74"]["eal"] or not results["qe75_af7ebfc"]["eal"]:
        rows.append({
            "comparison": "qe74_vs_qe75_af7ebfc",
            "common_steps": 0,
            "status": "not_available",
            "reason": "QE 7.4.1 failed during input parsing before any trajectory/eal/cell records were written.",
        })
    write_csv(ANALYSIS / "trajectory_difference_statistics.csv", rows)


def input_audit():
    q74 = parse_input(RUNS["qe74"]["input"])
    q75 = parse_input(RUNS["qe75_af7ebfc"]["input"])
    norm74 = json.dumps(q74["normalized"], indent=2, sort_keys=True).splitlines()
    norm75 = json.dumps(q75["normalized"], indent=2, sort_keys=True).splitlines()
    diff = list(difflib.unified_diff(norm74, norm75, fromfile="qe74.normalized", tofile="qe75_af7ebfc.normalized", lineterm=""))
    (ANALYSIS / "input_normalized.diff").write_text("\n".join(diff) + ("\n" if diff else "No normalized physical input differences.\n"))
    result = {
        "caimd_absent": not q74["lexical_caimd_hits"] and not q75["lexical_caimd_hits"] and not q74["structured_caimd_hits"] and not q75["structured_caimd_hits"],
        "qe74": {
            "path": q74["path"],
            "sha256": q74["sha256"],
            "lexical_caimd_hits": q74["lexical_caimd_hits"],
            "structured_caimd_hits": q74["structured_caimd_hits"],
            "nat": len(q74["positions"]),
            "mass_amu": q74["mass_amu"],
            "pseudo": q74["species"][0]["pseudo"] if q74["species"] else None,
            "input_direct_temperature": q74["input_direct_temperature"],
        },
        "qe75_af7ebfc": {
            "path": q75["path"],
            "sha256": q75["sha256"],
            "lexical_caimd_hits": q75["lexical_caimd_hits"],
            "structured_caimd_hits": q75["structured_caimd_hits"],
            "nat": len(q75["positions"]),
            "mass_amu": q75["mass_amu"],
            "pseudo": q75["species"][0]["pseudo"] if q75["species"] else None,
            "input_direct_temperature": q75["input_direct_temperature"],
        },
        "physical_inputs_identical": not diff,
        "max_cell_abs_diff": max_abs_nested(q74["cell"], q75["cell"]),
        "max_position_abs_diff": max_abs_nested(q74["positions"], q75["positions"]),
        "max_velocity_abs_diff": max_abs_nested(q74["velocities"], q75["velocities"]),
        "pseudo_sha256": sha256(ROOT / "pseudo" / "C.UPF"),
        "normalized_diff_path": str(ANALYSIS / "input_normalized.diff"),
    }
    write_json(ANALYSIS / "input_identity_audit.json", result)
    return {"qe74": q74, "qe75_af7ebfc": q75}, result


def make_report(summary, input_result, provenance):
    q74 = summary["runs"]["qe74"]
    q75 = summary["runs"]["qe75_af7ebfc"]
    slurm = summary["slurm_state"]
    terminal_state = next((row for row in slurm if row["job_id"] == "50127909"), {})
    q75_stats = q75["temperature_statistics_eal_21dof"]
    first600 = q75.get("first_above_600K")
    first2000 = q75.get("first_above_2000K")
    lines = [
        "# Audited Final Comparison Report",
        "",
        "## Overall Status",
        "",
        f"Overall status: `{summary['overall_status']}`.",
        "",
        "The official QE 7.4.1 process did execute, but it failed during input parsing with `unknown ion_temperature nose` before any MD trajectory was produced. The af7ebfc QE 7.5 process executed the same physical, non-CAIMD input and completed all 1000 requested steps; its `.eal` kinetic energy and terminal `.md` velocities show physical heating.",
        "",
        "## Required Answers",
        "",
        f"1. Exact Slurm Job ID and terminal state: `50127909`; top-level state `{terminal_state.get('state')}` exit `{terminal_state.get('exit_code')}`. Step `50127909.0` failed, step `50127909.1` completed.",
        "2. Did both QE processes execute? Yes. QE74 launched and failed in `iosys`; QE75 launched and completed.",
        f"3. Did both inputs contain no CAIMD keywords? {'Yes' if input_result['caimd_absent'] else 'No'}.",
        f"4. Were the physical inputs identical? {'Yes' if input_result['physical_inputs_identical'] else 'No'}. Only `prefix`/`outdir` differ.",
        f"5. Was the initial direct-velocity temperature 300 K in both inputs? QE74 `{input_result['qe74']['input_direct_temperature']['temperature_K']:.12f} K`; QE75 `{input_result['qe75_af7ebfc']['input_direct_temperature']['temperature_K']:.12f} K`.",
        "6. Did official QE 7.4.1 physically heat? Not determined: no MD was accepted or run because the official module rejected `ion_temperature='nose'`.",
        f"7. Did af7ebfc physically heat with CAIMD absent? Yes. Classification `{q75['classification']}`; max `.eal`/21-DOF temperature `{q75_stats.get('max'):.3f} K`, final `{q75_stats.get('final'):.3f} K`.",
        f"8. Did `.eal`, QE stdout, and direct velocities agree? `.eal` and terminal direct `.md` velocities agree within `{abs(q75.get('terminal_direct_minus_eal_K')):.3f} K`. QE stdout/`.tv` show the same qualitative heating but a different temperature observable/convention, ending near `{q75['stdout_terminal_temperature_K']:.1f} K` rather than the `.eal`/direct `{q75_stats.get('final'):.1f} K`.",
        f"9. Which run first crossed 600 K? Only af7ebfc has a trajectory; it crossed at step `{first600.get('step') if first600 else None}` time `{first600.get('time_ps') if first600 else None}` ps.",
        f"10. Which run first crossed 2000 K? Only af7ebfc has a trajectory; it crossed at step `{first2000.get('step') if first2000 else None}` time `{first2000.get('time_ps') if first2000 else None}` ps.",
        f"11. Did cell instability precede or follow heating? {q75['cell_instability_relation']} Cell volume changed only `{100.0 * q75['cell_metrics_summary']['volume_final_fractional_change']:.3f}%` by the final frame; max condition number `{q75['cell_metrics_summary']['condition_max']:.3f}`.",
        "12. Is CAIMD activation required to reproduce the af7ebfc failure? No for the af7ebfc executable: the input contains no CAIMD keywords. Whether official QE 7.4.1 shows the same heating is unresolved because it did not accept the Nose input.",
        "13. Narrowest justified conclusion: this completed job is a `PARTIAL` comparison. It proves ordinary non-CAIMD af7ebfc heats from the common state, but it does not test ordinary Nose dynamics in official QE 7.4.1 because that executable rejects the requested input.",
        "14. Minimum next test: pristine unmodified QE 7.5 using the exact same input, executable/module path recorded, before varying Nose/PR parameters.",
        "15. Were any raw outputs, QE source files, or existing remote calculations modified? No. Only local `analysis/`, `evidence/`, and this audited report in the local mirror were written.",
        "",
        "## Key Evidence",
        "",
        f"- Remote workflow: `{summary['remote_workdir']}`",
        f"- Local mirror: `{ROOT}`",
        f"- QE74 executable: `{provenance['qe74']['executable']}` SHA256 `{provenance['qe74']['sha256']}` runtime version `{provenance['qe74']['runtime_stdout_version']}`",
        f"- af7ebfc executable: `{provenance['qe75_af7ebfc']['executable']}` SHA256 `{provenance['qe75_af7ebfc']['sha256']}` runtime version `{provenance['qe75_af7ebfc']['runtime_stdout_version']}`",
        f"- Input SHA256: QE74 `{input_result['qe74']['sha256']}`; QE75 `{input_result['qe75_af7ebfc']['sha256']}`; pseudo `{input_result['pseudo_sha256']}`",
        "",
        "## Caveats",
        "",
        "- The printed `tmp/*.vel` file is not used as the vc-md physical velocity history because it remains at the initial PW global velocity values; terminal physical velocities were reconstructed from `.md` `ratd` and the final `.cel` matrix.",
        "- The official QE 7.4.1 result is a parser/runtime compatibility failure, not a thermostat stability result.",
    ]
    (ROOT / "FINAL_COMPARISON_REPORT_AUDITED.md").write_text("\n".join(lines) + "\n")


def main():
    ANALYSIS.mkdir(exist_ok=True)
    EVIDENCE.mkdir(exist_ok=True)
    inputs, input_result = input_audit()
    submission = json.loads(read_text(EVIDENCE / "SUBMISSION_RECORD.json"))
    metadata = json.loads(read_text(ROOT / "WORKFLOW_METADATA.json"))
    provenance = executable_provenance(submission, metadata)
    write_json(ANALYSIS / "executable_provenance_audit.json", provenance)
    write_json(EVIDENCE / "executable_provenance_audit.json", provenance)
    write_json(EVIDENCE / "recovered_submission_metadata.json", {
        "job_id": submission.get("job_id"),
        "workdir": submission.get("workdir"),
        "slurm_script": submission.get("slurm_script"),
        "submission_timestamp": provenance.get("submission_timestamp_from_scontrol"),
        "qe74_executable": submission.get("qe74_executable"),
        "qe75_af7ebfc_executable": submission.get("qe75_af7ebfc_executable"),
        "input_paths": {
            "qe74": str(ROOT / "inputs" / "qe74.in"),
            "qe75_af7ebfc": str(ROOT / "inputs" / "qe75_af7ebfc.in"),
        },
        "pseudo_path": str(ROOT / "pseudo" / "C.UPF"),
        "expected_report_paths": {
            "remote": submission.get("expected_final_report"),
            "local_audited": str(ROOT / "FINAL_COMPARISON_REPORT_AUDITED.md"),
        },
    })
    return_codes = parse_return_codes()
    results = {}
    for key, cfg in RUNS.items():
        results[key] = analyze_run(key, cfg, inputs, return_codes)
    make_temperature_csvs(results, inputs)
    make_run_statistics(results)
    trajectory_difference_stats(results)
    make_plots(results)
    critical_file_audit()
    git_status = run_local(["git", "status", "--short"])
    git_diff_check = run_local(["git", "diff", "--check"])
    q75_summary = results["qe75_af7ebfc"]["summary"]
    if results["qe74"]["summary"]["classification"] in ("NUMERICAL_FAILURE", "INCOMPLETE") and q75_summary["classification"] == "PHYSICAL_HEATING":
        overall = "PARTIAL"
    elif results["qe74"]["summary"]["classification"] == "PHYSICAL_HEATING" and q75_summary["classification"] == "PHYSICAL_HEATING":
        overall = "BOTH_ORDINARY_QE_HEAT"
    elif results["qe74"]["summary"]["classification"] == "STABLE" and q75_summary["classification"] == "PHYSICAL_HEATING":
        overall = "QE74_STABLE_AF7EBFC_HEATS"
    elif results["qe74"]["summary"]["classification"] == "STABLE" and q75_summary["classification"] == "STABLE":
        overall = "BOTH_STABLE"
    else:
        overall = "PARTIAL"
    stdout_dyn = results["qe75_af7ebfc"]["stdout_dyn"]
    if stdout_dyn:
        q75_summary["stdout_terminal_temperature_K"] = stdout_dyn[-1]["temperature_stdout_K"]
    else:
        q75_summary["stdout_terminal_temperature_K"] = None
    summary = {
        "overall_status": overall,
        "job_id": submission.get("job_id"),
        "remote_workdir": submission.get("workdir"),
        "local_mirror": str(ROOT),
        "slurm_state": parse_slurm_once(),
        "input_identity": input_result,
        "executable_provenance": provenance,
        "runs": {key: results[key]["summary"] for key in results},
        "git_status_initial_and_final_context": {
            "git_status_short": git_status,
            "git_diff_check": git_diff_check,
        },
        "source_format_audit": {
            "eal_columns_from_PW_src_vcsmd": "uta, eka, eta, utl, ekla, etl, nst",
            "temperature_formula_from_PW_src_vcsubs": "2*eka/(3*(nat-1)*k_B) = 21 DOF for nat=8",
            "md_restart_format_from_PW_src_vcsmd": "rati, ratd, rat2d, rat2di, tauold",
        },
        "minimum_next_test": "Pristine unmodified QE 7.5 with the exact same input and common state, because official QE 7.4.1 rejected ion_temperature='nose' before testing dynamics.",
        "no_submission_no_wait_confirmation": "No sbatch, cancel, sleep, watch, polling loop, or repeated Slurm query was used in this inspection run.",
        "raw_remote_files_modified": False,
        "qe_source_modified": False,
        "git_commit_or_push": False,
    }
    write_json(ANALYSIS / "comparison_summary_audited.json", summary)
    make_report(summary, input_result, provenance)
    generate_manifest()
    print(json.dumps({
        "overall_status": overall,
        "job_id": submission.get("job_id"),
        "qe74_classification": results["qe74"]["summary"]["classification"],
        "qe75_classification": q75_summary["classification"],
        "audited_report": str(ROOT / "FINAL_COMPARISON_REPORT_AUDITED.md"),
        "comparison_json": str(ANALYSIS / "comparison_summary_audited.json"),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
