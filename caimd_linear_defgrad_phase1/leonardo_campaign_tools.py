#!/usr/bin/env python3
"""Runtime tools for the Leonardo Diamond64 linear-defgrad campaign."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

CAMPAIGN_ROOT = Path(os.environ.get("CAMPAIGN_ROOT", "__CAMPAIGN_ROOT__"))
SOFTWARE_ROOT = Path(os.environ.get("SOFTWARE_ROOT", "__SOFTWARE_ROOT__"))
PW_EXE = SOFTWARE_ROOT / "install-current-gpu/bin/pw.x"
OLD_PW_EXE = Path("/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc/install-current-gpu/bin/pw.x")
PSEUDO = CAMPAIGN_ROOT / "pseudo/C.UPF"
ORIENTED = CAMPAIGN_ROOT / "source_structure/diamond64_oriented_111.json"

BOHR_ANG = 0.529177210903
RY_TO_EV = 13.605693122994
K_BOLTZMANN_RY = 6.333630300744e-6
DT_QE_2FS = 41.3413733366

EXPECTED = {
    "pw_sha": "__PW_SHA__",
    "old_pw_sha": "89d4fda827d2cc22a3b4c1670e90dc40ce629af42096c785b55c70c355a10e7d",
    "pseudo_sha": "9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639",
    "source_patch_sha": "ac5006766f0196d8750363bd6d802ec74f29b5ecc161ea143f8951172777fd9b",
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


def inv3(m):
    d = det3(m)
    if abs(d) < 1.0e-30:
        raise ValueError("singular 3x3 matrix")
    c = [
        [
            m[1][1] * m[2][2] - m[1][2] * m[2][1],
            -(m[0][1] * m[2][2] - m[0][2] * m[2][1]),
            m[0][1] * m[1][2] - m[0][2] * m[1][1],
        ],
        [
            -(m[1][0] * m[2][2] - m[1][2] * m[2][0]),
            m[0][0] * m[2][2] - m[0][2] * m[2][0],
            -(m[0][0] * m[1][2] - m[0][1] * m[1][0]),
        ],
        [
            m[1][0] * m[2][1] - m[1][1] * m[2][0],
            -(m[0][0] * m[2][1] - m[0][1] * m[2][0]),
            m[0][0] * m[1][1] - m[0][1] * m[1][0],
        ],
    ]
    return [[c[i][j] / d for j in range(3)] for i in range(3)]


def matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def matvec_row(fr, h):
    return [sum(fr[i] * h[i][j] for i in range(3)) for j in range(3)]


def eye():
    return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]


def fload(event_index, stride, dmat):
    f = eye()
    for i in range(3):
        for j in range(3):
            f[i][j] += event_index * stride * dmat[i][j]
    return f


def flatten(m):
    return [m[i][j] for i in range(3) for j in range(3)]


def read_structure(path: Path | str):
    return json.loads(Path(path).read_text())


def min_cc_ang(struct):
    h = struct["cell_bohr_rows"]
    pos = struct["fractional_positions"]
    best = 1.0e99
    for i in range(len(pos)):
        for j in range(i + 1, len(pos)):
            fi = pos[i][1]
            fj = pos[j][1]
            for sx in (-1, 0, 1):
                for sy in (-1, 0, 1):
                    for sz in (-1, 0, 1):
                        df = [fj[0] - fi[0] + sx, fj[1] - fi[1] + sy, fj[2] - fi[2] + sz]
                        r = math.sqrt(sum(x * x for x in matvec_row(df, h)))
                        if 1.0e-8 < r < best:
                            best = r
    return best * BOHR_ANG


def write_qe_input(
    path,
    prefix,
    struct,
    calculation="vc-md",
    nstep=1,
    tempw=3000.0,
    ecutwfc=80.0,
    ecutrho=320.0,
    kpoints="gamma",
    caimd=False,
    defgrad=None,
    legacy_strain=None,
    velocities=None,
    tprnfor=False,
    iprint=10,
    conv_thr="1.0d-7",
):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"! generated for {prefix}; calculation={calculation}",
        "&CONTROL",
        f"  calculation = '{calculation}',",
        "  restart_mode = 'from_scratch',",
        f"  prefix = '{prefix}',",
        f"  pseudo_dir = '{CAMPAIGN_ROOT}/pseudo/',",
        "  outdir = './tmp',",
    ]
    if calculation in ("md", "vc-md"):
        lines += [f"  nstep = {int(nstep)},", f"  dt = {DT_QE_2FS:.10f},"]
    lines += [
        f"  iprint = {int(iprint)},",
        "  tstress = .true.,",
        f"  tprnfor = {'.true.' if tprnfor else '.false.'},",
        "  disk_io = 'low',",
        "  verbosity = 'low',",
        "/",
        "&SYSTEM",
        "  nosym = .true.,",
        "  ibrav = 0,",
        f"  nat = {len(struct['fractional_positions'])},",
        "  ntyp = 1,",
        f"  ecutwfc = {float(ecutwfc):.8f},",
        f"  ecutrho = {float(ecutrho):.8f},",
        "  occupations = 'fixed',",
        "/",
        "&ELECTRONS",
        f"  conv_thr = {conv_thr},",
        "/",
    ]
    if calculation in ("md", "vc-md"):
        lines += [
            "&IONS",
            "  ion_dynamics = 'beeman',",
            "  ion_temperature = 'nose',",
        ]
        if velocities is not None:
            lines.append("  ion_velocities = 'from_input',")
        lines += [
            f"  tempw = {float(tempw):.8f},",
            "  fnosep(1) = 5.0,",
            "  nhpcl = 3,",
            "/",
        ]
    if calculation == "vc-md":
        lines += [
            "&CELL",
            "  cell_dynamics = 'pr',",
            "  cell_temperature = 'not_controlled',",
            "  press = 0.0,",
            "  cell_dofree = 'all',",
        ]
        if caimd:
            lines += ["  caimd = .true.,", "  caimd_stride = 10,"]
            if defgrad is not None:
                lines.append("  caimd_loading_mode = 'linear_defgrad',")
                for i in range(3):
                    for j in range(3):
                        if abs(defgrad[i][j]) > 0.0:
                            lines.append(f"  caimd_step_defgrad({i + 1},{j + 1}) = {defgrad[i][j]:.12e},")
            if legacy_strain is not None:
                for idx, value in enumerate(legacy_strain, start=1):
                    lines.append(f"  caimd_mask({idx}) = {'.true.' if abs(value) > 0 else '.false.'},")
                for idx, value in enumerate(legacy_strain, start=1):
                    lines.append(f"  caimd_step_strain({idx}) = {value:.12e},")
        lines.append("/")
    lines += ["ATOMIC_SPECIES", "C 12.01070000 C.UPF", "CELL_PARAMETERS bohr"]
    for row in struct["cell_bohr_rows"]:
        lines.append("  " + "  ".join(f"{x:.12f}" for x in row))
    lines.append("ATOMIC_POSITIONS crystal")
    for sym, pos in struct["fractional_positions"]:
        lines.append(f"{sym}  {pos[0]:.12f}  {pos[1]:.12f}  {pos[2]:.12f}")
    if velocities is not None:
        lines.append("ATOMIC_VELOCITIES a.u.")
        for (sym, _), vel in zip(struct["fractional_positions"], velocities):
            lines.append(f"{sym}  {vel[0]:.12e}  {vel[1]:.12e}  {vel[2]:.12e}")
    if kpoints == "gamma":
        lines.append("K_POINTS gamma")
    elif kpoints == "2x2x2":
        lines += ["K_POINTS automatic", "2 2 2 0 0 0"]
    else:
        lines += ["K_POINTS automatic", str(kpoints)]
    path.write_text("\n".join(lines) + "\n")
    return path


def run_cmd(cmd, cwd, out, err=None, env=None):
    cwd = Path(cwd)
    cwd.mkdir(parents=True, exist_ok=True)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    err = Path(err) if err else out.with_suffix(out.suffix + ".err")
    start = time.time()
    with open(out, "w") as fout, open(err, "w") as ferr:
        fout.write("COMMAND " + " ".join(map(str, cmd)) + "\n")
        fout.flush()
        proc = subprocess.run(list(map(str, cmd)), cwd=str(cwd), stdout=fout, stderr=ferr, env=env)
    return {"returncode": proc.returncode, "elapsed_seconds": time.time() - start, "stdout": str(out), "stderr": str(err)}


def srun_pw(input_path, run_dir, ranks=4, omp=8):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(omp)
    env["OMPI_MCA_coll_hcoll_enable"] = "0"
    cmd = ["srun", "--mpi=pmix_v3", "-n", str(ranks), str(PW_EXE), "-in", str(input_path)]
    return run_cmd(cmd, run_dir, Path(run_dir) / "qe.out", Path(run_dir) / "qe.err", env)


def parse_defgrad(path):
    rows = []
    path = Path(path)
    if not path.exists():
        return rows
    for line in path.read_text(errors="ignore").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        values = [float(x) for x in line.split()]
        if len(values) < 41:
            continue
        rows.append(
            {
                "ionic_step": int(values[0]),
                "event_index": int(values[1]),
                "det_f_before": values[2],
                "det_f_after": values[3],
                "det_a_event": values[4],
                "f_before": values[5:14],
                "f_after": values[14:23],
                "a_event": values[23:32],
                "green_lagrange": values[32:41],
                "velocity_jump_max": values[41] if len(values) > 41 else 0.0,
                "velocity_jump_rms": values[42] if len(values) > 42 else 0.0,
            }
        )
    return rows


def parse_caimd(path):
    samples = []
    path = Path(path)
    if not path.exists():
        return samples
    for line in path.read_text(errors="ignore").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        try:
            samples.append(
                {
                    "step": int(fields[0]),
                    "temp": float(fields[22]),
                    "volume": float(fields[23]),
                    "stress": [float(x) for x in fields[16:22]],
                }
            )
        except Exception:
            pass
    return samples


def parse_stdout_temperature(path):
    temps = []
    rx = re.compile(r"^\s*temperature\s+=\s*([-+0-9.Ee]+)\s*K")
    path = Path(path)
    if not path.exists():
        return temps
    for line in path.read_text(errors="ignore").splitlines():
        match = rx.search(line)
        if match:
            temps.append(float(match.group(1)))
    return temps


def temp_stats(temps, late=2000):
    if not temps:
        return {"count": 0, "finite": False}
    late_vals = temps[-min(late, len(temps)) :]
    mean = sum(late_vals) / len(late_vals)
    var = sum((x - mean) ** 2 for x in late_vals) / len(late_vals)
    half = max(1, len(late_vals) // 2)
    slope = sum(late_vals[-half:]) / half - sum(late_vals[:half]) / half
    return {
        "count": len(temps),
        "finite": all(math.isfinite(x) for x in temps),
        "min": min(temps),
        "max": max(temps),
        "late_count": len(late_vals),
        "late_mean": mean,
        "late_std": math.sqrt(var),
        "late_half_mean_delta": slope,
    }


def strip_tag(tag):
    return tag.split("}")[-1]


def extract_final_structure(xml_path, out_path, source):
    xml_path = Path(xml_path)
    root = ET.parse(xml_path).getroot()
    outputs = [e for e in root.iter() if strip_tag(e.tag) == "output"]
    scope = outputs[-1] if outputs else root
    atomic_structure = None
    for elem in scope.iter():
        if strip_tag(elem.tag) == "atomic_structure":
            atomic_structure = elem
    if atomic_structure is None:
        raise ValueError("no atomic_structure in QE XML")
    cell_node = None
    for elem in atomic_structure.iter():
        if strip_tag(elem.tag) == "cell":
            cell_node = elem
    rows = []
    for name in ("a1", "a2", "a3"):
        found = None
        for child in cell_node:
            if strip_tag(child.tag) == name:
                found = [float(x) for x in child.text.split()]
        if found is None:
            raise ValueError(f"missing {name} in cell")
        rows.append(found)
    atoms = []
    for elem in atomic_structure.iter():
        if strip_tag(elem.tag) == "atom":
            atoms.append((elem.attrib.get("name", "C").strip(), [float(x) for x in elem.text.split()]))
    hinv = inv3(rows)
    frac = []
    for sym, xyz in atoms:
        fr = [sum(xyz[j] * hinv[j][k] for j in range(3)) for k in range(3)]
        frac.append([sym, fr])
    struct = {
        "cell_bohr_rows": rows,
        "fractional_positions": frac,
        "n_atoms": len(frac),
        "determinant_bohr3": det3(rows),
        "source": source,
    }
    struct["min_CC_ang"] = min_cc_ang(struct)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(struct, indent=2) + "\n")
    return struct


def find_xml(run_dir, prefix):
    candidates = [
        Path(run_dir) / "tmp" / f"{prefix}.xml",
        Path(run_dir) / "tmp" / f"{prefix}.save/data-file-schema.xml",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(f"no QE XML for {prefix} in {run_dir}")


def validate_defgrad(prefix, run_dir, dmat, nstep, stride, tol=1.0e-10):
    rows = parse_defgrad(Path(run_dir) / "tmp" / f"{prefix}.caimd_defgrad")
    expected_events = (nstep - 1) // stride + 1
    errors = []
    if len(rows) != expected_events:
        errors.append(f"event count {len(rows)} != {expected_events}")
    for index, row in enumerate(rows, start=1):
        if row["event_index"] != index:
            errors.append(f"bad event_index at {index}")
        if row["ionic_step"] != 1 + (index - 1) * stride:
            errors.append(f"bad ionic_step at event {index}")
        f_after = fload(index, stride, dmat)
        f_before = fload(index - 1, stride, dmat)
        a_event = matmul(f_after, inv3(f_before))
        if max(abs(a - b) for a, b in zip(row["f_after"], flatten(f_after))) > tol:
            errors.append(f"F_after mismatch event {index}")
        if max(abs(a - b) for a, b in zip(row["a_event"], flatten(a_event))) > tol:
            errors.append(f"A_event mismatch event {index}")
        if row["det_f_after"] <= 0.0 or row["det_a_event"] <= 0.0:
            errors.append(f"non-positive determinant event {index}")
        if abs(row["velocity_jump_max"]) > 1.0e-8:
            errors.append(f"velocity jump too large event {index}")
    return {
        "pass": not errors,
        "errors": errors,
        "events": len(rows),
        "expected_events": expected_events,
        "last": rows[-1] if rows else None,
    }


def analyze_md_gate(stage, prefix, run_dir, nstep, temp_target, caimd=False, strict_3000=False):
    run_dir = Path(run_dir)
    rc = int((run_dir / "QE_RETURN_CODE.txt").read_text().strip()) if (run_dir / "QE_RETURN_CODE.txt").exists() else 999
    stdout = run_dir / "qe.out"
    text = stdout.read_text(errors="ignore") if stdout.exists() else ""
    temps = [x["temp"] for x in parse_caimd(run_dir / "tmp" / f"{prefix}.caimd")] if caimd else parse_stdout_temperature(stdout)
    stats = temp_stats(temps)
    gate = {"stage": stage, "returncode": rc, "job_done": "JOB DONE" in text, "temperature": stats, "errors": []}
    if rc != 0:
        gate["errors"].append(f"QE return code {rc}")
    if "JOB DONE" not in text:
        gate["errors"].append("JOB DONE missing")
    if stats.get("count", 0) < max(1, min(nstep, nstep // 2)):
        gate["errors"].append("too few temperature samples")
    if not stats.get("finite", False):
        gate["errors"].append("non-finite or missing temperature")
    if stats.get("max", 0.0) > max(6000.0, 2.5 * temp_target):
        gate["errors"].append("temperature runaway threshold exceeded")
    if strict_3000 and not (2850.0 <= stats.get("late_mean", -1.0e9) <= 3150.0):
        gate["errors"].append("T3000 late mean outside 2850-3150 K")
    try:
        final = extract_final_structure(find_xml(run_dir, prefix), run_dir / "final_structure.json", stage)
        gate["final_structure"] = {
            "n_atoms": final["n_atoms"],
            "determinant_bohr3": final["determinant_bohr3"],
            "min_CC_ang": final["min_CC_ang"],
        }
        if final["n_atoms"] != 64:
            gate["errors"].append("final atom count is not 64")
        if final["determinant_bohr3"] <= 0.0:
            gate["errors"].append("non-positive final cell determinant")
        if final["min_CC_ang"] < 0.7:
            gate["errors"].append("unphysical C-C distance")
    except Exception as exc:
        gate["errors"].append(f"final structure extraction failed: {exc}")
    gate["pass"] = not gate["errors"]
    (run_dir / "gate_report.json").write_text(json.dumps(gate, indent=2) + "\n")
    return gate


def parse_static(stdout):
    text = Path(stdout).read_text(errors="ignore") if Path(stdout).exists() else ""
    energy = None
    for match in re.finditer(r"!\s+total energy\s+=\s+([-+0-9.Ee]+)\s+Ry", text):
        energy = float(match.group(1))
    forces = [[float(a), float(b), float(c)] for a, b, c in re.findall(r"force\s+=\s*([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)", text)]
    stress = None
    lines = text.splitlines()
    for idx, line in enumerate(lines):
        if "total   stress" in line and idx + 3 < len(lines):
            values = []
            for off in (1, 2, 3):
                try:
                    values.extend(float(x) for x in lines[idx + off].split()[:3])
                except Exception:
                    values = []
                    break
            if len(values) == 9:
                stress = values
    return {"energy_Ry": energy, "forces_Ry_bohr": forces, "stress_kbar_9": stress, "job_done": "JOB DONE" in text}


def compare_static(ref, cand, nat=64):
    errors = []
    metrics = {}
    if ref["energy_Ry"] is None or cand["energy_Ry"] is None:
        errors.append("missing energy")
    else:
        metrics["dE_meV_atom"] = abs(cand["energy_Ry"] - ref["energy_Ry"]) * RY_TO_EV * 1000.0 / nat
        if metrics["dE_meV_atom"] > 10.0:
            errors.append("energy difference > 10 meV/atom")
    if len(ref["forces_Ry_bohr"]) == nat and len(cand["forces_Ry_bohr"]) == nat:
        conv = RY_TO_EV / BOHR_ANG
        diffs = []
        for a, b in zip(ref["forces_Ry_bohr"], cand["forces_Ry_bohr"]):
            diffs.extend((b[i] - a[i]) * conv for i in range(3))
        metrics["force_rms_eV_A"] = math.sqrt(sum(x * x for x in diffs) / len(diffs))
        metrics["force_max_eV_A"] = max(abs(x) for x in diffs)
        if metrics["force_rms_eV_A"] > 0.10:
            errors.append("force RMS > 0.10 eV/A")
        if metrics["force_max_eV_A"] > 0.30:
            errors.append("force max > 0.30 eV/A")
    else:
        errors.append("missing forces")
    if ref["stress_kbar_9"] and cand["stress_kbar_9"]:
        idx = [0, 4, 8, 5, 2, 1]
        diffs = [(cand["stress_kbar_9"][i] - ref["stress_kbar_9"][i]) * 0.1 for i in idx]
        metrics["stress_rms_GPa"] = math.sqrt(sum(x * x for x in diffs) / 6.0)
        if metrics["stress_rms_GPa"] > 2.0:
            errors.append("stress RMS > 2 GPa")
    else:
        errors.append("missing stress")
    if not cand["job_done"]:
        errors.append("JOB DONE missing")
    return {"pass": not errors, "errors": errors, "metrics": metrics}


def deform_structure(struct, fmat):
    h = struct["cell_bohr_rows"]
    new_h = [[sum(fmat[i][k] * h[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    return {
        "cell_bohr_rows": new_h,
        "fractional_positions": struct["fractional_positions"],
        "n_atoms": struct["n_atoms"],
        "source": struct.get("source", "") + " deformed",
    }


def read_selected_env():
    pass_file = CAMPAIGN_ROOT / "validation/VALIDATION_PASS"
    env_file = CAMPAIGN_ROOT / "validation/selected_settings.env"
    sha_file = CAMPAIGN_ROOT / "validation/selected_settings.env.sha256"
    if not pass_file.exists():
        raise SystemExit("VALIDATION_PASS missing")
    if not env_file.exists() or not sha_file.exists():
        raise SystemExit("selected_settings files missing")
    expected = sha_file.read_text().split()[0]
    actual = sha256(env_file)
    if expected != actual:
        raise SystemExit(f"selected_settings.env SHA mismatch: {actual} != {expected}")
    settings = {}
    for line in env_file.read_text().splitlines():
        if line and not line.startswith("#") and "=" in line:
            key, value = line.split("=", 1)
            settings[key] = value.strip().strip("'\"")
    if sha256(settings["PW_EXE"]) != settings["PW_SHA256"]:
        raise SystemExit("PW_EXE SHA mismatch")
    return settings


def selected_runtime_settings():
    data = read_selected_env()
    return float(data["ECUTWFC"]), float(data["ECUTRHO"]), data.get("KPOINT_MODE", "gamma"), int(data["MPI_RANKS"]), int(data["OMP_THREADS"])


def deterministic_velocities(struct, temp=3000.0, seed=1113000):
    rng = random.Random(seed)
    nat = len(struct["fractional_positions"])
    mass = 12.0107
    sigma = math.sqrt(K_BOLTZMANN_RY * temp / mass)
    velocities = [[rng.gauss(0.0, sigma) for _ in range(3)] for _ in range(nat)]
    center = [sum(v[j] for v in velocities) / nat for j in range(3)]
    for vel in velocities:
        for j in range(3):
            vel[j] -= center[j]
    ekin = sum(0.5 * mass * sum(x * x for x in vel) for vel in velocities)
    tcalc = 2.0 * ekin / (3.0 * (nat - 1) * K_BOLTZMANN_RY)
    factor = math.sqrt(temp / tcalc)
    for vel in velocities:
        for j in range(3):
            vel[j] *= factor
    ekin = sum(0.5 * mass * sum(x * x for x in vel) for vel in velocities)
    return velocities, {
        "seed": seed,
        "temperature_formula_K": 2.0 * ekin / (3.0 * (nat - 1) * K_BOLTZMANN_RY),
        "center_of_mass_velocity": [sum(v[j] for v in velocities) / nat for j in range(3)],
    }


def validation_main():
    valdir = CAMPAIGN_ROOT / "validation"
    runs = valdir / "runs"
    benches = valdir / "benchmarks"
    status = {"status": "RUNNING", "timestamp_start": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "tests": {}, "errors": []}

    def save():
        (valdir / "VALIDATION_STATUS.json").write_text(json.dumps(status, indent=2) + "\n")

    try:
        v0 = {
            "pw_exe": str(PW_EXE),
            "pw_sha256": sha256(PW_EXE),
            "old_pw_sha256": sha256(OLD_PW_EXE),
            "pseudo_sha256": sha256(PSEUDO),
            "source_patch_sha256": sha256(SOFTWARE_ROOT / "src-current/caimd_linear_defgrad_phase1/evidence/linear_defgrad_phase1.patch"),
        }
        v0["pass"] = (
            v0["pw_sha256"] == EXPECTED["pw_sha"]
            and v0["old_pw_sha256"] == EXPECTED["old_pw_sha"]
            and v0["pseudo_sha256"] == EXPECTED["pseudo_sha"]
            and v0["source_patch_sha256"] == EXPECTED["source_patch_sha"]
        )
        status["tests"]["V0_executable_identity"] = v0
        save()
        if not v0["pass"]:
            raise RuntimeError("V0 identity failed")

        struct = read_structure(ORIENTED)
        d0 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        parser_cases = [
            ("parser_d13", [[0.0, 0.0, 2.0e-6], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            ("parser_d33", [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, -2.0e-6]]),
        ]
        v1_cases = []
        for name, dmat in parser_cases:
            run_dir = runs / name
            inp = CAMPAIGN_ROOT / "inputs/generated" / f"validation_{name}.in"
            write_qe_input(inp, f"v_{name}", struct, nstep=1, tempw=300.0, ecutwfc=45.0, ecutrho=180.0, caimd=True, defgrad=dmat, iprint=1)
            result = srun_pw(inp, run_dir, 1, 8)
            text = (run_dir / "qe.out").read_text(errors="ignore")
            ok = result["returncode"] == 0 and "caimd_step_defgrad matrix" in text and ("2.000000E-06" in text or "-2.000000E-06" in text)
            v1_cases.append({"case": name, "run": result, "pass": ok})
        status["tests"]["V1_parser"] = {"pass": all(x["pass"] for x in v1_cases), "cases": v1_cases}
        save()
        if not status["tests"]["V1_parser"]["pass"]:
            raise RuntimeError("V1 parser failed")

        run_dir = runs / "identity_zero"
        inp = CAMPAIGN_ROOT / "inputs/generated/validation_identity_zero.in"
        write_qe_input(inp, "v_identity_zero", struct, nstep=20, tempw=300.0, ecutwfc=45.0, ecutrho=180.0, caimd=True, defgrad=d0, iprint=1)
        result = srun_pw(inp, run_dir, 4, 8)
        dg = validate_defgrad("v_identity_zero", run_dir, d0, 20, 10)
        status["tests"]["V2_identity_path"] = {"run": result, "defgrad": dg, "pass": result["returncode"] == 0 and dg["pass"]}
        save()
        if not status["tests"]["V2_identity_path"]["pass"]:
            raise RuntimeError("V2 identity failed")

        loading_cases = [
            ("V3_short_compression", "v_short_comp", [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, -2.0e-6]], {"F33": 0.998}),
            ("V4_short_shear", "v_short_shear", [[0.0, 0.0, 2.0e-6], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], {"F13": 0.002}),
        ]
        for key, prefix, dmat, endpoint in loading_cases:
            run_dir = runs / key
            inp = CAMPAIGN_ROOT / "inputs/generated" / f"{key}.in"
            write_qe_input(inp, prefix, struct, nstep=1000, tempw=3000.0, ecutwfc=80.0, ecutrho=320.0, caimd=True, defgrad=dmat, iprint=10)
            result = srun_pw(inp, run_dir, 4, 8)
            (run_dir / "QE_RETURN_CODE.txt").write_text(str(result["returncode"]) + "\n")
            dg = validate_defgrad(prefix, run_dir, dmat, 1000, 10)
            gate = analyze_md_gate(key, prefix, run_dir, 1000, 3000.0, caimd=True)
            status["tests"][key] = {"run": result, "defgrad": dg, "md_gate": gate, "endpoint_expected": endpoint, "pass": result["returncode"] == 0 and dg["pass"] and gate["pass"]}
            save()
            if not status["tests"][key]["pass"]:
                raise RuntimeError(f"{key} failed")

        run_dir = runs / "legacy_regression"
        inp = CAMPAIGN_ROOT / "inputs/generated/validation_legacy_regression.in"
        write_qe_input(inp, "v_legacy", struct, nstep=30, tempw=300.0, ecutwfc=45.0, ecutrho=180.0, caimd=True, legacy_strain=[1.0e-6, 0.0, 0.0, 0.0, 0.0, 0.0], iprint=1)
        result = srun_pw(inp, run_dir, 1, 8)
        status["tests"]["V5_legacy_regression"] = {"run": result, "caimd_samples": len(parse_caimd(run_dir / "tmp/v_legacy.caimd")), "pass": result["returncode"] == 0}
        save()
        if not status["tests"]["V5_legacy_regression"]["pass"]:
            raise RuntimeError("V5 legacy regression failed")

        source = (SOFTWARE_ROOT / "src-current/PW/src/caimd_module.f90").read_text()
        phrases = ["h_pr = h_candidate", "h_new = MATMUL( a_event, h_pr )", "h_candidate = h_new", "TRIM( caimd_loading_mode_ ) == 'linear_defgrad' ) RETURN"]
        status["tests"]["V6_PR_candidate_preservation"] = {
            "source_phrases_present": {phrase: phrase in source for phrase in phrases},
            "direct_H_PR_candidate_columns_available": False,
            "explanation": "Runtime output has no H_PR_candidate columns; source-level gate verifies left multiplication of the current PR candidate and tangent-projection skip for linear_defgrad.",
        }
        status["tests"]["V6_PR_candidate_preservation"]["pass"] = all(phrase in source for phrase in phrases)
        save()
        if not status["tests"]["V6_PR_candidate_preservation"]["pass"]:
            raise RuntimeError("V6 source gate failed")

        comp = status["tests"]["V3_short_compression"]["defgrad"]["last"]
        shear = status["tests"]["V4_short_shear"]["defgrad"]["last"]
        status["tests"]["V7_loading_path_sanity"] = {
            "compression_F33": comp["f_after"][8] if comp else None,
            "shear_F13": shear["f_after"][2] if shear else None,
            "pass": comp is not None and shear is not None and abs(comp["f_after"][8] - 0.998) < 1.0e-12 and abs(shear["f_after"][2] - 0.002) < 1.0e-12,
        }
        save()
        if not status["tests"]["V7_loading_path_sanity"]["pass"]:
            raise RuntimeError("V7 path sanity failed")

        run_dir = runs / "T3000_preflight"
        inp = CAMPAIGN_ROOT / "inputs/generated/validation_T3000_preflight.in"
        write_qe_input(inp, "v_T3000_pre", struct, nstep=2000, tempw=3000.0, ecutwfc=80.0, ecutrho=320.0, caimd=True, defgrad=d0, iprint=10)
        result = srun_pw(inp, run_dir, 4, 8)
        (run_dir / "QE_RETURN_CODE.txt").write_text(str(result["returncode"]) + "\n")
        gate = analyze_md_gate("T3000_preflight", "v_T3000_pre", run_dir, 2000, 3000.0, caimd=True)
        stats = gate["temperature"]
        thermal_ok = result["returncode"] == 0 and gate["pass"] and stats.get("max", 99999.0) < 6000.0 and 2700.0 <= stats.get("late_mean", -1.0) <= 3300.0
        status["tests"]["V8_3000K_preflight"] = {"run": result, "md_gate": gate, "physical_ionic_dof": 189, "pass": thermal_ok}
        save()
        if not thermal_ok:
            raise RuntimeError("V8 thermal preflight failed")

        ref = ("REF_GAMMA", 80.0, 320.0)
        candidates = [("C1", 60.0, 240.0), ("C2", 50.0, 200.0), ("C3", 45.0, 180.0)]
        geometries = {
            "G0": struct,
            "GC": deform_structure(struct, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.85]]),
            "GS": deform_structure(struct, [[1.0, 0.0, 0.20], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        }
        static = {"allowed_candidates": [], "comparisons": {}, "fallback_to_REF_GAMMA": False}
        refs = {}
        for geom, stc in geometries.items():
            run_dir = benches / f"static_REF_GAMMA_{geom}"
            inp = CAMPAIGN_ROOT / "inputs/generated" / f"static_REF_GAMMA_{geom}.in"
            write_qe_input(inp, f"sp_ref_{geom.lower()}", stc, calculation="scf", ecutwfc=80.0, ecutrho=320.0, kpoints="gamma", tprnfor=True, conv_thr="1.0d-8")
            rr = srun_pw(inp, run_dir, 4, 8)
            refs[geom] = parse_static(run_dir / "qe.out")
            refs[geom]["run"] = rr
        for geom in ("G0", "GC"):
            run_dir = benches / f"static_REF_2X2X2_{geom}"
            inp = CAMPAIGN_ROOT / "inputs/generated" / f"static_REF_2X2X2_{geom}.in"
            write_qe_input(inp, f"sp_ref2_{geom.lower()}", geometries[geom], calculation="scf", ecutwfc=80.0, ecutrho=320.0, kpoints="2x2x2", tprnfor=True, conv_thr="1.0d-8")
            rr = srun_pw(inp, run_dir, 4, 8)
            static.setdefault("ref_2x2x2", {})[geom] = {"run": rr, "parsed": parse_static(run_dir / "qe.out")}
        for cname, ew, er in candidates:
            candidate_ok = True
            data = {}
            for geom, stc in geometries.items():
                run_dir = benches / f"static_{cname}_{geom}"
                inp = CAMPAIGN_ROOT / "inputs/generated" / f"static_{cname}_{geom}.in"
                write_qe_input(inp, f"sp_{cname.lower()}_{geom.lower()}", stc, calculation="scf", ecutwfc=ew, ecutrho=er, kpoints="gamma", tprnfor=True, conv_thr="1.0d-8")
                rr = srun_pw(inp, run_dir, 4, 8)
                parsed = parse_static(run_dir / "qe.out")
                cmp = compare_static(refs[geom], parsed)
                data[geom] = {"run": rr, "parsed": parsed, "comparison_to_REF_GAMMA": cmp}
                candidate_ok = candidate_ok and rr["returncode"] == 0 and cmp["pass"]
            static["comparisons"][cname] = {"ecutwfc": ew, "ecutrho": er, "pass": candidate_ok, "geometries": data}
            if candidate_ok:
                static["allowed_candidates"].append((cname, ew, er))
        if not static["allowed_candidates"]:
            static["allowed_candidates"].append(ref)
            static["fallback_to_REF_GAMMA"] = True
        status["tests"]["candidate_static_accuracy"] = static
        save()

        perf = {"results": [], "eligible": []}
        bench_structs = {"early": struct, "compression_final": geometries["GC"], "shear_final": geometries["GS"]}
        for cname, ew, er in static["allowed_candidates"]:
            worst = 0.0
            for geom, stc in bench_structs.items():
                run_dir = benches / f"perf_{cname}_L4_{geom}"
                inp = CAMPAIGN_ROOT / "inputs/generated" / f"perf_{cname}_L4_{geom}.in"
                write_qe_input(inp, f"perf_{cname.lower()}_{geom}", stc, nstep=320, tempw=3000.0, ecutwfc=ew, ecutrho=er, caimd=True, defgrad=d0, iprint=10)
                rr = srun_pw(inp, run_dir, 4, 8)
                seconds_per_step = rr["elapsed_seconds"] / 300.0
                perf["results"].append({"candidate": cname, "layout": "L4", "geometry": geom, "run": rr, "steady_seconds_per_step_estimate": seconds_per_step})
                worst = max(worst, seconds_per_step if rr["returncode"] == 0 else 1.0e9)
            compression_raw = worst * 75000.0 / 3600.0
            shear_raw = worst * 100000.0 / 3600.0
            entry = {
                "candidate": cname,
                "ecutwfc": ew,
                "ecutrho": er,
                "layout": "L4",
                "worst_seconds_per_step": worst,
                "compression_raw_hours": compression_raw,
                "compression_safe_hours": 1.30 * compression_raw + 1.0,
                "shear_raw_hours": shear_raw,
                "shear_safe_hours": 1.30 * shear_raw + 1.0,
            }
            entry["pass"] = entry["compression_safe_hours"] <= 23.0 and entry["shear_safe_hours"] <= 23.0
            if entry["pass"]:
                perf["eligible"].append(entry)
        status["tests"]["performance_benchmark"] = perf
        save()
        if not perf["eligible"]:
            status["status"] = "NO_24H_FEASIBLE_CONFIGURATION"
            save()
            raise RuntimeError("NO_24H_FEASIBLE_CONFIGURATION")

        selected = sorted(perf["eligible"], key=lambda x: x["worst_seconds_per_step"])[0]
        selected_json = {
            "selection": selected,
            "KPOINT_MODE": "gamma",
            "MPI_RANKS": 4,
            "OMP_THREADS": 8,
            "GPU_COUNT": 4,
            "QUALITATIVE_ONLY": True,
            "label": "GAMMA_ONLY_QUALITATIVE_PILOT",
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }
        env_lines = [
            f"PW_EXE={PW_EXE}",
            f"PW_SHA256={EXPECTED['pw_sha']}",
            f"ECUTWFC={selected['ecutwfc']}",
            f"ECUTRHO={selected['ecutrho']}",
            "KPOINT_MODE=gamma",
            "MPI_RANKS=4",
            "OMP_THREADS=8",
            "GPU_COUNT=4",
            f"SECONDS_PER_STEP_WORST={selected['worst_seconds_per_step']}",
            f"COMPRESSION_RAW_HOURS={selected['compression_raw_hours']}",
            f"COMPRESSION_SAFE_HOURS={selected['compression_safe_hours']}",
            f"SHEAR_RAW_HOURS={selected['shear_raw_hours']}",
            f"SHEAR_SAFE_HOURS={selected['shear_safe_hours']}",
            "QUALITATIVE_ONLY=true",
            f"VALIDATION_TIMESTAMP={selected_json['timestamp']}",
        ]
        env_path = valdir / "selected_settings.env"
        env_path.write_text("\n".join(env_lines) + "\n")
        (valdir / "selected_settings.env.sha256").write_text(f"{sha256(env_path)}  selected_settings.env\n")
        (valdir / "selected_settings.json").write_text(json.dumps(selected_json, indent=2) + "\n")
        (valdir / "VALIDATION_PASS").write_text("VALIDATION_PASS\n")
        status["status"] = "VALIDATION_PASS"
        status["selected_settings"] = selected_json
        save()
    except Exception as exc:
        status["errors"].append(str(exc))
        if status["status"] == "RUNNING":
            status["status"] = "VALIDATION_FAIL"
        save()
        raise


def run_eq(stage, temp, prev_stage=None):
    ecutwfc, ecutrho, kpoints, ranks, omp = selected_runtime_settings()
    run_dir = CAMPAIGN_ROOT / "equilibration" / stage
    run_dir.mkdir(parents=True, exist_ok=True)
    if prev_stage:
        prev_gate = CAMPAIGN_ROOT / "equilibration" / prev_stage / "gate_report.json"
        if not prev_gate.exists() or not json.loads(prev_gate.read_text()).get("pass"):
            raise SystemExit(f"predecessor gate missing or failed: {prev_gate}")
        struct = read_structure(CAMPAIGN_ROOT / "equilibration" / prev_stage / "final_structure.json")
    else:
        struct = read_structure(ORIENTED)
    prefix = f"d64_{stage.lower()}"
    inp = CAMPAIGN_ROOT / "inputs/generated" / f"{stage}.in"
    write_qe_input(inp, prefix, struct, nstep=3000, tempw=temp, ecutwfc=ecutwfc, ecutrho=ecutrho, kpoints=kpoints, caimd=False, iprint=10)
    result = srun_pw(inp, run_dir, ranks, omp)
    (run_dir / "QE_RETURN_CODE.txt").write_text(str(result["returncode"]) + "\n")
    gate = analyze_md_gate(stage, prefix, run_dir, 3000, temp, caimd=False, strict_3000=(stage == "T3000"))
    if not gate["pass"]:
        raise SystemExit(f"{stage} gate failed")


def prepare_loading_start():
    read_selected_env()
    gate_path = CAMPAIGN_ROOT / "equilibration/T3000/gate_report.json"
    if not gate_path.exists() or not json.loads(gate_path.read_text()).get("pass"):
        raise SystemExit("T3000 gate missing or failed")
    struct = read_structure(CAMPAIGN_ROOT / "equilibration/T3000/final_structure.json")
    velocities, metadata = deterministic_velocities(struct, 3000.0, 1113000)
    data = {
        "structure": struct,
        "velocities_au": velocities,
        "velocity_metadata": metadata,
        "same_state_for_compression_and_shear": True,
    }
    out = CAMPAIGN_ROOT / "prepare_loading_start/loading_start.json"
    out.write_text(json.dumps(data, indent=2) + "\n")
    lines = ["CELL_PARAMETERS bohr"]
    for row in struct["cell_bohr_rows"]:
        lines.append("  " + "  ".join(f"{x:.12f}" for x in row))
    lines.append("ATOMIC_POSITIONS crystal")
    for sym, pos in struct["fractional_positions"]:
        lines.append(f"{sym}  {pos[0]:.12f}  {pos[1]:.12f}  {pos[2]:.12f}")
    lines.append("ATOMIC_VELOCITIES a.u.")
    for (sym, _), vel in zip(struct["fractional_positions"], velocities):
        lines.append(f"{sym}  {vel[0]:.12e}  {vel[1]:.12e}  {vel[2]:.12e}")
    (CAMPAIGN_ROOT / "prepare_loading_start/loading_start.in.fragment").write_text("\n".join(lines) + "\n")
    (CAMPAIGN_ROOT / "prepare_loading_start/loading_start_report.md").write_text(
        "# Loading Start\n\n"
        f"Velocity handoff from QE XML was not used. A deterministic 3000 K velocity set was generated with seed {metadata['seed']} "
        "using the QE kinetic-temperature formula and is shared by compression and shear.\n\n"
        f"loading_start_sha256={sha256(out)}\n"
    )
    (CAMPAIGN_ROOT / "prepare_loading_start/PREPARE_LOADING_START_PASS").write_text("PREPARE_LOADING_START_PASS\n")


def run_loading(kind):
    ecutwfc, ecutrho, kpoints, ranks, omp = selected_runtime_settings()
    if not (CAMPAIGN_ROOT / "prepare_loading_start/PREPARE_LOADING_START_PASS").exists():
        raise SystemExit("prepare_loading_start gate missing")
    data = json.loads((CAMPAIGN_ROOT / "prepare_loading_start/loading_start.json").read_text())
    struct = data["structure"]
    velocities = data["velocities_au"]
    if kind == "compression":
        run_dir = CAMPAIGN_ROOT / "compression_111"
        prefix = "d64_compression_111"
        nstep = 75000
        dmat = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, -2.0e-6]]
    else:
        run_dir = CAMPAIGN_ROOT / "shear_111_11m2"
        prefix = "d64_shear_111_11m2"
        nstep = 100000
        dmat = [[0.0, 0.0, 2.0e-6], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    run_dir.mkdir(parents=True, exist_ok=True)
    inp = CAMPAIGN_ROOT / "inputs/generated" / f"{prefix}.in"
    write_qe_input(inp, prefix, struct, nstep=nstep, tempw=3000.0, ecutwfc=ecutwfc, ecutrho=ecutrho, kpoints=kpoints, caimd=True, defgrad=dmat, velocities=velocities, iprint=10)
    result = srun_pw(inp, run_dir, ranks, omp)
    (run_dir / "QE_RETURN_CODE.txt").write_text(str(result["returncode"]) + "\n")
    gate = analyze_md_gate(kind, prefix, run_dir, nstep, 3000.0, caimd=True)
    defgrad = validate_defgrad(prefix, run_dir, dmat, nstep, 10)
    gate["defgrad_gate"] = defgrad
    if kind == "compression":
        gate["expected"] = {"events": 7500, "final_F33": 0.85, "Ezz": -0.13875}
        gate["pass"] = gate["pass"] and defgrad["pass"] and defgrad["events"] == 7500 and abs(defgrad["last"]["f_after"][8] - 0.85) < 1.0e-10
    else:
        gate["expected"] = {"events": 10000, "final_F13": 0.20, "detF": 1.0, "Exz": 0.10, "Ezx": 0.10, "Ezz": 0.02}
        gate["pass"] = gate["pass"] and defgrad["pass"] and defgrad["events"] == 10000 and abs(defgrad["last"]["f_after"][2] - 0.20) < 1.0e-10 and abs(defgrad["last"]["det_f_after"] - 1.0) < 1.0e-10
    (run_dir / "gate_report.json").write_text(json.dumps(gate, indent=2) + "\n")
    if not gate["pass"]:
        raise SystemExit(f"{kind} loading gate failed")


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("validate")
    p_eq = sub.add_parser("run-eq")
    p_eq.add_argument("stage")
    p_eq.add_argument("temp", type=float)
    p_eq.add_argument("--prev")
    sub.add_parser("prepare-loading-start")
    p_load = sub.add_parser("run-loading")
    p_load.add_argument("kind", choices=["compression", "shear"])
    args = parser.parse_args()
    if args.cmd == "validate":
        validation_main()
    elif args.cmd == "run-eq":
        run_eq(args.stage, args.temp, args.prev)
    elif args.cmd == "prepare-loading-start":
        prepare_loading_start()
    elif args.cmd == "run-loading":
        run_loading(args.kind)


if __name__ == "__main__":
    main()
