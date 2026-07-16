#!/usr/bin/env python3
"""CAIMD V1 software validation harness for local serial and Leonardo runs.

The harness generates its own small QE inputs, refuses to run without a local
pseudopotential, records executable checksums and command logs, and writes a
machine-readable summary for expected successes and expected rejections.
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
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LOCAL_PSEUDO = REPO_ROOT / "pseudo" / "C.UPF"
DEFAULT_REMOTE_ROOT = Path(
    os.environ.get(
        "CAIMD_REMOTE_ROOT",
        "/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919",
    )
)
FLOAT_RE = re.compile(
    r"[-+]?(?:(?:\d+\.\d*|\.\d+)(?:[EeDd]?[-+]?\d+)?|\d+(?:[EeDd][-+]?\d+)?)"
)
PSEUDO_NAME = "C.UPF"


@dataclass(frozen=True)
class CaseSpec:
    name: str
    group: str
    code: str
    input_key: str
    executable: str = "current"
    ranks: int = 1
    expect: str = "success"
    tokens: tuple[str, ...] = ()
    monitor_gpu: bool = False
    compare_to: str | None = None
    description: str = ""


@dataclass
class CaseResult:
    name: str
    group: str
    code: str
    executable: str | None
    executable_sha256: str | None
    ranks: int
    input_path: str | None
    run_dir: str | None
    output_path: str | None
    command: list[str]
    returncode: int | None
    elapsed_s: float
    expected: str
    success: bool
    skipped: bool = False
    reason: str | None = None
    job_done: bool = False
    mpi_abort: bool = False
    missing_library: bool = False
    timeout: bool = False
    metrics: dict[str, Any] = field(default_factory=dict)


def utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run_capture(
    cmd: list[str], cwd: Path | None = None, timeout: int | None = None
) -> tuple[int, str, bool]:
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
        )
        return proc.returncode, proc.stdout, False
    except subprocess.TimeoutExpired as exc:
        out = exc.stdout or ""
        if isinstance(out, bytes):
            out = out.decode(errors="replace")
        return 124, out + "\nTIMEOUT\n", True


def extract_floats(text: str) -> list[float]:
    vals: list[float] = []
    for raw in FLOAT_RE.findall(text):
        token = raw.replace("D", "E").replace("d", "e")
        if "e" not in token.lower():
            compact = re.match(r"^([+-]?(?:\d+\.\d*|\.\d+))([-+]\d+)$", token)
            if compact:
                token = compact.group(1) + "E" + compact.group(2)
        try:
            vals.append(float(token))
        except ValueError:
            continue
    return vals


def read_float_file(path: Path) -> list[float]:
    if not path.exists():
        return []
    return extract_floats(path.read_text(errors="ignore"))


def compare_float_files(a: Path, b: Path) -> dict[str, Any]:
    av = read_float_file(a)
    bv = read_float_file(b)
    n = min(len(av), len(bv))
    if n == 0:
        return {"count_a": len(av), "count_b": len(bv), "max_abs": None, "rms": None}
    diff = [av[i] - bv[i] for i in range(n)]
    max_abs = max(abs(x) for x in diff)
    rms = math.sqrt(sum(x * x for x in diff) / n)
    denom = max(max(abs(x) for x in av[:n]), max(abs(x) for x in bv[:n]), 1.0)
    return {
        "count_a": len(av),
        "count_b": len(bv),
        "same_count": len(av) == len(bv),
        "max_abs": max_abs,
        "rms": rms,
        "relative_max": max_abs / denom,
    }


def parse_energy(output: str) -> list[float]:
    vals: list[float] = []
    for line in output.splitlines():
        match = re.search(r"!\s+total energy\s+=\s+([-+0-9.EeDd]+)", line)
        if match:
            vals.append(float(match.group(1).replace("D", "E").replace("d", "e")))
    return vals


def parse_caimd(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"exists": False}
    rows: list[list[Any]] = []
    headers: list[str] = []
    for line in path.read_text(errors="ignore").splitlines():
        if not line.strip():
            continue
        if line.lstrip().startswith("#"):
            headers.append(line.strip())
            continue
        parts = line.split()
        parsed: list[Any] = []
        for idx, part in enumerate(parts):
            if idx == 3:
                parsed.append(part.upper().startswith("T"))
            else:
                parsed.append(float(part.replace("D", "E").replace("d", "e")))
        rows.append(parsed)
    column_counts = sorted({len(row) for row in rows})
    return {
        "exists": True,
        "path": str(path),
        "rows": len(rows),
        "header_lines": headers,
        "column_counts": column_counts,
        "single_column_count_26": column_counts == [26],
        "event_flags": [bool(row[3]) for row in rows] if rows else [],
        "event_counters": [int(round(row[2])) for row in rows] if rows else [],
        "temperature_finite": all(math.isfinite(float(row[22])) for row in rows) if rows else False,
        "volume_finite": all(math.isfinite(float(row[23])) for row in rows) if rows else False,
        "velocity_residual_finite": all(
            math.isfinite(float(row[24])) and math.isfinite(float(row[25])) for row in rows
        )
        if rows
        else False,
        "velocity_residual_max": max((float(row[24]) for row in rows), default=None),
        "velocity_residual_rms_max": max((float(row[25]) for row in rows), default=None),
        "nominal_final": rows[-1][4:10] if rows else None,
        "actual_final": rows[-1][10:16] if rows else None,
    }


def det3(a: list[list[float]]) -> float:
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def inv3(a: list[list[float]]) -> list[list[float]]:
    det = det3(a)
    if abs(det) < 1.0e-18:
        raise ValueError("singular 3x3 matrix")
    return [
        [
            (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det,
            -(a[0][1] * a[2][2] - a[0][2] * a[2][1]) / det,
            (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det,
        ],
        [
            -(a[1][0] * a[2][2] - a[1][2] * a[2][0]) / det,
            (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det,
            -(a[0][0] * a[1][2] - a[0][2] * a[1][0]) / det,
        ],
        [
            (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det,
            -(a[0][0] * a[2][1] - a[0][1] * a[2][0]) / det,
            (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det,
        ],
    ]


def matmul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def frob(a: list[list[float]]) -> float:
    return math.sqrt(sum(a[i][j] * a[i][j] for i in range(3) for j in range(3)))


def parse_cell_blocks(path: Path) -> list[list[list[float]]]:
    if not path.exists():
        return []
    blocks: list[list[list[float]]] = []
    rows: list[list[float]] = []
    for line in path.read_text(errors="ignore").splitlines():
        vals = extract_floats(line)
        if len(vals) == 3:
            rows.append(vals[:3])
            if len(rows) == 3:
                blocks.append(rows)
                rows = []
        elif rows and vals:
            rows = []
    return blocks


def rotation_metrics(cel_path: Path) -> dict[str, Any]:
    blocks = parse_cell_blocks(cel_path)
    if len(blocks) < 2:
        return {"cell_blocks": len(blocks)}
    e_norms: list[float] = []
    w_norms: list[float] = []
    ratios: list[float] = []
    for h0, h1 in zip(blocks[:-1], blocks[1:]):
        try:
            f = matmul(h1, inv3(h0))
        except ValueError:
            continue
        l = [[f[i][j] - (1.0 if i == j else 0.0) for j in range(3)] for i in range(3)]
        e = [[0.5 * (l[i][j] + l[j][i]) for j in range(3)] for i in range(3)]
        w = [[0.5 * (l[i][j] - l[j][i]) for j in range(3)] for i in range(3)]
        en = frob(e)
        wn = frob(w)
        e_norms.append(en)
        w_norms.append(wn)
        ratios.append(wn / en if en > 0.0 else 0.0)
    if not e_norms:
        return {"cell_blocks": len(blocks), "valid_cell_pairs": 0}
    return {
        "cell_blocks": len(blocks),
        "valid_cell_pairs": len(e_norms),
        "max_E_fro": max(e_norms),
        "rms_E_fro": math.sqrt(sum(x * x for x in e_norms) / len(e_norms)),
        "max_W_fro": max(w_norms),
        "rms_W_fro": math.sqrt(sum(x * x for x in w_norms) / len(w_norms)),
        "max_skew_to_strain": max(ratios),
    }


def parse_gpu_monitor(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"exists": False}
    samples: list[dict[str, Any]] = []
    for row in csv.reader(path.read_text(errors="ignore").splitlines()):
        if len(row) < 5:
            continue
        try:
            samples.append(
                {
                    "timestamp": row[0].strip(),
                    "index": int(row[1]),
                    "utilization_gpu_pct": float(row[2]),
                    "memory_used_mib": float(row[3]),
                    "power_w": float(row[4]),
                }
            )
        except ValueError:
            continue
    if not samples:
        return {"exists": True, "samples": 0}
    return {
        "exists": True,
        "samples": len(samples),
        "max_utilization_gpu_pct": max(s["utilization_gpu_pct"] for s in samples),
        "max_memory_used_mib": max(s["memory_used_mib"] for s in samples),
        "max_power_w": max(s["power_w"] for s in samples),
    }


def cell_text(kind: str = "orthogonal") -> str:
    if kind == "nonorthogonal":
        rows = [(6.80, 0.00, 0.00), (0.75, 6.95, 0.00), (0.30, 0.45, 7.10)]
    elif kind == "anisotropic":
        rows = [(6.20, 0.00, 0.00), (0.00, 7.30, 0.00), (0.00, 0.00, 8.10)]
    else:
        rows = [(6.74, 0.00, 0.00), (0.00, 6.74, 0.00), (0.00, 0.00, 6.74)]
    return "\n".join(f" {a:14.8f} {b:14.8f} {c:14.8f}" for a, b, c in rows)


def positions_text(nat: int = 2, fixed_flags: bool = False) -> str:
    pts = [(0.00, 0.00, 0.00), (0.25, 0.25, 0.25)]
    rows = []
    for idx in range(nat):
        x, y, z = pts[idx % len(pts)]
        flags = " 0 1 1" if fixed_flags and idx == 0 else ""
        rows.append(f" C {x:10.6f} {y:10.6f} {z:10.6f}{flags}")
    return "\n".join(rows)


def velocities_text(nat: int = 2) -> str:
    base = [(0.00020, 0.00010, -0.00015), (-0.00016, 0.00013, 0.00011)]
    vals = [base[i % len(base)] for i in range(nat)]
    return "\n".join(f" C {x: .8e} {y: .8e} {z: .8e}" for x, y, z in vals)


def qe_input(
    prefix: str,
    calculation: str,
    pseudo_dir: Path,
    *,
    nat: int = 2,
    nstep: int = 4,
    cell_kind: str = "orthogonal",
    caimd: str = "",
    include_cell: bool = True,
    velocities: bool = False,
    fixed_flags: bool = False,
    extra_control: str = "",
    extra_ions: str = "",
    extra_system: str = "",
    extra_cell: str = "",
    after_cards: str = "",
) -> str:
    ions = ""
    if calculation in {"md", "vc-md", "relax", "vc-relax"}:
        if calculation == "md":
            dyn = "ion_dynamics='verlet',"
        elif calculation == "vc-md":
            dyn = "ion_dynamics='beeman',"
        else:
            dyn = "ion_dynamics='bfgs',"
        vel = "ion_velocities='from_input'," if velocities else ""
        ion_temp = "" if "ion_temperature" in extra_ions.lower() else "ion_temperature='not_controlled',"
        ions = f"""
&IONS
 {dyn}
 {ion_temp}
 {vel}
 {extra_ions}
/
"""
    cell = ""
    if include_cell:
        cdyn = "cell_dynamics='pr'," if calculation in {"vc-md", "vc-relax"} else ""
        if calculation == "vc-relax":
            cdyn = "cell_dynamics='bfgs',"
        cell = f"""
&CELL
 {cdyn}
 press=0.0,
 wmass=20.0,
 cell_dofree='all',
 {extra_cell}
 {caimd}
/
"""
    cards = f"""
ATOMIC_SPECIES
 C 12.011 {PSEUDO_NAME}
CELL_PARAMETERS bohr
{cell_text(cell_kind)}
ATOMIC_POSITIONS crystal
{positions_text(nat, fixed_flags=fixed_flags)}
K_POINTS gamma
"""
    if velocities:
        cards += f"ATOMIC_VELOCITIES {{ a.u }}\n{velocities_text(nat)}\n"
    return f"""
&CONTROL
 calculation='{calculation}',
 prefix='{prefix}',
 outdir='./tmp',
 pseudo_dir='{pseudo_dir}',
 tstress=.true.,
 tprnfor=.true.,
 nstep={nstep},
 dt=5.0,
 verbosity='high',
 {extra_control}
/
&SYSTEM
 ibrav=0,
 nat={nat},
 ntyp=1,
 ecutwfc=12,
 ecutrho=48,
 occupations='smearing',
 smearing='gaussian',
 degauss=0.01,
 {extra_system}
/
&ELECTRONS
 conv_thr=1.0d-4,
 electron_maxstep=30,
 mixing_beta=0.30,
/
{ions}
{cell}
{cards}
{after_cards}
"""


def cp_input(prefix: str, pseudo_dir: Path, *, caimd: bool) -> str:
    caimd_line = "caimd=.true., caimd_mask(1)=.true.," if caimd else "caimd=.false.,"
    return f"""
&CONTROL
 calculation='vc-cp',
 restart_mode='from_scratch',
 prefix='{prefix}',
 outdir='./tmp',
 pseudo_dir='{pseudo_dir}',
 nstep=1,
 dt=5.0,
/
&SYSTEM
 ibrav=0,
 nat=1,
 ntyp=1,
 nbnd=4,
 ecutwfc=12,
 ecutrho=48,
/
&ELECTRONS
 electron_dynamics='verlet',
 emass=400.d0,
 emass_cutoff=2.5d0,
/
&IONS
 ion_dynamics='verlet',
/
&CELL
 cell_dynamics='pr',
 wmass=20.d0,
 {caimd_line}
/
ATOMIC_SPECIES
 C 12.011 {PSEUDO_NAME}
CELL_PARAMETERS bohr
{cell_text()}
ATOMIC_POSITIONS bohr
 C 0.0 0.0 0.0
"""


def input_templates(pseudo_dir: Path) -> dict[str, str]:
    caimd_false = "caimd=.false.,"
    all_false = (
        "caimd=.true., caimd_stride=1, "
        "caimd_mask(1)=.false., caimd_mask(2)=.false., caimd_mask(3)=.false., "
        "caimd_mask(4)=.false., caimd_mask(5)=.false., caimd_mask(6)=.false.,"
    )
    fixed_xx_zero = "caimd=.true., caimd_stride=1, caimd_step_strain(1)=0.0, caimd_mask(1)=.true.,"
    stride2_xx = "caimd=.true., caimd_stride=2, caimd_step_strain(1)=1.0d-4, caimd_mask(1)=.true.,"
    shear_xy = "caimd=.true., caimd_stride=1, caimd_step_strain(6)=1.0d-4, caimd_mask(6)=.true.,"
    cumulative = "caimd=.true., caimd_stride=1, caimd_step_strain(1)=6.0d-3, caimd_mask(1)=.true.,"
    active_system = "nosym=.true., noinv=.true.,"
    return {
        "startup_scf": qe_input("startup_scf", "scf", pseudo_dir, include_cell=False),
        "scf_ignored_cell_false": qe_input("scf_ignored_cell_false", "scf", pseudo_dir, caimd=caimd_false),
        "scf_cell_true": qe_input("scf_cell_true", "scf", pseudo_dir, caimd="caimd=.true.,"),
        "scf_malformed_cell": qe_input(
            "scf_malformed_cell", "scf", pseudo_dir, caimd="caimd=.false., bogus_caimd_key=1,"
        ),
        "scf_cell_after_cards": qe_input(
            "scf_cell_after_cards",
            "scf",
            pseudo_dir,
            include_cell=False,
            after_cards="\n&CELL\n caimd=.true.,\n/\n",
        ),
        "md_fixed": qe_input(
            "md_fixed",
            "md",
            pseudo_dir,
            include_cell=False,
            velocities=True,
            extra_system=active_system,
        ),
        "vcmd_default": qe_input("vcmd_default", "vc-md", pseudo_dir, velocities=True),
        "vcmd_false": qe_input("vcmd_false", "vc-md", pseudo_dir, caimd=caimd_false, velocities=True),
        "vcmd_allfalse": qe_input("vcmd_allfalse", "vc-md", pseudo_dir, caimd=all_false, velocities=True),
        "active_fixed_xx_zero": qe_input(
            "active_fixed_xx_zero",
            "vc-md",
            pseudo_dir,
            caimd=fixed_xx_zero,
            velocities=True,
            extra_system=active_system,
        ),
        "active_stride2_xx": qe_input(
            "active_stride2_xx",
            "vc-md",
            pseudo_dir,
            nstep=6,
            caimd=stride2_xx,
            velocities=True,
            extra_system=active_system,
        ),
        "active_shear_xy": qe_input(
            "active_shear_xy",
            "vc-md",
            pseudo_dir,
            caimd=shear_xy,
            velocities=True,
            extra_system=active_system,
        ),
        "active_nonorthogonal": qe_input(
            "active_nonorthogonal",
            "vc-md",
            pseudo_dir,
            caimd=fixed_xx_zero,
            cell_kind="nonorthogonal",
            velocities=True,
            extra_system=active_system,
        ),
        "active_cumulative_warning": qe_input(
            "active_cumulative_warning",
            "vc-md",
            pseudo_dir,
            nstep=4,
            caimd=cumulative,
            velocities=True,
            extra_system=active_system,
        ),
        "reject_fixed_atom": qe_input(
            "reject_fixed_atom",
            "vc-md",
            pseudo_dir,
            caimd=fixed_xx_zero,
            velocities=True,
            fixed_flags=True,
            extra_system=active_system,
        ),
        "reject_ionic_thermostat": qe_input(
            "reject_ionic_thermostat",
            "vc-md",
            pseudo_dir,
            caimd=fixed_xx_zero,
            velocities=True,
            extra_ions="ion_temperature='rescale-v',",
            extra_system=active_system,
        ),
        "reject_cell_thermostat": qe_input(
            "reject_cell_thermostat",
            "vc-md",
            pseudo_dir,
            caimd=fixed_xx_zero,
            velocities=True,
            extra_cell="cell_temperature='nose', temph=300.0,",
            extra_system=active_system,
        ),
        "reject_non_vcmd": qe_input("reject_non_vcmd", "scf", pseudo_dir, caimd="caimd=.true.,"),
        "reject_cell_w": qe_input(
            "reject_cell_w",
            "vc-md",
            pseudo_dir,
            caimd="caimd=.true., caimd_mask(1)=.true.,",
            velocities=True,
        ).replace("cell_dynamics='pr',", "cell_dynamics='w',"),
        "reject_dofree_xy": qe_input(
            "reject_dofree_xy",
            "vc-md",
            pseudo_dir,
            caimd="caimd=.true., caimd_mask(1)=.true.,",
            velocities=True,
        ).replace("cell_dofree='all',", "cell_dofree='xy',"),
        "reject_restart": qe_input(
            "reject_restart",
            "vc-md",
            pseudo_dir,
            caimd="caimd=.true., caimd_mask(1)=.true.,",
            velocities=True,
            extra_control="restart_mode='restart',",
        ),
        "reject_stride_zero": qe_input(
            "reject_stride_zero",
            "vc-md",
            pseudo_dir,
            caimd="caimd=.true., caimd_stride=0, caimd_mask(1)=.true.,",
            velocities=True,
        ),
        "reject_event_large": qe_input(
            "reject_event_large",
            "vc-md",
            pseudo_dir,
            caimd="caimd=.true., caimd_stride=1, caimd_step_strain(1)=1.1d-2, caimd_mask(1)=.true.,",
            velocities=True,
        ),
        "cp_ordinary": cp_input("cp_ordinary", pseudo_dir, caimd=False),
        "cp_active_reject": cp_input("cp_active_reject", pseudo_dir, caimd=True),
    }


def case_specs() -> dict[str, CaseSpec]:
    specs = [
        CaseSpec("startup_scf", "A_startup_gpu", "pw", "startup_scf", monitor_gpu=True),
        CaseSpec("baseline_vcmd_default", "B_noop", "pw", "vcmd_default", executable="baseline"),
        CaseSpec("current_vcmd_default", "B_noop", "pw", "vcmd_default"),
        CaseSpec("current_vcmd_false", "B_noop", "pw", "vcmd_false", compare_to="current_vcmd_default"),
        CaseSpec("current_vcmd_allfalse", "B_noop", "pw", "vcmd_allfalse", compare_to="current_vcmd_default"),
        CaseSpec("ordinary_scf_no_cell", "C_ordinary", "pw", "startup_scf"),
        CaseSpec("ordinary_scf_ignored_cell_false", "C_ordinary", "pw", "scf_ignored_cell_false"),
        CaseSpec("ordinary_md_fixed", "C_ordinary", "pw", "md_fixed"),
        CaseSpec(
            "parser_scf_cell_true",
            "D_parser",
            "pw",
            "scf_cell_true",
            expect="failure",
            tokens=("caimd", "vc-md"),
        ),
        CaseSpec(
            "parser_scf_malformed_cell",
            "D_parser",
            "pw",
            "scf_malformed_cell",
            expect="failure",
            tokens=("error",),
        ),
        CaseSpec(
            "parser_after_card_cell",
            "D_parser",
            "pw",
            "scf_cell_after_cards",
            expect="observe",
            tokens=("cell",),
        ),
        CaseSpec("cp_ordinary", "E_cp", "cp", "cp_ordinary", executable="cp"),
        CaseSpec(
            "cp_active_reject",
            "E_cp",
            "cp",
            "cp_active_reject",
            executable="cp",
            expect="failure",
            tokens=("caimd", "cp"),
        ),
        CaseSpec("active_fixed_xx_zero", "F_active", "pw", "active_fixed_xx_zero"),
        CaseSpec("active_stride2_xx", "F_active", "pw", "active_stride2_xx"),
        CaseSpec("active_shear_xy", "F_active", "pw", "active_shear_xy"),
        CaseSpec("active_nonorthogonal", "F_active", "pw", "active_nonorthogonal"),
        CaseSpec(
            "policy_cumulative_warning",
            "G_policies",
            "pw",
            "active_cumulative_warning",
            tokens=("cumulative masked nominal strain",),
        ),
        CaseSpec(
            "policy_fixed_atom_reject",
            "G_policies",
            "pw",
            "reject_fixed_atom",
            expect="failure",
            tokens=("fixed atomic coordinates",),
        ),
        CaseSpec(
            "policy_ionic_thermostat_reject",
            "G_policies",
            "pw",
            "reject_ionic_thermostat",
            expect="failure",
            tokens=("ionic thermostats",),
        ),
        CaseSpec(
            "policy_cell_thermostat_reject",
            "G_policies",
            "pw",
            "reject_cell_thermostat",
            expect="failure",
            tokens=("cell thermostats",),
        ),
        CaseSpec("caimd_output_schema", "H_caimd_output", "pw", "active_stride2_xx"),
        CaseSpec(
            "reject_non_vcmd",
            "I_rejections",
            "pw",
            "reject_non_vcmd",
            expect="failure",
            tokens=("caimd", "vc-md"),
        ),
        CaseSpec(
            "reject_cell_w",
            "I_rejections",
            "pw",
            "reject_cell_w",
            expect="failure",
            tokens=("cell_dynamics", "pr"),
        ),
        CaseSpec(
            "reject_dofree_xy",
            "I_rejections",
            "pw",
            "reject_dofree_xy",
            expect="failure",
            tokens=("cell_dofree", "all"),
        ),
        CaseSpec(
            "reject_restart",
            "I_rejections",
            "pw",
            "reject_restart",
            expect="failure",
            tokens=("restart",),
        ),
        CaseSpec(
            "reject_stride_zero",
            "I_rejections",
            "pw",
            "reject_stride_zero",
            expect="failure",
            tokens=("stride",),
        ),
        CaseSpec(
            "reject_event_large",
            "I_rejections",
            "pw",
            "reject_event_large",
            expect="failure",
            tokens=("1.0d-2",),
        ),
        CaseSpec("rotation_continuity_nonorthogonal", "J_rotation", "pw", "active_nonorthogonal"),
    ]
    return {spec.name: spec for spec in specs}


def resolve_cases(args: argparse.Namespace, specs: dict[str, CaseSpec]) -> list[CaseSpec]:
    if args.list_cases:
        for name in sorted(specs):
            spec = specs[name]
            print(f"{name:34s} {spec.group:16s} {spec.expect:8s} {spec.description}")
        raise SystemExit(0)
    selected: set[str] = set()
    if args.groups:
        groups = {g.strip() for part in args.groups for g in part.split(",") if g.strip()}
        selected.update(name for name, spec in specs.items() if spec.group in groups)
    if args.cases:
        for part in args.cases:
            for name in part.split(","):
                name = name.strip()
                if name:
                    if name not in specs:
                        raise SystemExit(f"unknown case: {name}")
                    selected.add(name)
    if not selected:
        selected = set(specs)
    return [specs[name] for name in sorted(selected, key=lambda n: list(specs).index(n))]


def prepare_run_root(args: argparse.Namespace) -> dict[str, Path]:
    run_id = args.run_id or f"run-{utc_stamp()}"
    root = Path(args.work_root).expanduser().resolve() / run_id
    paths = {
        "root": root,
        "inputs": root / "inputs",
        "runs": root / "runs",
        "analysis": root / "analysis",
        "pseudos": root / "pseudopotentials",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def executable_for(spec: CaseSpec, args: argparse.Namespace) -> Path | None:
    if spec.executable == "baseline":
        return Path(args.baseline_exe).expanduser().resolve() if args.baseline_exe else None
    if spec.executable == "cp":
        return Path(args.cp_exe).expanduser().resolve() if args.cp_exe else None
    return Path(args.current_exe).expanduser().resolve() if args.current_exe else None


def write_inputs(paths: dict[str, Path]) -> dict[str, Path]:
    templates = input_templates(paths["pseudos"])
    written: dict[str, Path] = {}
    for name, text in templates.items():
        path = paths["inputs"] / f"{name}.in"
        path.write_text(text.strip() + "\n")
        written[name] = path
    return written


def copy_pseudo(args: argparse.Namespace, paths: dict[str, Path]) -> dict[str, Any]:
    source = Path(args.pseudo).expanduser().resolve()
    if not source.exists():
        raise SystemExit(f"pseudopotential not found; refusing to run: {source}")
    if not source.is_file():
        raise SystemExit(f"pseudopotential is not a file; refusing to run: {source}")
    target = paths["pseudos"] / PSEUDO_NAME
    shutil.copy2(source, target)
    return {
        "source": str(source),
        "copied": str(target),
        "filename": PSEUDO_NAME,
        "sha256": sha256(target),
        "policy": "local file copied; no download attempted",
    }


def build_command(
    args: argparse.Namespace,
    spec: CaseSpec,
    exe: Path,
    input_path: Path,
    run_dir: Path,
) -> tuple[list[str], Path | None]:
    monitor_path: Path | None = None
    if args.mode == "local":
        if spec.ranks != 1 and not args.mpi_launcher:
            raise RuntimeError("local mode supports ranks>1 only with --mpi-launcher")
        if spec.ranks != 1:
            return [args.mpi_launcher, "-n", str(spec.ranks), str(exe), "-inp", input_path.name], None
        return [str(exe), "-inp", input_path.name], None

    if not args.allocation:
        raise RuntimeError("Leonardo mode requires --allocation")
    loader = Path(args.loader).expanduser().resolve()
    monitor = ""
    cleanup = "status=$?"
    if spec.monitor_gpu:
        monitor_path = run_dir / "nvidia_smi.csv"
        monitor = (
            "nvidia-smi --query-gpu=timestamp,index,utilization.gpu,memory.used,power.draw "
            f"--format=csv,noheader,nounits -lms 200 > {shlex.quote(monitor_path.name)} 2>&1 & mon=$!; "
        )
        cleanup = "status=$?; kill $mon 2>/dev/null || true; wait $mon 2>/dev/null || true"
    launcher = f"{shlex.quote(str(exe))} -inp {shlex.quote(input_path.name)}"
    outer_ranks = spec.ranks
    cpus = max(1, args.cpus_per_task)
    if spec.ranks > 1:
        launcher = f"mpiexec --oversubscribe -n {spec.ranks} {shlex.quote(str(exe))} -inp {shlex.quote(input_path.name)}"
        outer_ranks = 1
        cpus = max(cpus, 8)
    inner = (
        f"source {shlex.quote(str(loader))}; export OMP_NUM_THREADS={args.omp_threads}; "
        "export OMPI_MCA_coll_hcoll_enable=0; "
        f"cd {shlex.quote(str(run_dir))}; {monitor}{launcher}; {cleanup}; exit $status"
    )
    return [
        "srun",
        f"--jobid={args.allocation}",
        "--overlap",
        "--nodes=1",
        f"--ntasks={outer_ranks}",
        f"--cpus-per-task={cpus}",
        "--mpi=pmix_v3",
        "bash",
        "-lc",
        inner,
    ], monitor_path


def enrich_result(result: CaseResult, spec: CaseSpec) -> None:
    if result.skipped or not result.output_path or not result.run_dir or not result.input_path:
        return
    output_path = Path(result.output_path)
    run_dir = Path(result.run_dir)
    prefix = Path(result.input_path).stem
    text = output_path.read_text(errors="ignore") if output_path.exists() else ""
    result.job_done = "JOB DONE" in text
    result.mpi_abort = "MPI_ABORT" in text or "MPI_ABORT was invoked" in text
    result.missing_library = "not found" in text and "lib" in text
    result.metrics["energies_Ry"] = parse_energy(text)
    for ext in ("cel", "pos", "vel", "for", "ave", "avec", "e", "eal", "p", "tv"):
        path = run_dir / f"{prefix}.{ext}"
        if path.exists():
            result.metrics[f"{ext}_path"] = str(path)
            result.metrics[f"{ext}_float_count"] = len(read_float_file(path))
    md_path = run_dir / "tmp" / f"{prefix}.md"
    if md_path.exists():
        result.metrics["md_path"] = str(md_path)
        result.metrics["md_float_count"] = len(read_float_file(md_path))
    caimd_path = run_dir / "tmp" / f"{prefix}.caimd"
    if caimd_path.exists():
        result.metrics["caimd"] = parse_caimd(caimd_path)
    cel_path = run_dir / f"{prefix}.cel"
    if cel_path.exists():
        result.metrics["rotation"] = rotation_metrics(cel_path)
    if spec.monitor_gpu:
        monitor_path = run_dir / "nvidia_smi.csv"
        result.metrics["gpu_monitor"] = parse_gpu_monitor(monitor_path)


def interpret_success(spec: CaseSpec, rc: int, out: str, timeout: bool) -> bool:
    lower = out.lower()
    tokens_ok = all(token.lower() in lower for token in spec.tokens)
    if spec.expect == "success":
        return rc == 0 and "JOB DONE" in out and not timeout and tokens_ok
    if spec.expect == "failure":
        return rc != 0 and not timeout and tokens_ok
    if spec.expect == "observe":
        return not timeout
    raise ValueError(f"unknown expectation: {spec.expect}")


def run_case(
    args: argparse.Namespace,
    spec: CaseSpec,
    inputs: dict[str, Path],
    paths: dict[str, Path],
) -> CaseResult:
    exe = executable_for(spec, args)
    if exe is None:
        return CaseResult(
            spec.name,
            spec.group,
            spec.code,
            None,
            None,
            spec.ranks,
            None,
            None,
            None,
            [],
            None,
            0.0,
            spec.expect,
            True,
            skipped=True,
            reason=f"no executable configured for {spec.executable}",
        )
    if not exe.exists():
        return CaseResult(
            spec.name,
            spec.group,
            spec.code,
            str(exe),
            None,
            spec.ranks,
            None,
            None,
            None,
            [],
            None,
            0.0,
            spec.expect,
            False,
            skipped=True,
            reason=f"executable does not exist: {exe}",
        )
    input_path = inputs[spec.input_key]
    run_dir = paths["runs"] / spec.name
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    shutil.copy2(input_path, run_dir / input_path.name)
    try:
        cmd, _monitor = build_command(args, spec, exe, run_dir / input_path.name, run_dir)
    except RuntimeError as exc:
        return CaseResult(
            spec.name,
            spec.group,
            spec.code,
            str(exe),
            sha256(exe),
            spec.ranks,
            str(input_path),
            str(run_dir),
            None,
            [],
            None,
            0.0,
            spec.expect,
            False,
            skipped=True,
            reason=str(exc),
        )
    start = time.time()
    rc, out, timeout = run_capture(cmd, cwd=run_dir if args.mode == "local" else None, timeout=args.timeout)
    elapsed = time.time() - start
    output_path = run_dir / f"{spec.name}.out"
    output_path.write_text(out)
    success = interpret_success(spec, rc, out, timeout)
    result = CaseResult(
        name=spec.name,
        group=spec.group,
        code=spec.code,
        executable=str(exe),
        executable_sha256=sha256(exe),
        ranks=spec.ranks,
        input_path=str(input_path),
        run_dir=str(run_dir),
        output_path=str(output_path),
        command=cmd,
        returncode=rc,
        elapsed_s=elapsed,
        expected=spec.expect,
        success=success,
        timeout=timeout,
    )
    enrich_result(result, spec)
    return result


def result_to_json(result: CaseResult) -> dict[str, Any]:
    return {
        "name": result.name,
        "group": result.group,
        "code": result.code,
        "executable": result.executable,
        "executable_sha256": result.executable_sha256,
        "ranks": result.ranks,
        "input_path": result.input_path,
        "run_dir": result.run_dir,
        "output_path": result.output_path,
        "command": result.command,
        "returncode": result.returncode,
        "elapsed_s": result.elapsed_s,
        "expected": result.expected,
        "success": result.success,
        "skipped": result.skipped,
        "reason": result.reason,
        "job_done": result.job_done,
        "mpi_abort": result.mpi_abort,
        "missing_library": result.missing_library,
        "timeout": result.timeout,
        "metrics": result.metrics,
    }


def write_summary(paths: dict[str, Path], summary: dict[str, Any]) -> None:
    (paths["analysis"] / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def add_comparisons(summary: dict[str, Any]) -> None:
    results = summary["cases"]
    comparisons: dict[str, Any] = {}
    if "baseline_vcmd_default" in results and "current_vcmd_default" in results:
        comparisons["baseline_vs_current_default"] = compare_result_files(
            results["baseline_vcmd_default"], results["current_vcmd_default"]
        )
    if "current_vcmd_default" in results and "current_vcmd_false" in results:
        comparisons["current_default_vs_false"] = compare_result_files(
            results["current_vcmd_default"], results["current_vcmd_false"]
        )
    if "current_vcmd_default" in results and "current_vcmd_allfalse" in results:
        comparisons["current_default_vs_allfalse"] = compare_result_files(
            results["current_vcmd_default"], results["current_vcmd_allfalse"]
        )
    summary["comparisons"] = comparisons


def compare_result_files(a: dict[str, Any], b: dict[str, Any]) -> dict[str, Any]:
    if a.get("skipped") or b.get("skipped") or not a.get("run_dir") or not b.get("run_dir"):
        return {"skipped": True}
    ar = Path(a["run_dir"])
    br = Path(b["run_dir"])
    ap = Path(a["input_path"]).stem
    bp = Path(b["input_path"]).stem
    result = {"stdout": compare_float_files(Path(a["output_path"]), Path(b["output_path"]))}
    for ext in ("cel", "pos", "vel", "for", "ave", "avec", "e", "eal", "p", "tv"):
        result[ext] = compare_float_files(ar / f"{ap}.{ext}", br / f"{bp}.{ext}")
    result["md"] = compare_float_files(ar / "tmp" / f"{ap}.md", br / "tmp" / f"{bp}.md")
    ea = a.get("metrics", {}).get("energies_Ry", [])
    eb = b.get("metrics", {}).get("energies_Ry", [])
    n = min(len(ea), len(eb))
    if n:
        diff = [ea[i] - eb[i] for i in range(n)]
        result["total_energy_Ry"] = {
            "count_a": len(ea),
            "count_b": len(eb),
            "max_abs": max(abs(x) for x in diff),
            "rms": math.sqrt(sum(x * x for x in diff) / n),
        }
    return result


def default_paths_for_mode(mode: str) -> dict[str, str]:
    if mode == "leonardo":
        return {
            "current_exe": str(DEFAULT_REMOTE_ROOT / "install-current-gpu" / "bin" / "pw.x"),
            "baseline_exe": str(DEFAULT_REMOTE_ROOT / "build-baseline-gpu" / "bin" / "pw.x"),
            "pseudo": str(DEFAULT_REMOTE_ROOT / "src-current" / "pseudo" / PSEUDO_NAME),
            "loader": str(DEFAULT_REMOTE_ROOT / "scripts" / "load_build_environment.sh"),
            "work_root": str(DEFAULT_REMOTE_ROOT / "validation"),
        }
    return {
        "current_exe": "",
        "baseline_exe": "",
        "pseudo": str(DEFAULT_LOCAL_PSEUDO),
        "loader": "",
        "work_root": "/private/tmp/qe_caimd_validation",
    }


def parse_args() -> argparse.Namespace:
    defaults = default_paths_for_mode("local")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("local", "leonardo"), default="local")
    parser.add_argument("--current-exe", default=None, help="current pw.x executable")
    parser.add_argument("--baseline-exe", default=None, help="optional baseline pw.x executable")
    parser.add_argument("--cp-exe", default=None, help="optional cp.x executable for CP boundary cases")
    parser.add_argument("--pseudo", default=None, help="local C.UPF pseudopotential path")
    parser.add_argument("--work-root", default=None, help="directory under which a timestamped run root is made")
    parser.add_argument("--run-id", default=None, help="optional deterministic run id")
    parser.add_argument("--cases", action="append", help="comma-separated case names to run")
    parser.add_argument("--groups", action="append", help="comma-separated group names to run")
    parser.add_argument("--list-cases", action="store_true")
    parser.add_argument("--allocation", help="Slurm allocation/job id for Leonardo mode")
    parser.add_argument("--loader", default=None, help="Leonardo module loader script")
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--cpus-per-task", type=int, default=8)
    parser.add_argument("--omp-threads", type=int, default=1)
    parser.add_argument("--mpi-launcher", default=None, help="local MPI launcher for ranks>1 cases")
    args = parser.parse_args()
    mode_defaults = default_paths_for_mode(args.mode)
    for key, value in mode_defaults.items():
        arg_name = key.replace("_", "-")
        if getattr(args, key, None) in (None, ""):
            setattr(args, key, value)
        if arg_name == "current-exe" and not args.current_exe:
            args.current_exe = defaults["current_exe"]
    if args.mode == "local" and not args.current_exe and not args.list_cases:
        parser.error("--current-exe is required in local mode")
    return args


def main() -> int:
    args = parse_args()
    specs = case_specs()
    selected = resolve_cases(args, specs)
    paths = prepare_run_root(args)
    pseudo = copy_pseudo(args, paths)
    inputs = write_inputs(paths)
    summary: dict[str, Any] = {
        "schema": "caimd_validation_v1",
        "created_utc": utc_stamp(),
        "mode": args.mode,
        "run_root": str(paths["root"]),
        "python": sys.version,
        "platform": platform.platform(),
        "repo_root": str(REPO_ROOT),
        "pseudopotential": pseudo,
        "cases": {},
        "comparisons": {},
        "failures": [],
        "skipped": [],
    }
    write_summary(paths, summary)
    for spec in selected:
        result = run_case(args, spec, inputs, paths)
        summary["cases"][spec.name] = result_to_json(result)
        if result.skipped:
            summary["skipped"].append(spec.name)
        elif not result.success and spec.expect != "observe":
            summary["failures"].append(spec.name)
        write_summary(paths, summary)
    add_comparisons(summary)
    for name, case in summary["cases"].items():
        if name.startswith("current_vcmd_") or name.startswith("active_") or name.startswith("caimd_"):
            caimd = case.get("metrics", {}).get("caimd")
            if caimd and caimd.get("exists") and not caimd.get("single_column_count_26"):
                summary["failures"].append(f"{name}:caimd_columns")
    write_summary(paths, summary)
    print(json.dumps({"summary": str(paths["analysis"] / "summary.json"), "failures": summary["failures"]}, indent=2))
    return 0 if not summary["failures"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
