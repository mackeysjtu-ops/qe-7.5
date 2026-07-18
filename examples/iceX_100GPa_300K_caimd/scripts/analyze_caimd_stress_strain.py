#!/usr/bin/env python3
"""Analyze CAIMD stress-strain diagnostic files."""

from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path


EXAMPLE_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = EXAMPLE_ROOT / "convergence" / "analysis"
PILOT_STRAIN_BOUND = 5.0e-3
TEMP_SPIKE_WARN_K = 50.0
VELOCITY_JUMP_RMS_WARN = 1.0e-6
COLUMNS = [
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


def parse_caimd(path: Path) -> list[dict[str, float | int | bool]]:
    rows: list[dict[str, float | int | bool]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) < len(COLUMNS):
            continue
        row: dict[str, float | int | bool] = {}
        for name, value in zip(COLUMNS, fields):
            if name in {"ionic_step", "imposed_events_done"}:
                row[name] = int(value)
            elif name == "following_move_is_event":
                row[name] = value.upper().startswith("T")
            else:
                row[name] = float(value.replace("D", "E").replace("d", "e"))
        rows.append(row)
    return rows


def case_kind(path: Path) -> str:
    parts = path.relative_to(EXAMPLE_ROOT).parts
    joined = "/".join(parts)
    if "compression" in parts:
        return "compression_xx"
    if "shear" in parts:
        return "shear_xy"
    if "tension" in parts:
        return "tension_xx"
    if "zero_strain" in parts:
        return "zero_strain"
    if "parameter_screen" in parts:
        return "parameter_screen"
    if joined.startswith("deformation/"):
        return "deformation"
    return "other"


def finite_row(row: dict[str, float | int | bool]) -> bool:
    for value in row.values():
        if isinstance(value, bool) or isinstance(value, int):
            continue
        if not math.isfinite(float(value)):
            return False
    return True


def row_for_csv(path: Path, row: dict[str, float | int | bool]) -> dict[str, object]:
    out: dict[str, object] = {"path": str(path.relative_to(EXAMPLE_ROOT)), "case_kind": case_kind(path)}
    out.update(row)
    out["engineering_gamma_xy"] = 2.0 * float(row["actual_xy"])
    out["mean_pressure_from_stress_GPa"] = -(
        float(row["stress_xx_GPa"]) + float(row["stress_yy_GPa"]) + float(row["stress_zz_GPa"])
    ) / 3.0
    return out


def collect_flags(path: Path, rows: list[dict[str, float | int | bool]]) -> list[str]:
    if not rows:
        return ["missing_rows"]

    flags: list[str] = []
    steps = [int(row["ionic_step"]) for row in rows]
    if len(set(steps)) != len(steps):
        flags.append("duplicate_ionic_steps")
    if any(not finite_row(row) for row in rows):
        flags.append("non_finite_values")

    events = [int(row["imposed_events_done"]) for row in rows]
    if any(later < earlier for earlier, later in zip(events, events[1:])):
        flags.append("non_monotonic_event_counter")

    kind = case_kind(path)
    final = rows[-1]
    final_xx = float(final["final_actual_xx"] if "final_actual_xx" in final else final["actual_xx"])
    final_xy = float(final["final_actual_xy"] if "final_actual_xy" in final else final["actual_xy"])
    final_nominal_xx = float(final["nominal_xx"])
    final_nominal_xy = float(final["nominal_xy"])
    if kind == "tension_xx" and (final_nominal_xx < -1.0e-12 or final_xx < -1.0e-12):
        flags.append("unexpected_tension_xx_sign")
    if kind == "compression_xx" and (final_nominal_xx > 1.0e-12 or final_xx > 1.0e-12):
        flags.append("unexpected_compression_xx_sign")
    if kind == "shear_xy" and (final_nominal_xy < -1.0e-12 or final_xy < -1.0e-12):
        flags.append("unexpected_shear_xy_sign")

    max_abs_strain = max(
        max(abs(float(row[name])) for name in ("nominal_xx", "nominal_yy", "nominal_zz", "nominal_yz", "nominal_xz", "nominal_xy"))
        for row in rows
    )
    if max_abs_strain > PILOT_STRAIN_BOUND + 1.0e-10:
        flags.append("cumulative_strain_above_pilot_bound")

    temps = [float(row["temperature_K"]) for row in rows]
    if any(later - earlier > TEMP_SPIKE_WARN_K for earlier, later in zip(temps, temps[1:])):
        flags.append("large_temperature_spike")

    if max(float(row["velocity_jump_rms"]) for row in rows) > VELOCITY_JUMP_RMS_WARN:
        flags.append("large_velocity_transform_residual")
    if any(float(row["volume_bohr3"]) <= 0.0 for row in rows):
        flags.append("non_positive_volume")
    return flags


def summarize(path: Path, rows: list[dict[str, float | int | bool]]) -> dict[str, object]:
    if not rows:
        return {"path": str(path.relative_to(EXAMPLE_ROOT)), "case_kind": case_kind(path), "rows": 0, "status": "empty", "flags": ["missing_rows"]}
    temps = [float(row["temperature_K"]) for row in rows]
    jumps = [float(row["velocity_jump_rms"]) for row in rows]
    pressures = [
        -(float(row["stress_xx_GPa"]) + float(row["stress_yy_GPa"]) + float(row["stress_zz_GPa"])) / 3.0
        for row in rows
    ]
    final = rows[-1]
    first = rows[0]
    flags = collect_flags(path, rows)
    return {
        "path": str(path.relative_to(EXAMPLE_ROOT)),
        "case_kind": case_kind(path),
        "rows": len(rows),
        "status": "ok" if not flags else "warn",
        "flags": flags,
        "first_step": first["ionic_step"],
        "last_step": final["ionic_step"],
        "events_done": final["imposed_events_done"],
        "initial_temperature_K": temps[0],
        "final_temperature_K": temps[-1],
        "min_temperature_K": min(temps),
        "max_temperature_K": max(temps),
        "mean_temperature_K": statistics.fmean(temps),
        "max_temperature_step_jump_K": max((later - earlier for earlier, later in zip(temps, temps[1:])), default=0.0),
        "max_velocity_jump_rms": max(jumps),
        "mean_pressure_from_stress_GPa": statistics.fmean(pressures),
        "final_mean_pressure_from_stress_GPa": pressures[-1],
        "final_nominal_xx": final["nominal_xx"],
        "final_nominal_yy": final["nominal_yy"],
        "final_nominal_zz": final["nominal_zz"],
        "final_nominal_yz": final["nominal_yz"],
        "final_nominal_xz": final["nominal_xz"],
        "final_nominal_xy": final["nominal_xy"],
        "final_actual_xx": final["actual_xx"],
        "final_actual_yy": final["actual_yy"],
        "final_actual_zz": final["actual_zz"],
        "final_actual_yz": final["actual_yz"],
        "final_actual_xz": final["actual_xz"],
        "final_actual_xy": final["actual_xy"],
        "final_engineering_gamma_xy": 2.0 * float(final["actual_xy"]),
        "final_stress_xx_GPa": final["stress_xx_GPa"],
        "final_stress_yy_GPa": final["stress_yy_GPa"],
        "final_stress_zz_GPa": final["stress_zz_GPa"],
        "final_stress_yz_GPa": final["stress_yz_GPa"],
        "final_stress_xz_GPa": final["stress_xz_GPa"],
        "final_stress_xy_GPa": final["stress_xy_GPa"],
    }


def discover_caimd_files() -> list[Path]:
    files: list[Path] = []
    for path in sorted(EXAMPLE_ROOT.glob("**/*.caimd")):
        if path.parent.name == "tmp" and (path.parent.parent / path.name).exists():
            continue
        files.append(path)
    return files


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    caimd_files = discover_caimd_files()
    parsed = [(path, parse_caimd(path)) for path in caimd_files]
    summaries = [summarize(path, rows) for path, rows in parsed]
    detail_rows = [row_for_csv(path, row) for path, rows in parsed for row in rows]
    detail_csv_path = OUTPUT_DIR / "stress_strain.csv"
    detail_json_path = OUTPUT_DIR / "stress_strain_summary.json"
    json_path = OUTPUT_DIR / "caimd_stress_strain_summary.json"
    csv_path = OUTPUT_DIR / "caimd_stress_strain_summary.csv"
    md_path = OUTPUT_DIR / "caimd_stress_strain_summary.md"
    payload = {"files": summaries}
    detail_json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    detail_fields = ["path", "case_kind", *COLUMNS, "engineering_gamma_xy", "mean_pressure_from_stress_GPa"]
    with detail_csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=detail_fields, lineterminator="\n")
        writer.writeheader()
        for row in detail_rows:
            writer.writerow(row)
    fieldnames = sorted({key for row in summaries for key in row.keys()}) if summaries else ["path", "rows", "status"]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in summaries:
            writer.writerow(row)
    with md_path.open("w", encoding="utf-8") as handle:
        handle.write("# CAIMD stress-strain summary\n\n")
        if not summaries:
            handle.write("No `.caimd` files found yet.\n")
        else:
            handle.write("| File | Rows | Events | Final nominal xx | Final nominal xy | Final T K |\n")
            handle.write("| --- | ---: | ---: | ---: | ---: | ---: |\n")
            for row in summaries:
                handle.write(
                    f"| {row['path']} | {row['rows']} | {row.get('events_done', '')} | "
                    f"{row.get('final_nominal_xx', '')} | {row.get('final_nominal_xy', '')} | "
                    f"{row.get('final_temperature_K', '')} |\n"
                )
            handle.write("\n## Flags\n\n")
            for row in summaries:
                flags = row.get("flags") or []
                handle.write(f"- `{row['path']}`: {', '.join(flags) if flags else 'none'}\n")
    print(json.dumps({"caimd_files": len(caimd_files)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
